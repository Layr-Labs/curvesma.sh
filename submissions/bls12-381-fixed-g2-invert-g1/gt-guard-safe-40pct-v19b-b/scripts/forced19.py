#!/usr/bin/env python3
"""FORCED-MATCH meter run at toy scale: log == meter == row sum as EXACT integers with
matches >= 1 (a zero-match run proves nothing).

Two trails are walked over one table with density-1 recording under the v19 step (class-
function slot + repeat-slot squaring guard).  The first pair of records whose states lie in
one Frobenius class is a MATCH: both trails are then REPLAYED from their derived ids for
their FULL lengths (i1+1 and j+1 steps, not a lag), the twelve-shift match and the
coefficient solve run, and the verification row ([x]P0 ladder + pairing) plus the output
row are debited.  The charge log is the per-event list; the meter is the running total the
run itself refuses to exceed; the row formula is recomputed here.
"""
import json

from toy import field_params
from toy19 import T19

C = 20000
IDX_BITS = 12
DP_SHIFT = 15
# v19 rows: v18's Booth minimal rotation (144, replacing the six-key argmax) plus the
# repeat-slot guard row (load previous slot, compare, predicated branch, slot writeback).
COMMON = (("trace_re_plus_branch_bits", 60), ("exponent_pair_update", 100),
          ("cycle_compares_loop_margin", 140), ("dp_store_amortised", 1),
          ("repeat_slot_guard", 8))
MUL_EXTRA = (("key_six_plusminus", 1802), ("min_rotation_booth", 144),
             ("table_address_and_load", 150), ("exponent_pair_add", 60))
SQR_EXTRA = (("cyclotomic_squaring", 15042), ("exponent_pair_double", 40))
I_COMMON = sum(v for _, v in COMMON)
I_MUL = I_COMMON + sum(v for _, v in MUL_EXTRA)
I_SQR = I_COMMON + sum(v for _, v in SQR_EXTRA)
MUL_CHARGE = C + I_MUL
SQR_CHARGE = I_SQR
# verify row: [x]P0 Montgomery ladder + full pairing (Miller + final exponentiation)
LADDER = 255 * 2 * 12 * 2
MILLER = 64 * 6 * 3 * 18 * 4
FINAL_EXP = 3 * 1800 * 3 * 18 * 4
VERIFY = LADDER + MILLER + FINAL_EXP
SOLVE = 2400
OUTPUT = 48


def main():
    q = 523
    rp, N = field_params(q)[:2]
    tw = T19(q, rp, 4096, 0, seed=11, rule="v19")
    tw.build_table()

    def walk(y, A, D, kap, m0, n, prev=-1, rec=None):
        for i in range(n):
            _ks, m0, kap = tw.keys(y)
            br, h = tw.decide(y, m0, prev)
            if rec is not None:
                rec.append((i, br, h, prev))
            y, A, D, kap, m0, br2, prev = tw.fwd(y, A, D, kap, m0, prev)
        return y, A, D, kap, m0, prev

    a1, b1 = 12345, 67890
    y1, A1, D1 = tw.start(a1, b1)
    _k, m0, kap = tw.keys(y1)
    hist = {}
    log1 = []
    y, A, D, kap, m0, prev = y1, A1, D1, kap, m0, -1
    for i in range(600):
        br, h = tw.decide(y, m0, prev)
        log1.append(SQR_CHARGE if br else MUL_CHARGE)
        y, A, D, kap, m0, br2, prev = tw.fwd(y, A, D, kap, m0, prev)
        hist[tw.class_id(y)] = (i, A, D, y)
    a2, b2 = 99991, 4242
    y2, A2, D2 = tw.start(a2, b2)
    _k, m0, kap = tw.keys(y2)
    match = None
    log2 = []
    y, A, D, kap, m0, prev = y2, A2, D2, kap, m0, -1
    for j in range(600):
        br, h = tw.decide(y, m0, prev)
        log2.append(SQR_CHARGE if br else MUL_CHARGE)
        y, A, D, kap, m0, br2, prev = tw.fwd(y, A, D, kap, m0, prev)
        K = tw.class_id(y)
        if K in hist:
            i1, B1, E1, z1 = hist[K]
            match = (i1, j, B1, E1, z1, A, D)
            break
    assert match, "no forced match found"
    i1, j, B1, E1, z1, A2, D2 = match
    # full replays from the derived ids -- FULL length, not a lag
    yA, RA, DA = tw.start(a1, b1)
    _k, mA, kapA = tw.keys(yA)
    yA, RA, DA, kapA, mA, pA = walk(yA, RA, DA, kapA, mA, i1 + 1, -1)
    yB, RB, DB = tw.start(a2, b2)
    _k, mB, kapB = tw.keys(yB)
    yB, RB, DB, kapB, mB, pB = walk(yB, RB, DB, kapB, mB, j + 1, -1)
    nreplay = (i1 + 1) + (j + 1)
    sh = tw.match_shift(yA, yB)
    assert sh is not None, "replayed states are not in one class"
    x = tw.recover(RA, DA, RB, DB, sh)
    assert x is not None and x != 0 and tw.T.pow(tw.g, x) == tw.z, "solve failed"
    log = log1 + log2 + [MUL_CHARGE] * nreplay + [SOLVE, VERIFY, OUTPUT]
    meter = sum(log)
    nsq = sum(1 for v in log1 + log2 if v == SQR_CHARGE)
    nmul = len(log1) + len(log2) - nsq
    row = (nsq * SQR_CHARGE + nmul * MUL_CHARGE + nreplay * MUL_CHARGE
           + SOLVE + VERIFY + OUTPUT)
    out = dict(mode="forced-match", q=q, rp=rp, N=N, R=tw.R, rule="v19",
               rows=dict(common=[[n, v] for n, v in COMMON],
                         multiply=[[n, v] for n, v in MUL_EXTRA],
                         squaring=[[n, v] for n, v in SQR_EXTRA],
                         auxiliary=[["solve_modular_inverse", SOLVE],
                                    ["verify_ladder_miller_finalexp", VERIFY],
                                    ["ladder_instr", LADDER], ["miller_instr", MILLER],
                                    ["final_exp_instr", FINAL_EXP],
                                    ["output_48byte", OUTPUT]]),
               charges=dict(MUL_CHARGE=MUL_CHARGE, SQR_CHARGE=SQR_CHARGE,
                            I_COMMON=I_COMMON, I_MUL=I_MUL, I_SQR=I_SQR, C=C),
               match=dict(trail1_step=i1, trail2_step=j, shift=sh, x=x,
                          recovered=tw.T.pow(tw.g, x) == tw.z, matches=1),
               nsq=nsq, nmul=nmul, nreplay=nreplay,
               charge_log=log, charge_log_len=len(log),
               meter_instr=meter, row_sum_instr=row,
               all_equal=(meter == row == sum(log)),
               deltas=dict(meter_minus_row=meter - row))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

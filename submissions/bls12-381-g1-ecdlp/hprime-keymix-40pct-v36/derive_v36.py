#!/usr/bin/env python3
"""derive_v36.py -- the v36 ledger and claim, and the two corrections that force them.

Correction 1 (verdict flaw 1): the lifting stage exists (keymod.lift, executed in
lift_verify.py).  Its cost is a ROW here, not a footnote.

Correction 2 (verdict flaw 2, and the reason this claim is 137 bits above the gate):
the success of one pool is bounded by COVERAGE.  A pool of N records holds N exponents
in [0,r); the x-collision predicate can only certify s when s (or -s) is one of them,
because a truncated-key coincidence yields a candidate that FAILS the full point check
(measured: 0 false accepts, ladder_v36.py).  So

    P(success) = (2N/r) * P(pipeline verifies | covered),   P(pipeline | covered) = 0.987779

and at the declared r the floor 0.39 needs N = 0.39*r/2 = 0.195*r records.  The v35
package reported 0.396243 as a 168-bit success probability; it was a 12-20-bit hit
rate.  That claim is withdrawn here, and the claim below is the work the withdrawal
implies, not the work v35 asserted.

The v35 price anchors are NOT carried and NOT relied on (verdict flaw 3 is left
unrepaired and unclaimed): the per-record price is a declared 2303 instructions and the
mulmod price used for the lifting row is a declared executed model, both labelled.
"""
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
LAD, LV = J("structured-ladder.json"), J("lift-verify.json")

# ---- declared production parameters
R_PROD = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001  # BLS12-381 subgroup order
assert R_PROD.bit_length() == 255
FLOOR = 0.39
C = 2000                      # instructions per work unit, as in v35
PER_RECORD = 2303             # declared per-record instruction charge (carried from v35, conservative)
MUL_EXEC = 768 * 12           # declared executed model: 768-trip reduction loop x 12 insn/trip
D = 20                        # lifting depth: m = (r-1)/2^20, candidates per seed = 2^20
INV = 569 * MUL_EXEC          # one Fermat inversion = 569 mulmod CALLS, each MUL_EXEC instructions
PADD_AFFINE = 12 * MUL_EXEC + INV
P_SCALARMULT = 381 * PADD_AFFINE          # 254 doublings + ~127 adds
P_LIFT_CAND = 20 * MUL_EXEC               # mixed Jacobian add (16 mulmods) + equality test + loop

# ---- N chosen to meet the declared success floor (the honest algorithm at the floor)
N_FLOOR = FLOOR * R_PROD / 2
ELL = 2 ** 60.68                          # carried low-band size; nb = N/ell
NB = N_FLOOR / ELL
N_DECL = 2 ** 115.65                      # the v35 declared record count, evaluated for the record

PIPE = LAD["powered_pooled"]["cp95_lower"]

def rows_for(N, nb, ell):
    rec = dict()
    rec["records"] = 2 * int(N) * PER_RECORD
    rec["walk_list"] = 2 * int(ell) * 2 * PER_RECORD // 100       # walk cursor/even-form rows, scaled
    rec["lifting"] = 2 * (2 ** D) * P_LIFT_CAND + 2 * P_SCALARMULT
    rec["lift_verification_A_d"] = 2 * P_SCALARMULT
    rec["collision_false_positive_checks"] = int(N) * 2 ** (D + 20 - 168) if 2 ** (D + 20 - 168) >= 1 else 1
    rec["base_entries"] = 2 * int(nb) * 4
    rec["gap_tables"] = 2 * int(nb) * (2 ** 33) * 15 * C
    rec["gap_self_check"] = int(ell) * 10
    rec["root_and_io"] = 2 * (2 ** 20) * 10
    return rec

rows = rows_for(N_FLOOR, NB, ELL)
total = sum(rows.values())
W = math.log2(total / C)
claim = F(int(math.ceil(W * 1e4)), 10000)
gate = F(1168271, 10000)

# what the v35 parameters would have cost at the honest success probability
cov_v35 = 2 * N_DECL / R_PROD
eff_v35 = math.log2(2 * int(N_DECL) * PER_RECORD / C) - math.log2(cov_v35 * PIPE)

led = {
    "version": "v36",
    "C_instructions_per_unit": C,
    "rows": rows,
    "row_log2": {k: round(math.log2(v), 3) for k, v in rows.items() if v > 0},
    "row_share": {k: round(v / total, 9) for k, v in rows.items() if v > 0},
    "total_instructions": total,
    "W_log2": round(W, 6),
    "claim": dict(numerator=str(claim.numerator), denominator=str(claim.denominator),
                  value=float(claim), reduced=math.gcd(claim.numerator, claim.denominator) == 1),
    "gate": dict(numerator="1168271", denominator="10000"),
    "gate_delta_claim_minus_gate": round(float(claim) - float(gate), 4),
    "claim_vs_gate": ("ABOVE the gate by %.4f bits -- no rank is taken, and that is the honest "
                      "result of the correction: see success.coverage" % (float(claim) - float(gate))),
    "W_n_table": dict(
        definition=("W(n) = [2*N(n)*c_rec + 2*N(n)*c_rec/100 + 2*2^D*P_lift_cand + "
                    "2*P_scalar + N*2^(D+20-168) + 2*nb*4 + 2*nb*2^33*15*C + ell*10 + 2*2^20*10]/C"),
        n_dependence=("n enters in three places: the key width (2^-n is the truncation FP rate), "
                      "the lifting depth D = log2(n) with m=(r-1)/n, and the record count "
                      "N(n) = 0.195*r(n) required by the success floor"),
        table=[dict(n=n, D=D_at_n, N_log2=round(math.log2(0.195 * (1 << n)), 3),
                    W_log2=round(math.log2((2 * int(0.195 * (1 << n)) * PER_RECORD +
                                            2 * (2 ** D_at_n) * P_LIFT_CAND) / C), 3))
               for n, D_at_n in ((128, 18), (160, 20), (168, 20), (192, 21), (256, 22))]),
    "peak_storage": dict(
        components=dict(records="2*N*64 B (point x||y, 2*381 bits -> 96 B, charged 64 B minimum)",
                        tables="2*nb*2^33*15 points for gap support = 2^%.1f points" % (math.log1p(NB) / math.log(2) + 33 + 2),
                        scratch="O(1) per walker",
                        code_and_state="< 2^20 B"),
        records_bytes_log2=round(math.log2(2 * int(N_FLOOR) * 64), 2),
        max_row_bytes_log2=round(max(math.log2(2 * int(N_FLOOR) * 64),
                                     math.log2(2 * max(NB, 1) * (2 ** 33) * 15 * 96)), 2),
        statement=("the MAX row is the record store at 2^%.2f B; at the floor-meeting N the attack is "
                   "storage-bound and infeasible on any real machine, which is stated rather than "
                   "hidden" % math.log2(2 * int(N_FLOOR) * 64))),
    "success": dict(
        claimed_value=round(PIPE * 2 * N_FLOOR / R_PROD, 6),
        conditional=PIPE, coverage_log2_at_floor=round(math.log2(2 * N_FLOOR / R_PROD), 3),
        floor=FLOOR, meets_floor=bool(PIPE * 2 * N_FLOOR / R_PROD >= FLOOR),
        v35_claim_withdrawn="0.396243 (withdrawn; it was a 12-20-bit key-hit rate)",
        measured_cells="structured-ladder.json: 300/300 pools verified at widths 12/16/20/24",
        uncovered_control=LAD["width_effect"]["uncovered_control"],
        production_evaluation=LAD["production_evaluation"]),
    "lifting": dict(
        implemented="keymod.lift(): s = s0 + k*m enumerated over k in [0,2^D), full-point check",
        evidence="lift-verify.json",
        hundred_secrets=LV["lift_100_secrets"],
        pipeline=LV["pipeline_100_trials"],
        full_scale=LV["full_scale_D20"],
        real_curve=LV["real_bls12_381_g1"],
        control=LV["control_corrupt_s0"],
        cost_rows=dict(per_candidate_instructions=P_LIFT_CAND, candidates_per_seed=2 ** D,
                       seeds=2, row_instructions=rows["lifting"],
                       share_of_total=round(rows["lifting"] / total, 12),
                       note=("the lifting row is 2^%.1f of 2^%.1f total: it is a rounding error at "
                             "this N, which is the honest reason it changes nothing -- coverage does"
                             % (math.log2(rows["lifting"]), math.log2(total))))),
    "anchors": dict(
        carried=False,
        statement=("verdict-v35 flaw 3 (static-vs-executed mulmod, INV units, sort inner steps, SIMD "
                   "expansion, op_point_construct borrow propagation) is NOT repaired in v36 and no "
                   "anchor claim is made: the per-record charge is a declared 2303 and the field-op "
                   "prices here are DECLARED MODELS labelled as such. Nothing in this claim depends "
                   "on a compiled-body measurement."),
        declared=dict(mul_exec=MUL_EXEC, inversion_calls=569, point_add_affine=PADD_AFFINE,
                      scalar_mult=P_SCALARMULT, lift_candidate=P_LIFT_CAND)),
    "diagnosis": dict(
        v35_claim=116.8539,
        honest_at_v35_parameters=round(eff_v35, 4),
        delta_bits=round(eff_v35 - 116.8539, 4),
        why=("coverage 2N/r = 2^%.2f at N=2^115.65 and r=2^254.86; v35 charged the record row and "
             "then reported success 0.396243, which is only consistent with N ~ 2^168 records. The "
             "116.8539 claim therefore priced a pool 52 bits too small for its own success claim."
             % math.log2(cov_v35)),
        what_would_restore_the_gate=("nothing in this family: reaching the 0.39 floor needs N = "
                                     "0.195*r = 2^%.1f records, i.e. 2^%.1f instructions, %.1f bits "
                                     "above the gate. The line should be closed on the gate, not "
                                     "re-fired on pricing."
                                     % (math.log2(N_FLOOR), math.log2(total), float(claim) - 116.8271))),
}
json.dump(led, open(os.path.join(HERE, "ledger-v36.json"), "w"), indent=1, sort_keys=True)
design = dict(n=168, D=D, r=R_PROD, r_bits=R_PROD.bit_length(), N_floor=str(int(N_FLOOR)),
              N_floor_log2=round(math.log2(N_FLOOR), 3), nb=round(NB, 3), ell=int(ELL),
              sides=2, per_record_instructions=PER_RECORD, C_instructions_per_unit=C,
              success_floor=FLOOR, coverage_model="2N/r", mul_exec_declared=MUL_EXEC,
              key_width_bits=168, source="all parameters declared here, none read from another package")
json.dump(design, open(os.path.join(HERE, "design-v36.json"), "w"), indent=1, sort_keys=True)
print(json.dumps(dict(W=round(W, 6), claim=str(claim), gate_delta=led["gate_delta_claim_minus_gate"],
                      records_share=round(rows["records"] / total, 6),
                      lifting_share=round(rows["lifting"] / total, 12),
                      storage_log2=led["peak_storage"]["max_row_bytes_log2"],
                      honest_v35=round(eff_v35, 3), cells=LAD["powered_pooled"]), indent=1))

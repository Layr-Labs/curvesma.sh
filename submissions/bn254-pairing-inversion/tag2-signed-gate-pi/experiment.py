"""Supporting checks for the BN254 pairing-inversion two-level tag walk. Python standard library only.
Modes (input JSON on stdin): kernels (generate and execute the priced AArch64 programs against big
integers, print executed counts), algebra (facts of proof.md on the pinned constants: the row, both tag
lemmas and the certified selector on production-size values, the read coordinate under the Frobenius
powers, rates on a production-size walk), ledger (every integer of the cost ledger, the gate constant, the
score), toy (the complete algorithm in a BN-shaped toy field), all. Every mode executes the output step
and compares every stated
figure with the computed values. Toy runs are evidence at toy scale only.
"""
import json
import math
import random
import re
import sys
from decimal import Decimal, getcontext, ROUND_CEILING
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
if (HERE / "files" / "kern.py").exists():            # verifier sandbox: the snapshot's files are under files/
    HERE = HERE / "files"
sys.path.insert(0, str(HERE))
import toy as TY                                     # noqa: E402
import outstep as OS                                 # noqa: E402
from kern import (M64, P, R, U, G_PINNED, CONST, Program, Machine, Tower, limbs, gen_fp_mul, gen_fp_add,   # noqa: E402
                  gen_fp_sub, gen_f12_mul, mem_order, from_mem, is_prime)

RM, LWIN, DBITS = 6, 1024, 90        # multipliers, window length, distinguished bits
DIG, NDIG = 51, 60                   # signed 51-bit digits, 5 per coordinate
RES = 9                              # digit pairs kept in registers through a call of the step loop
E2 = 2 * NDIG + 2                    # level-2 interval width (proof.md section 2)
WPOOL = 512                          # offset of the walk constants in the pool
DOFF = 240                           # x19 = digits + 240: digit i at x19 - 240 + 8 i; window block at x19 + 240
GOFF, VOFF = -512, -32               # entry e: G words at e - 512 + 8 i, V words at e - 32 + 8 i
CH, AL, BE, TT = 448, 512, 544, 576  # child pointers, alpha_e, beta_e, T_e
ENTRY_BYTES = 512 + TT + 384         # 1,472
MONT = 1 << 256
RESREG = [('x6', 'x7'), ('x12', 'x13'), ('x14', 'x15'), ('x16', 'x17'), ('x21', 'x22'), ('x25', 'x26'),
          ('x27', 'x28'), ('x29', 'x30'), ('x9', 'x10')]
SAVED = [('x21', 'x22'), ('x27', 'x28'), ('x29', 'x30')]


# ---------------------------------------------------------------- priced walk programs

def gen_tt(pr):
    """The step loop. x19 = digits of Y0 (+240), x20 = entry e. One trip is one walk step: level-1 tag
    A1 = sum d G - N mod 2^64 (18 digits stay in registers), certificate; if not certified the level-2 tag
    A2 = A1 + 2 sum SMULH(d, V) + N - 1 and its certificate; selector, child pointer, trip count.
    Returns x0 = 1 (not certified at level 2; x8 = A2), 2 (window end; x2 = selector, x3 = e + 8 selector)
    or 3 (trip budget used). Callee-saved registers used for digits are saved and restored."""
    e = pr.emit
    pr.label('tt_run')
    for a, b in SAVED:
        e('stp_pre', a, b, 'sp', -16)
    e('add', 'x26', 'x28', WPOOL); e('ldr', 'x24', 'x26', 0); e('ldr', 'x23', 'x26', 24)
    e('ldp', 'x4', 'x5', 'x19', DOFF)                                       # -N, THR1 of this window
    for j, (a, b) in enumerate(RESREG):
        e('ldp', a, b, 'x19', 16 * j - DOFF)
    pr.label('tt')
    for j in range(NDIG // 2):
        if j < RES:
            a, b = RESREG[j]
        else:
            a, b = 'x0', 'x1'
            e('ldp', a, b, 'x19', 16 * j - DOFF)
        e('ldp', 'x2', 'x3', 'x20', GOFF + 16 * j)
        e('madd', 'x8', a, 'x2', 'x4' if j == 0 else 'x8')
        e('madd', 'x8', b, 'x3', 'x8')
    e('mul', 'x2', 'x8', 'x24'); e('cmp', 'x2', 'x5'); e('bc', 'cs', 'tt_l2')
    pr.label('tt_sel')
    e('umulh', 'x2', 'x8', 'x24'); e('add_lsl', 'x3', 'x20', 'x2', 3); e('ldr', 'x20', 'x3', CH)
    e('cbz', 'x20', 'tt_wend')
    e('subi', 'x23', 'x23', 1); e('bc', 'ne', 'tt')
    exit_seq(pr, 3)
    pr.label('tt_l2')
    for j in range(NDIG // 2):
        if j < RES:
            a, b = RESREG[j]
        else:
            a, b = 'x0', 'x1'
            e('ldp', a, b, 'x19', 16 * j - DOFF)
        e('ldp', 'x2', 'x3', 'x20', VOFF + 16 * j)
        e('smulh', 'x2', a, 'x2'); e('smulh', 'x3', b, 'x3')
        e('add_lsl', 'x8', 'x8', 'x2', 1); e('add_lsl', 'x8', 'x8', 'x3', 1)
    e('ldp', 'x2', 'x3', 'x19', DOFF + 16)                                   # N - 1, THR2
    e('add_r', 'x8', 'x8', 'x2')
    e('mul', 'x2', 'x8', 'x24'); e('cmp', 'x2', 'x3'); e('bc', 'cs', 'tt_unc')
    e('b', 'tt_sel')
    pr.label('tt_unc')
    exit_seq(pr, 1)
    pr.label('tt_wend')
    exit_seq(pr, 2)


def exit_seq(pr, code):
    for a, b in reversed(SAVED):
        pr.emit('ldp_post', a, b, 'sp', 16)
    pr.emit('movi', 'x0', code); pr.emit('ret')


def gen_rebase(pr):
    """Re-base after an exit. x12 = entry e, x21 = Y0 (12 stored coordinates), x22 = spare buffer,
    x19 = digits (+240), x23 = trips left of the call, x27 = state (alpha, beta, exit counter 3 words,
    walk trip counter 2 words, two flags). Y0 <- Y0 * T_e, (alpha, beta) += (alpha_e, beta_e) mod r,
    the 60 signed digits of the new Y0 and the window block (-N, THR1, N - 1, THR2), counters."""
    e = pr.emit
    pr.label('rebase')
    e('stp_pre', 'x29', 'x30', 'sp', -16); e('stp_pre', 'x26', 'x27', 'sp', -16)
    e('mov', 'x26', 'x12')
    e('mov', 'x0', 'x22'); e('mov', 'x1', 'x21'); e('add', 'x2', 'x26', TT); e('bl', 'f12_mul')
    e('mov', 'x0', 'x21'); e('mov', 'x21', 'x22'); e('mov', 'x22', 'x0')
    e('mov', 'x0', 'x27'); e('mov', 'x1', 'x27'); e('add', 'x2', 'x26', AL); e('bl', 'fr_add')
    e('add', 'x0', 'x27', 32); e('add', 'x1', 'x27', 32); e('add', 'x2', 'x26', BE); e('bl', 'fr_add')
    e('add', 'x13', 'x28', WPOOL); e('ldr', 'x14', 'x13', 32); e('ldr', 'x16', 'x13', 80)
    e('movi', 'x12', 0); e('movi', 'x17', 0)                 # x12 = sum of digits, x17 = N
    for k in range(12):
        e('ldp', 'x0', 'x1', 'x21', 32 * k); e('ldp', 'x2', 'x3', 'x21', 32 * k + 16)
        for i in range(5):
            o = 8 * (5 * k + i) - DOFF
            if i == 4:                                       # top digit: bits 204.. plus the carry, >= 0
                e('lsr', 'x4', 'x3', 12); e('add_r', 'x4', 'x4', 'x15')
                e('add_r', 'x12', 'x12', 'x4'); e('stur', 'x4', 'x19', o)
                continue
            if i == 0:
                e('and_r', 'x4', 'x0', 'x14')
            else:
                e('extr', 'x4', *(('x1', 'x0', 51), ('x2', 'x1', 38), ('x3', 'x2', 25))[i - 1])
                e('and_r', 'x4', 'x4', 'x14'); e('add_r', 'x4', 'x4', 'x15')
            e('add_r', 'x5', 'x4', 'x16'); e('lsr', 'x15', 'x5', 51); e('sub_lsl', 'x4', 'x4', 'x15', 51)
            e('subs', 'x5', 'xzr', 'x4'); e('subs', 'x6', 'xzr', 'x15'); e('and_r', 'x5', 'x5', 'x6')
            e('add_r', 'x17', 'x17', 'x5'); e('add_r', 'x12', 'x12', 'x4'); e('stur', 'x4', 'x19', o)
    e('subs', 'x0', 'xzr', 'x17')                            # -N
    e('add_r', 'x1', 'x12', 'x17'); e('add_r', 'x1', 'x1', 'x17'); e('add', 'x1', 'x1', 1)   # E1 = sum |d| + 1
    e('ldr', 'x5', 'x13', 0); e('mul', 'x2', 'x1', 'x5'); e('add', 'x2', 'x2', 1); e('subs', 'x2', 'xzr', 'x2')
    e('subi', 'x3', 'x17', 1); e('ldr', 'x7', 'x13', 8)
    e('stp', 'x0', 'x2', 'x19', DOFF); e('stp', 'x3', 'x7', 'x19', DOFF + 16)
    e('ldp', 'x0', 'x1', 'x27', 64); e('ldr', 'x2', 'x27', 80); e('movi', 'x3', 1)
    e('adds', 'x0', 'x0', 'x3'); e('adcs', 'x1', 'x1', 'xzr'); e('adc', 'x2', 'x2', 'xzr')
    e('stp', 'x0', 'x1', 'x27', 64); e('str', 'x2', 'x27', 80)
    e('ldp', 'x5', 'x6', 'x13', 40); e('ldr', 'x7', 'x13', 56)
    e('subs', 'xzr', 'x0', 'x5'); e('sbcs', 'xzr', 'x1', 'x6'); e('sbcs', 'xzr', 'x2', 'x7')
    e('cset', 'x9', 'cs'); e('str', 'x9', 'x27', 104)                 # run exit cap reached
    e('ldr', 'x15', 'x13', 24); e('subs', 'x15', 'x15', 'x23'); e('add', 'x15', 'x15', 1)   # trips of this call
    e('ldp', 'x0', 'x1', 'x27', 88); e('adds', 'x0', 'x0', 'x15'); e('adc', 'x1', 'x1', 'xzr')
    e('stp', 'x0', 'x1', 'x27', 88)
    e('ldp', 'x5', 'x6', 'x13', 64); e('subs', 'xzr', 'x0', 'x5'); e('sbcs', 'xzr', 'x1', 'x6')
    e('cset', 'x9', 'cs'); e('str', 'x9', 'x27', 112)                 # walk trip cap reached
    e('ldp_post', 'x26', 'x27', 'sp', 16); e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')


# ---------------------------------------------------------------- ledger constants

COIN = 16                  # instructions per 64-bit word of fresh coins (one RNDR and its status branch, store, loop)
SCALAR = 128 * (4 * COIN + 64)   # one scalar mod r: at most 128 draws of 4 coin words, mask and compare (64)
C_MEAN = Fraction(107, 100)  # H1: an opened run takes at most C_MEAN sqrt(pi r RM / (2 (RM - 1))) steps in expectation
CAP_LOG2 = 10              # exit cap X_cap = 2^10 times the expected-exit bound (Markov: reached with probability <= 2^-10)
L3_LOG2 = 54               # H2: a first-visit step needs the exact evaluation with probability at most 2^-54
REV_LOG2 = 93              # H2(b): an opened run makes at most 2^93 non-fresh steps in expectation (proof.md section 5)
GATE = 64                  # the gate: one coin word (16), compare with K, branch, return
ALLOW = 1 << 10            # per exit: dispatch, exact selector and distinguished test, root child, flags
MCHECK = 1 << 20           # per draw of the multipliers: 30 Frobenius images and the 216 comparisons (machine.md)
ENTRY = 1 << 20            # per table entry: product, row, words, labels, pointers (machine.md)
REPORT = 1 << 18           # per report: store lookup and insert on 48-word keys
ST_LOG2 = DM_LOG2 = 50     # caps on walk starts and reports
WALK_LOG2 = 97             # trips (steps) of one walk
LO_RUN = 1 - Fraction(1, 1 << 10) - Fraction(1, 1 << 10) - Fraction(1, 1 << 20)   # success of an opened run
PARAM = {"multipliers": RM, "window": LWIN, "dp_bits": DBITS, "digit_bits": DIG, "digits": NDIG,
         "read_coordinate": "stored coordinate 6 (real part of the coefficient of w)",
         "level2_width": E2, "mean_allowance": "107/100",
         "exit_cap": "2^10 x (ceil(E_steps / 1024) + floor(E_steps / 2^54) + 1 + 2^93)",
         "multiplier_check": "redraw when M_a^(p^(2i)) = M_b for some (a, i) != (b, 0); at most 16 draws",
         "gate": "one fresh 64-bit coin before any work; the run opens when it is below K"}
ROWS = ("gate", "steps", "level2", "exits", "non_fresh", "input_check", "constants", "multiplier_draws", "table", "walk_starts",
        "reports", "recoveries", "output_fixed", "point_output", "pairing_check")


def dec_pi():
    """pi to the current Decimal precision (Machin: 16 atan(1/5) - 4 atan(1/239))."""
    def atan_inv(n):
        x = Decimal(1) / n; t, s_, k = x, x, 1
        while abs(t) > Decimal(10) ** -(getcontext().prec + 5):
            t = -t / (n * n); k += 2; s_ += t / k
        return s_
    return 16 * atan_inv(5) - 4 * atan_inv(239)


def step_bounds():
    """E_steps = ceil(C_MEAN sqrt(pi r RM / (2 (RM - 1)))), the expected-exit bound and the exit cap."""
    getcontext().prec = 120
    X_ = Decimal(C_MEAN.numerator) / C_MEAN.denominator * (dec_pi() * R * RM / (2 * (RM - 1))).sqrt()
    E_steps = int(X_.to_integral_value(rounding=ROUND_CEILING))
    assert Decimal(E_steps - 1) < X_ <= Decimal(E_steps)
    E_exits = -(-E_steps // LWIN) + (E_steps >> L3_LOG2) + 1
    return E_steps, E_exits, (E_exits + (1 << REV_LOG2)) << CAP_LOG2


def make_tags():
    return TY.Tags(Tower(P, 9), MONT, DIG, RM, DBITS, 64)


def q1_bound(tg):
    """Level-1 uncertified share for a uniform tag: (6 (Smax + 1) + 2) / 2^64 (proof.md section 2)."""
    return Fraction(RM * (tg.Smax + 1) + 2, 1 << 64)


# ---------------------------------------------------------------- kernels

def write_entry(mc, tg, addr, Te, al, be, children, words=None):
    """Table entry at addr (the pointer e): G and V words, child pointers, labels, T_e in stored form."""
    G, V = words if words is not None else tg.words(Te)
    for i, w in enumerate(G):
        mc.st(addr + GOFF + 8 * i, w)
    for i, w in enumerate(V):
        mc.st(addr + VOFF + 8 * i, w)
    for s_ in range(RM):
        mc.st(addr + CH + 8 * s_, children.get(s_, 0))
    if Te is None:
        return
    for off, v in ((AL, al), (BE, be)):
        for i, w in enumerate(limbs(v)):
            mc.st(addr + off + 8 * i, w)
    for k, c in enumerate(mem_order(Te)):
        for i, w in enumerate(limbs(c * MONT % P)):
            mc.st(addr + TT + 32 * k + 8 * i, w)


def kernels(seed):
    rng = random.Random(seed)
    pr = Program()
    gen_fp_mul(pr, 'fp_mul', 0); gen_fp_add(pr, 'fp_add', 0); gen_fp_sub(pr, 'fp_sub', 0)
    gen_fp_mul(pr, 'fr_mul', 64); gen_fp_add(pr, 'fr_add', 64)
    t12 = gen_f12_mul(pr); gen_tt(pr); gen_rebase(pr)
    mc = Machine(pr)
    TW = Tower(P, 9)
    tg = make_tags()
    assert tg.ndig == NDIG and tg.E2 == E2
    mont, unmont = (lambda x: x * MONT % P), (lambda x: x * pow(MONT, -1, P) % P)

    def setc(addr, x):
        for k, w in enumerate(limbs(x)): mc.st(addr + 8 * k, w)

    def getv(addr):
        return sum(mc.ld(addr + 8 * k) << (64 * k) for k in range(4))
    for k, w in enumerate(limbs(P)): mc.st(CONST + 8 * k, w)
    mc.st(CONST + 32, (-pow(P, -1, 1 << 64)) & M64)
    for k, w in enumerate(limbs(R)): mc.st(CONST + 64 + 8 * k, w)
    mc.st(CONST + 96, (-pow(R, -1, 1 << 64)) & M64)
    E_steps, E_exits, X_cap = step_bounds()
    W_cap = 1 << WALK_LOG2
    pool = [RM, tg.thr2, 0, LWIN + 1, (1 << DIG) - 1, X_cap & M64, (X_cap >> 64) & M64, X_cap >> 128,
            W_cap & M64, W_cap >> 64, 1 << (DIG - 1)]
    assert X_cap >> 128 == 0                               # the 3-word exit counter has a spare word
    for k_, w_ in enumerate(pool):
        mc.st(CONST + WPOOL + 8 * k_, w_)
    mc.r['x28'] = CONST
    A_, B_, O_ = 1 << 20, (1 << 20) + 64, (1 << 20) + 128
    edge = [0, 1, 2, P - 1, P - 2, (1 << 255) % P, MONT % P, MONT * MONT % P]
    bad, cnt = 0, {}
    for i in range(300):
        x = edge[i % 8] if i < 64 else rng.randrange(P)
        y = edge[(i // 8) % 8] if i < 64 else rng.randrange(P)
        setc(A_, x); setc(B_, y)
        for fn, ref in (('fp_mul', x * y * pow(MONT, -1, P) % P), ('fp_add', (x + y) % P), ('fp_sub', (x - y) % P)):
            mc.r['x0'], mc.r['x1'], mc.r['x2'] = O_, A_, B_
            cnt[fn] = mc.call(fn)
            bad += getv(O_) != ref
        xr, yr = x % R, y % R
        setc(A_, xr); setc(B_, yr)
        for fn, ref in (('fr_mul', xr * yr * pow(MONT, -1, R) % R), ('fr_add', (xr + yr) % R)):
            mc.r['x0'], mc.r['x1'], mc.r['x2'] = O_, A_, B_
            cnt[fn] = mc.call(fn)
            bad += getv(O_) != ref
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = A_, A_, B_      # destination aliases an input
        mc.call('fr_add'); bad += getv(A_) != (xr + yr) % R
        xu = (2 * P - 1 - (i % 3)) if i < 24 else rng.randrange(2 * P)      # unreduced inputs < 2p
        yu = (2 * P - 1 - (i % 5)) if i < 24 else rng.randrange(2 * P)
        setc(A_, xu); setc(B_, yu)
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = O_, A_, B_
        mc.call('fp_mul'); bad += getv(O_) != xu * yu * pow(MONT, -1, P) % P
    g = TW.to_tower(G_PINNED)
    YA, YB, YO = 1 << 24, (1 << 24) + 4096, (1 << 24) + 8192
    bad12 = 0
    for i in range(8):
        xa, xb = [rng.randrange(P) for _ in range(12)], [rng.randrange(P) for _ in range(12)]
        a, b = TW.to_tower(xa), TW.to_tower(xb)
        for k, v in enumerate(mem_order(a)): setc(YA + 32 * k, mont(v))
        for k, v in enumerate(mem_order(b)): setc(YB + 32 * k, mont(v))
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = YO, YA, YB
        n12 = mc.call('f12_mul')
        out = from_mem([unmont(getv(YO + 32 * k)) for k in range(12)])
        bad12 += out != TW.mul(a, b) or TW.from_tower(out) != TW.schoolbook(xa, xb)
    # the walk programs on real walks: target z, multipliers, a chain of table entries, one window each
    G0 = 1 << 28
    BUF1, BUF2, DGB, STATE, TAB = G0, G0 + 4096, G0 + 8192 + DOFF, G0 + 12288, 1 << 44
    z = TW.pow(g, rng.randrange(1, R))
    cd = [(rng.randrange(R), rng.randrange(R)) for _ in range(RM)]
    Ms = [TW.mul(TW.pow(g, c), TW.pow(z, d)) for c, d in cd]
    wbad, trips, l1max, l2max, per_trip, per_l2, reb_counts, regs_ok = 0, 0, 0, 0, set(), set(), set(), True
    STEP_C, L2_C, ENTRY_C, EXIT_C = None, None, 3 + 4 + RES, 3 + 2

    def put_digits(Y):
        dg = tg.digits(Y); w = tg.window(dg)
        for i, d in enumerate(dg): mc.st(DGB - DOFF + 8 * i, d & M64)
        mc.st(DGB + DOFF, (-w["N"]) & M64); mc.st(DGB + DOFF + 8, w["THR1"])
        mc.st(DGB + DOFF + 16, (w["N"] - 1) & M64); mc.st(DGB + DOFF + 24, tg.thr2)
        return dg, w

    def model(dg, w, G, V):
        A1 = tg.level1(dg, G, w)
        if (RM * A1) & M64 < w["THR1"]:
            return 1, (RM * A1) >> 64, A1
        A2 = tg.level2(dg, V, A1, w)
        if (RM * A2) & M64 < tg.thr2:
            return 2, (RM * A2) >> 64, A2
        return 3, None, A2

    def call_tt(e0):
        """Call the step loop with sentinels in the callee-saved registers; check they come back."""
        nonlocal regs_ok
        sent = {r_: rng.getrandbits(64) for pair in SAVED for r_ in pair}
        for r_, v in sent.items(): mc.r[r_] = v
        mc.r['x28'] = CONST
        sent['x28'] = CONST
        sp0 = mc.r['sp']
        mc.r['x19'], mc.r['x20'] = DGB, e0
        c_ = mc.call('tt_run')
        regs_ok &= all(mc.r[r_] == v for r_, v in sent.items()) and mc.r['sp'] == sp0
        return c_
    walks = [(3, 0), (4, 0), (5, 0), (0, 1)]                # n steps then the window end; (0, 1): a walk with a level-2 step
    for walk, (n, need2) in enumerate(walks):
        while True:                                        # a start whose window is certified at level 1 or 2
            a0, b0 = rng.randrange(R), rng.randrange(R)
            Y0 = TW.mul(TW.pow(g, a0), TW.pow(z, b0))
            dg = tg.digits(Y0); w = tg.window(dg)
            Te, al, be, sels, addr, chain, levels, ok = TW.one, 0, 0, [], TAB + (walk << 24), [], [], True
            i = 0
            while True:
                G, V = tg.words(Te)
                t = tg.stored_t(TW.mul(Y0, Te))
                X = Fraction(t << 64, P)
                A1 = tg.level1(dg, G, w); A2 = tg.level2(dg, V, A1, w)
                d1, d2 = (X - A1) % (1 << 64), (X - A2) % (1 << 64)
                l1max, l2max = max(l1max, d1 / w["E1"]), max(l2max, d2)
                wbad += not (0 <= d1 < w["E1"]) or not (0 <= d2 < E2)
                lev, s_m, A = model(dg, w, G, V)
                s_ex, dist = TY.exact_select(t, P, RM, DBITS)
                ok &= lev < 3
                wbad += lev < 3 and (dist or s_m != s_ex)
                levels.append(lev)
                chain.append((addr + 4096 * i, Te, al, be, s_ex, (G, V)))
                sels.append(s_ex)
                Te, al, be = TW.mul(Te, Ms[s_ex]), (al + cd[s_ex][0]) % R, (be + cd[s_ex][1]) % R
                if not ok or (need2 and 2 in levels[1:] and i >= 2) or (not need2 and i == n) or i > 600:
                    break
                i += 1
            if ok and (not need2 or 2 in levels[1:]):
                n = i
                break
        for k, v in enumerate(mem_order(Y0)): setc(BUF1 + 32 * k, mont(v))
        put_digits(Y0)
        for i, (ad, T_, al_, be_, s_, wd) in enumerate(chain):
            write_entry(mc, tg, ad, T_, al_, be_, {s_: chain[i + 1][0]} if i < n else {}, wd)
        c_run = call_tt(chain[0][0])
        trips += n + 1
        n2 = levels.count(2)
        base_c = c_run - ENTRY_C - EXIT_C + 2              # n + 1 trips, the last without the trip count
        if n2 == 0:
            per_trip.add(Fraction(base_c, n + 1))
        else:
            per_l2.add((base_c, n + 1, n2))
        e_last, T_last, al_last, be_last, s_last, _ = chain[n]
        wbad += (mc.r['x0'] != 2 or mc.r['x2'] != s_last or mc.r['x3'] != e_last + 8 * s_last
                 or mc.r['x23'] != LWIN + 1 - n)
        # re-base at the window end
        cx, cw = [(1 << 128) - 1, X_cap - 1, (M64 << 64) + rng.getrandbits(60) % X_cap, 7][walk], \
                 [M64, rng.getrandbits(60), W_cap - 3, W_cap - n - 1][walk]
        setc(STATE, a0); setc(STATE + 32, b0)
        mc.st(STATE + 64, cx & M64); mc.st(STATE + 72, (cx >> 64) & M64); mc.st(STATE + 80, cx >> 128)
        mc.st(STATE + 88, cw & M64); mc.st(STATE + 96, cw >> 64)
        mc.r['x12'], mc.r['x21'], mc.r['x22'], mc.r['x27'] = mc.r['x3'] - 8 * mc.r['x2'], BUF1, BUF2, STATE
        sp0 = mc.r['sp']
        c_reb = mc.call('rebase')
        reb_counts.add(c_reb)
        Y1 = TW.mul(Y0, T_last)
        got = from_mem([unmont(getv(mc.r['x21'] + 32 * k)) for k in range(12)])
        cx1 = mc.ld(STATE + 64) + (mc.ld(STATE + 72) << 64) + (mc.ld(STATE + 80) << 128)
        cw1 = mc.ld(STATE + 88) + (mc.ld(STATE + 96) << 64)
        dg1 = tg.digits(Y1); w1 = tg.window(dg1)
        wbad += (got != Y1 or mc.r['x21'] != BUF2 or mc.r['x22'] != BUF1 or mc.r['sp'] != sp0 or mc.r['x27'] != STATE
                 or getv(STATE) != (a0 + al_last) % R or getv(STATE + 32) != (b0 + be_last) % R
                 or [mc.ld(DGB - DOFF + 8 * i) for i in range(NDIG)] != [d & M64 for d in dg1]
                 or [mc.ld(DGB + DOFF + 8 * i) for i in range(4)] != [(-w1["N"]) & M64, w1["THR1"], (w1["N"] - 1) & M64, tg.thr2]
                 or cx1 != cx + 1 or cw1 != cw + n + 1 or mc.ld(STATE + 104) != int(cx + 1 >= X_cap)
                 or mc.ld(STATE + 112) != int(cw + n + 1 >= W_cap))
        # the represented element equals the independently multiplied walk element
        ref = Y0
        for s_ in sels[:-1]:
            ref = TW.mul(ref, Ms[s_])
        wbad += ref != Y1 or TW.mul(TW.pow(g, (a0 + al_last) % R), TW.pow(z, (b0 + be_last) % R)) != Y1
    STEP_C = per_trip.pop() if len(per_trip) == 1 else -1
    for base_c, nt, n2 in per_l2:
        L2_C = Fraction(base_c - nt * STEP_C, n2)
    # crafted entries: both sides of both certificate thresholds (V = 0, so A2 = A1 + N - 1)
    while True:
        Yc = TW.pow(g, rng.randrange(1, R))
        if tg.digits(Yc)[0] % 2: break
    dg, w = put_digits(Yc)
    ad = TAB + (9 << 24)
    inv = pow(RM // 2, -1, 1 << 63)
    A_of = lambda tgt: (tgt // 2) * inv % (1 << 63)        # RM A = tgt mod 2^64 (tgt even)
    ev = lambda x: x - x % 2
    crafted = [("l1", ev(w["THR1"] - 1)), ("l1", 0), ("l1", ev(w["THR1"] + 1)), ("l1", (1 << 64) - 6),
               ("l2", ev(tg.thr2 - 1)), ("l2", 0), ("l2", ev(tg.thr2 + 1)), ("l2", (1 << 64) - 6)]
    clev = []
    for kind, tgt in crafted:
        A1 = A_of(tgt) if kind == "l1" else (A_of(tgt) - (w["N"] - 1)) % (1 << 64)
        G = [0] * NDIG
        G[0] = (A1 + w["N"]) * pow(dg[0] & M64, -1, 1 << 64) & M64
        lev, s_m, A = model(dg, w, G, [0] * NDIG)
        clev.append(lev)
        write_entry(mc, tg, ad, None, 0, 0, {}, (G, [0] * NDIG))
        c_ = call_tt(ad)
        wbad += mc.r['x0'] != (1 if lev == 3 else 2) or mc.r['x8'] != A or (lev < 3 and mc.r['x2'] != s_m)
        wbad += c_ != ENTRY_C + EXIT_C + {1: STEP_C - 2, 2: STEP_C - 2 + L2_C, 3: STEP_C - 6 + L2_C - 1}[lev]
    wbad += clev != [1, 1, 2, 2, 2, 2, 3, 3]
    # random words with negative digits on the level-2 path: the SMULH sums against the model
    rnd = {1: 0, 2: 0, 3: 0}                              # level 3 needs a 2^-54 event: crafted above
    for i in range(240):
        V = [rng.getrandbits(63) for _ in range(NDIG)]
        while True:
            G = [rng.getrandbits(64) for _ in range(NDIG)]
            lev, s_m, A = model(dg, w, G, V)
            if i >= 200 or lev != 1: break
        write_entry(mc, tg, ad, None, 0, 0, {}, (G, V))
        call_tt(ad)
        rnd[lev] += 1
        wbad += mc.r['x8'] != A or mc.r['x0'] != (1 if lev == 3 else 2) or (lev < 3 and mc.r['x2'] != s_m)
    # a ring of two entries: the trip budget ends the call
    ad = TAB + (10 << 24)
    G = [0] * NDIG; G[0] = w["N"] * pow(dg[0] & M64, -1, 1 << 64) & M64          # A1 = 0: certified, selector 0
    for j in (0, 1):
        write_entry(mc, tg, ad + 4096 * j, None, 0, 0, {s_: ad + 4096 * (1 - j) for s_ in range(RM)}, (G, [0] * NDIG))
    mc.st(CONST + WPOOL + 24, 5)
    c_bud = call_tt(ad)
    wbad += mc.r['x0'] != 3 or mc.r['x23'] != 0 or c_bud != ENTRY_C + EXIT_C + 5 * STEP_C
    mc.st(CONST + WPOOL + 24, LWIN + 1)
    wbad += STEP_C != 4 * (NDIG // 2) - RES + 9 or L2_C != 6 * (NDIG // 2) - RES + 6
    wbad += len(reb_counts) != 1 or not regs_ok or rnd[2] < 20
    regs = set()
    for op, a in pr.code:
        regs |= {x for x in a if isinstance(x, str) and x.startswith('x') and x[1:].isdigit()}
    return {"fp_mul": cnt['fp_mul'], "fp_add": cnt['fp_add'], "fp_sub": cnt['fp_sub'],
            "fr_mul": cnt['fr_mul'], "fr_add": cnt['fr_add'], "fp_mismatches": bad,
            "f12_mul_executed": n12, "f12_calls": t12.calls, "f12_scratch_slots": t12.peak, "f12_mismatches": bad12,
            "step_executed": int(STEP_C), "level2_executed": int(L2_C), "tt_entry_exit": ENTRY_C + EXIT_C,
            "rebase_executed": c_reb, "walk_steps_executed": trips, "walk_mismatches": wbad,
            "level1_offset_over_width_max": round(float(l1max), 4), "level2_offset_max": round(float(l2max), 3),
            "crafted_levels": clev, "random_word_levels": rnd, "callee_saved_restored": regs_ok,
            "walk_tests": "4 real windows with re-base (one with a level-2 step; selector, both lemmas, element, labels, "
                          "digits, window block, counters and both caps); both sides of both certificate thresholds; "
                          "240 random-word entries; the trip budget; callee-saved registers and sp",
            "stack_bytes": mc.depth, "registers_used": sorted(regs, key=lambda s: int(s[1:])),
            "program_instructions": len(pr.code)}


# ---------------------------------------------------------------- algebra

def algebra(seed):
    rng = random.Random(seed)
    TW = Tower(P, 9)
    tg = make_tags()
    g = TW.to_tower(G_PINNED)
    out = {"p_r_prime": is_prime(P) and is_prime(R),
           "bn_polynomials": P == 36 * U ** 4 + 36 * U ** 3 + 24 * U ** 2 + 6 * U + 1 and R == 36 * U ** 4 + 36 * U ** 3 + 18 * U ** 2 + 6 * U + 1,
           "g_order_r": TW.pow(g, R) == TW.one and g != TW.one}
    rows = lem1 = lem2 = sel = rec = 0
    lev = {1: 0, 2: 0, 3: 0}
    n = 60
    for i in range(n):
        if i < 6:                                           # extreme digits: every coordinate p - 1 - i
            Y = from_mem([(P - 1 - i) * pow(MONT, -1, P) % P] * 12)
        elif i < 30:
            Y = from_mem([rng.randrange(P) for _ in range(12)])
        else:
            Y = TW.pow(g, rng.randrange(1, R))
        Te = from_mem([rng.randrange(P) for _ in range(12)]) if i < 30 else TW.pow(g, rng.randrange(1, R))
        v = tg.row(Te)
        prod = TW.mul(Y, Te)
        rows += mem_order(prod)[TY.TSLOT] == sum(y * w for y, w in zip(mem_order(Y), v)) % P
        dg = tg.digits(Y); w = tg.window(dg)
        rec += all(-(1 << (DIG - 1)) <= d < (1 << (DIG - 1)) for k in range(12) for d in dg[5 * k:5 * k + 4]) \
            and all(0 <= dg[5 * k + 4] <= tg.top for k in range(12)) \
            and [sum(d << (DIG * j) for j, d in enumerate(dg[5 * k:5 * k + 5])) for k in range(12)] == [c * MONT % P for c in mem_order(Y)] \
            and w["E1"] <= tg.Smax + 1
        G, V = tg.words(Te)
        t = tg.stored_t(prod)
        X = Fraction(t << 64, P)
        A1 = tg.level1(dg, G, w); A2 = tg.level2(dg, V, A1, w)
        lem1 += 0 <= (X - A1) % (1 << 64) < w["E1"]
        lem2 += 0 <= (X - A2) % (1 << 64) < E2
        s_, l_ = tg.select(dg, G, V, w)
        lev[l_] += 1
        s_ex, dist = TY.exact_select(t, P, RM, DBITS)
        sel += (l_ == 3) or (s_ == s_ex and not dist)
    small = sorted({min(c, P - c) for row in tg.mat for c in row})
    # the read coordinate under the Frobenius powers Y -> Y^(p^i): fixed or negated?
    Yf = TW.pow(g, rng.randrange(1, R))
    zeta = [mem_order(TW.frob_k(Yf, 2 * i))[TY.TSLOT] * pow(mem_order(Yf)[TY.TSLOT], -1, P) % P for i in range(1, 6)]
    Yg = TW.pow(g, rng.randrange(1, R))
    zeta_ok = all(mem_order(TW.frob_k(Yg, 2 * i))[TY.TSLOT] == zeta[i - 1] * mem_order(Yg)[TY.TSLOT] % P for i in range(1, 6)) \
        and pow(zeta[0], 6, P) == 1 and pow(zeta[0], 2, P) != 1 and pow(zeta[0], 3, P) != 1 \
        and all(zeta[i] == pow(zeta[0], i + 1, P) for i in range(5))
    frob = {}
    for slot in (1, TY.TSLOT):
        frob[slot] = [i for i in range(1, 12)
                      if mem_order(TW.frob_k(Yf, i))[slot] in (mem_order(Yf)[slot], (P - mem_order(Yf)[slot]) % P)]
    # rates on a production-size walk (the tag walk of section 3 evaluated in Python, 2,000 steps)
    z = TW.pow(g, rng.randrange(1, R))
    Ms = [TW.mul(TW.pow(g, rng.randrange(R)), TW.pow(z, rng.randrange(R))) for _ in range(RM)]
    Y0 = TW.pow(g, rng.randrange(1, R)); dg = tg.digits(Y0); w = tg.window(dg); Te = TW.one
    wl = {1: 0, 2: 0, 3: 0}; buckets = [0] * RM; e1max = 0
    for st in range(2000):
        if st % LWIN == 0 and st:
            Y0 = TW.mul(Y0, Te); Te = TW.one; dg = tg.digits(Y0); w = tg.window(dg)
        e1max = max(e1max, w["E1"])
        G, V = tg.words(Te)
        s_, l_ = tg.select(dg, G, V, w)
        wl[l_] += 1
        if l_ == 3:
            s_, dist = TY.exact_select(tg.stored_t(TW.mul(Y0, Te)), P, RM, DBITS)
        buckets[s_] += 1
        Te = TW.mul(Te, Ms[s_])
    chi2 = sum((b_ - 2000 / RM) ** 2 / (2000 / RM) for b_ in buckets)
    q1 = q1_bound(tg)
    ok = (rows == n and lem1 == n and lem2 == n and sel == n and rec == n and small == [0, 1, 9]
          and frob[TY.TSLOT] == [6] and 2 in frob[1] and zeta_ok and tg.nd * DIG >= P.bit_length() + 1)
    out.update({"row_constants_up_to_sign": small, "row_identity": f"{rows}/{n}", "recoding": f"{rec}/{n}",
                "level1_lemma": f"{lem1}/{n}", "level2_lemma": f"{lem2}/{n}",
                "certified_selector_exact": f"{sel}/{n}", "levels_of_the_60_pairs": lev,
                "frobenius_powers_fixing_or_negating_coordinate_1": frob[1],
                "frobenius_powers_fixing_or_negating_read_coordinate_6": frob[TY.TSLOT],
                "read_coordinate_scaled_by_primitive_6th_root_under_p^2_frobenius": zeta_ok, "zeta": zeta[0],
                "Smax": tg.Smax, "Smax_log2": round(math.log2(tg.Smax), 4), "THR2": tg.thr2,
                "level1_uncertified_bound": round(float(q1), 6), "level1_uncertified_bound_log2": round(math.log2(q1), 4),
                "level3_bound_log2": round(math.log2((RM * E2 + 2) / 2 ** 64 + 2 ** -90), 2),
                "walk_2000_steps_levels": wl, "walk_level1_uncertified_share": wl[2] / 2000 + wl[3] / 2000,
                "walk_selector_counts": buckets, "walk_selector_chi2_5df": round(chi2, 2),
                "walk_largest_E1_log2": round(math.log2(e1max), 4),
                "table_entries": math.comb(LWIN + RM, RM), "algebra_ok": ok})
    return out


# ---------------------------------------------------------------- ledger

def ledger(k):
    getcontext().prec = 120
    C = 20000
    tg = make_tags()
    F12, FRM = k["f12_mul_executed"], k["fr_mul"] + 4
    STEP, L2 = k["step_executed"], k["level2_executed"]
    EXIT = k["rebase_executed"] + k["tt_entry_exit"] + ALLOW
    E_steps, E_exits, X_cap = step_bounds()
    Q1 = q1_bound(tg)
    E_l2 = -(-(Q1 * E_steps).numerator // (Q1 * E_steps).denominator)     # ceil(q1 E_steps)
    REV = 1 << REV_LOG2
    ent = math.comb(LWIN + RM, RM)
    ST, DM, REC = 1 << ST_LOG2, 1 << DM_LOG2, 16
    ob = OS.bounds(k["fp_mul"], F12)
    other = {"input_check": 508 * F12 + (1 << 16),
             "constants": 1 << 22,
             "multiplier_draws": 16 * (RM * (1017 * F12 + 2 * SCALAR) + MCHECK),
             "table": ent * ENTRY + RM * (1017 * F12 + 2 * SCALAR),
             "walk_starts": ST * (1017 * F12 + 2 * SCALAR + (1 << 10)),
             "reports": DM * REPORT,
             "recoveries": REC * (508 * F12 + 512 * FRM + (1 << 10)),
             "output_fixed": (1 << 16) + (1 << 32),
             "point_output": 1 << 22, "pairing_check": 1 << 28}
    trips = (LWIN + 1) * X_cap
    exp_rows = dict(other, steps=STEP * E_steps, level2=L2 * E_l2, exits=EXIT * E_exits, non_fresh=(STEP + L2 + EXIT) * REV)
    cap_rows = dict(other, steps=STEP * trips, level2=L2 * trips, exits=EXIT * X_cap, non_fresh=0)
    E_run, cap_run = sum(exp_rows.values()), sum(cap_rows.values())
    lim = Fraction(39, 100) * (1 << 64) / LO_RUN
    K = lim.numerator // lim.denominator + (lim.numerator % lim.denominator != 0)   # least K with (K/2^64) LO_RUN >= 0.39
    assert Fraction(K, 1 << 64) * LO_RUN >= Fraction(39, 100) > Fraction(K - 1, 1 << 64) * LO_RUN and K < 1 << 64
    exp_total = GATE + Fraction(K * E_run, 1 << 64)        # expected total work: gate + P(open) x E_run
    every_run = GATE + cap_run
    rep_mean = Fraction(E_steps, 1 << DBITS) + 64          # expected reports, and expected walk starts
    assert rep_mean / ST + rep_mean / DM <= Fraction(1, 1 << 10)
    assert Fraction(E_exits + REV, X_cap) <= Fraction(1, 1 << CAP_LOG2)
    pre = other["constants"] + other["input_check"] + other["table"]
    mem = ent * ENTRY_BYTES + DM * 480 + (1 << 20)
    ln2 = Decimal(2).ln()
    lg = lambda v: (Decimal(v.numerator) / Decimal(v.denominator) if isinstance(v, Fraction) else Decimal(v)).ln() / ln2
    up = lambda v: (lg(v) * 10 ** 6).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** 6)
    up4 = lambda v: "%.4f" % (math.ceil(float(lg(v) - lg(C)) * 10 ** 4) / 10 ** 4) if v else "-"
    score = up(exp_total / C)
    sc = Fraction(int(score * 10 ** 6), 10 ** 6)
    exp_rows["gate"] = cap_rows["gate"] = GATE
    fmt = lambda n: f"{int(n):,}"
    table = [f"| {key} | {fmt(exp_rows[key])} | {up4(exp_rows[key])} | {fmt(cap_rows[key])} | {up4(cap_rows[key])} |" for key in ROWS]
    table.append(f"| **opened run** | **{fmt(GATE + E_run)}** | {up4(GATE + E_run)} | **{fmt(every_run)}** | {up4(every_run)} |")
    return {"step_instructions": STEP, "level2_instructions": L2, "exit_instructions": EXIT,
            "per_step_average": math.ceil(float(Fraction(exp_rows["steps"] + exp_rows["level2"] + exp_rows["exits"], E_steps)) * 100) / 100,
            "E_steps": str(E_steps), "E_steps_log2": float(lg(E_steps)), "E_exits": str(E_exits),
            "E_level2": str(E_l2), "level1_uncertified_bound": float(Q1),
            "X_cap": str(X_cap), "X_cap_log2": float(lg(X_cap)), "table_entries": str(ent),
            "expectation_rows": {k2: str(v) for k2, v in exp_rows.items()},
            "every_run_rows": {k2: str(v) for k2, v in cap_rows.items()},
            "table": table,
            "opened_run_expected_instructions": str(GATE + E_run), "every_run_instructions": str(every_run),
            "every_run_log2": str(up(Fraction(every_run, C))),
            "opened_run_expected_log2": str(up(Fraction(GATE + E_run, C))),
            "gate_K": str(K), "gate_open_probability": float(Fraction(K, 1 << 64)),
            "opened_run_success_lower": "%d/%d" % (LO_RUN.numerator, LO_RUN.denominator),
            "success_lower": str(Decimal(K) / Decimal(1 << 64) * Decimal(LO_RUN.numerator) / Decimal(LO_RUN.denominator))[:22],
            "expected_reports_log2": float(lg(rep_mean)),
            "cap_failure_starts_reports_log2": float(lg(rep_mean / ST + rep_mean / DM)),
            "log2_W_unrounded": str(lg(exp_total / C))[:24],
            "time_log2": str(score), "score_rational": [str(sc.numerator), str(sc.denominator)],
            "preprocessing_log2": str(up(Fraction(pre, C))), "peak_bytes": str(mem), "memory_log2_bytes": str(up(mem)),
            "memory_below_2^63": mem < (1 << 63),
            "output_bounds": {"point": ob["point"], "pairing": ob["pairing"]},
            "output_allowances_ok": ob["point"] <= other["point_output"] and ob["pairing"] <= other["pairing_check"]}


up4f = lambda v: "%.4f" % (math.ceil(v * 10 ** 4) / 10 ** 4)       # stated bounds are rounded up
up2f = lambda v: "%.2f" % (math.ceil(v * 100) / 100)


def toy_stats(toy):
    """Per-row and pooled statistics of toy_offline.json: steps at the stop over the model mean."""
    rows, pool = [], []
    for fam in ("bn", "generic"):
        F_ = toy[fam]
        for key, v in F_["runs"].items():
            m0 = TY.model_steps(F_["r_prime"], v["multipliers"])
            x = [s_ / m0 for s_ in v["steps"]]
            n_ = len(x); mean = sum(x) / n_
            se = math.sqrt(sum((t - mean) ** 2 for t in x) / (n_ - 1) / n_)
            rows.append((fam, key, n_, m0, mean, se, max(x), v))
            pool += x
    n_ = len(pool); mean = sum(pool) / n_
    se = math.sqrt(sum((t - mean) ** 2 for t in pool) / (n_ - 1) / n_)
    return rows, (n_, mean, se, max(pool))


def toy_rows(toy):
    """The proof.md toy table, recomputed from toy_offline.json."""
    rows, (n_, mean, se, mx) = toy_stats(toy)
    out = []
    for fam, key, n1, m0, m1, s1, x1, v in rows:
        tot = sum(v["steps"])
        out.append(f"| {fam} | {v['multipliers']} | {v['window']} | {n1} | {m0:,.0f} | {m1:.4f} | {s1:.4f} | {x1:.2f} | "
                   f"{v['correct']} | {v['wrong']} | {v['at_exit_cap']} | {v['exits'] / tot:.4f} | "
                   f"{v['level1_fail'] / tot:.5f} | {v['level2'] / tot:.5f} | {v['level3'] / tot:.5f} | "
                   f"{v['reports'] / tot:.5f} | {v['check_failures']} of {v['checked']:,} |")
    out.append(f"Pooled over {n_:,} runs: mean {mean:.4f}, standard error {se:.4f}, largest {mx:.2f}.")
    plain = [s_ / math.sqrt(math.pi * toy[fam]["r_prime"] / 2) for fam, key, n1, m0, m1, s1, x1, v in rows for s_ in v["steps"]]
    out.append(f"Negative control, the same runs over sqrt(pi r'/2): mean {sum(plain) / n_:.4f}.")
    return out


def stated_figures(L, k):
    """Compare every figure stated in submission.json, proof.md and machine.md with the computed values."""
    man = json.loads((HERE / "submission.json").read_text())
    st = man["claims"][0]["statement"]
    pr_ = (HERE / "proof.md").read_text()
    mach = (HERE / "machine.md").read_text()
    sc = man["claims"][0]["score"]
    getcontext().prec = 120
    claimed = Decimal(sc["numerator"]) / Decimal(sc["denominator"])
    tg = make_tags()
    fmt = lambda n: f"{int(n):,}"
    need_claim = ["2^" + L["time_log2"], "2^" + L["every_run_log2"], "K = " + L["gate_K"],
                  "%d instructions per walk step" % L["step_instructions"], "%d more" % L["level2_instructions"],
                  fmt(L["exit_instructions"]), fmt(k["f12_mul_executed"]), "2^" + L["memory_log2_bytes"],
                  "2^" + L["preprocessing_log2"], "at least 0.39"]
    need_proof = ["**time_log2 = " + L["time_log2"] + "**", "(exact rational %s/%s)" % tuple(L["score_rational"]),
                  "K = " + L["gate_K"], "E_steps = " + L["E_steps"], "X_cap = " + L["X_cap"],
                  "E_level2 = " + L["E_level2"], "THR2 = %d" % tg.thr2, "Smax = %d" % tg.Smax,
                  "(below 2^%s)" % up4f(L["E_steps_log2"]), "(below 2^%s)" % up4f(L["X_cap_log2"]),
                  "2^" + L["every_run_log2"], "2^" + L["opened_run_expected_log2"],
                  "Preprocessing 2^" + L["preprocessing_log2"], L["opened_run_success_lower"],
                  "**%d** executed instructions" % L["step_instructions"], "**%d** more" % L["level2_instructions"],
                  "**%s** instructions" % fmt(L["exit_instructions"]),
                  "2^" + L["memory_log2_bytes"], fmt(SCALAR), "%.6f" % L["level1_uncertified_bound"],
                  "%.6f" % (math.ceil(L["gate_open_probability"] * 10 ** 6) / 10 ** 6),
                  "2^%s" % up2f(L["expected_reports_log2"]), "2^%s" % up2f(L["cap_failure_starts_reports_log2"]),
                  "%s entries" % fmt(L["table_entries"]), "average %.2f instructions" % L["per_step_average"]]
    need_proof += L["table"]
    op, om = OS.doc_lines(k["fp_mul"], k["f12_mul_executed"], 22, 28)
    need_proof += op
    need_mach = [fmt(k["f12_mul_executed"]), "**%d**" % L["step_instructions"], "**%d**" % L["level2_instructions"],
                 fmt(k["rebase_executed"]), fmt(L["exit_instructions"]), fmt(k["stack_bytes"]) + " bytes", str(k["fp_mul"]),
                 fmt(L["peak_bytes"]) + " = 2^" + L["memory_log2_bytes"], "%s entries" % fmt(L["table_entries"])] + om
    toy = json.loads((HERE / "toy_offline.json").read_text())
    need_proof += toy_rows(toy)
    rows_, (n_, mean, se, mx) = toy_stats(toy)
    plain = [s_ / math.sqrt(math.pi * toy[r_[0]]["r_prime"] / 2) for r_ in rows_ for s_ in r_[7]["steps"]]
    h1_ok = mean + 5 * se <= float(C_MEAN) < sum(plain) / n_ - 5 * se     # H1 support and its negative control
    h2_ok = all(r_[7]["check_failures"] == 0 and r_[7]["wrong"] == 0 and r_[7]["at_exit_cap"] == 0 for r_ in rows_)
    miss = ([x for x in need_claim if x not in st] + [x for x in need_proof if x not in pr_]
            + [x for x in need_mach if x not in mach])
    allowed = {L["time_log2"], L["every_run_log2"], L["opened_run_expected_log2"], L["memory_log2_bytes"], L["preprocessing_log2"]}
    for name, txt in (("claim", st), ("proof.md", pr_), ("machine.md", mach), ("README.md", (HERE / "README.md").read_text())):
        miss += ["%s: 2^%s" % (name, t) for t in re.findall(r"2\^(\d+\.\d{6})", txt) if t not in allowed]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"\d+/1048576", txt) if t != L["opened_run_success_lower"]]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"time_log2 = (\d+\.\d+)", txt) if t != L["time_log2"]]
    return {"claim_score_equals_ledger": claimed == Decimal(L["time_log2"]), "missing_figures": miss,
            "toy_pooled_mean": round(mean, 4), "toy_standard_error": round(se, 4),
            "toy_mean_plus_5se_within_allowance_and_control_outside": h1_ok, "toy_selector_checks_clean": h2_ok,
            "stated_figures_match": claimed == Decimal(L["time_log2"]) and not miss and h1_ok and h2_ok}


def main():
    args = json.load(sys.stdin)
    mode, seed = args.get("mode", "all"), args.get("seed", 1)
    out = {}
    k = kernels(seed)
    L = ledger(k)
    out["stated_figures"] = stated_figures(L, k)          # every mode: declared figures against the computed ledger
    rx = random.Random(seed).randrange(1, R)
    out["output_step"] = OS.check([1, R - 1, rx], k["fp_mul"])   # every mode: P = [x]P0 executed and checked
    alg = None
    if mode in ("all", "kernels"):
        out["kernels"] = k
    if mode in ("all", "algebra"):
        alg = out["algebra"] = algebra(seed)
    if mode in ("all", "ledger"):
        out["parameters"] = PARAM
        out["ledger"] = L
    toy_bad = False
    if mode in ("all", "toy"):
        t_ = out["toy"] = TY.toy(seed, args.get("toy_instances", 32), RM, args.get("toy_window", 4))
        cen = out["toy_census"] = TY.census(-7)
        q_, r_ = cen["q"], cen["r_prime"]
        c1, c6 = cen["slots"]["1"]["distinct_values"], cen["slots"][str(TY.TSLOT)]["distinct_values"]
        rnd_map = q_ * (1 - math.exp(-r_ / q_))                # distinct values of a random map of r' points to q values
        cen["coordinate_1_at_most_1_plus_(r-1)/6"] = c1 <= 1 + (r_ - 1) // 6
        cen["read_coordinate_distinct_over_random_map"] = round(c6 / rnd_map, 4)
        toy_bad = bool(t_["wrong_outputs"] or t_["selector_check_failures"] or t_["reference_walk_mismatches"]
                       or not cen["coordinate_1_at_most_1_plus_(r-1)/6"] or abs(c6 / rnd_map - 1) > 0.02)
    out["scope"] = ("kernel counts are executed instructions of the generated programs; algebra on the pinned "
                    "constants; toy runs are evidence at toy scale only; no production attack is run")
    print(json.dumps(out, indent=1))
    bad = [n for n in ("fp_mismatches", "f12_mismatches", "walk_mismatches") if k[n]]
    bad += ["output_step"] * (not out["output_step"]["output_step_ok"]) + ["output_allowances"] * (not L["output_allowances_ok"])
    bad += ["algebra"] * (alg is not None and not alg["algebra_ok"]) + ["toy"] * toy_bad
    if bad or not out["stated_figures"]["stated_figures_match"]:
        print("FAIL: " + ", ".join(bad + ["stated_figures"] * (not out["stated_figures"]["stated_figures_match"])), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

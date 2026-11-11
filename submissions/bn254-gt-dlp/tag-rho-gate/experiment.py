"""Supporting checks for the BN254 GT discrete-logarithm tag-walk claim. Python standard library only.
Modes (input JSON on stdin): kernels (generate and execute the priced AArch64 programs against big
integers, print executed counts), algebra (facts of proof.md on the pinned constants, the tag lemma on
production-size values), ledger (every integer of the cost ledger, the gate constant, the score), toy
(the complete algorithm in a BN-shaped toy field), all. Every mode compares
every stated figure with the computed values. Toy runs are evidence at toy scale only.
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
from kern import (M64, P, R, U, G_PINNED, CONST, Program, Machine, Tower, limbs, gen_fp_mul, gen_fp_add,   # noqa: E402
                  gen_fp_sub, gen_f12_mul, mem_order, from_mem, is_prime)

RM, LWIN, DBITS = 6, 1024, 90        # multipliers, window length, distinguished bits
DIG, NDIG = TY.DIG, 72               # 43-bit digits, 6 per coordinate
ERR = NDIG << DIG                    # tag error bound E of the lemma
THR = TY.threshold(RM, ERR, DBITS)   # a tag is certified when RM * A mod 2^64 is below THR
WPOOL = 512                          # offset of the walk constants in the pool
EOFF, CH, AL, BE, TT = 288, 288, 352, 384, 416   # entry: words at e-288..e+287, children, alpha, beta, T_e
ENTRY_BYTES = EOFF + TT + 384
MONT = 1 << 256


# ---------------------------------------------------------------- priced walk programs

HOLD = ['x4', 'x5', 'x6', 'x7', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17', 'x21', 'x22', 'x26', 'x27', 'x29', 'x30']
HOLD_SAVE = (('x21', 'x22'), ('x26', 'x27'), ('x29', 'x30'))   # callee-saved registers among them
TT_IO = 3 + len(HOLD_SAVE) + len(HOLD) // 2 + len(HOLD_SAVE) + 2   # entry and exit instructions of one call
TRIP = 2 * NDIG - len(HOLD) // 2 + 9                       # instructions of one full trip


def gen_tt(pr):
    """The step loop. x19 = digits of Y0 (base + 288), x20 = entry e. One trip is one walk step:
    tag A = sum of digit * word mod 2^64, certificate, selector, child pointer, trip count.
    Returns x0 = 1 (not certified; x8 = A), 2 (window end; x10 = selector, x11 = e + 8 selector) or
    3 (trip budget used)."""
    e = pr.emit
    pr.label('tt_run')
    e('add', 'x9', 'x28', WPOOL); e('ldp', 'x24', 'x25', 'x9', 0); e('ldr', 'x23', 'x9', 24)
    for a_, b_ in HOLD_SAVE:
        e('stp_pre', a_, b_, 'sp', -16)
    for j in range(len(HOLD) // 2):                        # the first 16 digits stay in registers for the call
        e('ldp', HOLD[2 * j], HOLD[2 * j + 1], 'x19', 16 * j - EOFF)
    pr.label('tt')
    for j in range(NDIG // 2):
        if 2 * j < len(HOLD):
            d0, d1 = HOLD[2 * j], HOLD[2 * j + 1]
        else:
            d0, d1 = 'x0', 'x1'
            e('ldp', 'x0', 'x1', 'x19', 16 * j - EOFF)
        e('ldp', 'x2', 'x3', 'x20', 16 * j - EOFF)
        if j == 0:
            e('mul', 'x8', d0, 'x2')
        else:
            e('madd', 'x8', d0, 'x2', 'x8')
        e('madd', 'x8', d1, 'x3', 'x8')
    e('mul', 'x9', 'x8', 'x24'); e('cmp', 'x9', 'x25'); e('bc', 'cs', 'tt_unc')
    e('umulh', 'x10', 'x8', 'x24'); e('add_lsl', 'x11', 'x20', 'x10', 3); e('ldr', 'x20', 'x11', CH)
    e('cbz', 'x20', 'tt_wend')
    e('subi', 'x23', 'x23', 1); e('bc', 'ne', 'tt')
    for lab, code in ((None, 3), ('tt_unc', 1), ('tt_wend', 2)):
        if lab:
            pr.label(lab)
        for a_, b_ in HOLD_SAVE[::-1]:
            e('ldp_post', a_, b_, 'sp', 16)
        e('movi', 'x0', code); e('ret')


def gen_rebase(pr):
    """Re-base after an exit. x12 = entry e, x21 = Y0 (12 stored coordinates), x22 = spare buffer,
    x19 = digits (base + 288), x23 = trips left of the call, x27 = state (alpha, beta, exit counter 3 words,
    walk trip counter 2 words, two flags). Y0 <- Y0 * T_e, (alpha, beta) += (alpha_e, beta_e) mod r, digits of the new Y0, counters."""
    e = pr.emit
    pr.label('rebase')
    e('stp_pre', 'x29', 'x30', 'sp', -16); e('stp_pre', 'x26', 'x27', 'sp', -16)
    e('mov', 'x26', 'x12')
    e('mov', 'x0', 'x22'); e('mov', 'x1', 'x21'); e('add', 'x2', 'x26', TT); e('bl', 'f12_mul')
    e('mov', 'x0', 'x21'); e('mov', 'x21', 'x22'); e('mov', 'x22', 'x0')
    e('mov', 'x0', 'x27'); e('mov', 'x1', 'x27'); e('add', 'x2', 'x26', AL); e('bl', 'fr_add')
    e('add', 'x0', 'x27', 32); e('add', 'x1', 'x27', 32); e('add', 'x2', 'x26', BE); e('bl', 'fr_add')
    e('add', 'x13', 'x28', WPOOL); e('ldr', 'x14', 'x13', 32)
    for k in range(12):
        e('ldp', 'x0', 'x1', 'x21', 32 * k); e('ldp', 'x2', 'x3', 'x21', 32 * k + 16)
        e('and_r', 'x4', 'x0', 'x14')
        e('extr', 'x5', 'x1', 'x0', 43); e('and_r', 'x5', 'x5', 'x14')
        e('extr', 'x6', 'x2', 'x1', 22); e('and_r', 'x6', 'x6', 'x14')
        e('lsr', 'x7', 'x2', 1); e('and_r', 'x7', 'x7', 'x14')
        e('extr', 'x8', 'x3', 'x2', 44); e('and_r', 'x8', 'x8', 'x14')
        e('lsr', 'x9', 'x3', 23)
        o = 48 * k - EOFF
        e('stp', 'x4', 'x5', 'x19', o); e('stp', 'x6', 'x7', 'x19', o + 16); e('stp', 'x8', 'x9', 'x19', o + 32)
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
CAP_LOG2 = 8               # exit cap X_cap = 2^8 times the expected-exit bound (Markov: reached with probability <= 2^-8)
UNC_LOG2 = 12              # H2: a step is not certified with probability at most 2^-12
GATE = 64                  # the gate: one coin word (16), compare with K, branch, return
ALLOW = 1 << 10            # per exit: dispatch, exact selector and distinguished test, root child, flags
ENTRY = 1 << 20            # per table entry: product, row, words, labels, pointers (machine.md)
REPORT = 1 << 18           # per report: store lookup and insert on 48-word keys
ST_LOG2 = DM_LOG2 = 50     # caps on walk starts and reports
WALK_LOG2 = 97             # trips (steps) of one walk
REVISIT = 1 << 98          # H3: expected evaluations of elements the run has already visited (2^8 mean walk lengths)
LO_RUN = 1 - Fraction(1, 1 << 8) - Fraction(1, 1 << 10) - Fraction(1, 1 << 20)   # success of an opened run
PARAM = {"multipliers": RM, "window": LWIN, "dp_bits": DBITS, "digit_bits": DIG, "digits": NDIG,
         "tag_error_bound": "72 x 2^43", "mean_allowance": "107/100",
         "exit_cap": "2^8 x (ceil(E_steps / 1024) + floor(E_steps / 2^12) + 1 + 2^98)",
         "gate": "one fresh 64-bit coin before any work; the run opens when it is below K"}
ROWS = ("gate", "steps", "exits", "input_check", "constants", "table", "walk_starts", "reports", "recoveries",
        "output_fixed")


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
    E_exits = -(-E_steps // LWIN) + (E_steps >> UNC_LOG2) + 1 + REVISIT   # every revisit is charged an exit
    return E_steps, E_exits, E_exits << CAP_LOG2


# ---------------------------------------------------------------- kernels

def write_entry(mc, tg, addr, Te, al, be, children):
    """Table entry at addr (the pointer e): words, child pointers, labels, T_e in stored form."""
    for i, w in enumerate(tg.words(Te)):
        mc.st(addr - EOFF + 8 * i, w)
    for s_ in range(RM):
        mc.st(addr + CH + 8 * s_, children.get(s_, 0))
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
    tg = TY.Tags(TW, MONT)
    assert tg.ndig == NDIG and tg.E == ERR
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
    pool = [RM, THR, 0, LWIN + 1, (1 << DIG) - 1, X_cap & M64, (X_cap >> 64) & M64, X_cap >> 128, W_cap & M64, W_cap >> 64]
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
    # the walk programs on a real walk: target z, multipliers, a chain of table entries, one window
    G0 = 1 << 28
    BUF1, BUF2, DGB, STATE, TAB = G0, G0 + 4096, G0 + 8192 + EOFF, G0 + 12288, 1 << 44
    z = TW.pow(g, rng.randrange(1, R))
    cd = [(rng.randrange(R), rng.randrange(R)) for _ in range(RM)]
    Ms = [TW.mul(TW.pow(g, c), TW.pow(z, d)) for c, d in cd]
    wbad, trips, lemma_max, trip_counts, reb_counts = 0, 0, 0, set(), set()

    def put_digits(Y):
        for i, d in enumerate(tg.digits(Y)): mc.st(DGB - EOFF + 8 * i, d)

    def lemma(Y, Te, A):
        """(2^64 t / p - A) mod 2^64 lies in [0, E), with t the stored coordinate of Y * Te."""
        t = tg.stored_t(TW.mul(Y, Te))
        d = (Fraction(t << 64, P) - A) % (1 << 64)
        return t, d
    for walk in range(3):
        n = 3 + walk                                       # window of n steps, then the window end
        while True:                                        # a start whose first n + 1 tags are all certified
            a0, b0 = rng.randrange(R), rng.randrange(R)
            Y0 = TW.mul(TW.pow(g, a0), TW.pow(z, b0))
            dg = tg.digits(Y0)
            Te, al, be, sels, addr = TW.one, 0, 0, [], TAB + (walk << 20)
            chain, allcert = [], True
            for i in range(n + 1):                         # entries along the walk's own path; others null
                A = tg.tag(dg, tg.words(Te))
                t, d = lemma(Y0, Te, A)
                lemma_max = max(lemma_max, d)
                s_ex, dist = TY.exact_select(t, P, RM, DBITS)
                cert = (RM * A) & M64 < THR
                allcert &= cert
                wbad += not (0 <= d < ERR) or (cert and (dist or (RM * A) >> 64 != s_ex))
                chain.append((addr + 2048 * i, Te, al, be, s_ex))
                sels.append(s_ex)
                Te, al, be = TW.mul(Te, Ms[s_ex]), (al + cd[s_ex][0]) % R, (be + cd[s_ex][1]) % R
            if allcert:
                break
        for k, v in enumerate(mem_order(Y0)): setc(BUF1 + 32 * k, mont(v))
        put_digits(Y0)
        for i, (ad, T_, al_, be_, s_) in enumerate(chain):
            write_entry(mc, tg, ad, T_, al_, be_, {s_: chain[i + 1][0]} if i < n else {})
        mc.r['x19'], mc.r['x20'] = DGB, chain[0][0]
        keep = {reg: rng.getrandbits(64) for pair in HOLD_SAVE for reg in pair if reg != 'x30'}
        mc.r.update(keep); sp_in = mc.r['sp']
        c_run = mc.call('tt_run')
        wbad += any(mc.r[reg] != v for reg, v in keep.items()) or mc.r['sp'] != sp_in   # saved registers restored
        trips += n + 1
        trip_counts.add((c_run - TT_IO - (TRIP - 2)) // n if (c_run - TT_IO - (TRIP - 2)) % n == 0 else -1)
        e_last, T_last, al_last, be_last, s_last = chain[n]
        wbad += (mc.r['x0'] != 2 or mc.r['x10'] != s_last or mc.r['x11'] != e_last + 8 * s_last
                 or mc.r['x23'] != LWIN + 1 - n)
        # re-base at the window end
        cx, cw = [(1 << 128) - 1, X_cap - 1, (M64 << 64) + rng.getrandbits(60) % X_cap][walk], [M64, rng.getrandbits(60), W_cap - 3][walk]
        setc(STATE, a0); setc(STATE + 32, b0)
        mc.st(STATE + 64, cx & M64); mc.st(STATE + 72, (cx >> 64) & M64); mc.st(STATE + 80, cx >> 128)
        mc.st(STATE + 88, cw & M64); mc.st(STATE + 96, cw >> 64)
        mc.r['x12'], mc.r['x21'], mc.r['x22'], mc.r['x27'] = mc.r['x11'] - 8 * mc.r['x10'], BUF1, BUF2, STATE
        sp0 = mc.r['sp']
        c_reb = mc.call('rebase')
        reb_counts.add(c_reb)
        Y1 = TW.mul(Y0, T_last)
        got = from_mem([unmont(getv(mc.r['x21'] + 32 * k)) for k in range(12)])
        cx1 = mc.ld(STATE + 64) + (mc.ld(STATE + 72) << 64) + (mc.ld(STATE + 80) << 128)
        cw1 = mc.ld(STATE + 88) + (mc.ld(STATE + 96) << 64)
        wbad += (got != Y1 or mc.r['x21'] != BUF2 or mc.r['x22'] != BUF1 or mc.r['sp'] != sp0 or mc.r['x27'] != STATE
                 or getv(STATE) != (a0 + al_last) % R or getv(STATE + 32) != (b0 + be_last) % R
                 or [mc.ld(DGB - EOFF + 8 * i) for i in range(NDIG)] != tg.digits(Y1)
                 or cx1 != cx + 1 or cw1 != cw + n + 1 or mc.ld(STATE + 104) != int(cx + 1 >= X_cap)
                 or mc.ld(STATE + 112) != int(cw + n + 1 >= W_cap))
        # the represented element equals the independently multiplied walk element
        ref = Y0
        for s_ in sels[:-1]:
            ref = TW.mul(ref, Ms[s_])
        wbad += ref != Y1 or TW.mul(TW.pow(g, (a0 + al_last) % R), TW.pow(z, (b0 + be_last) % R)) != Y1
    # crafted entries: both sides of the certificate threshold, and the trip budget
    Yc = TW.pow(g, rng.randrange(1, R))
    while tg.digits(Yc)[0] % 2 == 0:
        Yc = TW.mul(Yc, g)
    put_digits(Yc)
    d0 = tg.digits(Yc)[0]
    half = pow(RM // 2, -1, 1 << 63)
    for lo, code in ((THR - 1, 2), (THR + 1, 1), ((1 << 64) - RM, 1), (0, 2)):
        A = (lo // 2) * half % (1 << 63)                   # RM * A mod 2^64 = lo (lo even)
        assert (RM * A) & M64 == lo
        ad = TAB + (9 << 20)
        for i in range(NDIG): mc.st(ad - EOFF + 8 * i, 0)
        mc.st(ad - EOFF, A * pow(d0, -1, 1 << 64) & M64)
        for s_ in range(RM): mc.st(ad + CH + 8 * s_, 0)
        mc.r['x19'], mc.r['x20'] = DGB, ad
        c_ = mc.call('tt_run')
        wbad += mc.r['x0'] != code or (code == 1 and (mc.r['x8'] != A or c_ != TT_IO + TRIP - 6))
    ad = TAB + (10 << 20)                                  # a ring of two entries: the budget ends the call
    for j in (0, 1):                                       # all words zero: tag 0, certified, selector 0
        for i in range(NDIG): mc.st(ad + 2048 * j - EOFF + 8 * i, 0)
        for s_ in range(RM): mc.st(ad + 2048 * j + CH + 8 * s_, ad + 2048 * (1 - j))
    mc.st(CONST + WPOOL + 24, 5)
    mc.r['x19'], mc.r['x20'] = DGB, ad
    c_bud = mc.call('tt_run')
    wbad += mc.r['x0'] != 3 or mc.r['x23'] != 0
    mc.st(CONST + WPOOL + 24, LWIN + 1)
    step = trip_counts.pop() if len(trip_counts) == 1 else -1
    wbad += step != TRIP or c_bud != TT_IO + 5 * step or len(reb_counts) != 1
    regs = set()
    for op, a in pr.code:
        regs |= {x for x in a if isinstance(x, str) and x.startswith('x') and x[1:].isdigit()}
    return {"fp_mul": cnt['fp_mul'], "fp_add": cnt['fp_add'], "fp_sub": cnt['fp_sub'],
            "fr_mul": cnt['fr_mul'], "fr_add": cnt['fr_add'], "fp_mismatches": bad,
            "f12_mul_executed": n12, "f12_calls": t12.calls, "f12_scratch_slots": t12.peak, "f12_mismatches": bad12,
            "step_executed": step, "tt_entry_exit": TT_IO, "rebase_executed": c_reb,
            "walk_steps_executed": trips, "walk_mismatches": wbad,
            "largest_tag_offset_log2": round(math.log2(max(1, lemma_max)), 2),
            "walk_tests": "3 real windows with re-base (selector, lemma, element, labels, digits, counters and both caps); "
                          "both sides of the certificate threshold; the trip budget",
            "stack_bytes": mc.depth, "registers_used": sorted(regs, key=lambda s: int(s[1:])),
            "program_instructions": len(pr.code)}


# ---------------------------------------------------------------- algebra

def algebra(seed):
    rng = random.Random(seed)
    TW = Tower(P, 9)
    tg = TY.Tags(TW, MONT)
    g = TW.to_tower(G_PINNED)
    out = {"p_r_prime": is_prime(P) and is_prime(R),
           "bn_polynomials": P == 36 * U ** 4 + 36 * U ** 3 + 24 * U ** 2 + 6 * U + 1 and R == 36 * U ** 4 + 36 * U ** 3 + 18 * U ** 2 + 6 * U + 1,
           "g_order_r": TW.pow(g, R) == TW.one and g != TW.one}
    rows = lem = sel = unc = 0
    worst = 0
    n = 60
    for i in range(n):
        if i < 6:                                           # extreme digits: every coordinate p - 1 or 2^254 - like
            Y = from_mem([(P - 1 - i) * pow(MONT, -1, P) % P] * 12)
        elif i < 30:
            Y = from_mem([rng.randrange(P) for _ in range(12)])
        else:
            Y = TW.pow(g, rng.randrange(1, R))
        Te = from_mem([rng.randrange(P) for _ in range(12)]) if i < 30 else TW.pow(g, rng.randrange(1, R))
        v = tg.row(Te)
        prod = TW.mul(Y, Te)
        rows += mem_order(prod)[TY.TSLOT] == sum(y * w for y, w in zip(mem_order(Y), v)) % P
        A = tg.tag(tg.digits(Y), tg.words(Te))
        t = tg.stored_t(prod)
        d = (Fraction(t << 64, P) - A) % (1 << 64)
        lem += 0 <= d < ERR
        worst = max(worst, d)
        s_ex, dist = TY.exact_select(t, P, RM, DBITS)
        if (RM * A) & M64 < THR:
            sel += (RM * A) >> 64 == s_ex and not dist
        else:
            unc += 1; sel += 1
    conj = 0                                               # no Frobenius power fixes the walk's coordinate on GT
    for _ in range(4):
        y = TW.pow(g, rng.randrange(1, R))
        conj += len({tg.stored_t(TW.frob_k(y, k_)) for k_ in range(12)}) == 12
    ent = math.comb(LWIN + RM, RM)
    small = all(c in (0, 1, P - 1, 9, P - 9) for row in tg.mat for c in row)   # structure constants of the row
    out.update({"coordinate_distinct_on_the_12_conjugates": f"{conj}/4", "row_constants_in_0_1_9": small, "row_identity": f"{rows}/{n}", "tag_lemma": f"{lem}/{n}", "certified_selector_exact": f"{sel}/{n}",
                "uncertified_cases": unc, "largest_tag_offset_log2": round(math.log2(max(1, worst)), 2),
                "tag_error_bound_log2": round(math.log2(ERR), 4),
                "threshold": THR, "uncertified_fraction_log2": round(math.log2(((1 << 64) - THR) / 2 ** 64), 4),
                "table_entries": ent, "digits_cover_p": DIG * (NDIG // 12) >= P.bit_length(),
                "algebra_ok": rows == n and lem == n and sel == n and small and conj == 4 and DIG * (NDIG // 12) >= P.bit_length()
                and ((1 << 64) - THR) <= 1 << (64 - UNC_LOG2)})
    return out


# ---------------------------------------------------------------- ledger

def ledger(k):
    getcontext().prec = 120
    C = 20000
    F12, FRM = k["f12_mul_executed"], k["fr_mul"] + 4
    STEP = k["step_executed"]
    EXIT = k["rebase_executed"] + k["tt_entry_exit"] + ALLOW
    E_steps, E_exits, X_cap = step_bounds()
    ent = math.comb(LWIN + RM, RM)
    ST, DM, REC = 1 << ST_LOG2, 1 << DM_LOG2, 16
    other = {"input_check": 508 * F12 + (1 << 16),
             "constants": 1 << 22,
             "table": ent * ENTRY + RM * (1017 * F12 + 2 * SCALAR),
             "walk_starts": ST * (1017 * F12 + 2 * SCALAR + (1 << 10)),
             "reports": DM * REPORT,
             "recoveries": REC * (508 * F12 + 512 * FRM + (1 << 10)),
             "output_fixed": (1 << 16) + (1 << 32)}
    exp_rows = dict(other, steps=STEP * E_steps, exits=EXIT * E_exits)
    cap_rows = dict(other, steps=STEP * (LWIN + 1) * X_cap, exits=EXIT * X_cap)
    E_run, cap_run = sum(exp_rows.values()), sum(cap_rows.values())
    lim = Fraction(39, 100) * (1 << 64) / LO_RUN
    K = lim.numerator // lim.denominator + (lim.numerator % lim.denominator != 0)   # least K with (K/2^64) LO_RUN >= 0.39
    assert Fraction(K, 1 << 64) * LO_RUN >= Fraction(39, 100) > Fraction(K - 1, 1 << 64) * LO_RUN and K < 1 << 64
    exp_total = GATE + Fraction(K * E_run, 1 << 64)        # expected total work: gate + P(open) x E_run
    every_run = GATE + cap_run
    rep_mean = Fraction(2 * E_steps, 1 << DBITS) + (1 << 18)   # expected reports, and expected walk starts (F2)
    assert rep_mean / ST + rep_mean / DM <= Fraction(1, 1 << 10)
    pre = other["constants"] + other["input_check"] + other["table"]
    mem = ent * ENTRY_BYTES + DM * 480 + (1 << 20)
    ln2 = Decimal(2).ln()
    lg = lambda v: (Decimal(v.numerator) / Decimal(v.denominator) if isinstance(v, Fraction) else Decimal(v)).ln() / ln2
    up = lambda v: (lg(v) * 10 ** 6).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** 6)
    up4 = lambda v: "%.4f" % (math.ceil(float(lg(v) - lg(C)) * 10 ** 4) / 10 ** 4)
    score = up(exp_total / C)
    sc = Fraction(int(score * 10 ** 6), 10 ** 6)
    exp_rows["gate"] = cap_rows["gate"] = GATE
    fmt = lambda n: f"{int(n):,}"
    table = [f"| {key} | {fmt(exp_rows[key])} | {up4(exp_rows[key])} | {fmt(cap_rows[key])} | {up4(cap_rows[key])} |" for key in ROWS]
    table.append(f"| **opened run** | **{fmt(GATE + E_run)}** | {up4(GATE + E_run)} | **{fmt(every_run)}** | {up4(every_run)} |")
    return {"step_instructions": STEP, "exit_instructions": EXIT,
            "per_step_average": math.ceil(float(Fraction(exp_rows["steps"] + exp_rows["exits"], E_steps)) * 100) / 100,
            "E_steps": str(E_steps), "E_steps_log2": float(lg(E_steps)), "E_exits": str(E_exits),
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
            "memory_below_2^63": mem < (1 << 63)}


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
        wl = 1 << toy[fam]["dp_bits"]                        # the toy's mean walk length
        out.append(f"| {fam} | {v['multipliers']} | {v['window']} | {n1} | {m0:,.0f} | {m1:.4f} | {s1:.4f} | {x1:.2f} | "
                   f"{v['correct']} | {v['wrong']} | {v['at_exit_cap']} | {v['exits'] / tot:.4f} | "
                   f"{v['uncertified'] / tot:.5f} | {v['reports'] / tot:.5f} | {v['check_failures']} of {v['checked']:,} | "
                   f"{sum(v['revisits']) / n1 / wl:.2f} | {max(v['revisits']) / wl:.1f} |")
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
    fmt = lambda n: f"{int(n):,}"
    need_claim = ["2^" + L["time_log2"], "2^" + L["every_run_log2"], "K = " + L["gate_K"],
                  "%d instructions per walk step" % L["step_instructions"], fmt(L["exit_instructions"]),
                  fmt(k["f12_mul_executed"]), "2^" + L["memory_log2_bytes"], "2^" + L["preprocessing_log2"], "at least 0.39"]
    need_proof = ["**time_log2 = " + L["time_log2"] + "**", "(exact rational %s/%s)" % tuple(L["score_rational"]),
                  "K = " + L["gate_K"], "E_steps = " + L["E_steps"], "X_cap = " + L["X_cap"],
                  "(below 2^%s)" % up4f(L["E_steps_log2"]), "(below 2^%s)" % up4f(L["X_cap_log2"]),
                  "2^" + L["every_run_log2"], "2^" + L["opened_run_expected_log2"],
                  "Preprocessing 2^" + L["preprocessing_log2"], L["opened_run_success_lower"],
                  "**%d** executed instructions" % L["step_instructions"], "**%s** instructions" % fmt(L["exit_instructions"]),
                  "2^" + L["memory_log2_bytes"], fmt(SCALAR),
                  "%.6f" % (math.ceil(L["gate_open_probability"] * 10 ** 6) / 10 ** 6),
                  "2^%s" % up2f(L["expected_reports_log2"]), "2^%s" % up2f(L["cap_failure_starts_reports_log2"]),
                  "THR = %d" % THR, "%s entries" % fmt(L["table_entries"]), "average %.2f instructions" % L["per_step_average"]]
    need_proof += L["table"]
    need_mach = [fmt(k["f12_mul_executed"]), "**%d**" % L["step_instructions"], fmt(k["rebase_executed"]),
                 fmt(L["exit_instructions"]), fmt(k["stack_bytes"]) + " bytes", str(k["fp_mul"]),
                 fmt(L["peak_bytes"]) + " = 2^" + L["memory_log2_bytes"], "%s entries" % fmt(L["table_entries"])]
    toy = json.loads((HERE / "toy_offline.json").read_text())
    need_proof += toy_rows(toy)
    rows_, (n_, mean, se, mx) = toy_stats(toy)
    plain = [s_ / math.sqrt(math.pi * toy[r_[0]]["r_prime"] / 2) for r_ in rows_ for s_ in r_[7]["steps"]]
    h1_ok = mean + 5 * se <= float(C_MEAN) < sum(plain) / n_ - 5 * se     # H1 support and its negative control
    h2_ok = all(r_[7]["check_failures"] == 0 and r_[7]["wrong"] == 0 for r_ in rows_)
    rv = [x / (1 << toy[r_[0]]["dp_bits"]) for r_ in rows_ for x in r_[7]["revisits"]]   # H3: revisits in walk lengths
    rvm = sum(rv) / len(rv)
    h3_ok = rvm + 5 * math.sqrt(sum((x - rvm) ** 2 for x in rv) / (len(rv) - 1) / len(rv)) <= 16   # 2^4, against 2^8 allowed
    miss = ([x for x in need_claim if x not in st] + [x for x in need_proof if x not in pr_]
            + [x for x in need_mach if x not in mach])
    allowed = {L["time_log2"], L["every_run_log2"], L["opened_run_expected_log2"], L["memory_log2_bytes"], L["preprocessing_log2"]}
    for name, txt in (("claim", st), ("proof.md", pr_), ("machine.md", mach), ("README.md", (HERE / "README.md").read_text())):
        miss += ["%s: 2^%s" % (name, t) for t in re.findall(r"2\^(\d+\.\d{6})", txt) if t not in allowed]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"\d+/1048576", txt) if t != L["opened_run_success_lower"]]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"time_log2 = (\d+\.\d+)", txt) if t != L["time_log2"]]
    return {"claim_score_equals_ledger": claimed == Decimal(L["time_log2"]), "missing_figures": miss,
            "toy_mean_plus_5se_within_allowance_and_control_outside": h1_ok, "toy_selector_checks_clean": h2_ok,
            "toy_revisits_within_allowance": h3_ok,
            "stated_figures_match": claimed == Decimal(L["time_log2"]) and not miss and h1_ok and h2_ok and h3_ok}


def main():
    args = json.load(sys.stdin)
    mode, seed = args.get("mode", "all"), args.get("seed", 1)
    out = {}
    k = kernels(seed)
    L = ledger(k)
    out["stated_figures"] = stated_figures(L, k)          # every mode: declared figures against the computed ledger
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
        toy_bad = bool(t_["wrong_outputs"] or t_["selector_check_failures"])
    out["scope"] = ("kernel counts are executed instructions of the generated programs; algebra on the pinned "
                    "constants; toy runs are evidence at toy scale only; no production attack is run")
    print(json.dumps(out, indent=1))
    bad = [n for n in ("fp_mismatches", "f12_mismatches", "walk_mismatches") if k[n]]
    bad += ["algebra"] * (alg is not None and not alg["algebra_ok"]) + ["toy"] * toy_bad
    if bad or not out["stated_figures"]["stated_figures_match"]:
        print("FAIL: " + ", ".join(bad + ["stated_figures"] * (not out["stated_figures"]["stated_figures_match"])), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

"""Supporting checks for the BN254 pairing-inversion tag-walk claim. Python standard library only.
Modes (input JSON on stdin): kernels (generate and execute the priced AArch64 programs against big
integers, print executed counts), algebra (facts of proof.md on the pinned constants, the tag lemma on
production-size values), ledger (every integer of the cost ledger, the gate constant, the score), toy
(the complete algorithm in a BN-shaped toy field), all. Every mode executes the output step and compares
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
import outstep as OS                                 # noqa: E402
import toy as TY                                     # noqa: E402
from kern import (M64, P, R, U, G_PINNED, CONST, Program, Machine, Tower, limbs, gen_fp_mul, gen_fp_add,   # noqa: E402
                  gen_fp_sub, gen_f12_mul, mem_order, from_mem, is_prime)

RM, LWIN, DBITS = 5, 3500, 90        # multipliers, window length, distinguished bits
DIG, NDIG = TY.DIG, 60               # 51-bit digits, 5 per coordinate
ERR = NDIG << DIG                    # first-tag error bound E1 of lemma 1
ERR2 = 2 * NDIG                      # second-tag error bound E2 of lemma 2
THR = TY.threshold(RM, ERR, DBITS)   # the first tag is certified when RM * A mod 2^64 is below THR
THR2 = TY.threshold(RM, ERR2, DBITS)  # the second tag is certified when RM * A2 mod 2^64 is below THR2
WPOOL = 512                          # offset of the walk constants in the pool
EOFF, CH, AL, BE, TT = 480, 480, 544, 576, 608   # entry: first words at e-480..e-8, second words at e..e+472, children, alpha, beta, T_e
DOFF = 240                           # digits at x19-240..x19+232, THR2 at x19+240
ENTRY_BYTES = EOFF + TT + 384
MONT = 1 << 256


# ---------------------------------------------------------------- priced walk programs

RES_PAIRS = tuple(range(20, 30))     # digit pairs held in registers for a whole call of the step loop
RES_REGS = (('x4', 'x5'), ('x6', 'x7'), ('x9', 'x10'), ('x11', 'x12'), ('x13', 'x14'), ('x15', 'x16'),
            ('x17', 'x26'), ('x21', 'x22'), ('x27', 'x28'), ('x29', 'x30'))
SAVED = (('x21', 'x22'), ('x27', 'x28'), ('x29', 'x30'))     # callee-saved registers the loop borrows


def gen_tt(pr):
    """The step loop. x19 = digits of Y0 (base + 240; THR2 at x19 + 240), x20 = entry e. The call saves x21, x22
    and x27 to x30, loads the constants and the 20 digits of pairs 20..29 into registers. One trip is one walk
    step: first tag A = sum of digit * first word mod 2^64 (pairs 0..19 load their digits, pairs 20..29 use the
    registers) and its certificate; if it is not certified, the second tag A2 = A + sum of UMULH(digit, second
    word) and its certificate; then selector, child pointer, trip count. Returns x0 = 1 (second tag not
    certified; x8 = A2), 2 (window end; x10 = selector, x11 = e + 8 selector) or 3 (trip budget used), with the
    borrowed registers restored."""
    e = pr.emit
    res = dict(zip(RES_PAIRS, RES_REGS))
    n0 = len(pr.code)
    pr.label('tt_run')
    for a, b in SAVED:
        e('stp_pre', a, b, 'sp', -16)
    e('add', 'x26', 'x28', WPOOL); e('ldp', 'x24', 'x25', 'x26', 0); e('ldr', 'x23', 'x26', 24)
    for j in RES_PAIRS:
        e('ldp', res[j][0], res[j][1], 'x19', 16 * j - DOFF)
    n_entry = len(pr.code) - n0
    pr.label('tt')
    for j in range(NDIG // 2):
        if j in res:
            ra, rb = res[j]
        else:
            ra, rb = 'x0', 'x1'
            e('ldp', 'x0', 'x1', 'x19', 16 * j - DOFF)
        e('ldp', 'x2', 'x3', 'x20', 16 * j - EOFF)
        if j == 0:
            e('mul', 'x8', ra, 'x2')
        else:
            e('madd', 'x8', ra, 'x2', 'x8')
        e('madd', 'x8', rb, 'x3', 'x8')
    e('mul', 'x0', 'x8', 'x24'); e('cmp', 'x0', 'x25'); e('bc', 'cs', 'tt_t2')
    pr.label('tt_sel')
    e('umulh', 'x1', 'x8', 'x24'); e('add_lsl', 'x2', 'x20', 'x1', 3); e('ldr', 'x20', 'x2', CH)
    e('cbz', 'x20', 'tt_wend')
    e('subi', 'x23', 'x23', 1); e('bc', 'ne', 'tt')

    def leave(code):
        for a, b in reversed(SAVED):
            e('ldp_post', a, b, 'sp', 16)
        e('movi', 'x0', code); e('ret')
    leave(3)
    pr.label('tt_t2')
    n1 = len(pr.code)
    for j in range(NDIG // 2):
        if j in res:
            ra, rb = res[j]
            e('ldp', 'x2', 'x3', 'x20', 16 * j)
            e('umulh', 'x2', ra, 'x2'); e('add_r', 'x8', 'x8', 'x2')
            e('umulh', 'x3', rb, 'x3'); e('add_r', 'x8', 'x8', 'x3')
        else:
            e('ldp', 'x0', 'x1', 'x19', 16 * j - DOFF); e('ldp', 'x2', 'x3', 'x20', 16 * j)
            e('umulh', 'x0', 'x0', 'x2'); e('add_r', 'x8', 'x8', 'x0')
            e('umulh', 'x1', 'x1', 'x3'); e('add_r', 'x8', 'x8', 'x1')
    e('ldr', 'x1', 'x19', DOFF); e('mul', 'x0', 'x8', 'x24'); e('cmp', 'x0', 'x1'); e('bc', 'cs', 'tt_unc')
    e('b', 'tt_sel')
    n_t2 = len(pr.code) - n1
    pr.label('tt_unc')
    leave(1)
    pr.label('tt_wend')
    e('mov', 'x10', 'x1'); e('mov', 'x11', 'x2')
    leave(2)
    k = len(SAVED)
    return {"entry": n_entry, "exit_budget": k + 2, "exit_uncertified": k + 2, "exit_window_end": k + 4,
            "tag": 4 * (NDIG // 2) - len(RES_PAIRS), "second": n_t2,
            "second_to_exit": n_t2 - 1}


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
    for k in range(0, 12, 2):                              # two coordinates, ten 51-bit digits, five STP
        regs = (('x4', 'x5', 'x6', 'x7', 'x8'), ('x9', 'x10', 'x11', 'x12', 'x16'))
        for h in (0, 1):
            d0, d1, d2, d3, d4 = regs[h]
            e('ldp', 'x0', 'x1', 'x21', 32 * (k + h)); e('ldp', 'x2', 'x3', 'x21', 32 * (k + h) + 16)
            e('and_r', d0, 'x0', 'x14')
            e('extr', d1, 'x1', 'x0', 51); e('and_r', d1, d1, 'x14')
            e('extr', d2, 'x2', 'x1', 38); e('and_r', d2, d2, 'x14')
            e('extr', d3, 'x3', 'x2', 25); e('and_r', d3, d3, 'x14')
            e('lsr', d4, 'x3', 12)
        o = 40 * k - DOFF
        flat = regs[0] + regs[1]
        for i in range(5):
            e('stp', flat[2 * i], flat[2 * i + 1], 'x19', o + 16 * i)
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
CAP_LOG2 = 12              # exit cap X_cap = 2^12 times the expected-exit bound (Markov: reached with probability <= 2^-12)
UNC_LOG2 = 50              # H2 (a2): steps whose second tag is not certified, at most E_steps / 2^50 in expectation
SECOND = Fraction(1, 26)   # H2 (a1): steps whose first tag is not certified, at most E_steps / 26 in expectation
GATE = 64                  # the gate: one coin word (16), compare with K, branch, return
ALLOW = 1 << 10            # per exit: dispatch, exact selector and distinguished test, root child, flags
ENTRY = 1 << 20            # per table entry: product, row, words, labels, pointers (machine.md)
REPORT = 1 << 18           # per report: store lookup and insert on 48-word keys
ST_LOG2 = DM_LOG2 = 50     # caps on walk starts and reports
WALK_LOG2 = 97             # trips (steps) of one walk
LO_RUN = 1 - Fraction(1, 1 << CAP_LOG2) - Fraction(1, 1 << 10) - Fraction(1, 1 << 20)   # success of an opened run
PARAM = {"multipliers": RM, "window": LWIN, "dp_bits": DBITS, "digit_bits": DIG, "digits": NDIG,
         "tag_error_bounds": "60 x 2^51 (first tag), 120 (second tag)", "second_tag_rate": "1/26", "mean_allowance": "107/100",
         "exit_cap": "2^12 x (ceil(E_steps / 3500) + floor(E_steps / 2^12) + 1)",
         "gate": "one fresh 64-bit coin before any work; the run opens when it is below K"}
ROWS = ("gate", "steps", "second_tags", "exits", "input_check", "constants", "table", "walk_starts", "reports", "recoveries",
        "output_fixed", "point_output", "pairing_check")


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
    E_exits = -(-E_steps // LWIN) + (E_steps >> UNC_LOG2) + 1
    return E_steps, E_exits, E_exits << CAP_LOG2


# ---------------------------------------------------------------- kernels

def write_entry(mc, tg, addr, Te, al, be, children):
    """Table entry at addr (the pointer e): first and second words, child pointers, labels, T_e in stored form."""
    for i, w in enumerate(tg.words(Te)):
        mc.st(addr - EOFF + 8 * i, w)
    for i, w in enumerate(tg.words2(Te)):
        mc.st(addr + 8 * i, w)
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
    t12 = gen_f12_mul(pr); tt = gen_tt(pr); gen_rebase(pr)
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
    assert X_cap >> 192 == 0                               # the cap fits the 3-word exit counter
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
        for i, d in enumerate(tg.digits(Y)): mc.st(DGB - DOFF + 8 * i, d)
        mc.st(DGB + DOFF, THR2)

    def lemma(Y, Te, A):
        """(2^64 t / p - A) mod 2^64, which lies in [0, E1) for the first tag and in [0, E2) for the second."""
        t = tg.stored_t(TW.mul(Y, Te))
        d = (Fraction(t << 64, P) - A) % (1 << 64)
        return t, d
    second_trips, lemma2_max = 0, 0
    for walk in range(3):
        n = 3 + walk                                       # window of n steps, then the window end
        while True:                                        # a start whose first n + 1 tags are all certified
            a0, b0 = rng.randrange(R), rng.randrange(R)
            Y0 = TW.mul(TW.pow(g, a0), TW.pow(z, b0))
            dg = tg.digits(Y0)
            Te, al, be, sels, addr = TW.one, 0, 0, [], TAB + (walk << 20)
            chain, allcert, sec = [], True, []
            for i in range(n + 1):                         # entries along the walk's own path; others null
                A = tg.tag(dg, tg.words(Te))
                t, d = lemma(Y0, Te, A)
                lemma_max = max(lemma_max, d)
                s_ex, dist = TY.exact_select(t, P, RM, DBITS)
                cert = (RM * A) & M64 < THR
                wbad += not (0 <= d < ERR) or (cert and (dist or (RM * A) >> 64 != s_ex))
                sec.append(not cert)
                if not cert:                               # the second tag
                    A = tg.tag2(dg, tg.words2(Te), A)
                    t, d = lemma(Y0, Te, A)
                    lemma2_max = max(lemma2_max, d)
                    cert = (RM * A) & M64 < THR2
                    wbad += not (0 <= d < ERR2) or (cert and (dist or (RM * A) >> 64 != s_ex))
                allcert &= cert
                chain.append((addr + 2048 * i, Te, al, be, s_ex))
                sels.append(s_ex)
                Te, al, be = TW.mul(Te, Ms[s_ex]), (al + cd[s_ex][0]) % R, (be + cd[s_ex][1]) % R
            if allcert and (walk > 0 or any(sec)):         # the first window uses the second tag at least once
                break
        for k, v in enumerate(mem_order(Y0)): setc(BUF1 + 32 * k, mont(v))
        put_digits(Y0)
        for i, (ad, T_, al_, be_, s_) in enumerate(chain):
            write_entry(mc, tg, ad, T_, al_, be_, {s_: chain[i + 1][0]} if i < n else {})
        mc.r['x19'], mc.r['x20'] = DGB, chain[0][0]
        keep = {r_: rng.getrandbits(64) for r_ in ('x21', 'x22', 'x27', 'x29', 'x30')}
        mc.r.update(keep); sp0 = mc.r['sp']
        c_run = mc.call('tt_run')
        trips += n + 1
        second_trips += sum(sec)
        rest = c_run - tt["entry"] - (tt["tag"] + 3 + 4) - tt["exit_window_end"] - tt["second"] * sum(sec)   # n full trips remain
        trip_counts.add(rest // n if rest % n == 0 else -1)
        e_last, T_last, al_last, be_last, s_last = chain[n]
        wbad += (mc.r['x0'] != 2 or mc.r['x10'] != s_last or mc.r['x11'] != e_last + 8 * s_last
                 or mc.r['x23'] != LWIN + 1 - n or any(mc.r[r_] != v_ for r_, v_ in keep.items())
                 or mc.r['x19'] != DGB or mc.r['x28'] != CONST or mc.r['sp'] != sp0)
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
                 or [mc.ld(DGB - DOFF + 8 * i) for i in range(NDIG)] != tg.digits(Y1) or mc.ld(DGB + DOFF) != THR2
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
    # (RM * A mod 2^64, increment of the second tag, expected code, expected selector or None)
    crafted = ((THR - 1, 0, 2, None), (0, 0, 2, 0), (THR, 0, 2, None), (THR + 1, 0, 2, None), (THR2 - 1, 0, 2, None),
               (THR2, 0, 1, None), ((1 << 64) - RM, 0, 1, None), ((1 << 64) - RM, 1, 2, 0))
    c_cert = {}
    for lo, inc, code, sel in crafted:
        A = lo * pow(RM, -1, 1 << 64) & M64                # RM * A mod 2^64 = lo (RM odd)
        assert (RM * A) & M64 == lo
        ad = TAB + (9 << 20)
        for i in range(NDIG): mc.st(ad - EOFF + 8 * i, 0); mc.st(ad + 8 * i, 0)
        mc.st(ad - EOFF, A * pow(d0, -1, 1 << 64) & M64)
        w2 = -(-(inc << 64) // d0)                         # UMULH(d0, w2) = inc
        assert (d0 * w2) >> 64 == inc and w2 < 1 << 64
        mc.st(ad, w2)
        for s_ in range(RM): mc.st(ad + CH + 8 * s_, 0)
        mc.r['x19'], mc.r['x20'] = DGB, ad
        c_ = mc.call('tt_run')
        A2 = (A + inc) & M64
        first = lo < THR
        wbad += mc.r['x0'] != code or mc.r['x28'] != CONST
        if code == 1:
            wbad += mc.r['x8'] != A2 or c_ != tt["entry"] + tt["tag"] + 3 + tt["second_to_exit"] + tt["exit_uncertified"]
        else:
            wbad += sel is not None and mc.r['x10'] != sel
            c_cert[first] = c_
    second = c_cert[False] - c_cert[True]                  # executed instructions of the second tag
    ad = TAB + (10 << 20)                                  # a ring of two entries: the budget ends the call
    for j in (0, 1):                                       # all words zero: tag 0, certified, selector 0
        for i in range(NDIG): mc.st(ad + 2048 * j - EOFF + 8 * i, 0); mc.st(ad + 2048 * j + 8 * i, 0)
        for s_ in range(RM): mc.st(ad + 2048 * j + CH + 8 * s_, ad + 2048 * (1 - j))
    mc.st(CONST + WPOOL + 24, 5)
    mc.r['x19'], mc.r['x20'] = DGB, ad
    c_bud = mc.call('tt_run')
    wbad += mc.r['x0'] != 3 or mc.r['x23'] != 0
    mc.st(CONST + WPOOL + 24, LWIN + 1)
    step = trip_counts.pop() if len(trip_counts) == 1 else -1
    wbad += (step != tt["tag"] + 9 or c_bud != tt["entry"] + 5 * step + tt["exit_budget"] or len(reb_counts) != 1
             or second != tt["second"] or second_trips == 0)
    regs = set()
    for op, a in pr.code:
        regs |= {x for x in a if isinstance(x, str) and x.startswith('x') and x[1:].isdigit()}
    return {"fp_mul": cnt['fp_mul'], "fp_add": cnt['fp_add'], "fp_sub": cnt['fp_sub'],
            "fr_mul": cnt['fr_mul'], "fr_add": cnt['fr_add'], "fp_mismatches": bad,
            "f12_mul_executed": n12, "f12_calls": t12.calls, "f12_scratch_slots": t12.peak, "f12_mismatches": bad12,
            "step_executed": step, "second_tag_executed": second, "second_tag_trips_in_windows": second_trips,
            "largest_second_tag_offset": int(lemma2_max) + 1, "tt_entry_exit": tt["entry"] + max(tt["exit_budget"], tt["exit_uncertified"], tt["exit_window_end"]),
            "tt_entry": tt["entry"], "tt_exits": [tt["exit_uncertified"], tt["exit_window_end"], tt["exit_budget"]],
            "rebase_executed": c_reb,
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
    rows = lem = sel = unc = lem2 = sel2 = unc2 = 0
    worst = worst2 = 0
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
        A2 = tg.tag2(tg.digits(Y), tg.words2(Te), A)        # the second tag, on every pair
        d2 = (Fraction(t << 64, P) - A2) % (1 << 64)
        lem2 += 0 <= d2 < ERR2
        worst2 = max(worst2, d2)
        if (RM * A2) & M64 < THR2:
            sel2 += (RM * A2) >> 64 == s_ex and not dist
        else:
            unc2 += 1; sel2 += 1
    # Frobenius structure of the selector coordinate (proof.md section 7, H2): on GT, Y -> Y^(p^2) has order 6
    # and fixes the constant Fp2 coefficient, so t(Y^(p^2)) = t(Y); the same holds in the toy fields.
    fro, nfro = 0, 3
    for i in range(nfro):
        Y = TW.pow(g, rng.randrange(1, R))
        Yf = TW.pow(Y, P * P % R)
        fro += Yf != Y and TW.pow(Yf, R) == TW.one and tg.stored_t(Yf) == tg.stored_t(Y) and pow(P * P % R, 6, R) == 1
    q_, rp_, a0_, T_, g_ = TY.toy_field()
    tq = TY.Tags(T_)
    Yt = T_.pow(g_, rng.randrange(2, rp_))
    Ytf = T_.pow(Yt, q_ * q_ % rp_)
    fro_toy = Ytf != Yt and tq.stored_t(Ytf) == tq.stored_t(Yt)
    ent = math.comb(LWIN + RM, RM)
    small = all(c in (0, 1, P - 1, 9, P - 9) for row in tg.mat for c in row)   # structure constants of the row
    out.update({"row_constants_in_0_1_9": small, "row_identity": f"{rows}/{n}", "tag_lemma": f"{lem}/{n}", "certified_selector_exact": f"{sel}/{n}",
                "uncertified_cases": unc, "largest_tag_offset_log2": round(math.log2(max(1, worst)), 2),
                "tag_error_bound_log2": round(math.log2(ERR), 4),
                "threshold": THR, "uncertified_fraction_log2": round(math.log2(((1 << 64) - THR) / 2 ** 64), 4),
                "second_tag_lemma": f"{lem2}/{n}", "second_certified_selector_exact": f"{sel2}/{n}",
                "second_uncertified_cases": unc2, "largest_second_tag_offset": int(worst2) + 1, "second_tag_error_bound": ERR2,
                "second_threshold": THR2, "second_uncertified_fraction_log2": round(math.log2(((1 << 64) - THR2) / 2 ** 64), 4),
                "table_entries": ent, "digits_cover_p": DIG * (NDIG // 12) >= P.bit_length(),
                "frobenius_p2_fixes_t": f"{fro}/{nfro}", "frobenius_q2_fixes_t_toy": fro_toy,
                "algebra_ok": rows == n and lem == n and sel == n and small and DIG * (NDIG // 12) >= P.bit_length()
                and fro == nfro and fro_toy and lem2 == n and sel2 == n
                and Fraction((1 << 64) - THR, 1 << 64) <= SECOND and ((1 << 64) - THR2) <= 1 << (64 - UNC_LOG2)})
    return out


# ---------------------------------------------------------------- ledger

def ledger(k):
    getcontext().prec = 120
    C = 20000
    F12, FRM = k["f12_mul_executed"], k["fr_mul"] + 4
    STEP, SEC = k["step_executed"], k["second_tag_executed"]
    EXIT = k["rebase_executed"] + k["tt_entry_exit"] + ALLOW
    E_steps, E_exits, X_cap = step_bounds()
    ent = math.comb(LWIN + RM, RM)
    ST, DM, REC = 1 << ST_LOG2, 1 << DM_LOG2, 16
    ob = OS.bounds(k["fp_mul"], F12)
    other = {"input_check": 508 * F12 + (1 << 16),
             "constants": 1 << 22,
             "table": ent * ENTRY + RM * (1017 * F12 + 2 * SCALAR),
             "walk_starts": ST * (1017 * F12 + 2 * SCALAR + (1 << 10)),
             "reports": DM * REPORT,
             "recoveries": REC * (508 * F12 + 512 * FRM + (1 << 10)),
             "output_fixed": (1 << 16) + (1 << 32),
             "point_output": 1 << 22, "pairing_check": 1 << 28}
    E_second = -(-E_steps * SECOND.numerator // SECOND.denominator)    # H2 (a1): expected second tags, rounded up
    exp_rows = dict(other, steps=STEP * E_steps, second_tags=SEC * E_second, exits=EXIT * E_exits)
    cap_rows = dict(other, steps=STEP * (LWIN + 1) * X_cap, second_tags=SEC * (LWIN + 1) * X_cap, exits=EXIT * X_cap)
    E_run, cap_run = sum(exp_rows.values()), sum(cap_rows.values())
    lim = Fraction(39, 100) * (1 << 64) / LO_RUN
    K = lim.numerator // lim.denominator + (lim.numerator % lim.denominator != 0)   # least K with (K/2^64) LO_RUN >= 0.39
    assert Fraction(K, 1 << 64) * LO_RUN >= Fraction(39, 100) > Fraction(K - 1, 1 << 64) * LO_RUN and K < 1 << 64
    exp_total = GATE + Fraction(K * E_run, 1 << 64)        # expected total work: gate + P(open) x E_run
    every_run = GATE + cap_run
    rep_mean = Fraction(E_steps, 1 << (DBITS - 1))          # H2 (b): expected reports
    start_mean = rep_mean + Fraction(E_steps, 1 << (WALK_LOG2 - 1)) + 2   # walks end at a report, the walk cap or a start equal to 1
    assert start_mean / ST + rep_mean / DM <= Fraction(1, 1 << 10)
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
            "second_tag_instructions": SEC, "E_second": str(E_second),
            "per_step_average": math.ceil(float(Fraction(exp_rows["steps"] + exp_rows["second_tags"] + exp_rows["exits"], E_steps)) * 100) / 100,
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
            "cap_failure_starts_reports_log2": float(lg(start_mean / ST + rep_mean / DM)),
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
                   f"{v['correct']} | {v['wrong']} | {v['at_exit_cap']} | {v['exits'] / tot:.4f} | {v['second_tags'] / tot:.5f} | "
                   f"{v['uncertified'] / tot:.5f} | {v['reports'] / tot:.5f} | {v['check_failures']} of {v['checked']:,} |")
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
                  "%d instructions per walk step" % L["step_instructions"], "%d more per second tag" % L["second_tag_instructions"], fmt(L["exit_instructions"]),
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
                  "THR = %d" % THR, "THR2 = %d" % THR2, "**%d** more executed instructions" % L["second_tag_instructions"],
                  "%s entries" % fmt(L["table_entries"]), "average %.2f instructions" % L["per_step_average"]]
    need_proof += L["table"]
    op, om = OS.doc_lines(k["fp_mul"], k["f12_mul_executed"], 22, 28)
    need_proof += op
    need_mach = [fmt(k["f12_mul_executed"]), "**%d**" % L["step_instructions"], "**%d**" % L["second_tag_instructions"], fmt(k["rebase_executed"]),
                 fmt(L["exit_instructions"]), fmt(k["stack_bytes"]) + " bytes", str(k["fp_mul"]),
                 fmt(L["peak_bytes"]) + " = 2^" + L["memory_log2_bytes"], "%s entries" % fmt(L["table_entries"])] + om
    toy = json.loads((HERE / "toy_offline.json").read_text())
    need_proof += toy_rows(toy)
    rows_, (n_, mean, se, mx) = toy_stats(toy)
    plain = [s_ / math.sqrt(math.pi * toy[r_[0]]["r_prime"] / 2) for r_ in rows_ for s_ in r_[7]["steps"]]
    h1_ok = mean + 5 * se <= float(C_MEAN) < sum(plain) / n_ - 5 * se     # H1 support and its negative control
    h2_ok = all(r_[7]["check_failures"] == 0 and r_[7]["wrong"] == 0 for r_ in rows_)
    miss = ([x for x in need_claim if x not in st] + [x for x in need_proof if x not in pr_]
            + [x for x in need_mach if x not in mach])
    allowed = {L["time_log2"], L["every_run_log2"], L["opened_run_expected_log2"], L["memory_log2_bytes"], L["preprocessing_log2"]}
    for name, txt in (("claim", st), ("proof.md", pr_), ("machine.md", mach), ("README.md", (HERE / "README.md").read_text())):
        miss += ["%s: 2^%s" % (name, t) for t in re.findall(r"2\^(\d+\.\d{6})", txt) if t not in allowed]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"\d+/1048576", txt) if t != L["opened_run_success_lower"]]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"time_log2 = (\d+\.\d+)", txt) if t != L["time_log2"]]
    return {"claim_score_equals_ledger": claimed == Decimal(L["time_log2"]), "missing_figures": miss,
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
        toy_bad = bool(t_["wrong_outputs"] or t_["selector_check_failures"])
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

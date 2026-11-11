"""Single declared experiment for proof.md and machine.md (standard library only).

Prints one JSON object with four parts:
  kernel     - the emitted AArch64 sequences of a64.py executed on the emulator:
               Montgomery product, add and subtract on edge and random inputs,
               the tower Fp12 product against the pinned-basis product, and
               the instruction count of every priced sequence.
  production - the emitted walk step executed on production BN254 values with a
               1024-record table, checked against the reference walk (partition,
               table address, product, labels, invariant Y = g^a z^b), plus the
               distinguished-point and counter branches and the record append.
  ledger     - exact recomputation of the cap, the work and the score, and a
               comparison of every figure the claim declares.
  toy        - the identical walk rule on a small BN-shape target group (same
               field polynomial, tower and Montgomery encoding), fresh targets
               and coins per trial, all outcomes reported.  Toy scale only.
"""
import importlib.util
import json
import math
import random
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
spec = importlib.util.spec_from_file_location("a64", base / "a64.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)

X_BN = 0x44e992b44a6909f1
p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
r = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
G = [int(c, 16) for c in (
    "28c6e04df059260df7d2d2a1f9b5f77676d1939847852c4ed50d2318744c1d5f",
    "17bb74adab1705c26133af1dac87044a3833ac011018e8158da48382bbd2dcd6",
    "0d3bd72f54d742f78ea9e6015c8ea2f2e7fbb728c9c905ec531dcf7de5b246f0",
    "090cb8ee97e091a667af03882b06c3ecb4e437993cbd1b05b98c7f9dfcfe9c40",
    "016b6d855b5cbf76f9829a309db52f5c442f65ae29f996af59d65f85f4afe78a",
    "0a0272204db51dadc0342bd318b9302a44faec12ff500bdd4d4b012ffe45f36f",
    "084f330485b09e866bc2f2ea2b897394deaf3f12aa31f28cb0552990967d4704",
    "27ed208e7a0b55ae6e710bbfbd2fd922669c026360e37cc5b2ab862411536104",
    "2067586885c3318eeffa1938c754fe3c60224ee5ae15e66af6b5104c47c8c5d8",
    "279db296f9d479292532c7c493d8e0722b6efae42158387564889c79fc038ee3",
    "2b03614464f04dd772d86df88674c270ffc8747ea13e72da95e3594468f222c4",
    "108c19d15f9446f744d0f110405d3856d6cc3bda6c4d537663729f5257628417")]
RM = 1 << 256
ONE = [1] + [0] * 11

# ---- ledger constants (machine.md, proof.md)
STEP_CHARGE = 13312
DP_CHARGE = 128
D_MAX = 2 ** 36
FIXED = 2 ** 52
DELTA = 2 ** 102
ADDERS = 1024
DP_BITS = 94
C_UNIT = 20000


def pmul(a, b, m=p):
    t = [0] * 23
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                t[i + j] += x * y
    for k in range(22, 11, -1):
        t[k - 6] += 18 * t[k]
        t[k - 12] -= 82 * t[k]
        t[k] = 0
    return [x % m for x in t[:12]]


def ppow(a, e, m=p):
    o = ONE
    while e:
        if e & 1:
            o = pmul(o, a, m)
        a = pmul(a, a, m)
        e >>= 1
    return o


def to_tower(c, m=p):
    return [((c[i] + 9 * c[i + 6]) % m, c[i + 6] % m) for i in range(6)]


def from_tower(t, m=p):
    c = [0] * 12
    for i, (a, b) in enumerate(t):
        c[i + 6] = b
        c[i] = (a - 9 * b) % m
    return c


def mont(v, m=p):
    return v * RM % m


def partition(c, m=p, adders=ADDERS):
    """Low 10 bits of word 0 of Montgomery(A0.im) = Montgomery(c6)."""
    return mont(c[6] % m, m) & (adders - 1)


def is_dp(c, bits=DP_BITS, m=p):
    """Low `bits` bits of Montgomery(A0.re) = Montgomery(c0 + 9 c6) are zero."""
    return mont((c[0] + 9 * c[6]) % m, m) % (1 << bits) == 0


# ---- emulator plumbing
def machine():
    X = [0] * 31
    for i, reg in enumerate((A.CUR, A.TAB, A.NXT, A.S0, A.S1, A.CTL)):
        X[reg] = 0x100000 * (i + 1) + (0x1000 if reg in (A.S0, A.S1) else 0)
    for j in range(4):
        X[A.P_REGS[j]] = (p >> (64 * j)) & A.MASK
    X[A.N0] = (-pow(p, -1, 1 << 64)) % (1 << 64)
    return X, {}


def put_words(mem, addr, v):
    for k in range(4):
        mem[addr + 8 * k] = (v >> (64 * k)) & A.MASK


def get_words(mem, addr):
    return sum(mem[addr + 8 * k] << (64 * k) for k in range(4))


def put_f12(mem, addr, c):
    t = to_tower(c)
    for k in range(2):
        for m in range(3):
            a, b = t[k + 2 * m]
            put_words(mem, addr + 192 * k + 64 * m, mont(a))
            put_words(mem, addr + 192 * k + 64 * m + 32, mont(b))


def get_f12(mem, addr):
    inv = pow(RM, -1, p)
    t = [None] * 6
    for k in range(2):
        for m in range(3):
            a = get_words(mem, addr + 192 * k + 64 * m)
            b = get_words(mem, addr + 192 * k + 64 * m + 32)
            assert a < p and b < p, "non-canonical output"
            t[k + 2 * m] = (a * inv % p, b * inv % p)
    return from_tower(t)


def kernel(rng):
    out = {}
    cases = [(a, b) for a in (0, 1, 2, p - 2, p - 1) for b in (0, 1, 2, p - 2, p - 1)]
    cases += [(rng.randrange(p), rng.randrange(p)) for _ in range(400)]
    for name, emit, ref, mods in (
            ('fp_mul', A.fp_mul, lambda a, b, m: a * b * pow(RM, -1, m) % m, (p,)),
            ('fp_add', A.fp_add, lambda a, b, m: (a + b) % m, (p, r)),
            ('fp_sub', A.fp_sub, lambda a, b, m: (a - b) % m, (p,))):
        n = 0
        for m in mods:
            mregs = A.P_REGS if m == p else (9, 10, 11, 12)
            P = A.Prog()
            emit(P, (A.S0, 0), (A.S0, 32), (A.S0, 64), mod=mregs)
            for a, b in cases:
                a, b = a % m, b % m
                X, mem = machine()
                for j in range(4):
                    X[9 + j] = (r >> (64 * j)) & A.MASK
                put_words(mem, X[A.S0] + 32, a)
                put_words(mem, X[A.S0] + 64, b)
                A.run(P.ins, X, mem)
                assert get_words(mem, X[A.S0]) == ref(a, b, m), (name, a, b)
                n += 1
        out[name + '_cases_passed'] = n
    P, S = A.Prog(), A.Scratch()
    A.f12_mul(P, S, A.f12_loc(A.NXT), A.f12_loc(A.CUR), A.f12_loc(A.TAB))
    tests = [(ONE, ONE), ([p - 1] * 12, [p - 1] * 12), (G, G)]
    tests += [([rng.randrange(p) for _ in range(12)], [rng.randrange(p) for _ in range(12)]) for _ in range(5)]
    n = 0
    for x, y in tests:
        X, mem = machine()
        put_f12(mem, X[A.CUR], x)
        put_f12(mem, X[A.TAB], y)
        A.run(P.ins, X, mem)
        assert get_f12(mem, X[A.NXT]) == pmul(x, y)
        n += 1
    out['fp12_mul_cases_passed_vs_pinned_basis'] = n
    out['fp12_scratch_slots_peak_of_64'] = S.peak
    out['instruction_counts'] = A.counts()
    guards = []
    for bad in (('ldp', 0, 1, A.S0, 512), ('ldp', 0, 0, A.S0, 0), ('ldr', 0, A.S0, 4),
                ('and_imm', 0, 0, 0x5), ('and_imm', 0, 0, (1 << 64) - 1),
                ('and_imm', 0, 0, (1 << 65) - 1), ('and_imm', 0, 0, -1)):
        try:
            A.Prog().e(*bad)
            guards.append(False)
        except AssertionError:
            guards.append(True)
    out['encoding_guards_reject_bad_forms'] = all(guards)
    prof = json.loads((base / 'target-profile.json').read_text())
    out['constants_match_profile'] = (
        int(prof['field_modulus'], 16) == p and int(prof['subgroup_order'], 16) == r and
        int(prof['parameter_x'], 16) == X_BN and
        [int(c, 16) for c in prof['generator_g']['coefficients']] == G)
    return out


def production(rng):
    out = {}
    g = G
    out['g_order_r'] = g != ONE and ppow(g, r) == ONE
    out['p_r_from_x'] = (36 * X_BN**4 + 36 * X_BN**3 + 24 * X_BN**2 + 6 * X_BN + 1 == p and
                         36 * X_BN**4 + 36 * X_BN**3 + 18 * X_BN**2 + 6 * X_BN + 1 == r)
    x = rng.randrange(1, r)
    z = ppow(g, x)
    # table: valid labels, cheap construction M_j = M_0 * S^j (mechanics test only)
    c0, d0, gm, dl = (rng.randrange(r) for _ in range(4))
    m0 = pmul(ppow(g, c0), ppow(z, d0))
    st = pmul(ppow(g, gm), ppow(z, dl))
    X, mem = machine()
    tab = 0x40000000
    ctl = X[A.CTL]
    put_words(mem, ctl + A.CTL_R, r)
    put_words(mem, ctl + A.CTL_R + 16, r >> 128)
    mem[ctl + A.CTL_TAB] = tab
    labels = []
    mj = m0
    for j in range(ADDERS):
        cj, dj = (c0 + j * gm) % r, (d0 + j * dl) % r
        put_f12(mem, tab + 448 * j, mj)
        put_words(mem, tab + 448 * j + A.REC_C, cj)
        put_words(mem, tab + 448 * j + A.REC_D, dj)
        labels.append((cj, dj))
        mj = pmul(mj, st)
    a, b = rng.randrange(r), rng.randrange(r)
    y = pmul(ppow(g, a), ppow(z, b))
    put_f12(mem, X[A.CUR], y)
    put_words(mem, ctl + A.CTL_A, a)
    put_words(mem, ctl + A.CTL_B, b)
    cap = (1 << 127) + 12345
    mem[ctl + A.CTL_CNT] = 0
    mem[ctl + A.CTL_CNT + 8] = 0
    mem[ctl + A.CTL_CAP] = cap & A.MASK
    mem[ctl + A.CTL_CAP + 8] = cap >> 64
    progs = [A.step_program(A.CUR, A.NXT)[0].ins, A.step_program(A.NXT, A.CUR)[0].ins]
    regs = [A.CUR, A.NXT]
    steps, used = 12, set()
    for s in range(steps):
        cur = regs[s % 2]
        j = partition(y)
        used.add(j)
        br = A.run(progs[s % 2], X, mem)
        assert X[A.TAB] == tab + 448 * j
        y = pmul(y, get_f12(mem, tab + 448 * j))
        a, b = (a + labels[j][0]) % r, (b + labels[j][1]) % r
        nxt = regs[(s + 1) % 2]
        assert get_f12(mem, X[nxt]) == y
        assert get_words(mem, ctl + A.CTL_A) == a and get_words(mem, ctl + A.CTL_B) == b
        assert [t for _, t in br] == [is_dp(y), False, True]
    out['emitted_steps_match_reference'] = steps
    out['distinct_partitions_used'] = len(used)
    out['invariant_Y_eq_g^a_z^b'] = y == pmul(ppow(g, a), ppow(z, b))
    out['counter_after_steps'] = mem[ctl + A.CTL_CNT] == steps
    # counter exit branch
    mem[ctl + A.CTL_CNT] = (cap - 1) & A.MASK
    mem[ctl + A.CTL_CNT + 8] = (cap - 1) >> 64
    br = A.run(progs[steps % 2], X, mem)
    out['cap_exit_branch_taken'] = br[-2][1] is True
    # distinguished-point predicate and record append on crafted words
    hits = []
    for w0, w1 in ((0, 0), (0, 1 << 30), (0, 1), (1, 0), (0, (1 << 64) - (1 << 30))):
        P = A.Prog()
        A.dp_test(P, A.NXT)
        X2, mem2 = machine()
        mem2[X2[A.NXT]] = w0
        mem2[X2[A.NXT] + 8] = w1
        hits.append(A.run(P.ins, X2, mem2)[0][1])
    out['dp_predicate_cases'] = hits == [True, True, False, False, True]
    X2, mem2 = machine()
    put_f12(mem2, X2[A.NXT], y)
    put_words(mem2, X2[A.CTL] + A.CTL_A, a)
    put_words(mem2, X2[A.CTL] + A.CTL_B, b)
    mem2[X2[A.CTL] + A.CTL_NDP] = 5
    mem2[X2[A.CTL] + A.CTL_DPB] = 0x50000000
    mem2[X2[A.CTL] + A.CTL_DPMAX] = D_MAX
    A.run(A.dp_program(A.NXT).ins, X2, mem2)
    rec = 0x50000000 + 5 * 448
    out['dp_record_append'] = (get_f12(mem2, rec) == y and get_words(mem2, rec + 384) == a and
                               get_words(mem2, rec + 416) == b and mem2[X2[A.CTL] + A.CTL_NDP] == 6)
    return out


CAP_NUM, CAP_DEN = 27, 25      # N0 = ceil((27/25) * ceil(sqrt(r)))
PRED_STATED = "0.441572"       # proof.md's lower figure for the H1 prediction at N0


def declared():
    """The shipped claim: its exact score and the bounds stated in its text."""
    claim = json.loads((base / "submission.json").read_text())["claims"][0]
    st = claim["statement"]
    grab = lambda pat: re.search(pat, st).group(1)
    return {
        "score": Fraction(int(claim["score"]["numerator"]), int(claim["score"]["denominator"])),
        "score_text": grab(r"Every run uses at most 2\^(\d+\.\d+) units"),
        "step": int(grab(r"the walk \(([\d,]+) instructions per step").replace(",", "")),
        "emitted": int(grab(r"emitted ([\d,]+)-instruction sequence").replace(",", "")),
        "memory_log2": int(grab(r"Peak memory is below 2\^(\d+) bytes")),
        "preprocessing_log2": int(grab(r"preprocessing is below 2\^(\d+) units")),
        "success": grab(r"success at least (0\.\d+)"),
    }


def ledger():
    getcontext().prec = 90
    dc = declared()
    s = math.isqrt(r - 1) + 1
    n0 = -(-CAP_NUM * s // CAP_DEN)
    cap = n0 + DELTA
    counts = A.counts()
    r_run = cap * STEP_CHARGE + D_MAX * DP_CHARGE + FIXED
    lg = (Decimal(r_run).ln() - Decimal(C_UNIT).ln()) / Decimal(2).ln()
    num = int((lg * 10**6).to_integral_value(rounding='ROUND_CEILING'))
    rho = 1 - Fraction(1, ADDERS)
    pred = 1 - math.exp(-float(rho) * n0 * n0 / (2 * r))
    fail_terms = 2.0 ** -23 + 2.0 ** -255 + 2.0 ** -1000 + 1 / r + 2.0 ** -1012
    floor_ = float(PRED_STATED) - fail_terms
    return {
        'ceil_sqrt_r': str(s), 'N0': str(n0), 'cap': str(cap),
        'step_emitted': counts['step'], 'step_charged': STEP_CHARGE,
        'dp_path_emitted': counts['dp_path'], 'dp_path_charged': DP_CHARGE,
        'charges_cover_emitted': counts['step'] <= STEP_CHARGE and counts['dp_path'] <= DP_CHARGE,
        'R': str(r_run), 'log2_W': float(round(lg, 9)),
        'submitted_score': [str(num), '1000000'],
        'submitted_score_decimal': '%d.%06d' % divmod(num, 10**6),
        'score_is_upper_bound': Decimal(num) / Decimal(10**6) >= lg,
        'heuristic_prediction_at_N0_floor6': '%.6f' % (math.floor(pred * 10**6) / 10**6),
        'proof_prediction_figure_is_lower_bound': pred > float(PRED_STATED),
        'success_floor': '%.6f' % (math.floor(floor_ * 10**6) / 10**6),
        'success_floor_ge_0.39': floor_ >= 0.39,
        'declared_score': '%d/%d' % (dc['score'].numerator, dc['score'].denominator),
        'declared_score_equals_ledger': dc['score'] == Fraction(num, 10**6),
        'declared_score_is_upper_bound': Decimal(dc['score'].numerator) / Decimal(dc['score'].denominator) >= lg,
        'declared_score_text_matches': Fraction(dc['score_text']) == dc['score'],
        'declared_success_is_0.39': dc['success'] == '0.39',
        'declared_step_charge_matches': dc['step'] == STEP_CHARGE and dc['emitted'] == counts['step'],
        'declared_preprocessing_holds': (2**24 + 2**36 + 2**35 < 2**37 and
                                         Fraction(2**37, C_UNIT) <= 2**dc['preprocessing_log2']),
        'declared_memory_holds': (D_MAX * 448 + 2 * D_MAX * 8 + ADDERS * 448 + 2**20
                                  < 2**dc['memory_log2']),
        'stated_fixed_items_below_F': (2**24 + 2**36 + 2**35 + 2**51 + 2**29 + 2**23 +
                                       2**28 + 2**32) < FIXED,
        'proof_success_floor_holds': floor_ > 0.4415,
        'proof_dp_mean_bound_holds': math.log2(cap / 2**DP_BITS) < 32.91,
        'proof_capacity_ratio_holds': D_MAX / (cap / 2**DP_BITS) > 8.5,
        'proof_short_cycle_holds': DELTA * math.sqrt(math.pi * float(rho) / (2 * r)) < 2.0 ** -23,
        'expected_dp_count_log2': round(math.log2(cap / 2**DP_BITS), 4),
        'dp_capacity_over_mean': round(D_MAX / (cap / 2**DP_BITS), 4),
        'dp_table_bytes_log2': round(math.log2(D_MAX * 448), 4),
    }


TOYS = {5: (27631, 27481), -7: (75223, 74929)}


def toy(xt, trials, rng, adders):
    pt, rt = TOYS[xt]
    assert 36 * xt**4 + 36 * xt**3 + 24 * xt**2 + 6 * xt + 1 == pt
    assert 36 * xt**4 + 36 * xt**3 + 18 * xt**2 + 6 * xt + 1 == rt
    e = (pt**12 - 1) // rt
    gen = next(h for h in (ppow([k, 1] + [0] * 10, e, pt) for k in range(2, 50)) if h != ONE)
    assert ppow(gen, rt, pt) == ONE
    dbits = 2
    st = math.isqrt(rt - 1) + 1
    n0 = -(-CAP_NUM * st // CAP_DEN)
    delta = 8 << dbits
    cap = n0 + delta
    # walk rule needs only partition/DP of g^k; tabulate them (exact, bijective k <-> g^k)
    part, dp = [0] * rt, [False] * rt
    y = ONE
    for k in range(rt):
        part[k] = partition(y, pt, adders)
        dp[k] = is_dp(y, dbits, pt)
        y = pmul(y, gen, pt)
    outcomes = {'success': 0, 'cap_or_no_duplicate': 0, 'degenerate': 0, 'wrong': 0}
    for _ in range(trials):
        x = rng.randrange(rt)
        cs = [rng.randrange(rt) for _ in range(adders)]
        ds = [rng.randrange(rt) for _ in range(adders)]
        a, b = rng.randrange(rt), rng.randrange(rt)
        seen, got = {}, None
        for _ in range(cap):
            k = (a + x * b) % rt
            j = part[k]
            a, b = (a + cs[j]) % rt, (b + ds[j]) % rt
            k = (a + x * b) % rt
            if dp[k]:
                if k in seen and got is None:
                    a2, b2 = seen[k]
                    got = (a, b, a2, b2)
                seen.setdefault(k, (a, b))
        if got is None:
            outcomes['cap_or_no_duplicate'] += 1
            continue
        a, b, a2, b2 = got
        if (b - b2) % rt == 0:
            outcomes['degenerate'] += 1
            continue
        xr = (a2 - a) * pow(b - b2, -1, rt) % rt
        if xr == x:
            outcomes['success'] += 1
        else:
            outcomes['wrong'] += 1
    rho = 1 - 1 / adders
    return {'toy_x': xt, 'adders': adders, 'p': pt, 'r': rt, 'dp_bits': dbits, 'N0': n0, 'delta': delta, 'cap': cap,
            'trials': trials, 'outcomes': outcomes,
            'run_success_fraction': round(outcomes['success'] / trials, 4),
            'prediction_at_N0': round(1 - math.exp(-rho * n0 * n0 / (2 * rt)), 4),
            'prediction_at_cap': round(1 - math.exp(-rho * cap * cap / (2 * rt)), 4)}


def main():
    req = json.load(sys.stdin)
    rng = random.Random(req['seed'])
    print(json.dumps({'kernel': kernel(rng), 'production': production(rng), 'ledger': ledger(),
                      'toy': toy(req['toy_x'], req['trials'], rng, req['adders'])}))


if __name__ == '__main__':
    main()

"""Exact Cheon reduction with baby-step giant-step: production ledger and an exhaustive toy run.

Part 1 recomputes every production number quoted in proof.md, README.md, submission.json and claim.json
(group orders, the step counts m and G, the exact success fraction, the itemised instruction ledger, W(2^20),
the score, memory and preprocessing) and checks the production facts the proof uses (factorisation of r-1,
the generator 7, the cube-root pairing of the G1 endomorphism on P0).

Part 2 runs the same algorithm, step for step, on a small curve y^2 = x^3 + b with a prime-order subgroup of
order rt, d a power of two dividing rt-1 and 6 dividing (rt-1)/d: three-window fixed-base tables for P0 and A_d,
incomplete affine additions that raise on any exceptional case, the canonical x-key under the six units,
an LSD radix sort of all records on the full key, one adjacent-pair scan, the six-unit check against A_d and the
d-candidate search against T. It runs every secret s in 1..rt-1 and checks that the number of recovered secrets
equals d * 6 * m * G exactly, that every output satisfies [s]P0 = T, and that no exceptional addition occurs.
Small experiments are evidence at their own scope only. Standard library only."""
import json
import math
import sys
from fractions import Fraction

# ---------------------------------------------------------------- production constants (target-profile.json)
P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16)
R = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
GX = int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16)
GY = int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16)
D = 2 ** 20
FACTORS_R_MINUS_1 = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2,
                     2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}
ZETA = 7            # generator of F_r^*, checked below
UNIT = 2000         # instructions per reference operation (C for bls12-381-g1-ecdlp)
FP_MUL = 383        # measured blst_fp_mul instructions (COST_MODEL.md)
WIN = 85            # window width: three windows cover 255 bits, and r < 2^255
LANES = 2 ** 16     # lanes per batch (Montgomery batch inversion across lanes)
DIGIT = 16          # radix-sort digit width
KEY_BITS = 384      # sort key width: canonical x < p < 2^381, padded to 6 words
REC_BYTES = 64      # record: 6-word key + 2-word index with the list tag in the top bit
TAB_BYTES = 96      # packed affine table entry (x, y), 12 words, no padding

# Per step (baby or giant), instructions outside the two reference additions; worst case, no early exit.
STEP_ITEMS = [
    ("scalar update: one F_r Montgomery multiplication by the stored constant kappa*2^256 mod r", 383),
    ("split the 4-word scalar into three 85-bit two-word digits, 8 per digit", 24),
    ("zero test per digit (orr, compare, branch), 3 per digit", 9),
    ("table address per digit: two-word index times 96 plus window base, with carries, 10 per digit", 30),
    ("table reads: 3 entries x 12 words x 2 per large-store access", 72),
    ("operand staging into and out of the batch buffers, 60 per addition, 2 additions", 120),
    ("beta*x: one F_p Montgomery multiplication", 383),
    ("beta^2*x = -(x + beta*x) mod p: modular addition 19, modular negation 21", 40),
    ("canonical reduction of x and beta*x to [0,p) (6 subtractions, 6 selects, 1 compare each), kept even if already reduced", 26),
    ("minimum of three 6-word values: two compare-and-select of 12", 24),
    ("record write: 8 words x 2 per large-store access 16, two-word output pointer increment 2", 18),
    ("loop control: two-word counter increment 2, two-word compare 2, branch 1", 5),
]
# Per record per radix pass (LSD counting sort, stable), instructions.
HIST_ITEMS = [
    ("read the key word holding the digit (large store)", 2),
    ("extract the 16-bit digit (shift, and)", 2),
    ("counter address", 1),
    ("two-word counter: load 2, add with carry 2, store 2", 6),
    ("loop control (two-word pointer increment, two-word compare, branch)", 5),
]
SCAT_ITEMS = [
    ("read the key word holding the digit (large store)", 2),
    ("extract the digit", 2),
    ("bucket pointer address", 1),
    ("two-word bucket pointer load", 2),
    ("copy the 8-word record: 8 reads and 8 writes at 2 per large-store access", 32),
    ("two-word bucket pointer increment", 2),
    ("two-word bucket pointer store", 2),
    ("loop control", 5),
]
SCAN_ITEMS = [
    ("read the 6-word key (large store)", 12),
    ("compare with the previous key: 6 eor, 5 orr, 1 branch", 12),
    ("read the tag word and test the tag", 4),
    ("keep the previous record pointer (two words)", 2),
    ("loop control", 5),
]
BUCKET_FIXED = 8    # per bucket per pass: clear 2, two-word prefix add 2, store 2, loop 2


def total(items):
    return sum(c for _, c in items)


def log2q(x):
    x = Fraction(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


def ceil_div(a, b):
    return -(-a // b)


def is_prime_small(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


# ---------------------------------------------------------------- elliptic-curve arithmetic (affine, y^2 = x^3 + b)
class Exceptional(Exception):
    pass


def add_incomplete(p, A, B):
    """Affine addition A + B for A, B finite with x(A) != x(B); anything else is an exceptional case and raises."""
    if A is None or B is None or A[0] == B[0]:
        raise Exceptional((A, B))
    lam = (B[1] - A[1]) * pow(B[0] - A[0], -1, p) % p
    x = (lam * lam - A[0] - B[0]) % p
    return (x, (lam * (A[0] - x) - A[1]) % p)


def dbl(p, A):
    """Affine doubling of a finite point with y != 0 (odd prime order subgroup, so y != 0)."""
    if A is None or A[1] == 0:
        raise Exceptional((A, A))
    lam = 3 * A[0] * A[0] * pow(2 * A[1], -1, p) % p
    x = (lam * lam - 2 * A[0]) % p
    return (x, (lam * (A[0] - x) - A[1]) % p)


def add_complete(p, A, B):
    """Setup-only addition (input generation and reference checks); never used inside the attack's steps."""
    if A is None:
        return B
    if B is None:
        return A
    if A[0] == B[0]:
        return dbl(p, A) if A[1] == B[1] else None
    return add_incomplete(p, A, B)


def mul_ref(p, k, A):
    out = None
    for bit in bin(k)[2:] if k > 0 else "":
        out = add_complete(p, out, out)
        if bit == "1":
            out = add_complete(p, out, A)
    return out


# ---------------------------------------------------------------- production facts and ledger
def production():
    out = {}
    n = 1
    for q, e in FACTORS_R_MINUS_1.items():
        assert is_prime_small(q)
        n *= q ** e
    assert n == R - 1, "factorisation of r-1"
    assert FACTORS_R_MINUS_1[2] == 32 and (R - 1) % D == 0
    assert all(pow(ZETA, (R - 1) // q, R) != 1 for q in FACTORS_R_MINUS_1), "7 generates F_r^*"
    N = (R - 1) // D
    assert N % 6 == 0
    M = N // 6
    assert R < 2 ** 255 and 3 * WIN >= 255
    # cube roots: beta in F_p, lambda in F_r, paired by phi(P0) = [lambda]P0
    lam = pow(ZETA, (R - 1) // 3, R)
    c = 2
    while pow(c, (P - 1) // 3, P) == 1:
        c += 1
    beta = pow(c, (P - 1) // 3, P)
    P0 = (GX, GY)
    assert (GY * GY - GX ** 3 - 4) % P == 0
    if mul_ref(P, lam, P0) != (beta * GX % P, GY):
        beta = beta * beta % P
    assert mul_ref(P, lam, P0) == (beta * GX % P, GY), "endomorphism pairing on P0"
    assert mul_ref(P, R, P0) is None, "P0 has order r"
    # step counts: m = ceil(sqrt(2M/5)), G = ceil(2M/(5m)), success m*G/M >= 2/5 exactly
    m = math.isqrt(2 * M // 5)
    while m * m * 5 < 2 * M:
        m += 1
    G = ceil_div(2 * M, 5 * m)
    assert m * G <= M
    succ = Fraction(m * G, M)
    assert succ >= Fraction(2, 5)
    step_ins = total(STEP_ITEMS)
    passes = ceil_div(KEY_BITS, DIGIT)
    sort_ins = passes * (total(HIST_ITEMS) + total(SCAT_ITEMS)) + total(SCAN_ITEMS)
    per_record = Fraction(2) + Fraction(step_ins + sort_ins, UNIT)
    inv_ins = (P.bit_length() - 1 + bin(P - 2).count("1") - 1) * FP_MUL + 50   # Fermat inversion in F_p
    dbl_ins = inv_ins + 6 * FP_MUL + 100                                         # unbatched affine doubling or addition
    exp_r_ins = 2 * 256 * FP_MUL                                                 # any F_r exponentiation
    rounds = ceil_div(m, LANES) + ceil_div(G, LANES)
    entries = 2 ** WIN - 1
    tables = 6
    t_entry_ins = 60 + 24 + 6                                                    # staging, 12-word write at 2, address/loop
    items = {}
    items["steps"] = (m + G) * (Fraction(2) + Fraction(step_ins, UNIT))
    items["sort_and_scan"] = (m + G) * Fraction(sort_ins, UNIT) + Fraction(passes * (2 ** DIGIT) * BUCKET_FIXED, UNIT)
    items["batch_inversions"] = Fraction(2 * rounds * inv_ins, UNIT)
    items["lane_starts"] = Fraction(2 * LANES * (exp_r_ins + 64), UNIT)
    items["constants"] = Fraction(8 * exp_r_ins + 12 + 2 * FP_MUL, UNIT)
    prep = (tables * entries * (1 + Fraction(t_entry_ins, UNIT))
            + Fraction(tables * (LANES * 170 * dbl_ins + ceil_div(entries - 2, LANES) * inv_ins + dbl_ins), UNIT)
            + Fraction(2 * 170 * dbl_ins, UNIT))
    items["preprocessing_tables"] = prep
    cand_ins = step_ins + 25
    items["phase2"] = (6 * (2 + Fraction(cand_ins + FP_MUL, UNIT)) + D * (2 + Fraction(cand_ins, UNIT))
                       + Fraction(60 + 3 * exp_r_ins + 10, UNIT))
    W = sum(items.values())
    lw = log2q(W)
    score = Fraction(math.ceil(lw * 100), 100)
    assert log2q(W) < float(score) and float(score) - lw > 1e-6
    mem = tables * entries * TAB_BYTES + 2 * (m + G) * REC_BYTES + 2 ** 30
    out.update({
        "log2_N": round(log2q(N), 6), "log2_M": round(log2q(M), 6), "M": str(M), "m": str(m), "G": str(G),
        "log2_m": round(log2q(m), 6), "log2_steps": round(log2q(m + G), 6),
        "success": str(succ), "success_decimal": "%.9f" % float(succ),
        "step_instructions": step_ins, "sort_passes": passes, "sort_scan_instructions": sort_ins,
        "per_record_units": str(per_record), "per_record_units_decimal": "%.6f" % float(per_record),
        "fp_inversion_instructions": inv_ins, "unbatched_group_op_instructions": dbl_ins,
        "items_log2": {k: round(log2q(v), 6) for k, v in items.items()},
        "log2_W": "%.6f" % lw, "score": "%d/%d" % (score.numerator, score.denominator),
        "preprocessing_log2": round(log2q(prep), 6),
        "table_entries": tables * entries, "table_bytes_log2": round(log2q(tables * entries * TAB_BYTES), 6),
        "record_bytes_log2": round(log2q(2 * (m + G) * REC_BYTES), 6),
        "memory_bytes_log2": round(log2q(mem), 6),
        "beta": hex(beta), "lambda": hex(lam),
    })
    return out


# ---------------------------------------------------------------- toy run of the same algorithm
TOY_CURVES = [  # (p, b, cofactor h, subgroup order rt, d)
    (2731, 3, 1, 2689, 8),
    (3187, 2, 1, 3217, 8),
    (3391, 7, 3, 1153, 8),
]


def prime_factors(n):
    fs, q = [], 2
    while q * q <= n:
        if n % q == 0:
            fs.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        fs.append(n)
    return fs


def toy(curve_index, target_permille):
    p, b, h, rt, d = TOY_CURVES[curve_index]
    assert is_prime_small(p) and is_prime_small(rt) and p % 3 == 1 and (rt - 1) % d == 0
    N = (rt - 1) // d
    assert N % 6 == 0
    M = N // 6
    x = 0
    while True:
        rhs = (x ** 3 + b) % p
        ys = [y for y in range(p) if y * y % p == rhs]
        if ys:
            P0 = mul_ref(p, h, (x, min(ys)))
            if P0 is not None and mul_ref(p, rt, P0) is None:
                break
        x += 1
    zeta = 2
    while any(pow(zeta, (rt - 1) // q, rt) == 1 for q in prime_factors(rt - 1)):
        zeta += 1
    g = pow(zeta, d, rt)                     # generator of E, order N
    lam = pow(zeta, (rt - 1) // 3, rt)
    beta = next(bb for bb in range(2, p) if pow(bb, 3, p) == 1 and mul_ref(p, lam, P0) == (bb * P0[0] % p, P0[1]))
    w = -(-rt.bit_length() // 3)             # three windows cover every scalar below rt
    m = math.isqrt(target_permille * M // 1000)
    while m * m * 1000 < target_permille * M:
        m += 1
    G = min(-(-target_permille * M // (1000 * m)), M // m)
    assert m * G <= M

    def key(Q):
        x1 = Q[0]
        x2 = beta * x1 % p
        return min(x1, x2, (-x1 - x2) % p)

    # key lemma, checked on every subgroup point: equal keys <=> same class under the six units
    pts, Q = [], P0
    for _ in range(rt - 1):
        pts.append(Q)
        Q = add_complete(p, Q, P0)
    groups = {}
    for Q in pts:
        groups.setdefault(key(Q), set()).add(Q)
    for grp in groups.values():
        Q = next(iter(grp))
        cls = {Q, (Q[0], -Q[1] % p), (beta * Q[0] % p, Q[1]), (beta * Q[0] % p, -Q[1] % p),
               (beta * beta * Q[0] % p, Q[1]), (beta * beta * Q[0] % p, -Q[1] % p)}
        assert grp == cls and len(cls) == 6
    assert len(groups) == (rt - 1) // 6

    def tables(B):
        tab = []
        Bk = B
        for k in range(3):
            if k:
                for _ in range(w):
                    Bk = dbl(p, Bk)
            row = [None, Bk, dbl(p, Bk)]
            for c in range(3, 2 ** w):
                row.append(add_incomplete(p, row[-1], Bk))   # c-1 >= 2, so never +-Bk
            tab.append(row)
        return tab

    def fixed_base(tab, e):
        assert 0 < e < rt
        acc = None
        for k in range(3):
            c = (e >> (w * k)) & (2 ** w - 1)
            if c == 0:
                continue
            acc = tab[k][c] if acc is None else add_incomplete(p, acc, tab[k][c])
        return acc

    def radix_sort(recs, bits, dg):
        for shift in range(0, bits, dg):
            buckets = [[] for _ in range(2 ** dg)]
            for rec in recs:
                buckets[(rec[0] >> shift) & (2 ** dg - 1)].append(rec)
            recs = [rec for bk in buckets for rec in bk]
        return recs

    g_inv_m = pow(g, -m, rt)
    omega = pow(zeta, N, rt)
    kbits = p.bit_length()

    def attack(T, Ad):
        tP, tA = tables(P0), tables(Ad)
        recs, e = [], 1
        for j in range(m):
            recs.append((key(fixed_base(tP, e)), j, 0))
            e = e * g % rt
        e = 1
        for i in range(G):
            recs.append((key(fixed_base(tA, e)), i, 1))
            e = e * g_inv_m % rt
        recs = radix_sort(recs, kbits, 4)
        match = None
        for a, bb in zip(recs, recs[1:]):
            if a[0] == bb[0]:
                assert a[2] != bb[2] and match is None     # babies distinct, giants distinct, at most one pair
                match = (a, bb)
        if match is None:
            return None, None
        jj = match[0][1] if match[0][2] == 0 else match[1][1]
        ii = match[0][1] if match[0][2] == 1 else match[1][1]
        u0 = jj + ii * m
        hits = [k for k in range(6) if fixed_base(tP, pow(g, u0 + k * M, rt)) == Ad]
        assert len(hits) == 1
        u = (u0 + hits[0] * M) % N
        s0 = pow(zeta, u, rt)
        found = [l for l in range(d) if fixed_base(tP, s0 * pow(omega, l, rt) % rt) == T]
        assert len(found) == 1
        return s0 * pow(omega, found[0], rt) % rt, u0

    dlog = {}
    e = 1
    for u in range(N):
        dlog[e] = u
        e = e * g % rt
    ok = fail = exc = 0
    predicted_ok = 0
    for s in range(1, rt):
        T = mul_ref(p, s, P0)
        Ad = mul_ref(p, pow(s, d, rt), P0)
        u_true = dlog[pow(s, d, rt)]
        predicted = (u_true % M) < m * G
        predicted_ok += predicted
        try:
            out, _ = attack(T, Ad)
        except Exceptional:
            exc += 1
            continue
        if out is None:
            fail += 1
            assert not predicted
        else:
            assert out == s and mul_ref(p, out, P0) == T and predicted
            ok += 1
    formula = d * 6 * m * G
    return {
        "curve": {"p": p, "b": b, "cofactor": h, "rt": rt, "d": d, "N": N, "M": M, "window_bits": w},
        "m": m, "G": G, "secrets_run": rt - 1, "recovered": ok, "not_covered": fail, "exceptional_additions": exc,
        "formula_d_6_m_G": formula, "recovered_equals_formula": ok == formula and ok == predicted_ok,
        "success_fraction": str(Fraction(ok, rt - 1)), "m_G_over_M": str(Fraction(m * G, M)),
    }


def main():
    raw = sys.stdin.read()
    args = json.loads(raw) if raw.strip() else {}
    ci = int(args.get("curve", 0))
    tp = int(args.get("target_permille", 400))
    if not (0 <= ci < len(TOY_CURVES)) or not (100 <= tp <= 800):
        raise SystemExit("curve in 0..2, target_permille in 100..800")
    res = {"production": production(), "toy": toy(ci, tp)}
    print(json.dumps(res, indent=1))
    if not res["toy"]["recovered_equals_formula"] or res["toy"]["exceptional_additions"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

"""Cheon-type auxiliary-powers attack on BLS12-381 G1, with a GLV-priced scalar multiplication.

Stdlib only. Three separately labelled parts:
  (1) exact production ledger for W(2^20) and the declared score;
  (2) production number-theory checks (r-1 factorisation, generator, GLV eigenvalue
      and short-lattice decomposition bound);
  (3) an end-to-end Cheon attack on a small elliptic-curve group and on a small
      prime-order group, at their actual toy scope.
Toy runs are evidence for the algebra and the implementation only; they are not a
production cost measurement or a production success estimate.
"""
import json, math, random, sys
from decimal import Decimal, getcontext, ROUND_CEILING
from fractions import Fraction

# ---------- number theory helpers ----------
def is_prime(n):
    if n < 2:
        return False
    small = (2,3,5,7,11,13,17,19,23,29,31,37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def factor_trial(n, bound=10**7):
    f = {}
    d = 2
    while d * d <= n and d <= bound:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def element_order_check(g, r, primes):
    return all(pow(g, (r - 1) // q, r) != 1 for q in primes)

# ---------- fixed-base comb (used by the toy end-to-end searches) ----------
def comb_precompute(P, bits, w, add, neg, zero):
    """T[i][d] = [d*2^(i*w)]P for i < ceil(bits/w), d < 2^w."""
    digits = (bits + w - 1) // w
    blocks = []
    Pi = P
    for _ in range(digits):
        T = [zero, Pi]
        for _ in range(2, 1 << w):
            T.append(add(T[-1], Pi))
        blocks.append(T)
        for _ in range(w):
            Pi = add(Pi, Pi)
    return (blocks, digits, w, neg, zero)

def comb_eval(k, tables, add):
    blocks, digits, w, neg, zero = tables
    flip = k < 0
    if flip:
        k = -k
    R = zero
    for i in range(digits):
        d = (k >> (i * w)) & ((1 << w) - 1)
        if d:
            R = add(R, blocks[i][d])
    if flip and R is not zero and R is not None:
        R = neg(R)
    return R

# ---------- toy elliptic curve y^2 = x^3 + 4 ----------
class ToyEC:
    def __init__(self, p, a=0, b=4):
        self.p, self.a, self.b = p, a, b
    def add(self, P, Q):
        p = self.p
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0: return None
        if P == Q:
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)
    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)
    def mul(self, k, P):
        R = None
        if k < 0:
            k = -k; P = self.neg(P)
        while k:
            if k & 1: R = self.add(R, P)
            P = self.add(P, P); k >>= 1
        return R

def cheon_toy_ec(seed, p=100183, r=4759, P=(43374, 99049), d=2, w=4):
    ec = ToyEC(p)
    rng = random.Random(seed)
    primes = list(factor_trial(r - 1))
    zeta = next(g for g in range(2, r) if element_order_check(g, r, primes))
    c = pow(zeta, d, r)
    M = (r - 1) // d
    assert pow(c, M, r) == 1 and all(pow(c, M // q, r) != 1 for q in factor_trial(M))
    zcheck = pow(zeta, M, r)
    assert pow(zcheck, d, r) == 1
    m = math.isqrt(M)
    if m * m < M: m += 1
    mh = (M + m - 1) // m
    bits = r.bit_length()
    results = []
    for _ in range(6):
        s = rng.randrange(1, r)
        T = ec.mul(s, P)
        Ad = ec.mul(pow(s, d, r), P)
        # fixed-base tables for the two phase-1 bases
        tab_Ad = comb_precompute(Ad, bits, w, ec.add, ec.neg, None)
        tab_P = comb_precompute(P, bits, w, ec.add, ec.neg, None)
        baby = {}
        cu = 1
        cinv = pow(c, -1, r)
        for u in range(m):
            X = comb_eval(cu, tab_Ad, ec.add)
            baby.setdefault(X, u)
            cu = cu * cinv % r
        ccur = 1; cm = pow(c, m, r); found = None
        for v in range(mh + 1):
            Y = comb_eval(ccur, tab_P, ec.add)
            if Y in baby:
                found = baby[Y] + m * v
                break
            ccur = ccur * cm % r
        assert found is not None and pow(c, found, r) == pow(s, d, r)
        k0 = found
        # phase 2: S is a runtime base, so its table is built here
        S = ec.mul(pow(zeta, -k0, r), T)
        tab_S = comb_precompute(S, bits, w, ec.add, ec.neg, None)
        m2 = math.isqrt(d)
        if m2 * m2 < d: m2 += 1
        m2h = (d + m2 - 1) // m2
        zinv = pow(zcheck, -1, r)
        baby2 = {}
        zu = 1
        for u in range(m2):
            X = comb_eval(zu, tab_S, ec.add)
            baby2.setdefault(X, u)
            zu = zu * zinv % r
        zcur = 1; zm = pow(zcheck, m2, r); k1 = None
        for v in range(m2h + 1):
            Y = comb_eval(zcur, tab_P, ec.add)
            if Y in baby2:
                k1 = baby2[Y] + m2 * v
                break
            zcur = zcur * zm % r
        assert k1 is not None
        rec = pow(zeta, k0, r) * pow(zcheck, k1, r) % r
        results.append({
            "s_verified": ec.mul(rec, P) == T,
            "s_matches": rec == s,
            "k0_ok": pow(c, k0, r) == pow(s, d, r),
            "phase2_runtime_table": True,
        })
    return {
        "scope": "toy 29-bit curve y^2=x^3+4 over F_100183, order 4759, d=2, comb w=%d" % w,
        "trials": len(results),
        "all_verified": all(x["s_verified"] and x["s_matches"] and x["k0_ok"] for x in results),
        "details": results,
    }

def comb_selftest_ec(seed, p=100183, r=4759, P=(43374, 99049), w=4, trials=32):
    ec = ToyEC(p); rng = random.Random(seed)
    bits = r.bit_length()
    tab = comb_precompute(P, bits, w, ec.add, ec.neg, None)
    ok = 0
    for _ in range(trials):
        k = rng.randrange(1, r)
        ok += int(comb_eval(k, tab, ec.add) == ec.mul(k, P))
    k = rng.randrange(1, r); Q = ec.mul(k, P)
    tabQ = comb_precompute(Q, bits, w, ec.add, ec.neg, None)
    ok2 = 0
    for _ in range(trials):
        k2 = rng.randrange(1, r)
        ok2 += int(comb_eval(k2, tabQ, ec.add) == ec.mul(k2, Q))
    return {"trials": trials, "fixed_base_ok": ok, "runtime_base_ok": ok2,
            "all_verified": ok == trials and ok2 == trials}

# ---------- toy prime-order group (additive model of the exponent relation) ----------
def cheon_toy_group(seed, r=1032193, d=1024, w=4):
    k = 1
    while True:
        q = k * r + 1
        if is_prime(q): break
        k += 1
    g = 2
    P = pow(g, (q - 1) // r, q)
    assert pow(P, r, q) == 1 and P != 1
    add = lambda A, B: (A * B) % q
    neg = lambda A: pow(A, -1, q)
    zero = 1
    rng = random.Random(seed)
    primes = list(factor_trial(r - 1))
    zeta = next(g for g in range(2, r) if element_order_check(g, r, primes))
    c = pow(zeta, d, r); M = (r - 1) // d
    zcheck = pow(zeta, M, r)
    m = math.isqrt(M)
    if m * m < M: m += 1
    mh = (M + m - 1) // m
    bits = r.bit_length()
    ok = 0
    for _ in range(4):
        s = rng.randrange(1, r)
        Tp = pow(P, s, q); Ad = pow(P, pow(s, d, r), q)
        tab_Ad = comb_precompute(Ad, bits, w, add, neg, zero)
        tab_P = comb_precompute(P, bits, w, add, neg, zero)
        baby = {}; cu = 1; cinv = pow(c, -1, r)
        for u in range(m):
            baby.setdefault(comb_eval(cu, tab_Ad, add), u)
            cu = cu * cinv % r
        ccur = 1; cm = pow(c, m, r); k0 = None
        for v in range(mh + 1):
            Y = comb_eval(ccur, tab_P, add)
            if Y in baby:
                k0 = baby[Y] + m * v; break
            ccur = ccur * cm % r
        assert k0 is not None and pow(c, k0, r) == pow(s, d, r)
        S = pow(Tp, pow(zeta, -k0, r), q)
        tab_S = comb_precompute(S, bits, w, add, neg, zero)
        zinv = pow(zcheck, -1, r); m2 = math.isqrt(d)
        if m2 * m2 < d: m2 += 1
        m2h = (d + m2 - 1) // m2
        baby2 = {}; zu = 1
        for u in range(m2):
            baby2.setdefault(comb_eval(zu, tab_S, add), u)
            zu = zu * zinv % r
        zcur = 1; zm = pow(zcheck, m2, r); k1 = None
        for v in range(m2h + 1):
            Y = comb_eval(zcur, tab_P, add)
            if Y in baby2:
                k1 = baby2[Y] + m2 * v; break
            zcur = zcur * zm % r
        assert k1 is not None
        rec = pow(zeta, k0, r) * pow(zcheck, k1, r) % r
        ok += int(rec == s and pow(P, rec, q) == Tp)
    return {
        "scope": "toy prime-order group of order %d, d=%d, q=%d, comb w=%d" % (r, d, q, w),
        "trials": 4, "recovered": ok, "all_verified": ok == 4,
        "phase2_runtime_table": True,
    }

# ---------- production ledger ----------
def production(seed):
    r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
    p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
    d = 2 ** 20
    assert (r - 1) % d == 0
    M = (r - 1) // d
    m = math.isqrt(M)
    if m * m < M: m += 1
    mh = (M + m - 1) // m
    # factorisation of r-1 (public, verified by product and primality)
    fac = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1,
           906349: 2, 2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}
    prod = 1
    for q, e in fac.items():
        prod *= q ** e
    assert prod == r - 1
    assert all(is_prime(q) for q in fac)
    zeta = 7
    assert element_order_check(zeta, r, list(fac))
    c = pow(zeta, d, r)
    assert pow(c, M, r) == 1
    assert all(pow(c, M // q, r) != 1 for q in fac)
    # endomorphism: beta^3=1 mod p, lambda^2+lambda+1=0 mod r
    s3 = pow((-3) % p, (p + 1) // 4, p)
    assert s3 * s3 % p == (-3) % p
    beta = ((-1 + s3) * pow(2, -1, p)) % p
    assert pow(beta, 3, p) == 1 and beta != 1 and (beta * beta + beta + 1) % p == 0
    lam = pow(7, (r - 1) // 3, r)
    assert (lam * lam + lam + 1) % r == 0 and lam != 1
    # production curve: y^2=x^3+4, generator P0; verify [r]P0=O and the matching
    # nontrivial cube root beta with phi_beta(P0)=[lam]P0.
    P0 = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
          0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
    def ec_add(A, B):
        if A is None: return B
        if B is None: return A
        x1, y1 = A; x2, y2 = B
        if x1 == x2 and (y1 + y2) % p == 0: return None
        if A == B:
            lamd = (3 * x1 * x1) * pow(2 * y1, -1, p) % p
        else:
            lamd = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lamd * lamd - x1 - x2) % p
        return (x3, (lamd * (x1 - x3) - y1) % p)
    def ec_mul(k, A):
        R = None
        while k:
            if k & 1: R = ec_add(R, A)
            A = ec_add(A, A); k >>= 1
        return R
    assert ec_mul(r, P0) is None
    beta_match = None
    for btest in (beta, beta * beta % p):
        if (btest * P0[0] % p, P0[1]) == ec_mul(lam, P0):
            beta_match = btest
            break
    assert beta_match is not None
    # short lattice basis for GLV, Gauss-reduced
    def rdiv(a, b):
        # nearest integer to a/b, b>0, exact for integers
        return (a + b // 2) // b if a >= 0 else -((-a + b // 2) // b)

    def gauss(b1, b2):
        while True:
            n1 = b1[0]*b1[0] + b1[1]*b1[1]
            n2 = b2[0]*b2[0] + b2[1]*b2[1]
            if n1 > n2: b1, b2 = b2, b1; n1, n2 = n2, n1
            mu = rdiv(b1[0]*b2[0] + b1[1]*b2[1], n1)
            if mu == 0: return b1, b2
            b2 = (b2[0] - mu*b1[0], b2[1] - mu*b1[1])
    b1, b2 = gauss((r, 0), (-lam, 1))
    assert (b1[0] + lam*b1[1]) % r == 0 and (b2[0] + lam*b2[1]) % r == 0
    det = b1[0]*b2[1] - b1[1]*b2[0]
    rng = random.Random(seed)
    max_bits = 0
    for _ in range(64):
        k = rng.randrange(1, r)
        best = None
        for e1 in (-1, 0, 1):
            for e2 in (-1, 0, 1):
                c1 = rdiv(k*b2[1], det) + e1
                c2 = rdiv(-k*b1[1], det) + e2
                a = c1*b1[0] + c2*b2[0]
                bb = c1*b1[1] + c2*b2[1]
                k1 = k - a; k2 = -bb
                assert (k1 + lam*k2 - k) % r == 0
                mb = max(abs(k1).bit_length(), abs(k2).bit_length())
                if best is None or mb < best: best = mb
        max_bits = max(max_bits, best)
    assert max_bits <= 130
    # Fixed-base comb/window scalar multiplication.  Each 127-bit GLV half is
    # written in 8 base-2^16 digits; a precomputed table per base gives one
    # addition per digit, so the online phase contains only reference additions.
    W_BITS = 16
    DIGITS = 8                             # ceil(127/16) per half
    ADDITIONS = 2 * DIGITS + 1             # two halves plus final combination
    CHARGE_ADDITIONS = 33                  # worst-case charge (> ADDITIONS)
    ADD_UNITS = Fraction(1)                # protected batched incomplete affine addition
    TABLE_LOOKUPS = 33                     # one per digit plus combination
    LOOKUP_UNITS = Fraction(2)             # 112-byte random load, address, compare
    OVERHEAD = 16                          # itemised GLV/digits/sign/scalar/loop allowance
    E = CHARGE_ADDITIONS * ADD_UNITS + TABLE_LOOKUPS * LOOKUP_UNITS + OVERHEAD
    # One-time comb precomputation: 4 bases, 8 blocks, 2^16 entries per block.
    # Fully charged, with a complete instruction allowance:
    COMB_TABLES = 4
    ENTRIES_PER_TABLE = DIGITS * (1 << W_BITS)
    PRECOMB = Fraction(2 ** 30)            # phase-1 comb tables (A_d, phi(A_d), P0, phi(P0))
    PHASE2_PRECOMB = Fraction(2 ** 29)     # runtime comb tables for S and phi(S), charged
    # Merge sort + binary search: deterministic, no hashing assumptions.
    L = (m - 1).bit_length()
    SORT = Fraction(5)
    SEARCH = Fraction(4)
    PHASE2 = 2048 * E + 8 + Fraction(1024 * 12 * 48, 2000)
    FIXED = Fraction(2 ** 40, 2000)
    W = (m + mh) * E + m * SORT + mh * SEARCH + PHASE2 + PRECOMB + PHASE2_PRECOMB + FIXED
    # rounded-up log2 at 6 decimals, high precision (Decimal.ln, 60 digits)
    getcontext().prec = 60
    lnW = Decimal(W.numerator).ln() - Decimal(W.denominator).ln()
    lg = lnW / Decimal(2).ln()
    q = int((lg * Decimal(10**6)).to_integral_value(rounding=ROUND_CEILING))
    score = Fraction(q, 10**6)
    assert Decimal(q) / Decimal(10**6) >= lg
    assert Decimal(q - 1) / Decimal(10**6) < lg
    assert score == Fraction(31332353, 250000), score
    # two merge-sort buffers of m records (96-byte point + 16-byte index),
    # four comb tables, a B=2^64-chain batch workspace, fixed reserve.
    BATCH_LOG2 = 64
    record_bytes = 96 + 16
    comb_bytes = COMB_TABLES * ENTRIES_PER_TABLE * record_bytes
    peak_bytes = 2 * m * record_bytes + comb_bytes + 2 ** BATCH_LOG2 * 240 + 2 ** 22
    pre_units = m * (E + SORT) + PRECOMB + PHASE2_PRECOMB
    return {
        "r_bits": r.bit_length(), "d": d, "M_bits": M.bit_length(),
        "m": str(m), "m_hat": str(mh),
        "log2_m_plus_mhat": math.log2(m + mh),
        "addition_units": str(ADD_UNITS), "additions_per_scalar": ADDITIONS,
        "charged_additions": CHARGE_ADDITIONS, "table_lookups": TABLE_LOOKUPS,
        "lookup_units": str(LOOKUP_UNITS), "E_units": str(E),
        "comb_window_bits": W_BITS, "comb_digits_per_half": DIGITS,
        "comb_tables": COMB_TABLES, "comb_entries_per_table": ENTRIES_PER_TABLE,
        "precomb_units": str(PRECOMB), "phase2_precomb_units": str(PHASE2_PRECOMB), "sort_units_per_entry": str(SORT),
        "search_units_per_lookup": str(SEARCH), "sort_passes": L,
        "record_bytes": record_bytes, "batch_size_log2": BATCH_LOG2,
        "W_log2_high_precision": str(lg),
        "score": "%d/%d" % (score.numerator, score.denominator),
        "score_decimal": float(score),
        "success_probability": 1,
        "memory_log2_bytes": math.log2(peak_bytes),
        "preprocessing_log2_units": math.log2(pre_units),
        "r_minus_1_factorisation_verified": True,
        "zeta_generator": zeta, "zeta_d_order_M": True,
        "beta_cube_root_verified": True, "lambda_cubic_root_verified": True,
        "production_glv_pair_verified": True, "production_rP0_identity": True,
        "lambda_bits": lam.bit_length(),
        "glv_basis": [list(b1), list(b2)],
        "glv_max_half_bits": max_bits,
        "phase1_steps": str(m + mh), "phase2_scalar_mults": 2048,
    }

def main():
    cfg = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    seed = int(cfg.get("seed", 20261024))
    out = {
        "comb_selftest_ec": comb_selftest_ec(seed),
        "toy_ec_cheon": cheon_toy_ec(seed),
        "toy_group_cheon": cheon_toy_group(seed + 1),
        "production_ledger": production(seed),
    }
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()

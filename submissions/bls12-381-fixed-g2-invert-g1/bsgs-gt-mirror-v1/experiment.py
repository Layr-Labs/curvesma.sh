"""Finite-scope checks for the bls12-381-fixed-g2-invert-g1 plain-BSGS package.

Evidence only at the stated scope. Nothing here executes production attack
work, allocates the exponential table, times anything, or proves the
production cost bound; that bound is analytic (proof.md).

Modes:
  additive - plain-BSGS coverage/reconstruction logic on Z/nZ.
  group    - the full algorithm on small prime-order multiplicative groups
             with canonical keys: sort, giant stepping, reconstruction and
             the D^(a') = z verification.
  kernel   - the exact 6x64-bit CIOS Montgomery multiplication of proof.md over
             the real BLS12-381 base field, checked against Python integers,
             with the per-instruction tally; the Fermat chain for p-2; and the
             protected generator P0 on the curve with [r]P0 = infinity; then
             the Fp12 schoolbook multiply/reduce kernel, Fp2 arithmetic and
             the G2 generator Q0 on the twist with [r]Q0 = infinity.
  ledger   - exact integer checks of full coverage vs the real subgroup order r
             and of the work ledger against 20000*2^143.5.
"""
import json
import sys

BLS_P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
            "1eabfffeb153ffffb9feffffffffaaab", 16)
BLS_R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BLS_B = 4
P0X = int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac58"
          "6c55e83ff97a1aeffb3af00adb22c6bb", 16)
P0Y = int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3ed"
          "d03cc744a2888ae40caa232946c5e7e1", 16)
MASK = (1 << 64) - 1
LIMBS = 6


def plan(n, m):
    h = (n + m - 1) // m
    return h


def additive(n, m):
    # Plain BSGS on Z/nZ: table {0..m-1}, stride m, covered = m*h >= n.
    h = plan(n, m)
    table = {j: j for j in range(m)}
    covered = min(m * h, n)
    recovered = 0
    for x in range(n):
        out = None
        for i in range(h):
            t = (x - i * m) % n
            if t in table:
                out = (i * m + table[t]) % n
                break
        assert out == x, ("scalar not recovered", n, m, x)
        recovered += 1
    return {"order": n, "m": m, "h": h, "covered": covered,
            "recovered": recovered, "exact_coverage": recovered == n == covered}


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def prime_modulus(q, k):
    # smallest prime p = k*q + 1 (k even trial), so Fp* has a subgroup of order q
    while True:
        p = k * q + 1
        if is_prime(p):
            return p, k
        k += 2


def subgroup_gen(p, q):
    # generator of the order-q subgroup of Fp*
    k = (p - 1) // q
    for g in range(2, p):
        D = pow(g, k, p)
        if D != 1 and pow(D, q, p) == 1:
            return D
    raise AssertionError("no subgroup generator")


def group(q, m):
    # Full plain-BSGS algorithm over the order-q subgroup of Fp*, canonical
    # integer keys: baby table, sort, giant stepping, reconstruction a'=im+j
    # and the D^(a') = z verification, for every challenge in the subgroup.
    assert is_prime(q)
    p, k = prime_modulus(q, 2)
    D = subgroup_gen(p, q)
    h = (q + m - 1) // m
    E = pow(D, q - (m % q), p)  # D^(-m)
    table = {}
    R = 1
    for j in range(m):
        key = R  # canonical: Fp element, distinct while j < q
        assert key not in table
        table[key] = j
        R = R * D % p
    recs = sorted(table.items())  # deterministic sort by canonical key
    keys = [k for k, _ in recs]
    assert len(keys) == m
    import bisect
    recovered = 0
    for a in range(q):
        z = pow(D, a, p)
        assert z != 0 and pow(z, q, p) == 1  # subgroup check
        gamma = z
        out = None
        for _ in range(h):
            t = bisect.bisect_left(keys, gamma)
            if t < m and keys[t] == gamma:
                a_ = (_ * m + recs[t][1]) % q
                if pow(D, a_, p) == z:
                    out = a_
                break
            gamma = gamma * E % p
        assert out is not None, ("challenge not recovered", q, m, a)
        assert pow(D, out, p) == z
        recovered += 1
    return {"subgroup_order": q, "field": p, "cofactor": k, "m": m, "h": h,
            "recovered": recovered, "exact_coverage": recovered == q}


def groups():
    runs = []
    for q in (101, 257, 521, 1009):
        for m in (8, 16, 32):
            runs.append(group(q, m))
    return runs


class Tally:
    def __init__(self):
        self.n = 0

    def __call__(self, k=1):
        self.n += k


def mont_mul(a, b, p, n0, t_):
    """6x64-bit CIOS Montgomery multiplication with the proof.md tally."""
    L = LIMBS
    A = [(a >> (64 * i)) & MASK for i in range(L)]
    B = [(b >> (64 * i)) & MASK for i in range(L)]
    P = [(p >> (64 * i)) & MASK for i in range(L)]
    t_(2 * L)  # LDR a, b limbs
    t_(L)      # materialize p limbs
    T = [0] * (L + 2)
    for i in range(L):
        C = 0
        for j in range(L):
            prod = A[j] * B[i]
            t_(2)                     # MUL, UMULH
            s = T[j] + (prod & MASK) + C
            t_(3)                     # ADDS lo, ADC hi, ADDS carry-in
            T[j] = s & MASK
            C = (prod >> 64) + (s >> 64)
            t_(1)                     # ADC carry-out
        s = T[L] + C
        T[L] = s & MASK
        T[L + 1] = s >> 64
        t_(2)
        mq = (T[0] * n0) & MASK
        t_(1)                         # MUL
        C = 0
        for j in range(L):
            prod = mq * P[j]
            t_(2)
            s = T[j] + (prod & MASK) + C
            t_(3)
            T[j] = s & MASK
            C = (prod >> 64) + (s >> 64)
            t_(1)
        s = T[L] + C
        T[L] = s & MASK
        T[L + 1] += s >> 64
        t_(2)
        assert T[0] == 0
        T = T[1:] + [0]
        t_(L + 1)                     # limb shift moves
    val = sum(T[i] << (64 * i) for i in range(L + 1))
    t_(L + 1)                         # SUBS/SBCS chain incl. top
    if val >= p:
        val -= p
    t_(L)                             # CSEL
    t_(L)                             # STR
    t_(24)                            # call/ret, prologue/epilogue
    return val


def kernel(trials=300):
    p = BLS_P
    n0 = (-pow(p, -1, 1 << 64)) % (1 << 64)
    Rm = 1 << (64 * LIMBS)
    Rinv = pow(Rm, -1, p)
    seed = 0x5eed
    worst = 0
    ok = True
    for k in range(trials):
        seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        a = pow(seed, 7, p)
        b = pow(seed + 1, 11, p)
        if k == 0:
            a, b = p - 1, p - 1
        t = Tally()
        got = mont_mul(a, b, p, n0, t)
        ok &= (got == a * b * Rinv % p) and got < p
        worst = max(worst, t.n)
    e = p - 2
    t = Tally()
    x = (P0X * Rm) % p
    acc = x
    for bit in bin(e)[3:]:
        acc = mont_mul(acc, acc, p, n0, t)
        if bit == "1":
            acc = mont_mul(acc, x, p, n0, t)
    inv_ok = (acc * Rinv % p) * P0X % p == 1
    on_curve = (P0Y ** 2 - P0X ** 3 - BLS_B) % p == 0

    def g1_add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
            return None
        if P == Q:
            l = (3 * P[0] * P[0]) * pow(2 * P[1], p - 2, p) % p
        else:
            l = (Q[1] - P[1]) * pow(Q[0] - P[0], p - 2, p) % p
        x3 = (l * l - P[0] - Q[0]) % p
        return (x3, (l * (P[0] - x3) - P[1]) % p)

    def g1_mul(k, P):
        R = None
        while k:
            if k & 1:
                R = g1_add(R, P)
            P = g1_add(P, P)
            k >>= 1
        return R
    order_ok = g1_mul(BLS_R, (P0X, P0Y)) is None

    # Fp12 = Fp[w]/(w^12 - 2w^6 + 2): schoolbook multiply + reduction.
    def fp12_mul(A, Bc):
        assert len(A) == 12 and len(Bc) == 12
        prod = [0] * 23
        nmuls = 0
        for i in range(12):
            for j in range(12):
                prod[i + j] = (prod[i + j] + A[i] * Bc[j]) % p
                nmuls += 1
        for k in range(22, 11, -1):
            c = prod[k]
            if c:
                prod[k - 12] = (prod[k - 12] - 2 * c) % p
                prod[k - 6] = (prod[k - 6] + 2 * c) % p
        return prod[:12], nmuls
    f12_ok = True
    seed = 0x1f12
    for _ in range(200):
        seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        A = [pow(seed + i, 5, p) for i in range(12)]
        Bc = [pow(seed + 100 + i, 7, p) for i in range(12)]
        got, nm = fp12_mul(A, Bc)
        # independent check: naive polynomial product reduced mod w^12-2w^6+2
        naive = [0] * 23
        for i in range(12):
            for j in range(12):
                naive[i + j] = (naive[i + j] + A[i] * Bc[j]) % p
        for k in range(22, 11, -1):
            c = naive[k] % p
            naive[k - 12] = (naive[k - 12] - 2 * c) % p
            naive[k - 6] = (naive[k - 6] + 2 * c) % p
        f12_ok &= (got == [x % p for x in naive[:12]]) and nm == 144

    # Fp2 = Fp[u]/(u^2+1) and the G2 twist y^2 = x^3 + 4(1+u).
    def fp2_add(A, Bc):
        return ((A[0] + Bc[0]) % p, (A[1] + Bc[1]) % p)

    def fp2_mul(A, Bc):
        return ((A[0] * Bc[0] - A[1] * Bc[1]) % p,
                (A[0] * Bc[1] + A[1] * Bc[0]) % p)

    def fp2_inv(A):
        n = (A[0] * A[0] + A[1] * A[1]) % p
        ni = pow(n, p - 2, p)
        return (A[0] * ni % p, (-A[1]) * ni % p)
    fp2_ok = True
    for t in range(50):
        A = (pow(3 * t + 1, 5, p), pow(5 * t + 2, 5, p))
        Bc = (pow(7 * t + 3, 5, p), pow(11 * t + 4, 5, p))
        if A == (0, 0) or Bc == (0, 0):
            continue
        C = fp2_mul(A, Bc)
        fp2_ok &= fp2_mul(C, fp2_inv(Bc)) == A and fp2_mul(A, fp2_inv(A)) == (1, 0)
    B2 = (4 % p, 4 % p)  # twist constant 4(1+u)

    def g2_add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        negqy = ((-Q[1][0]) % p, (-Q[1][1]) % p)
        if P[0] == Q[0] and P[1] == negqy:
            return None
        if P == Q:
            xsq = fp2_mul(P[0], P[0])
            num = ((3 * xsq[0]) % p, (3 * xsq[1]) % p)
            den = ((2 * P[1][0]) % p, (2 * P[1][1]) % p)
        else:
            num = ((Q[1][0] - P[1][0]) % p, (Q[1][1] - P[1][1]) % p)
            den = ((Q[0][0] - P[0][0]) % p, (Q[0][1] - P[0][1]) % p)
        l = fp2_mul(num, fp2_inv(den))
        x3 = (l[0] * l[0] - l[1] * l[1] - P[0][0] - Q[0][0]) % p, \
             (2 * l[0] * l[1] - P[0][1] - Q[0][1]) % p
        dx = ((P[0][0] - x3[0]) % p, (P[0][1] - x3[1]) % p)
        ly = fp2_mul(l, dx)
        y3 = ((ly[0] - P[1][0]) % p, (ly[1] - P[1][1]) % p)
        return (x3, y3)

    def g2_mul(k, P):
        R = None
        while k:
            if k & 1:
                R = g2_add(R, P)
            P = g2_add(P, P)
            k >>= 1
        return R
    # Q0 exactly as pinned in the protected target profile ([c0, c1] = c0+c1*u).
    Q0X = (int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d17"
               "70bac0326a805bbefd48056c8c121bdb8", 16) % p,
           int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f50"
               "49334cf11213945d57e5ac7d055d042b7e", 16) % p)
    Q0Y = (int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d1"
               "2c923ac9cc3baca289e193548608b82801", 16) % p,
           int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99"
               "ab3f370d275cec1da1aaa9075ff05f79be", 16) % p)
    Q0 = (Q0X, Q0Y)
    lhs = fp2_mul(Q0Y, Q0Y)
    rhs = fp2_add(fp2_mul(fp2_mul(Q0X, Q0X), Q0X), B2)
    q0_curve = lhs == rhs
    q0_order = g2_mul(BLS_R, Q0) is None
    return {"field": "BLS12-381 Fp", "limbs": LIMBS, "n0_prime": hex(n0),
            "trials": trials, "mont_mul_correct": bool(ok),
            "mont_mul_worst_tally_instructions": worst,
            "mont_mul_budget_charged": 8192,
            "fermat_squarings": e.bit_length() - 1,
            "fermat_multiplies": bin(e).count("1") - 1,
            "fermat_total_mults": e.bit_length() - 1 + bin(e).count("1") - 1,
            "fermat_budget_mults": 1024, "fermat_chain_inverts": bool(inv_ok),
            "generator_on_curve": bool(on_curve),
            "generator_order_r": bool(order_ok),
            "fp12_kernel_correct_144muls": bool(f12_ok),
            "fp2_inverse_correct": bool(fp2_ok),
            "g2_generator_on_twist": bool(q0_curve),
            "g2_generator_order_r": bool(q0_order)}


def ledger():
    r = BLS_R
    m = h = 1 << 128
    cov = m * h
    R_pre = (1 << 149) + m * (1 << 29) + (1 << 48)
    R_total = R_pre + h * (1 << 28)
    fp_inv = 1024 * 8192 + (1 << 14)
    fp2_inv = 1 << 24
    g2step = fp2_inv + 4 * 65536 + 10 * 1024 + (1 << 12)
    rec = (1 << 21) + (1 << 20) + 128 * 2 * (1 << 20)
    it = 129 * (1 << 20) + (1 << 21)
    # W < 2^(287/2)  <=>  R^2 < 20000^2 * 2^287 (exact integers)
    bound_sq = (20000 ** 2) * (1 << 287)
    return {"coverage": str(cov), "coverage_above_r": cov > r,
            "full_success_probability_one": True,
            "r_bits": r.bit_length(), "r_below_2_255": r < (1 << 255),
            "R_pre": str(R_pre), "R_total": str(R_total),
            "R_total_sq_below_20000sq_2_287": R_total * R_total < bound_sq,
            "R_pre_sq_below_20000sq_2_287": R_pre * R_pre < bound_sq,
            "g2_step_bound": g2step, "g2_step_below_2_25": g2step < (1 << 25),
            "per_record": rec, "per_record_below_2_29": rec < (1 << 29),
            "per_iteration": it, "per_iteration_below_2_28": it < (1 << 28)}


def main():
    cfg = json.load(sys.stdin)
    mode = cfg.get("mode", "all")
    res = {"scope": "toy groups, one real-field kernel model and exact integer "
                    "checks; does not validate production cost bounds"}
    if mode in ("all", "additive"):
        runs = []
        for n in (101, 257, 521, 1009):
            for m in (2, 4, 8, 16):
                runs.append(additive(n, m))
        res["additive_runs"] = len(runs)
        res["additive_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "group"):
        runs = groups()
        res["group"] = runs
        res["group_runs"] = len(runs)
        res["group_all_exact"] = all(x["exact_coverage"] for x in runs)
    if mode in ("all", "kernel"):
        res["kernel"] = kernel()
    if mode in ("all", "ledger"):
        res["ledger"] = ledger()
    print(json.dumps(res))


if __name__ == "__main__":
    main()

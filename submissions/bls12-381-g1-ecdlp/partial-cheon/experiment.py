"""Evidence experiment for the partial-coverage Cheon submission
(bls12-381-g1-ecdlp, auxiliary powers, scored at n=2^20).

Part A verifies, at full production BLS12-381 parameters with pure
standard-library arithmetic: 2^20 divides r-1; the complete prime
factorization of r-1 (each listed factor is proved prime by Miller-Rabin and
the product is checked to equal r-1, which is a polynomial-time verification
of the factorization the attack needs); a generator zeta of Z_r* from that
factorization; eta = zeta^(2^20) of exact order N = (r-1)/2^20; and
zeta_e = zeta^N of exact order 2^20. No production-size discrete logarithm
is solved here.

Part B runs the complete two-phase Cheon attack end to end on a small toy
group (j=0 curve of prime order n=1047379 with d=6 | n-1, N=(n-1)/6=174563),
with the auxiliary powers A_1=[s]P0 and A_d=[s^d]P0, at both full coverage
(deterministic success 1) and the partial coverage used in the claim
(expected success m^2/N ~ 0.40). It reports measured success, wrong
extractions and the scalar-multiplication count. Toy evidence is evidence at
toy scope only.

Stdlib only; deterministic seeds; reads JSON input on stdin.
"""
import json
import math
import random
import sys

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def is_prime_mr(v):
    if v < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if v % sp == 0:
            return v == sp
    d = v - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, v)
        if x in (1, v - 1):
            continue
        for _ in range(s - 1):
            x = x * x % v
            if x == v - 1:
                break
        else:
            return False
    return True


def is_prime_trial(v):
    if v < 2:
        return False
    for d in range(2, math.isqrt(v) + 1):
        if v % d == 0:
            return False
    return True


def ec_add(Pt, Qt, mod, b):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    x1, y1 = Pt
    x2, y2 = Qt
    if x1 == x2 and (y1 + y2) % mod == 0:
        return None
    if Pt == Qt:
        m = (3 * x1 * x1) * pow(2 * y1, -1, mod) % mod
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, mod) % mod
    x3 = (m * m - x1 - x2) % mod
    return (x3, (m * (x1 - x3) - y1) % mod)


def ec_mul(k, Pt, mod, b):
    if k % mod == 0:
        return None
    acc = None
    base = Pt
    k %= mod
    while k:
        if k & 1:
            acc = ec_add(acc, base, mod, b)
        base = ec_add(base, base, mod, b)
        k >>= 1
    return acc


# ---------------------------------------------------------------------------
# Part A: production BLS12-381 parameter checks
# ---------------------------------------------------------------------------
BLS_R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D = 2 ** 20  # scoring n; d | r-1
# prime factorization of r-1 (hex): the attack re-derives it at runtime by
# cheap factoring (all prime factors are below 2^28); here it is *verified*.
R_M1_FACTORS = [
    (0x2, 32), (0x3, 1), (0xB, 1), (0x13, 1), (0x27C1, 1), (0x1EA57, 1),
    (0xD1C83, 1), (0xDD46D, 2), (0x264679, 1), (0x26987B, 1),
    (0x320238B, 1), (0xF2F5565, 2),
]


def production_checks():
    r = BLS_R
    out = {}
    assert (r - 1) % D == 0
    out["d_power_of_two"] = D
    out["d_divides_r_minus_1"] = True
    prod = 1
    for q, e in R_M1_FACTORS:
        assert is_prime_mr(q), hex(q)
        prod *= q ** e
    out["factorization_product_equals_r_minus_1"] = prod == r - 1
    out["all_prime_factors_below_2_pow_28"] = all(
        q < 2 ** 28 for q, _ in R_M1_FACTORS)
    primes = [q for q, _ in R_M1_FACTORS]

    def is_generator(g):
        if math.gcd(g, r) != 1:
            return False
        return all(pow(g, (r - 1) // q, r) != 1 for q in primes)

    zeta = next(g for g in range(2, 1000) if is_generator(g))
    out["zeta_generator_small"] = zeta
    N = (r - 1) // D
    eta = pow(zeta, D, r)
    assert pow(eta, N, r) == 1
    nfac = dict(R_M1_FACTORS)
    nfac[2] = 32 - 20
    primes_N = [q for q, e in nfac.items() if e > 0 and q != 2]
    out["eta_order_is_N"] = all(pow(eta, N // q, r) != 1 for q in primes_N)
    ze = pow(zeta, N, r)
    assert pow(ze, D, r) == 1
    out["zeta_e_order_is_d"] = pow(ze, D // 2, r) != 1
    out["N_decimal"] = N
    out["log2_N"] = math.log2(N)
    # d-th powers land in the order-N subgroup (Fermat)
    rng = random.Random(20261022)
    ok = all(pow(pow(rng.randrange(1, r), D, r), N, r) == 1 for _ in range(4))
    out["dth_powers_in_order_N_subgroup"] = ok
    # A_d index rule for the input encoding: A[i-1] encodes A_i
    out["A_d_index"] = D - 1
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale end-to-end Cheon attack
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379   # prime order; d = 6 | n-1
TOY_G = (1, 245324)
TOY_D = 6


def toy_cheon(trials, mode):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    assert is_prime_trial(n) and (n - 1) % TOY_D == 0
    assert ec_mul(n, G, p, b) is None
    N = (n - 1) // TOY_D
    fac = {}
    m = n - 1
    for q in range(2, math.isqrt(m) + 1):
        while m % q == 0:
            fac[q] = fac.get(q, 0) + 1
            m //= q
    if m > 1:
        fac[m] = fac.get(m, 0) + 1
    primes = list(fac)

    def is_generator(g):
        if math.gcd(g, n) != 1:
            return False
        return all(pow(g, (n - 1) // q, n) != 1 for q in primes)

    zeta = next(g for g in range(2, 1000) if is_generator(g))
    eta = pow(zeta, TOY_D, n)          # order N
    assert pow(eta, N, n) == 1
    ze = pow(zeta, N, n)               # order d

    m1 = math.isqrt(N) + 1
    M = (N + m1 - 1) // m1
    if mode == "full":
        m = m1
        predicted = 1.0
    else:
        m = math.isqrt(int(0.40 * N)) + 1
        predicted = (m / m1) * (m / M)

    def attack(s, seed):
        rng = random.Random(seed)
        A1 = ec_mul(s, G, p, b)
        Ad = ec_mul(pow(s, TOY_D, n), G, p, b)   # auxiliary power A_d
        if mode == "full":
            J = range(M)
        else:
            J = rng.sample(range(M), min(m, M))
        table = {}
        e = 1
        mults = 0
        for i in range(m):
            pt = ec_mul(e, G, p, b)
            mults += 1
            if pt not in table:
                table[pt] = i
            e = e * eta % n
        for j in J:
            Y = ec_mul(pow(eta, (-m1 * j) % N, n), Ad, p, b)
            mults += 1
            if Y in table:
                i = table[Y]
                k1 = i + m1 * j
                if k1 >= N:
                    continue
                Qp = ec_mul(pow(zeta, (-k1) % (n - 1), n), A1, p, b)
                mults += 1
                for t in range(TOY_D):
                    if ec_mul(pow(zeta, (t * N) % (n - 1), n), G, p, b) == Qp:
                        s_rec = pow(zeta, (k1 + t * N) % (n - 1), n)
                        mults += TOY_D // 2
                        return (s_rec == s), mults
        return None, mults

    rng = random.Random(2026)
    succ = 0
    wrong = 0
    totm = 0
    for ti in range(trials):
        s = rng.randrange(1, n)
        got, mults = attack(s, 7000 + ti)
        totm += mults
        if got is True:
            succ += 1
        elif got is False:
            wrong += 1
    return {
        "curve": "y^2 = x^3 + 29 over F_1048609, prime order 1047379",
        "d": TOY_D,
        "N": N,
        "m1": m1,
        "M": M,
        "m": m,
        "mode": mode,
        "predicted_success": predicted,
        "trials": trials,
        "successes": succ,
        "empirical_success": succ / trials,
        "wrong_extractions": wrong,
        "mean_scalar_mults": totm / trials,
    }


def main():
    data = json.load(sys.stdin)
    out = {
        "production_checks": production_checks(),
        "toy_cheon": toy_cheon(data["trials"], data["mode"]),
        "scope": (
            "Part A verifies the Cheon preconditions and subgroup structure "
            "at full BLS12-381 parameters (d=2^20 | r-1, generator, eta of "
            "order N=(r-1)/2^20, zeta_e of order 2^20). Part B demonstrates "
            "the complete two-phase attack only on a 20-bit toy group with "
            "d=6; it is not evidence about production costs."
        ),
    }
    print(json.dumps(out))


main()

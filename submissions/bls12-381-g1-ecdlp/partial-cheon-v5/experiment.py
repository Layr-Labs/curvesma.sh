"""Evidence experiment for the partial-cheon-v5 submission (BLS12-381 G1 ECDLP
with auxiliary powers, at n = 2^20).

v5 replaces the v4 hash table with a SORTED ARRAY keyed by (x, y mod 2) -
a key that uniquely determines an affine point on an odd-prime short
Weierstrass curve (the two points sharing an x-coordinate have opposite
y-parities) - merged by an explicit bottom-up merge sort with an auxiliary
buffer, and searched by an explicit binary search. The toy attack below
runs exactly that machinery (merge sort + binary search), not a
dictionary, so the measured zero-wrong-extraction record exercises the
claimed implementation.

Part A (production, exact): verifies d = 2^20 | r-1, the complete prime
factorization of r-1 (Miller-Rabin + product equality, all factors
< 2^28), zeta = 7 generates Z_r*, eta has exact order N = (r-1)/2^20,
zeta_e has exact order 2^20, d-th powers lie in the order-N subgroup, the
parity-key uniqueness property on a random sample, and the per-offset
minimum coverage lower bound (m-1)*m/N >= 39/100 in exact integers; it
also cross-checks the closed-form coverage bound against an exhaustive
every-offset brute force on five smaller parameter sets.

Part B (toy, complete two-phase attack): d = 6 | n-1 on a 20-bit prime
order toy curve, production-shaped contiguous-window sequential giants,
parity-keyed sorted array with the explicit merge sort and binary search,
wrap-inclusive exact coverage expectation, and brute-force verification
of every extraction.

Stdlib only; deterministic (fixed seeds); single JSON object on stdout.
"""
import json
import math
import random
import sys
from fractions import Fraction


def is_prime_mr(v):
    if v < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
            if v % sp == 0:
                return v == sp
    d = v - 1
    a = 0
    while d % 2 == 0:
        d //= 2
        a += 1
    for w in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(w, d, v)
        if x in (1, v - 1):
            continue
        ok = False
        for _ in range(a):
            x = x * x % v
            if x == v - 1:
                ok = True
                break
            if x == 1:
                return False
        if not ok:
            return False
    return True


def is_prime_trial(v):
    if v < 2:
        return False
    for q in range(2, math.isqrt(v) + 1):
        if v % q == 0:
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
    add = Pt
    while k:
        if k & 1:
            acc = ec_add(acc, add, mod, b)
        add = ec_add(add, add, mod, b)
        k >>= 1
    return acc


# ---------------------------------------------------------------------------
# Explicit sorted-array machinery: bottom-up merge sort with an auxiliary
# buffer, and binary search. These are the exact primitives charged in the
# v5 ledger (one key compare and one record move per record per pass;
# ceil(log2(2m)) probes per lookup).
# ---------------------------------------------------------------------------
def merge_sort(arr):
    n = len(arr)
    buf = [None] * n
    src = arr
    dst = buf
    width = 1
    compares = 0
    moves = 0
    while width < n:
        for lo in range(0, n, 2 * width):
            mid = min(lo + width, n)
            hi = min(lo + 2 * width, n)
            i = lo
            j = mid
            k = lo
            while i < mid and j < hi:
                if src[i][0] < src[j][0]:
                    dst[k] = src[i]
                    i += 1
                else:
                    dst[k] = src[j]
                    j += 1
                compares += 1
                moves += 1
                k += 1
            while i < mid:
                dst[k] = src[i]
                i += 1
                moves += 1
                k += 1
            while j < hi:
                dst[k] = src[j]
                j += 1
                moves += 1
                k += 1
        src, dst = dst, src
        width *= 2
    if src is not arr:
        arr[:] = src
    assert all(arr[i][0] < arr[i + 1][0] for i in range(n - 1))
    return compares, moves


def binary_search(arr, key):
    """Return the index of `key` in the strictly increasing key array, or -1.

    Probes exactly ceil(log2(2*len(arr))) slots in the worst case.
    """
    n = len(arr)
    lo = 0
    hi = n - 1
    probes = 0
    budget = math.ceil(math.log2(2 * n))
    while lo <= hi:
        mid = (lo + hi) // 2
        probes += 1
        if probes > budget:
            raise AssertionError("binary search exceeded probe budget")
        if arr[mid][0] == key:
            return mid
        if arr[mid][0] < key:
            lo = mid + 1
        else:
            hi = mid - 1
    return -1


# ---------------------------------------------------------------------------
# Part A: production BLS12-381 parameter checks
# ---------------------------------------------------------------------------
BLS_R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D = 2 ** 20
R_M1_FACTORS = [
    (0x2, 32), (0x3, 1), (0xB, 1), (0x13, 1), (0x27C1, 1), (0x1EA57, 1),
    (0xD1C83, 1), (0xDD46D, 2), (0x264679, 1), (0x26987B, 1),
    (0x320238B, 1), (0xF2F5565, 2),
]


def coverage_exact(N, m1, M, m, brute=False):
    """Exact coverage analysis of the contiguous-window partial phase 1.

    A window of m consecutive blocks at offset o covers k1 iff some baby
    index i < m and block j in the window satisfy k1 = (i + m1*j) mod N -
    the wrap-around hits near the end of the k1 range are included.
    per_offset_min is the exact minimum over offsets of the covered count
    (brute-forced when M is small); the closed-form return uses the
    disjoint-arc bound (m-1)*m, valid for every offset.
    """
    if brute and M <= 4000:
        covs = []
        for o in range(M):
            cov = set()
            for t in range(m):
                j = o + t
                if j >= M:
                    j -= M
                start = (m1 * j) % N
                for i in range(m):
                    cov.add((start + i) % N)
            covs.append(len(cov))
        per_min = Fraction(min(covs), N)
        expct = Fraction(sum(covs), M * N)
        return per_min, expct
    # closed-form: each block's arc is m k1-values {(m1*j + i) mod N :
    # i < m} (start spacing m1 > m keeps consecutive arcs disjoint);
    # dropping any one potentially-overlapping block pair leaves (m-1)*m
    # covered values for every offset.
    per_min = Fraction((m - 1) * m, N)
    cov = (m - 1) * m
    return per_min, Fraction(cov, N)


def coverage_crosscheck():
    """Exhaustively validate the closed-form (m-1)*m bound on every offset
    of five smaller parameter sets (same construction as production:
    m1 = isqrt(N)+1, m = isqrt(floor(2N/5))+1)."""
    rng = random.Random(20261023)
    checked = []
    for _ in range(5):
        N = rng.randrange(500, 9000)
        m1 = math.isqrt(N) + 1
        M = -(-N // m1)
        m = math.isqrt((2 * N) // 5) + 1
        per_min, _ = coverage_exact(N, m1, M, m, brute=True)
        assert per_min >= Fraction((m - 1) * m, N), (N, m1, m, per_min)
        checked.append({"N": N, "m1": m1, "M": M, "m": m,
                        "brute_min_over_all_offsets": str(per_min),
                        "closed_form": f"{(m - 1) * m}/{N}"})
    return checked


def parity_key_check():
    """The (x, y mod 2) key uniquely determines an affine point on an odd
    prime short Weierstrass curve: the two points sharing an x have
    parities y and p - y, which differ for odd p. Verified on a random
    sample of toy-curve points: every sampled point and its negation are
    sampled as a pair (they share x by construction), each pair yields two
    distinct keys, and every key collision among the samples comes from the
    same pair of group elements (an x-coordinate can carry at most the two
    points P and -P)."""
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    rng = random.Random(20261024)
    keys = {}
    pairs_ok = 0
    collisions_same_point = 0
    for _ in range(1000):
        k = rng.randrange(1, n)
        P = ec_mul(k, G, p, b)
        nP = ec_mul(n - k, G, p, b)
        assert P is not None and nP is not None
        assert P[0] == nP[0]
        k1, k2 = (P[0], P[1] % 2), (nP[0], nP[1] % 2)
        assert k1 != k2
        for key, pt in ((k1, P), (k2, nP)):
            if key in keys:
                assert keys[key] == pt
                collisions_same_point += 1
            keys[key] = pt
        pairs_ok += 1
    return {"sampled_points": 2 * pairs_ok,
            "distinct_keys": len(keys),
            "key_collisions_all_same_point": True,
            "collision_count": collisions_same_point}

def production_checks():
    r = BLS_R
    out = {}
    assert (r - 1) % D == 0
    out["d_power_of_two"] = D
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
    m1 = math.isqrt(N) + 1
    M = -(-(N) // m1)
    m = math.isqrt((2 * N) // 5) + 1
    per_min, expct = coverage_exact(N, m1, M, m, brute=False)
    out["m1"] = m1
    out["M"] = M
    out["m"] = m
    out["per_offset_min_rational"] = f"{per_min.numerator}/{per_min.denominator}"
    out["per_offset_min_float"] = float(per_min)
    out["expectation_lower_bound_rational"] = f"{expct.numerator}/{expct.denominator}"
    out["expectation_lower_bound_float"] = float(expct)
    assert per_min >= Fraction(39, 100)
    out["coverage_at_least_0_39_every_offset"] = True
    out["coverage_closed_form_crosscheck"] = coverage_crosscheck()
    out["parity_key_uniqueness_sample"] = parity_key_check()
    rng = random.Random(20261022)
    out["dth_powers_in_order_N_subgroup"] = all(
        pow(pow(rng.randrange(1, r), D, r), N, r) == 1 for _ in range(4))
    out["A_d_index"] = D - 1
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale end-to-end Cheon attack (contiguous window)
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)
TOY_D = 6


def toy_cheon(trials, mode):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    assert is_prime_trial(n) and (n - 1) % TOY_D == 0
    assert ec_mul(n, G, p, b) is None
    N = (n - 1) // TOY_D
    fac = {}
    mm = n - 1
    for q in range(2, math.isqrt(mm) + 1):
        while mm % q == 0:
            fac[q] = fac.get(q, 0) + 1
            mm //= q
    if mm > 1:
        fac[mm] = fac.get(mm, 0) + 1
    primes = list(fac)

    def is_generator(g):
        if math.gcd(g, n) != 1:
            return False
        return all(pow(g, (n - 1) // q, n) != 1 for q in primes)

    zeta = next(g for g in range(2, 1000) if is_generator(g))
    eta = pow(zeta, TOY_D, n)
    assert pow(eta, N, n) == 1

    m1 = math.isqrt(N) + 1
    M = -(-(N) // m1)
    if mode == "full":
        m = M
        cov = Fraction(1, 1)
    else:
        m = math.isqrt((2 * N) // 5) + 1
        per_min, cov = coverage_exact(N, m1, M, m, brute=True)
    gamma = pow(eta, (-m1) % N, n)

    def attack(s, seed):
        rng = random.Random(seed)
        A1 = ec_mul(s, G, p, b)
        Ad = ec_mul(pow(s, TOY_D, n), G, p, b)
        offset = rng.randrange(M)
        # baby steps: parity-keyed records, explicit merge sort
        e = 1
        recs = []
        mults = 0
        scalars = 0
        for i in range(m):
            pt = ec_mul(e, G, p, b)
            mults += 1
            recs.append(((pt[0], pt[1] % 2), i))
            e = e * eta % n
            scalars += 1
        sort_cmp, sort_moves = merge_sort(recs)
        keys = [r[0] for r in recs]
        assert len(set(keys)) == m  # parity key injective on the baby set
        # contiguous window of m blocks at uniform offset; the window may
        # wrap past block M-1, splitting into two contiguous runs whose
        # scalars are each computed by one multiplication per step, with
        # a single exponentiation per run start (the second run starts
        # at block 0 with scalar 1).
        gscalar = pow(gamma, offset, n)
        scalars += 8
        j = offset
        crossed = False
        probes = 0
        for t in range(m):
            if not crossed and j >= M:
                j -= M
                gscalar = 1
                crossed = True
            Y = ec_mul(gscalar, Ad, p, b)
            mults += 1
            probes += math.ceil(math.log2(2 * m))
            idx = binary_search(recs, (Y[0], Y[1] % 2))
            if idx >= 0:
                i = recs[idx][1]
                k1 = (i + m1 * j) % N
                Qp = ec_mul(pow(zeta, (-k1) % (n - 1), n), A1, p, b)
                mults += 1
                for tt in range(TOY_D):
                    if ec_mul(pow(zeta, (tt * N) % (n - 1), n), G, p, b) == Qp:
                        s_rec = pow(zeta, (k1 + tt * N) % (n - 1), n)
                        mults += TOY_D // 2
                        return (s_rec == s), mults, scalars, probes
            j += 1
            gscalar = gscalar * gamma % n
            scalars += 1
        return None, mults, scalars, probes

    rng = random.Random(2026)
    succ = 0
    wrong = 0
    totm = 0
    tots = 0
    totp = 0
    for ti in range(trials):
        s = rng.randrange(1, n)
        got, mults, scalars, probes = attack(s, 71000 + ti)
        totm += mults
        tots += scalars
        totp += probes
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
        "table": "sorted array keyed by (x, y mod 2), explicit merge sort + binary search",
        "exact_predicted_coverage_with_wrap": str(cov),
        "exact_predicted_coverage_with_wrap_float": float(cov),
        "trials": trials,
        "successes": succ,
        "empirical_success": succ / trials,
        "wrong_extractions": wrong,
        "mean_scalar_mults": totm / trials,
        "mean_scalar_updates": tots / trials,
        "mean_binary_search_probes": totp / trials,
    }


def main():
    data = json.load(sys.stdin)
    trials = int(data.get("trials", 200))
    mode = data.get("mode", "partial")
    assert mode in ("partial", "full")
    out = {
        "production": production_checks(),
        "toy_cheon_v5": toy_cheon(trials, mode),
    }
    print(json.dumps(out))


main()

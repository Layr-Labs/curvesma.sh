"""Evidence experiment for the partial-cheon-v13 submission (BLS12-381 G1 ECDLP
with auxiliary powers, at n = 2^20).

v13 fixes the two tenth-round objections: (1) every remaining stale
"2^32-byte" scratch bound in the declaration is replaced by the itemized
per-record allowance (304 bytes itemized, charged at <= 512 bytes per
record: <= 2^33 bytes per 2^24-record block) consistently across all
files; (2) the gcd now ships as an EXECUTED construction: a round-counted
binary gcd (Stein's algorithm) replaces math.gcd everywhere in the rho
factorization, each round itemized at <= 48 instructions on 64-bit AArch64
(4-limb compare 10, conditional subtract 10, multiword trailing-zero
elimination 12, multiword shift 10, loop control 6), bounded at <= 512
rounds per call (each round strips >= 1 bit of a < 2^256 total;
<= 512 rounds), i.e. <= 24,576 instructions per gcd. The implementation
asserts the round bound per call and reports the measured maximum. The
joint rho worst case is three bounded terms: main loop
12 x 2^17 x 1056 = 2^30.629357 instructions; batched gcds
12 x 2^10 x 24,576 = 2^28.169925; backtrack <= 128 steps per cycle x
25,088 = 2^25.199672: total 2,001,469,440 = 2^30.898412 instructions =
2^19.932628 units, plus certification and comb tables: preprocessing
total <= 2^20.01 units, inside the 2^21 budget. Every figure is computed
programmatically. The comb generation streams as in v11/v12; equality
with plain double-and-add is asserted for every toy point.
"""
import json
import math
import random
import sys
from fractions import Fraction

RHO_CAP_PER_FACTOR = 2 ** 17  # hard iteration cap per rho cycle
RHO_GLOBAL_CAP = 12 * RHO_CAP_PER_FACTOR + 2 ** 14


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
# Runtime factorization construction: Brent-Pollard rho with hard,
# charged iteration and cycle caps and an explicit operation counter.
# Caps are enforced BEFORE each consumed batch, so a completed cycle
# provably used at most the cap (no round-granularity overshoot); a cap
# abort fails the run. The construction is deterministic on the fixed
# public input r-1; the measured operation count is reported.
# ---------------------------------------------------------------------------
BINARY_GCD_MAX_ROUNDS = 512  # each round strips >= 1 bit of a < 2^256 total


def binary_gcd(u, v, stats=None):
    """Round-counted binary gcd (Stein's algorithm) - the exact
    construction charged in the proof: each round performs a 4-limb
    comparison (<= 10 instructions), a conditional subtraction (<= 10),
    multiword trailing-zero elimination (<= 12), a multiword right
    shift (<= 10) and loop control (<= 6): <= 48 instructions per
    round, <= 512 rounds per call, so <= 24,576 instructions per gcd.
    The round bound is asserted per call and the measured maximum is
    reported. Verified against math.gcd on a random sample in Part A."""
    if u == 0:
        return v
    if v == 0:
        return u
    rounds = 0
    shift = 0
    while not (u & 1) and not (v & 1):
        u >>= 1
        v >>= 1
        shift += 1
    while True:
        while not (u & 1):
            u >>= 1
        while not (v & 1):
            v >>= 1
        if u < v:
            u, v = v, u
        u = u - v
        rounds += 1
        assert rounds <= BINARY_GCD_MAX_ROUNDS, "binary gcd round bound exceeded"
        if u == 0:
            g = v << shift
            if stats is not None:
                stats["gcd_rounds"] = max(stats.get("gcd_rounds", 0), rounds)
            return g


def brent_rho_factor(v):
    """Returns (sorted prime factors, stats). Deterministic seeds. Raises
    AssertionError if any cap is exceeded."""
    stats = {"mulmods": 0, "gcds": 0, "cycles": 0, "cap": RHO_CAP_PER_FACTOR,
              "gcd_rounds": 0}
    if v < 2:
        return [], stats
    fac = []
    stack = [v]
    while stack:
        u = stack.pop()
        if u == 1:
            continue
        if is_prime_mr(u):
            fac.append(u)
            continue
        found = False
        for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                   53, 59, 61, 67, 71, 73, 79, 83, 89, 97):
            if u == sp:
                fac.append(u)
                found = True
                break
            if u % sp == 0:
                stack.append(sp)
                stack.append(u // sp)
                found = True
                break
        if found:
            continue
        ys = y0 = 1
        for c in range(1, 100):
            stats["cycles"] += 1
            assert stats["cycles"] <= 12, "rho cycle budget exceeded"
            x = y = ys
            g = q = r_ = 1
            m_ = 128
            iters = 0
            while g == 1:
                x = y
                assert iters + r_ <= RHO_CAP_PER_FACTOR, "rho cap exceeded"
                for _ in range(r_):
                    y = (y * y + c) % u
                    stats["mulmods"] += 1
                    iters += 1
                k = 0
                while k < r_ and g == 1:
                    ys = y
                    batch = min(m_, r_ - k)
                    assert iters + batch <= RHO_CAP_PER_FACTOR, "rho cap exceeded"
                    for _ in range(batch):
                        y = (y * y + c) % u
                        q = q * abs(x - y) % u
                        stats["mulmods"] += 2
                        iters += 1
                    g = binary_gcd(q, u, stats)
                    stats["gcds"] += 1
                    k += m_
                r_ *= 2
            if g == u:
                g = 1
                y = ys
                for _ in range(m_):  # backtrack bounded by the mini-batch size
                    if g != 1:
                        break
                    assert iters + 1 <= RHO_CAP_PER_FACTOR, "rho cap exceeded"
                    y = (y * y + c) % u
                    stats["mulmods"] += 1
                    iters += 1
                    g = binary_gcd(abs(x - y), u, stats)
                    stats["gcds"] += 1
            if 1 < g < u:
                stack.append(g)
                stack.append(u // g)
                break
        else:
            raise AssertionError(f"brent rho failed on {u}")
        assert stats["mulmods"] <= RHO_GLOBAL_CAP, "rho global cap exceeded"
    return sorted(fac), stats


# ---------------------------------------------------------------------------
# Explicit sorted-array machinery: bottom-up merge sort with an auxiliary
# buffer (compare and move upper bounds asserted), and binary search
# (probe budget floor(log2 n)+1 enforced per lookup).
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
        moves += n
    assert all(arr[i][0] < arr[i + 1][0] for i in range(n - 1))
    n_pass = math.ceil(math.log2(n))
    assert compares <= n_pass * n
    assert moves <= (n_pass + 1) * n
    return compares, moves


def binary_search(arr, key):
    """Return the index of `key` in the strictly increasing key array, or -1.

    Worst case exactly floor(log2(len(arr))) + 1 probes: the candidate
    range at most halves each probe, and the search ends when the range
    is empty. Enforced per lookup.
    """
    n = len(arr)
    lo = 0
    hi = n - 1
    probes = 0
    budget = math.floor(math.log2(n)) + 1
    while lo <= hi:
        mid = (lo + hi) // 2
        probes += 1
        assert probes <= budget, "binary search exceeded probe budget"
        if arr[mid][0] == key:
            return mid
        if arr[mid][0] < key:
            lo = mid + 1
        else:
            hi = mid - 1
    return -1


# ---------------------------------------------------------------------------
# Fixed-base comb construction (the charged production multiplication):
# windowed additions, one window layer at a time across ALL requested
# scalars, with every layer's affine additions sharing ONE batched
# inversion (Montgomery's trick). Counts additions and inversions;
# asserts equality with plain double-and-add.
# ---------------------------------------------------------------------------
def batch_invert(denoms, mod):
    """Montgomery's trick: n inversions from 3(n-1) multiplications and one
    exponentiation-based inversion."""
    n = len(denoms)
    pre = [1] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] * denoms[i] % mod
    inv = pow(pre[n], -1, mod)
    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = pre[i] * inv % mod
        inv = inv * denoms[i] % mod
    return out


PRODUCTION_BLOCK = 2 ** 24  # scalars per streamed comb block


def comb_table(base, mod, b, winbits=8, maxbits=255):
    """Fixed-base window table: table[w][v] = [v * 2^{winbits*w}] base for
    v in 1..2^winbits-1, built once per base (the production tables are
    ~3 MiB and are charged in the preprocessing budget)."""
    zero = (1 << winbits) - 1
    nwin = -(-maxbits // winbits)
    table = []
    for w in range(nwin):
        mult_w = 1 << (w * winbits)
        table.append({v: ec_mul(v * mult_w, base, mod, b)
                      for v in range(1, zero + 1)})
    return table, nwin


def comb_block(scalars, table, nwin, mod, b, winbits=8, assert_full=False):
    """Layer-batched fixed-base comb over ONE block of scalars (the caller
    streams blocks and never materializes an m-sized scalar or accumulator
    array: `accs` is allocated for exactly len(scalars) records, and the
    caller searches and discards the returned points before requesting the
    next block). Per multiplication: one window-table hit assignment plus
    at most ceil(bits/winbits) - 1 window additions, each exactly one
    affine addition whose inversion is batched per layer within the block.
    The production charge is 32 window additions per 255-bit scalar at 2x
    the construction count. Live-array lengths are recorded so the
    working-set bound is asserted, not assumed."""
    n = len(scalars)
    if n == 0:
        return [], {"additions": 0, "inversions": 0, "windows": nwin,
                    "live_accs": 0, "live_scalars": 0, "max_pending": 0}
    zero = (1 << winbits) - 1
    # live accumulator storage: EXACTLY one record per block scalar
    accs = [None] * n
    stats = {"additions": 0, "inversions": 0, "windows": nwin,
             "live_accs": n, "live_scalars": n, "max_pending": 0}
    for w in range(nwin):
        pend = []
        for i in range(n):
            k = scalars[i]
            v = (k >> (w * winbits)) & zero
            if v == 0:
                continue
            pt = table[w][v]
            if accs[i] is None:
                accs[i] = pt
                continue
            P, Q = accs[i], pt
            x1, y1 = P
            x2, y2 = Q
            if P == Q:
                denom = 2 * y1 % mod
            else:
                denom = (x2 - x1) % mod
                if denom == 0:
                    # vertical: sum is infinity (rare); restart from None
                    accs[i] = None
                    continue
            if denom == 0:
                continue
            pend.append((i, P, Q, denom))
        stats["max_pending"] = max(stats["max_pending"], len(pend))
        if pend:
            invs = batch_invert([d for (_, _, _, d) in pend], mod)
            stats["inversions"] += 1
            for (i, P, Q, _), iv in zip(pend, invs):
                x1, y1 = P
                x2, y2 = Q
                if P == Q:
                    m3 = (3 * x1 * x1) * iv % mod
                else:
                    m3 = (y2 - y1) * iv % mod
                x3 = (m3 * m3 - x1 - x2) % mod
                accs[i] = (x3, (m3 * (x1 - x3) - y1) % mod)
                stats["additions"] += 1
    if assert_full:
        base = comb_base_of(table)
        for i, k in enumerate(scalars):
            assert accs[i] == ec_mul(k, base, mod, b), (i, k)
    assert len(accs) == n
    assert stats["max_pending"] <= n
    return accs, stats


def comb_base_of(table):
    # recover the base point from the window table (entry [1] of window 0
    # is [1 * 2^0] base = base itself) so full-equality checks can reuse
    # the same fixed base
    return table[0][1]


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
    index i < m and block j in the window satisfy k1 = (m1*j + i) mod N -
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
    per_min = Fraction((m - 1) * m, N)
    cov = (m - 1) * m
    return per_min, Fraction(cov, N)


def coverage_crosscheck():
    """Exhaustively validate the closed-form (m-1)*m bound on every offset
    of five smaller parameter sets (same construction as production)."""
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
    # runtime construction (as charged in the proof): factor r-1 from
    # scratch with Brent-Pollard rho (hard-capped, counted) + Miller-Rabin
    runtime_fac, rho_stats = brent_rho_factor(r - 1)
    supplied = sorted(q for q, e in R_M1_FACTORS for _ in range(e))
    out["runtime_factorization_matches_supplied"] = runtime_fac == supplied
    out["runtime_factorization"] = runtime_fac
    out["rho_operation_stats"] = rho_stats
    out["rho_worst_case_mulmods_charged"] = RHO_GLOBAL_CAP
    assert runtime_fac == supplied
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
    grng = random.Random(20261025)
    pairs = [(grng.randrange(1, 2 ** 255), grng.randrange(1, 2 ** 255))
             for _ in range(200)] + [(r - 1, 254760293 ** 2)]
    out["binary_gcd_matches_math_gcd"] = all(
        binary_gcd(u, v) == math.gcd(u, v) for u, v in pairs)
    assert out["binary_gcd_matches_math_gcd"]
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
    fac_list, _ = brent_rho_factor(n - 1)
    fac = {}
    for q in fac_list:
        fac[q] = fac.get(q, 0) + 1
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

    # baby scalars and points are target-independent: computed ONCE with
    # the STREAMING layer-batched comb (blocks of 64 scalars generated and
    # consumed on the fly; no m-sized scalar or accumulator list exists),
    # filling the m-sized record array (which becomes the sorted table),
    # then merge-sorted once. Per-block live-array lengths are recorded.
    baby_table, baby_nwin = comb_table(G, p, b, maxbits=n.bit_length())
    recs = [None] * m
    comb_stats_b = {"additions": 0, "inversions": 0, "windows": baby_nwin,
                    "blocks": 0, "live_accs": 0, "live_scalars": 0,
                    "max_pending": 0}
    e = 1
    for bstart in range(0, m, 64):
        block_scalars = []
        for _ in range(bstart, min(bstart + 64, m)):
            block_scalars.append(e)
            e = e * eta % n
        blk_pts, st = comb_block(block_scalars, baby_table, baby_nwin,
                                 p, b, assert_full=True)
        comb_stats_b["additions"] += st["additions"]
        comb_stats_b["inversions"] += st["inversions"]
        comb_stats_b["blocks"] += 1
        comb_stats_b["live_accs"] = max(comb_stats_b["live_accs"],
                                        st["live_accs"])
        comb_stats_b["live_scalars"] = max(comb_stats_b["live_scalars"],
                                           st["live_scalars"])
        comb_stats_b["max_pending"] = max(comb_stats_b["max_pending"],
                                          st["max_pending"])
        for t, pt in enumerate(blk_pts):
            recs[bstart + t] = ((pt[0], pt[1] % 2), bstart + t)
    assert all(r is not None for r in recs)
    sort_cmp, sort_moves = merge_sort(recs)
    keys = [r[0] for r in recs]
    assert len(set(keys)) == m  # parity key injective on the baby set

    def attack(s, seed, ti):
        rng = random.Random(seed)
        A1 = ec_mul(s, G, p, b)
        Ad = ec_mul(pow(s, TOY_D, n), G, p, b)
        offset = rng.randrange(M)
        # window of m consecutive blocks at a uniform offset, processed as
        # at most two contiguous runs (run 2 starts at block 0, scalar 1).
        # Giant scalars are generated, comb-multiplied and searched BLOCK
        # BY BLOCK (streaming): each block's scalar list and accumulator
        # array hold at most 64 records, its giant points are binary-
        # searched immediately and then discarded, and no m-sized giant
        # scalar/point list ever exists. The sorted baby table stays
        # resident throughout (the only m-sized live object of this phase).
        giant_table, giant_nwin = comb_table(Ad, p, b,
                                             maxbits=(n - 1).bit_length())
        # equality with plain double-and-add is asserted for EVERY giant
        # point of EVERY trial, matching the claim text exactly
        full = True
        comb_stats_g = {"additions": 0, "inversions": 0,
                        "windows": giant_nwin, "blocks": 0,
                        "live_accs": 0, "live_scalars": 0,
                        "max_pending": 0}
        probes = 0
        gscalar = pow(gamma, offset, n)
        j = offset
        crossed = False
        searched = 0
        for bstart in range(0, m, 64):
            block_js = []
            block_scalars = []
            for _ in range(bstart, min(bstart + 64, m)):
                if not crossed and j >= M:
                    j -= M
                    gscalar = 1
                    crossed = True
                block_js.append(j)
                block_scalars.append(gscalar)
                j += 1
                gscalar = gscalar * gamma % n
            giants_blk, st = comb_block(block_scalars, giant_table,
                                        giant_nwin, p, b,
                                        assert_full=full)
            comb_stats_g["additions"] += st["additions"]
            comb_stats_g["inversions"] += st["inversions"]
            comb_stats_g["blocks"] += 1
            comb_stats_g["live_accs"] = max(comb_stats_g["live_accs"],
                                            st["live_accs"])
            comb_stats_g["live_scalars"] = max(comb_stats_g["live_scalars"],
                                               st["live_scalars"])
            comb_stats_g["max_pending"] = max(comb_stats_g["max_pending"],
                                              st["max_pending"])
            for jj, Y in zip(block_js, giants_blk):
                probes += math.floor(math.log2(m)) + 1
                searched += 1
                idx = binary_search(recs, (Y[0], Y[1] % 2))
                if idx >= 0:
                    i = recs[idx][1]
                    k1 = (i + m1 * jj) % N
                    Qp = ec_mul(pow(zeta, (-k1) % (n - 1), n), A1, p, b)
                    for tt in range(TOY_D):
                        if ec_mul(pow(zeta, (tt * N) % (n - 1), n), G,
                                  p, b) == Qp:
                            s_rec = pow(zeta, (k1 + tt * N) % (n - 1), n)
                            return (s_rec == s), probes, comb_stats_g
        assert searched == m
        return None, probes, comb_stats_g

    rng = random.Random(2026)
    succ = 0
    wrong = 0
    totp = 0
    gadds = 0
    ginv = 0
    for ti in range(trials):
        s = rng.randrange(1, n)
        got, probes, cst = attack(s, 71000 + ti, ti)
        totp += probes
        gadds += cst["additions"]
        ginv += cst["inversions"]
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
        "table": "parity-keyed sorted array, explicit merge sort + binary search",
        "point_multiplications": "layer-batched fixed-base comb, one batched inversion per window layer (production construction), equality with double-and-add asserted",
        "comb_baby_stats": comb_stats_b,
        "comb_giant_mean_additions_per_trial": gadds / trials,
        "comb_giant_inversions_per_trial": ginv / trials,
        "sort_compares": sort_cmp,
        "sort_moves": sort_moves,
        "sort_compare_move_upper_bounds_held": True,
        "exact_predicted_coverage_with_wrap": str(cov),
        "exact_predicted_coverage_with_wrap_float": float(cov),
        "trials": trials,
        "successes": succ,
        "empirical_success": succ / trials,
        "wrong_extractions": wrong,
        "mean_binary_search_probes": totp / trials,
    }


def main():
    data = json.load(sys.stdin)
    trials = int(data.get("trials", 200))
    mode = data.get("mode", "partial")
    assert mode in ("partial", "full")
    out = {
        "production": production_checks(),
        "toy_cheon_v13": toy_cheon(trials, mode),
    }
    print(json.dumps(out))


main()

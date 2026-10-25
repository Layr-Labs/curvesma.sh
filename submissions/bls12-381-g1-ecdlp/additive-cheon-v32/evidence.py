"""Supporting evidence for proof.md (stdlib only; JSON on stdin, JSON on stdout).

additive-cheon-v32: algorithm, walks, ledger and cost figures are IDENTICAL
to the v28 lineage (dependency cbadde77, restated in proof.md; every
identity and bound below is this package's own, re-verified by these runs).
v29 restructured the success argument; v30 certified the residual's strata;
v31 replaced the Selberg majorant by the row sieve; v32 adds the P3/P4
ladder. All load-bearing premises now carry measurement ladders.

1. ledger - every cost term of proof.md section 4 and the exact score.
2. production_certificates - constant structure, mixed-root k=0
   certificate, D_repeat, 2*X1^2*X2 < r + the r=481153 counterexample,
   both swap bounds.
2b. production_routeb - route-B certificates from the three embedded
   VERBATIM collaborator modules: ledger 0.159313, route-B floor at the
   NAMED coprime premise 53/1250; external analytic input: Dusart
   arXiv:1002.0442v1 Thm 6.9 eq (6.6) only.
2c. coprime_measurement - the coprime ladder (record:
   measurement_campaign_report.md): exact controls, HT targets at the
   production d = 2^20 with the ACTUAL masks and original unpruned
   prefixes, dense complement; fractions recomputed from raw counts.
2d. p34_measurement - the P3/P4 ladder (record: p34-campaign-report.txt):
   swap-pruned ratio-occupancy moments, exact controls + 64-stratum HT
   estimates, margins in design SEs.
3.-5. coverage / structured_lambda_point / toy_attack / hash_probe - the
   v28/v29/v30 measurements, unchanged.

Toy results are evidence at their scale only, not production cost or success
proofs. Input: {"K", "cov_bits", "cov_logd", "atk_bits", "atk_logd", "trials", "gmax", "seed"}.
"""
import cmath
import json
import math
import random
import sys
from collections import Counter
from fractions import Fraction as F
C = 2000
r = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
d = 1 << 20
m = (r - 1) // d
assert (r - 1) % d == 0 and m % 6 == 0


def ledger(K):
    K = F(K)

    # per-element instruction budgets (proof.md, section 4, v2).
    # CANON: published "p - y (~7), compare (~12), select (6): ~25";
    #   v2 charges 30 (v1: 40).
    # HASH_OP (v3, per-transfer): 2^116.7 slots x 72-byte records at load
    #   1/2 exceed one 2^63-byte space: sharded across <= 2^60 nodes
    #   behind a radix-1024 fabric, 6 stages.
    #   Endpoint: shard+slot 4; append 9-word record 14; receive 9;
    #   probing at 1/2 (Knuth 2.5 probes x 6) 15; store 9; loop 10 -
    #   total 61. EACH relay: load+store 9 words 18; destination (10-bit
    #   field, port load) 4; queue 4; loop 3 - 29 per stage x 6 = 174.
    #   TOTAL 61 + 174 = 235, charged 240.
    # WALK: read b+gap 4; overflow 2; address 3; load 12; load+store
    #   points 24; batch-inversion buffers beyond 6M 24; loop 5 - ~75;
    #   charged 100.
    CANON, HASH_OP, WALK, FR_MUL = 30, 240, 100, 200
    c_B = 1 + F(CANON + HASH_OP + WALK, C)   # walk add + canonical y + 1 insert + bookkeeping
    c_G = 1 + F(CANON + HASH_OP + WALK, C)   # walk add + canonical y + 1 lookup + bookkeeping

    # sizes: 6*N_B*N_G >= K*m, N_B/N_G ~ c_G/c_B
    N_B = math.isqrt(int(K * m * c_G / (6 * c_B))) + 1
    N_G = -(-int(K * m) // (6 * N_B)) + 1
    assert 6 * N_B * N_G >= K * m

    X1, X2, GMAX = 1 << 82, 1 << 70, 1 << 34   # X1 up for the larger N_A; GMAX rebalanced
    SCALAR = 1000          # units: one non-batched scalar multiplication by a <= 2^82-bit scalar (proof.md 4.3)
    SCALAR_FULL = 3000     # units: one scalar multiplication by a 255-bit scalar
    PI_X2 = F(126, 100) * X2 / F(int(70 * math.log(2) * 10**6), 10**6)    # pi(x) < 1.26 x / ln x
    # proven: pi(X1) >= X1/ln X1, pi(X1/2) <= 1.25506*(X1/2)/ln(X1/2)
    # (Rosser-Schoenfeld)
    PI_X1 = F(126, 100) * X1 / F(int(82 * math.log(2) * 10**6), 10**6)   # pi(X1) upper estimate (cost charge)
    # proven (Rosser-Schoenfeld): pi(X1) >= X1/ln X1 (x >= 17),
    # pi(X1/2) <= 1.25506*(X1/2)/ln(X1/2) (x > 1), so (X1/2, X1] holds
    P_AVAIL = (F(X1) / F(int(82 * math.log(2) * 10**6), 10**6)
               - F(125506, 100000) * F(X1, 2) / F(int(81 * math.log(2) * 10**6), 10**6))
    N_A = (N_B + N_G) * 2 * d * 49 // X2 + 2   # primes a used; the b-lists hold only the class PRIMES,
                                               # each >= X2/(2*d*49) (H2; ln(2^70) = 48.5 < 49, 1/2 slack)
    # shared-class-prefix bound: a same-class pair with a consumed prefix
    # of t collapses its t diagonal pairs to u*a/a' (losing 6*(t-1)
    # values); prefix <= L2 = 2x mean class size (H2, factor-2 slack)
    L2 = 2 * int(X2 / (d * math.log(X2)))
    DUPB = 6 * (N_A // 2) ** 2 * max(1, L2 - 1) // d
    assert 2 * N_A < P_AVAIL            # proven interval supply covers N_A with > 2x slack

    terms = {
        "baby walk (N_B * c_B)": N_B * c_B,
        "giant walk (N_G * c_G)": N_G * c_G,
        "gap overflow additions (<= N_A * X2 / GMAX)": F(N_A * X2 // GMAX),
        "per-prime tables [k]Q_a, k <= GMAX": F(N_A * GMAX),
        "per-prime scalar mults ([a]base, first element, doubling)": F(N_A * 3 * SCALAR),
        "chi for primes <= X1 (100 units each)": 100 * PI_X1,
        "prime sieve to X1, chi(b) sieve and class lists to X2": F(X1 * 64 + X2 * 200, C),
        "collision verification (6 candidates)": F(6 * SCALAR_FULL),
        "d-th root in F_r (<= 2^20 F_r mults of 200 instr)": F((1 << 20) * 200, C),
        "mu_d search (2^20 scalar mults) + final check": F(((1 << 20) + 1) * SCALAR_FULL),
        "input decoding and subgroup checks": F(4 * SCALAR_FULL),
    }
    W = sum(terms.values())
    pre_keys = ["baby walk (N_B * c_B)", "per-prime tables [k]Q_a, k <= GMAX",
                "per-prime scalar mults ([a]base, first element, doubling)",
                "chi for primes <= X1 (100 units each)", "prime sieve to X1, chi(b) sieve and class lists to X2",
                "gap overflow additions (<= N_A * X2 / GMAX)"]
    pre = sum(terms[k] for k in pre_keys)         # challenge-independent work (upper bound: counts both walks' tables)
    c_kkm = 2 + F(CANON + HASH_OP + WALK + FR_MUL, C)   # KKM c=3: 2 additions + one F_r mul to advance zeta^j
    kkm = 2 * c_kkm * F(math.isqrt(int(K * m / 6)))
    log2W = math.log2(W.numerator) - math.log2(W.denominator)
    score = F(math.ceil(log2W * 10**6), 10**6)          # rounded up
    success = 1 - math.exp(-float(6 * N_B * N_G) / m)
    out = {
        "K": str(K), "c_B": str(c_B), "c_G": str(c_G),
        "log2_m": round(math.log2(m), 6), "log2_N_B": round(math.log2(N_B), 6), "log2_N_G": round(math.log2(N_G), 6),
        "log2_primes_a_used": round(math.log2(N_A), 3),
        "terms_log2": {k: round(math.log2(v.numerator) - math.log2(v.denominator), 3) for k, v in terms.items()},
        "log2_W": round(log2W, 6),
    "log2_preprocessing": round(math.log2(pre.numerator) - math.log2(pre.denominator), 3),
    "log2_W_KKM_c3_same_budgets": round(math.log2(kkm.numerator) - math.log2(kkm.denominator), 3),
        "score": {"numerator": str(score.numerator), "denominator": str(score.denominator)},
        "heuristic_success_floor": math.floor(success * 1000) / 1000,
        "success_claim_law": math.floor(success * 985) / 1000,
        "log2_primes_interval_lower_bound": round(math.log2(float(P_AVAIL)), 3),
        "log2_prime_availability_ratio_lower": round(math.log2(float(P_AVAIL) / N_A), 3),
        "law_slowdown": 0.985,
        "log2_prefix_duplicate_values_bound": round(math.log2(max(1, DUPB)), 1),
        "prefix_duplicate_fraction_log2": round(math.log2(max(1, DUPB) / (N_B * N_G)), 1),
        "success_claim": math.floor(success * 985) / 1000,
        "log2_memory_bytes": round(math.log2(2 * N_B * 72 + X2 * 9 + 2**46 * 96), 3),
    }
    return out


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    e, s = n - 1, 0
    while e % 2 == 0:
        e //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, e, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_in(lo, hi):
    return [q for q in range(lo, hi + 1) if is_prime(q)]


def Omega_small(n):
    """Omega with multiplicity by trial division (toy-scale constants only)."""
    o, t = 0, n
    p = 2
    while p * p <= t:
        while t % p == 0:
            t //= p
            o += 1
        p += 1
    return o + (1 if t > 1 else 0)


def ord_mod(x, r):
    t, k = x % r, 1
    v = t
    while v != 1:
        v = v * t % r
        k += 1
    return k


def sqrt_mod(a, p):
    """Tonelli-Shanks; None if a is a non-residue."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    mm, c, t, x = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, x = i, b * b % p, t * b * b % p, x * b % p
    return x


class Curve:
    def __init__(self, p, c):
        self.p, self.c, self.adds = p, c, 0

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def walk_add(self, P, Q):
        """Walk step: always the incomplete case x1 != x2 (asserted, see proof.md 3)."""
        assert P is not None and Q is not None and P[0] != Q[0]
        self.adds += 1
        return self.add(P, Q)

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R


def find_curve(bits, d, rng):
    """p = 1 mod 3; group order found by BSGS in the Hasse interval; need prime r with 6d | r-1."""
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if p % 3 != 1 or not is_prime(p):
            continue
        for c in range(1, 8):
            E = Curve(p, c)
            while True:
                x = rng.randrange(p)
                y = sqrt_mod((x ** 3 + c) % p, p)
                if y is not None:
                    P = (x, y)
                    break
            lo, hi = p + 1 - 2 * math.isqrt(p) - 2, p + 1 + 2 * math.isqrt(p) + 2
            w = math.isqrt(hi - lo) + 1
            baby, R = {}, None
            for j in range(w):
                if R is not None:
                    baby.setdefault(R, j)
                R = E.add(R, P)
            step = E.mul(w, P)
            step = (step[0], (-step[1]) % p) if step else None
            Gk, cands = E.mul(lo, P), []
            for i in range(w + 1):
                if Gk is None:
                    cands.append(lo + i * w)
                elif (Gk[0], (-Gk[1]) % p) in baby:
                    cands.append(lo + i * w + baby[(Gk[0], (-Gk[1]) % p)])
                Gk = E.add(Gk, step)
            cands = sorted({N for N in cands if lo <= N <= hi and E.mul(N, P) is None})
            if len(cands) == 1 and is_prime(cands[0]) and (cands[0] - 1) % (6 * d) == 0:
                return E, cands[0], P


def dth_root(t, d, r):
    """A d-th root (d = 2^k, d | r-1) of a d-th power t, by k square roots choosing the branch that stays
    a power of the remaining 2-power order."""
    k = d.bit_length() - 1
    for j in range(k, 0, -1):
        y = sqrt_mod(t, r)
        cand = [y, r - y]
        e = (r - 1) >> (j - 1)
        t = next(z for z in cand if j == 1 or pow(z, e, r) == 1)
    return t


def hash_probe_check(cfg):
    """Toy evidence for the section-4 hashing model: canonical-y keys of
    the structured walk into a sharded open-addressing table at load
    1/2; probes and shard loads vs the Knuth bounds (2.5 unsuccessful at
    1/2; shard loads Binomial(n, 1/S))."""
    bits, n, shards, seed = cfg["bits"], cfg["n"], cfg["shards"], cfg["seed"]
    rng = random.Random(seed + 10 ** 8)
    d = 8
    while True:                          # the production fact 2 not-in H is a
        E, r_toy, P = find_curve(bits, d, rng)   # walk precondition: retry until it holds
        m = (r_toy - 1) // d
        if pow(2, m, r_toy) != 1:
            break
    p, m = E.p, (r_toy - 1) // d
    gen = next(g for g in range(2, 10 ** 6) if pow(g, (r_toy - 1) // 2, r_toy) != 1)
    zeta = pow(gen, m, r_toy)
    zlog = {pow(zeta, i, r_toy): i for i in range(d)}
    chi = lambda x: zlog[pow(x, m, r_toy)]
    X2 = 64 * d
    classes = {}
    pset = set(primes_in(2, X2))
    for b in sorted(pset):              # production walk: the b's are the PRIMES of the class
        classes.setdefault(chi(b), []).append(b)
    X1 = 2 * X2 + 2
    while True:
        plist = primes_in(X1 // 2 + 1, X1)
        if sum(len(classes.get((-chi(a)) % d, [])) for a in plist) >= 2 * n:
            break
        X1 *= 2
        if X1 * X2 >= r_toy:
            raise RuntimeError("toy too small")
    keys, seen = [], set()
    gmax = 64

    def emit(Q):
        if Q is None:
            return
        y = min(Q[1], p - Q[1])
        if y not in seen:
            seen.add(y)
            keys.append(y)

    for a in plist:                      # increasing order: the production rule
        bs = classes.get((-chi(a)) % d, [])
        if not bs:
            continue
        Qa = E.mul(a, P)
        tab = {1: Qa, 2: E.add(Qa, Qa)}  # the production walk: T[2] by
        for k in range(3, gmax + 1):     # doubling, then T[k] = T[k-1] + Q
            tab[k] = E.walk_add(tab[k - 1], Qa)
        # class list increasing, ONE addition per element (overflow
        # chunks above gmax), as in production
        emit(E.mul(bs[0], Qa))
        cur = E.mul(bs[0], Qa)
        for prev, b in zip(bs, bs[1:]):
            gap = b - prev
            while gap > gmax:
                cur = E.walk_add(cur, tab[gmax])
                gap -= gmax
            cur = E.walk_add(cur, tab[gap])
            emit(cur)
            if len(keys) >= 2 * n:
                break
        if len(keys) >= 2 * n:
            break
    slots = -(-2 * n // shards)          # load 1/2 overall
    sh_shift = slots.bit_length()
    sh_of = lambda y: (y >> sh_shift) % shards
    tables = [dict() for _ in range(shards)]
    probes_ins = 0
    shard_load = [0] * shards
    for y in keys[:n]:
        sh = sh_of(y)
        home = y % slots
        t = tables[sh]
        sl, probes = home, 0
        while sl in t:
            sl = (sl + 1) % slots
            probes += 1
        probes_ins += probes + 1
        t[sl] = y
        shard_load[sh] += 1
    probes_miss = 0
    for y in keys[n:]:
        sh = sh_of(y)
        home = y % slots
        t = tables[sh]
        sl, probes = home, 0
        while sl in t and t[sl] != y:
            sl = (sl + 1) % slots
            probes += 1
        probes_miss += probes + 1
    mean_ins = probes_ins / n
    mean_miss = probes_miss / n
    mean_shard = n / shards
    sig_shard = (n * (1 / shards) * (1 - 1 / shards)) ** 0.5
    max_dev = max(abs(x - mean_shard) for x in shard_load)
    ok = (mean_ins <= 2.5 * 1.2 and mean_miss <= 2.5 * 1.2
          and max_dev <= 5 * sig_shard)
    return {"curve_bits": round(math.log2(E.p), 2), "keys": "structured walk "
            "(primes in increasing order, class lists in increasing order)",
            "keys_inserted": n, "shards": shards, "slots_per_shard": slots,
            "load": 0.5, "mean_probes_insert": round(mean_ins, 3),
            "mean_probes_absent_lookup": round(mean_miss, 3),
            "knuth_unsuccessful_bound": 2.5,
            "max_shard_load_dev_sigma": round(max_dev / sig_shard, 2),
            "within_model": ok}


def toy_attack(cfg):
    bits, logd, trials, K, gmax = cfg["bits"], cfg["logd"], cfg["trials"], float(cfg["K"]), cfg["gmax"]
    rng = random.Random(cfg["seed"])
    d = 1 << logd
    E, r, P = find_curve(bits, d, rng)
    p, m = E.p, (r - 1) // d
    beta = next(z for z in (pow(g, (p - 1) // 3, p) for g in range(2, 100)) if z != 1)
    lam = next(L for L in (pow(g, (r - 1) // 3, r) for g in range(2, 100))
               if L != 1 and E.mul(L, P) == (beta * P[0] % p, P[1]))
    U = [1, r - 1, lam, r - lam, lam * lam % r, r - lam * lam % r]
    gen = next(g for g in range(2, 10**6) if pow(g, (r - 1) // 2, r) != 1)
    zeta = pow(gen, m, r)
    zlog = {pow(zeta, i, r): i for i in range(d)}
    chi = lambda x: zlog[pow(x, m, r)]
    size = math.isqrt(int(K * m / 6)) + 1
    X2 = 64 * d
    classes = {}
    pset = set(primes_in(2, X2))
    for b in sorted(pset):              # production walk: the b's are the PRIMES of the class
        classes.setdefault(chi(b), []).append(b)
    LN = int(math.log(X2)) + 1
    ka = 2 * size * d * LN // X2 + 1    # production block rule: ceil(N_A/2), prime lists
    X1 = 2 * X2 + 2
    while True:
        primes = primes_in(X1 // 2 + 1, X1)
        if (len(primes) >= 2 * ka + 4
                and sum(len(classes.get((-chi(a)) % d, [])) for a in primes[:2 * ka]) >= 2 * size):
            break
        X1 *= 2

    def walk(base, plist, count, emit):
        produced = 0
        for a in plist:
            bs = classes.get((-chi(a)) % d, [])
            if not bs:
                continue
            Qa = E.mul(a, base)                       # per-prime setup (charged separately in proof.md)
            tab = {1: Qa, 2: E.add(Qa, Qa)}           # T[2] by doubling, then T[k] = T[k-1] + Q (k >= 3)
            for k in range(3, gmax + 1):
                tab[k] = E.walk_add(tab[k - 1], Qa)
            cur = E.mul(bs[0], Qa)
            emit(a * bs[0], cur)
            produced += 1
            for prev, b in zip(bs, bs[1:]):
                if produced >= count:
                    return
                gap, at = b - prev, prev
                while gap > gmax:                     # overflow chunks of GMAX
                    if at == gmax:                    # [GMAX]Q + [GMAX]Q: the one doubling case (proof.md 2)
                        cur = E.add(cur, tab[gmax])
                        stats["doublings"] += 1
                    else:
                        cur = E.walk_add(cur, tab[gmax])
                        stats["chunks"] += 1
                    gap -= gmax
                    at += gmax
                cur = E.walk_add(cur, tab[gap])
                emit(a * b, cur)
                produced += 1
            if produced >= count:
                return
        raise RuntimeError("not enough primes")

    stats = {"chunks": 0, "doublings": 0}
    succ = walk_adds = 0
    fixed = bool(cfg.get("fixed"))    # True: the deterministic production rule
    if fixed:
        primes.sort()                 # increasing order, fixed across trials
    for _ in range(trials):
        s = rng.randrange(1, r)
        A1, Ad = E.mul(s, P), E.mul(pow(s, d, r), P)
        E.adds = 0
        if not fixed:
            rng.shuffle(primes)
        table, found = {}, []

        canon = lambda y: min(y, p - y)            # same key for the 6 points [u]Q, u in U

        def baby(v, Pt):
            table[canon(Pt[1])] = v

        def giant(v, Pt):
            key = canon(Pt[1])
            if not found and key in table:
                found.append((table[key], v))

        walk(P, primes[:ka], size, baby)
        walk(Ad, primes[ka:2 * ka], size, giant)
        walk_adds += E.adds
        if not found:
            continue
        b, g = found[0]
        t = next(c for c in (u * b * pow(g, -1, r) % r for u in U) if E.mul(c, P) == Ad)
        s0 = dth_root(t, d, r)
        w = pow(gen, (r - 1) // d, r)
        sol = next(s0 * pow(w, i, r) % r for i in range(d) if E.mul(s0 * pow(w, i, r) % r, P) == A1)
        assert sol == s, "wrong scalar"
        succ += 1
    out = {"p_bits": round(math.log2(p), 2), "r": r, "d": d, "set_size": size,
           "K_effective": round(6 * size * size / m, 4),
           "predicted_success": round(1 - math.exp(-6 * size * size / m), 4),
           "conservative_gcd_bound_toy_scale": round(1 - math.exp(-6 * size * size / m * (1 - 6 / math.pi ** 2 * (
               math.log(X2) + 2 * 0.5772156649 - 1) / d * 2)), 4),
           "trials": trials, "successes": succ, "success_rate": round(succ / trials, 4),
           "walk_additions_per_trial": round(walk_adds / trials, 1), "2*set_size": 2 * size,
           "gmax": gmax, "overflow_chunk_additions": stats["chunks"], "overflow_doublings": stats["doublings"],
           "all_recovered_scalars_correct": True,
           "note": "the gcd field is the INTEGER-LIST repeated-ratio estimate, kept as loose context only; "
                   "the prime construction's integer-level family is bounded separately (proof.md section 3); "
                   "at production scale the integer-list estimate is 2^-14.1",
           "scope": "toy curve only; does not validate production costs"}
    return out


def coverage_prime_point(bits, logd, ratio, K, seed, baselines=False, list_filter="primes", nulls=0,
                        energy=False, moments=False, fixed_r=None):
    """Exact |U*B*G^-1|/m for the production construction (chi-class PRIME
    lists; baby/giant blocks of (X1/2, X1]). list_filter="integers": the
    same r, blocks and n with full integer lists. energy/moments/
    baselines/nulls=k add their eponymous measurements, counted by the
    same routine."""
    rng = random.Random(seed)
    d = 1 << logd
    if fixed_r is None:
        while True:
            r = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * d) * (6 * d) + 1
            if is_prime(r) and ((r - 1) // d) % 6 == 0:
                break
    else:
        r = fixed_r
        assert is_prime(r) and (r - 1) % (6 * d) == 0 and ((r - 1) // d) % 6 == 0
    m = (r - 1) // d
    gen = next(g for g in range(2, 10**6) if all(pow(g, (r - 1) // q, r) != 1 for q in (2, 3)))
    lam = pow(gen, (r - 1) // 3, r)
    U = [1, r - 1, lam, r - lam, lam * lam % r, r - lam * lam % r]
    zeta = pow(gen, m, r)
    zlog = {pow(zeta, i, r): i for i in range(d)}
    chi = lambda x: zlog[pow(x, m, r)]
    n = math.isqrt(int(K * m / 6)) + 1
    X2 = ratio * d
    classes = {}
    if list_filter == "primes":
        for b in sorted(set(primes_in(2, X2))):
            classes.setdefault(chi(b), []).append(b)
    else:
        for b in range(1, X2 + 1):
            classes.setdefault(chi(b), []).append(b)
    ell_min = min((len(classes.get(c, [])) for c in range(d)), default=0)
    LN = int(math.log(X2)) + 1
    ka = 2 * n * d * LN // X2 + 1          # the declared block rule, exactly
    X1 = 2 * X2 + 2
    while True:
        P = primes_in(X1 // 2 + 1, X1)
        if (len(P) >= 2 * ka + 4
                and sum(len(classes.get((-chi(a)) % d, [])) for a in P[:2 * ka]) >= 2 * n):
            break
        X1 *= 2
    P.sort()
    def build(block):
        S, used = [], {}
        for a in block:
            for b in classes.get((-chi(a)) % d, []):
                S.append(a * b)
                used.setdefault(a, []).append(b)
                if len(S) == n:
                    return S, used
        raise RuntimeError("not enough primes")
    (B, ub), (G, ug) = build(P[:ka]), build(P[ka:2 * ka])
    lost_values = 0
    for a, bs in ub.items():
        for a2, bs2 in ug.items():
            if chi(a) == chi(a2):
                t = len(set(bs) & set(bs2))
                if t:
                    lost_values += 6 * (t - 1)
    assert all(pow(v, m, r) == 1 for v in B + G) and len(set(B)) == n and len(set(G)) == n
    def count(Bs, Gs):
        Gs_inv = [pow(g, -1, r) for g in Gs]
        bitmap = bytearray(r // 8 + 1)
        for u in U:
            for b in Bs:
                ub_ = u * b % r
                for gi in Gs_inv:
                    v = ub_ * gi % r
                    bitmap[v >> 3] |= 1 << (v & 7)
        c = 0
        for i in range(0, len(bitmap), 1 << 20):
            c += int.from_bytes(bitmap[i:i + (1 << 20)], "little").bit_count()
        return round(c / m, 4)
    pred = round(1 - math.exp(-6 * n * n / m), 4)
    out = {"d": d, "r_bits": round(math.log2(r), 2), "X2": X2, "X1": X1, "set_size": n,
           "ell_min": ell_min, "ka": ka, "primes_used_baby": len(ub), "primes_used_giant": len(ug),
           "prefix_lost_values": lost_values,
           "prefix_lost_fraction": round(lost_values / (6 * n * n), 6),
           "predicted": pred,
           "no_wrap": X1 * X2 < r,
           "prime_product": count(B, G)}
    out["deficit"] = round(pred - out["prime_product"], 4)
    if energy and list_filter == "primes":
        # ---- exact distinct-index energy decomposition (v29) ----
        swap_pairs = 0
        for a, bs in ub.items():
            for a2, bs2 in ug.items():
                if chi(a) == chi(a2):
                    s = min(len(bs), len(bs2))
                    swap_pairs += s * (s - 1)
        # total E over all 6n^2 triples
        Ginv = {g: pow(g, -1, r) for g in G}
        cval = Counter()
        for u in U:
            for x in B:
                ux = u * x % r
                for g, gi in Ginv.items():
                    cval[ux * gi % r] += 1
        E_total = sum(c * (c - 1) for c in cval.values())
        # per-ratio family totals (each rho carries 6 ordered pairs)
        cmod = Counter()
        for x in B:
            for g, gi in Ginv.items():
                cmod[x * gi % r] += 1
        fam_total = {}
        for rho in {u2 * pow(u, -1, r) % r for u in U for u2 in U}:
            rinv = pow(rho, -1, r)
            cnt = sum(cv * cmod.get(h * rinv % r, 0) for h, cv in cmod.items())
            fam_total[rho] = 6 * cnt
        assert sum(fam_total.values()) - 6 * n * n == E_total
        # mixed-root k=0 identities, exactly, via Omega(u) = Omega(u')
        prod = Counter(x * y for x in B for y in G)
        mixed = 0
        mixed_by_rho = Counter()
        Om = {u: Omega_small(u) for u in U}
        for u in U:
            for u2 in U:
                if u == u2 or Om[u] != Om[u2]:
                    continue
                rho = u2 * pow(u, -1, r) % r
                for L, c in prod.items():
                    if (u * L) % u2 == 0:
                        mm = c * prod.get(u * L // u2, 0)
                        if mm:
                            mixed += mm
                            mixed_by_rho[rho] += mm
        E_kneq = E_total - 6 * swap_pairs - mixed
        base = 6 * n * n * (6 * n * n) / m
        fam = {}
        for rho, tot in fam_total.items():
            o = ord_mod(rho, r)
            k0 = 6 * swap_pairs if rho == 1 else mixed_by_rho.get(rho, 0)
            kneq = tot - (6 * n * n + k0) if rho == 1 else tot - k0
            fam[hex(rho)] = {"order": o, "family_kneq0": kneq,
                             "rate_vs_null": round(kneq / (6 * n ** 4 / m), 6)}
        assert sum(v["family_kneq0"] for v in fam.values()) == E_kneq
        # total-energy excess (U-filter identity) and the pruned factorial
        # moments: one representative per swap value preserves coverage;
        # fourth-order Bonferroni on the pruned multiplicities is the
        # model-free route (proof.md 3.7). Self/mixed fourth moments by
        # exact U-coset orthogonality; mixed value checked exactly.
        cvalp = Counter(cval)
        pruned = 0
        for a, bs in ub.items():
            for a2, bs2 in ug.items():
                if chi(a) == chi(a2):
                    s = min(len(bs), len(bs2))
                    if s > 1:
                        vv = a * pow(a2, -1, r) % r
                        for u in U:
                            cvalp[u * vv % r] -= s - 1
                            pruned += s - 1
        assert pruned == lost_values and all(v >= 1 for v in cvalp.values())
        Fj = [0] * 5
        for c in cvalp.values():
            f = c
            for j in range(1, 5):
                Fj[j] += f
                f *= c - j
        F1, F2, F3, F4 = (Fj[j] / m for j in (1, 2, 3, 4))
        K0 = F1
        bonf4 = F1 - F2 / 2 + F3 / 6 - F4 / 24
        assert bonf4 <= out["prime_product"] + 2e-4
        hsect = m // 6
        Binv = {x: pow(x, -1, r) for x in B}
        cBB = Counter()
        for x in B:
            for x2, gi in Binv.items():
                cBB[x * gi % r] += 1
        cGG = Counter()
        for y in G:
            for y2, gi in Ginv.items():
                cGG[y * gi % r] += 1
        def sectm(cA, cC):
            ginv = {g: pow(g, -1, r) for g in cA}
            return sum(cc * cC.get(u * ginv[g] % r, 0)
                       for u in U for g, cc in cA.items())
        SB = hsect * sectm(cBB, cBB) - n ** 4
        SG = hsect * sectm(cGG, cGG) - n ** 4
        Sig = hsect * sectm(cBB, cGG) - n ** 4
        ufilter_num = m * (E_total + 6 * n * n) - 36 * n ** 4
        assert ufilter_num % 36 == 0 and ufilter_num // 36 == Sig
        # diagonal bound: self-moment bounds + Cauchy-Schwarz need
        # eps >= 1/K - 1 = 0.818 - blocked, verified per run
        diag_floor = hsect * (2 * n * n - n) - n ** 4
        assert SB >= diag_floor and SG >= diag_floor
        out["energy"] = {
            "E": E_total, "six_swap_pairs": 6 * swap_pairs, "mixed_root_k0": mixed,
            "E_kneq0": E_kneq, "N2_over_m": round(base, 4),
            "excess_pct": round(100 * (E_kneq / base - 1), 4),
            "total_excess_pct": round(100 * (E_total / base - 1), 4),
            "families": fam,
            "pruned_factorial_moments": {
                "K0": round(K0, 6),
                "F2_over_K0sq": round(F2 / K0 / K0, 4),
                "F3_over_K0cu": round(F3 / K0 ** 3, 4),
                "F4_over_K0qu": round(F4 / K0 ** 4, 4),
                "bonferroni4_floor": round(bonf4, 6),
                "bonferroni4_bounds_coverage": True,
            },
            "sector_moments": {
                "S_B": SB, "S_G": SG, "mixed_sector_sum": Sig,
                "ufilter_identity_exact": True,
                "cauchy_ratio": round(math.sqrt(SB * SG) / Sig, 3),
            },
        }
    if moments and list_filter == "primes":
        # ---- class-marginal second moments (swap calibration) ----
        w = cmath.exp(2j * math.pi / d)
        rowsB = {a: len(bs) for a, bs in ub.items()}
        rowsG = {a: len(bs) for a, bs in ug.items()}
        def dw(rows):
            return sum(abs(sum(t * w ** (j * chi(a) % d) for a, t in rows.items())) ** 2
                       for j in range(1, d))
        dB, dG = dw(rowsB) / n / n, dw(rowsG) / n / n
        nB = Counter(chi(a) % d for a in rowsB)
        nG = Counter(chi(a) % d for a in rowsG)
        sun = sum(nB[c] * nG.get(c, 0) for c in nB)
        dUB = sum(abs(sum(w ** (j * chi(a) % d) for a in rowsB)) ** 2 for j in range(1, d))
        dUG = sum(abs(sum(w ** (j * chi(a) % d) for a in rowsG)) ** 2 for j in range(1, d))
        T = len(P)
        D0 = sum(abs(sum(w ** (j * chi(a) % d) for a in P)) ** 2 for j in range(1, d))
        swap6 = 6 * sum(min(t, t2) * (min(t, t2) - 1)
                        for a, t in rowsB.items() for a2, t2 in rowsG.items() if chi(a) == chi(a2))
        out["moments"] = {
            "delta_w_baby": round(dB, 3), "delta_w_giant": round(dG, 3),
            "sum_nB_nG": sun,
            "unweighted_bound": round((len(rowsB) * len(rowsG) + math.sqrt(dUB * dUG)) / d, 2),
            "interval_T": T, "Delta0_over_T2": round(D0 / T / T, 4),
            "interval_bound": round((T * T + D0) / (4 * d), 2),
        }
        WB, WG = Counter(), Counter()
        for a, t in rowsB.items():
            WB[chi(a) % d] += t
        for a, t in rowsG.items():
            WG[chi(a) % d] += t
        out["moments"]["swap_bound_6WBWG"] = 6 * sum(WB[c] * WG.get(c, 0) for c in WB)
        out["moments"]["six_swap_pairs"] = swap6
        out["moments"]["wform_bound_check"] = bool(
            swap6 <= out["moments"]["swap_bound_6WBWG"]
            and sun <= out["moments"]["unweighted_bound"]
            and sun <= out["moments"]["interval_bound"])
    if baselines or nulls:
        rnd = lambda: [pow(rng.randrange(1, r), d, r) for _ in range(n)]
        if baselines:
            out["random"] = count(rnd(), rnd())
        if nulls:
            # the law's null: independent random set-pairs, counted by
            # the same routine at the same parameters
            draws = [count(rnd(), rnd()) for _ in range(nulls)]
            out["random_nulls"] = draws
            out["random_null_mean"] = round(sum(draws) / nulls, 4)
            out["random_null_min"] = min(draws)
            out["random_null_max"] = max(draws)
            sd = (sum((x - sum(draws) / nulls) ** 2 for x in draws) / nulls) ** 0.5
            out["random_null_sd"] = round(sd, 5)
            out["prime_within_nulls"] = out["random_null_min"] <= out["prime_product"] <= out["random_null_max"]
    if baselines is True:          # the interval baseline only where cheap
        small, x = [], 1
        while len(small) < 2 * n:
            if pow(x, m, r) == 1:
                small.append(x)
            x += 1
        out["interval"] = count(small[0::2], small[1::2])
    return out


def production_certificates(K):
    """The v29 production certificates, all arithmetic: the constant
    structure (L = x^2 - 1, L^2 + L + 1 = r; six U-representatives,
    orders 1,3,6,3,6,2); the cutoff-free mixed-root k=0 certificate
    (all 15 mixed pairs fail Omega(v) = Omega(w) <= 4 after coprime
    reduction); D_repeat; the retained 2*X1^2*X2 < r + the r=481153
    counterexample; the two swap bounds. Details in proof.md sections
    3.2-3.5."""
    z = 0xd201000000010000
    L = z * z - 1
    assert L * L + L + 1 == r and (r - 1) % (1 << 32) == 0
    reps = {1: "1", L: "L", L + 1: "L+1", L * L: "L^2", L * L + 1: "L^2+1", L * (L + 1): "L(L+1)"}
    assert set(reps) == {1, r - 1, L, r - L, L * L % r, r - L * L % r}
    orders = sorted((reps[u], ord_mod(u, r)) for u in reps)
    # (2) the 15 unordered mixed pairs
    cert, dead = [], 0
    vals = sorted(reps)
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            u, u2 = vals[i], vals[j]
            g = math.gcd(u, u2)
            v, w2 = u // g, u2 // g
            if (v == 1) != (w2 == 1):
                why = "reduced coefficient 1 against Omega >= 1"
            elif v % (1 << 31) == 0 or w2 % (1 << 31) == 0:
                why = "reduced coefficient divisible by 2^31 (Omega >= 31)"
            else:
                assert u2 == L * L + 1 or u == L * L + 1
                assert (L * L + 1) % (2 * 5 * 37 * 1637) == 0 and (L * L + 1) // (2 * 5 * 37 * 1637) > 1
                why = "Omega(L^2+1) >= 5 (2*5*37*1637 | L^2+1, quotient > 1)"
            cert.append({"pair": "%s vs %s" % (reps[u], reps[u2]),
                         "gcd_log2": round(math.log2(g), 2) if g > 1 else 0,
                         "v": reps[v] if v in reps else str(v)[:12],
                         "w": reps[w2] if w2 in reps else str(w2)[:12],
                         "dead_by": why})
            dead += 1
    assert dead == 15
    # (3) D_repeat at production
    X1, X2 = 1 << 82, 1 << 70
    Kf = F(K)
    NB = math.isqrt(int(Kf * m / 6)) + 1
    Z = X1 * X2
    Dr = 8 * (F(2 * Z * Z, r) + 1) / (Kf * NB)
    # (4) the retained inequality + counterexample
    slack2 = math.log2(r) - math.log2(2 * X1 * X1 * X2)
    slack3 = math.log2(r) - math.log2(3 * X1 * X1)
    rc = 481153
    Pc = (227 * 167) * (347 * 197)
    Qc = (241 * 73) * (347 * 151)
    mc_ = (rc - 1) // 16
    counter = {
        "r": rc, "prime": is_prime(rc), "r_mod_96": rc % 96,
        "P_minus_Q_is_3470r": Pc - Qc == 3470 * rc, "gcd_PQ": math.gcd(Pc, Qc),
        "all_four_semiprimes_in_H": all(pow(s, mc_, rc) == 1 for s in
                                        [227 * 167, 241 * 73, 347 * 197, 347 * 151]),
        "tiers_disjoint": True, "X1": 400, "X2": 197,
        "3X1_sq_lt_r": 3 * 400 * 400 < rc, "X1_sq_X2_gt_r": 400 * 400 * 197 > rc,
        "note": "a class-matched shared-prime wrapped collision (P - Q = 3470r, "
                "gcd(P,Q) = 347) exists with only 3*X1^2 < r: the shared-prime "
                "exclusion needs the stronger 2*X1^2*X2 < r, which production "
                "retains with 19.857 bits of slack",
    }
    # (5) the two swap bounds
    baseline = F(36 * NB ** 4, m)
    L2v = 2 * int(X2 / (d * math.log(X2)))
    NA = (2 * NB) * 2 * d * 49 // X2 + 2
    T = X1 // 81
    swap_wform = F(501, 1) / (Kf * d)
    swap_interval = F(6 * L2v * (L2v - 1) * 2 * T * T, 4 * d) / baseline
    eps_swap = float(swap_wform)
    # P1's swap charge for the ORIGINAL unpruned total: E_swap/B0 <=
    # 501/(K*d) = 8.69e-4. Route floors: routeb_certificates()
    # (production_routeb).
    beta = swap_wform
    return {
        "constant_structure": {
            "z": hex(z), "L": str(L), "L_squared_plus_L_plus_1_eq_r": True,
            "L_log2": round(math.log2(L), 3),
            "six_representatives": ["%s (order %d)" % s for s in orders],
        },
        "mixed_root_k0_certificate": {"pairs": cert, "all_15_dead": True,
                                      "cutoff_inequalities_needed": "none"},
        "d_repeat": {"formula": "8*(2*Z^2/r + 1)/(K*n)", "Z_log2": math.log2(Z),
                     "value": float(Dr), "log2": round(math.log2(float(Dr)), 3)},
        "size_inequalities": {"2X1sqX2_lt_r": True, "slack_bits": round(slack2, 3),
                              "3X1sq_lt_r": True, "slack_bits_3": round(slack3, 3),
                              "counterexample_r481153": counter},
        "swap_bounds": {
            "wform_fraction_of_baseline": round(eps_swap, 6),
            "wform_premise": "delta_B^w, delta_G^w <= 500",
            "interval_fraction_of_baseline": round(float(swap_interval), 6),
            "interval_premise": "Delta0 <= T^2 over I0 = (X1/2, X1], T <= X1/81",
            "l2_form_values_bound_log2": round(math.log2(6 * (NA // 2) ** 2 * max(1, L2v - 1) // d), 1),
        },
        "residual_premise": {
            "decomposition": "pruned E0/B0 <= 1 + known + Delta0/B0, known <= "
                             "0.159313 CERTIFIED under P1 (sieve 0.156980, inner + "
                             "multi 1.09e-8, shared-large 0.002332911916573675, "
                             "two-large 1.64e-17, repeated 2.671e-19; exact "
                             "rationals in production_routeb); Delta0/B0 = coprime "
                             "discrepancy of the original unpruned tuples, "
                             "UNPROVED, at 53/1250",
            "character_count": m // 6 - 1,
            "character_count_log2": round(math.log2(m / 6 - 1), 3),
            "u_filter_identity": "E + N = N^2/m + (36/m)*sum over the ~m/6 - 1 "
                                 "U-invariant nontrivial characters (verified "
                                 "exactly at both energy points)",
            "conditional_floors": "production_routeb.route_floors: premise "
                                  "Delta0/B0 <= 53/1250 gives route B "
                                  "(P3/P4) >= 0.3900083; route A retired",
            "original_total_caveat": "PRUNED coverage implication only; NOT "
                                     "original total P2(0.057) under any "
                                     "premise; no original-total claim",
            "joint_defect": "UNPROVED",
        },
    }


def structured_lambda_point(K):
    """The cyclotomic toy r = 48^4 - 48^2 + 1 = 5306113 (prime), d = 16,
    X2 = 1024, K = 0.55: L_c = 48^2 - 1 = 2303 with L_c^2 + L_c + 1 = r
    exactly - the production constants in miniature. Energy
    decomposition, family rates, class moments, exact coverage at
    n = 175, ka = 39."""
    r_c = 48 ** 4 - 48 * 48 + 1
    assert r_c == 5306113 and is_prime(r_c)
    Lc = 48 * 48 - 1
    assert Lc * Lc + Lc + 1 == r_c
    out = coverage_prime_point(22, 4, 64, K, 0, energy=True, moments=True, fixed_r=r_c)
    gen = next(g for g in range(2, 10**6) if all(pow(g, (r_c - 1) // q, r_c) != 1 for q in (2, 3)))
    lam = pow(gen, (r_c - 1) // 3, r_c)
    reps = {1, r_c - 1, lam, r_c - lam, lam * lam % r_c, r_c - lam * lam % r_c}
    out["constant_structure"] = {
        "L_c": Lc,
        "L_c_squared_plus_L_c_plus_1_eq_r": True,
        "six_representatives_are_polynomials":
            reps == {1, Lc, Lc + 1, Lc * Lc, Lc * Lc + 1, Lc * (Lc + 1)},
        "representatives": sorted(reps),
    }
    return out


# The collaborator's three certificate modules (shared_large_certificate.py,
# direct_row_sieve_verify.py, norm_dual_joint_verify.py - the routeb pair
# with standalone-output tails removed, as the v30 strata modules) as ONE
# compressed VERBATIM image: decompress, assert in-folder byte-identity
# when present, report the executed texts' sha256s, execute FRESH;
# prime_count_inputs.py ships in-folder (stdlib, runnable) for the Dusart
# transcription spot-checks; the hosted runner stages only this file.
_MODULES_B64 = "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"


import base64
import zlib
_DATA_B64 = "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"
_DATA = json.loads(zlib.decompress(base64.b64decode(_DATA_B64)))

def routeb_certificates():
    """The v31 route-B certificates: the three embedded VERBATIM collaborator
    modules execute FRESH (hosted stages only this file; in-folder
    identity asserted when present, sha256s reported). Under P1 and the
    .57/.64 endpoints, NO H2, no bound on the Selberg remainder (the
    sieve REPLACES it): shared-small 0.156980, inner+multi 1.09e-8,
    shared-large 0.002332911916573675, two-large 1.64e-17, repeated
    2.671e-19, pruning 4.79e-14, ledger 0.159313, route-B floor at the
    NAMED coprime premise 53/1250 (UNPROVED; see coprime_measurement)."""
    import contextlib
    import io
    import tempfile
    from pathlib import Path
    root = Path(__file__).resolve().parent
    import base64 as _b64mod
    import hashlib as _sha
    import lzma as _lzma
    _SEP = "\n### MODULE BOUNDARY ###\n"
    _PARTS = _lzma.decompress(_b64mod.b64decode(_MODULES_B64)).decode("utf-8").split(_SEP)
    assert len(_PARTS) == 3
    _SCRATCH = tempfile.mkdtemp(prefix="cheon-v31-")
    def fresh(name, part):
        # Hosted stages only this file; an in-folder module must match
        # the embedded image. Modules write a .json beside __file__, so
        # __file__ is a scratch copy.
        digest = _sha.sha256(part.encode("utf-8")).hexdigest()
        if (root / name).exists():
            assert (root / name).read_text(encoding="utf-8") == part, name
        namespace = {"__name__": name[:-3],
                     "__file__": str(Path(_SCRATCH) / name)}
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(part, name, "exec"), namespace)
        namespace["_module_sha256"] = digest
        return namespace
    large = fresh("shared_large_certificate.py", _PARTS[0])
    row = fresh("direct_row_sieve_verify.py", _PARTS[1])
    geom = fresh("norm_dual_joint_verify.py", _PARTS[2])
    rr = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    dd = 1 << 20
    assert rr == r and (rr - 1) % dd == 0
    mm = (rr - 1) // dd
    assert mm == m
    nn = math.isqrt(11 * mm // 120)
    if 120 * nn * nn < 11 * mm:
        nn += 1
    assert 120 * (nn - 1) ** 2 < 11 * mm <= 120 * nn * nn
    NA_ = 196 * nn * dd // (1 << 70) + 2
    MM = (NA_ + 1) // 2
    NN = 6 * nn * nn
    B0 = F(36 * nn ** 4, mm)
    KK = F(NN, mm)
    X1_, X2_ = 1 << 82, 1 << 70
    ZZ = X1_ * X2_
    for ns in [large, row, geom]:
        assert ns["r"] == rr and ns["m"] == mm
        assert ns["n"] == nn and ns["M"] == MM
    assert large["baseline"] == B0 and geom["base"] == B0 and row["B0"] == B0
    Lcal = large["ratio"]
    ss = F(row["best"]["total_exact"])
    assert row["best"]["pmax"] == 31 and row["out"]["safe_bound"] == "0.156980"
    assert ss < F("0.156980")
    inner = F(geom["inner_tight"], 1) / B0
    multi = F(geom["twoL"] + geom["oneLoneS_tight"] + geom["twoS_tight"], 1) / B0
    two_large = F(96 * MM * MM * (2 * X2_ ** 4 // rr + 1), 1) / B0
    assert two_large == large["Etwo_large"]
    repeated = F(48 * nn * (2 * ZZ * ZZ // rr + 1), 1) / B0
    assert repeated == large["Erepeat_ratio"]
    loss = F(6 * MM * nn, mm)
    assert loss == large["pruning_mass_over_m"]
    # the multi cover already contains the repeats; the tiny repeated
    # charge is kept for safety
    known = ss + inner + multi + Lcal + repeated
    assert known <= F("0.1593122871411313")
    A = F(159313, 1000000)
    assert known < A
    kappa = KK - loss
    assert F(".549") < kappa < KK < F(".551")
    def gain3and4(x):
        return F(19, 20) * x ** 3 / 6 - F(6, 5) * x ** 4 / 24
    assert KK < F(19, 8)
    Bterm = kappa + gain3and4(kappa)
    premise = F(53, 1250)
    goal = F(".39")
    floorB = Bterm - (1 + A + premise) * KK * KK / 2
    floorB_x = Bterm - (1 + known + premise) * KK * KK / 2
    floorA0 = kappa - (1 + A) * KK * KK / 2
    max_d = 2 * (Bterm - goal) / (KK * KK) - 1 - A
    assert floorB > goal and floorB_x > goal and floorA0 < goal
    beta_ = F(501) * mm / (6 * dd * nn * nn)
    assert beta_ < F(".000869")
    return {
        "status": "ROUTE_B_COMPONENTS_PASS_COPRIME_PREMISE_UNPROVED",
        "external_theorem": "Dusart arXiv:1002.0442v1 Thm 6.9 eq 6.6 (stated in "
                            "prime-count-dependency.txt); stated, not reproved; every "
                            "consuming figure recomputed locally",
        "same_side_one_small_exact": str(ss),
        "same_side_one_small_safe_bound": "0.156980",
        "inner_bound_exact": str(inner),
        "all_multicommon_bound_exact": str(multi),
        "shared_large_exact": str(Lcal),
        "shared_large_ceiling": "0.002332911916573675",
        "two_large_correction_exact": str(two_large),
        "repeated_coordinate_exact": str(repeated),
        "pruning_loss_over_m_exact": str(loss),
        "module_sha256": {"shared_large_certificate.py": large["_module_sha256"],
                          "direct_row_sieve_verify.py": row["_module_sha256"],
                          "norm_dual_joint_verify.py": geom["_module_sha256"]},
        "known_pruned_excess_exact": str(known),
        "ledger_bound_exact": str(A),
        "conditional_inequality": "pruned E0/B0 <= 1 + known + Delta0/B0; "
                                  "Delta0/B0 <= 53/1250 NAMED, UNPROVED",
        "route_floors": {
            "kappa_exact": str(kappa),
            "premise": "53/1250",
            "routeB_floor_at_premise_exact": str(floorB),
            "routeB_floor_at_premise_decimal": str(float(floorB)),
            "routeB_floor_premise_exact_excess_known": str(floorB_x),
            "max_coprime_for_0.39_exact": str(max_d),
            "routeA_floor_at_zero_defect_exact": str(floorA0),
            "routeA_retired": True,
            "swap_beta_exact": str(beta_),
            "coprime_premise": "UNPROVED; supported by the "
                               "coprime_measurement ladder",
        },
    }


def coprime_measurement():
    """The v31 coprime ladder (record: measurement_campaign_report.md): seven
    EXACT controls (d = 2^4..2^20), three whole-row HT targets at the
    production d = 2^20 with the ACTUAL masks and the original unpruned
    full-critical prefixes (design-unbiasedness by exhaustive oracles in
    the collaborator's bundle), one dense complement. Every consuming
    fraction RECOMPUTED from embedded raw counts. SEs descriptive only;
    the 35-bit X2/d, ~182.3-bit field-size and tier/block-width gaps stay
    explicit. Premise consumed: the NAMED fixed-instance extrapolation
    E_coprime <= 36*M0/m + (53/1250)*B0."""
    allowance = F(53, 1250)
    _L = _DATA["ladder"]
    controls, ht = _L["controls"], _L["ht"]
    out_c = []
    for name, rr, dd, nn, E, M0, rep, x2r in controls:
        mm = (rr - 1) // dd
        val = F(mm * E - 36 * M0, 36 * nn ** 4)
        assert val == F(rep), name
        assert val < allowance, name
        out_c.append({"case": name, "x2_over_d": x2r, "d": dd,
                      "delta0_over_b0_exact": str(val),
                      "theta_exact": str(F(M0, nn ** 4)),
                      "below_allowance": True})
    out_h = []
    for name, ratio, nfull, seeds, rep_mean, rep_se in [tuple(x) for x in ht]:
        frs = [F(s) for s in seeds]
        mean = sum(frs) / len(frs)
        assert mean == F(rep_mean), name
        var = sum((f - mean) ** 2 for f in frs) / (len(frs) - 1)
        se = math.sqrt(float(var) / len(frs))
        assert abs(se - float(rep_se)) < 1e-9, name
        margin = float((allowance - mean)) / se
        assert mean < allowance
        out_h.append({"case": name, "x2_over_d": ratio, "n_full": nfull,
                      "replicates": len(frs), "mean_ht_delta0_over_b0_exact": str(mean),
                      "empirical_se": str(se),
                      "margin_in_SEs_below_allowance": str(margin)})
    dn = _L["dense"]
    dmm = (dn["r"] - 1) // dn["d"]
    dval = F(dmm * dn["E_coprime"] - 36 * dn["M0"], 36 * dn["n"] ** 4)
    assert dval == F(dn["delta0_over_b0_exact"])
    assert dval < allowance
    return {
        "status": "LADDER_RECOMPUTED_FROM_RAW_RECORDS_PREMISE_UNPROVED",
        "premise": "fixed-instance extrapolation E_coprime <= 36*M0/m + "
                   "(53/1250)*B0 - NAMED, ladder-supported, NOT a theorem",
        "exact_controls": out_c,
        "ht_targets": out_h,
        "dense_lowload_exact": {"r": dn["r"], "d": dn["d"], "n": dn["n"],
                                "delta0_over_b0_exact": str(dval)},
        "gaps_disclosure": "35-bit X2/d gap (2^15 vs 2^50), ~182.3-bit "
                           "field-size gap, tier/block-width gaps; SEs "
                           "descriptive only; no confidence interval, no "
                           "uniform bound",
        "record_note": "raw records embedded as exact integers; the full "
                       "closure (registrations, hashes, variance batches) is "
                       "the collaborator's bundle and the in-folder report",
    }


def p34_measurement():
    """The v32 P3/P4 ladder (record: p34-campaign-report.txt): swap-pruned
    RATIO occupancy moments normalized by the actual retained mean
    kappa = T0/h. Seven exact controls (recomputed from embedded
    histograms) + three 64-stratum output-character HT estimates at the
    production d = 2^20 with the exact original prefixes; SEs are design
    sampling uncertainty only; P3/P4 stay premises with ladder support."""
    import math as _math3
    _L = _DATA["p34"]
    out_c = []
    for name, hh, hist in _L["controls"]:
        hist = {int(k): v for k, v in hist.items()}
        f1 = F(sum(k * v for k, v in hist.items()), hh)
        f3 = F(sum(k * (k - 1) * (k - 2) * v for k, v in hist.items()), hh) / f1 ** 3
        f4 = F(sum(k * (k - 1) * (k - 2) * (k - 3) * v for k, v in hist.items()), hh) / f1 ** 4
        weak = name.endswith("low-K-screen")
        assert float(f3) >= 0.95 and (weak or float(f4) <= 1.20), name
        out_c.append({"case": name, "k0": float(f1), "F3_over_kappa3": str(f3),
                      "F4_over_kappa4": str(f4),
                      "F4_note": "weak: ~0.73 expected null quadruple events" if weak else ""})
    out_t = []
    for name, k0e, f3e, f4e, c33, c44, j34e, jv34e in _L["targets"]:
        k0 = F(k0e)
        f3, f4 = F(f3e) / k0 ** 3, F(f4e) / k0 ** 4
        se3 = _math3.sqrt(float(F(c33) / k0 ** 6))
        se4 = _math3.sqrt(float(F(c44) / k0 ** 8))
        j34, vj = F(j34e), F(jv34e)
        bound = F(95, 100) * k0 ** 3 / 6 - F(6, 5) * k0 ** 4 / 24
        out_t.append({"case": name, "F3_over_kappa3": str(f3), "F3_design_se": str(se3),
                      "F4_over_kappa4": str(f4), "F4_design_se": str(se4),
                      "P3_margin_in_SEs": str((float(f3) - 0.95) / se3),
                      "P4_margin_in_SEs": str((1.20 - float(f4)) / se4),
                      "joint_gain_J34_exact": str(j34), "J34_variance_exact": str(vj),
                      "J34_margin_in_SEs": str(float(j34 - bound) / _math3.sqrt(float(vj)))})
    return {
        "status": "P34_LADDER_RECOMPUTED_FROM_RAW_RECORDS_PREMISES_RETAINED",
        "normalizer": "the actual retained mean kappa = T0/h after the "
                      "declared U-consistent swap pruning (ratio occupancy)",
        "exact_controls": out_c,
        "ht_targets": out_t,
        "null_reference": "the hypergeometric Cartesian null predicts "
                          "normalized moments 1 - ~8.9e-35 at production; a "
                          "model reference, not a theorem",
        "gaps_disclosure": "SEs are design sampling uncertainty at the "
                           "fixed targets only; the 35-bit X2/d gap, "
                           "geometry gaps and PRNG idealization remain; "
                           "P3/P4 stay premises; Delta0 <= 53/1250 stays "
                           "the separate named coprime premise",
    }

cfg = json.load(sys.stdin)
K = cfg["K"]
print(json.dumps({
    "ledger": ledger(K),
    "production_certificates": production_certificates(K),
    "production_routeb": routeb_certificates(),
    "coprime_measurement": coprime_measurement(),
    "p34_measurement": p34_measurement(),
    "hash_probe": hash_probe_check({"bits": 28, "n": 1 << 14,
                                    "shards": 16, "seed": cfg["seed"]}),
    "coverage_prime_sweep": [
        coverage_prime_point(25, 4, 64, float(K), cfg["seed"] + 1, baselines=True, nulls=8, moments=True),
        coverage_prime_point(25, 4, 64, float(K), cfg["seed"] + 1, list_filter="integers"),
        coverage_prime_point(cfg["cov_bits"] - 3, cfg["cov_logd"], 16, float(K), cfg["seed"] + 2,
                             nulls=6, energy=True, moments=True),
    ],
    "structured_lambda_point": structured_lambda_point(float(K)),
    "toy_attack_fixed_order": toy_attack({"bits": cfg["atk_bits"],
                                          "logd": cfg["atk_logd"],
                                          "trials": 16, "K": K,
                                          "seed": cfg["seed"], "gmax": cfg["gmax"],
                                          "fixed": True}),
    "toy_attack": toy_attack({"bits": cfg["atk_bits"], "logd": cfg["atk_logd"], "trials": cfg["trials"],
                              "K": K, "seed": cfg["seed"], "gmax": cfg["gmax"]}),
}, indent=1))

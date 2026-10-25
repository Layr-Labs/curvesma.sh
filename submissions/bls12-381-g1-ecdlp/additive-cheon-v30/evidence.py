"""Supporting evidence for proof.md (stdlib only; JSON on stdin, JSON on stdout).

additive-cheon-v30: the algorithm, walks, ledger and every cost figure are
IDENTICAL to the v28 lineage (and through it to the declared dependency
cbadde77, restated in full in proof.md; every identity, itemisation and
bound below is this package's own, re-verified by these runs). What v29
restructured is the SUCCESS argument, and what v30 adds is the residual's
CERTIFIED STRATA (item 2b below): every structured stratum proved at
production under P1 with no H2, leaving one named joint-defect bound.

1. ledger      - unchanged: every cost term of proof.md section 4 and the exact score.
2. production_certificates - the fixed target's constant structure (L = x^2 - 1,
                 L^2 + L + 1 = r exactly; the six U-representatives are the
                 polynomials {1, L, L+1, L^2, L^2+1, L(L+1)} with orders
                 1, 3, 6, 3, 6, 2), the cutoff-free mixed-root k=0 certificate
                 (all 15 unordered mixed pairs die by the Omega-equality
                 necessary condition: eight have a reduced coefficient 1 against
                 >1, five have a reduced coefficient divisible by 2^31, two die
                 by Omega(L^2+1) >= 5 from 2*5*37*1637 | L^2+1 with quotient >1),
                 the repeated-coordinate bound D_repeat = 8(2Z^2/r + 1)/(K n),
                 the retained size inequality 2*X1^2*X2 < r with its slack, the
                 r=481153 class-matched counterexample showing the weaker
                 3*X1^2 < r does NOT exclude shared-prime wrapped collisions,
                 and the two production swap bounds.
2b. production_strata - the residual's stratum certificates, executed FRESH
                 (runpy) from the three in-folder source modules
                 (prime_count_inputs.py, shared_large_certificate.py,
                 selberg_certificate.py): the exactly-one-shared-large
                 lattice theorem (0.002332911916573675 of the birthday
                 baseline, the mod-30 wheel's 16/225, the actual P1 class
                 masses, Dusart-certified block endpoints a_B <= 0.57*X1,
                 a_G <= 0.64*X1), the shared-small Selberg main
                 (0.038976659226098 with the integer-certified hyperbola
                 denominator G(2^35) > 322.7629594) and small-prime
                 exception (0.002641842), the two-large-only (1.633e-17)
                 and repeated-coordinate (2.671e-19) corrections, the
                 coprime main to 3.492e-13, the coverage-preserving
                 pruning (loss 4.79e-14, no H2), the ledger bound
                 alpha = 0.043951414, and the conditional route floors as
                 exact functions of the UNPROVED joint defect
                 t = (Delta_0 + R_SS)/B0. The only external analytic input
                 is Dusart arXiv:1002.0442v1 Thm 6.9 eq (6.6), stated in
                 prime-count-dependency.txt; stated, not reproved; every
                 consuming figure recomputed locally.
3. coverage    - as v28 (calibrated nulls, the matched integer-list control,
                 the smallest-integers baseline), extended with the exact
                 distinct-index energy decomposition E = 6*sum s(s-1) + E_{k!=0}
                 (the diagonal belongs to E+N), the per-ratio-order family
                 rates, and the class-marginal second moments delta^w and
                 Delta0/T^2 that the v29 swap bound and its interval
                 corroboration consume.
4. structured_lambda_point - the cyclotomic toy r = 48^4 - 48^2 + 1
                 (L_c = 48^2 - 1 with L_c^2 + L_c + 1 = r exactly: the
                 production constant structure in miniature, three small and
                 three large representatives), full energy decomposition and
                 family rates at d = 16, X2 = 1024, n = 175, ka = 39.
5. toy_attack / hash_probe - unchanged.

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
    # CANON: published package itemised "p - y (~7), 6-limb compare (~12),
    #   select (6): ~25"; v2 charges 30 (v1 charged 40).
    # HASH_OP (v3, complete per-transfer accounting): the table of 2^116.7
    #   slots x 72-byte records at load <= 1/2 cannot sit in one 2^63-byte
    #   address space, so it is sharded across <= 2^60 nodes and reached
    #   through a radix-1024 switch fabric: ceil(60/10) = 6 relay stages.
    #   Endpoint work: compute shard and slot from key bits 4; append the
    #   9-word record to the outbound buffer 14; read at the receiver 9;
    #   linear probing at load <= 1/2 (Knuth unsuccessful-search bound
    #   (1+1/(1-a)^2)/2 = 2.5 expected probes of 6 instructions) 15; store
    #   the record (inserts) 9; loop and control 10 - total 61. EACH relay
    #   stage: load and store each of the 9 record words 18; destination
    #   selection (extract the 10-bit next-stage field, port-table load) 4;
    #   queue head/tail updates 4; loop and branch control 3 - 29 per
    #   stage, 6 stages: 174. TOTAL 61 + 174 = 235, charged 240.
    # WALK: itemised - read b and the gap 4; overflow test 2; address T[gap]
    #   3; load T[gap] 12; load and store points 24; batch-inversion
    #   buffers beyond the 6M 24; loop 5 - ~75; charged 100.
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
    # proven interval lower bound: pi(X1) >= X1/ln X1 (Rosser-Schoenfeld, X1 >= 17)
    # and pi(X1/2) <= 1.25506 * (X1/2)/ln(X1/2) (Rosser-Schoenfeld, x > 1)
    PI_X1 = F(126, 100) * X1 / F(int(82 * math.log(2) * 10**6), 10**6)   # pi(X1) upper estimate (cost charge)
    # proven interval lower bound (Rosser-Schoenfeld): pi(X1) >= X1/ln X1 (x >= 17) and
    # pi(X1/2) <= 1.25506 * (X1/2)/ln(X1/2) (x > 1), so (X1/2, X1] holds at least
    P_AVAIL = (F(X1) / F(int(82 * math.log(2) * 10**6), 10**6)
               - F(125506, 100000) * F(X1, 2) / F(int(81 * math.log(2) * 10**6), 10**6))
    N_A = (N_B + N_G) * 2 * d * 49 // X2 + 2   # primes a used; the b-lists hold only the class PRIMES,
                                               # each >= X2/(2*d*49) (H2; ln(2^70) = 48.5 < 49, 1/2 slack)
    # shared-class-prefix bound: a same-class (baby, giant) prime pair with a shared
    # consumed prefix of t elements collapses its t diagonal pairs (x = a*b, y = a'*b)
    # to the single value u*a/a' for each u - losing 6*(t-1) values; the shared prefix
    # is at most L2 = 2x the mean class size (H2 concentration, factor-2 slack)
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
    """Toy evidence for the section-4 hashing model with keys from the
    STRUCTURED production-style walk: for primes a of a production-style
    block, the b's of the matching character class are taken in increasing
    order (the production walk order) and the keys are the canonical-y of
    the resulting points [b]Q_a. The keys are inserted into a sharded
    open-addressing table at load 1/2 (shard field disjoint from the slot
    field, as in production); probe counts and per-shard loads are measured
    against the standard uniform-hashing bounds (Knuth: unsuccessful search
    at load 1/2 expects (1+1/(1-a)^2)/2 = 2.5 probes; shard loads are
    Binomial(n, 1/S))."""
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
        # walk the class list in increasing order at ONE addition per element
        # (with overflow chunks for gaps above gmax), as in production
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
    """Exact |U*B*G^-1|/m under the PRODUCTION construction: the b's are the
    PRIMES of the chi-class lists (in increasing order), the a's are the
    first ka = ceil(N_A/2) primes of (X1/2, X1] (baby block) and the next ka
    (giant block), ascending. With list_filter="primes" (production) every
    element is a product a*b of two primes from disjoint ranges (a > X2 >= b);
    with list_filter="integers" the SAME r, blocks and n use the full integer
    class lists, isolating the list construction alone.

    v29 extensions (energy=True): the exact distinct-index energy
    decomposition. E counts ordered pairs of distinct triples (u,x,y) with
    equal coverage values u*x*y^-1 mod r; the diagonal (x = x', y = y')
    belongs to E + N, not to E, and the k = 0 integer identities split as
        E = 6*sum_{same-class row pairs} s(s-1)  +  mixed  +  E_{k != 0},
    with s = min(t_a, t_a') the shared consumed prefix of a same-class
    (baby row a, giant row a') pair. The mixed-root count is computed
    exactly (it is 0 at every measured point) using the necessary condition
    Omega(u) = Omega(u') from unique factorisation over 4-prime products.
    Per-ratio-order family rates group the 30 ordered mixed (u,u') pairs by
    the order of rho = u'/u (six pairs per rho value; null per family n^4/m).
    moments=True adds the class-marginal second moments consumed by the v29
    swap bound: delta_s^w = (1/n^2) sum_j |sum_a t_a chi(a)^j|^2 over the
    emitted row masses, and Delta0/T^2 for the full interval I0."""
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
        # ---- the exact distinct-index energy decomposition (v29) ----
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
        # per-ratio family totals (each rho value carries exactly 6 ordered pairs)
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
        # mixed-root k=0 integer identities, exactly, via Omega(u) = Omega(u')
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
        # total-energy excess (what the P2 box prices via the U-filter identity)
        # and the pruned factorial moments: removing all but one representative
        # of each same-root swap value preserves coverage exactly; fourth-order
        # Bonferroni on the pruned multiplicities is the model-free alternative
        # to the law (proof.md 3.7). The sector self fourth-moments and the
        # mixed moment are computed by exact U-coset orthogonality, and the
        # mixed value is checked against the U-filter identity exactly.
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
        # the friend's diagonal lower bound: any proof via separate self
        # fourth-moment bounds plus Cauchy-Schwarz needs eps >= h/n^2*(1-1/n)-1
        # ~ 1/K - 1 (0.818 at K = 0.55) - the route is blocked, verified per run
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
        # ---- class-marginal second moments (v29 swap premise calibration) ----
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
            # the law's null, calibrated: nulls independent random set-pairs,
            # counted by the same routine at the same parameters
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
    """The v29 production certificates: all arithmetic, no sampling.
    (1) the constant structure L = x^2 - 1, L^2 + L + 1 = r exactly, the six
    U-representatives as polynomials in L with their multiplicative orders;
    (2) the cutoff-free mixed-root k=0 certificate: a mixed integer identity
    u*P = u'*Q with P, Q four-prime products forces, after reducing
    u/u' = v/w (coprime), v | Q, w | P and Omega(v) = Omega(w) <= 4; all 15
    unordered mixed pairs of the six constants fail this - eight with a
    reduced coefficient 1 against >1, five with a reduced coefficient
    divisible by 2^31 (Omega >= 31), two by Omega(L^2+1) >= 5 (divisors
    2*5*37*1637 with quotient > 1), so no size hypothesis is needed;
    (3) the repeated-coordinate bound D_repeat = 8(2Z^2/r + 1)/(K n) with
    Z = X1*X2, the amplification mechanism of the small-d toys, at the
    production parameters;
    (4) the retained size inequality 2*X1^2*X2 < r (19.857 bits of slack) and
    the class-matched counterexample r = 481153 showing the weaker
    3*X1^2 < r does NOT exclude shared-prime wrapped collisions;
    (5) the two production swap bounds: the emitted-mass form
    E_swap/(N^2/m) <= (1 + sqrt(delta_B^w delta_G^w))/(K d) (premise
    delta^w <= 500 gives 0.0869% of baseline) and the full-interval form
    (premise Delta0 <= T^2 gives 0.1452%); the v28 L2-form is kept as
    corroboration by the ledger."""
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
    # (4) the retained inequality and the counterexample
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
    # P1's swap charge beta for the ORIGINAL unpruned total (the comparison
    # the pruned budgets must not be conflated with): E_swap/B0 <= 501/(K*d)
    # = 8.69e-4. The route floors are computed exactly, as functions of the
    # unproved joint defect t, by strata_certificates() below (the
    # production_strata block of the output).
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
            "decomposition": "the pruned distinct-index energy obeys E0/B0 <= 1 + alpha + t "
                             "with alpha = 0.043951414 CERTIFIED at production under P1 "
                             "(shared-large 0.002332911916573675, Selberg main "
                             "0.03897665922609721, small-prime exception 0.00264184166416664, "
                             "two-large 1.64e-17, repeated-coordinate 2.671e-19 - the exact "
                             "rationals and certificates are in production_strata) and "
                             "t = (Delta_0 + R_SS)/B0 UNPROVED: the coprime modular "
                             "discrepancy plus the signed Selberg remainder",
            "character_count": m // 6 - 1,
            "character_count_log2": round(math.log2(m / 6 - 1), 3),
            "u_filter_identity": "E + N = N^2/m + (36/m)*sum over the ~m/6 - 1 U-invariant "
                                 "nontrivial characters (verified exactly at both energy points)",
            "conditional_floors": "production_strata.route_floors: t <= 0.0138998 gives "
                                  "route A > 0.39 (the claim); t <= 0.157816 gives route B "
                                  "> 0.39 additionally under P3/P4",
            "original_total_caveat": "the pruned-route budgets must not be described as "
                                     "proving original total P2(0.057); at the stricter "
                                     "t = 0.0121794 the original unpruned total follows "
                                     "conditionally by adding P1's swap beta < 0.000869",
            "joint_defect": "UNPROVED",
        },
    }


def structured_lambda_point(K):
    """The cyclotomic toy: r = 48^4 - 48^2 + 1 = 5306113 (prime), d = 16,
    X2 = 1024, K = 0.55. Its small cube root L_c = 48^2 - 1 = 2303 satisfies
    L_c^2 + L_c + 1 = r EXACTLY, so the six U-representatives are the
    polynomials {1, L_c, L_c+1, L_c^2, L_c^2+1, L_c^2+L_c} - the production
    constant structure in miniature (three small, three large). Full energy
    decomposition, per-order family rates, class moments and the exact
    coverage count at n = 175, ka = 39, rows 15/19."""
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


def strata_certificates():
    """The v30 production strata certificates. The three source modules
    (prime_count_inputs.py, shared_large_certificate.py,
    selberg_certificate.py) ship in this folder and are re-executed FRESH
    here (runpy), exactly as the standalone checker does; no prior-run
    result is ever an input. The only external analytic input is Dusart
    (arXiv:1002.0442v1, Theorem 6.9, equation 6.6), stated with its
    validity domains in prime-count-dependency.txt; the lemma is stated,
    not reproved, and every consuming numerical figure is recomputed
    locally, with independent sieved transcription spot-checks. No H2
    clause enters any bound or the pruning. Output: the exact shared-large,
    Selberg-main and small-prime-exception rationals with directed decimal
    ceilings, the two-large and repeated-coordinate corrections, the
    coverage-preserving pruning loss, the ledger bound alpha, and the
    conditional route floors as exact functions of the UNPROVED joint
    defect t = (Delta_0 + R_SS)/B0."""
    import contextlib
    import io
    import runpy
    from pathlib import Path
    root = Path(__file__).resolve().parent
    def fresh(name):
        with contextlib.redirect_stdout(io.StringIO()):
            return runpy.run_path(str(root / name))
    pc = fresh("prime_count_inputs.py")
    prime_inputs = pc["finite_consuming_checks"]()
    prime_spotchecks = pc["small_transcription_spotchecks"]()
    large = fresh("shared_large_certificate.py")
    small = fresh("selberg_certificate.py")
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
    for ns in [large, small]:
        assert ns["r"] == rr and ns["d"] == dd and ns["m"] == mm
        assert ns["n"] == nn and ns["M"] == MM and ns["baseline"] == B0
    assert prime_inputs["M"] == MM
    Lcal = large["ratio"]
    main_, exception_, Glo = small["main"], small["exception"], small["Dlo"]
    two_large = F(96 * MM * MM * (2 * X2_ ** 4 // rr + 1), 1) / B0
    assert two_large == large["Etwo_large"]
    repeated = F(48 * nn * (2 * ZZ * ZZ // rr + 1), 1) / B0
    assert repeated == large["Erepeat_ratio"]
    loss = F(6 * MM * nn, mm)
    assert loss == large["pruning_mass_over_m"]
    known = Lcal + main_ + exception_ + two_large + repeated
    # the package's conservative 9-place ledger bound: directed ceiling
    # plus one final decimal unit
    scale_alpha = 10 ** 9
    alpha_units = (known.numerator * scale_alpha + known.denominator - 1) // known.denominator + 1
    alpha = F(alpha_units, scale_alpha)
    assert alpha == F(".043951414")
    assert known < alpha
    assert F(".002332911916573674") < Lcal < F(".002332911916573675")
    assert main_ < F(".03897665922609721")
    assert exception_ < F(".00264184166416664")
    assert Glo > F("322.7629594")
    kappa = KK - loss
    assert F(".549") < kappa < KK < F(".551")
    def gain3and4(x):
        return F(19, 20) * x ** 3 / 6 - F(6, 5) * x ** 4 / 24
    # x + gain3and4(x) is increasing on [0, KK]: derivative 1 + x^2(19/40 - x/5) > 0
    assert KK < F(19, 8)
    def floor_A(t):
        return kappa - (1 + alpha + F(t)) * KK * KK / 2
    def floor_B(t):
        return kappa + gain3and4(kappa) - (1 + alpha + F(t)) * KK * KK / 2
    goal = F(".39")
    max_A = 2 * (kappa - goal) / (KK * KK) - 1 - alpha
    max_B = 2 * (kappa + gain3and4(kappa) - goal) / (KK * KK) - 1 - alpha
    max_B_398 = 2 * (kappa + gain3and4(kappa) - F(".398")) / (KK * KK) - 1 - alpha
    safe_A, safe_B = F(".0138998"), F(".157816")
    assert floor_A(safe_A) > goal and floor_B(safe_B) > goal
    assert floor_B(F(".158")) < goal
    assert floor_B(safe_B) < F(".398")
    original_t = F(".0121794")
    beta_ = F(501) * mm / (6 * dd * nn * nn)
    assert beta_ < F(".000869")
    assert alpha + original_t + beta_ < F(".057")
    assert floor_A(original_t) > F(".390260")
    return {
        "status": "COMPONENT_CERTIFICATES_PASS_JOINT_DEFECT_UNPROVED",
        "external_theorem": "Dusart arXiv:1002.0442v1 Thm 6.9 eq 6.6 (stated in "
                            "prime-count-dependency.txt); stated, not reproved; every "
                            "consuming figure recomputed locally",
        "prime_count_inputs": prime_inputs,
        "prime_count_spotchecks": prime_spotchecks,
        "shared_large_exact": str(Lcal),
        "shared_large_ceiling": "0.002332911916573675",
        "selberg_main_exact": str(main_),
        "selberg_main_ceiling": "0.03897665922609721",
        "small_prime_exception_exact": str(exception_),
        "small_prime_exception_ceiling": "0.00264184166416664",
        "two_large_correction_exact": str(two_large),
        "repeated_coordinate_exact": str(repeated),
        "pruning_loss_over_m_exact": str(loss),
        "all_known_exact": str(known),
        "alpha_ledger_bound_exact": str(alpha),
        "conditional_inequality": "pruned E0/B0 <= 1 + alpha + t, t = (Delta_0 + R_SS)/B0 UNPROVED",
        "route_floors": {
            "kappa_exact": str(kappa),
            "routeA_floor_at_safe_gap_exact": str(floor_A(safe_A)),
            "routeA_safe_gap": "0.0138998",
            "routeA_max_gap_exact": str(max_A),
            "routeB_floor_at_safe_gap_exact": str(floor_B(safe_B)),
            "routeB_safe_gap": "0.157816",
            "routeB_max_gap_exact": str(max_B),
            "routeB_gap_for_0.398_exact": str(max_B_398),
            "at_original_target_t": "0.0121794: route A > 0.390260; with P1's swap beta "
                                     "< 0.000869 the ORIGINAL unpruned total P2(0.057) "
                                     "follows conditionally",
            "swap_beta_exact": str(beta_),
            "joint_defect": "UNPROVED",
        },
    }

cfg = json.load(sys.stdin)
K = cfg["K"]
print(json.dumps({
    "ledger": ledger(K),
    "production_certificates": production_certificates(K),
    "production_strata": strata_certificates(),
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

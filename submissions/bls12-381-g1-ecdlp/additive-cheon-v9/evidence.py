"""Supporting evidence for proof.md (stdlib only; JSON on stdin, JSON on stdout). Three parts.
This is additive-cheon-v3: the construction of the declared dependency
cbadde77 (provenance only; every identity, itemisation and bound below is
this package's own, re-verified by these runs) with tightened per-element
instruction budgets, a complete per-stage hash accounting (relay word moves
PLUS destination selection, queue updates and control), a floor-rounded
success claim, and a measured probe/shard-distribution check for the
linear-probing model.

1. ledger      - recomputes every cost term of proof.md section 4 and the exact score (arithmetic check only).
2. coverage    - toy evidence for H1/H2: exact |U*B*G^-1|/m for random sets, smallest integers in H and the
                 production product construction a*b (chi(a)+chi(b)=0), vs 1-exp(-K) and the gcd-corrected
                 prediction; class-list spread.
3. toy_attack  - the full algorithm on a toy j=0 curve (walks at one addition per element, keyed on
                 canonical y = min(y, p-y), which is invariant under U = {±1, ±λ, ±λ²}; collision, d-th root, mu_d search), success rate vs prediction, scalars checked.
Toy results are evidence at their scale only, not production cost or success proofs.
Input: {"K", "cov_bits", "cov_logd", "atk_bits", "atk_logd", "trials", "gmax", "seed"}.
"""
import json
import math
import random
import sys
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

    X1, X2, GMAX = 1 << 76, 1 << 70, 1 << 30
    SCALAR = 1000          # units: one non-batched scalar multiplication by a <= 2^76-bit scalar (proof.md 4.3)
    SCALAR_FULL = 3000     # units: one scalar multiplication by a 255-bit scalar
    PI_X2 = F(126, 100) * X2 / F(int(70 * math.log(2) * 10**6), 10**6)    # pi(x) < 1.26 x / ln x
    PI_X1 = F(126, 100) * X1 / F(int(76 * math.log(2) * 10**6), 10**6)
    N_A = (N_B + N_G) * 2 * d // X2 + 2                    # primes a used, heuristic H2 (lists >= X2/(2d))
    assert N_A < PI_X1 / 4                                  # enough primes in (X1/2, X1] (>= pi(X1)/2 - pi(X1/2)...)

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
        "structural_loss_bound_corroboration": 2 * 6 / math.pi ** 2
            * (math.log(1 << 70) + 2 * 0.5772156649 - 1) / d,
        "success_claim": math.floor(success * 1000) / 1000 - 0.004,
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


def coverage(cfg):
    bits, logd, K, seed = cfg["bits"], cfg["logd"], float(cfg["K"]), cfg["seed"]
    fixed, ab_below_r = cfg.get("fixed", False), cfg.get("ab_below_r", False)
    products_only = cfg.get("products_only", False)
    rng = random.Random(seed)
    d = 1 << logd
    while True:
        r = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * d) * (6 * d) + 1
        if is_prime(r) and ((r - 1) // d) % 6 == 0:
            break
    m = (r - 1) // d
    gen = next(g for g in range(2, 10**6) if all(pow(g, (r - 1) // q, r) != 1 for q in (2, 3)))
    lam = pow(gen, (r - 1) // 3, r)
    U = [1, r - 1, lam, r - lam, lam * lam % r, r - lam * lam % r]
    assert all(pow(u, m, r) == 1 for u in U)
    zeta = pow(gen, m, r)                       # generator of mu_d
    zlog = {pow(zeta, i, r): i for i in range(d)}
    chi = lambda x: zlog[pow(x, m, r)]          # character of order d, kernel H
    n = math.isqrt(int(K * m / 6)) + 1

    def coverage(B, G):
        Ginv = [pow(g, -1, r) for g in G]
        bitmap = bytearray(r // 8 + 1)          # exact |U*B*G^-1| without a large set
        for u in U:
            for b in B:
                ub = u * b % r
                for gi in Ginv:
                    v = ub * gi % r
                    bitmap[v >> 3] |= 1 << (v & 7)
        return int.from_bytes(bitmap, "little").bit_count() / m

    out = {"r_bits": round(math.log2(r), 2), "d": d, "m_log2": round(math.log2(m), 2), "set_size": n,
           "K_effective": round(6 * n * n / m, 4), "predicted": round(1 - math.exp(-6 * n * n / m), 4)}
    # random elements of H
    rnd = lambda: [pow(rng.randrange(1, r), d, r) for _ in range(n)]
    if not products_only:
        out["random"] = round(coverage(rnd(), rnd()), 4)
    # smallest integers in H, split alternately between B and G
    small, x = [], 1
    while len(small) < 2 * n:
        if pow(x, m, r) == 1:
            small.append(x)
        x += 1
    if not products_only:
        out["interval"] = round(coverage(small[0::2], small[1::2]), 4)
    # production product construction, at several X2/d (production: 2^50)
    out["product_by_X2_over_d"] = []
    for ratio in (cfg["ratios"] if cfg.get("ratios") else ((4,) if cfg.get("only_ratio4") else (4, 16, 64, 256))):
        out["product_by_X2_over_d"].append(product_run(ratio, d, n, m, r, chi, coverage, rng, fixed, ab_below_r))
    out["scope"] = "toy prime only; evidence for H1/H2 at this scale, not a production proof"
    return out


def product_run(ratio, d, n, m, r, chi, coverage, rng, fixed=False, ab_below_r=False):
    out = {"X2_over_d": ratio}
    X2 = ratio * d
    classes = {}
    for b in range(1, X2 + 1):
        classes.setdefault(chi(b), []).append(b)
    lens = [len(classes.get(c, [])) for c in range(d)]
    X1 = 2 * X2 + 2
    while True:
        P = primes_in(X1 // 2 + 1, X1)
        if len(P) * X2 // d >= 3 * n:
            break
        X1 *= 2
    if ab_below_r:
        # production regime: every product a*b < r strictly (a <= X1, b <= X2)
        assert X1 * X2 < r, (X1, X2, r)
    if not fixed:
        rng.shuffle(P)
    else:
        P.sort()  # deterministic ascending prime sets, as in the production construction

    def build(plist):
        S = []
        for a in plist:
            for b in classes.get((-chi(a)) % d, []):
                S.append(a * b)
                if len(S) == n:
                    return S
        raise RuntimeError("not enough primes")

    B, G = build(P[: len(P) // 2]), build(P[len(P) // 2:])
    assert all(pow(v, m, r) == 1 for v in B + G) and len(set(B)) == n and len(set(G)) == n
    out["product"] = round(coverage(B, G), 4)
    out["predicted"] = round(1 - math.exp(-6 * n * n / m), 4)
    # repeated rationals b/b' = c/c' (shared factor g, cofactor k with chi(k) = chi(g)): expected repeats per
    # pair <~ 2 * E[gcd(b, b')]/d ~ 2 (6/pi^2)(ln X2 + 2*gamma - 1)/d; production: X2 = 2^70, d = 2^20 -> 2^-14.1
    delta = 2 * 6 / math.pi ** 2 * (math.log(X2) + 2 * 0.5772156649 - 1) / d
    out["conservative_gcd_bound_toy_scale"] = round(1 - math.exp(-6 * n * n / m * (1 - delta)), 4)
    out["H2_class_lengths_over_X2/d"] = {"min": round(min(lens) * d / X2, 3), "max": round(max(lens) * d / X2, 3)}
    out["X1"] = X1
    return out


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
    for b in range(1, X2 + 1):
        classes.setdefault(chi(b), []).append(b)
    X1 = 2 * X2 + 2
    while True:
        plist = primes_in(X1 // 2 + 1, X1)
        if sum(len(classes.get((-chi(a)) % d, [])) for a in plist) >= 2 * n:
            break
        X1 *= 2
        if X1 * X2 >= r_toy:            # keep a*b < r as in production
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
    X2 = 16 * d
    classes = {}
    for b in range(1, X2 + 1):
        classes.setdefault(chi(b), []).append(b)
    X1 = 2 * X2 + 2
    while len(primes_in(X1 // 2 + 1, X1)) * X2 // d < 3 * size:
        X1 *= 2
    primes = primes_in(X1 // 2 + 1, X1)

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

        walk(P, primes[: len(primes) // 2], size, baby)
        walk(Ad, primes[len(primes) // 2:], size, giant)
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
           "note": "the gcd bound is a conservative upper bound on repeated-rational loss; it grows with ln X2 "
                   "but is 2^-14.1 at production scale (X2/d = 2^50)",
           "scope": "toy curve only; does not validate production costs"}
    return out


def coverage_ladder_point(bits, logd, K, seed):
    """Exact |S|/m at X2/d = 4 with deterministic ascending prime sets (no
    shuffle) and X1*X2 < r asserted for every product, counted by direct set
    construction (no bitmap), at large d close to production's d = 2^20."""
    rng = random.Random(seed)
    d = 1 << logd
    while True:
        r = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * d) * (6 * d) + 1
        if is_prime(r) and ((r - 1) // d) % 6 == 0:
            break
    m = (r - 1) // d
    gen = next(g for g in range(2, 10**6) if all(pow(g, (r - 1) // q, r) != 1 for q in (2, 3)))
    lam = pow(gen, (r - 1) // 3, r)
    U = [1, r - 1, lam, r - lam, lam * lam % r, r - lam * lam % r]
    zeta = pow(gen, m, r)
    zlog = {pow(zeta, i, r): i for i in range(d)}
    chi = lambda x: zlog[pow(x, m, r)]
    n = math.isqrt(int(K * m / 6)) + 1
    X2 = 4 * d
    classes = {}
    for b in range(1, X2 + 1):
        classes.setdefault(chi(b), []).append(b)
    X1 = 2 * X2 + 2
    while True:
        P = primes_in(X1 // 2 + 1, X1)
        if len(P) * X2 // d >= 3 * n:
            break
        X1 *= 2
    assert X1 * X2 < r, (X1, X2, r)
    P.sort()
    def build(plist):
        S = []
        for a in plist:
            for b in classes.get((-chi(a)) % d, []):
                S.append(a * b)
                if len(S) == n:
                    return S
        raise RuntimeError("not enough primes")
    B, G = build(P[: len(P) // 2]), build(P[len(P) // 2:])
    assert all(pow(v, m, r) == 1 for v in B + G) and len(set(B)) == n and len(set(G)) == n
    vals = set()
    for u in U:
        for b in B:
            ub = u * b % r
            for g in G:
                vals.add(ub * g % r)
    return {"d": d, "r_bits": round(math.log2(r), 2), "set_size": n,
            "predicted": round(1 - math.exp(-6 * n * n / m), 4),
            "product": round(len(vals) / m, 4), "X1": X1,
            "X1_times_X2": X1 * X2, "method": "exact set construction"}


cfg = json.load(sys.stdin)
K = cfg["K"]
print(json.dumps({
    "ledger": ledger(K),
    "hash_probe": hash_probe_check({"bits": 28, "n": 1 << 14,
                                    "shards": 16, "seed": cfg["seed"]}),
    "coverage_ab_below_r": coverage({"bits": 25, "logd": 4, "K": K,
                                     "seed": cfg["seed"] + 1, "ratios": [4, 16, 64],
                                     "products_only": True, "fixed": True, "ab_below_r": True}),
    "coverage_ab_below_r_trend": [
        coverage({"bits": 29, "logd": 8, "K": K, "seed": cfg["seed"] + 3,
                  "ratios": [4], "products_only": True, "fixed": True, "ab_below_r": True}),
        coverage({"bits": 29, "logd": 10, "K": K, "seed": cfg["seed"] + 4,
                  "ratios": [4], "products_only": True, "fixed": True, "ab_below_r": True}),
    ],
    "coverage_d_ladder": [
        coverage_ladder_point(31, 12, float(K), cfg["seed"] + 5),
        coverage_ladder_point(35, 14, float(K), cfg["seed"] + 6),
    ],
    "toy_attack_fixed_order": toy_attack({"bits": cfg["atk_bits"],
                                          "logd": cfg["atk_logd"],
                                          "trials": 120, "K": K,
                                          "seed": cfg["seed"], "gmax": cfg["gmax"],
                                          "fixed": True}),
    "coverage": coverage({"bits": cfg["cov_bits"], "logd": cfg["cov_logd"], "K": K, "seed": cfg["seed"]}),
    "toy_attack": toy_attack({"bits": cfg["atk_bits"], "logd": cfg["atk_logd"], "trials": cfg["trials"],
                              "K": K, "seed": cfg["seed"], "gmax": cfg["gmax"]}),
}, indent=1))

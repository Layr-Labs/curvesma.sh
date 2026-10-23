"""Supporting evidence for proof.md (stdlib only; JSON on stdin, JSON on stdout). Four parts:

1. ledger      - recomputes every cost term of proof.md section 4, the itemised per-element instruction budgets of the
                 fixed-depth key store (two-word page-table entries at both levels and a two-word record address), the
                 entry widths against the heap size, the slot count, the exact score and the exact integer inequality
                 that proves the round-up (arithmetic check only).
2. store       - measures the fixed-depth window rule itself: with the production load factor N_B/S, a window of D slots
                 probed in order inside a block of 1024 slots and M overflow records per block, the fraction of inserts
                 whose window is full and the fraction discarded. This is a measurement of the store rule at production
                 load, at a reduced slot count.
3. coverage    - toy evidence for H1/H2: exact |U*B*G^-1|/m for random sets, smallest integers in H and the production
                 product construction a*b (chi(a)+chi(b)=0), vs 1-exp(-K) and the gcd-corrected prediction; class spread.
4. toy_attack  - the full algorithm on a toy j=0 curve (walks at one addition per element, keyed on canonical y =
                 min(y, p-y), invariant under U = {±1, ±λ, ±λ²}; collision, d-th root, mu_d search), success rate vs
                 prediction, scalars checked. Its key table is an exact dictionary, so it exercises the walk and the
                 coverage argument, not the store's window rule (part 2 does that).
Toy and reduced-size results are evidence at their scale only, not production cost or success proofs.
Input: {"K", "cov_bits", "cov_logd", "atk_bits", "atk_logd", "trials", "gmax", "seed", "store_bits", "store_seeds"}.
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F
C = 2000
r = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
p = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
d = 1 << 20
m = (r - 1) // d
assert (r - 1) % d == 0 and m % 6 == 0

PROBE = 8        # D: fixed probe depth, slots h, h+1, ..., h+7 inside the block, wrapping at the block's end
OVER = 4         # M: overflow records per block of 1024 primary slots
DISCARD = F(1, 100)   # charged allowance for keys the store discards (proof.md section 3, H3)
BLOCK = 1024     # primary slots per block
ENTRY = 16       # bytes per page-table entry: two 64-bit words, at both levels (proof.md section 4.2)
L1BITS, L2BITS = 54, 53   # index bits per level, 54 + 53 = 107 >= the block number's 107 bits


def floor_log2(frac):
    """Exact floor(log2(frac)) for a positive Fraction, by integer comparison: 2^e <= frac < 2^(e+1)."""
    e = frac.numerator.bit_length() - frac.denominator.bit_length()
    c = frac / F(1 << e) if e >= 0 else frac * F(1 << -e)     # frac / 2^e, a Fraction near 1
    while c >= 2:
        c, e = c / 2, e + 1
    while c < 1:
        c, e = c * 2, e - 1
    return e


def exact_rounding(W, claimed):
    """Prove claimed = ceil(log2(W) * 10^6) with exact integers and no floating point. W = 2^E * f with 1 <= f < 2, so
    W^10^6 = 2^(E*10^6) * f^10^6 and the claim is 2^(m-1) < f^10^6 <= 2^m with m = claimed - E*10^6. f^1000 is
    bracketed by g/2^40 <= f^1000 < (g+1)/2^40 with the integer g, and raising to the 1000th power keeps the bracket:
    the two integer comparisons below are exact."""
    e, sh = 1000, 40
    E = floor_log2(W)
    m = claimed - E * 10**6
    fn = W.numerator ** e
    fd = W.denominator ** e * (1 << (e * E))
    g = (fn << sh) // fd
    lo = g ** e > (1 << (m - 1 + sh * e))
    hi = (g + 1) ** e <= (1 << (m + sh * e))
    return {"exponent_E": E, "mantissa_ceiling_m": m, "scale_bits": sh, "root_power": e, "g_bits": g.bit_length(),
            "2^%d < g^1000" % (m - 1 + sh * e): lo, "(g+1)^1000 <= 2^%d" % (m + sh * e): hi, "proved": bool(lo and hi),
            "inequality": "2^(%d) < W^1000000 < 2^(%d) after scaling by 2^E and by 2^40 per 1000th root"
                          % (E * 10**6 + m - 1, E * 10**6 + m)}


def ledger(K):
    K = F(K)

    # Itemised per-element instruction budgets (proof.md section 4.3). CANON and WALK are shared by both walks. ADDR is
    # the shared access prologue: the two-word slot index and block number, the two-level page-table walk whose entries
    # are two words at both levels, and the two-word record address. INSERT and LOOKUP are the twelve slot visits. Every
    # budget is the itemised count plus a slack of SLACK.
    CANON, WALK, SLACK = 26, 80, 2
    ADDR = {
        "slot index j (118 bits, two words: shift/shift/or per word)": 6,
        "block number P = j >> 10 (two words)": 4,
        "level-1 index P >> 53 (two words: shift, shift, or)": 3,
        "level-1 entry address T1base + 16*index (shift, add)": 2,
        "level-1 entry load (two words, 16-byte entry)": 2,
        "level-2 index P mod 2^53 (and)": 1,
        "level-2 entry address: two-word table base + 16*index (shift, add, adc)": 3,
        "level-2 entry load (two words, 16-byte entry)": 2,
        "in-block slot offset o = j mod 1024 (and)": 1,
        "in-block byte offset 72*o (shift, shift, add)": 3,
        "record address: two-word block base + 72*o (add, adc)": 2,
    }
    INSERT = {
        "8 window probes x [two-word advance and in-block wrap (5), flag load (1), flag test (1), branch (1)]":
            8 * 8,
        "overflow base: two-word block base + 73,728 (add, adc)": 2,
        "4 overflow probes x [two-word advance (2), flag load (1), flag test (1), branch (1)]": 4 * 5,
        "store the nine record words from the two-word address": 9,
        "class-list advance and gap handling": 9,
    }
    LOOKUP = {
        "8 window probes x [two-word advance and in-block wrap (5), flag test (1), 6 key words loaded and compared "
        "with branches (15)]": 8 * 21,
        "overflow base: two-word block base + 73,728 (add, adc)": 2,
        "4 overflow records x [two-word advance (2), flag test (1), 6 key words compared (15)]": 4 * 18,
        "class-list advance and gap handling": 9,
    }
    addr, ins, look = sum(ADDR.values()), sum(INSERT.values()), sum(LOOKUP.values())
    ACC_B, ACC_G = addr + ins + SLACK, addr + look + SLACK
    c_B = 1 + F(CANON + ACC_B + WALK, C)   # walk add + canonical y + 1 store insert + bookkeeping
    c_G = 1 + F(CANON + ACC_G + WALK, C)   # walk add + canonical y + 1 store lookup + bookkeeping

    # K is the coverage exponent the claim needs. The store discards a fraction at most DISCARD of the stored keys, so
    # the charge is K/(1-DISCARD) and the effective exponent is K again.
    K_charged = K / (1 - DISCARD)

    # sizes: 6*N_B*N_G*(1-DISCARD) >= K*m, N_B/N_G = c_G/c_B
    N_B = math.isqrt(int(K_charged * m * c_G / (6 * c_B))) + 1
    N_G = -(-int(K_charged * m) // (6 * N_B)) + 1
    assert 6 * N_B * N_G >= K_charged * m

    X1, X2, GMAX = 1 << 76, 1 << 70, 1 << 30
    SCALAR = 1000          # units: one non-batched scalar multiplication by a <= 2^76-bit scalar (proof.md 4.4)
    SCALAR_FULL = 3000     # units: one scalar multiplication by a 255-bit scalar
    PI_X2 = F(126, 100) * X2 / F(int(70 * math.log(2) * 10**6), 10**6)    # pi(x) < 1.26 x / ln x
    PI_X1 = F(126, 100) * X1 / F(int(76 * math.log(2) * 10**6), 10**6)
    N_A = (N_B + N_G) * 2 * d // X2 + 2                    # primes a used, heuristic H2 (lists >= X2/(2d))
    assert N_A < PI_X1 / 4                                  # enough primes in (X1/2, X1]

    # fixed-depth store (proof.md 4.2): S slots of 72 bytes indexed by the 118-bit prefix of the canonical y, laid out
    # 1024 primary slots plus M overflow records per block. S+1 possible prefix values, so the load is N_B/(S+1).
    S = (((p - 1) // 2) >> 263) + 1
    blocks = -(-S // BLOCK)
    records = blocks * (BLOCK + OVER)
    tables2 = -(-blocks >> L2BITS)                          # level-2 tables, each covering 2^53 block numbers
    assert L1BITS + L2BITS >= blocks.bit_length()
    l1_bytes = (1 << L1BITS) * ENTRY                        # 16-byte entries: two 64-bit words per entry
    l2_bytes = tables2 * (1 << L2BITS) * ENTRY
    pt_bytes = l1_bytes + l2_bytes
    pt_words = pt_bytes // 8
    # the widths the two words of each entry have to hold: a block's byte address in the heap, and a level-2 table's
    # byte address in the page-table region
    block_bytes = (BLOCK + OVER) * 72
    block_base_bits = (block_bytes * (blocks - 1)).bit_length()
    l2_base_bits = (tables2 * (1 << L2BITS) * ENTRY).bit_length()
    assert block_base_bits <= ENTRY * 8 and l2_base_bits <= ENTRY * 8

    terms = {
        "baby walk (N_B * c_B)": N_B * c_B,
        "giant walk (N_G * c_G)": N_G * c_G,
        "gap overflow additions (<= N_A * X2 / GMAX)": F(N_A * X2 // GMAX),
        "per-prime tables [k]Q_a, k <= GMAX": F(N_A * GMAX),
        "per-prime scalar mults ([a]base, first element, doubling)": F(N_A * 3 * SCALAR),
        "chi for primes <= X1 (100 units each)": 100 * PI_X1,
        "prime sieve to X1, chi(b) sieve and class lists to X2": F(X1 * 64 + X2 * 200, C),
        "store and page-table initialisation (records*9 + page-table words)": F(records * 9 + pt_words, C),
        "collision verification (6 candidates)": F(6 * SCALAR_FULL),
        "d-th root in F_r (<= 2^20 F_r mults of 200 instr)": F((1 << 20) * 200, C),
        "mu_d search (2^20 scalar mults) + final check": F(((1 << 20) + 1) * SCALAR_FULL),
        "input decoding and subgroup checks": F(4 * SCALAR_FULL),
    }
    W = sum(terms.values())
    pre_keys = ["baby walk (N_B * c_B)", "per-prime tables [k]Q_a, k <= GMAX",
                "per-prime scalar mults ([a]base, first element, doubling)",
                "chi for primes <= X1 (100 units each)", "prime sieve to X1, chi(b) sieve and class lists to X2",
                "gap overflow additions (<= N_A * X2 / GMAX)",
                "store and page-table initialisation (records*9 + page-table words)"]
    pre = sum(terms[k] for k in pre_keys)         # challenge-independent work (upper bound: counts both walks' tables)
    c_kkm = 2 + F(CANON + ACC_G + WALK + 200, C)   # KKM c=3: 2 additions + one F_r mul to advance zeta^j
    kkm = 2 * c_kkm * F(math.isqrt(int(K * m / 6)))
    # exact score: ceil(log2(W) * 10^6) / 10^6, from the exact rational W and an exact integer proof of the rounding
    getcontext().prec = 80
    log2W = (Decimal(W.numerator).ln() - Decimal(W.denominator).ln()) / Decimal(2).ln()
    claimed = int((log2W * 10**6).to_integral_value(rounding="ROUND_CEILING"))
    score = F(claimed, 10**6)
    rounding = exact_rounding(W, claimed)
    assert rounding["proved"] and score == F(116889745, 10**6)   # the claim proved by rounding
    success = 1 - math.exp(-float(6 * N_B * N_G) / m * float(1 - DISCARD))   # store keeps at least 1-DISCARD

    def lg(v):
        return round(math.log2(v.numerator) - math.log2(v.denominator), 3)

    out = {
        "K_effective": str(K), "K_charged": str(K_charged), "discard_allowance": str(DISCARD),
        "c_B": str(c_B), "c_G": str(c_G), "probe_depth_D": PROBE, "overflow_M": OVER,
        "element_budget": {"canonical y": CANON, "store insert (baby), itemised": addr + ins,
                           "store lookup (giant), itemised": addr + look,
                           "budget slack per access": SLACK,
                           "store insert (baby), budgeted": ACC_B, "store lookup (giant), budgeted": ACC_G,
                           "walk bookkeeping": WALK,
                           "per-element total (baby)": CANON + ACC_B + WALK,
                           "per-element total (giant)": CANON + ACC_G + WALK,
                           "store probes per access (D window + M overflow)": PROBE + OVER},
        "access_itemisation": {"shared addressing (instructions)": addr, "insert-only rows (instructions)": ins,
                               "lookup-only rows (instructions)": look,
                               "shared": ADDR, "insert": INSERT, "lookup": LOOKUP,
                               "two-word advance and in-block wrap per slot visit": 5},
        "page_table": {"entry_bytes": ENTRY, "entry_words": ENTRY // 8,
                       "level1_index_bits": L1BITS, "level2_index_bits": L2BITS,
                       "level1_entries": 1 << L1BITS, "level1_bytes_log2": round(math.log2(l1_bytes), 3),
                       "level2_tables_log2": round(math.log2(tables2), 6),
                       "level2_bytes_log2": round(math.log2(l2_bytes), 6),
                       "level2_words": l2_bytes // 8, "page_table_bytes_log2": round(math.log2(pt_bytes), 3),
                       "block_base_bits_held": block_base_bits, "level2_table_base_bits_held": l2_base_bits,
                       "block_base_spare_bits": ENTRY * 8 - block_base_bits,
                       "level2_table_base_spare_bits": ENTRY * 8 - l2_base_bits,
                       "block_base_fits_one_64bit_word": block_base_bits <= 64,
                       "level2_table_base_fits_one_64bit_word": l2_base_bits <= 64,
                       "fits_entry_bits": ENTRY * 8},
        "slot_count_log2": round(math.log2(S), 6), "block_log2": round(math.log2(blocks), 6),
        "block_bytes": block_bytes, "load_factor": round(N_B / S, 5),
        "log2_m": round(math.log2(m), 6), "log2_N_B": round(math.log2(N_B), 6), "log2_N_G": round(math.log2(N_G), 6),
        "log2_primes_a_used": round(math.log2(N_A), 3),
        "terms_log2": {k: lg(v) for k, v in terms.items()},
        "log2_W": round(float(log2W), 9),
        "log2_preprocessing": lg(pre),
        "log2_W_KKM_c3_same_budgets": lg(kkm),
        "score": {"numerator": str(score.numerator), "denominator": str(score.denominator)},
        "score_exact_check": rounding,
        "heuristic_success": round(success, 6),
        "store_bytes_log2": round(math.log2(records * 72), 4),
        "page_table_bytes_log2": round(math.log2(pt_bytes), 4),
        "log2_memory_bytes": round(math.log2(records * 72 + pt_bytes + X2 * 9 + 2**46 * 96), 4),
        "memory_claim_log2": 122.88,
    }
    return out


def store(load, bits, seeds):
    """Measure the fixed-depth window rule at the production load factor. Slots h..h+D-1 inside a block of 1024 primary
    slots (wrapping at the block's end) are probed in order; a key that finds none free takes one of the M overflow
    records of its block; a key that finds neither is discarded. Uniform homes, which proof.md section 1 establishes for
    the production index. Reduced slot count; the rule and the load are the production ones."""
    out = {"slots": 1 << bits, "load": round(load, 5), "D": PROBE, "blocks_of": BLOCK, "by_seeds": []}
    for M in (0, 2, OVER):
        win = dis = 0.0
        for seed in range(1, seeds + 1):
            rng = random.Random(seed)
            S = 1 << bits
            n = int(load * S)
            occ = bytearray(S)
            cnt = bytearray(S // BLOCK)
            w = dcount = 0
            for _ in range(n):
                j = rng.randrange(S)
                b, o = divmod(j, BLOCK)
                placed = False
                for t in range(PROBE):
                    s = b * BLOCK + (o + t) % BLOCK
                    if not occ[s]:
                        occ[s] = 1
                        placed = True
                        break
                if not placed:
                    w += 1
                    if cnt[b] < M:
                        cnt[b] += 1
                    else:
                        dcount += 1
            win += w / n
            dis += dcount / n
        out["by_seeds"].append({"M": M, "window_full_fraction": round(win / seeds, 6),
                                "discard_fraction": round(dis / seeds, 6),
                                "discard_log2": round(math.log2(max(dis / seeds, 2.0 ** -60)), 3)})
    best = [x for x in out["by_seeds"] if x["M"] == OVER][0]
    out["charged_discard_allowance"] = str(DISCARD)
    out["charge_over_measured_discard"] = round(float(DISCARD) / max(best["discard_fraction"], 1e-12), 1)
    out["charge_over_window_full"] = round(float(DISCARD) / max(best["window_full_fraction"], 1e-12), 2)
    out["scope"] = ("reduced slot count with the production load factor; the rule is local, so the fraction depends on "
                    "the load and on D and M, not on the slot count")
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
    out["random"] = round(coverage(rnd(), rnd()), 4)
    # smallest integers in H, split alternately between B and G
    small, x = [], 1
    while len(small) < 2 * n:
        if pow(x, m, r) == 1:
            small.append(x)
        x += 1
    out["interval"] = round(coverage(small[0::2], small[1::2]), 4)
    # production product construction, at several X2/d (production: 2^50)
    out["product_by_X2_over_d"] = []
    for ratio in (4, 16, 64, 256):
        out["product_by_X2_over_d"].append(product_run(ratio, d, n, m, r, chi, coverage, rng))
    out["scope"] = "toy prime only; evidence for H1/H2 at this scale, not a production proof"
    return out


def product_run(ratio, d, n, m, r, chi, coverage, rng):
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
    rng.shuffle(P)

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
    # repeated rationals b/b' = c/c' (shared factor g, cofactor k with chi(k) = chi(g)): repeats per pair below
    # 2 * E[gcd(b, b')]/d ~ 2 (6/pi^2)(ln X2 + 2*gamma - 1)/d; production: X2 = 2^70, d = 2^20 -> 2^-14.1
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
    for _ in range(trials):
        s = rng.randrange(1, r)
        A1, Ad = E.mul(s, P), E.mul(pow(s, d, r), P)
        E.adds = 0
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
           "key_table": "exact dictionary, so the walk and the coverage argument are tested, not the store's window rule",
           "scope": "toy curve only; does not validate production costs"}
    return out


cfg = json.load(sys.stdin)
K = cfg["K"]
_led = ledger(K)
print(json.dumps({
    "ledger": _led,
    "store": store(_led["load_factor"], cfg["store_bits"], cfg["store_seeds"]),
    "coverage": coverage({"bits": cfg["cov_bits"], "logd": cfg["cov_logd"], "K": K, "seed": cfg["seed"]}),
    "toy_attack": toy_attack({"bits": cfg["atk_bits"], "logd": cfg["atk_logd"], "trials": cfg["trials"],
                              "K": K, "seed": cfg["seed"], "gmax": cfg["gmax"]}),
}, indent=1))

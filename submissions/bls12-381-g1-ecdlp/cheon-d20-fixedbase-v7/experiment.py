"""Toy two-phase Cheon check plus production parameter/ledger checks.

The toy group is the additive group of integers modulo a toy prime. It has the
same prime-order cyclic structure as G1, but the experiment is evidence only at
toy scope. It does not execute or validate the production-size attack.
"""
import json
import math
import random
import sys
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60

R = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
D = 1 << 20
BITS = 85
WINDOWS = 3
C = 2000
HASH_CHARGE = 560          # P2: 206 itemised + 13 per extra link, worst case 284, charged 1.97x
NODE_BYTES = 128           # 96-byte affine point + 128-bit tag + 128-bit next reference
HEAD_BYTES = 16            # bucket head is a 128-bit reference: 8 bytes cannot address the arena
N_BUCKETS = 1 << 119       # chained matching table: 2^119 bucket heads
P2_BASE = 206              # itemised per-operation base, 16-byte head + 128-byte node + y-sign fold
P2_LINK = 13               # itemised cost of one extra chain link
LINK_BUDGET = 6            # GLOBAL budget of 6*N links per attempt: no per-operation cap
P2_ALLOW = P2_BASE + LINK_BUDGET * P2_LINK   # 284 = 206 + 6*13, inside the 560 charge
TAG_SHIFT = 137            # t(x) = bits 137..264 of p, immediately BELOW h(x) = p >> 265
TAG_MASK = (1 << 128) - 1  # so (h,t) is the TOP 247 bits of p: 2/2^247, and a != 0 mod 2^128
# ---- distributed arena: how a 128-bit logical address becomes a real access ----
# The matching table and the six comb tables do not fit any single machine, so the
# arena is distributed: 2^77 nodes with 2^48 bytes (256 TiB) of local DRAM each,
# total 2^125 bytes >= the 2^124.88-byte arena.  A 128-bit logical address is a
# 77-bit node field followed by a 51-bit byte offset inside that node, so no
# AArch64 LDR/STR ever takes a register pair as an effective address: the address
# is decomposed once, the request is routed to the owning node, the owning node
# performs a 64-bit local load/store, and the payload is routed back.
ARENA_NODE_BITS = 77       # ceil(log2(2^124.88 / 2^48)) = 77
ARENA_OFF_BITS = 51        # 77 + 51 = 128-bit logical address
ARENA_LEVELS = 11          # flattened butterfly, 7 node-id bits per level: 11*7 = 77
ARENA_HOP = 4              # per level: LDR the (level, 7-bit group) entry, EXTR/AND, ORR splice
ARENA_ADDR_SPLIT = 8       # UBFX/AND both fields out of the 128-bit address
ARENA_REQ_BUILD = 6        # serialise the request header (op, node, offset, nonce)
ARENA_DEST_FIX = 4         # destination: base+offset arithmetic and bounds check
ARENA_RESP_FIX = 14        # tag compare, counter update, branch, chain-state store


def arena_access(words):
    """Itemised instructions for one remote read of `words` 8-byte words."""
    return (ARENA_ADDR_SPLIT + ARENA_REQ_BUILD + ARENA_LEVELS * ARENA_HOP
            + ARENA_DEST_FIX + words + ARENA_LEVELS * 2 * words + ARENA_RESP_FIX)


CMP_BUDGET = 1 << 20       # GLOBAL budget of false full-point comparisons per attempt
CMP_CHARGE = 24            # instructions for one full 96-byte point comparison


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
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


def distinct_prime_factors(n):
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors.append(n)
    return factors


def primitive_root(p):
    factors = distinct_prime_factors(p - 1)
    g = 2
    while True:
        if all(pow(g, (p - 1) // q, p) != 1 for q in factors):
            return g
        g += 1


def toy_prime(seed, d_bits):
    k = d_bits
    q = ((1 << (30 - k)) | 1) + 2 * (seed % 64)
    while True:
        r0 = q * (1 << k) + 1
        if r0.bit_length() >= 30 and is_prime(r0):
            return r0, q
        q += 2


def window_tables(base, mod, w, windows):
    tables = []
    for j in range(windows):
        step = (base * pow(2, w * j, mod)) % mod
        tab = [0] * (1 << w)
        acc = 0
        for c in range(1 << w):
            tab[c] = acc
            acc = (acc + step) % mod
        tables.append(tab)
    return tables


def window_mul(x, tables, w, windows, mod):
    top = (x >> (w * (windows - 1))) & ((1 << w) - 1)
    acc = tables[windows - 1][top]
    for j in range(windows - 2, -1, -1):
        digit = (x >> (w * j)) & ((1 << w) - 1)
        acc = (acc + tables[j][digit]) % mod
    return acc


def naive_mul(k, base, mod):
    acc = 0
    addend = base % mod
    while k:
        if k & 1:
            acc = (acc + addend) % mod
        addend = (addend + addend) % mod
        k >>= 1
    return acc


class ChainedTable:
    """Separate chaining with the declared 2-universal multiply-shift hash.

    Production uses h(x) = ((a*x + b) mod 2^384) >> 265 with a odd and a, b
    drawn uniformly per attempt; the toy uses the same family at the toy key
    width w = r0.bit_length() and l bucket bits. Every chain link walked by an
    insert or a lookup increments `links`. When `budget` is set and the counter
    would exceed it, `abort` is raised: this is the single GLOBAL counter of the
    proof. There is no per-operation cap and no chain is ever truncated.
    """

    def __init__(self, w, l, rng, budget=None):
        self.shift = w - l
        self.mask = (1 << l) - 1
        self.a = rng.getrandbits(w) | 1
        self.b = rng.getrandbits(w)
        self.buckets = [[] for _ in range(1 << l)]
        self.links = 0
        self.budget = budget
        self.abort = False
        self.max_chain = 0

    def _bucket(self, key):
        return ((self.a * key + self.b) >> self.shift) & self.mask

    def _charge(self, n):
        self.links += n
        if self.budget is not None and self.links > self.budget:
            self.abort = True

    def insert(self, key, value):
        chain = self.buckets[self._bucket(key)]
        self._charge(len(chain))
        chain.append((key, value))
        if len(chain) > self.max_chain:
            self.max_chain = len(chain)

    def lookup(self, key):
        chain = self.buckets[self._bucket(key)]
        self._charge(len(chain))
        for k, v in chain:
            if k == key:
                return v
        return None

    def stats(self):
        return {
            "links": self.links,
            "budget": self.budget,
            "max_chain": self.max_chain,
            "abort": self.abort,
        }


def cheon_toy(r0, zeta, d0, secret, seed):

    M1 = (r0 - 1) // d0
    m = math.isqrt(M1)
    if m * m < M1:
        m += 1
    g1 = (M1 + m - 1) // m
    zeta_d = pow(zeta, d0, r0)
    sd = pow(secret, d0, r0)

    w = max(1, math.ceil(r0.bit_length() / WINDOWS))
    baby_tabs = window_tables(sd, r0, w, WINDOWS)
    giant_tabs = window_tables(1, r0, w, WINDOWS)

    rng = random.Random(seed)
    w_bits = r0.bit_length()
    L_BITS = max(6, m.bit_length() + 2)   # mirrors the production bucket/node ratio 2^119/m
    link_budget = LINK_BUDGET * (m + g1)
    lookup = ChainedTable(w_bits, L_BITS, rng, budget=link_budget)
    scalar = 1
    inv = pow(zeta_d, -1, r0)
    for u in range(m):
        point = window_mul(scalar, baby_tabs, w, WINDOWS, r0)
        lookup.insert(point, u)
        scalar = scalar * inv % r0

    omega = pow(zeta, M1, r0)
    m2 = math.isqrt(d0)
    if m2 * m2 < d0:
        m2 += 1
    g2 = (d0 + m2 - 1) // m2

    target = secret
    scalar = 1
    step = pow(zeta_d, m, r0)
    for v in range(g1):
        giant = window_mul(scalar, giant_tabs, w, WINDOWS, r0)
        u_hit = lookup.lookup(giant)
        if u_hit is not None:
            k0 = (u_hit + m * v) % M1
            tprime = pow(zeta, -k0, r0) * target % r0
            baby2 = ChainedTable(w_bits, L_BITS, rng, budget=link_budget)
            s2 = 1
            inv2 = pow(omega, -1, r0)
            for up in range(m2):
                baby2.insert(naive_mul(s2, tprime, r0), up)
                s2 = s2 * inv2 % r0
            s3 = 1
            step2 = pow(omega, m2, r0)
            for vp in range(g2):
                w_hit = baby2.lookup(s3)
                if w_hit is not None:
                    k1 = w_hit + m2 * vp
                    recovered = pow(zeta, (k0 + M1 * k1) % (r0 - 1), r0)
                    if recovered == target and pow(recovered, d0, r0) == sd:
                        st = lookup.stats()
                        st["buckets"] = 1 << L_BITS
                        st["phase2_links"] = baby2.links
                        st["phase2_max_chain"] = baby2.max_chain
                        return True, k0, k1, st
                s3 = s3 * step2 % r0
        scalar = scalar * step % r0
    st = lookup.stats()
    st["buckets"] = 1 << L_BITS
    st["phase2_links"] = 0
    st["phase2_max_chain"] = 0
    return False, None, None, st


def digit_checks():
    """Check on production integers that the 3x85-bit comb has no exception."""
    B = 1 << BITS
    xs = {0, 1, 2, B - 1, B, B + 1, R - 1, R - 2, (R - 1) // 2, (R - 1) // 3}
    seed = 123456789
    for _ in range(64):
        seed = (seed * 6364136223846793005 + 1442695040888963407) & ((1 << 255) - 1)
        xs.add(seed % R)
    for x in xs:
        x0 = x % B
        a = x // B
        x1 = a % B
        x2 = a // B
        assert x0 < B and x1 < B and x2 < B
        # first addition is safe unless both table digits are zero
        assert (B * x2 == x1 or B * x2 + x1 == 0) == (x1 == 0 and x2 == 0)
        # second addition: S = [x - x0]B equals +/- T[B,0][x0] only for x = 0
        if x > 0:
            assert x != 2 * x0
            assert x0 != B * (B * x2 + x1)
    return True


def production_checks():
    out = {}
    assert (R - 1) % D == 0
    assert R < (1 << 255)
    assert digit_checks()
    M1 = (R - 1) // D
    assert M1 == 50006747412801924209068050869165387952509453297164571592906626415194112
    m = math.isqrt(M1)
    if m * m < M1:
        m += 1
    g1 = (M1 + m - 1) // m
    N = m + g1
    assert m == 223621884914696853578958571250384896
    assert g1 == m
    assert N == 447243769829393707157917142500769792
    assert (m - 1) * (m - 1) < M1 <= m * m
    assert m * g1 - M1 == 218380746987008646073201729736704
    assert 0 < m * g1 - M1 < m

    h_pre = 3 * 6 * ((1 << BITS) - 1) + 510
    h_main = 2 * N
    h_p2 = 1025 * 384 + 1025 * 2
    h_exc = 1024
    h_total = h_pre + h_main + h_p2 + h_exc
    r_lookup = N * 3 * 200
    r_hash = N * HASH_CHARGE
    r_scal = N * 300
    r_digits = N * 100
    c16 = arena_access(2)        # 16-byte record: a bucket head, or (tag, next reference)
    c96 = arena_access(12)       # 96-byte affine comb entry
    # Every arena access is a remote access: 3 comb lookups per multiplication and
    # 1 + LINK_BUDGET table operations per multiplication (one bucket head plus the
    # links the global counter can afford).  Hard cap, no expectation.
    r_arena = N * (3 * c96 + (1 + LINK_BUDGET) * c16)
    r_stores = 6 * ((1 << BITS) - 1) * c96
    r_init = N_BUCKETS * c16
    r_setup = 1 << 40
    r_p2 = 1 << 22
    r_out = 1 << 12
    r_cmp = CMP_BUDGET * CMP_CHARGE      # hard allowance of false full-point comparisons
    r_total = (r_lookup + r_hash + r_scal + r_digits + r_stores + r_init + r_setup
               + r_p2 + r_out + r_cmp + r_arena)

    W = Decimal(h_total) + Decimal(r_total) / C
    log2_w = W.ln() / Decimal(2).ln()
    claimed = Decimal(120364249) / Decimal(1000000)
    assert log2_w <= claimed, (log2_w, claimed)

    prep = Decimal(h_pre) + Decimal(r_stores + r_setup + r_init) / C
    log2_prep = prep.ln() / Decimal(2).ln()
    assert log2_prep <= Decimal(115100000) / Decimal(1000000)   # declared 2^115.1

    mem = N_BUCKETS * HEAD_BYTES + NODE_BYTES * m + 6 * ((1 << BITS) - 1) * 96 + (1 << 27)
    log2_mem = Decimal(mem).ln() / Decimal(2).ln()
    assert log2_mem <= Decimal(125000000) / Decimal(1000000)   # declared 2^125
    # The distributed layout must actually hold the arena, and the address split must
    # be total: 2^77 nodes x 2^48 bytes >= mem, node field + offset = 128 bits, and
    # ARENA_LEVELS levels of 7 node-id bits cover the 77-bit node field.
    assert (1 << ARENA_NODE_BITS) * (1 << 48) >= mem, "arena does not fit the machine"
    assert ARENA_NODE_BITS + ARENA_OFF_BITS == 128, "address split is not total"
    assert ARENA_LEVELS * 7 >= ARENA_NODE_BITS, "network too shallow"
    assert ARENA_LEVELS == -(-ARENA_NODE_BITS // 7), "extra levels are not charged"
    assert arena_access(2) == 122 and arena_access(12) == 352, "itemisation changed"
    # The joint (bucket, tag) event is an event on the TOP 247 bits of p, so the
    # 2-universality factor is 2/2^247.  The low 128 bits would NOT give this: for
    # odd a, x = y (mod 2^128) forces equal low tags for every a.
    assert TAG_SHIFT + 128 <= 265, (TAG_SHIFT, 265)

    # Matching cost model: separate chaining with the multiply-shift hash.
    # Every pair of distinct keys collides with probability at most 2/2^119, and
    # the total link count is a sum of such indicators, so no independence
    # between operations is needed (only pairwise independence + Markov).
    e_insert = F(m * (m - 1), N_BUCKETS)        # C(m,2) * 2/2^119
    e_lookup = F(2 * m * g1, N_BUCKETS)         # m*g1 * 2/2^119
    e_probes = N + e_insert + e_lookup
    ratio = e_probes / N
    abort = e_probes / (LINK_BUDGET * N)
    # Second abort event: false full-point comparisons.  The first moment is taken
    # over the m*g1 POTENTIAL pairs (stored baby node, walked giant key), which are
    # fixed before any coin is drawn -- not over walked links, whose set is a
    # function of the coins (review 23face18, defect 1).  Each such pair collides in
    # the 247-bit (bucket, tag) window with probability at most 2/2^247.
    # Both pair classes are counted: stored baby vs walked giant, and stored baby vs
    # inserted baby (insertion walks the chain and compares too -- review 693583a6).
    e_false_cmp = F(2 * m * g1, 1 << 247) + F(m * (m - 1), 1 << 247)
    # The counter caps ALL full comparisons while this probability counts only FALSE
    # ones.  At most one comparison per attempt can be true (P and -P are separated
    # by the compressed encoding), so all <= false + 1 and
    # Pr[abort] <= Pr[false comparisons >= CMP_BUDGET - 1] <= E/(CMP_BUDGET - 1).
    cmp_abort = e_false_cmp / (CMP_BUDGET - 1)
    success = 1 - abort - cmp_abort
    assert ratio < F(151, 100), ratio
    assert abort < F(251, 1000), abort    # = E[probes]/(6N) = 0.2508
    assert cmp_abort < F(1, 1000000), cmp_abort
    assert success > F(74, 100), success
    assert success >= F(39, 100), success
    # The itemised P2 accounting must fit inside the per-operation charge.
    # 16-byte bucket head (+1 address, +1 for the second LDR word), a middle-window
    # tag extraction (+4) and the compressed-encoding y-sign fold (+3): the 2-word
    # next-reference load costs 5, the 16 x STR push 17.
    p2_itemised = 108 + 24 + 7 + 3 + 3 + 5 + 5 + 24 + 17 + 10
    p2_allow = p2_itemised + LINK_BUDGET * P2_LINK
    p2_work = (N * p2_itemised + LINK_BUDGET * N * P2_LINK
               + CMP_BUDGET * CMP_CHARGE)
    assert p2_itemised == P2_BASE, p2_itemised
    assert p2_allow == P2_ALLOW, p2_allow
    assert p2_allow <= HASH_CHARGE, p2_allow
    assert p2_work <= N * HASH_CHARGE, p2_work

    product = 1
    for q, e in ((2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1),
                 (859267, 1), (906349, 2), (2508409, 1), (2529403, 1),
                 (52437899, 1), (254760293, 2)):
        assert is_prime(q)
        product *= q ** e
    assert product == R - 1

    root_primes = (2, 3, 11, 19, 10177, 125527, 859267, 906349,
                   2508409, 2529403, 52437899, 254760293)
    for g0 in range(2, 8):
        is_root = all(pow(g0, (R - 1) // q, R) != 1 for q in root_primes)
        assert is_root == (g0 == 7), (g0, is_root)

    out.update({
        "M1_bits": M1.bit_length(),
        "m": str(m),
        "g1": str(g1),
        "N": str(N),
        "H": str(h_total),
        "R": str(r_total),
        "r_arena": str(r_arena),
        "log2_W": str(log2_w),
        "claimed_log2": "120.364249",
        "arena_access_16B": c16,
        "arena_access_96B": c96,
        "arena_accesses_per_N": 3 + 1 + LINK_BUDGET,
        "r_arena_log2": "%.6f" % math.log2(r_arena),
        "arena_bytes_log2": "%.6f" % math.log2((1 << ARENA_NODE_BITS) * (1 << 48)),
        "ledger_ok": True,
        "digit_checks_ok": digit_checks(),
        "preprocessing_log2": str(log2_prep),
        "hash_charge": HASH_CHARGE,
        "link_budget_per_op": LINK_BUDGET,
        "buckets": str(N_BUCKETS),
        "node_bytes": NODE_BYTES,
        "head_bytes": HEAD_BYTES,
        "tag_shift": TAG_SHIFT,
        "cmp_budget": CMP_BUDGET,
        "cmp_charge_instructions": CMP_BUDGET * CMP_CHARGE,
        "E_probes_over_N": "%.6f" % float(ratio),
        "abort_bound": "%.6f" % float(abort),
        "success_lower_bound": "%.6f" % float(success),
        "cmp_abort_log2": "%.4f" % (math.log2(float(cmp_abort))),
        "e_false_cmp_log2": "%.4f" % (math.log2(float(e_false_cmp))),
        "p2_itemised": p2_itemised,
        "p2_allow": p2_allow,
        "p2_work_ok": p2_work <= N * HASH_CHARGE,
        "memory_log2": str(log2_mem),
        "memory_declared_log2": "125.0",
        "factorization_ok": True,
        "primitive_root_ok": True,
    })
    return out



def sign_and_comparison_checks(r0, seed=20261022):
    """Model the two things review 23face18 said the toy did not model.

    (1) An ECDLP x-coordinate sign collision.  In the exponent ring the analogue
        of "P and -P share x" is "a and r0-a share the key", so the x-only key is
        min(a, r0-a) and the compressed key is a itself.
    (2) The production comparison counter.  A chained table whose nodes carry a
        tag: only a tag match triggers a full key comparison, every such
        comparison is counted, and one whose keys differ is a FALSE comparison.

    The x-only run must show the forced false comparison; the compressed run must
    show none; both must still find the true match.
    """
    w = r0.bit_length()
    tau = 20                                   # toy tag width (production: 128)
    l = max(6, min(16, w - tau - 2))
    a_coin = (seed * 2654435761 + 1) | 1
    b_coin = (seed * 40503 + 7) % (1 << w)

    def bucket_tag(key):
        p = (a_coin * key + b_coin) % (1 << w)
        return p >> (w - l), (p >> (w - l - tau)) & ((1 << tau) - 1)

    class Tab:
        def __init__(self):
            self.nodes = {}
            self.cmp = 0
            self.false_cmp = 0

        def insert(self, key, full):
            h, t = bucket_tag(key)
            self.nodes.setdefault(h, []).append((t, full))

        def lookup(self, key, full):
            h, t = bucket_tag(key)
            for tt, ff in self.nodes.get(h, []):
                if tt == t:
                    self.cmp += 1
                    if ff == full:
                        return True
                    self.false_cmp += 1
            return False

    half = r0 // 2
    x = seed % half + 1                        # an exponent in the "x < r0/2" class
    neg = r0 - x                               # its negation: same x, other sign
    assert 0 < x <= half and neg > half

    # coin-by-coin: under x-only keys the pair is the same hash input for every
    # coin (by construction), under compressed keys it is an ordinary distinct pair
    xonly_forced = 0
    signed_hits = 0
    coins = 64
    for i in range(coins):
        aa = ((seed + i) * 2654435761 + 1) | 1
        bb = ((seed + i) * 40503 + 7) % (1 << w)
        px = (aa * x + bb) % (1 << w)
        pn = (aa * neg + bb) % (1 << w)
        kx = min(x, neg)
        kn = min(neg, x)
        assert kx == kn
        p_forced = (aa * kx + bb) % (1 << w)
        if p_forced == (aa * kn + bb) % (1 << w):
            xonly_forced += 1
        if (px >> (w - l), (px >> (w - l - tau)) & ((1 << tau) - 1)) == (
                pn >> (w - l), (pn >> (w - l - tau)) & ((1 << tau) - 1)):
            signed_hits += 1

    t1 = Tab()
    t1.insert(min(x, neg), x)                  # baby node, x-only key
    hit1 = t1.lookup(min(neg, x), neg)         # giant = -baby, full compare on exponents
    t2 = Tab()
    t2.insert(x, x)                            # compressed key
    hit2 = t2.lookup(neg, neg)

    t3 = Tab()
    t3.insert(min(x, neg), x)
    hit3 = t3.lookup(min(x, neg), x)           # the TRUE match, x-only keys
    t4 = Tab()
    t4.insert(x, x)
    hit4 = t4.lookup(x, x)                     # the TRUE match, compressed keys

    ok = (
        xonly_forced == coins                  # forced, not probabilistic
        and hit1 is False and t1.false_cmp == 1
        and hit2 is False and t2.false_cmp == 0 and t2.cmp == 0
        and hit3 is True and t3.false_cmp == 0
        and hit4 is True and t4.false_cmp == 0
    )
    return {
        "sign_collision": {
            "coins": coins,
            "xonly_forced_tag_matches": xonly_forced,
            "xonly_comparisons_on_sign_pair": t1.cmp,
            "xonly_false_comparisons_on_sign_pair": t1.false_cmp,
            "xonly_lookup_returned_hit": hit1,
            "compressed_tag_matches_on_sign_pair": signed_hits,
            "compressed_comparisons_on_sign_pair": t2.cmp,
            "compressed_false_comparisons_on_sign_pair": t2.false_cmp,
            "compressed_lookup_returned_hit": hit2,
            "xonly_true_match_found": hit3,
            "compressed_true_match_found": hit4,
            "true_match_false_comparisons": t3.false_cmp + t4.false_cmp,
            "counter_ok": ok,
        },
    }



def main():
    params = json.load(sys.stdin)
    trials = int(params["trials"])
    seed = int(params["seed"])
    d_bits = int(params["d_bits"])
    if not 1 <= trials <= 8:
        raise ValueError("trials out of range")
    if not 0 <= seed <= 2147483647:
        raise ValueError("seed out of range")
    if not 8 <= d_bits <= 12:
        raise ValueError("d_bits out of range")

    r0, q = toy_prime(seed, d_bits)
    zeta = primitive_root(r0)
    d0 = 1 << d_bits
    rng = random.Random(seed)
    last_st = {}
    secrets = [1, r0 - 1, r0 - 2, 2]
    while len(secrets) < trials:
        secrets.append(rng.randrange(1, r0))
    recovered = 0
    for s in secrets[:trials]:
        ok, k0, k1, st = cheon_toy(r0, zeta, d0, s, seed)
        assert ok, (r0, d0, s)
        recovered += 1
        last_st = st

    toy = {
        "prime": str(r0),
        "cofactor_q": str(q),
        "d_bits": d_bits,
        "trials": trials,
        "recovered": recovered,
        "chained_table": last_st,
        "link_budget_aborts": int(last_st["abort"]),
        "window_bits": max(1, math.ceil(r0.bit_length() / WINDOWS)),
    }
    prod = production_checks()
    prod.update(sign_and_comparison_checks(r0, seed))
    print(json.dumps({
        "toy": toy,
        "production": prod,
        "scope": "toy algebra and parameter/ledger arithmetic only; no production-size run",
    }, sort_keys=True))


if __name__ == "__main__":
    main()

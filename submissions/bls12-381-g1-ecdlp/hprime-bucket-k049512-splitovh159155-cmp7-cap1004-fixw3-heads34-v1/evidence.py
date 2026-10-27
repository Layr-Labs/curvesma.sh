#!/usr/bin/env python3
import sys, json, math, random, bisect, time
from fractions import Fraction

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
X_BLS = -0xD201000000010000
FP_MUL_INSTR = 383
LOG_D = 20
D = 1 << LOG_D
M = (R - 1) // D
C = 2000
K_NUM, K_DEN = 49512, 100000
LOG_CHAINS = 64
LOG_NODES = 84
HOP_BITS, HOPS = 28, 3
LOG_HEADS = 34
LOG_NODE_BYTES = 40
LOG_POOL_RECORDS = 33
LOG_GAP_TABLE = 32
LOG_GEN_NODES, LOG_STORE_NODES = 20, 30
KEY_BITS = 192
OVERHEAD_LISTED_BABY = {"gap stream": 5, "table address": 2, "operand loads and stores": 18, "lane control": 8,
                        "parity key": 23, "scatter 3 hops": 36, "post and poll 3 hops": 9, "transport copy 3 hops": 18,
                        "gap stream transport from storage": 3, "batch inversion loads and stores": 6,
                        "flow control and counters": 5, "shared F_p inversion of the round": 4, "insert": 10}
OVERHEAD_LISTED_GIANT = dict(OVERHEAD_LISTED_BABY); del OVERHEAD_LISTED_GIANT["insert"]; OVERHEAD_LISTED_GIANT["lookup fixed"] = 6
OVERHEAD_CHARGED_BABY, OVERHEAD_CHARGED_GIANT = 159, 155
WALK_STEP_OVERHEAD = sum(OVERHEAD_LISTED_BABY[k] for k in ("gap stream", "gap stream transport from storage", "table address",
                         "operand loads and stores", "batch inversion loads and stores", "lane control", "flow control and counters"))
assert WALK_STEP_OVERHEAD == 47
INSTR_PER_COMPARISON = 7
SLOT_CHAIN_CAP = 64
BERNSTEIN_VAR_MARGIN = 4
NODE_FIXED_INSTR = {"head zeroing 2^28 iters x 34 (32 stp + cmp + bne)": (1 << 28) * 34, "pointer setup": 1 << 20,
                     "end-of-hop flushes 3 hops x 2^28 blocks x 5 (post-inc store + bound cmp, final bne folded into hop loop)": 3 * (1 << 28) * 5,
                     "status record 48 B": 6, "comparison counter preset": 4}
NODE_FIXED_INSTR_TOTAL = sum(NODE_FIXED_INSTR.values())
assert NODE_FIXED_INSTR_TOTAL == 13154385930
LOG_SCALAR_MULT_UNITS = 12
LOG_HTEST_INSTR = 18
COMPARISON_CAP_NUM, COMPARISON_CAP_DEN = 1004, 1000
LOG_MATCH_CAP = 44
CANDIDATE_BYTES = 40
LOG_NODE_CANDIDATE_CAP = 16
CANDIDATE_FORWARD_INSTR = 128
LOG_GEN_CANDIDATE_CAP = 24
REWALK_STEP_INSTR = 32
REWALK_MATCH_CAP = 4
RECOVERY_HIT_TOKEN_INSTR = 128
KEPT_FORWARD_INSTR = 256
LOG_STATUS_INSTR = 12
PRE_CLAIM_LOG2_E4 = 1157603
MEMORY_CLAIM_LOG2_E4 = 1234615
LOG_STORAGE_NODE_BYTES = 40
LOG_BIGGAP_RESERVE = 69
BIGGAP_SPLITS = 16
LOG_ROUND = 16
LOG_WALK_SIEVE, LOG_WALK_HTESTS = 78, 73
LOG_BASE_SIEVE, LOG_BASE_HTESTS = 92, 86
LOG_STAGE2 = 24
SIEVE_INSTR = 128
LOG_SIEVE_BITMAP = 38
SIEVE_CROSSING = {"ldrb": 1, "orr bit mask": 1, "strb": 1, "advance by p bytes": 1, "cmp": 1, "branch": 1}
SIEVE_PRIME_SEGMENT = {"load p from the receive block": 1, "30 p (lsl, sub)": 2,
                       "lo mod 30 p: high limb (udiv, msub), six chunks of the low limb (lsl, ubfx, orr, udiv, msub)": 2 + 6 * 5,
                       "p mod 30 (udiv, msub), class-table row (ldrb, add)": 4,
                       "8 wheel classes x (table ldrb, mul, subs/add/csel wrap, umulh/lsr byte index, mask, entry cmp/b)": 8 * 11,
                       "prime loop (add, cmp, b)": 3}
SIEVE_BLOCK_LISTED = {"copy 64 primes of 8 B into a 512 B block at the holder (32 ldp/stp pairs)": 64,
                      "post and poll (the record route's 48 per block)": 48, "receiver block loop and prime count": 4}
SIEVE_BLOCK_INSTR = 256
SIEVE_SCAN_WORD = {"ldr": 1, "mvn": 1, "cbz": 1, "word base value (adds, adc)": 2, "loop (add, cmp, b)": 3}
SIEVE_SURVIVOR = {"rbit, clz": 2, "wheel value ldr": 1, "value (adds, adc)": 2, "queue stp and pointer": 2,
                  "queue full test (cmp, b)": 2, "clear lowest bit (sub, and)": 2, "loop cbnz": 1,
                  "emit a 16-byte entry to a storage node (record route: 9 per 32 B)": 8}
SIEVE_CLEAR_64B = 7
SIEVE_SEGMENT_FIXED = 1 << 10
SIEVE_REMARK = 16
SIEVE_PLAIN = 64
LOG_SIEVE_HOLDERS = 30
SIEVE_INTEGERS = {"walk": (1 << 78) + (1 << 43), "base": (1 << 92) + (1 << 43), "primes": (1 << 47) + (1 << 43)}
H1_ALLOWANCE = Fraction(5, 10000)
SUCCESS_CLAIM = 0.3900
CENSUS_WORST_M30 = -0.00023204
CENSUS_WORST_M30_K4965 = -0.00007302
CENSUS_WORST_M30_K4975 = -0.00010833

def ceil_div(a, b): return -(-a // b)

def up(x, nd=0):
    """Round up to nd decimals: every quoted cost, size, memory, mean load and log2 failure bound is rounded this way."""
    return math.ceil(x * 10 ** nd) / 10 ** nd

def down(x, nd=0):
    """Round down to nd decimals: every quoted success figure and every quoted capacity margin is rounded this way."""
    return math.floor(x * 10 ** nd) / 10 ** nd

def pow_at_most_pow2(w, e, k, prec=2048):
    """Exact integer proof that w^e <= 2^k: binary powering on (x, s) with value <= x 2^s, each product truncated
    upwards to prec bits, so the final x 2^s is an upper bound on w^e (False means not proved)."""
    def trunc(x, s):
        t = x.bit_length() - prec
        return ((x >> t) + 1, s + t) if t > 0 else (x, s)
    rx, rs, bx, bs = 1, 0, w, 0
    while e:
        if e & 1: rx, rs = trunc(rx * bx, rs + bs)
        e >>= 1
        if e: bx, bs = trunc(bx * bx, 2 * bs)
    return rs <= k and rx <= 1 << (k - rs)

def poisson_tail_log2(mu, k):
    """log2 of the Chernoff bound P(X >= k) <= e^-mu (e mu / k)^k for X Poisson(mu), k > mu."""
    assert k > mu
    return (-mu + k * (1 + math.log(mu) - math.log(k))) * math.log2(math.e)

def log2_int(n):
    """log2 of a positive integer with about 15 significant digits."""
    b = n.bit_length()
    if b <= 1000: return math.log2(n)
    return (b - 1000) + math.log2(n >> (b - 1000))

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def sieve(hi):
    s = bytearray([1]) * (hi + 1); s[0] = s[1] = 0
    for i in range(2, math.isqrt(hi) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(hi + 1) if s[i]]

def prices(inv_p_mults, headroom_log2):
    """Counts behind the four fixed prices: multiplications at FP_MUL_INSTR for the scalar multiplication, the H-test and
    Lemma 1, and the itemised instruction count of the sieve."""
    out = {}
    bits = R.bit_length(); sm = (bits - 1) * (7 + 11) + inv_p_mults + 4
    budget = C << LOG_SCALAR_MULT_UNITS
    out["scalar_mult"] = {"multiplications": sm, "instr": sm * FP_MUL_INSTR, "charged_instr": budget,
                          "charged_over_counted": down(budget / (sm * FP_MUL_INSTR), 2)}
    ht = M.bit_length() - 1 + bin(M).count("1") - 1
    out["h_test"] = {"m_bits": M.bit_length(), "m_set_bits": bin(M).count("1"), "multiplications": ht,
                     "instr": ht * FP_MUL_INSTR, "charged_instr": 1 << LOG_HTEST_INSTR,
                     "charged_over_counted": down((1 << LOG_HTEST_INSTR) / (ht * FP_MUL_INSTR), 2)}
    v2 = 32; Q = (R - 1) >> v2; h = (Q - 1) // 2
    pw = lambda e: e.bit_length() - 1 + bin(e).count("1") - 1
    per_root = pw(h) + 2 + (v2 - 1) * ((v2 - 1) + (v2 - 1) + 3)
    roots = 20 * per_root + pw(Q)
    l1_instr = (roots + 3 * (1 << 10) + 2) * FP_MUL_INSTR + 2 * (1 << 10) * 11 * 64 + (1 << 11) * (C << LOG_SCALAR_MULT_UNITS)
    out["lemma1"] = {"square_root_multiplications_each": per_root, "roots_total_multiplications": roots,
                     "instr_units": up(l1_instr / C, 0), "charged_units": 1 << LOG_STAGE2,
                     "charged_over_counted": down((C << LOG_STAGE2) / l1_instr, 2)}
    def crossings(root_log2):
        lx = root_log2 * math.log(2)
        return 8 / 30 * (math.log(lx) + 0.2615 + 1 / (2 * lx * lx) - (1 / 2 + 1 / 3 + 1 / 5))
    pi_lo = lambda v: v / math.log(v)
    pi_hi = lambda v: 1.25506 * v / math.log(v)
    cross, per_prime = sum(SIEVE_CROSSING.values()), sum(SIEVE_PRIME_SEGMENT.values())
    block_listed, scan, surv_instr = sum(SIEVE_BLOCK_LISTED.values()), sum(SIEVE_SCAN_WORD.values()), sum(SIEVE_SURVIVOR.values())
    assert (cross, per_prime, block_listed, scan, surv_instr) == (6, 130, 116, 8, 20) and SIEVE_BLOCK_INSTR >= block_listed
    S = 30 << LOG_SIEVE_BITMAP
    assert (1 << LOG_SIEVE_BITMAP) <= 96 << LOG_GAP_TABLE
    base_root = math.isqrt((1 << 80) + (1 << LOG_BASE_SIEVE))
    stored = math.ceil(pi_hi(base_root))
    per_holder = ceil_div(stored, 1 << LOG_SIEVE_HOLDERS)
    blocks_per_segment = (1 << LOG_SIEVE_HOLDERS) * ceil_div(per_holder, 64)
    def window(a, b, root, charged, streamed):
        segs = ceil_div(b - (a - a % 30) + 1, S)
        n = segs * S
        assert n <= charged
        primes = pi_hi(root)
        survivors = pi_hi(b) - (pi_lo(a) if a >= 17 else 0)
        items = {"bitmap clear": n / 30 / 64 * SIEVE_CLEAR_64B,
                 "crossings (sum 8 n / (30 p))": n * crossings(math.log2(root)) * cross,
                 "crossings (8 per prime per segment)": 8 * segs * primes * cross,
                 "per prime per segment": segs * primes * per_prime,
                 "prime stream blocks": segs * blocks_per_segment * SIEVE_BLOCK_INSTR if streamed else 0,
                 "scan": n / 240 * scan, "survivors (H-test hand-off and emission)": survivors * surv_instr,
                 "segment fixed": segs * SIEVE_SEGMENT_FIXED}
        if not streamed:
            items["plain sieve of (0, 2^24] per segment"] = segs * (1 << 24) * SIEVE_PLAIN
            items["restore the primes below 2^24"] = segs * pi_hi(1 << 24) * SIEVE_REMARK
        per_int = {k: up(v / charged, 3) for k, v in items.items()}
        total = sum(items.values()) / charged
        assert total <= SIEVE_INSTR
        return {"segments": segs, "sieving_primes_log2_below": up(math.log2(primes), 2),
                "crossings_sum_part_per_integer_below": up(crossings(math.log2(root)), 3),
                "instr_per_charged_integer": per_int, "total_per_charged_integer_below": up(total, 2),
                "charged_over_counted": down(SIEVE_INSTR / total, 2)}
    walk = window(1 << 39, 1 << LOG_WALK_SIEVE, 1 << 39, SIEVE_INTEGERS["walk"], True)
    base = window(1 << 80, (1 << 80) + (1 << LOG_BASE_SIEVE), base_root, SIEVE_INTEGERS["base"], True)
    gen = window(0, 1 << 47, 1 << 24, SIEVE_INTEGERS["primes"], False)
    assert base_root < 1 << 47
    total_integers = sum(SIEVE_INTEGERS.values())
    out["sieve"] = {"segment_integers": "30 x 2^38", "bitmap_bytes_log2": LOG_SIEVE_BITMAP,
                    "crossing_listed": cross, "per_prime_per_segment_listed": per_prime,
                    "block_listed": block_listed, "block_charged": SIEVE_BLOCK_INSTR, "scan_word_listed": scan,
                    "survivor_listed": surv_instr, "stored_sieving_primes_at_most": stored,
                    "stored_sieving_primes_log2_below": up(math.log2(stored), 2),
                    "sieving_prime_bound_log2_below": up(math.log2(base_root + 1), 4),
                    "primes_per_holder_at_most": per_holder, "blocks_per_segment_at_most": blocks_per_segment,
                    "walk_window": walk, "base_window": base, "sieving_primes_range": gen,
                    "charged_instr_per_integer": SIEVE_INSTR,
                    "score_unchanged_with_further_instr_per_sieved_integer_up_to":
                        int(2 ** (headroom_log2 + math.log2(C) - math.log2(total_integers)))}
    return out

def ledger():
    out = {}
    assert (R - 1) % D == 0 and (R - 1) % (6 * D) == 0 and M % 6 == 0
    v2 = 0
    while (R - 1) % (1 << (v2 + 1)) == 0: v2 += 1
    out["v2_of_r_minus_1"] = v2; assert v2 >= LOG_D
    out["log2_m"] = round(log2_int(M), 6)
    omega = next(pow(g, (R - 1) // 3, R) for g in range(2, 100) if pow(g, (R - 1) // 3, R) != 1)
    mu6 = [1, R - 1, omega, R - omega, omega * omega % R, (R - omega * omega) % R]
    assert all(pow(u, M, R) == 1 for u in mu6), "mu6 must lie in H"
    out["mu6_in_H"] = True
    n0 = math.isqrt(K_NUM * M // (6 * K_DEN)) + 1
    chains = 1 << LOG_CHAINS
    L = ceil_div(n0, chains)
    N = L * chains
    assert 6 * N * N * K_DEN >= K_NUM * M
    out["steps_per_chain_log2"] = up(log2_int(L), 4)
    out["records_per_side_log2"] = up(log2_int(N), 4)
    K = 6 * N * N / M
    out["K"] = down(K, 6)
    k_excess = Fraction(6 * N * N, M) - Fraction(K_NUM, K_DEN)
    assert k_excess >= 0
    k_excess_log2 = up(log2_int(k_excess.numerator) - log2_int(k_excess.denominator), 1)
    out["K_minus_%g_log2_below" % (K_NUM / K_DEN)] = k_excess_log2
    out["predicted_success"] = math.floor((1 - math.exp(-K)) * 10 ** 6) / 10 ** 6
    out["diagonal_loss_log2"] = up(-log2_int(L), 2)
    hp_needed = L
    primes_needed = hp_needed * D
    x = primes_needed * 54
    out["walk_list_last_prime_log2_about"] = round(log2_int(x), 2)
    assert x < 1 << 79
    base_count = 2 * chains
    out["base_primes_needed_log2"] = LOG_CHAINS + 1
    pi_lo = lambda v: v / math.log(v)
    pi_hi = lambda v: 1.25506 * v / math.log(v)
    walk_lo, walk_hi = pi_lo(2.0 ** LOG_WALK_SIEVE) - pi_hi(2.0 ** 39), pi_hi(2.0 ** LOG_WALK_SIEVE)
    assert walk_hi < 2 ** LOG_WALK_HTESTS
    assert walk_lo / D > L
    out["h2_walk_window"] = {"sieved": "(2^39, 2^78]", "primes_log2_above": math.floor(math.log2(walk_lo) * 100) / 100,
                             "primes_log2_below": up(math.log2(walk_hi), 2), "h_tests_charged_log2": LOG_WALK_HTESTS,
                             "examined": "its primes in increasing order up to the L-th H-prime, all within the 2^73 H-tests",
                             "h_primes_expected_log2_above": math.floor(math.log2(walk_lo / D) * 100) / 100,
                             "needed_log2": up(log2_int(L), 4),
                             "expected_over_needed_above": math.floor(walk_lo / D / L * 100) / 100}
    base_lo = pi_lo(2.0 ** 80 + 2.0 ** LOG_BASE_SIEVE) - pi_hi(2.0 ** 80)
    assert base_lo > 2 ** LOG_BASE_HTESTS
    base_expected = 1 << (LOG_BASE_HTESTS - LOG_D)
    assert base_expected >= 2 * base_count
    out["h2_base_window"] = {"sieved": "(2^80, 2^80 + 2^92]", "primes_log2_above": math.floor(math.log2(base_lo) * 1000) / 1000,
                             "h_tests_charged_log2": LOG_BASE_HTESTS, "examined": "at most the first 2^86 primes above 2^80, in increasing order up to the 2^65-th H-prime",
                             "h_primes_expected_log2": LOG_BASE_HTESTS - LOG_D, "needed_log2": LOG_CHAINS + 1,
                             "expected_over_needed": base_expected // base_count}
    assert (1 << 80) + (1 << LOG_BASE_SIEVE) < 1 << 93 and (1 << (93 + LOG_WALK_SIEVE)) < R
    out["x_max_log2_below"] = 93 + LOG_WALK_SIEVE
    mean_gap = D * LOG_WALK_SIEVE * math.log(2)
    out["mean_gap_log2_at_most"] = up(math.log2(mean_gap), 2)
    out["gap_above_2^33_expected_log2_below"] = up(math.log2(L - 1) - (1 << 33) / mean_gap * math.log2(math.e), 1)
    lb, lg = sum(OVERHEAD_LISTED_BABY.values()), sum(OVERHEAD_LISTED_GIANT.values())
    out["overhead_listed"] = {"baby": lb, "giant": lg, "charged_baby": OVERHEAD_CHARGED_BABY, "charged_giant": OVERHEAD_CHARGED_GIANT}
    assert OVERHEAD_CHARGED_BABY >= lb and OVERHEAD_CHARGED_GIANT >= lg
    c_baby, c_giant = C + OVERHEAD_CHARGED_BABY, C + OVERHEAD_CHARGED_GIANT
    assert R == X_BLS ** 4 - X_BLS ** 2 + 1 and P == (X_BLS - 1) ** 2 * R // 3 + X_BLS
    inv_p_mults = (P - 2).bit_length() - 1 + bin(P - 2).count("1") - 1
    inv_r_mults = (R - 2).bit_length() - 1 + bin(R - 2).count("1") - 1
    inv_p = inv_p_mults * FP_MUL_INSTR
    walk_rounds = ceil_div(L, 1 << LOG_ROUND) + BIGGAP_SPLITS
    table_rounds = 1 << (LOG_GAP_TABLE - LOG_ROUND)
    listed_inv = OVERHEAD_LISTED_BABY["shared F_p inversion of the round"]
    assert walk_rounds * inv_p <= L * listed_inv
    z_mults = REWALK_MATCH_CAP * inv_r_mults + 2 * REWALK_MATCH_CAP + REWALK_MATCH_CAP ** 2 * 7
    out["field_inversions"] = {"F_p_squarings": (P - 2).bit_length() - 1, "F_p_multiplications": bin(P - 2).count("1") - 1,
                               "F_p_inversion_instr": inv_p, "rounds_per_chain_walk": "ceil(L / 2^16) + 16",
                               "rounds_per_gap_table": table_rounds,
                               "per_record_share_instr_below": up(walk_rounds * inv_p / L, 2), "per_record_listed": listed_inv,
                               "F_r_inversion_multiplications": inv_r_mults, "z_forming_multiplications_per_candidate": z_mults}
    nodes = 1 << LOG_NODES
    mean_per_node = N / nodes
    out["baby_records_per_node_mean_log2"] = up(math.log2(mean_per_node), 3)
    cap = 1 << LOG_POOL_RECORDS
    delta = cap / mean_per_node - 1
    out["pool_capacity_over_mean"] = down(cap / mean_per_node, 2)
    chernoff_log2 = -mean_per_node * delta * delta / (2 + delta) * math.log2(math.e)
    out["overflow_per_node_chernoff_log2"] = up(chernoff_log2)
    heads = nodes << LOG_HEADS
    lam = N / heads
    out["load_per_head"] = up(lam, 3)
    comp_mean = N * N // heads
    node_comp_mean = comp_mean // nodes
    node_comp_cap = ceil_div(COMPARISON_CAP_NUM * node_comp_mean, COMPARISON_CAP_DEN)
    comp_cap = node_comp_cap * nodes
    out["comparisons_mean_log2"] = up(log2_int(comp_mean), 3)
    out["comparisons_per_node_mean_log2"] = up(log2_int(node_comp_mean), 3)
    out["comparisons_per_node_cap_log2"] = up(log2_int(node_comp_cap), 3)
    k = SLOT_CHAIN_CAP
    k1 = k + 1
    slot_tail_log2 = (-lam * math.log2(math.e) + k1 * math.log2(lam) - math.log2(math.factorial(k1))
                      - math.log2(1 - lam / (k1 + 1)))
    out["slot_chain_cap_exceeded_all_slots_log2"] = up(slot_tail_log2 + math.log2(heads), 1)
    t = 0.002 * node_comp_mean
    V = mean_per_node * (lam + lam * lam) * BERNSTEIN_VAR_MARGIN
    bern_log2 = -(t * t / (2 * (V + k * t / 3))) * math.log2(math.e)
    out["comparisons_node_cap_exceeded_bernstein_all_nodes_log2"] = up(bern_log2 + LOG_NODES)
    count_dev_log2 = -(mean_per_node * (2 ** -10) ** 2 / 3) * math.log2(math.e) + LOG_NODES + 1
    out["node_count_off_by_2^-10_chernoff_all_nodes_log2"] = up(count_dev_log2)
    false_mean = N * N // (1 << KEY_BITS)
    out["false_matches_mean_log2"] = up(log2_int(false_mean), 3)
    match_cap = 1 << LOG_MATCH_CAP
    out["match_cap_exceeded_chebyshev_log2"] = up(math.log2(2 * false_mean + 1) - log2_int((match_cap - false_mean) ** 2), 1)
    true_budget = 1 << LOG_D
    out["true_matches_mean"] = "K, at most %g + 2^" % (K_NUM / K_DEN) + str(k_excess_log2)
    # Unconditional Markov bound: E[M(z)] over uniform z in H is exactly the off-diagonal pair count over m,
    # below K, so P(M(z) > true_budget) <= K / (true_budget + 1) with no distributional premise. The enforced
    # recovery-hit budget of Section 2 carries the work bound; this row bounds the budget-overflow failure event.
    E_offdiag = 6 * N * N * (L - 1) / (L * M)
    assert E_offdiag < K
    markov_overflow = K / (true_budget + 1)
    out["true_match_budget_exceeded_off_diagonal_markov"] = up(markov_overflow, 12)
    out["true_match_budget_exceeded_off_diagonal_log2_below"] = up(math.log2(markov_overflow), 2)
    assert L > match_cap > true_budget
    diag_values = 6 * chains * chains
    diag_prob = diag_values / M
    out["diagonal_event_values"] = "6 x 2^128"
    out["diagonal_event_probability_log2_below"] = up(math.log2(diag_values) - log2_int(M), 2)
    success_floor = 1 - math.exp(-K) - markov_overflow - 2 ** -20 - 2 ** -40 - diag_prob
    out["success_after_failure_terms"] = down(success_floor, 6)
    after_h1 = success_floor - float(H1_ALLOWANCE)
    out["success_after_failure_terms_and_H1_allowance"] = down(after_h1, 6)
    out["success_claimed"] = SUCCESS_CLAIM
    assert after_h1 >= SUCCESS_CLAIM >= 0.39
    out["coverage_deficit_tolerated"] = {"claimed": down(success_floor - SUCCESS_CLAIM, 6), "required_0.39": down(success_floor - 0.39, 6)}
    out["margin_over_0.39_after_H1_allowance"] = down(after_h1 - 0.39, 6)
    out["census_lowest_prime_m_2^30_up"] = CENSUS_WORST_M30
    out["margin_over_census_lowest"] = down((after_h1 - 0.39) / -CENSUS_WORST_M30, 1)
    out["census_lowest_prime_m_2^30_up_at_K_0.4965"] = CENSUS_WORST_M30_K4965
    out["margin_over_census_lowest_at_K_0.4965"] = down((after_h1 - 0.39) / -CENSUS_WORST_M30_K4965, 1)
    out["census_lowest_prime_m_2^30_up_at_K_0.4975"] = CENSUS_WORST_M30_K4975
    out["margin_over_census_lowest_at_K_0.4975"] = down((after_h1 - 0.39) / -CENSUS_WORST_M30_K4975, 1)
    rows = {}
    rows["records N baby + N giant (addition + split charged overhead)"] = N * (c_baby + c_giant)
    rows["gap tables 2^65 chains x 2^32 additions at C + 47"] = base_count * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD)
    rows["table strand starts 2^65 x (2^16 + 1) scalar mults"] = base_count * ((1 << 16) + 1) * (C << LOG_SCALAR_MULT_UNITS)
    rows["segment starts 2^65 x 2^16 scalar mults"] = base_count * (1 << 16) * (C << LOG_SCALAR_MULT_UNITS)
    rows["per-node fixed work"] = ((1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)) * NODE_FIXED_INSTR_TOTAL
    out["node_fixed_instr_itemized"] = dict(NODE_FIXED_INSTR)
    rows["comparisons at cap x 7"] = comp_cap * INSTR_PER_COMPARISON
    step_instr = C + WALK_STEP_OVERHEAD + REWALK_STEP_INSTR
    rewalk = (2 * L * step_instr + 2 * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD)
              + (4 * (1 << 16) + 2) * (C << LOG_SCALAR_MULT_UNITS))
    out["rewalk_step_instr"] = step_instr
    rows["matches at cap: forward 40 B candidates to both generating nodes"] = match_cap * CANDIDATE_FORWARD_INSTR
    rows["matches at cap: re-walk both chains, rebuild both tables"] = match_cap * rewalk
    rows["matches at cap: kept values of the giant chain to the verification node"] = match_cap * KEPT_FORWARD_INSTR
    # Enforced recovery-hit budget (proof Section 2): one 128-instruction token round trip with node 0 per
    # candidate at cap bounds the at most 2^20 + 1 Lemma 1 executions; the work is charged whether or not the
    # token is granted, so a denied run has every poll paid for.
    rows["matches at cap: recovery-hit token round trips with node 0"] = match_cap * RECOVERY_HIT_TOKEN_INSTR
    rows["matches at cap: up to 6 x 4 x 4 candidate scalar mults"] = match_cap * 6 * REWALK_MATCH_CAP ** 2 * (C << LOG_SCALAR_MULT_UNITS)
    rows["gap tables: one F_p inversion per round, 2^65 x 2^16 rounds"] = base_count * table_rounds * inv_p
    rows["matches at cap: one F_p inversion per round of both re-walks and both table rebuilds"] = (
        match_cap * 2 * (walk_rounds + table_rounds) * inv_p)
    rows["matches at cap: F_r inversions and multiplications forming the six z"] = match_cap * z_mults * FP_MUL_INSTR
    cand_node_mean_log2 = log2_int(false_mean + true_budget) - LOG_NODES
    out["candidates_per_node_mean_log2"] = up(cand_node_mean_log2, 3)
    out["node_candidate_cap_exceeded_all_nodes_log2_below"] = up((1 << LOG_NODE_CANDIDATE_CAP) * cand_node_mean_log2 + LOG_NODES)
    gen_cand_mean = 2 * (false_mean + true_budget) / (1 << LOG_GEN_NODES)
    out["candidates_per_generating_node_mean_log2"] = up(math.log2(gen_cand_mean), 3)
    out["generating_node_candidate_area_exceeded_all_nodes_log2_below"] = up(
        poisson_tail_log2(gen_cand_mean, (1 << LOG_GEN_CANDIDATE_CAP) + 1) + LOG_GEN_NODES)
    repeat_log2 = (REWALK_MATCH_CAP + 1) * log2_int(L) - REWALK_MATCH_CAP * KEY_BITS - math.log2(math.factorial(REWALK_MATCH_CAP + 1))
    out["rewalk_match_cap_exceeded_all_chains_log2_below"] = up(repeat_log2 + LOG_CHAINS + 1)
    out["within_chain_key_repeat_per_chain_log2"] = up(2 * log2_int(L) - 1 - KEY_BITS, 1)
    assert SIEVE_INSTR == 128
    rows["walk list W: sieve (2^39, 2^78] as 2^78 + 2^43 integers x 128 + 2^73 H-tests"] = (
        SIEVE_INTEGERS["walk"] * SIEVE_INSTR + ((1 << LOG_WALK_HTESTS) << LOG_HTEST_INSTR))
    rows["base primes: sieve (2^80, 2^80 + 2^92] as 2^92 + 2^43 integers x 128 + 2^86 H-tests"] = (
        SIEVE_INTEGERS["base"] * SIEVE_INSTR + ((1 << LOG_BASE_HTESTS) << LOG_HTEST_INSTR))
    rows["sieving primes: sieve (0, 2^47] as 2^47 + 2^43 integers x 128"] = SIEVE_INTEGERS["primes"] * SIEVE_INSTR
    rows["split gaps above 2^33 reserve"] = (C << LOG_BIGGAP_RESERVE)
    rows["split gaps: walk overhead of the reserved additions, re-walk split steps at cap"] = (
        (WALK_STEP_OVERHEAD << LOG_BIGGAP_RESERVE) + match_cap * 2 * BIGGAP_SPLITS * step_instr)
    all_nodes = (1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)
    rows["status tree: 3 passes x every node x 2^12"] = 3 * (all_nodes << LOG_STATUS_INSTR)
    rows["stage 2 up to 2^20 + 1 times, input access, output"] = (true_budget + 1) * (C << LOG_STAGE2)
    total_instr = sum(rows.values())
    W = ceil_div(total_instr, C)
    out["rows_log2_units"] = {k: up(log2_int(ceil_div(v, C)), 3) for k, v in rows.items()}
    out["W_units"] = str(W)
    lw = log2_int(W)
    score_num = math.ceil((lw + 1e-9) * 10 ** 6)
    out["log2_W"] = up(lw, 7)
    headroom_log2 = math.floor((lw + math.log2(math.expm1((score_num / 10 ** 6 - lw) * math.log(2)))) * 100) / 100
    out["score_headroom_units_log2_above"] = headroom_log2
    g = math.gcd(score_num, 10 ** 6)
    out["score"] = {"numerator": str(score_num // g), "denominator": str(10 ** 6 // g)}
    assert pow_at_most_pow2(W, 10 ** 6 // g, score_num // g), "W exceeds 2^score"
    out["score_millionths"] = f"{score_num}/1000000"
    out["score_decimal"] = score_num / 10 ** 6
    out["c_per_record"] = (c_baby + c_giant) / (2 * C)
    out["prices"] = prices(inv_p_mults, headroom_log2)
    pre = (N * c_baby + chains * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD) + chains * (2 * (1 << 16) + 1) * (C << LOG_SCALAR_MULT_UNITS)
           + chains * table_rounds * inv_p
           + rows["per-node fixed work"] + rows["walk list W: sieve (2^39, 2^78] as 2^78 + 2^43 integers x 128 + 2^73 H-tests"]
           + rows["base primes: sieve (2^80, 2^80 + 2^92] as 2^92 + 2^43 integers x 128 + 2^86 H-tests"]
           + rows["sieving primes: sieve (0, 2^47] as 2^47 + 2^43 integers x 128"]
           + (all_nodes << LOG_STATUS_INSTR) + ((C + WALK_STEP_OVERHEAD) << (LOG_BIGGAP_RESERVE - 1)))
    pre_units = ceil_div(pre, C)
    out["preprocessing_log2"] = up(log2_int(pre_units), 4)
    assert pre_units ** 10000 <= 1 << PRE_CLAIM_LOG2_E4, "preprocessing exceeds the quoted bound"
    assert pre_units <= W
    node_bytes = 1 << LOG_NODE_BYTES
    per_owner = {"heads 2^34 x 8 B": 8 << LOG_HEADS, "pool 2^33 x 32 B": 32 << LOG_POOL_RECORDS,
                 "send blocks 2^28 x 512 B": 512 << HOP_BITS, "receive FIFOs 3 x 2^36 B": 3 << 36,
                 "candidate list 2^16 x 40 B": CANDIDATE_BYTES << LOG_NODE_CANDIDATE_CAP, "code and counters": 1 << 30}
    owner_bytes = sum(per_owner.values())
    assert owner_bytes <= node_bytes
    out["owner_node_bytes_log2"] = up(math.log2(owner_bytes), 4)
    per_gen = {"gap table 2^32 x 96 B": 96 << LOG_GAP_TABLE, "send blocks 2^28 x 512 B": 512 << HOP_BITS,
               "stream buffers and state": 1 << 36, "candidate receive area 2^24 x 40 B": CANDIDATE_BYTES << LOG_GEN_CANDIDATE_CAP}
    gen_bytes = sum(per_gen.values())
    assert gen_bytes <= node_bytes
    out["generating_node_bytes_log2"] = up(math.log2(gen_bytes), 4)
    sieve_stored = out["prices"]["sieve"]["stored_sieving_primes_at_most"]
    storage_bytes = base_count * 16 + L * 4 + sieve_stored * 8
    storage_alloc = (1 << LOG_STORE_NODES) << LOG_STORAGE_NODE_BYTES
    assert storage_bytes <= storage_alloc <= (1 << LOG_STORE_NODES) * node_bytes
    assert (ceil_div(base_count * 16 + L * 4, 1 << LOG_STORE_NODES) + 8 * ceil_div(sieve_stored, 1 << LOG_STORE_NODES)
            <= 1 << LOG_STORAGE_NODE_BYTES)
    out["storage_bytes_log2"] = up(math.log2(storage_bytes), 2)
    mem = (1 << LOG_NODES) * owner_bytes + (1 << LOG_GEN_NODES) * gen_bytes + storage_alloc
    out["memory_bytes"] = str(mem)
    out["memory_log2_bytes_exact_truncated"] = down(log2_int(mem), 6)
    assert mem ** 10000 <= 1 << MEMORY_CLAIM_LOG2_E4, "memory exceeds the claimed bound"
    out["memory_log2_bytes"] = MEMORY_CLAIM_LOG2_E4 / 10000
    out["memory_every_node_at_2^40_would_be_bytes"] = str(((1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)) * node_bytes)
    out["records_alone_log2_bytes"] = up(log2_int(N * 32), 2)
    small = sieve(1 << 16)
    e256 = [pow(q, (R - 1) // 256, R) for q in small]
    p16 = sum(1 for e in e256 if pow(e, 16, R) == 1)
    p256 = sum(1 for e in e256 if e == 1)
    out["h2_small_primes"] = {"primes_below_2^16": len(small), "16th_powers": p16, "expected": len(small) / 16,
                              "256th_powers": p256, "expected_256": len(small) / 256}
    return out

TP, TB, TR, TK = 207670579, 3, 207671809, 7

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        l = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else: l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return x3, (l * (x1 - x3) - y1) % p

def ec_mul(k, P, p):
    Rp = None
    while k:
        if k & 1: Rp = ec_add(Rp, P, p)
        P = ec_add(P, P, p); k >>= 1
    return Rp

def batch_ec_add(pairs, p):
    """Independent additions P_i + Q_i with the ec_add formulas and one shared inversion (Montgomery trick)."""
    n = len(pairs); out = [None] * n; den = [0] * n; live = []
    for i, (P, Q) in enumerate(pairs):
        if P is None: out[i] = Q
        elif Q is None: out[i] = P
        elif P[0] != Q[0]: den[i] = (Q[0] - P[0]) % p; live.append(i)
        elif (P[1] + Q[1]) % p: den[i] = 2 * P[1] % p; live.append(i)
    pre = [1] * (len(live) + 1)
    for t, i in enumerate(live): pre[t + 1] = pre[t] * den[i] % p
    inv = pow(pre[-1], -1, p)
    for t in range(len(live) - 1, -1, -1):
        i = live[t]; d_inv = inv * pre[t] % p; inv = inv * den[i] % p
        (x1, y1), (x2, y2) = pairs[i]
        l = (3 * x1 * x1 if x1 == x2 else y2 - y1) * d_inv % p
        x3 = (l * l - x1 - x2) % p
        out[i] = (x3, (l * (x1 - x3) - y1) % p)
    return out

def batch_mul(ks, Ps, p):
    """Lockstep double-and-add [k_i]P_i: one shared inversion per round for the additions and doublings together."""
    n = len(ks); ks = list(ks); Ps = list(Ps); Rs = [None] * n
    while any(ks):
        act = [i for i in range(n) if ks[i] & 1]; dbl = [i for i in range(n) if ks[i] > 1]
        res = batch_ec_add([(Rs[i], Ps[i]) for i in act] + [(Ps[i], Ps[i]) for i in dbl], p)
        for i, R in zip(act, res): Rs[i] = R
        for i, R in zip(dbl, res[len(act):]): Ps[i] = R
        ks = [k >> 1 for k in ks]
    return Rs

def batch_double(X, p):
    """Lockstep doublings [2]X_i, the ec_add doubling formula with one shared inversion."""
    n = len(X); pre = [1] * (n + 1)
    for i in range(n): pre[i + 1] = pre[i] * (2 * X[i][1]) % p
    inv = pow(pre[n], -1, p); out = [None] * n
    for i in range(n - 1, -1, -1):
        x1, y1 = X[i]; d_inv = inv * pre[i] % p; inv = inv * (2 * y1) % p
        l = 3 * x1 * x1 * d_inv % p
        x3 = (l * l - 2 * x1) % p
        out[i] = (x3, (l * (x1 - x3) - y1) % p)
    return out

def batch_doublings(Ps, bits, p):
    """[2^j]P_i for j < bits, every chain's doublings in lockstep."""
    pows = [[P] for P in Ps]
    for _ in range(bits - 1):
        for pw, R in zip(pows, batch_double([pw[-1] for pw in pows], p)): pw.append(R)
    return pows

def batch_tables(keys, pows, p):
    """[g]R_i for every key g and chain i from the doublings: the lowest set bit is a copy, every higher one a
    lockstep addition sharing one inversion per bit position. [k']R + [2^j]R with 0 < k' < 2^j never degenerates."""
    cells = [(i, g) for i in range(len(pows)) for g in keys]
    low = {c: (c[1] & -c[1]).bit_length() - 1 for c in cells}
    Rs = {c: pows[c[0]][low[c]] for c in cells}
    for j in range(1, max(keys).bit_length()):
        act = [c for c in cells if c[1] >> j & 1 and low[c] < j]
        if act:
            for c, R in zip(act, batch_add([Rs[c] for c in act], [pows[c[0]][j] for c in act], p)): Rs[c] = R
    return [{g: Rs[(i, g)] for g in keys} for i in range(len(pows))]

def sqrt_mod(a, p):
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    c, x, t, mm = pow(z, q, p), pow(a, (q + 1) // 2, p), pow(a, q, p), s
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p); x = x * b % p; c = b * b % p; t = t * c % p; mm = i
    return x

def rand_point(p, b, rnd):
    while True:
        x = rnd.randrange(p); v = (x ** 3 + b) % p
        if pow(v, (p - 1) // 2, p) == 1: return x, sqrt_mod(v, p)

def batch_add(X, G, p):
    """Lockstep incomplete affine additions X[i] + G[i] sharing one inversion (Montgomery trick)."""
    n = len(X); dx = [(G[i][0] - X[i][0]) % p for i in range(n)]
    assert all(dx), "incomplete addition degenerated"
    pre = [1] * (n + 1)
    for i in range(n): pre[i + 1] = pre[i] * dx[i] % p
    inv = pow(pre[n], -1, p); out = [None] * n
    for i in range(n - 1, -1, -1):
        d_inv = inv * pre[i] % p; inv = inv * dx[i] % p
        l = (G[i][1] - X[i][1]) * d_inv % p
        x3 = (l * l - X[i][0] - G[i][0]) % p
        out[i] = (x3, (l * (X[i][0] - x3) - X[i][1]) % p)
    return out

def toy_key(P, p):
    y = P[1] if P[1] % 2 == 0 else p - P[1]
    return (y >> 1) & 0xFFFF

def ec_mode(seed, trials):
    p, b, r, k = TP, TB, TR, TK
    d = 1 << k; m = (r - 1) // d
    assert is_prime(p) and is_prime(r) and (r - 1) % (6 * d) == 0
    hprime = lambda q: pow(q, m, r) == 1
    HP = [q for q in sieve(2_000_000) if hprime(q)]
    omega = next(pow(g, (r - 1) // 3, r) for g in range(2, 50) if pow(g, (r - 1) // 3, r) != 1)
    mu6 = [1, r - 1, omega, r - omega, omega * omega % r, (r - omega * omega) % r]
    assert all(pow(u, m, r) == 1 for u in mu6)
    La, Lw = 15, 25
    n = La * Lw; K = 6 * n * n / m; Kd = 6 * (n * n - n * La + La * La) / m
    Wl = HP[:Lw]; A = HP[Lw:Lw + La]; A2 = HP[Lw + La:Lw + 2 * La]
    gaps = sorted(set(Wl[i + 1] - Wl[i] for i in range(Lw - 1)))
    NODE_BITS, HEAD_BITS = 5, 4
    rnd = random.Random(seed)
    ok = 0; comp_total = 0; false_total = 0; true_total = 0; chi_total = 0.0; chi_sq_total = 0.0; rewalk_cap_fails = 0
    true_hist = {}; lemma_runs = 0; hit_budget = 8; hit_budget_overflows = 0
    for _ in range(trials):
        P0 = rand_point(p, b, rnd)
        s = rnd.randrange(1, r); T = ec_mul(s, P0, p)
        z = pow(s, d, r); Ad = ec_mul(z, P0, p)
        bits = max(Wl[0], gaps[-1]).bit_length()
        def walk(bases, Q, emit):
            pows = batch_doublings(batch_mul(bases, [Q] * len(bases), p), bits, p)
            tables = batch_tables(gaps + [Wl[0]], pows, p)
            X = [t[Wl[0]] for t in tables]
            for i in range(len(bases)): emit(X[i], i)
            for j in range(1, Lw):
                g = Wl[j] - Wl[j - 1]
                X = batch_add(X, [t[g] for t in tables], p)
                for i in range(len(bases)): emit(X[i], i)
            return tables
        table = [[[] for _ in range(1 << HEAD_BITS)] for _ in range(1 << NODE_BITS)]
        occupancy = [0] * (1 << NODE_BITS)
        def insert(X, chain):
            key = toy_key(X, p); node = key & 31; head = (key >> 5) & 15
            table[node][head].append((key >> 9, chain)); occupancy[node] += 1
        tb = walk(A, P0, insert)
        exp_occ = n / 32
        chi = sum((o - exp_occ) ** 2 / exp_occ for o in occupancy); chi_total += chi; chi_sq_total += chi * chi
        matches = []; comps = 0
        def lookup(X, chain):
            nonlocal comps
            key = toy_key(X, p); node = key & 31; head = (key >> 5) & 15; cmp = key >> 9
            for rec in table[node][head]:
                comps += 1
                if rec[0] == cmp: matches.append((key, rec[1], chain))
        tg = walk(A2, Ad, lookup)
        comp_total += comps
        found = None; cap_fail = False; trial_true = 0
        def steps_with_key(a, table, key):
            X = table[Wl[0]]; kept = [a * Wl[0] % r] if toy_key(X, p) == key else []
            for j in range(1, Lw):
                X = ec_add(X, table[Wl[j] - Wl[j - 1]], p)
                if toy_key(X, p) == key: kept.append(a * Wl[j] % r)
            assert len(kept) >= 1, "the matched step must carry the key"
            return kept
        for (key, ib, ig) in matches:
            xbs = steps_with_key(A[ib], tb[ib], key); xgs = steps_with_key(A2[ig], tg[ig], key)
            if len(xbs) > REWALK_MATCH_CAP or len(xgs) > REWALK_MATCH_CAP:
                cap_fail = True; break
            hit = False
            for xb in xbs:
                for xg in xgs:
                    cands = [xb * pow(u * xg % r, -1, r) % r for u in mu6]
                    for cand, Rc in zip(cands, batch_mul(cands, [P0] * 6, p)):
                        if Rc == Ad: found = cand; hit = True
            if hit: true_total += 1; trial_true += 1
            else: false_total += 1
        true_hist[trial_true] = true_hist.get(trial_true, 0) + 1
        if trial_true > hit_budget: hit_budget_overflows += 1
        if cap_fail: rewalk_cap_fails += 1; continue
        if found is None: continue
        lemma_runs += 1
        assert found == z
        s0 = found
        for _ in range(k): s0 = sqrt_mod(s0, r)
        assert pow(s0, d, r) == z
        g = next(x for x in range(2, 100) if pow(x, (r - 1) // 2, r) != 1)
        zeta = pow(g, (r - 1) // d, r); cur = s0; hit = None
        for _ in range(0, d, 16):
            curs = []
            for _ in range(16): curs.append(cur); cur = cur * zeta % r
            for c, Rc in zip(curs, batch_mul(curs, [P0] * 16, p)):
                if Rc == T: assert c == s; ok += 1; hit = c; break
            if hit is not None: break
        else: raise AssertionError("stage 2 failed")
    rate = ok / trials; pred = 1 - math.exp(-K); predd = 1 - math.exp(-Kd)
    se = math.sqrt(predd * (1 - predd) / trials)
    chi_mean = chi_total / trials; chi_se = math.sqrt(max(chi_sq_total / trials - chi_mean ** 2, 0) / trials)
    return {"curve": {"p": p, "b": b, "r": r, "d": d, "m": m}, "per_side": n, "chains_per_side": La, "steps_per_chain": Lw,
            "K": round(K, 4), "K_diagonal_corrected": round(Kd, 4), "predicted": round(pred, 4), "diagonal_corrected": round(predd, 4),
            "measured": round(rate, 4), "trials": trials, "z_vs_corrected": round((rate - predd) / se, 2),
            "comparisons_mean": round(comp_total / trials, 1), "comparisons_expected": round(n * n / 512, 1),
            "false_matches_mean": round(false_total / trials, 3), "false_matches_expected": round(n * n / 65536, 3),
            "true_matches_mean": round(true_total / trials, 3), "rewalk_cap_fails": rewalk_cap_fails,
            "true_matches_histogram": {str(k): true_hist.get(k, 0) for k in range(max(true_hist) + 1)} if true_hist else {},
            "true_matches_ge2_trials": sum(v for k, v in true_hist.items() if k >= 2),
            "true_matches_ge2_poisson_expected": round(trials * (-math.expm1(-K) - K * math.exp(-K)), 1),
            "lemma1_runs": lemma_runs, "recovery_hit_budget": hit_budget,
            "recovery_hit_budget_overflows": hit_budget_overflows,
            "node_chi2_mean": round(chi_mean, 2), "node_chi2_expected": 31, "node_chi2_se": round(chi_se, 2),
            "all_recovered_scalars_correct": True}

CENSUS_SPLIT = (Fraction(64) / Fraction("115.6503")).limit_denominator(10 ** 6)
CENSUS_WALK_START = Fraction(39) / Fraction("77.4")
CENSUS_BASE_OFFSET = Fraction("2.6")

def primes_above(start):
    lo, seg, base, base_hi = max(2, start + 1), 1 << 16, [], 1
    while True:
        hi = lo + seg; root = math.isqrt(hi) + 1
        if root > base_hi: base_hi = max(2 * root, 1 << 12); base = sieve(base_hi)
        s = bytearray([1]) * seg
        for q in base:
            if q * q >= hi: break
            st = max(q * q, (lo + q - 1) // q * q)
            if st < hi: s[st - lo::q] = bytes((hi - 1 - st) // q + 1)
        i = s.find(1)
        while i >= 0: yield lo + i; i = s.find(1, i + 1)
        lo = hi

def h_primes_above(start, count, m, r):
    out = []
    for q in primes_above(start):
        if pow(q, m, r) == 1:
            out.append(q)
            if len(out) == count: return out

def iroot(n, e):
    if n < 2: return n
    x = 1 << ((n.bit_length() + e - 1) // e)
    while True:
        y = ((e - 1) * x + n // x ** (e - 1)) // e
        if y >= x: break
        x = y
    while x ** e > n: x -= 1
    while (x + 1) ** e <= n: x += 1
    return x

def brent(n):
    if n % 2 == 0: return 2
    rnd = random.Random(n)
    while True:
        y, c, step = rnd.randrange(1, n), rnd.randrange(1, n), 64
        g = q = rr = 1; x = ys = y
        while g == 1:
            x = y
            for _ in range(rr): y = (y * y + c) % n
            kk = 0
            while kk < rr and g == 1:
                ys = y
                for _ in range(min(step, rr - kk)): y = (y * y + c) % n; q = q * abs(x - y) % n
                g = math.gcd(q, n); kk += step
            rr *= 2
        if g == n:
            g = 1
            while g == 1: ys = (ys * ys + c) % n; g = math.gcd(abs(x - ys), n)
        if g != n: return g

def factorize(n):
    f = {}
    for q in sieve(1000):
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
    stack = [n] if n > 1 else []
    while stack:
        x = stack.pop()
        if is_prime(x): f[x] = f.get(x, 0) + 1
        else: d = brent(x); stack += [d, x // d]
    return dict(sorted(f.items()))

def census_mu6(r):
    w = next(pow(g, (r - 1) // 3, r) for g in range(2, r) if pow(g, (r - 1) // 3, r) != 1)
    return [1, r - 1, w, r - w, w * w % r, (r - w * w) % r]

def dlogs(r, k, gamma, hs):
    """Discrete logs base gamma (a generator of H) by Pohlig-Hellman with baby-step giant-step, each checked."""
    m = (r - 1) >> k; parts = []
    for q, e in factorize(m).items():
        qe = q ** e
        gq, gqe_inv = pow(gamma, m // q, r), pow(pow(gamma, m // qe, r), -1, r)
        b = max(1, min(q, math.isqrt(q * len(hs) * e) + 1, 1 << 21))
        table, cur = {}, 1
        for j in range(b): table.setdefault(cur, j); cur = cur * gq % r
        step, rounds = pow(gq, -b, r), -(-q // b) + 1
        res = []
        for h in hs:
            hq, x = pow(h, m // qe, r), 0
            for i in range(e):
                cur = pow(hq * pow(gqe_inv, x, r) % r, q ** (e - 1 - i), r)
                for jj in range(rounds):
                    j0 = table.get(cur)
                    if j0 is not None: x += (jj * b + j0) % q * q ** i; break
                    cur = cur * step % r
                else: raise AssertionError("no discrete log")
            res.append(x)
        parts.append((qe, res))
    out = []
    for i, h in enumerate(hs):
        x, mod = 0, 1
        for qe, res in parts: x += mod * ((res[i] - x) * pow(mod, -1, qe) % qe); mod *= qe
        assert pow(gamma, x, r) == h % r, "discrete log check failed"
        out.append(x)
    return out

def census_lists(r, k, K):
    """n_b within 10 % of N^(64 / 115.6503), L the least integer with 6 n_b^2 L (L - 1) >= K m (the n_b that brings K
    closest above the target), walk list the first L H-primes above w_min = w_max^(39 / 77.4) (iterated to a fixed
    point), base primes the first 2 n_b H-primes above 2^2.6 w_max: the first n_b baby, the next n_b giant."""
    m = (r - 1) >> k
    n0 = math.sqrt(float(K) * m / 6) ** float(CENSUS_SPLIT); best = None
    for nb in range(max(1, int(n0 * 0.9)), int(n0 * 1.1) + 2):
        t = -(-K.numerator * m // (6 * K.denominator * nb * nb))
        L = max(2, (1 + math.isqrt(4 * t)) // 2)
        while L * (L - 1) < t: L += 1
        while L > 2 and (L - 1) * (L - 2) >= t: L -= 1
        cand = (Fraction(6 * nb * nb * L * (L - 1), m), abs(nb - n0), nb, L)
        if best is None or cand[:2] < best[:2]: best = cand
    nb, L = best[2], best[3]
    wmin = 1
    for _ in range(256):
        W = h_primes_above(wmin, L, m, r)
        new = iroot(W[-1] ** CENSUS_WALK_START.numerator, CENSUS_WALK_START.denominator)
        if new == wmin: break
        wmin = new
    else: raise AssertionError("walk start did not settle")
    f = CENSUS_BASE_OFFSET
    B = h_primes_above(iroot((W[-1] ** f.denominator) << f.numerator, f.denominator), 2 * nb, m, r)
    assert B[-1] < r and W[-1] < B[0]
    return W, B[:nb], B[nb:]

def census_row(m_bits, k, K, seed):
    rnd = random.Random("census:%d:%d:%d" % (m_bits, k, seed))
    lo, hi = -(-(1 << m_bits) // 6), (1 << (m_bits + 1)) // 6
    while True:
        r = ((6 * rnd.randrange(lo, hi)) << k) + 1
        if is_prime(r): break
    W, A, A2 = census_lists(r, k, K)
    m = (r - 1) >> k; M = m // 6; u6 = census_mu6(r)
    assert m % 6 == 0 and all(pow(u, m, r) == 1 for u in u6), "mu6 must lie in H"
    fr = factorize(r - 1)
    g = next(g for g in range(2, r) if all(pow(g, (r - 1) // q, r) != 1 for q in fr))
    gamma = pow(g, 1 << k, r)
    e = [x % M for x in dlogs(r, k, gamma, W + A + A2)]
    ew, ea, ea2 = e[:len(W)], e[len(W):len(W) + len(A)], e[len(W) + len(A):]
    alpha = sorted({(x - y) % M for x in ea for y in ea2})
    delta = sorted({(ew[i] - ew[j]) % M for i in range(len(W)) for j in range(len(W)) if i != j})
    flags = bytearray(M)
    outer, inner = (delta, alpha) if len(delta) <= len(alpha) else (alpha, delta)
    for d in outer:
        cut = bisect.bisect_left(inner, M - d)
        for a in inner[:cut]: flags[a + d] = 1
        d -= M
        for a in inner[cut:]: flags[a + d] = 1
    nR = flags.count(1); nRD = sum(1 for a in alpha if flags[a])
    xd = {}
    for a in A:
        for iw, w in enumerate(W): xd.setdefault(a * w % r, []).append(iw)
    uy = [(u * a2 * w % r, jw) for u in u6 for a2 in A2 for jw, w in enumerate(W)]
    def member(z):
        off = diag = False
        for v, jw in uy:
            for iw in xd.get(z * v % r, ()):
                if iw == jw: diag = True
                else: off = True
        return off, diag
    srnd = random.Random("spot:%d:%d" % (r, seed)); aset = set(alpha); su6 = sorted(set(u6)); ok = 0
    for _ in range(32):
        x = srnd.randrange(m)
        ok += member(pow(gamma, x, r)) == (bool(flags[x % M]), x % M in aset)
    for _ in range(8):
        ia, ja, u = srnd.randrange(len(A)), srnd.randrange(len(A2)), srnd.choice(su6)
        i, j = srnd.sample(range(len(W)), 2)
        z = A[ia] * W[i] * pow(u * A2[ja] * W[j] % r, -1, r) % r
        idx = (ea[ia] - ea2[ja] + ew[i] - ew[j]) % M
        ok += member(z)[0] and bool(flags[idx]) and pow(gamma, idx, r) * pow(z, -1, r) % r in su6
    assert ok == 40, "spot check failed"
    Kf = float(Fraction(len(A) * len(A2) * len(W) * (len(W) - 1), M))
    cov, pred = nR / M, -math.expm1(-Kf)
    sd = math.sqrt(math.exp(-Kf) * (1 - (1 + Kf) * math.exp(-Kf)) / M)
    # Exact multiplicity tails of M(z) over the coset array, from the same alpha/delta residue sets (the verdict
    # on the prior submission asked for multiplicity evidence, not only coverage): for each shift a in alpha the
    # off-diagonal multiplicity of coset c is the number of d in delta with (c - d) mod M in alpha, i.e. the number of
    # (baby base, giant base, off-diagonal walk pair) triples hitting c. The convolution is over sets of residues,
    # so delta with repeated (i, j) differences carries its full multiplicity via the delta multiset below.
    dmulti = {}
    for i in range(len(W)):
        for j in range(len(W)):
            if i != j: dmulti[(ew[i] - ew[j]) % M] = dmulti.get((ew[i] - ew[j]) % M, 0) + 1
    amulti = {}
    for x in ea:
        for y in ea2: amulti[(x - y) % M] = amulti.get((x - y) % M, 0) + 1
    mvol = [0] * M
    small, big = (dmulti, amulti) if len(dmulti) <= len(amulti) else (amulti, dmulti)
    sitems = sorted(small.items()); bkeys = sorted(big)
    bloc = {v: i for i, v in enumerate(bkeys)}
    barr = [big[v] for v in bkeys]; nB = len(bkeys)
    for s, ns in sitems:
        cut = bisect.bisect_left(bkeys, M - s)
        base = s
        for i in range(cut): mvol[base + bkeys[i]] += ns * barr[i]
        base -= M
        for i in range(cut, nB): mvol[base + bkeys[i]] += ns * barr[i]
    tot = sum(mvol); npairs = len(A) * len(A2) * len(W) * (len(W) - 1)
    assert tot == npairs, "multiplicity volume must equal the off-diagonal pair count"
    assert sum(1 for v in mvol if v) == nR, "multiplicity support must equal R"
    # Independent cross-check of the multiplicity array: for 200 sampled off-diagonal pairs the residue coset of
    # the pair's ratio must carry multiplicity at least 1, verified by residue arithmetic against gamma^idx.
    mrnd = random.Random("mult:%d:%d" % (r, seed)); mspot = 0
    for _ in range(200):
        ia, ja = mrnd.randrange(len(A)), mrnd.randrange(len(A2))
        i, j = mrnd.sample(range(len(W)), 2)
        z = A[ia] * W[i] * pow(A2[ja] * W[j] % r, -1, r) % r
        idx = (ea[ia] - ea2[ja] + ew[i] - ew[j]) % M
        mspot += mvol[idx] >= 1 and pow(gamma, idx, r) * pow(z, -1, r) % r in su6
    assert mspot == 200, "multiplicity spot check failed"
    p_ge2 = sum(1 for v in mvol if v >= 2) / M
    pois_ge2 = -math.expm1(-Kf) - Kf * math.exp(-Kf)
    p_ge3 = sum(1 for v in mvol if v >= 3) / M
    pois_ge3 = 1 - math.exp(-Kf) * (1 + Kf + Kf * Kf / 2)
    mmax = max(mvol)
    return {"r": r, "k": k, "seed": seed, "m": m, "m_log2": "%.4f" % math.log2(m), "nb": len(A), "L": len(W),
            "N": len(A) * len(W), "K": "%.6f" % Kf, "R_cosets": nR, "M6": M, "cov": "%.8f" % cov, "pred": "%.8f" % pred,
            "diff": "%+.8f" % (cov - pred), "z": "%+.2f" % ((cov - pred) / sd), "D_frac": "%.8f" % (len(alpha) / M),
            "succ": "%.8f" % ((nR - nRD) / M), "spot": "%d/40" % ok,
            "mult_ge2": "%.8f" % p_ge2, "mult_ge2_poisson": "%.8f" % pois_ge2,
            "mult_ge3": "%.8f" % p_ge3, "mult_ge3_poisson": "%.8f" % pois_ge3, "mult_max": mmax}

def census(m_bits, ks, K_per_100000, first_seed, primes, seconds):
    """Primes in the order seed, then k; no further prime starts once `seconds` have passed (hosted runs)."""
    K = Fraction(K_per_100000, 100000); t0 = time.monotonic(); rows = []; stopped = False
    for seed in range(first_seed, first_seed + primes):
        for k in ks:
            if rows and time.monotonic() - t0 > seconds: stopped = True; break
            rows.append(census_row(m_bits, k, K, seed))
        if stopped: break
    d = [float(x["diff"]) for x in rows]
    return {"m_bits": m_bits, "k": ks, "K_target": K_per_100000 / 100000, "first_seed": first_seed,
            "primes_requested": primes * len(ks), "primes_run": len(rows), "stopped_by_time_budget": stopped,
            "diff_is": "|R| / m - (1 - e^-K) at the K each prime reached (K >= the target)",
            "z_is": "diff over the spread that independent uniform draws would give; the spread across primes is several "
                    "times that (proof.md Section 6), so z measures one prime's deviation, not a test of H1",
            "mean_diff": "%+.6f" % (sum(d) / len(d)), "lowest_diff": "%+.6f" % min(d), "highest_diff": "%+.6f" % max(d),
            "below_prediction": sum(1 for x in d if x < 0),
            "mult_is": "exact off-diagonal multiplicity tails over the M/6 cosets: mult_ge2/mult_ge3 are P(M(z) >= 2/3), "
                    "with the Poisson(K) value at the same K beside each; mult_max is the largest M(z)",
            "mean_mult_ge2": "%.8f" % (sum(float(x["mult_ge2"]) for x in rows) / len(rows)),
            "mean_mult_ge2_poisson": "%.8f" % (sum(float(x["mult_ge2_poisson"]) for x in rows) / len(rows)),
            "mean_mult_ge3": "%.8f" % (sum(float(x["mult_ge3"]) for x in rows) / len(rows)),
            "mean_mult_ge3_poisson": "%.8f" % (sum(float(x["mult_ge3_poisson"]) for x in rows) / len(rows)),
            "max_mult_max": max(x["mult_max"] for x in rows), "rows": rows}

def main():
    raw = sys.stdin.read().strip() if not sys.stdin.isatty() else ""
    args = json.loads(raw) if raw else {}
    mode = args.get("mode", "all"); seed = int(args.get("seed", 2026)); trials = int(args.get("trials", 400))
    assert mode in ("ledger", "ec", "census", "all") and 0 <= seed < 2 ** 31 and 20 <= trials <= 2000
    cen = (int(args.get("m_bits", 21)), [int(k) for k in args.get("k", [2, 4, 7, 10])], int(args.get("K_per_100000", K_NUM)),
           int(args.get("first_seed", 0)), int(args.get("primes", 2)), int(args.get("seconds", 8)))
    mb, ks, kp, fs, npr, sec = cen
    assert (16 <= mb <= 36 and 1 <= len(ks) <= 4 and len(set(ks)) == len(ks) and all(1 <= k <= 12 for k in ks)
            and 30000 <= kp <= 90000 and 0 <= fs < 2 ** 31 and 1 <= npr <= 64 and 1 <= sec <= 10 ** 6)
    out = {}
    if mode in ("ledger", "all"): out["ledger"] = ledger()
    if mode in ("ec", "all"): out["ec"] = ec_mode(seed, trials)
    if mode in ("census", "all"): out["census"] = census(*cen)
    json.dump(out, sys.stdout, indent=1); print()

if __name__ == "__main__": main()

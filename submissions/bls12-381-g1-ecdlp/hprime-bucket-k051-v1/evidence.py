#!/usr/bin/env python3
# Evidence for hprime-bucket-k051-v1 (bls12-381-g1-ecdlp). Standard library only, no imports from the package.
# Input: one JSON object on stdin, {"mode": "ledger" | "ec" | "census" | "all", "seed": int, "trials": int,
# "m_bits": int, "k": [int], "K_percent": int, "first_seed": int, "primes": int, "seconds": int}; a missing key takes
# its default ("all", 2026, 400, 21, [2, 4, 7, 10], 51, 0, 2, 8). Output: one JSON object on stdout.
#   ledger  exact integer recomputation of every number in proof.md (sizes, caps, memory, score)
#   ec      Algorithm A end to end on a genuine j=0 curve of prime order: lockstep affine walks with batch
#           inversion, truncated parity key, bucket routing to pseudo-nodes, chained head table, counted
#           comparisons and false matches, re-walk on match, six-candidate check, d-th root and root-of-unity scan
#   census  exact count of the ratio set R on toy primes built by the production rules (H1, proof.md Section 6)
import sys, json, math, random, bisect, time
from fractions import Fraction

# ---------------------------------------------------------------- production parameters (proof.md Section 2)
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001   # BLS12-381 scalar field order
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB  # base field order
X_BLS = -0xD201000000010000     # BLS parameter: r = x^4 - x^2 + 1, p = (x - 1)^2 r / 3 + x
FP_MUL_INSTR = 383              # one blst_fp_mul as measured by the cost model (S = M); also charged per F_r multiplication
LOG_D = 20                      # d = 2^20 divides r-1
D = 1 << LOG_D
M = (R - 1) // D                # |H|, H = d-th powers in F_r^*
C = 2000                        # instructions per reference unit (batched incomplete affine G1 addition)
K_NUM, K_DEN = 51, 100          # coverage parameter K = 6 N^2 / m >= 0.51
LOG_CHAINS = 64                 # chains per side; each chain is one base H-prime a walking the shared list W
LOG_NODES = 84                  # owner nodes, 3 routing hops of 28 bits
HOP_BITS, HOPS = 28, 3
LOG_HEADS = 33                  # direct-mapped head array per owner node (8-byte slots)
LOG_NODE_BYTES = 40             # every node is counted at its full 2^40 bytes
LOG_POOL_RECORDS = 33           # chained pool capacity per owner node (32-byte records)
LOG_GAP_TABLE = 32              # even gaps g <= 2^33, entry [a g]Q, 96 bytes each
LOG_GEN_NODES, LOG_STORE_NODES = 20, 30
KEY_BITS = 192
OVERHEAD_LISTED_BABY = {"gap stream": 5, "table address": 2, "operand loads and stores": 18, "lane control": 8,
                        "parity key": 23, "scatter 3 hops": 36, "post and poll 3 hops": 9, "transport copy 3 hops": 18,
                        "gap stream transport from storage": 3, "batch inversion loads and stores": 6,
                        "flow control and counters": 5, "shared F_p inversion of the round": 4, "insert": 10}
OVERHEAD_LISTED_GIANT = dict(OVERHEAD_LISTED_BABY); del OVERHEAD_LISTED_GIANT["insert"]; OVERHEAD_LISTED_GIANT["lookup fixed"] = 6
OVERHEAD_CHARGED = 180          # instructions charged per record on top of the reference addition
# per-step overhead of a walk that produces no record (gap-table build, re-walk): the stream, address, operand,
# inversion, lane and flow items of the listing; the table build stores a 96 B entry (8) instead of reading the gap
# stream (5 + 3), so the same 47 covers both
WALK_STEP_OVERHEAD = sum(OVERHEAD_LISTED_BABY[k] for k in ("gap stream", "gap stream transport from storage", "table address",
                         "operand loads and stores", "batch inversion loads and stores", "lane control", "flow control and counters"))
assert WALK_STEP_OVERHEAD == 47
INSTR_PER_COMPARISON = 8        # load 24 B, compare, branch, pointer advance, counter decrement
SLOT_CHAIN_CAP = 64             # baby records per head slot; more -> fail
BERNSTEIN_VAR_MARGIN = 4
LOG_NODE_FIXED_INSTR = 35       # head zeroing 2^33, bucket pointer setup and end-of-hop flushes 3 x 2^33, per node
LOG_SCALAR_MULT_UNITS = 12      # one 255-bit scalar multiplication, Jacobian double-and-add, one final inversion (about 2^9.5)
LOG_HTEST_INSTR = 18            # one exponentiation q^m mod r at 383 instructions per F_r multiplication (about 2^17.1)
COMPARISON_CAP_NUM, COMPARISON_CAP_DEN = 101, 100
LOG_MATCH_CAP = 44
CANDIDATE_BYTES = 40            # a candidate pair keeps the full 192-bit key and both chain ids
LOG_NODE_CANDIDATE_CAP = 16     # candidate list capacity per owner node; more -> fail
CANDIDATE_FORWARD_INSTR = 128   # copy one 40-byte candidate into two 512-byte blocks, post and poll, at both generating nodes
LOG_GEN_CANDIDATE_CAP = 24      # candidate receive area per generating node, 40-byte slots; more -> fail
REWALK_STEP_INSTR = 32          # per re-walk step: parity key 23, compare with the retained key 8, position bookkeeping 1
REWALK_MATCH_CAP = 4            # steps of one chain allowed to carry the retained key; more -> fail
KEPT_FORWARD_INSTR = 256        # one 168-byte block (40 B header + 4 x 32 B kept values) from g(id') to g(id): copy, post, poll
LOG_STATUS_INSTR = 12           # per node per pass of the fan-in-256 status tree (48 B records, 3 passes)
PRE_CLAIM_LOG2_E4 = 1157815    # quoted preprocessing bound 2^115.7815 units, checked exactly below
MEMORY_CLAIM_LOG2_E4 = 1233242  # claimed memory bound 2^123.3242 bytes, checked exactly below (exact value 2^123.324186...)
LOG_STORAGE_NODE_BYTES = 40     # walk list and base primes (2^69.0 B) spread over 2^30 storage nodes at their full 2^40 B
LOG_BIGGAP_RESERVE = 69         # 16 split gaps x 2^65 chains x one extra addition
BIGGAP_SPLITS = 16              # split gaps allowed in the walk list; more, or one above 2^34, stops the run
LOG_ROUND = 16                  # a chain walks as 2^16 lockstep segments and a gap table is built by 2^16 strands
LOG_WALK_SIEVE, LOG_WALK_HTESTS = 78, 73   # walk list: sieve (2^39, 2^78], H-test every prime there
LOG_BASE_SIEVE, LOG_BASE_HTESTS = 92, 86   # base primes: sieve (2^80, 2^80 + 2^92], H-test at most its first 2^86 primes
LOG_STAGE2 = 24

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

# ---------------------------------------------------------------- ledger
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
    # sizes: smallest N with 6 N^2 / m >= K, rounded up to a whole number of chains
    n0 = math.isqrt(K_NUM * M // (6 * K_DEN)) + 1
    chains = 1 << LOG_CHAINS
    L = ceil_div(n0, chains)                      # steps per chain
    N = L * chains                                # records per side
    assert 6 * N * N * K_DEN >= K_NUM * M
    out["steps_per_chain_log2"] = up(log2_int(L), 4)
    out["records_per_side_log2"] = up(log2_int(N), 4)
    K = 6 * N * N / M
    out["K"] = down(K, 6)
    k_excess = Fraction(6 * N * N, M) - Fraction(K_NUM, K_DEN)          # exact K - 0.51 >= 0
    assert k_excess >= 0
    k_excess_log2 = up(log2_int(k_excess.numerator) - log2_int(k_excess.denominator), 1)
    out["K_minus_0.51_log2_below"] = k_excess_log2
    out["predicted_success"] = math.floor((1 - math.exp(-K)) * 10 ** 6) / 10 ** 6   # rounded down
    out["diagonal_loss_log2"] = up(-log2_int(L), 2)
    # walk list and base primes (H2 sizes)
    hp_needed = L
    primes_needed = hp_needed * D                 # one prime in 2^20 is an H-prime
    x = primes_needed * 54                        # pi(x) ~ x / ln x, ln 2^78 = 54.1
    out["walk_list_last_prime_log2_about"] = round(log2_int(x), 2)
    assert x < 1 << 79
    base_count = 2 * chains
    out["base_primes_needed_log2"] = LOG_CHAINS + 1
    # H2 windows: the primes each capped search examines. Rosser-Schoenfeld: pi(x) > x / ln x for x >= 17 and
    # pi(x) < 1.25506 x / ln x for x > 1, so a window (a, b] holds more than lo(b) - hi(a) and fewer than hi(b) primes
    pi_lo = lambda v: v / math.log(v)
    pi_hi = lambda v: 1.25506 * v / math.log(v)
    walk_lo, walk_hi = pi_lo(2.0 ** LOG_WALK_SIEVE) - pi_hi(2.0 ** 39), pi_hi(2.0 ** LOG_WALK_SIEVE)
    assert walk_hi < 2 ** LOG_WALK_HTESTS          # every prime of (2^39, 2^78] is H-tested within the charge
    assert walk_lo / D > L                         # expected H-primes in the window exceed the L needed
    out["h2_walk_window"] = {"sieved": "(2^39, 2^78]", "primes_log2_above": math.floor(math.log2(walk_lo) * 100) / 100,
                             "primes_log2_below": up(math.log2(walk_hi), 2), "h_tests_charged_log2": LOG_WALK_HTESTS,
                             "examined": "its primes in increasing order up to the L-th H-prime, all within the 2^73 H-tests",
                             "h_primes_expected_log2_above": math.floor(math.log2(walk_lo / D) * 100) / 100,
                             "needed_log2": up(log2_int(L), 4),
                             "expected_over_needed_above": math.floor(walk_lo / D / L * 100) / 100}
    base_lo = pi_lo(2.0 ** 80 + 2.0 ** LOG_BASE_SIEVE) - pi_hi(2.0 ** 80)
    assert base_lo > 2 ** LOG_BASE_HTESTS          # the 2^86 H-test cap ends the search inside the sieved window
    base_expected = 1 << (LOG_BASE_HTESTS - LOG_D) # H-primes expected among the first 2^86 primes above 2^80
    assert base_expected >= 2 * base_count
    out["h2_base_window"] = {"sieved": "(2^80, 2^80 + 2^92]", "primes_log2_above": math.floor(math.log2(base_lo) * 1000) / 1000,
                             "h_tests_charged_log2": LOG_BASE_HTESTS, "examined": "at most the first 2^86 primes above 2^80, in increasing order up to the 2^65-th H-prime",
                             "h_primes_expected_log2": LOG_BASE_HTESTS - LOG_D, "needed_log2": LOG_CHAINS + 1,
                             "expected_over_needed": base_expected // base_count}
    assert (1 << 80) + (1 << LOG_BASE_SIEVE) < 1 << 93 and (1 << (93 + LOG_WALK_SIEVE)) < R
    out["x_max_log2_below"] = 93 + LOG_WALK_SIEVE  # a < 2^93 by the base window, w <= 2^78 by the walk window
    mean_gap = D * LOG_WALK_SIEVE * math.log(2)    # 2^20 ln q at the top of the walk window
    out["mean_gap_log2_at_most"] = up(math.log2(mean_gap), 2)
    out["gap_above_2^33_expected_log2_below"] = up(math.log2(L - 1) - (1 << 33) / mean_gap * math.log2(math.e), 1)
    # per-record charge
    lb, lg = sum(OVERHEAD_LISTED_BABY.values()), sum(OVERHEAD_LISTED_GIANT.values())
    out["overhead_listed"] = {"baby": lb, "giant": lg, "charged": OVERHEAD_CHARGED}
    assert OVERHEAD_CHARGED >= max(lb, lg)
    c_num = C + OVERHEAD_CHARGED
    # field inversions: Fermat powers x^(p-2) and x^(r-2) by binary exponentiation, every squaring and
    # multiplication charged at FP_MUL_INSTR
    assert R == X_BLS ** 4 - X_BLS ** 2 + 1 and P == (X_BLS - 1) ** 2 * R // 3 + X_BLS
    inv_p_mults = (P - 2).bit_length() - 1 + bin(P - 2).count("1") - 1
    inv_r_mults = (R - 2).bit_length() - 1 + bin(R - 2).count("1") - 1
    inv_p = inv_p_mults * FP_MUL_INSTR
    walk_rounds = ceil_div(L, 1 << LOG_ROUND) + BIGGAP_SPLITS   # rounds of one chain walk; a split delays one segment one round
    table_rounds = 1 << (LOG_GAP_TABLE - LOG_ROUND)               # 2^16 strands of at most 2^16 additions each
    listed_inv = OVERHEAD_LISTED_BABY["shared F_p inversion of the round"]
    assert walk_rounds * inv_p <= L * listed_inv                   # the listed 4 per record pays every record round's inversion
    # six z per kept pair: one F_r inversion per kept giant value, the 2 x 4 products a w and a' w', and per pair one
    # multiplication for a w / (a' w') and six for the z
    z_mults = REWALK_MATCH_CAP * inv_r_mults + 2 * REWALK_MATCH_CAP + REWALK_MATCH_CAP ** 2 * 7
    out["field_inversions"] = {"F_p_squarings": (P - 2).bit_length() - 1, "F_p_multiplications": bin(P - 2).count("1") - 1,
                               "F_p_inversion_instr": inv_p, "rounds_per_chain_walk": "ceil(L / 2^16) + 16",
                               "rounds_per_gap_table": table_rounds,
                               "per_record_share_instr_below": up(walk_rounds * inv_p / L, 2), "per_record_listed": listed_inv,
                               "F_r_inversion_multiplications": inv_r_mults, "z_forming_multiplications_per_candidate": z_mults}
    # node occupancy and caps
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
    comp_cap = node_comp_cap * nodes                      # every node stops at its own cap, so the total is bounded
    out["comparisons_mean_log2"] = up(log2_int(comp_mean), 3)
    out["comparisons_per_node_mean_log2"] = up(log2_int(node_comp_mean), 3)
    out["comparisons_per_node_cap_log2"] = up(log2_int(node_comp_cap), 3)
    # per-slot chain cap: Poisson(lam) tail at SLOT_CHAIN_CAP + 1 (the fail event), times all head slots
    k = SLOT_CHAIN_CAP
    k1 = k + 1
    slot_tail_log2 = (-lam * math.log2(math.e) + k1 * math.log2(lam) - math.log2(math.factorial(k1))
                      - math.log2(1 - lam / (k1 + 1)))
    out["slot_chain_cap_exceeded_all_slots_log2"] = up(slot_tail_log2 + math.log2(heads), 1)
    # per-node comparisons: sum over the node's giant records of a slot count bounded by SLOT_CHAIN_CAP;
    # Bernstein with M = cap, variance (lam + lam^2) per term times a margin; 0.2 % of the 1 % cap margin is
    # reserved for the baby and giant counts of the node each deviating up to 2^-10 (0.098 %) from their means
    t = 0.008 * node_comp_mean
    V = mean_per_node * (lam + lam * lam) * BERNSTEIN_VAR_MARGIN
    bern_log2 = -(t * t / (2 * (V + k * t / 3))) * math.log2(math.e)
    out["comparisons_node_cap_exceeded_bernstein_all_nodes_log2"] = up(bern_log2 + LOG_NODES)
    count_dev_log2 = -(mean_per_node * (2 ** -10) ** 2 / 3) * math.log2(math.e) + LOG_NODES + 1
    out["node_count_off_by_2^-10_chernoff_all_nodes_log2"] = up(count_dev_log2)
    false_mean = N * N // (1 << KEY_BITS)
    out["false_matches_mean_log2"] = up(log2_int(false_mean), 3)
    match_cap = 1 << LOG_MATCH_CAP
    out["match_cap_exceeded_chebyshev_log2"] = up(math.log2(2 * false_mean + 1) - log2_int((match_cap - false_mean) ** 2), 1)
    # true matches (pairs with a w = u a' w' z); the ledger budgets 2^20 of them. Off the diagonal (w != w') their
    # number is Poisson with mean at most K under H1. On the diagonal, z = a / (u a') for one baby base a, one giant
    # base a' and one u in mu6 makes all L positions w = w' match, so the global candidate cap fails the run; z is
    # uniform in H and these z form a set of at most 6 x 2^64 x 2^64 elements, so the event has probability at most
    # 6 n_b^2 / m, a separate failure term
    true_budget = 1 << LOG_D
    out["true_matches_mean"] = "K, at most 0.51 + 2^" + str(k_excess_log2)
    out["true_match_budget_exceeded_off_diagonal_log2_below"] = up(poisson_tail_log2(K, true_budget + 1))
    assert L > match_cap > true_budget
    diag_values = 6 * chains * chains
    diag_prob = diag_values / M
    out["diagonal_event_values"] = "6 x 2^128"
    out["diagonal_event_probability_log2_below"] = up(math.log2(diag_values) - log2_int(M), 2)
    success_floor = 1 - math.exp(-K) - 2 ** -20 - 2 ** -40 - diag_prob
    out["success_after_failure_terms"] = down(success_floor, 6)
    # H1 as stated (proof.md Section 6): coverage |R| / m at least 1 - e^-K - 0.0005, the census deviation allowance
    out["success_after_failure_terms_and_H1_allowance"] = down(success_floor - 0.0005, 6)
    out["success_claimed"] = 0.399
    assert success_floor > 0.399
    # how far below 1 - e^-K the exact coverage |R| / m may fall (H1, mode census) before the claimed 0.399 and the
    # required 0.39 fail
    assert success_floor - 0.0005 > 0.399 and success_floor - 0.0095 > 0.39
    out["coverage_deficit_tolerated"] = {"claimed_0.399": 0.0005, "required_0.39": 0.0095}
    # ledger rows, all in instructions, divided by C at the end
    rows = {}
    rows["records 2N (addition + charged overhead)"] = 2 * N * c_num
    rows["gap tables 2^65 chains x 2^32 additions at C + 47"] = base_count * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD)
    rows["table strand starts 2^65 x (2^16 + 1) scalar mults"] = base_count * ((1 << 16) + 1) * (C << LOG_SCALAR_MULT_UNITS)
    rows["segment starts 2^65 x 2^16 scalar mults"] = base_count * (1 << 16) * (C << LOG_SCALAR_MULT_UNITS)
    rows["per-node fixed work"] = ((1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)) << LOG_NODE_FIXED_INSTR
    rows["comparisons at cap x 8"] = comp_cap * INSTR_PER_COMPARISON
    step_instr = C + WALK_STEP_OVERHEAD + REWALK_STEP_INSTR                        # 2079 per re-walk step
    rewalk = (2 * L * step_instr + 2 * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD)   # both chains re-walked to their ends, both tables rebuilt
              + (4 * (1 << 16) + 2) * (C << LOG_SCALAR_MULT_UNITS))             # strand and segment starts for both
    out["rewalk_step_instr"] = step_instr
    rows["matches at cap: forward 40 B candidates to both generating nodes"] = match_cap * CANDIDATE_FORWARD_INSTR
    rows["matches at cap: re-walk both chains, rebuild both tables"] = match_cap * rewalk
    rows["matches at cap: kept values of the giant chain to the verification node"] = match_cap * KEPT_FORWARD_INSTR
    rows["matches at cap: up to 6 x 4 x 4 candidate scalar mults"] = match_cap * 6 * REWALK_MATCH_CAP ** 2 * (C << LOG_SCALAR_MULT_UNITS)
    rows["gap tables: one F_p inversion per round, 2^65 x 2^16 rounds"] = base_count * table_rounds * inv_p
    rows["matches at cap: one F_p inversion per round of both re-walks and both table rebuilds"] = (
        match_cap * 2 * (walk_rounds + table_rounds) * inv_p)
    rows["matches at cap: F_r inversions and multiplications forming the six z"] = match_cap * z_mults * FP_MUL_INSTR
    # a node's candidate count is Poisson with mean (false + budgeted true) / nodes; the fail event is count > cap
    cand_node_mean_log2 = log2_int(false_mean + true_budget) - LOG_NODES
    out["candidates_per_node_mean_log2"] = up(cand_node_mean_log2, 3)
    out["node_candidate_cap_exceeded_all_nodes_log2_below"] = up((1 << LOG_NODE_CANDIDATE_CAP) * cand_node_mean_log2 + LOG_NODES)
    # each candidate is copied to two generating nodes; a generating node's receive count is Poisson with mean
    # 2 (false + budgeted true) / 2^20 and fails above 2^24 slots
    gen_cand_mean = 2 * (false_mean + true_budget) / (1 << LOG_GEN_NODES)
    out["candidates_per_generating_node_mean_log2"] = up(math.log2(gen_cand_mean), 3)
    out["generating_node_candidate_area_exceeded_all_nodes_log2_below"] = up(
        poisson_tail_log2(gen_cand_mean, (1 << LOG_GEN_CANDIDATE_CAP) + 1) + LOG_GEN_NODES)
    # within-chain key repeats (H3): a chain whose steps carry one key more than REWALK_MATCH_CAP times
    repeat_log2 = (REWALK_MATCH_CAP + 1) * log2_int(L) - REWALK_MATCH_CAP * KEY_BITS - math.log2(math.factorial(REWALK_MATCH_CAP + 1))
    out["rewalk_match_cap_exceeded_all_chains_log2_below"] = up(repeat_log2 + LOG_CHAINS + 1)
    out["within_chain_key_repeat_per_chain_log2"] = up(2 * log2_int(L) - 1 - KEY_BITS, 1)
    rows["walk list W: sieve (2^39, 2^78] x 4 + 2^73 H-tests"] = (
        (1 << LOG_WALK_SIEVE) * 4 + ((1 << LOG_WALK_HTESTS) << LOG_HTEST_INSTR))
    rows["base primes: sieve (2^80, 2^80 + 2^92] x 4 + 2^86 H-tests"] = (
        (1 << LOG_BASE_SIEVE) * 4 + ((1 << LOG_BASE_HTESTS) << LOG_HTEST_INSTR))
    rows["split gaps above 2^33 reserve"] = (C << LOG_BIGGAP_RESERVE)
    # the reserve pays the split additions at C; their walk overhead, and the re-walks' split steps (each re-walked chain
    # may meet the same 16 splits), are charged here at their maximum although the path has expected count below 2^-166.9
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
    g = math.gcd(score_num, 10 ** 6)                                     # reduced form, as the manifest states it
    out["score"] = {"numerator": str(score_num // g), "denominator": str(10 ** 6 // g)}
    assert pow_at_most_pow2(W, 10 ** 6 // g, score_num // g), "W exceeds 2^score"   # W^200000 <= 2^23355819, exactly
    out["score_millionths"] = f"{score_num}/1000000"
    out["score_decimal"] = score_num / 10 ** 6
    out["c_per_record"] = c_num / C
    # preprocessing: everything that does not read A_d (baby side, lists, tables, node setup, the status pass after the
    # baby side, and the baby half of the split-gap reserve with its walk overhead)
    pre = (N * c_num + chains * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD) + chains * (2 * (1 << 16) + 1) * (C << LOG_SCALAR_MULT_UNITS)
           + chains * table_rounds * inv_p
           + rows["per-node fixed work"] + rows["walk list W: sieve (2^39, 2^78] x 4 + 2^73 H-tests"]
           + rows["base primes: sieve (2^80, 2^80 + 2^92] x 4 + 2^86 H-tests"]
           + (all_nodes << LOG_STATUS_INSTR) + ((C + WALK_STEP_OVERHEAD) << (LOG_BIGGAP_RESERVE - 1)))
    pre_units = ceil_div(pre, C)
    out["preprocessing_log2"] = up(log2_int(pre_units), 4)
    assert pre_units ** 10000 <= 1 << PRE_CLAIM_LOG2_E4, "preprocessing exceeds the quoted bound"
    assert pre_units <= W
    # memory: every node at its full itemised allocation, summed in exact integers
    node_bytes = 1 << LOG_NODE_BYTES
    per_owner = {"heads 2^33 x 8 B": 8 << LOG_HEADS, "pool 2^33 x 32 B": 32 << LOG_POOL_RECORDS,
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
    storage_bytes = base_count * 16 + L * 4
    storage_alloc = (1 << LOG_STORE_NODES) << LOG_STORAGE_NODE_BYTES
    assert storage_bytes <= storage_alloc <= (1 << LOG_STORE_NODES) * node_bytes
    out["storage_bytes_log2"] = up(math.log2(storage_bytes), 2)
    mem = (1 << LOG_NODES) * owner_bytes + (1 << LOG_GEN_NODES) * gen_bytes + storage_alloc
    out["memory_bytes"] = str(mem)
    out["memory_log2_bytes_exact_truncated"] = down(log2_int(mem), 6)
    assert mem ** 10000 <= 1 << MEMORY_CLAIM_LOG2_E4, "memory exceeds the claimed bound"
    out["memory_log2_bytes"] = MEMORY_CLAIM_LOG2_E4 / 10000
    out["memory_every_node_at_2^40_would_be_bytes"] = str(((1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)) * node_bytes)
    out["records_alone_log2_bytes"] = up(log2_int(N * 32), 2)
    # H2 sanity on this r: fraction of small primes that are 16th and 256th powers
    small = sieve(1 << 16)
    e256 = [pow(q, (R - 1) // 256, R) for q in small]      # q^((R-1)/16) = (q^((R-1)/256))^16, one exponentiation per prime
    p16 = sum(1 for e in e256 if pow(e, 16, R) == 1)
    p256 = sum(1 for e in e256 if e == 1)
    out["h2_small_primes"] = {"primes_below_2^16": len(small), "16th_powers": p16, "expected": len(small) / 16,
                              "256th_powers": p256, "expected_256": len(small) / 256}
    return out

# ---------------------------------------------------------------- toy curve (j = 0, prime order, 6 * 2^7 | r - 1)
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
    return (y >> 1) & 0xFFFF           # bits 1..16 of the even one of y, p - y

def ec_mode(seed, trials):
    p, b, r, k = TP, TB, TR, TK
    d = 1 << k; m = (r - 1) // d
    assert is_prime(p) and is_prime(r) and (r - 1) % (6 * d) == 0
    hprime = lambda q: pow(q, m, r) == 1
    HP = [q for q in sieve(2_000_000) if hprime(q)]
    omega = next(pow(g, (r - 1) // 3, r) for g in range(2, 50) if pow(g, (r - 1) // 3, r) != 1)
    mu6 = [1, r - 1, omega, r - omega, omega * omega % r, (r - omega * omega) % r]
    assert all(pow(u, m, r) == 1 for u in mu6)
    La, Lw = 15, 25                     # chains per side, steps per chain: n = 375, K = 6 n^2 / m
    n = La * Lw; K = 6 * n * n / m; Kd = 6 * (n * n - n * La + La * La) / m
    Wl = HP[:Lw]; A = HP[Lw:Lw + La]; A2 = HP[Lw + La:Lw + 2 * La]
    gaps = sorted(set(Wl[i + 1] - Wl[i] for i in range(Lw - 1)))
    NODE_BITS, HEAD_BITS = 5, 4
    rnd = random.Random(seed)
    ok = 0; comp_total = 0; false_total = 0; true_total = 0; chi_total = 0.0; chi_sq_total = 0.0; rewalk_cap_fails = 0
    for _ in range(trials):
        P0 = rand_point(p, b, rnd)
        s = rnd.randrange(1, r); T = ec_mul(s, P0, p)
        z = pow(s, d, r); Ad = ec_mul(z, P0, p)
        bits = max(Wl[0], gaps[-1]).bit_length()
        def walk(bases, Q, emit):
            # [a g]Q = [g]([a]Q): one scalar multiplication per chain, then small fixed-base multiples; all chains in lockstep
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
                if rec[0] == cmp: matches.append((key, rec[1], chain))      # the candidate keeps the key and both ids
        tg = walk(A2, Ad, lookup)
        comp_total += comps
        # every candidate (key, id, id'): re-walk each chain once from its start, keep the steps whose key equals the
        # retained key (at most REWALK_MATCH_CAP, else fail), then test the six z for every pair of kept steps
        found = None; cap_fail = False
        def steps_with_key(a, table, key):
            # one linear pass over the chain from its gap table (the full-size ledger charges rebuilding that table too)
            X = table[Wl[0]]; kept = [a * Wl[0] % r] if toy_key(X, p) == key else []
            for j in range(1, Lw):
                X = ec_add(X, table[Wl[j] - Wl[j - 1]], p)
                if toy_key(X, p) == key: kept.append(a * Wl[j] % r)
            assert len(kept) >= 1, "the matched step must carry the key"
            return kept
        for (key, ib, ig) in matches:                  # in the order found, duplicates included, as in production
            xbs = steps_with_key(A[ib], tb[ib], key); xgs = steps_with_key(A2[ig], tg[ig], key)
            if len(xbs) > REWALK_MATCH_CAP or len(xgs) > REWALK_MATCH_CAP:
                cap_fail = True; break                 # the production run fails here; the trial counts as a miss
            hit = False
            for xb in xbs:
                for xg in xgs:
                    cands = [xb * pow(u * xg % r, -1, r) % r for u in mu6]
                    for cand, Rc in zip(cands, batch_mul(cands, [P0] * 6, p)):
                        if Rc == Ad: found = cand; hit = True
            if hit: true_total += 1
            else: false_total += 1
        if cap_fail: rewalk_cap_fails += 1; continue
        if found is None: continue
        assert found == z
        s0 = found
        for _ in range(k): s0 = sqrt_mod(s0, r)
        assert pow(s0, d, r) == z
        g = next(x for x in range(2, 100) if pow(x, (r - 1) // 2, r) != 1)
        zeta = pow(g, (r - 1) // d, r); cur = s0; hit = None
        for _ in range(0, d, 16):                       # the 2^k roots of unity in order, 16 scalar multiplications per batch
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
            "node_chi2_mean": round(chi_mean, 2), "node_chi2_expected": 31, "node_chi2_se": round(chi_se, 2),
            "all_recovered_scalars_correct": True}

# ---------------------------------------------------------------- census: exact |R| / m on toy primes (H1)
# z = s^d is uniform in H and independent of the lists, so the probability that z is an off-diagonal candidate is
# exactly |R| / m, R the set of mu6 multiples of the ratios a w / (a' w') with w != w'. A toy prime r = 6 c 2^k + 1
# (m = (r - 1) / 2^k in [2^b, 2^(b+1)), c from a generator seeded by (b, k, seed)) gets its lists by the production
# rules, and R is counted exactly: the mu6 coset of a w / (a' w') is e_a - e_a' + e_w - e_w' mod m / 6, e the
# discrete logs in H, and each coset is one byte of an array of m / 6 bytes.
CENSUS_SPLIT = (Fraction(64) / Fraction("115.6503")).limit_denominator(10 ** 6)   # log n_b / log N as in production
CENSUS_WALK_START = Fraction(39) / Fraction("77.4")      # log w_min / log w_max as in production (2^39, last 2^77.4)
CENSUS_BASE_OFFSET = Fraction("2.6")                     # first base prime 2^2.6 above the last walk prime (2^80 / 2^77.4)

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
    gamma = pow(g, 1 << k, r)                                   # a generator of H
    e = [x % M for x in dlogs(r, k, gamma, W + A + A2)]
    ew, ea, ea2 = e[:len(W)], e[len(W):len(W) + len(A)], e[len(W) + len(A):]
    alpha = sorted({(x - y) % M for x in ea for y in ea2})     # cosets of a / a'
    delta = sorted({(ew[i] - ew[j]) % M for i in range(len(W)) for j in range(len(W)) if i != j})   # of w / w', w != w'
    flags = bytearray(M)
    outer, inner = (delta, alpha) if len(delta) <= len(alpha) else (alpha, delta)
    for d in outer:
        cut = bisect.bisect_left(inner, M - d)
        for a in inner[:cut]: flags[a + d] = 1
        d -= M
        for a in inner[cut:]: flags[a + d] = 1
    nR = flags.count(1); nRD = sum(1 for a in alpha if flags[a])
    # spot checks: random z in H and random candidates, membership by residue arithmetic against the byte array
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
    Kf = float(Fraction(len(A) * len(A2) * len(W) * (len(W) - 1), M))   # 6 n_b^2 L (L - 1) / m, off-diagonal pairs
    cov, pred = nR / M, -math.expm1(-Kf)
    sd = math.sqrt(math.exp(-Kf) * (1 - (1 + Kf) * math.exp(-Kf)) / M)
    return {"r": r, "k": k, "seed": seed, "m": m, "m_log2": "%.4f" % math.log2(m), "nb": len(A), "L": len(W),
            "N": len(A) * len(W), "K": "%.6f" % Kf, "R_cosets": nR, "M6": M, "cov": "%.8f" % cov, "pred": "%.8f" % pred,
            "diff": "%+.8f" % (cov - pred), "z": "%+.2f" % ((cov - pred) / sd), "D_frac": "%.8f" % (len(alpha) / M),
            "succ": "%.8f" % ((nR - nRD) / M), "spot": "%d/40" % ok}

def census(m_bits, ks, K_percent, first_seed, primes, seconds):
    """Primes in the order seed, then k; no further prime starts once `seconds` have passed (hosted runs)."""
    K = Fraction(K_percent, 100); t0 = time.monotonic(); rows = []; stopped = False
    for seed in range(first_seed, first_seed + primes):
        for k in ks:
            if rows and time.monotonic() - t0 > seconds: stopped = True; break
            rows.append(census_row(m_bits, k, K, seed))
        if stopped: break
    d = [float(x["diff"]) for x in rows]
    return {"m_bits": m_bits, "k": ks, "K_target": K_percent / 100, "first_seed": first_seed,
            "primes_requested": primes * len(ks), "primes_run": len(rows), "stopped_by_time_budget": stopped,
            "diff_is": "|R| / m - (1 - e^-K) at the K each prime reached (K >= the target)",
            "z_is": "diff over the spread that independent uniform draws would give; the spread across primes is several "
                    "times that (proof.md Section 6), so z measures one prime's deviation, not a test of H1",
            "mean_diff": "%+.6f" % (sum(d) / len(d)), "lowest_diff": "%+.6f" % min(d), "highest_diff": "%+.6f" % max(d),
            "below_prediction": sum(1 for x in d if x < 0), "rows": rows}

def main():
    raw = sys.stdin.read().strip() if not sys.stdin.isatty() else ""
    args = json.loads(raw) if raw else {}
    mode = args.get("mode", "all"); seed = int(args.get("seed", 2026)); trials = int(args.get("trials", 400))
    assert mode in ("ledger", "ec", "census", "all") and 0 <= seed < 2 ** 31 and 20 <= trials <= 2000
    cen = (int(args.get("m_bits", 21)), [int(k) for k in args.get("k", [2, 4, 7, 10])], int(args.get("K_percent", 51)),
           int(args.get("first_seed", 0)), int(args.get("primes", 2)), int(args.get("seconds", 8)))
    mb, ks, kp, fs, npr, sec = cen
    assert (16 <= mb <= 36 and 1 <= len(ks) <= 4 and len(set(ks)) == len(ks) and all(1 <= k <= 12 for k in ks)
            and 30 <= kp <= 90 and 0 <= fs < 2 ** 31 and 1 <= npr <= 64 and 1 <= sec <= 10 ** 6)
    out = {}
    if mode in ("ledger", "all"): out["ledger"] = ledger()
    if mode in ("ec", "all"): out["ec"] = ec_mode(seed, trials)
    if mode in ("census", "all"): out["census"] = census(*cen)
    json.dump(out, sys.stdout, indent=1); print()

if __name__ == "__main__": main()

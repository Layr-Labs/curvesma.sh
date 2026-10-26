#!/usr/bin/env python3
import sys, json, math, random

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
LOG_D = 20
D = 1 << LOG_D
M = (R - 1) // D
C = 2000
K_NUM, K_DEN = 52, 100
LOG_CHAINS = 64
LOG_NODES = 84
HOP_BITS, HOPS = 28, 3
LOG_HEADS = 33
LOG_NODE_BYTES = 40
LOG_POOL_RECORDS = 33
LOG_GAP_TABLE = 32
LOG_GEN_NODES, LOG_STORE_NODES = 20, 30
KEY_BITS = 192
OVERHEAD_LISTED_BABY = {"gap stream": 5, "table address": 2, "operand loads and stores": 18, "lane control": 8,
                        "parity key": 23, "scatter 3 hops": 36, "post and poll 3 hops": 9, "transport copy 3 hops": 18,
                        "gap stream transport from storage": 3, "batch inversion loads and stores": 6,
                        "flow control and counters": 5, "insert": 10}
OVERHEAD_LISTED_GIANT = dict(OVERHEAD_LISTED_BABY); del OVERHEAD_LISTED_GIANT["insert"]; OVERHEAD_LISTED_GIANT["lookup fixed"] = 6
OVERHEAD_CHARGED = 180
WALK_STEP_OVERHEAD = sum(OVERHEAD_LISTED_BABY[k] for k in ("gap stream", "gap stream transport from storage", "table address",
                         "operand loads and stores", "batch inversion loads and stores", "lane control", "flow control and counters"))
assert WALK_STEP_OVERHEAD == 47
INSTR_PER_COMPARISON = 8
SLOT_CHAIN_CAP = 64
BERNSTEIN_VAR_MARGIN = 4
LOG_NODE_FIXED_INSTR = 35
LOG_SCALAR_MULT_UNITS = 12
LOG_HTEST_INSTR = 18
COMPARISON_CAP_NUM, COMPARISON_CAP_DEN = 1002, 1000
LOG_MATCH_CAP = 44
CANDIDATE_BYTES = 40
LOG_NODE_CANDIDATE_CAP = 16
CANDIDATE_FORWARD_INSTR = 128
LOG_GEN_CANDIDATE_CAP = 24
REWALK_STEP_INSTR = 32
REWALK_MATCH_CAP = 4
KEPT_FORWARD_INSTR = 256
LOG_STATUS_INSTR = 12
PRE_CLAIM_LOG2_E4 = 1157954
MEMORY_CLAIM_LOG2_E4 = 1233242
LOG_STORAGE_NODE_BYTES = 40
LOG_BIGGAP_RESERVE = 69
LOG_STAGE2 = 24

def ceil_div(a, b): return -(-a // b)

def up(x, nd=0):
    """Round a log2 probability bound towards zero (the safe side): the text says 'below 2^-k' with k = -up(x)."""
    return math.ceil(x * 10 ** nd) / 10 ** nd

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
    out["steps_per_chain_log2"] = round(log2_int(L), 4)
    out["records_per_side_log2"] = round(log2_int(N), 4)
    K = 6 * N * N / M
    out["K"] = round(K, 6)
    out["predicted_success"] = round(1 - math.exp(-K), 4)
    out["diagonal_loss_log2"] = -round(log2_int(L), 2)
    hp_needed = L
    primes_needed = hp_needed * D
    x = primes_needed * 54
    out["walk_list_last_prime_log2_about"] = round(log2_int(x), 2)
    assert x < 1 << 79
    base_count = 2 * chains
    out["base_primes_needed_log2"] = LOG_CHAINS + 1
    out["x_max_log2_below"] = 99 + 79; assert (1 << 178) < R
    out["mean_gap_log2"] = round(math.log2(D * 54), 2)
    out["gap_above_2^33_expected_log2"] = round(log2_int(L) - (1 << 33) / (D * 54) * math.log2(math.e), 1)
    lb, lg = sum(OVERHEAD_LISTED_BABY.values()), sum(OVERHEAD_LISTED_GIANT.values())
    out["overhead_listed"] = {"baby": lb, "giant": lg, "charged": OVERHEAD_CHARGED}
    assert OVERHEAD_CHARGED >= max(lb, lg)
    c_num = C + OVERHEAD_CHARGED
    nodes = 1 << LOG_NODES
    mean_per_node = N / nodes
    out["baby_records_per_node_mean_log2"] = round(math.log2(mean_per_node), 3)
    cap = 1 << LOG_POOL_RECORDS
    delta = cap / mean_per_node - 1
    out["pool_capacity_over_mean"] = round(cap / mean_per_node, 2)
    chernoff_log2 = -mean_per_node * delta * delta / (2 + delta) * math.log2(math.e)
    out["overflow_per_node_chernoff_log2"] = up(chernoff_log2)
    heads = nodes << LOG_HEADS
    lam = N / heads
    out["load_per_head"] = round(lam, 3)
    comp_mean = N * N // heads
    node_comp_mean = comp_mean // nodes
    node_comp_cap = ceil_div(COMPARISON_CAP_NUM * node_comp_mean, COMPARISON_CAP_DEN)
    comp_cap = node_comp_cap * nodes
    out["comparisons_mean_log2"] = round(log2_int(comp_mean), 3)
    out["comparisons_per_node_cap_log2"] = round(log2_int(node_comp_cap), 3)
    k = SLOT_CHAIN_CAP
    k1 = k + 1
    slot_tail_log2 = (-lam * math.log2(math.e) + k1 * math.log2(lam) - math.log2(math.factorial(k1))
                      - math.log2(1 - lam / (k1 + 1)))
    out["slot_chain_cap_exceeded_all_slots_log2"] = up(slot_tail_log2 + math.log2(heads), 1)
    t = 0.001 * node_comp_mean
    V = mean_per_node * (lam + lam * lam) * BERNSTEIN_VAR_MARGIN
    bern_log2 = -(t * t / (2 * (V + k * t / 3))) * math.log2(math.e)
    out["comparisons_node_cap_exceeded_bernstein_all_nodes_log2"] = up(bern_log2 + LOG_NODES)
    count_dev_log2 = -(mean_per_node * (2 ** -11) ** 2 / 3) * math.log2(math.e) + LOG_NODES + 1
    out["node_count_off_by_2^-11_chernoff_all_nodes_log2"] = up(count_dev_log2)
    false_mean = N * N // (1 << KEY_BITS)
    out["false_matches_mean_log2"] = round(log2_int(false_mean), 2)
    match_cap = 1 << LOG_MATCH_CAP
    out["match_cap_exceeded_chebyshev_log2"] = up(math.log2(2 * false_mean + 1) - log2_int((match_cap - false_mean) ** 2), 1)
    true_budget = 1 << LOG_D
    out["true_matches_mean"] = round(K, 4)
    out["true_match_budget_exceeded_log2_below"] = up(poisson_tail_log2(K, true_budget + 1))
    out["success_claimed"] = 0.40
    assert 1 - math.exp(-K) - 2 ** -20 - 2 ** -40 > 0.40
    rows = {}
    rows["records 2N (addition + charged overhead)"] = 2 * N * c_num
    rows["gap tables 2^65 chains x 2^32 additions at C + 47"] = base_count * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD)
    rows["table strand starts 2^65 x (2^16 + 1) scalar mults"] = base_count * ((1 << 16) + 1) * (C << LOG_SCALAR_MULT_UNITS)
    rows["segment starts 2^65 x 2^16 scalar mults"] = base_count * (1 << 16) * (C << LOG_SCALAR_MULT_UNITS)
    rows["per-node fixed work"] = ((1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)) << LOG_NODE_FIXED_INSTR
    rows["comparisons at cap x 8"] = comp_cap * INSTR_PER_COMPARISON
    step_instr = C + WALK_STEP_OVERHEAD + REWALK_STEP_INSTR
    rewalk = (2 * L * step_instr + 2 * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD)
              + (4 * (1 << 16) + 2) * (C << LOG_SCALAR_MULT_UNITS))
    out["rewalk_step_instr"] = step_instr
    rows["matches at cap: forward 40 B candidates to both generating nodes"] = match_cap * CANDIDATE_FORWARD_INSTR
    rows["matches at cap: re-walk both chains, rebuild both tables"] = match_cap * rewalk
    rows["matches at cap: kept values of the giant chain to the verification node"] = match_cap * KEPT_FORWARD_INSTR
    rows["matches at cap: up to 6 x 4 x 4 candidate scalar mults"] = match_cap * 6 * REWALK_MATCH_CAP ** 2 * (C << LOG_SCALAR_MULT_UNITS)
    cand_node_mean_log2 = log2_int(false_mean + true_budget) - LOG_NODES
    out["candidates_per_node_mean_log2"] = round(cand_node_mean_log2, 1)
    out["node_candidate_cap_exceeded_all_nodes_log2_below"] = up((1 << LOG_NODE_CANDIDATE_CAP) * cand_node_mean_log2 + LOG_NODES)
    gen_cand_mean = 2 * (false_mean + true_budget) / (1 << LOG_GEN_NODES)
    out["candidates_per_generating_node_mean_log2"] = round(math.log2(gen_cand_mean), 2)
    out["generating_node_candidate_area_exceeded_all_nodes_log2_below"] = up(
        poisson_tail_log2(gen_cand_mean, (1 << LOG_GEN_CANDIDATE_CAP) + 1) + LOG_GEN_NODES)
    repeat_log2 = (REWALK_MATCH_CAP + 1) * log2_int(L) - REWALK_MATCH_CAP * KEY_BITS - math.log2(math.factorial(REWALK_MATCH_CAP + 1))
    out["rewalk_match_cap_exceeded_all_chains_log2_below"] = up(repeat_log2 + LOG_CHAINS + 1)
    out["within_chain_key_repeat_per_chain_log2"] = round(2 * log2_int(L) - 1 - KEY_BITS, 1)
    rows["walk list W: sieve 2^78 x 4 + 2^73 H-tests"] = (1 << 78) * 4 + ((1 << 73) << LOG_HTEST_INSTR)
    rows["base primes: sieve 2^92 x 4 + 2^86 H-tests"] = (1 << 92) * 4 + ((1 << 86) << LOG_HTEST_INSTR)
    rows["split gaps above 2^33 reserve"] = (C << LOG_BIGGAP_RESERVE)
    all_nodes = (1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)
    rows["status tree: 3 passes x every node x 2^12"] = 3 * (all_nodes << LOG_STATUS_INSTR)
    rows["stage 2 up to 2^20 + 1 times, input access, output"] = (true_budget + 1) * (C << LOG_STAGE2)
    total_instr = sum(rows.values())
    W = ceil_div(total_instr, C)
    out["rows_log2_units"] = {k: round(log2_int(ceil_div(v, C)), 3) for k, v in rows.items()}
    out["W_units"] = str(W)
    lw = log2_int(W)
    score_num = math.ceil((lw + 1e-9) * 10 ** 6)
    out["log2_W"] = round(lw, 7)
    g = math.gcd(score_num, 10 ** 6)
    out["score"] = {"numerator": str(score_num // g), "denominator": str(10 ** 6 // g)}
    out["score_millionths"] = f"{score_num}/1000000"
    out["score_decimal"] = score_num / 10 ** 6
    out["c_per_record"] = c_num / C
    pre = (N * c_num + chains * (1 << LOG_GAP_TABLE) * (C + WALK_STEP_OVERHEAD) + chains * (2 * (1 << 16) + 1) * (C << LOG_SCALAR_MULT_UNITS)
           + rows["per-node fixed work"] + rows["walk list W: sieve 2^78 x 4 + 2^73 H-tests"]
           + rows["base primes: sieve 2^92 x 4 + 2^86 H-tests"])
    pre_units = ceil_div(pre, C)
    out["preprocessing_log2"] = math.ceil(log2_int(pre_units) * 10 ** 4) / 10 ** 4
    assert out["preprocessing_log2"] * 10000 <= PRE_CLAIM_LOG2_E4, "preprocessing exceeds the quoted bound"
    assert pre_units ** 10000 <= 1 << PRE_CLAIM_LOG2_E4, "preprocessing exceeds the quoted bound (integer check)"
    assert pre_units <= W
    node_bytes = 1 << LOG_NODE_BYTES
    per_owner = {"heads 2^33 x 8 B": 8 << LOG_HEADS, "pool 2^33 x 32 B": 32 << LOG_POOL_RECORDS,
                 "send blocks 2^28 x 512 B": 512 << HOP_BITS, "receive FIFOs 3 x 2^36 B": 3 << 36,
                 "candidate list 2^16 x 40 B": CANDIDATE_BYTES << LOG_NODE_CANDIDATE_CAP, "code and counters": 1 << 30}
    owner_bytes = sum(per_owner.values())
    assert owner_bytes <= node_bytes
    out["owner_node_bytes_log2"] = round(math.log2(owner_bytes), 4)
    per_gen = {"gap table 2^32 x 96 B": 96 << LOG_GAP_TABLE, "send blocks 2^28 x 512 B": 512 << HOP_BITS,
               "stream buffers and state": 1 << 36, "candidate receive area 2^24 x 40 B": CANDIDATE_BYTES << LOG_GEN_CANDIDATE_CAP}
    gen_bytes = sum(per_gen.values())
    assert gen_bytes <= node_bytes
    out["generating_node_bytes_log2"] = round(math.log2(gen_bytes), 4)
    storage_bytes = base_count * 16 + L * 4
    storage_alloc = (1 << LOG_STORE_NODES) << LOG_STORAGE_NODE_BYTES
    assert storage_bytes <= storage_alloc <= (1 << LOG_STORE_NODES) * node_bytes
    out["storage_bytes_log2"] = round(math.log2(storage_bytes), 2)
    mem = (1 << LOG_NODES) * owner_bytes + (1 << LOG_GEN_NODES) * gen_bytes + storage_alloc
    out["memory_bytes"] = str(mem)
    out["memory_log2_bytes_exact_about"] = round(log2_int(mem), 6)
    assert mem ** 10000 <= 1 << MEMORY_CLAIM_LOG2_E4, "memory exceeds the claimed bound"
    out["memory_log2_bytes"] = MEMORY_CLAIM_LOG2_E4 / 10000
    out["memory_every_node_at_2^40_would_be_bytes"] = str(((1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)) * node_bytes)
    out["records_alone_log2_bytes"] = round(log2_int(N * 32), 2)
    small = sieve(1 << 16)
    e256 = [pow(q, (R - 1) // 256, R) for q in small]
    p16 = sum(1 for e in e256 if pow(e, 16, R) == 1)
    p256 = sum(1 for e in e256 if e == 1)
    out["h2_small_primes"] = {"primes_below_2^16": len(small), "16th_powers": p16, "expected": round(len(small) / 16, 1),
                              "256th_powers": p256, "expected_256": round(len(small) / 256, 1)}
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
        found = None; cap_fail = False
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
            "K": round(K, 4), "predicted": round(pred, 4), "diagonal_corrected": round(predd, 4),
            "measured": round(rate, 4), "trials": trials, "z_vs_corrected": round((rate - predd) / se, 2),
            "comparisons_mean": round(comp_total / trials, 1), "comparisons_expected": round(n * n / 512, 1),
            "false_matches_mean": round(false_total / trials, 3), "false_matches_expected": round(n * n / 65536, 3),
            "true_matches_mean": round(true_total / trials, 3), "rewalk_cap_fails": rewalk_cap_fails,
            "node_chi2_mean": round(chi_mean, 2), "node_chi2_expected": 31, "node_chi2_se": round(chi_se, 2),
            "all_recovered_scalars_correct": True}

def main():
    raw = sys.stdin.read().strip() if not sys.stdin.isatty() else ""
    args = json.loads(raw) if raw else {}
    mode = args.get("mode", "all"); seed = int(args.get("seed", 2026)); trials = int(args.get("trials", 400))
    assert mode in ("ledger", "ec", "all") and 0 <= seed < 2 ** 31 and 20 <= trials <= 2000
    out = {}
    if mode in ("ledger", "all"): out["ledger"] = ledger()
    if mode in ("ec", "all"): out["ec"] = ec_mode(seed, trials)
    json.dump(out, sys.stdout, indent=1); print()

if __name__ == "__main__": main()

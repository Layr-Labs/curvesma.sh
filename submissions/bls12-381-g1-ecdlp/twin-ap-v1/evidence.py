#!/usr/bin/env python3
import sys, json, math, random, decimal

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
LOG_D = 20
D = 1 << LOG_D
M = (R - 1) // D
C = 2000
K_NUM, K_DEN = 4944, 10000
LOG_CHAINS = 56
LOG_X0 = 210
LOG_SPAN_CAP = 99
LOG_T = 44
LOG_MULT_K0, LOG_MULT_RANGE = 103, 85
LOG_MULT_LO = 211
LOG_NODES = 80
HOP_BITS, HOPS = 40, 2
LOG_HEADS = 37
LOG_NODE_BYTES = 43
LOG_STORE_NODE_BYTES = 40
LOG_GEN_NODE_BYTES = 52
LOG_POOL_RECORDS = 37
LOG_GEN_NODES, LOG_STORE_NODES = 20, 30
KEY_BITS = 192
DELIVERY_LISTED = {"sender: credit check and send-buffer reclaim": 6, "sender: descriptor (3 stores), doorbell, addressing": 8,
                   "interface copy of the 576-byte slot charged as 36 ldp + 36 stp plus loop control": 80,
                   "interface: arrival-queue append and event at the receiver": 3,
                   "receiver: wake on the data put (at most one per put)": 8,
                   "receiver: arrival-queue read, slot address, record count from the slot header": 7,
                   "receiver: credit return, an 8-byte put (descriptor, doorbell, interface copy charged)": 10,
                   "sender: wake on the credit put (at most one per put)": 8, "sender: credit update and resume check": 6}
DELIVERY_PER_BLOCK = 144
SUSPEND_PER_RECORD_HOP = 8
assert sum(DELIVERY_LISTED.values()) <= DELIVERY_PER_BLOCK
RECORDS_PER_BLOCK, STEPS_PER_STREAM_BLOCK = 16, 85
SLOT_BYTES = 576
LOG_RELAY_NODES = 40
LOG_RELAY_FIXED_INSTR = 50
WALK_OVERHEAD_LISTED = {"gap stream read and pointer advance": 5, "gap stream transport (6 B per step, blocked)": 4,
                        "gap stream delivery (144 per block of 85 steps, 1.7)": 2, "stream request (144 per 85-step block, 1.7)": 2,
                        "table address": 2, "operand loads and stores": 18, "batch inversion prefix/suffix loads and stores": 6,
                        "lane control and loop": 8, "flow control and counters": 5}
SIDE_ARITH_LISTED = {"2 multiplications + 1 squaring at 383 (blst_fp_mul, measured)": 3 * 383,
                     "1 addition + 4 subtractions mod p at 48": 5 * 48,
                     "argument setup for 8 calls at 4": 32, "keep the lane inverse: 3 stp + 3 ldp": 6}
RECORD_LISTED_BABY = {"parity key from the stored Montgomery residue": 23, "full reduction of the residue to [0, p) (6-limb compare, conditional subtract)": 12,
                      "scatter 2 hops": 24, "block bookkeeping 2 hops": 6, "transport copy 2 hops": 12, "insert": 10,
                      "delivery 2 hops (144 per block per hop, 16 records)": HOPS * DELIVERY_PER_BLOCK // RECORDS_PER_BLOCK,
                      "suspend/resume 2 hops (at most one episode per record per hop)": HOPS * SUSPEND_PER_RECORD_HOP}
RECORD_LISTED_GIANT = dict(RECORD_LISTED_BABY); del RECORD_LISTED_GIANT["insert"]; RECORD_LISTED_GIANT["lookup fixed"] = 6
WALK_STEP_OVERHEAD = sum(WALK_OVERHEAD_LISTED.values())
SIDE_ARITH = sum(SIDE_ARITH_LISTED.values())
assert WALK_STEP_OVERHEAD == 52 and SIDE_ARITH == 1427
STEP_EXTRA_CHARGED = 1760
FALLBACK_STEP_INSTR = 2200
assert FALLBACK_STEP_INSTR >= C + WALK_STEP_OVERHEAD + sum(RECORD_LISTED_BABY.values())
INSTR_PER_COMPARISON = 8
SLOT_CHAIN_CAP = 64
BERNSTEIN_VAR_MARGIN = 4
LOG_NODE_FIXED_INSTR = 37
LOG_GEN_FIXED_INSTR = 50
LOG_SCALAR_MULT_UNITS = 12
LOG_HTEST_INSTR = 17
LOG_POCK_INSTR = 24
LOG_POCK_CAP = 66
LOG_CU_CAP = 81
PRE_CLAIM_LOG2_E4 = 1155456
LOG_SEGMENTS = 16
COMPARISON_CAP_NUM, COMPARISON_CAP_DEN = 101, 100
LOG_MATCH_CAP = 44
CANDIDATE_BYTES = 40
LOG_NODE_CANDIDATE_CAP = 16
CANDIDATE_FORWARD_INSTR = 1024
LOG_GEN_CANDIDATE_CAP = 24
REWALK_MATCH_CAP = 4
KEPT_FORWARD_INSTR = 512
LOG_STATUS_INSTR = 12
MEMORY_CLAIM_LOG2_E4 = 1223223
LOG_FALLBACK_RESERVE = 40
LIST_ENTRY_BYTES = 6
LOG_STAGE2 = 24

FIXED_ROW = "per-node fixed work: owner and storage 2^37, relay roles 2^50, generating 2^50"
MULT_ROW = "multipliers: H-test 2^85 candidates 1 + 2kF, Pocklington certificate on at most 2^66 survivors"
LIST_ROW = "twin list S: (2^99 + 2^44 + 2^81 + 2^84) H-tests at 2^17 (window duplicate checks included)"

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

def trial_prime(n):
    """Primality by trial division (the bottom of a certificate chain and the toy curve constants, all below 2^30)."""
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def pocklington(N, F):
    """Pocklington-Lehmer certificate. Preconditions: F prime, F | N - 1, F^2 > N. If some base b has b^(N-1) = 1 mod N
    and gcd(b^((N-1)/F) - 1, N) = 1, then every prime factor q of N has q = 1 mod F, so q > F > sqrt(N) and N is
    prime. Returns True only with such a proof (never for a composite N); a prime N fails all four bases only if
    each base is an F-th power residue mod N (probability about F^-4)."""
    assert (N - 1) % F == 0 and F * F > N
    for b in (2, 3, 5, 7):
        if pow(b, N - 1, N) != 1: return False
        g = math.gcd(pow(b, (N - 1) // F, N) - 1, N)
        if g == 1: return True
        if g != N: return False
    return False

def certified_prime_chain(bits0, bits):
    """p_0 = least prime >= 2^bits0 (trial division); p_{i+1} = least N >= 2^bits[i] with N = 1 + 2k p_i and a
    Pocklington certificate over p_i (p_i^2 > N asserted). Every element is proven prime."""
    p = 1 << bits0
    while not trial_prime(p): p += 1
    chain = [p]
    for b in bits:
        k = -(-(1 << b) // (2 * p))
        while not pocklington(1 + 2 * k * p, p): k += 1
        p = 1 + 2 * k * p; chain.append(p)
    return chain

def build_twin_list(start, steps, inU, T, lo):
    """The list rule of proof.md Section 3. From u_0 = start, each step picks the least c in (u, u + T] with c in U,
    2u - c in U, 2u - c >= lo and 2u - c not among the list elements in the window [u - T, u); it emits the main
    u' = c and the side v = 2u - c. If no c qualifies (a fallback step) it takes the least c in (u, u + T] with c in U
    and emits only the main. Returns the gap list [(g, twin)] and the element list S. Only the elements in the
    current window are kept for the duplicate check (every element of S below u - T can never equal a later side)."""
    u = start; window = [u]; gaps = []; S = []
    for _ in range(steps):
        c = u + 1; pick = None; single = None
        while c <= u + T:
            if inU(c):
                if single is None: single = c
                v = 2 * u - c
                if v >= lo and v not in window and inU(v): pick = c; break
            c += 1
        if pick is not None:
            v = 2 * u - pick; gaps.append((pick - u, True)); S += [pick, v]; window += [pick, v]; u = pick
        else:
            assert single is not None, "no element of U within T"
            gaps.append((single - u, False)); S.append(single); window.append(single); u = single
        window = [w for w in window if w >= u - T]
    return gaps, S

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
    L = ceil_div(n0, 2 * chains)
    N = 2 * L * chains
    assert 6 * N * N * K_DEN >= K_NUM * M
    out["twin_steps_per_chain_log2"] = round(log2_int(L), 4)
    out["records_per_side_log2"] = round(log2_int(N), 4)
    K = 6 * N * N / M
    out["K"] = round(K, 6)
    Kd = K * (1 - 1 / (2 * L))
    out["predicted_success"] = round(1 - math.exp(-Kd), 5)
    T = 1 << LOG_T
    span_mean = L << 40
    out["list_span_mean_log2"] = round(log2_int(span_mean), 3)
    span_cap = 1 << LOG_SPAN_CAP
    out["list_span_cap_log2"] = LOG_SPAN_CAP
    delta = 2 * T + span_cap
    assert delta * delta < (1 << LOG_X0)
    out["delta_log2"] = round(log2_int(delta), 4)
    eps = span_cap / span_mean - 1
    span_tail_log2 = -L * (eps - math.log(1 + eps)) * math.log2(math.e)
    out["list_span_over_cap_chernoff_log2_below"] = up(span_tail_log2) if span_tail_log2 > -1e300 else "-inf"
    p_fallback_log2 = -T / 2 ** 40 * math.log2(math.e)
    out["fallback_per_step_log2"] = round(p_fallback_log2, 2)
    fallback_mean_log2 = log2_int(L) + p_fallback_log2
    out["fallback_steps_mean_log2"] = round(fallback_mean_log2, 2)
    fallback_cap = 1 << LOG_FALLBACK_RESERVE
    out["fallback_over_reserve_log2_below"] = up(poisson_tail_log2(2 ** fallback_mean_log2, fallback_cap + 1))
    chain = certified_prime_chain(28, [55, 107])
    F = chain[-1]
    out["certificate_chain"] = [str(q) for q in chain]
    out["F_log2"] = round(math.log2(F), 4)
    k_hi = (1 << LOG_MULT_K0) + (1 << LOG_MULT_RANGE)
    assert 1 + 2 * (1 << LOG_MULT_K0) * F > (1 << LOG_MULT_LO) and F * F > 1 + 2 * k_hi * F
    demo, k = [], 1 << LOG_MULT_K0
    while len(demo) < 3:
        if pocklington(1 + 2 * k * F, F): demo.append(k - (1 << LOG_MULT_K0))
        k += 1
    out["first_certified_multiplier_candidates_k_minus_2^103"] = demo
    assert (1 << LOG_MULT_LO) > (1 << LOG_X0) + delta
    mult_needed = 2 * chains
    hprime_density_log2 = -LOG_D + math.log2(2 / ((LOG_MULT_LO + 1) * math.log(2)))
    out["multiplier_hprimes_expected_log2"] = round(LOG_MULT_RANGE + hprime_density_log2, 2)
    assert LOG_MULT_RANGE + hprime_density_log2 > LOG_CHAINS + 1 + 1.5
    lb = WALK_STEP_OVERHEAD + SIDE_ARITH + 2 * sum(RECORD_LISTED_BABY.values())
    lg = WALK_STEP_OVERHEAD + SIDE_ARITH + 2 * sum(RECORD_LISTED_GIANT.values())
    out["twin_step_extra_listed"] = {"baby": lb, "giant": lg, "charged": STEP_EXTRA_CHARGED,
                                     "walk_overhead": WALK_STEP_OVERHEAD, "side_arithmetic": SIDE_ARITH}
    assert STEP_EXTRA_CHARGED >= max(lb, lg)
    step_instr = C + STEP_EXTRA_CHARGED
    out["instr_per_record"] = step_instr / 2
    nodes = 1 << LOG_NODES
    mean_per_node = N / nodes
    out["baby_records_per_node_mean_log2"] = round(math.log2(mean_per_node), 3)
    cap = 1 << LOG_POOL_RECORDS
    dlt = cap / mean_per_node - 1
    out["pool_capacity_over_mean"] = round(cap / mean_per_node, 2)
    out["overflow_any_node_chernoff_log2"] = up(-mean_per_node * dlt * dlt / (2 + dlt) * math.log2(math.e) + LOG_NODES)
    heads = nodes << LOG_HEADS
    lam = N / heads
    out["load_per_head"] = round(lam, 3)
    comp_mean = N * N // heads
    node_comp_mean = comp_mean // nodes
    node_comp_cap = ceil_div(COMPARISON_CAP_NUM * node_comp_mean, COMPARISON_CAP_DEN)
    comp_cap = node_comp_cap * nodes
    out["comparisons_mean_log2"] = round(log2_int(comp_mean), 3)
    out["comparisons_per_node_cap_log2"] = round(log2_int(node_comp_cap), 3)
    k = SLOT_CHAIN_CAP; k1 = k + 1
    slot_tail_log2 = (-lam * math.log2(math.e) + k1 * math.log2(lam) - math.log2(math.factorial(k1))
                      - math.log2(1 - lam / (k1 + 1)))
    out["slot_chain_cap_exceeded_all_slots_log2"] = up(slot_tail_log2 + math.log2(heads), 1)
    t = 0.008 * node_comp_mean
    V = mean_per_node * (lam + lam * lam) * BERNSTEIN_VAR_MARGIN
    bern_log2 = -(t * t / (2 * (V + k * t / 3))) * math.log2(math.e)
    out["comparisons_node_cap_exceeded_bernstein_all_nodes_log2"] = up(bern_log2 + LOG_NODES)
    count_dev_log2 = -(mean_per_node * (2 ** -10) ** 2 / 3) * math.log2(math.e) + LOG_NODES + 1
    out["node_count_off_by_2^-10_chernoff_all_nodes_log2"] = up(count_dev_log2)
    false_mean = N * N // (1 << KEY_BITS)
    out["false_matches_mean_log2"] = round(log2_int(false_mean), 2)
    match_cap = 1 << LOG_MATCH_CAP
    out["match_cap_exceeded_chebyshev_log2"] = up(math.log2(2 * false_mean + 1) - log2_int((match_cap - (1 << LOG_D) - false_mean) ** 2), 1)
    true_budget = 1 << LOG_D
    out["true_matches_mean"] = round(K, 4)
    out["true_match_budget_exceeded_log2_below"] = up(poisson_tail_log2(K, true_budget + 1))
    cand_node_mean_log2 = log2_int(false_mean + true_budget) - LOG_NODES
    out["candidates_per_node_mean_log2"] = round(cand_node_mean_log2, 1)
    out["node_candidate_cap_exceeded_all_nodes_log2_below"] = up((1 << LOG_NODE_CANDIDATE_CAP) * cand_node_mean_log2 + LOG_NODES)
    gen_cand_mean = 2 * (false_mean + true_budget) / (1 << LOG_GEN_NODES)
    out["candidates_per_generating_node_mean_log2"] = round(math.log2(gen_cand_mean), 2)
    out["generating_node_candidate_area_exceeded_all_nodes_log2_below"] = up(
        poisson_tail_log2(gen_cand_mean, (1 << LOG_GEN_CANDIDATE_CAP) + 1) + LOG_GEN_NODES)
    outputs_per_chain = 2 * L
    repeat_log2 = ((REWALK_MATCH_CAP + 1) * log2_int(outputs_per_chain) - REWALK_MATCH_CAP * KEY_BITS
                   - math.log2(math.factorial(REWALK_MATCH_CAP + 1)))
    out["rewalk_match_cap_exceeded_all_chains_log2_below"] = up(repeat_log2 + LOG_CHAINS + 1)
    out["within_chain_key_repeat_per_chain_log2"] = round(2 * log2_int(outputs_per_chain) - 1 - KEY_BITS, 1)
    tails = [out["slot_chain_cap_exceeded_all_slots_log2"], out["comparisons_node_cap_exceeded_bernstein_all_nodes_log2"],
             out["node_count_off_by_2^-10_chernoff_all_nodes_log2"], out["match_cap_exceeded_chebyshev_log2"],
             out["true_match_budget_exceeded_log2_below"], out["node_candidate_cap_exceeded_all_nodes_log2_below"],
             out["generating_node_candidate_area_exceeded_all_nodes_log2_below"], out["rewalk_match_cap_exceeded_all_chains_log2_below"]]
    assert max(tails) < -47 and sum(2.0 ** t for t in tails) <= 2 ** -40
    assert out["overflow_any_node_chernoff_log2"] < -20
    out["cap_tails_sum_log2_below"] = up(math.log2(sum(2.0 ** t for t in tails)), 1)
    fail_terms = 2 ** -20 + 2 ** -40
    success = 1 - math.exp(-Kd) - fail_terms
    out["success_after_failure_terms"] = round(success, 5)
    out["success_claimed"] = 0.390
    assert success >= 0.390
    S = 1 << LOG_SEGMENTS
    smul = C << LOG_SCALAR_MULT_UNITS
    rows = {}
    rows[f"records: 2 x 2^56 chains x L twin steps (one reference addition + {STEP_EXTRA_CHARGED} instr, two records each)"] = 2 * chains * L * step_instr
    rows["fallback steps: reserve 2^40 per list x 2^57 chains, one single-record step charge each"] = (2 * chains) * fallback_cap * FALLBACK_STEP_INSTR
    rows[f"gap tables: 2^57 chains x 2^44 additions at C + {WALK_STEP_OVERHEAD}"] = 2 * chains * T * (C + WALK_STEP_OVERHEAD)
    rows["table strand starts: 2^57 x 2^16 scalar mults"] = 2 * chains * S * smul
    rows["segment starts: 2^57 x 2^16 scalar mults"] = 2 * chains * S * smul
    rows["chain bases [a]Q: 2^57 scalar mults"] = 2 * chains * smul
    rows[LIST_ROW] = (span_cap + T + (1 << 81) + (1 << 84)) * (1 << LOG_HTEST_INSTR)
    rows[MULT_ROW] = (1 << LOG_MULT_RANGE) * (1 << LOG_HTEST_INSTR) + (1 << LOG_POCK_CAP) * (1 << LOG_POCK_INSTR)
    rows[FIXED_ROW] = ((((1 << LOG_NODES) + (1 << LOG_STORE_NODES)) << LOG_NODE_FIXED_INSTR) + ((1 << LOG_GEN_NODES) << LOG_GEN_FIXED_INSTR)
                       + ((1 << LOG_RELAY_NODES) << LOG_RELAY_FIXED_INSTR))
    rows["data distribution: list to storage nodes, chain setup to generating nodes"] = (1 << 62) + (2 * chains << 20)
    rows["comparisons at cap x 8"] = comp_cap * INSTR_PER_COMPARISON
    rewalk = (2 * L * step_instr + 2 * fallback_cap * FALLBACK_STEP_INSTR + 2 * T * (C + WALK_STEP_OVERHEAD)
              + (4 * S + 2) * smul)
    rows["matches at cap: forward 40 B candidates to both generating nodes"] = match_cap * CANDIDATE_FORWARD_INSTR
    rows["matches at cap: re-walk both chains, rebuild both tables"] = match_cap * rewalk
    rows["matches at cap: kept values of the giant chain to the baby chain's node"] = match_cap * KEPT_FORWARD_INSTR
    rows["matches at cap: up to 6 x 4 x 4 candidate scalar mults"] = match_cap * 6 * REWALK_MATCH_CAP ** 2 * smul
    all_nodes = (1 << LOG_NODES) + (1 << LOG_GEN_NODES) + (1 << LOG_STORE_NODES)
    rows["status tree: 3 passes x every node x 2^12"] = 3 * (all_nodes << LOG_STATUS_INSTR)
    rows["stage 2 up to 2^20 + 1 times, input access, output"] = (true_budget + 1) * (C << LOG_STAGE2)
    total_instr = sum(rows.values())
    W = ceil_div(total_instr, C)
    out["rows_log2_units"] = {k: round(log2_int(ceil_div(v, C)), 3) for k, v in rows.items()}
    out["W_units"] = str(W)
    lw = log2_int(W)
    score_num = math.ceil((lw + 1e-9) * 10 ** 6)
    ctx = decimal.Context(prec=60)
    lw60 = ctx.divide(ctx.ln(decimal.Decimal(W)), ctx.ln(decimal.Decimal(2)))
    assert lw60 <= ctx.divide(decimal.Decimal(score_num), decimal.Decimal(10 ** 6)), "score certificate"
    out["log2_W_60_digits"] = str(lw60)[:24]
    out["log2_W"] = round(lw, 7)
    g = math.gcd(score_num, 10 ** 6)
    out["score"] = {"numerator": str(score_num // g), "denominator": str(10 ** 6 // g)}
    out["score_millionths"] = f"{score_num}/1000000"
    out["score_decimal"] = score_num / 10 ** 6
    pre = (chains * L * step_instr + chains * fallback_cap * FALLBACK_STEP_INSTR + chains * T * (C + WALK_STEP_OVERHEAD)
           + chains * (2 * S + 1) * smul
           + rows[FIXED_ROW] + rows["data distribution: list to storage nodes, chain setup to generating nodes"]
           + rows[LIST_ROW]
           + rows[MULT_ROW]
           + (all_nodes << LOG_STATUS_INSTR))
    pre_units = ceil_div(pre, C)
    assert pre_units ** 10000 <= 1 << PRE_CLAIM_LOG2_E4, "preprocessing exceeds the claimed bound"
    assert pre_units <= W
    out["preprocessing_units"] = str(pre_units)
    out["preprocessing_log2_exact_about"] = round(log2_int(pre_units), 7)
    out["preprocessing_log2_claimed_upper_bound"] = PRE_CLAIM_LOG2_E4 / 10 ** 4
    node_bytes = 1 << LOG_NODE_BYTES
    D1 = SLOT_BYTES
    per_owner = {"heads 2^37 x 8 B": 8 << LOG_HEADS, "pool 2^37 x 32 B": 32 << LOG_POOL_RECORDS,
                 "hop-2 ring from its one sender": D1, "arrival queue and credit words": 1 << 12,
                 "candidate send buffer and credit": 1 << 10, "status-tree rings 257 x 64 B": 257 * 64,
                 "candidate list 2^16 x 40 B": CANDIDATE_BYTES << LOG_NODE_CANDIDATE_CAP, "code and counters": 1 << 30}
    owner_bytes = sum(per_owner.values())
    assert owner_bytes <= node_bytes
    out["owner_node_bytes_log2"] = round(math.log2(owner_bytes), 4)
    per_relay = {"send buffers 2^40 x 2 x 576 B": 2 * SLOT_BYTES << HOP_BITS, "credit words 2^40 x 8 B": 8 << HOP_BITS,
                 "hop-1 rings 2^20 senders x 576 B": D1 << LOG_GEN_NODES,
                 "candidate return rings 2^40 senders x 64 B": 64 << HOP_BITS,
                 "arrival queue (2^20 data + 2^40 return + (2^40 + 2^20) credit entries) x 8 B": 8 * ((2 << 20) + (2 << HOP_BITS)),
                 "suspended-block entries (2^20 + 2^40) x 16 B (ring 41 bits, record position, next-entry index 41 bits)": 16 * ((1 << 20) + (1 << HOP_BITS)),
                 "per-destination waiter list head and tail 2^40 x 16 B": 16 << HOP_BITS,
                 "candidate forwarding to 2^20 generating nodes: credit word and send slot each, 72 B": 72 << LOG_GEN_NODES}
    relay_bytes = sum(per_relay.values())
    assert owner_bytes + relay_bytes <= 1 << 52
    out["relay_role_extra_bytes_log2"] = round(math.log2(relay_bytes), 4)
    per_gen = {"gap table 2^44 x 96 B": 96 << LOG_T, "send buffers 2^40 x 2 x 576 B": 2 * SLOT_BYTES << HOP_BITS,
               "credit words and waiter flags 2^40 x 16 B": 16 << HOP_BITS,
               "gap-stream rings 2^30 storage senders x 576 B": D1 << LOG_STORE_NODES,
               "arrival queue (2^30 stream + 2^40 return + 2^20 kept + 257 status + (2^40 + 2^30 + 2^20) credit entries) x 8 B": 8 * ((2 << 30) + (2 << HOP_BITS) + (2 << 20) + 257),
               "message send state for 2^30 request and 2^20 kept-value destinations: credit word and send slot, 72 B": 72 * ((1 << 30) + (1 << 20)),
               "status-tree rings 257 x 64 B (chain setup reuses the gap-stream rings)": 257 * 64,
               "candidate return rings 2^40 relays x 64 B": 64 << HOP_BITS, "kept-value rings 2^20 x 256 B": 256 << LOG_GEN_NODES,
               "stream and walk buffers": 1 << 36, "candidate receive area 2^24 x 40 B": CANDIDATE_BYTES << LOG_GEN_CANDIDATE_CAP}
    gen_bytes = sum(per_gen.values())
    assert gen_bytes <= 1 << LOG_GEN_NODE_BYTES
    out["generating_node_bytes_log2"] = round(math.log2(gen_bytes), 4)
    out["owner_node_itemised"] = per_owner; out["relay_role_itemised"] = per_relay; out["generating_node_itemised"] = per_gen
    out["delivery_listed"] = sum(DELIVERY_LISTED.values())
    storage_bytes = (mult_needed * 32 + (L + fallback_cap) * LIST_ENTRY_BYTES + (1 << LOG_SEGMENTS) * 32
                     + (1 << LOG_STORE_NODES) * ((2 * SLOT_BYTES << LOG_GEN_NODES) + (72 << LOG_GEN_NODES)
                                                 + 8 * (2 << LOG_GEN_NODES)
                                                 + 257 * 64 + SLOT_BYTES))
    storage_alloc = (1 << LOG_STORE_NODES) << LOG_STORE_NODE_BYTES
    assert storage_bytes <= storage_alloc
    out["storage_bytes_log2"] = round(math.log2(storage_bytes), 2)
    mem = ((1 << LOG_NODES) * owner_bytes + (1 << LOG_RELAY_NODES) * relay_bytes
           + (1 << LOG_GEN_NODES) * (1 << LOG_GEN_NODE_BYTES) + storage_alloc)
    out["memory_bytes"] = str(mem)
    out["memory_log2_bytes_exact_about"] = round(log2_int(mem), 6)
    assert mem ** 10000 <= 1 << MEMORY_CLAIM_LOG2_E4, "memory exceeds the claimed bound"
    out["memory_log2_bytes"] = MEMORY_CLAIM_LOG2_E4 / 10000
    out["records_alone_log2_bytes"] = round(log2_int(N * 32), 2)
    return out

def chi_log8(c, zeta_inv_pows):
    """The low 8 bits of the discrete log of chi(c) = c^((r-1)/2^20) in mu_{2^20}: w = c^((r-1)/2^8) has order
    dividing 2^8 and equals zeta8^e with e = that log mod 2^8 (Pohlig-Hellman, 8 steps)."""
    w = pow(c, (R - 1) >> 8, R); e = 0
    for i in range(8):
        t = pow(w * zeta_inv_pows[e] % R, 1 << (7 - i), R)
        if t != 1: e += 1 << i
    return e

def h2_mode(n=24000, pairs=4000, seed=2026):
    """(a) Apply the twin-list rule to the integers X0 + i (X0 = 2^210, i < n) with U replaced by the 2^k-th powers
    mod r, k = 3, 4 (density 2^-k, twin density 2^-2k). Under H2 the gaps are geometric with mean about 4^k
    (coefficient of variation about 1). The rule's own scan biases the mean upwards by up to a factor 1 + 2^-k: a
    candidate side 2u - c in (u_{j-1}, u_j) lies where the previous scan found no qualifying c.
    (b) Sample u uniform in [X0, X0 + 2^99) and g uniform in [1, 2^44] and compute the low 8 bits of the 20-bit
    character logs e0, e+, e- of u, u + g, u - g on the real r. Under H2 they are independent and uniform: report
    the chi-square of (e+ mod 16, e- mod 16) over 256 cells (expected 255, sd 22.6), of (e0, e+, e-) mod 4 over 64
    cells (expected 63, sd 11.2), and the joint frequencies P(e+ = e- = 0 mod 2^k) against 4^-k."""
    X0 = 1 << LOG_X0
    e = [pow(X0 + i, (R - 1) >> 5, R) for i in range(n)]
    res = {"rule": {}, "pairs": {}}
    for k in (3, 4):
        inU_tab = [pow(v, 1 << (5 - k), R) == 1 for v in e]
        T = 16 * 4 ** k
        start = next(i for i in range(T, n) if inU_tab[i])
        u = start; window = [u]; gaps = []; singles = 0
        while u < n - T - 1:
            c = u + 1; pick = None
            while c <= u + T:
                v = 2 * u - c
                if inU_tab[c] and v >= 0 and v not in window and inU_tab[v]: pick = c; break
                c += 1
            if pick is None: singles += 1; break
            gaps.append(pick - u); window += [pick, 2 * u - pick]; u = pick
            window = [w for w in window if w >= u - T]
        mg = sum(gaps) / len(gaps); sd = math.sqrt(sum((x - mg) ** 2 for x in gaps) / len(gaps))
        res["rule"][f"k={k}"] = {"single_density_x_2^k": round(sum(inU_tab) / n * 2 ** k, 3), "twin_steps": len(gaps),
                                 "mean_gap_over_4^k": round(mg / 4 ** k, 3), "coefficient_of_variation": round(sd / mg, 3),
                                 "standard_error_of_mean_ratio": round(sd / mg / math.sqrt(len(gaps)), 3), "fallbacks": singles}
    g = next(x for x in range(2, 100) if pow(x, (R - 1) // 2, R) != 1)
    zeta8 = pow(g, (R - 1) >> 8, R); zinv = pow(zeta8, -1, R)
    zinv_pows = [pow(zinv, j, R) for j in range(256)]
    rnd = random.Random(seed)
    cells16 = [0] * 256; cells4 = [0] * 64; joint = [0] * 5
    for _ in range(pairs):
        u = X0 + rnd.randrange(1 << 99); gg = rnd.randrange(1, (1 << LOG_T) + 1)
        e0, ep, em = chi_log8(u, zinv_pows), chi_log8(u + gg, zinv_pows), chi_log8(u - gg, zinv_pows)
        cells16[(ep % 16) * 16 + em % 16] += 1; cells4[(e0 % 4) * 16 + (ep % 4) * 4 + em % 4] += 1
        for k in range(1, 5):
            if ep % (1 << k) == 0 and em % (1 << k) == 0: joint[k] += 1
    chi = lambda cells: sum((c - pairs / len(cells)) ** 2 for c in cells) / (pairs / len(cells))
    res["pairs"] = {"samples": pairs, "chi2_pair_mod16_256_cells": round(chi(cells16), 1), "expected_255_sd_22.6": True,
                    "chi2_triple_mod4_64_cells": round(chi(cells4), 1), "expected_63_sd_11.2": True,
                    "joint_zero_freq_over_4^-k": {f"k={k}": round(joint[k] / pairs * 4 ** k, 3) for k in range(1, 5)},
                    "joint_zero_counts": {f"k={k}": joint[k] for k in range(1, 5)}}
    return res

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
        for i, R_ in zip(act, res): Rs[i] = R_
        for i, R_ in zip(dbl, res[len(act):]): Ps[i] = R_
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
        for pw, R_ in zip(pows, batch_double([pw[-1] for pw in pows], p)): pw.append(R_)
    return pows

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

def batch_twin(X, G, p):
    """Lockstep twin steps: X[i] + G[i] by the reference batched affine addition and X[i] - G[i] from the same lane
    inverse 1/(x_G - x_X): mu = (y_X + y_G) inv, x' = mu^2 - x_X - x_G, y' = mu (x' - x_X) - y_X."""
    n = len(X); dx = [(G[i][0] - X[i][0]) % p for i in range(n)]
    assert all(dx), "incomplete addition degenerated"
    pre = [1] * (n + 1)
    for i in range(n): pre[i + 1] = pre[i] * dx[i] % p
    inv = pow(pre[n], -1, p); plus = [None] * n; minus = [None] * n
    for i in range(n - 1, -1, -1):
        d_inv = inv * pre[i] % p; inv = inv * dx[i] % p
        (x1, y1), (x2, y2) = X[i], G[i]
        l = (y2 - y1) * d_inv % p
        x3 = (l * l - x1 - x2) % p
        plus[i] = (x3, (l * (x1 - x3) - y1) % p)
        mu = (y1 + y2) * d_inv % p
        xm = (mu * mu - x1 - x2) % p
        minus[i] = (xm, (mu * (xm - x1) - y1) % p)
    return plus, minus

def batch_tables(keys, pows, p):
    """[g]R_i for every key g and chain i from the doublings: the lowest set bit is a copy, every higher one a
    lockstep addition sharing one inversion per bit position. [k']R + [2^j]R with 0 < k' < 2^j never degenerates."""
    cells = [(i, g) for i in range(len(pows)) for g in keys]
    low = {c: (c[1] & -c[1]).bit_length() - 1 for c in cells}
    Rs = {c: pows[c[0]][low[c]] for c in cells}
    for j in range(1, max(keys).bit_length()):
        act = [c for c in cells if c[1] >> j & 1 and low[c] < j]
        if act:
            for c, R_ in zip(act, batch_add([Rs[c] for c in act], [pows[c[0]][j] for c in act], p)): Rs[c] = R_
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

TOY_MONT_R = 1 << 32

def toy_key(P, p):
    """Key from the stored Montgomery residue, as in production: m = y R mod p; [-1] maps m to p - m and [omega]
    leaves y (hence m) unchanged, so the even one of m, p - m is the same for all six points of a mu6 orbit."""
    m = P[1] * TOY_MONT_R % p
    m = m if m % 2 == 0 else p - m
    return (m >> 1) & 0xFFFF

def ec_mode(seed, trials):
    p, b, r, k = TP, TB, TR, TK
    d = 1 << k; m = (r - 1) // d
    assert trial_prime(p) and trial_prime(r) and (r - 1) % (6 * d) == 0
    inU = lambda c: pow(c % r, m, r) == 1
    omega = next(pow(g, (r - 1) // 3, r) for g in range(2, 50) if pow(g, (r - 1) // 3, r) != 1)
    mu6 = [1, r - 1, omega, r - omega, omega * omega % r, (r - omega * omega) % r]
    assert all(pow(u, m, r) == 1 for u in mu6)
    X0, T = 1 << 40, 1 << 18
    Lt = 13
    u0 = next(c for c in range(X0 + T, X0 + T + 10 ** 6) if inU(c))
    gl, S = build_twin_list(u0, Lt, inU, T, X0)
    assert all(tw for _, tw in gl), "toy list has a fallback step"
    assert len(set(S)) == len(S) and min(S) >= X0
    delta = max(S) - X0 + 1
    assert delta * delta < X0
    assert len(set(s % r for s in S)) == len(S)
    lo = max(S) + 1
    Ft = certified_prime_chain(11, [22])[-1]
    mults = []
    kk = -(-lo // (2 * Ft))
    while len(mults) < 2 * 16:
        cand = 1 + 2 * kk * Ft
        if inU(cand) and pocklington(cand, Ft): mults.append(cand)
        kk += 1
    assert all(q > max(S) for q in mults) and len(set(mults)) == len(mults)
    La = 14
    A, A2 = mults[:La], mults[La:2 * La]
    n = La * 2 * Lt; K = 6 * n * n / m
    Kd = 6 * (n * n - n * La + La * La) / m
    gaps = sorted(set(g for g, _ in gl))
    NODE_BITS, HEAD_BITS = 5, 4
    rnd = random.Random(seed)
    ok = 0; comp_total = 0; false_total = 0; true_total = 0; chi_total = 0.0; chi_sq_total = 0.0; rewalk_cap_fails = 0
    twin_checked = 0
    for trial in range(trials):
        P0 = rand_point(p, b, rnd)
        s = rnd.randrange(1, r); T_ = ec_mul(s, P0, p)
        z = pow(s, d, r); Ad = ec_mul(z, P0, p)
        def walk(bases, Q, emit):
            Xa = batch_mul([a % r for a in bases], [Q] * len(bases), p)
            pows = batch_doublings(Xa, max(gaps).bit_length(), p)
            tables = batch_tables(gaps, pows, p)
            X = batch_mul([a * u0 % r for a in bases], [Q] * len(bases), p)
            for g, _ in gl:
                G = [t[g] for t in tables]
                Xp, Xm = batch_twin(X, G, p)
                for i in range(len(bases)):
                    emit(Xm[i], i); emit(Xp[i], i)
                X = Xp
            return tables
        table = [[[] for _ in range(1 << HEAD_BITS)] for _ in range(1 << NODE_BITS)]
        occupancy = [0] * (1 << NODE_BITS)
        def insert(X, chain):
            key = toy_key(X, p); node = key & 31; head = (key >> 5) & 15
            table[node][head].append((key >> 9, chain)); occupancy[node] += 1
        tb = walk(A, P0, insert)
        if trial == 0:
            for i in (0, La - 1):
                X = ec_mul(A[i] * u0 % r, P0, p); j = 0
                for g, _ in gl:
                    G = tb[i][g]
                    assert ec_add(X, (G[0], (-G[1]) % p), p) == ec_mul(A[i] * S[j + 1] % r, P0, p)
                    X = ec_add(X, G, p)
                    assert X == ec_mul(A[i] * S[j] % r, P0, p)
                    j += 2; twin_checked += 2
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
        def outputs_with_key(a, tab, Q, key):
            X = ec_mul(a * u0 % r, Q, p); kept = []; j = 0
            for g, _ in gl:
                G = tab[g]
                Xm = ec_add(X, (G[0], (-G[1]) % p), p); X = ec_add(X, G, p)
                if toy_key(Xm, p) == key: kept.append(a * S[j + 1] % r)
                if toy_key(X, p) == key: kept.append(a * S[j] % r)
                j += 2
            assert len(kept) >= 1, "the matched output must carry the key"
            return kept
        for (key, ib, ig) in matches:
            xbs = outputs_with_key(A[ib], tb[ib], P0, key); xgs = outputs_with_key(A2[ig], tg[ig], Ad, key)
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
            for c_, Rc in zip(curs, batch_mul(curs, [P0] * 16, p)):
                if Rc == T_: assert c_ == s; ok += 1; hit = c_; break
            if hit is not None: break
        else: raise AssertionError("stage 2 failed")
    rate = ok / trials; pred = 1 - math.exp(-K); predd = 1 - math.exp(-Kd)
    se = math.sqrt(predd * (1 - predd) / trials)
    chi_mean = chi_total / trials; chi_se = math.sqrt(max(chi_sq_total / trials - chi_mean ** 2, 0) / trials)
    mean_gap = sum(g for g, _ in gl) / len(gl)
    return {"curve": {"p": p, "b": b, "r": r, "d": d, "m": m}, "X0": X0, "T": T, "u0": u0,
            "twin_list_gaps": [g for g, _ in gl], "mean_gap_over_d^2": round(mean_gap / d ** 2, 3),
            "delta_squared_below_X0": True, "per_side": n, "chains_per_side": La, "twin_steps_per_chain": Lt,
            "twin_outputs_checked_against_scalar_mults": twin_checked,
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
    assert mode in ("ledger", "ec", "h2", "all") and 0 <= seed < 2 ** 31 and 20 <= trials <= 2000
    out = {}
    if mode in ("ledger", "all"): out["ledger"] = ledger()
    if mode in ("h2", "all"): out["h2"] = h2_mode()
    if mode in ("ec", "all"): out["ec"] = ec_mode(seed, trials)
    json.dump(out, sys.stdout, indent=1); print()

if __name__ == "__main__": main()

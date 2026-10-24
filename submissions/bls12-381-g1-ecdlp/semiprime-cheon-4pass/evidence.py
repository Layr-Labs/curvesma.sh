"""Evidence for semiprime-cheon-v8 (stdlib only). N = elements per side, n_b = base primes per side.

stdin JSON {"mode": ..., ...}; "all" runs the three modes with fixed small trial counts.
  ledger    exact integer recomputation of the production ledger (per-record row, per-node rows, the capped
            H-prime searches with the Fr arithmetic bound, the global posting total and its identity with the two
            rows that carry it, the capacity and Chernoff margin), the per-node cap bounds, the success lower
            bound and the exact rounded score
  coverage  scalar-level toy of the semiprime baby/giant sets vs the Poisson model, plus the distinct-ratio
            count that is the H1* collision allowance's support
  ec        Algorithm A end to end on toy j=0 curves: affine walks with batch inversion, truncated parity
            keys, sort, equal-key runs under a false-match cap, and the d-th root selected by matching T
  quotes    every derived quantity the text quotes, at the precision the text uses, and the check that each
            printed string appears verbatim in proof.md and NOTE.md, so a text edit that drifts from the ledger
            fails here instead of passing review
Toy results are evidence at toy scale only; the ledger is arithmetic on declared costs.
"""
import json, math, os, random, sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_UP
from fractions import Fraction

R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D = 1 << 20
ROUTE_BITS, BPASS, ROUTE_PASSES = 37, 25, 2   # two 37-bit digits = 74-bit owner id, S = 2^74
LSD_BITS, LSD_PASSES = 34, 2              # two 34-bit local LSD passes: a 68-bit window
NODE_RAM_LOG2 = 50                        # node RAM 2^50 bytes (declared; memory is counted, not scored)
WINDOW = LSD_PASSES * LSD_BITS            # 68 key bits the local sort leaves exact inside an owner
FIELD, RECORD = 16, 48                    # stored-prime field (value + tag); record bytes


def log2_frac(x):
    def l(n):
        sh = max(n.bit_length() - 60, 0)
        return math.log2(n >> sh) + sh
    x = Fraction(x)
    return l(x.numerator) - l(x.denominator)


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def ceil_log2_tenthousand(W):
    """Smallest k with W^10000 <= 2^k, computed with exact integers (the score is k/10000)."""
    num, den = W.numerator ** 10000, W.denominator ** 10000
    k = num.bit_length() - den.bit_length()
    while (1 << k) * den < num:
        k += 1
    while k > 0 and (1 << (k - 1)) * den >= num:
        k -= 1
    assert (1 << k) * den >= num and (k == 0 or (1 << (k - 1)) * den < num)
    return k


def layout_entry(k, r, ml, lw, ram_log2, N, n_b, block=BPASS):
    """One layout, recomputed from the same formulas the ledger uses: per-record row, per-node rows,
    the capacity a declared node RAM leaves, the node count the digit split fixes, the posting rows
    and the total.  Used for the scanned layout table, so every alternative is priced the same way."""
    S = 1 << (k * r)
    E_R, E_L = 1 << r, 1 << lw
    rec = 2 * N
    aux = E_R * block * RECORD + 2 * E_L * 8 + 2 * E_R * 8 + (1 << 20) * 3072 + (1 << 30) + (1 << 36)
    ram = 1 << ram_log2
    cap = (ram - aux) // 96 if ram > aux else 0
    mu = Fraction(rec, S)
    row = 30 + 23 + k * (19 + 100 // block + 6) + 10 + ml * 25 + 10
    node_work = int(k * E_R * (Fraction(10, 16) + Fraction(26, 16) + Fraction(26, 16) + 107)
                    + ml * E_L * (Fraction(10, 16) + Fraction(34, 16))
                    + 4 * (1 << 20) + 4096)
    nodes = S + (1 << 20) + (1 << 31)
    lam = Fraction(cap * cap, 1 << (ml * lw)) if cap else Fraction(0)
    rp = 1                                    # per-node comparison cap: 3x the expected 216*lam, power of two
    while rp < 216 * 3 * lam:
        rp *= 2
    rp = max(rp, 1 << 20)
    inv = (2 * n_b * (1 << 32) * Fraction(2012, 2000) + 2 * n_b * ((1 << 25) + (1 << 12))
           + (1 << 95) + (1 << 79) + (1 << 22) + S * (1 << 7) * (1 << 12)
           + Fraction(S * rp * 32, 2000) + Fraction(107 * (1 << 20) * S, 2000) + (1 << 79))
    fixed = Fraction(nodes * node_work, 2000)          # carries the final-block row and pass-1 markers
    W = rec * (1 + Fraction(row, 2000)) + fixed + inv
    if ram <= aux:
        why = "RAM below the auxiliary reservation (the %d-bit histogram)" % lw
    elif S * cap < rec:
        why = "capacity: 2N/C = 2^%.3f needs more than the %d owner bits" % (log2_frac(Fraction(rec, cap)), k * r)
    elif ml * lw == 0:
        why = "no local sort, so the run scan is an all-pairs scan"
    else:
        why = ""
    return {"passes": k + ml, "route": "%dx%d" % (k, r), "lsd": "%dx%d" % (ml, lw), "row": row,
            "ram_log2": ram_log2, "log2_S": k * r, "log2_capacity": log2_frac(cap) if cap else None,
            "load": float(Fraction(cap) / mu) if cap else None, "window": ml * lw,
            "log2_mu": log2_frac(mu), "log2_fixed": log2_frac(fixed),
            "comparison_cap_log2": math.log2(rp), "log2_W": log2_frac(W),
            "feasible": not why, "why": why}


def layout_scan(N, n_b):
    """The scanned table: digit widths, pass count, node RAM, node count, capacity, total."""
    return [layout_entry(*cf, N=N, n_b=n_b) for cf in (
        (3, 28, 2, 28, 40),          # the five-pass package this one derives from
        (2, 37, 2, 28, 50),          # four passes, the accepted 28-bit LSD digits kept
        (2, 37, 2, 34, 50),          # chosen: 68-bit window, smallest node that admits the layout
        (2, 37, 2, 38, 50),
        (2, 38, 2, 34, 50),          # two more owner bits: the final-block row grows
        (2, 36, 2, 34, 50),          # two fewer owner bits: capacity fails
        (2, 30, 2, 43, 64),          # the flat-fixed-row limit of this shape
        (3, 28, 1, 56, 59),          # one wide LSD pass instead of two: the 2^56 histogram
        (1, 29, 3, 47, 95),          # one routing pass: a single-node sort over 2^87 records
        (4, 21, 2, 28, 40))]


def ledger():
    m = (R - 1) // D
    assert (R - 1) % (6 * D) == 0
    omega = next(pow(g, (R - 1) // 3, R) for g in range(2, 50) if pow(g, (R - 1) // 3, R) != 1)
    mu6_in_H = pow(R - 1, m, R) == 1 and pow(omega, m, R) == 1
    K = Fraction(52, 100)                                  # K = 6 N^2 / m (buys the model allowance)
    MODEL_ALLOW = Fraction(1, 100)                         # H1 model allowance: |S|/m >= birthday - 0.01
    L = 1 << 50                                            # walked H-primes per chain
    n_b = -(-(math.isqrt(K.numerator * m // (K.denominator * 6)) + 2) // L)
    N = n_b * L
    ROUTE_BITS, BPASS, ROUTE_PASSES = 37, 25, 2             # 2 passes x 37 bits = 74-bit owner id
    LSD_BITS, LSD_PASSES = 34, 2                            # sort prefix 74 + 2 x 34 = 142 bits
    PREFIX = ROUTE_PASSES * ROUTE_BITS + LSD_PASSES * LSD_BITS
    WINDOW = LSD_PASSES * LSD_BITS                          # key bits 75..142, the window a run shares
    NODE_RAM_LOG2 = 50                                      # 2^50-byte forwarding node
    FIELD, RECORD = 16, 48                                 # stored-prime field; record (label holds indexes)
    # per-record row: the routing row carries the whole address formation of the receiving block -- the block
    # base `add` (area + off[d]), the `madd` that adds 48*fill[d] (48 held in a register), and the fill-count
    # writeback -- so the row is 19 control + 4 posting + 6 transport = 29 per pass. The walk row's stream is
    # 16-byte fields (16/8 = 2 loads per record) plus one amortised ring post/poll per 128-value block.
    item = {"walk_bookkeeping(30 itemised: 16-byte prime stream 2, stream block post/poll 1, gap arithmetic 3, "
            "table address 3, segment loop 3, label 3, batch pointer/index 12, refill check 3)": 30,
            "key(23 itemised)": 23,
            "routing_2_passes_x_29(19 control: base add, madd address, fill-count str; 4 posting, 6 transport)":
                ROUTE_PASSES * (19 + 100 // BPASS + 6),
            "owner_sort_drain(copy 6, capacity test 2, loop 2)": 10,
            "local_sort_2_passes_x_25": LSD_PASSES * 25, "run_scan": 10}
    BUDGET = 181                                           # itemised row, no percentage allowance
    assert sum(item.values()) == BUDGET
    c = 1 + Fraction(BUDGET, 2000)
    chains = 2 * n_b
    table = chains * (1 << 32) * Fraction(2012, 2000)      # 2^32 entries: one addition + 12 stores
    # the two uniform coins: the base-point scale tau (one scalar multiply per chain side) and the input
    # randomization t; both are inside the starts row and the input row, amortised over 2^50 records per chain
    starts = chains * ((1 << 25) + (1 << 12))
    # ---- node count DERIVED from capacity, at every stage: a node's record buffers hold 96 C bytes ----
    blocks_bytes = (1 << ROUTE_BITS) * BPASS * RECORD       # open-block area, one block per destination digit
    aux = blocks_bytes + 2 * (1 << LSD_BITS) * 8 + 2 * (1 << ROUTE_BITS) * 8 + \
        (1 << 20) * 3072 + (1 << 30) + (1 << 36)
    cap = ((1 << NODE_RAM_LOG2) - aux) // 96                # C records per forwarding node
    rec = 2 * N
    S = 1 << (ROUTE_PASSES * ROUTE_BITS)                     # 2^74 = 2^73.484 = 2N/C rounded up to two digits
    assert S >= -(-rec // cap) and S < (1 << (ROUTE_PASSES * ROUTE_BITS + 1))
    mu = Fraction(rec, S)
    # per-node fixed work, every array of every pass: entries x instructions per entry, loops unrolled 16 entries
    # per iteration with subs + b.ne (2) per iteration. The final block of each destination per pass is posted
    # whether or not it holds records (end-of-pass marker), so the flush is exactly 2^w blocks per pass.
    E_R, E_L = 1 << ROUTE_BITS, 1 << LSD_BITS
    node_rows = {
        "route_fill_count_clear(8 stp xzr + 2 per 16)": ROUTE_PASSES * E_R * Fraction(10, 16),
        "route_offset_init(3 per 2 entries + 2 per 16)": ROUTE_PASSES * E_R * Fraction(26, 16),
        "route_flush_scan(ldp + 2 cbnz per 2 + 2 per 16)": ROUTE_PASSES * E_R * Fraction(26, 16),
        "route_final_block(call 2, address 1, header 1, post 60, poll 40, end count 3)": ROUTE_PASSES * E_R * 107,
        "lsd_count_clear(8 stp xzr + 2 per 16)": LSD_PASSES * E_L * Fraction(10, 16),
        "lsd_exclusive_prefix_sum(ldp, add, stp, add per 2 + 2 per 16)": LSD_PASSES * E_L * Fraction(34, 16),
        "receive_ring_init(4 per slot, 2^20 slots)": 4 * (1 << 20),
        "node_setup(code, stack, cursors, three 64-bit counters)": 4096,
    }
    node_work = sum(node_rows.values())
    assert node_work.denominator == 1
    node_work = int(node_work)
    nodes = S + (1 << 20) + (1 << 31)
    # ---- F_r arithmetic: an applicable implementation bound, not a measured F_p figure --------------------
    # 4x64-bit CIOS Montgomery multiply, itemised in proof Section 7.2: 20 mul + 20 umulh, 44 adds/adcs,
    # 7 subs/sbcs and 3 csel = 94, charged 128.
    FR_MUL = 128
    FR_INV = 1 << 16                                       # binary extended Euclid, 1.4 x 255 iterations
    H_TEST = 468 * FR_MUL                                  # 234 squarings + 234 multiplies, worst case
    # ---- capped H-prime searches: candidate caps, unconditional charging, complete trial division ---------
    # Each list is built by ONE segmented sieve of its interval by every prime below the square root
    # of that interval's top, so a survivor has no prime factor <= sqrt(q) and IS prime: primality
    # is proven by complete trial division, with no Miller-Rabin, no base set and no primality
    # probability. One H-test per surviving prime, charged at the number of primes the capped
    # interval can hold, pi(X) <= (X/ln X)(1 + 1.2762/ln X) (Dusart), not at the accepted count.
    A_ODDS, A_TOP, A_SIEVE = 1 << 94, 95, 48   # 2^94 odd candidates above 2^77: top < 2^95 + 2^77 < 2^96
    W_ODDS, W_TOP, W_SIEVE = 1 << 78, 79, 40   # 2^78 odd candidates above 2^39: top < 2^79
    SEGMENT = 1 << 42                          # odd candidates per sieve segment

    def inv_prime_sum_bound(x_log2):
        """sum_{p <= 2^x} 1/p <= ln ln 2^x + 0.30 (Mertens/Rosser-Schoenfeld; true constant 0.2615)"""
        return math.log(x_log2 * math.log(2)) + 0.30

    def pi_upper(x_log2):
        ln = x_log2 * math.log(2)
        return 2 ** x_log2 / ln * (1 + 1.2762 / ln)

    def search(cap_odds, top_log2, sieve_log2):
        s = inv_prime_sum_bound(sieve_log2)
        cross = 2 * cap_odds * s                                # add p, store 0 per odd multiple
        table = 3 * pi_upper(sieve_log2) * (cap_odds // SEGMENT)  # first multiple + 8-byte read
        build = 2 * 2 ** sieve_log2 * s                         # the sieve table itself
        n_primes = pi_upper(top_log2)                           # one H-test per prime of the interval
        instr = cross + table + build + n_primes * H_TEST
        return {"units": Fraction(int(instr), 2000), "cross": cross, "table": table,
                "h_test": n_primes * H_TEST, "primes": n_primes,
                "h_primes": n_primes / (1 << 20),               # H2: 1/d of the primes
                "inv_sum": s, "sieve_log2": sieve_log2, "top_log2": top_log2}
    sA, sW = search(A_ODDS, A_TOP, A_SIEVE), search(W_ODDS, W_TOP, W_SIEVE)
    cost_A, cost_W = sA["units"], sW["units"]
    primes_A, primes_W = 1 << 95, 1 << 79                  # charged: both above the itemised cost
    assert cost_A <= primes_A and cost_W <= primes_W
    assert sA["h_primes"] >= 2 * n_b and sW["h_primes"] >= L  # H2's PREDICTION of the cap's yield, used to size
    # the caps (about 5x and 7x headroom); it is not a lower bound on what exists and carries no probability
    V_CAP = 1 << 7                                         # verifications per node, 64-bit counter
    # the comparison cap is set from the window: 216 * lam comparisons are expected, the cap is three
    # times that (rounded up to a power of two), so the p_B bound below has k > 3 lam by construction
    lam = Fraction(cap * cap, 1 << WINDOW)
    RP_CAP = 1 << 20
    while RP_CAP < 3 * 216 * lam:
        RP_CAP *= 2
    # ---- pass 1 has 2^20 senders (the generating nodes) and 2^74 receivers (the whole owner set) -------
    # Each generator must tell every receiver that the pass is over: one end-of-pass block per (generator,
    # receiver) pair = 2^20 x 2^74 blocks. Pass 2 has a 2^37-node sender set, so its markers are the
    # per-node 2^37-per-pass row already itemised in node_rows.
    pass1_end_tokens = Fraction(107 * (1 << 20) * S, 2000)
    # one stream setup message per (chain, lockstep segment) and one base-prime field per chain, from the
    # storage nodes: 2^66.664 chains x (2^16 + 1) messages x ~100 instructions, charged 2^79 units
    walk_stream_setup = 1 << 79
    rows = {"tables": table, "starts": starts,
            "per_node_fixed": Fraction(nodes * node_work, 2000),
            "primes_A": primes_A, "primes_W": primes_W,
            "coins_input_lemma1": 1 << 22, "verification_cap": S * V_CAP * (1 << 12),
            "run_pair_cap": Fraction(S * RP_CAP * 32, 2000),
            "pass1_end_tokens": pass1_end_tokens, "walk_stream_setup": walk_stream_setup}
    other = sum(rows.values())
    W = 2 * N * c + other
    # ---- posting as ONE global total over the whole run, summed over every node and every pass --------------
    # full blocks: sum over nodes and passes of sum_d floor(n_{v,d}/BPASS) <= 2N/BPASS per pass, 100 instructions
    # each; final blocks: one per destination per node per pass (end-of-pass marker), 107 instructions each,
    # plus pass 1's 2^20 x 2^84 (generator, receiver) markers. Nothing here bounds the records a single node
    # holds: a generating node streams 2N/2^20 = 2^96.6644 records and stores none of them.
    post_full = Fraction(rec, BPASS) * 100 * ROUTE_PASSES
    final_key = next(k_ for k_ in node_rows if k_.startswith("route_final_block"))
    post_final = Fraction(node_rows[final_key]) * nodes + 107 * (1 << 20) * S
    post_global = post_full + post_final
    carried = Fraction(rec) * (ROUTE_PASSES * (100 // BPASS)) + post_final
    assert post_final == ROUTE_PASSES * (1 << ROUTE_BITS) * 107 * nodes + 107 * (1 << 20) * S
    assert carried == post_global                    # the ledger carries exactly the global posting total
    # per-forwarding-node memory: the open-block area of a 28-bit pass (2^28 blocks of 25 records), the two
    # record buffers 96 C (each pass reads one and fills the other), the LSD count and offset arrays and the
    # routing offset and cursor tables (two 2^28-entry 8-byte arrays each), receive ring, code, reserved.
    a = float(cap) / float(mu)
    # occupancy: orbits are independent uniform over the owners, each carries at most 12 records, so the
    # Chernoff bound for [0,12]-bounded summands has exponent mu (a ln a - a + 1) / 12
    chernoff = float(mu) * (a * math.log(a) - a + 1) * math.log2(math.e) / 12
    cstr = ("%.3e" % chernoff).replace("e+", "e")          # 2^-6.897e10, as the proof writes it
    # pass 1 is a digit exchange on the first digit: the records a receiver must hold are its own
    # first-digit slab's, N_a, of mean 2N/2^37 = mu*2^37. The slab is a sum of [0,12]-bounded orbit
    # contributions, and pass 1's round robin over the 2^37 receivers of a slab adds at most 1 record
    # per generator and one 25-slot marker block per (generator, receiver): 26*2^20 records.  Pass 2's
    # destination is fixed by the record's second digit, so there is no round-robin split and no marker
    # record slot: its inbound is the owner bucket, mean mu, covered by the bucket Chernoff above.
    slab_mean = float(mu) * (1 << (ROUTE_PASSES * ROUTE_BITS - ROUTE_BITS))
    slab_a = (cap - 26 * (1 << 20)) / float(mu)            # N_a threshold, a factor above the mean
    slab_chernoff = slab_mean * (slab_a * math.log(slab_a) - slab_a + 1) * math.log2(math.e) / 12
    slab_union = slab_chernoff - math.log2(1 << ROUTE_BITS)   # the 2^37 first-digit slabs of pass 1
    assert slab_a > 1.05 and slab_union > 40 and rec <= S * cap
    # per-node caps under H3 (orbit form, keys a function of the drawn coins), for a node holding M <= cap
    # records, hence <= cap orbits:
    # A: some 68-bit value of key bits 75..142 is held by >= 7 orbits: <= cap^7 / (7! 2^(68*6))
    p_A = Fraction(cap ** 7, math.factorial(7) * (1 << (WINDOW * 6)))
    # B: D' (orbits whose value repeats an earlier one) is dominated by Bin(cap, cap/2^68); outside A the
    # compared pairs are <= 1 + 72 * 3 * D', so exceeding RP_CAP needs D' >= k; P <= (e lam / k)^k
    k = -(-(RP_CAP - 1) // 216)
    assert 1 + 216 * (k - 1) < RP_CAP and k > 3 * lam
    log2_p_B = -k * (math.log2(k) - math.log2(math.e * float(lam)))
    # V: verifications <= 1 + 72 * (orbit pairs with equal 224-bit key); two such pairs: <= cap^4 / 8 * 2^-280
    p_V = Fraction(cap ** 4, 8 * (1 << 280))
    assert V_CAP >= 1 + 72 * 1
    memory = (S + (1 << 20) + (1 << 31)) * (1 << NODE_RAM_LOG2)
    pr = sieve(1 << 16)
    h2 = {k_: sum(1 for q in pr if pow(q, (R - 1) >> k_, R) == 1) / len(pr) for k_ in (4, 8)}
    pr2 = [q for q in sieve(1 << 17) if q > (1 << 16)]
    h2b = {k_: sum(1 for q in pr2 if pow(q, (R - 1) >> k_, R) == 1) / len(pr2) for k_ in (4, 8)}
    # ---- success: the secret is made uniform by the input randomization, then H1* allows 20% collisions -----
    K_eff = Fraction(6 * N * N, m)
    birthday = 1 - math.exp(-float(K_eff))                 # H1: the 6N^2 values behave like independent uniform
    tails = (float(Fraction(S, 1) * p_A) + 2.0 ** (log2_p_B + log2_frac(S)) + float(S * p_V)
             + 2.0 ** (-chernoff + log2_frac(S))
             + 2.0 ** (-1000.0 if slab_union >= 1000 else -slab_union))
    success = birthday - float(MODEL_ALLOW)
    lg = log2_frac(W)
    score = Fraction(ceil_log2_tenthousand(W), 10 ** 4)
    alternatives = layout_scan(N, n_b)
    chosen = [al for al in alternatives if al["row"] == BUDGET and al["log2_S"] == ROUTE_PASSES * ROUTE_BITS
              and al["ram_log2"] == NODE_RAM_LOG2 and al["window"] == WINDOW and al["feasible"]]
    assert len(chosen) == 1
    assert abs(chosen[0]["log2_W"] - lg) < 1e-9            # the chosen layout is the ledger's own total
    assert chosen[0]["log2_capacity"] == log2_frac(cap)
    # ---- the four stages, each with its node count derived from C and its per-node inbound load ----------
    stages = [
        {"stage": "generate (pass-1 senders)", "nodes": 1 << 20,
         "records_per_node": "streams 2^96.6644, stores none",
         "held_in": "open-block area, 2^%.3f B" % math.log2(blocks_bytes)},
        {"stage": "pass 1 (digit 1)", "nodes": S, "records_per_node": log2_frac(mu),
         "inbound_bound": "N_a/2^37 + 26*2^20 under H3; digit-slab Chernoff 2^-2^%.3f per slab, union over "
                          "2^37 slabs 2^-2^73" % math.log2(slab_chernoff),
         "held_in": "pass buffer A"},
        {"stage": "pass 2 (digit 2, owner bucket)", "nodes": S, "records_per_node": log2_frac(mu),
         "inbound_bound": "exactly the owner bucket under H3; bucket Chernoff 2^-%s per node, union over "
                          "2^74 nodes 2^-%s" % (cstr, cstr), "held_in": "pass buffer B"},
        {"stage": "local sort and scan", "nodes": S, "records_per_node": log2_frac(mu),
         "held_in": "sort buffer pair (A,B)"},
        {"stage": "storage (W and base lists)", "nodes": 1 << 31,
         "records_per_node": "2^39.664 B of 16-byte fields", "held_in": "local RAM, <= 2^40 B"},
    ]
    return {
        "log2_m": log2_frac(m), "mu6_subset_H": mu6_in_H, "K_eff": float(K_eff),
        "success_birthday_model": birthday, "success_model_allowance": float(MODEL_ALLOW),
        "success_lower_bound": success - tails,
        "success_tail_sum": tails, "log2_N": log2_frac(N), "log2_n_b": log2_frac(n_b),
        "overhead_budgeted": item, "overhead_sum": sum(item.values()), "c": float(c),
        "per_node_rows_instr": {k_: float(v) for k_, v in node_rows.items()},
        "per_node_instr": node_work, "log2_per_node_instr": math.log2(node_work),
        "log2_2Nc": log2_frac(2 * N * c), "log2_other": log2_frac(other),
        "log2_rows": {k_: log2_frac(v) for k_, v in rows.items()},
        "log2_W": lg, "claim_score": [score.numerator, score.denominator, float(score)],
        "rounding_slack_bits": float(score) - lg,
        "log2_preprocessing": log2_frac(N * c + other),
        "fr_mul_instr": FR_MUL, "fr_inv_instr": FR_INV, "h_test_instr": H_TEST,
        "prime_search": {
            "A_cap_odd_candidates": A_ODDS, "A_cost_units": log2_frac(cost_A),
            "A_charged_units": log2_frac(primes_A), "W_cap_odd_candidates": W_ODDS,
            "W_cost_units": log2_frac(cost_W), "W_charged_units": log2_frac(primes_W),
            "probe": "every candidate below 2^%d is trial-divided by every prime below 2^%d" % (A_TOP, A_SIEVE),
            "costs": {"A": {"crossing_out_instr_log2": math.log2(sA["cross"]),
                            "per_segment_table_instr_log2": math.log2(sA["table"]),
                            "h_test_instr_log2": math.log2(sA["h_test"]),
                            "h_tests": float(sA["primes"]),
                            "h_primes_available_log2": math.log2(sA["h_primes"]),
                            "inv_prime_sum_bound": sA["inv_sum"]},
                      "W": {"crossing_out_instr_log2": math.log2(sW["cross"]),
                            "per_segment_table_instr_log2": math.log2(sW["table"]),
                            "h_test_instr_log2": math.log2(sW["h_test"]),
                            "h_tests": float(sW["primes"]),
                            "h_primes_available_log2": math.log2(sW["h_primes"]),
                            "inv_prime_sum_bound": sW["inv_sum"]}},
            "primality": "complete trial division to sqrt(q): every survivor of the sieve is prime, so "
                         "there is no Miller-Rabin base set, no primality probability and no false accept",
            "narrowing_rejected": "the first 13 prime bases are deterministic below 3.317e24 = 2^81.46 "
                                  "(Sorenson-Webster); H2 (a density prediction, not a bound) puts only "
                                  "pi(2^81.46)/2^20 = 2^55.67 H-primes below that bound against the 2^66.664 needed",
            "cap_policy": "both caps are deterministic functions of the fixed modulus, built before t and tau are "
                          "drawn and charged in full in every run whatever they find; the sieves count the H-primes "
                          "they find and the construction proceeds on those, so the search carries no probability "
                          "term in the success line"},
        "digits": {"owner_bits": ROUTE_PASSES * ROUTE_BITS, "node_ram_log2": NODE_RAM_LOG2,
                   "routing_digit_bits": ROUTE_BITS, "routing_passes": ROUTE_PASSES,
                   "lsd_digit_bits": LSD_BITS, "lsd_passes": LSD_PASSES, "sort_prefix_bits": PREFIX,
                   "window_bits": WINDOW, "compared_key_bits": 224 - PREFIX,
                   "block_records": BPASS, "block_bytes": BPASS * RECORD, "record_bytes": RECORD},
        "stages": stages,
        "field_layout": {"field_bytes": FIELD, "value_bits": 97, "base_prime_bits": 96, "tag_bits": 16, "reserved_bits": FIELD * 8 - 97 - 16,
                         "record_bytes": RECORD,
                         "label_bits": "67-bit base index + 50-bit step + 1 type bit = 118",
                         "read_path": "a label is resolved by a request/reply to the owning storage node: "
                                      "two 16-byte fields per verification, 2 x ~250 instructions, inside the "
                                      "2^12 units charged per verification"},
        "posting": {"per_record_per_pass": 100 // BPASS,
                    "blocks_per_node_pass_bound": "sum_d ceil(n_d/BPASS) <= M_v/BPASS + 2^37 (valid for every M_v)",
                    "global_formula": "2 x [(2N/25) x 100 + (2^74+2^20+2^31) x 2^37 x 107] "
                                      "+ 2^20 x 2^74 x 107 instructions",
                    "global_instr_log2": log2_frac(post_global),
                    "global_units_log2": log2_frac(post_global / 2000),
                    "full_blocks_instr_log2": log2_frac(post_full),
                    "final_blocks_instr_log2": log2_frac(post_final),
                    "pass1_end_tokens_instr_log2": log2_frac(107 * (1 << 20) * S),
                    "global_equals_carried_rows": carried == post_global,
                    "per_node_mu_premise": False,
                    "generating_node_records": log2_frac(Fraction(rec, 1 << 20)),
                    "log2_blocks_bound_from_mean": log2_frac(mu / BPASS + (1 << ROUTE_BITS)),
                    "final_blocks_per_node_pass": 1 << ROUTE_BITS, "final_block_instr": 107,
                    "open_block_area_bytes": blocks_bytes},
        "per_node_caps": {"run_pairs": RP_CAP, "verifications": V_CAP, "counter_bits": 64,
                          "log2_p_A_per_node": log2_frac(p_A), "log2_p_A_union": log2_frac(S * p_A),
                          "mean_Dprime_bound": float(lam), "Dprime_threshold_k": k,
                          "log2_p_B_per_node": log2_p_B, "log2_p_B_union": log2_p_B + log2_frac(S),
                          "log2_p_V_per_node": log2_frac(p_V), "log2_p_V_union": log2_frac(S * p_V)},
        "run_pair_cap_log2_units": log2_frac(rows["run_pair_cap"]),
        "verification_cap_log2_units": log2_frac(rows["verification_cap"]),
        "log2_mean_records_per_node": log2_frac(mu), "log2_node_capacity": math.log2(cap), "capacity_over_mean": a,
        "node_memory_aux_bytes": aux, "record_buffers_bytes": 96 * cap,
        "overflow_union_bound_log2": "-%s + %.1f" % (cstr, log2_frac(S)),
        "pass_inbound_union_bound_log2": "-2^%.1f" % math.log2(slab_chernoff),
        "pass_inbound": {"assumption": "H3: orbit keys behave as independent uniform strings at the relevant width",
                         "slab_mean_log2": math.log2(slab_mean), "threshold_factor": slab_a,
                         "rounding_records": 26 * (1 << 20),
                         "slabs_in_union": 1 << ROUTE_BITS,
                         "per_slab_exponent_log2": math.log2(slab_chernoff),
                         "pass1_union_bound": "2^37 slabs x 2^-2^%.3f = 2^-2^73" % math.log2(slab_chernoff),
                         "pass2_inbound": "exactly the owner bucket under H3, mean mu, no round-robin split and "
                                          "no marker record slot",
                         "pass2_union_bound": "2^74 nodes x 2^-%s = 2^-%s" % (cstr, cstr),
                         "deterministic_claim": "withdrawn: no pass's inbound is deterministically bounded by a "
                                                "slab count over a node count, and no for-all-passes statement",
                         "deduction_in_success": "2^-%s + 2^-2^73, both far below 2^-36.27" % cstr},
        "log2_memory_bytes": log2_frac(memory),
        "min_owner_nodes_2N_over_C": math.log2(rec / cap), "log2_S": math.log2(S),
        "gap_storage_bytes": (1 << 50) * FIELD, "base_list_bytes": 2 * n_b * FIELD,
        "storage_bytes_per_node": math.log2(((1 << 50) + 2 * n_b) * FIELD) - 31,
        "layout_scan": alternatives,
        "H2_fraction_of_primes_below_2^16_that_are_2^k_th_powers(k:measured)": h2,
        "H2_fraction_of_primes_2^16_to_2^17_that_are_2^k_th_powers(k:measured)": h2b,
    }


TOY = {  # name: (p, r, b, P0, k): y^2 = x^3 + b, #E = r prime, d = 2^k, 6d | r-1
    "t28": (207670579, 207671809, 3, (1, 207670577), 5),
    "t30": (903270253, 903326209, 17, (2, 903270248), 6),
    "t32": (3445477687, 3445418497, 3, (1, 2), 7),
}


def toy_sets(name, K, wave=4):
    p, r, b, P0, k = TOY[name]
    d = 1 << k
    m = (r - 1) // d
    Hp = [q for q in sieve(1 << 21) if pow(q, m, r) == 1]
    wlo = 1 << 11
    while True:                                   # chain start must exceed every gap
        W = [q for q in Hp if q > wlo][:60]
        if W[0] > max(W[j + 1] - W[j] for j in range(len(W) - 1)):
            break
        wlo *= 2
    L = len(W)
    nb = -(-int(math.sqrt(K * m / 6) + 1) // L)
    Ap = [q for q in Hp if q > W[-1]][:2 * nb]
    assert len(Ap) == 2 * nb
    return dict(p=p, r=r, b=b, P0=P0, k=k, d=d, m=m, W=W, L=L, nb=nb, wave=wave,
                baby=Ap[0::2], giant=Ap[1::2], K=6.0 * (nb * L) ** 2 / m)


def mu6(r):
    w = next(pow(g, (r - 1) // 3, r) for g in range(2, 50) if pow(g, (r - 1) // 3, r) != 1)
    return [1, r - 1, w, r - w, w * w % r, r - w * w % r]


def coverage(name, trials, seed, K):
    P = toy_sets(name, K)
    r, d, L = P["r"], P["d"], P["L"]
    U = mu6(r)
    tab = {u * a * w % r for a in P["baby"] for w in P["W"] for u in U}
    giants = [a * w for a in P["giant"] for w in P["W"]]
    rnd = random.Random(seed)
    succ = 0
    for _ in range(trials):
        z = pow(rnd.randrange(1, r), d, r)
        succ += any(x * z % r in tab for x in giants)
    Kd = P["K"] * (1 - (L - 1) / L ** 2)            # diagonal pairs w = w' repeat the ratio a/a'
    pred = 1 - math.exp(-Kd)
    se = math.sqrt(pred * (1 - pred) / trials)
    return {"toy": name, "log2_m": math.log2(P["m"]), "N": len(giants), "L": L, "trials": trials,
            "success": succ / trials, "poisson_K": 1 - math.exp(-P["K"]), "poisson_Kdiag": pred,
            "z_vs_Kdiag": (succ / trials - pred) / se}


def distinct_ratios(name, K):
    """#distinct x/x' over all baby/giant pairs, a lower bound on |S| (the u = 1 members of S).
    H1* allows 20% of the 6N^2 values to coincide; this measures the pairwise collision fraction."""
    P = toy_sets(name, K)
    r = P["r"]
    xs = [a * w % r for a in P["baby"] for w in P["W"]]
    ys = [a * w % r for a in P["giant"] for w in P["W"]]
    inv = [pow(y, -1, r) for y in ys]
    vals = {x * i % r for x in xs for i in inv}
    return {"toy": name, "pairs": len(xs) * len(ys), "distinct_ratios": len(vals),
            "collision_fraction": 1 - len(vals) / (len(xs) * len(ys))}


def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        l = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)


def ec_mul(k, P, p):
    R_ = None
    while k:
        if k & 1: R_ = ec_add(R_, P, p)
        P = ec_add(P, P, p)
        k >>= 1
    return R_


def batch_inv(vals, p):
    pre = [1]
    for v in vals: pre.append(pre[-1] * v % p)
    inv = pow(pre[-1], -1, p)
    out = [0] * len(vals)
    for i in range(len(vals) - 1, -1, -1):
        out[i] = inv * pre[i] % p
        inv = inv * vals[i] % p
    return out


def sqrt_mod(a, p):
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    c, x, t, mm = pow(z, q, p), pow(a, (q + 1) // 2, p), pow(a, q, p), s
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        x, c = x * b % p, b * b % p
        t, mm = t * c % p, i
    return x


def canon_y(P, p):
    return P[1] if P[1] % 2 == 0 else p - P[1]       # the even one of y, p - y


def ec_run(name, trials, seed, K, key_bits):
    P = toy_sets(name, K)
    p, r, b, P0, k, d = P["p"], P["r"], P["b"], P["P0"], P["k"], P["d"]
    W, L, wave = P["W"], P["L"], P["wave"]
    assert ec_mul(r, P0, p) is None
    gaps = [W[j + 1] - W[j] for j in range(L - 1)]
    U = mu6(r)
    assert all(canon_y(ec_mul(u, P0, p), p) == canon_y(P0, p) for u in U)       # key constant on the mu6 orbit
    omega_d = next(pow(g, (r - 1) // d, r) for g in range(2, 100) if pow(pow(g, (r - 1) // d, r), d // 2, r) != 1)
    mask = (1 << key_bits) - 1
    CAP = 4096
    key = lambda Y: (canon_y(Y, p) >> 1) & mask      # bits 1..key_bits of the even representative

    def side(base_pt, bases):
        recs = []
        for w0 in range(0, len(bases), wave):
            chain = bases[w0:w0 + wave]
            Q = [ec_mul(a, base_pt, p) for a in chain]
            tabs = [{g: ec_mul(g, q, p) for g in set(gaps)} for q in Q]
            cur = [ec_mul(W[0], q, p) for q in Q]
            recs += [(key(c), a * W[0]) for a, c in zip(chain, cur)]
            for j in range(L - 1):
                dens = [(tabs[i][gaps[j]][0] - cur[i][0]) % p for i in range(len(chain))]
                assert all(dens)                       # incomplete additions never degenerate
                inv = batch_inv(dens, p)
                for i in range(len(chain)):
                    (x1, y1), (x2, y2) = cur[i], tabs[i][gaps[j]]
                    l = (y2 - y1) * inv[i] % p
                    x3 = (l * l - x1 - x2) % p
                    cur[i] = (x3, (l * (x1 - x3) - y1) % p)
                    recs.append((key(cur[i]), chain[i] * W[j + 1]))
        return recs

    baby = side(P0, P["baby"])
    tabb = {u * a * w % r for a in P["baby"] for w in W for u in U}
    giants_x = [a * w for a in P["giant"] for w in W]
    nb_bins = 32
    shift = key_bits - 5
    cnt_b = [0] * nb_bins
    for kk, _ in baby: cnt_b[kk >> shift] += 1
    e_b = len(baby) / nb_bins
    chi_b = sum((c - e_b) ** 2 / e_b for c in cnt_b)
    chi_g = 0.0
    rnd = random.Random(seed)
    found = ok = disagree = false_total = aborts = 0
    root_matches = 0
    for _ in range(trials):
        s = rnd.randrange(1, r)
        z = pow(s, d, r)
        Ad, T_pt = ec_mul(z, P0, p), ec_mul(s, P0, p)
        in_S = any(x * z % r in tabb for x in giants_x)
        giant = side(Ad, P["giant"])
        cnt = [0] * nb_bins
        for kk, _ in giant: cnt[kk >> shift] += 1
        chi_g += sum((c - e_b) ** 2 / e_b for c in cnt)
        recs = sorted([(kk, 0, x) for kk, x in baby] + [(kk, 1, x) for kk, x in giant])
        zc, false_n, i = None, 0, 0
        while i < len(recs) and zc is None and false_n <= CAP:
            j = i
            while j < len(recs) and recs[j][0] == recs[i][0]: j += 1
            bs = [x for _k, t, x in recs[i:j] if t == 0]
            gs = [x for _k, t, x in recs[i:j] if t == 1]
            for xb in bs:
                for xg in gs:
                    cand = [xb * pow(u * xg % r, -1, r) % r for u in U]
                    hit = [cz for cz in cand if ec_mul(cz, P0, p) == Ad]
                    if hit: zc = hit[0]; break
                    false_n += 1
                if zc is not None: break
            i = j
        false_total += false_n
        if false_n > CAP:
            aborts += 1
            continue
        disagree += (zc is not None) != in_S
        if zc is None: continue
        found += 1
        assert zc == z
        # d-th root: take one root, then SELECT the root whose point is T (as the algorithm does) and only then
        # compare the recovered scalar with the secret. The d roots are root*omega_d^j, tested one scalar
        # multiplication each against the supplied T (only successful trials pay this).
        root = z
        for _ in range(k): root = sqrt_mod(root, r)
        cand = root
        rec_s, sel = None, None
        for jj in range(d):
            if ec_mul(cand, P0, p) == T_pt:
                rec_s, sel = cand, jj
                break
            cand = cand * omega_d % r
        assert rec_s is not None                     # exactly one of the d roots has point T
        root_matches += 1
        assert rec_s == s                            # the selection, not the known secret, produced s
        ok += 1
    dof = nb_bins - 1
    return {"toy": name, "trials": trials, "key_bits": key_bits, "found": found, "recovered_correct": ok,
            "root_selected_from_T": root_matches, "disagree_with_scalar_model": disagree,
            "false_key_pairs_per_trial": false_total / trials, "false_cap_aborts": aborts, "K_eff": P["K"],
            "occupancy_chi2_z_baby": (chi_b - dof) / math.sqrt(2 * dof),
            "occupancy_chi2_z_giants": (chi_g - trials * dof) / math.sqrt(2 * dof * trials)}


def dp(x, places, rounding=ROUND_HALF_UP):
    """x to a fixed number of decimals: Decimal(float) is exact and a Fraction is divided at 60 digits,
    so the result is the value the text displays."""
    if isinstance(x, Fraction):
        from decimal import localcontext
        with localcontext() as ctx:
            ctx.prec = 60
            x = Decimal(x.numerator) / Decimal(x.denominator)
    else:
        x = Decimal(x)
    return str(x.quantize(Decimal(1).scaleb(-places), rounding=rounding))


def sci(x):
    """One significant decimal, exponent written without a zero pad (5.0e-5, 1.2e-11)."""
    return ("%.1e" % x).replace("e-0", "e-").replace("e+0", "e")


def document(name):
    """proof.md sits beside this file; NOTE.md is the top-level note of the submission tree, so it is
    searched for upward from here and from the working directory (it is outside the snapshot the
    verifier copies, and then the check simply does not run)."""
    starts = [os.path.dirname(os.path.abspath(__file__)), os.path.abspath(os.getcwd())]
    for start in starts:
        d = start
        for _ in range(5):
            p = os.path.join(d, name)
            if os.path.isfile(p):
                with open(p, encoding="utf-8") as fh:
                    return p, fh.read()
            up = os.path.dirname(d)
            if up == d:
                break
            d = up
    return None, None


def quotes(L=None, strict=True):
    """Every derived quantity the text quotes, printed at the precision the text uses (Section 7.2 row
    totals and the shape table in the note, the capacity and the load, the exact total and the score,
    the memory bound and the success bound), plus the assertion that each printed string appears
    verbatim in the documents. The ledger is the authority: a string that is not found means the text
    is wrong, not the arithmetic."""
    L = L or ledger()
    scan = L["layout_scan"]
    cap = ((1 << NODE_RAM_LOG2) - L["node_memory_aux_bytes"]) // 96
    assert "2^" + dp(math.log2(cap), 4) == "2^" + dp(L["log2_node_capacity"], 4)
    q, check = {}, []

    def add(name, s, files=("proof.md", "NOTE.md")):
        q[name] = s
        check.append((s, tuple(files)))

    # the total in the three pieces the text quotes: the leading 2N term, every other row, the exact sum
    add("leading_2N_term", "2^" + dp(L["log2_2Nc"], 4))
    add("other_rows", "2^" + dp(L["log2_other"], 4), ("proof.md",))
    add("W_exact", "2^" + dp(L["log2_W"], 6))
    add("claim_fraction", "%d/%d" % tuple(L["claim_score"][:2]))
    add("score", dp(L["claim_score"][2], 4))
    add("rounding_slack", sci(L["rounding_slack_bits"]), ("proof.md",))
    # the per-record factor c: one record of each side carries the same 181-instruction row
    add("c_B", dp(L["c"], 4))
    add("c_G", dp(L["c"], 4))
    add("c_row210", dp(1 + Fraction(210, 2000), 3), ("NOTE.md",))
    add("2Nc_row210", "2^" + dp(L["log2_N"] + 1 + math.log2(1 + Fraction(210, 2000)), 4), ("NOTE.md",))
    add("2Nc_saving", dp(math.log2(1 + Fraction(210, 2000)) - math.log2(L["c"]), 3))
    # capacity, owner set and the load factor they leave
    add("log2_C", "2^" + dp(L["log2_node_capacity"], 4))
    add("C_records", str(cap), ("proof.md",))
    add("S", "2^%d" % int(round(L["log2_S"])))
    add("C_over_mu", dp(L["capacity_over_mean"], 4))
    # the two places the memory is quoted: the figure and the bound claimed on it
    add("memory", "2^" + dp(L["log2_memory_bytes"], 1), ("proof.md",))
    add("memory_records", "2^" + dp(L["log2_N"] + 1 + math.log2(2 * RECORD), 2), ("proof.md",))
    add("memory_claim", "2^124.1")
    # the success line, floored at the precision the text writes (it is a lower bound)
    add("success_bound", dp(L["success_lower_bound"], 7, ROUND_FLOOR))
    add("success_birthday", dp(L["success_birthday_model"], 7, ROUND_FLOOR))
    add("success_tails", sci(L["success_tail_sum"]))
    # preprocessing is an upper bound, so its display is a ceiling
    add("preprocessing_3dp", "2^" + dp(L["log2_preprocessing"], 3, ROUND_CEILING), ("proof.md",))
    add("preprocessing_2dp", "2^" + dp(L["log2_preprocessing"], 2, ROUND_CEILING))
    # the fixed-row items the text quotes from the ledger
    add("per_node_instr", "2^" + dp(L["log2_per_node_instr"], 3))
    add("per_node_fixed", "2^" + dp(L["log2_rows"]["per_node_fixed"], 4))
    add("posting_global", "2^" + dp(L["posting"]["global_instr_log2"], 4))
    add("run_pair_cap_3dp", "2^" + dp(L["log2_rows"]["run_pair_cap"], 3), ("proof.md",))
    add("run_pair_cap", "2^" + dp(L["log2_rows"]["run_pair_cap"], 2))
    add("verification_cap", "2^" + dp(L["log2_rows"]["verification_cap"], 0))
    # every Section 7.2 row: its total at 4 decimals, its fixed rows at 3 and its load at 2.  The proof's
    # table prints the total bare (the column is time_log2) and the note's shape table prefixes it.
    rows = []
    for i, al in enumerate(scan):
        row = {"row": i + 1, "passes": al["passes"], "route": al["route"], "lsd": al["lsd"],
               "total_exact": round(al["log2_W"], 9), "total": dp(al["log2_W"], 4),
               "fixed": "2^" + dp(al["log2_fixed"], 3),
               "load_exact": None if al["load"] is None else round(al["load"], 9), "load": None}
        if al["load"] is not None:
            row["load"] = dp(al["load"], 6 if i == 4 else 2)
        rows.append(row)
        if i in (5, 7):                       # row 6 fails capacity, row 8 fails on memory: no total, no row
            continue
        add("row%d_total" % (i + 1), row["total"], ("proof.md",))
        add("row%d_fixed" % (i + 1), row["fixed"], ("proof.md",))
        if i in (0, 1, 2, 4, 5, 6, 7, 8):     # the shapes the note's table lists beside the proof's
            add("row%d_total_prefixed" % (i + 1), "2^" + row["total"], ("NOTE.md",))
        if row["load"] is not None:
            add("row%d_load" % (i + 1), row["load"], ("NOTE.md",) if i in (0, 1, 2, 4, 6, 8) else ("proof.md",))
    # the two row differences the text argues from
    add("penalty_2x28", dp(scan[1]["log2_W"] - scan[2]["log2_W"], 4))
    add("gain_row7", dp(scan[2]["log2_W"] - scan[6]["log2_W"], 4), ("proof.md",))
    add("cost_2x38", dp(scan[4]["log2_W"] - scan[2]["log2_W"], 3))
    # the check: every printed string verbatim in the documents
    docs, text, missing = {}, {}, []
    for name in ("proof.md", "NOTE.md"):
        path, body = document(name)
        docs[name] = path or "(absent from this checkout)"
        text[name] = body
        if body is None and name == "proof.md":
            raise SystemExit("quote check: proof.md not found beside evidence.py")
    for s, files in check:
        for f in files:
            if text[f] is not None and s not in text[f]:
                missing.append("%s is not in %s" % (s, f))
    if missing and strict:
        raise SystemExit("quote desync: " + "; ".join(sorted(set(missing))))
    return {"quotes": dict(sorted(q.items())), "rows": rows,
            "asserted": [{"text": s, "documents": list(f)} for s, f in check],
            "desync": sorted(set(missing)), "documents": docs}


def main():
    arg = json.load(sys.stdin)
    mode, seed = arg["mode"], arg.get("seed", 2026)
    K = arg.get("K_permille", 500) / 1000.0
    if mode == "all":
        out = {"ledger": ledger(),
               "coverage": [coverage(t, n, seed, K) for t, n in (("t28", 8000), ("t30", 5000), ("t32", 3000))]
                           + [dict(distinct_ratios("t28", K), note="lower bound on |S|, u = 1 members")],
               "ec": [ec_run(t, n, seed, K, kb) for t, n, kb in (("t28", 80, 18), ("t30", 30, 20), ("t32", 20, 22))]}
        out["quotes"] = quotes(out["ledger"])
    elif mode == "ledger":
        out = ledger()
    elif mode == "quotes":
        out = quotes()
    elif mode == "coverage":
        out = coverage(arg["toy"], arg["trials"], seed, K)
    elif mode == "distinct":
        out = distinct_ratios(arg["toy"], K)
    elif mode == "ec":
        out = ec_run(arg["toy"], arg["trials"], seed, K, arg.get("key_bits", 18))
    else:
        raise SystemExit("unknown mode")
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

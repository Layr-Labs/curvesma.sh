"""Evidence for semiprime-cheon-v7 (stdlib only). N = elements per side, n_b = base primes per side.

stdin JSON {"mode": ..., ...}; "all" runs the three modes with fixed small trial counts.
  ledger    exact integer recomputation of the production ledger, cap and overflow bounds, memory sum
  coverage  scalar-level toy of the semiprime baby/giant sets vs the Poisson model
  ec        Algorithm A end to end on toy j=0 curves: affine walks with batch inversion, truncated parity
            keys, sort, equal-key runs under a false-match cap, d-th root, node-occupancy test
Toy results are evidence at toy scale only; the ledger is arithmetic on declared costs.
"""
import json, math, random, sys
from fractions import Fraction

R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D = 1 << 20


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


def ledger():
    m = (R - 1) // D
    assert (R - 1) % (6 * D) == 0
    omega = next(pow(g, (R - 1) // 3, R) for g in range(2, 50) if pow(g, (R - 1) // 3, R) != 1)
    mu6_in_H = pow(R - 1, m, R) == 1 and pow(omega, m, R) == 1
    K = Fraction(1, 2)                                     # K = 6 N^2 / m
    L = 1 << 50                                            # walked H-primes per chain
    n_b = -(-(math.isqrt(K.numerator * m // (K.denominator * 6)) + 2) // L)
    N = n_b * L
    ROUTE_BITS, BPASS, ROUTE_PASSES = 28, 25, 3             # 3 passes x 28 bits = 84-bit owner id
    LSD_BITS, LSD_PASSES = 28, 2                            # sort prefix 84 + 2 x 28 = 140 bits
    PREFIX = 84 + LSD_PASSES * LSD_BITS
    item = {"walk_bookkeeping(28 itemised)": 28, "key(23 itemised)": 23,
            "routing_3_passes_x_27": ROUTE_PASSES * (17 + 100 // BPASS + 6),
            "owner_drain(copy 6, capacity test 2, loop 2)": 10,
            "local_sort_2_passes_x_25": LSD_PASSES * 25, "run_scan": 10}
    BUDGET = 202                                          # itemised row, no percentage allowance
    assert sum(item.values()) == BUDGET
    c = 1 + Fraction(BUDGET, 2000)
    chains = 2 * n_b
    table = chains * (1 << 32) * Fraction(2012, 2000)      # 2^32 entries: one addition + 12 stores
    S = 1 << 84                                            # owner nodes of 2^40 bytes
    mu = Fraction(2 * N, S)
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
        "node_setup(code, stack, counters)": 4096,
    }
    node_work = sum(node_rows.values())
    assert node_work.denominator == 1
    node_work = int(node_work)
    nodes = S + (1 << 20) + (1 << 31)
    RP_CAP, V_CAP = 1 << 20, 1 << 7                        # per-node caps, 64-bit counters
    rows = {"tables": table, "starts": chains * (1 << 25),
            "per_node_fixed": Fraction(nodes * node_work, 2000),
            "primes_A": 1 << 97, "primes_W": 1 << 84, "input": 1 << 20, "lemma1": 1 << 21,
            "verification_cap": S * V_CAP * (1 << 12), "run_pair_cap": Fraction(S * RP_CAP * 32, 2000)}
    other = sum(rows.values())
    W = 2 * N * c + other
    rec = 2 * N
    # per-owner-node memory: the open-block area of a 28-bit pass (2^28 blocks of 25 records), the two record
    # buffers 96 C (48 source + 48 destination), the LSD count and offset arrays and the routing offset and
    # cursor tables (two 2^28-entry 8-byte arrays each), receive ring, code, reserved. Every forwarding node has it.
    blocks_bytes = (1 << ROUTE_BITS) * BPASS * 48
    aux = blocks_bytes + 2 * (1 << LSD_BITS) * 8 + 2 * (1 << ROUTE_BITS) * 8 + \
        (1 << 20) * 3072 + (1 << 30) + (1 << 36)
    cap = ((1 << 40) - aux) // 96
    a = float(cap) / float(mu)
    # occupancy: orbits are independent uniform over the owners, each carries at most 12 records, so the
    # Chernoff bound for [0,12]-bounded summands has exponent mu (a ln a - a + 1) / 12
    chernoff = float(mu) * (a * math.log(a) - a + 1) * math.log2(math.e) / 12
    assert chernoff > 1000 and cap * 96 + aux <= (1 << 40) and rec <= S * cap
    # per-node caps under H3 (orbit form), for a node holding M <= cap records, hence <= cap orbits:
    # A: some 56-bit value of key bits 85..140 is held by >= 7 orbits: <= cap^7 / (7! 2^(56*6))
    p_A = Fraction(cap ** 7, math.factorial(7) * (1 << (56 * 6)))
    # B: D' (orbits whose value repeats an earlier one) is dominated by Bin(cap, cap/2^56); outside A the
    # compared pairs are <= 1 + 72 * 3 * D', so exceeding RP_CAP needs D' >= k; P <= (e lam / k)^k
    lam = Fraction(cap * cap, 1 << 56)
    k = -(-(RP_CAP - 1) // 216)
    assert 1 + 216 * (k - 1) < RP_CAP and k > 3 * lam
    log2_p_B = -k * (math.log2(k) - math.log2(math.e * float(lam)))
    # V: verifications <= 1 + 72 * (orbit pairs with equal 224-bit key); two such pairs: <= cap^4 / 8 * 2^-280
    p_V = Fraction(cap ** 4, 8 * (1 << 280))
    assert V_CAP >= 1 + 72 * 1
    memory = (S + (1 << 20) + (1 << 31)) * (1 << 40)
    pr = sieve(1 << 16)
    h2 = {k_: sum(1 for q in pr if pow(q, (R - 1) >> k_, R) == 1) / len(pr) for k_ in (4, 8)}
    lg = log2_frac(W)
    sc = Fraction(math.ceil(lg * 10 ** 4), 10 ** 4)
    assert float(sc) >= lg
    assert W.numerator ** sc.denominator <= W.denominator ** sc.denominator * 2 ** sc.numerator
    success = 1 - math.exp(-float(Fraction(6 * N * N, m)))
    loss = float(S * p_A) + 2.0 ** (log2_p_B + 84) + float(S * p_V) + 2.0 ** -20
    return {
        "log2_m": log2_frac(m), "mu6_subset_H": mu6_in_H, "K_eff": float(Fraction(6 * N * N, m)),
        "success_poisson": success, "log2_N": log2_frac(N), "log2_n_b": log2_frac(n_b),
        "overhead_budgeted": item, "overhead_sum": sum(item.values()), "c": float(c),
        "per_node_rows_instr": {k_: float(v) for k_, v in node_rows.items()},
        "per_node_instr": node_work, "log2_per_node_instr": math.log2(node_work),
        "log2_2Nc": log2_frac(2 * N * c), "log2_other": log2_frac(other),
        "log2_rows": {k_: log2_frac(v) for k_, v in rows.items()},
        "log2_W": lg, "log2_preprocessing": log2_frac(N * c + other), "claim_score": [sc.numerator, sc.denominator, float(sc)],
        "digits": {"owner_bits": 84, "routing_digit_bits": ROUTE_BITS, "routing_passes": ROUTE_PASSES,
                   "lsd_digit_bits": LSD_BITS, "lsd_passes": LSD_PASSES, "sort_prefix_bits": PREFIX,
                   "block_records": BPASS, "block_bytes": BPASS * 48},
        "posting": {"per_record_per_pass": 100 // BPASS,
                    "log2_blocks_per_node_pass_bound": log2_frac(mu / BPASS + (1 << ROUTE_BITS)),
                    "final_blocks_per_node_pass": 1 << ROUTE_BITS, "final_block_instr": 107,
                    "open_block_area_bytes": blocks_bytes},
        "per_node_caps": {"run_pairs": RP_CAP, "verifications": V_CAP, "counter_bits": 64,
                          "log2_p_A_per_node": log2_frac(p_A), "log2_p_A_union": log2_frac(S * p_A),
                          "mean_Dprime_bound": float(lam), "Dprime_threshold_k": k,
                          "log2_p_B_per_node": log2_p_B, "log2_p_B_union": log2_p_B + 84,
                          "log2_p_V_per_node": log2_frac(p_V), "log2_p_V_union": log2_frac(S * p_V)},
        "run_pair_cap_log2_units": log2_frac(rows["run_pair_cap"]),
        "verification_cap_log2_units": log2_frac(rows["verification_cap"]),
        "success_loss": loss, "success_lower_bound": success - loss,
        "log2_mean_records_per_node": log2_frac(mu), "log2_node_capacity": math.log2(cap), "capacity_over_mean": a,
        "node_memory_aux_bytes": aux, "record_buffers_bytes": 96 * cap,
        "overflow_union_bound_log2": "-%.3e + 84" % chernoff, "log2_memory_bytes": log2_frac(memory),
        "H2_fraction_of_primes_below_2^16_that_are_2^k_th_powers(k:measured)": h2,
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
        root = z
        for _ in range(k): root = sqrt_mod(root, r)
        sc = [root * pow(omega_d, j, r) % r for j in range(d)]
        assert s in sc and ec_mul(s, P0, p) == T_pt
        ok += 1
    dof = nb_bins - 1
    return {"toy": name, "trials": trials, "key_bits": key_bits, "found": found, "recovered_correct": ok,
            "disagree_with_scalar_model": disagree, "false_key_pairs_per_trial": false_total / trials,
            "false_cap_aborts": aborts, "K_eff": P["K"],
            "occupancy_chi2_z_baby": (chi_b - dof) / math.sqrt(2 * dof),
            "occupancy_chi2_z_giants": (chi_g - trials * dof) / math.sqrt(2 * dof * trials)}


def main():
    arg = json.load(sys.stdin)
    mode, seed = arg["mode"], arg.get("seed", 2026)
    K = arg.get("K_permille", 500) / 1000.0
    if mode == "all":
        out = {"ledger": ledger(),
               "coverage": [coverage(t, n, seed, K) for t, n in (("t28", 8000), ("t30", 5000), ("t32", 3000))],
               "ec": [ec_run(t, n, seed, K, kb) for t, n, kb in (("t28", 80, 18), ("t30", 30, 20), ("t32", 20, 22))]}
    elif mode == "ledger":
        out = ledger()
    elif mode == "coverage":
        out = coverage(arg["toy"], arg["trials"], seed, K)
    elif mode == "ec":
        out = ec_run(arg["toy"], arg["trials"], seed, K, arg.get("key_bits", 18))
    else:
        raise SystemExit("unknown mode")
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()

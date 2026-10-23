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
    K = Fraction(13, 25)                                   # K = 6 N^2 / m
    L = 1 << 50                                            # walked H-primes per chain
    n_b = -(-(math.isqrt(K.numerator * m // (K.denominator * 6)) + 2) // L)
    N = n_b * L
    item = {"walk_bookkeeping(27 listed)": 40, "key(23 listed)": 30, "routing_4_hops_x_19": 76,
            "local_sort_3_passes_x_25": 75, "run_scan": 10}
    assert sum(item.values()) <= 300
    c = 1 + Fraction(300, 2000)
    chains = 2 * n_b
    table = chains * (1 << 32) * Fraction(2012, 2000)      # 2^32 entries: one addition + 12 stores
    S = 1 << 85                                            # owner nodes of 2^40 bytes
    rows = {"tables": table, "starts": chains * (1 << 25), "per_node_fixed": Fraction(S * (1 << 31), 2000),
            "primes_A": 1 << 97, "primes_W": 1 << 84, "input": 1 << 20, "lemma1": 1 << 21,
            "false_match_cap": 1 << 32, "run_pair_cap": 1 << 84}
    other = sum(rows.values())
    W = 2 * N * c + other
    rec = 2 * N
    pairs_prefix = Fraction(rec * (rec - 1), 2) / (1 << 154)
    pairs_full = Fraction(rec * (rec - 1), 2) / (1 << 224)
    mu = Fraction(rec, S)
    reserved = 1 << 36
    cap = ((1 << 40) - reserved) // 96
    a = float(cap) / float(mu)
    chernoff = float(mu) * (a * math.log(a) - a + 1) * math.log2(math.e)
    assert chernoff > 1000 and cap * 96 + reserved <= (1 << 40) and rec <= S * cap
    memory = (S + (1 << 20) + (1 << 31)) * (1 << 40)
    pr = sieve(1 << 16)
    h2 = {k: sum(1 for q in pr if pow(q, (R - 1) >> k, R) == 1) / len(pr) for k in (4, 8)}
    lg = log2_frac(W)
    sc = Fraction(math.ceil(lg * 10 ** 4), 10 ** 4)
    assert float(sc) >= lg
    return {
        "log2_m": log2_frac(m), "mu6_subset_H": mu6_in_H, "K_eff": float(Fraction(6 * N * N, m)),
        "success_poisson": 1 - math.exp(-6 * N * N / m), "log2_N": log2_frac(N), "log2_n_b": log2_frac(n_b),
        "overhead_budgeted": item, "overhead_sum": sum(item.values()), "c": float(c),
        "log2_2Nc": log2_frac(2 * N * c), "log2_other": log2_frac(other),
        "log2_rows": {k: log2_frac(v) for k, v in rows.items()},
        "log2_W": lg, "log2_preprocessing": log2_frac(N * c + other), "claim_score": [sc.numerator, sc.denominator, float(sc)],
        "log2_prefix_pairs": log2_frac(pairs_prefix), "log2_p_pair_cap": log2_frac(pairs_prefix / (1 << 90)),
        "log2_false_pairs": log2_frac(pairs_full), "log2_p_false_cap": log2_frac(pairs_full / (1 << 20)),
        "log2_mean_records_per_node": log2_frac(mu), "log2_node_capacity": math.log2(cap), "capacity_over_mean": a,
        "overflow_union_bound_log2": "-%.3e + 85" % chernoff, "log2_memory_bytes": log2_frac(memory),
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
    K = arg.get("K_permille", 520) / 1000.0
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

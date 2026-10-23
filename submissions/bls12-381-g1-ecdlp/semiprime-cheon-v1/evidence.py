"""Evidence for the H-prime semiprime Cheon walk submissions (stdlib only).
(Here T below is the number of elements per side, called N in proof.md.)

modes (stdin JSON {"mode": ..., ...}); "all" runs the three below with fixed small trial counts:
  ledger   exact recomputation of the production cost ledger and algebra facts
  coverage scalar-level toy of the semiprime baby/giant walk: success and
           early-stopping statistics against the Poisson model
  ec       toy BLS-like j=0 curve end to end: batched incomplete affine walks,
           canonical-y keys, unified table, candidate z, d-th root, mu_d scan
Everything here is toy scale evidence or exact arithmetic on the public
production constants; it does not prove a production complexity bound.
"""
import json, math, random, sys
from fractions import Fraction

R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D = 1 << 20

# ---------------------------------------------------------------- ledger ---
def isqrt_ceil(n):
    s = math.isqrt(n)
    return s if s * s == n else s + 1

def _log2_int(n):
    sh = max(n.bit_length() - 60, 0)
    return math.log2(n >> sh) + sh

def log2_frac(x):
    """log2 of a positive Fraction/int as float, via integer bit length."""
    x = Fraction(x)
    return _log2_int(x.numerator) - _log2_int(x.denominator)

def trapezoid_free_I(K):
    # int_0^1 exp(-K x^2) dx = sqrt(pi)/2 * erf(sqrt K)/sqrt K
    return math.sqrt(math.pi) / 2 * math.erf(math.sqrt(K)) / math.sqrt(K)

def ledger():
    m = (R - 1) // D
    assert (R - 1) % D == 0 and (R - 1) % (6 * D) == 0
    v2 = ((R - 1) & -(R - 1)).bit_length() - 1
    # facts used by the attack
    two_in_H = pow(2, m, R) == 1
    omega = None
    for g in range(2, 50):
        w = pow(g, (R - 1) // 3, R)
        if w != 1:
            omega = w
            break
    mu6_in_H = pow(R - 1, m, R) == 1 and pow(omega, m, R) == 1
    K = Fraction(13, 25)                      # coverage parameter K = 6 T^2 / m
    L = 1 << 50                               # walked H-primes per chain
    T0 = isqrt_ceil(K.numerator * m // (K.denominator * 6)) + 1
    n = -(-T0 // L)                           # base primes per side
    T = n * L                                 # elements per side
    Keff = Fraction(6 * T * T, m)
    # per-element overhead in AArch64 instructions beyond the reference addition
    item = {"walk_bookkeeping": 40, "canonical_y_key": 30, "local_probe_insert": 61,
            "routing_6_hops_x_18": 6 * 18}
    itemised = sum(item.values())
    budget = 300
    assert itemised <= budget
    c = 1 + Fraction(budget, 2000)
    chains = 2 * n
    G = 1 << 38                               # table gap bound; only even gaps
    table_adds = chains * (G // 2)
    table_store = Fraction(12 * table_adds, 2000)
    chain_init = chains * 512                 # [a]P0, [b0]Q: <= 2*(94+36) ops < 512
    chunk = chains * (1 << 38)                # sum of gaps <= 2^76 => <= 2^76/G extra adds
    prime_gen_A = 1 << 96                     # H-prime search up to ~2^93: ~2^91 + 4x slack
    prime_gen_W = 1 << 84                     # H-primes below ~2^76: ~2^75 + slack
    input_access = 1 << 20                    # locate/decode A_d and T in the input file
    final_step = 1 << 23                      # false matches, z candidates, 20 sqrt, mu_d scan
    other = (table_adds + table_store + chain_init + chunk + prime_gen_A
             + prime_gen_W + input_access + final_step)
    elems = 2 * T * c
    W_worst = elems + other
    I = Fraction(trapezoid_free_I(float(Keff))).limit_denominator(10 ** 15)
    W_exp = elems * I + other
    pr = sieve(1 << 16)
    h2 = {}
    for k in (4, 8):
        e = (R - 1) >> k
        hits = sum(1 for q in pr if pow(q, e, R) == 1)
        h2["2^%d-th_power_fraction_of_primes_below_2^16" % k] = [hits, len(pr), hits / len(pr), 1.0 / (1 << k)]
    out = {
        "H2_sanity_equidistribution": h2,
        "log2_r": log2_frac(R), "log2_m": log2_frac(m), "v2_r_minus_1": v2,
        "d_divides_r_minus_1": True, "6d_divides_r_minus_1": True,
        "mu6_subset_H": mu6_in_H, "two_in_H": two_in_H,
        "K_target": float(K), "K_eff": float(Keff),
        "success_poisson": 1 - math.exp(-float(Keff)),
        "log2_N": log2_frac(T), "log2_n_b": log2_frac(n),
        "overhead_itemised": item, "overhead_itemised_sum": itemised,
        "overhead_budget": budget, "c": float(c),
        "log2_2Tc": log2_frac(elems), "log2_other": log2_frac(other),
        "I_K": float(I),
        "log2_W_worst": log2_frac(W_worst), "log2_W_expected": log2_frac(W_exp),
        "other_terms_log2": {
            "tables": log2_frac(table_adds + table_store),
            "chain_init": log2_frac(chain_init), "gap_chunk_bound": log2_frac(chunk),
            "prime_gen_A": 96, "prime_gen_W": 84, "input_access": 20, "final_step": 23},
        "memory_log2_bytes_static": log2_frac(Fraction(T * 48 * 10, 9)),   # bucketed table, load 0.9
        "memory_log2_bytes_unified": log2_frac(Fraction(2 * T * 48 * 10, 9)),
    }
    for key, w in (("worst", W_worst), ("expected", W_exp)):
        lg = log2_frac(w)
        sc = Fraction(math.ceil(lg * 10 ** 4), 10 ** 4)      # round up to 4 decimals
        assert float(sc) >= lg
        out["claim_score_" + key] = {"numerator": str(sc.numerator),
                                     "denominator": str(sc.denominator),
                                     "decimal": float(sc)}
    return out

# --------------------------------------------------------- toy parameters ---
TOY = {
    # name: (p, r, b, P0, k) with d = 2^k, y^2 = x^3 + b, #E = r prime, 6*2^k | r-1
    "t28": (207670579, 207671809, 3, (1, 207670577), 5),
    "t30": (903270253, 903326209, 17, (2, 903270248), 6),
    "t32": (3445477687, 3445418497, 3, (1, 2), 7),
}

def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def toy_sets(name, K=0.52, L=None, wave=4, wlo=None):
    p, r, b, P0, k = TOY[name]
    d = 1 << k
    m = (r - 1) // d
    Tt = math.sqrt(K * m / 6)
    primes = sieve(1 << 21)
    Hp = [q for q in primes if pow(q, m, r) == 1]
    wlo = wlo or (1 << 11)
    while True:                              # chain start must exceed every gap
        Wl = [q for q in Hp if q > wlo][: (L or 60)]
        if Wl[0] > max(Wl[j + 1] - Wl[j] for j in range(len(Wl) - 1)):
            break
        wlo *= 2
    L = len(Wl)
    amin = Wl[-1] + 1
    nb = -(-int(Tt + 1) // L)
    Ap = [q for q in Hp if q > amin][: 2 * nb]
    assert len(Ap) == 2 * nb, "increase prime range"
    baby_bases, giant_bases = Ap[0::2], Ap[1::2]
    return dict(p=p, r=r, b=b, P0=P0, k=k, d=d, m=m, W=Wl, L=L, nb=nb, wave=wave,
                baby_bases=baby_bases, giant_bases=giant_bases, T=nb * L,
                K=6.0 * (nb * L) ** 2 / m)

def order_stream(bases, W, wave):
    """yield (x=a*b, t) in lockstep-wave order; t = elements of this side generated
    up to and including the round."""
    out = []
    done = 0
    L = len(W)
    for w0 in range(0, len(bases), wave):
        chain = bases[w0:w0 + wave]
        for j in range(L):
            t = done + (j + 1) * len(chain)
            for a in chain:
                out.append((a * W[j], t))
        done += L * len(chain)
    return out

def mu6(r):
    for g in range(2, 50):
        w = pow(g, (r - 1) // 3, r)
        if w != 1:
            break
    return [1, r - 1, w, r - w, w * w % r, r - w * w % r]

def coverage(name, trials, seed, K):
    P = toy_sets(name, K=K)
    r, d, m, T = P["r"], P["d"], P["m"], P["T"]
    U = mu6(r)
    baby = order_stream(P["baby_bases"], P["W"], P["wave"])
    giant = order_stream(P["giant_bases"], P["W"], P["wave"])
    assert len({x for x, _ in baby}) == len(baby) == T
    tab = {}
    dup_orbit = 0
    for x, t in baby:
        for u in U:
            v = u * x % r
            if v in tab:
                if tab[v][0] != x:
                    dup_orbit += 1
                tab[v] = (tab[v][0], min(tab[v][1], t))
            else:
                tab[v] = (x, t)
    rnd = random.Random(seed)
    succ = 0
    tot = 0.0
    for _ in range(trials):
        s = rnd.randrange(1, r)
        z = pow(s, d, r)
        best = None
        for x2, tg in giant:
            if best is not None and tg >= best:
                break
            hit = tab.get(x2 * z % r)
            if hit is not None:
                cand = max(hit[1], tg)
                if best is None or cand < best:
                    best = cand
        if best is not None:
            succ += 1
            tot += best
        else:
            tot += T
    Kp = P["K"]
    L = P["L"]
    Kd = Kp * (1 - (L - 1) / L ** 2)     # diagonal pairs w = w' repeat the ratio a/a'
    pred_s = 1 - math.exp(-Kp)
    pred_d = 1 - math.exp(-Kd)
    ph = succ / trials
    se_s = math.sqrt(pred_d * (1 - pred_d) / trials)
    return {"toy": name, "r": r, "d": d, "log2_m": math.log2(m), "T_per_side": T,
            "chain_len_L": L, "base_primes_per_side": P["nb"], "K_eff": Kp,
            "K_diag_corrected": Kd, "trials": trials, "success_measured": ph,
            "success_poisson": pred_s, "success_poisson_diag_corrected": pred_d,
            "success_zscore_vs_diag_corrected": (ph - pred_d) / se_s,
            "mean_stop_over_T_measured": tot / trials / T,
            "mean_stop_over_T_poisson": trapezoid_free_I(Kp),
            "mean_stop_over_T_poisson_diag_corrected": trapezoid_free_I(Kd),
            "baby_orbit_duplicates": dup_orbit, "baby_elements": len(baby)}

# ----------------------------------------------------------------- toy EC ---
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
    for v in vals:
        pre.append(pre[-1] * v % p)
    inv = pow(pre[-1], -1, p)
    out = [0] * len(vals)
    for i in range(len(vals) - 1, -1, -1):
        out[i] = inv * pre[i] % p
        inv = inv * vals[i] % p
    return out

def sqrt_mod(a, p):
    a %= p
    q = p - 1; s = 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    c = pow(z, q, p); x = pow(a, (q + 1) // 2, p); t = pow(a, q, p); mm = s
    while t != 1:
        i = 0; tt = t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        x = x * b % p; c = b * b % p; t = t * c % p; mm = i
    return x

def dth_root(z, r, k):
    w = z
    for _ in range(k):
        w = sqrt_mod(w, r)
    return w

def canon_y(P, p):
    y = P[1]
    return y if y % 2 == 0 else p - y     # p odd: exactly one of y, p-y is even

def ec_run(name, trials, seed, K):
    P = toy_sets(name, K=K)
    p, r, b, P0, k, d, m = P["p"], P["r"], P["b"], P["P0"], P["k"], P["d"], P["m"]
    W, L, wave = P["W"], P["L"], P["wave"]
    assert (P0[1] ** 2 - P0[0] ** 3 - b) % p == 0 and ec_mul(r, P0, p) is None
    gaps = [W[j + 1] - W[j] for j in range(L - 1)]
    assert W[0] > max(gaps)                 # b_j > every gap => P != +-D, no degenerate add
    # GLV check: [omega]P0 = (beta*x, y) for some cube roots; key is orbit invariant
    beta = [pow(g, (p - 1) // 3, p) for g in range(2, 40) if pow(g, (p - 1) // 3, p) != 1][0]
    omegas = [w for w in mu6(r) if w not in (1, r - 1)]
    lam = None
    for w in omegas:
        Q = ec_mul(w, P0, p)
        if Q[1] == P0[1] and Q[0] == beta * P0[0] % p: lam = w
        if Q[1] == P0[1] and Q[0] == beta * beta % p * P0[0] % p: lam = w
    assert lam is not None
    for u in mu6(r):
        assert canon_y(ec_mul(u, P0, p), p) == canon_y(P0, p)
    Omega_d = None
    for g in range(2, 100):
        w = pow(g, (r - 1) // d, r)
        if pow(w, d // 2, r) != 1: Omega_d = w; break

    def side(base_pt, bases, mult_base):
        """lockstep waves; yields rounds of (label, key, t). base_pt scalar-multiplied by a."""
        rounds = []
        done = 0
        for w0 in range(0, len(bases), wave):
            chain = bases[w0:w0 + wave]
            Q = [ec_mul(a, base_pt, p) for a in chain]
            tabs = []
            for q in Q:
                tab = {}
                for g_ in set(gaps):
                    tab[g_] = ec_mul(g_, q, p)
                tabs.append(tab)
            cur = [ec_mul(W[0], q, p) for q in Q]
            keys0 = [canon_y(c_, p) for c_ in cur]
            rounds.append((done + len(chain), [(a, W[0], key) for a, key in zip(chain, keys0)]))
            for j in range(L - 1):
                g_ = gaps[j]
                dens = [(tabs[i][g_][0] - cur[i][0]) % p for i in range(len(chain))]
                assert all(dd != 0 for dd in dens)          # incomplete addition never degenerates
                invs = batch_inv(dens, p)
                row = []
                for i in range(len(chain)):
                    x1, y1 = cur[i]; x2, y2 = tabs[i][g_]
                    l = (y2 - y1) * invs[i] % p
                    x3 = (l * l - x1 - x2) % p
                    cur[i] = (x3, (l * (x1 - x3) - y1) % p)
                    row.append((chain[i], W[j + 1], canon_y(cur[i], p)))
                rounds.append((done + (j + 2) * len(chain), row))
            done += L * len(chain)
        return rounds

    baby_rounds = side(P0, P["baby_bases"], 1)
    rnd = random.Random(seed)
    ok = 0
    found = 0
    mism = 0
    scal_baby = order_stream(P["baby_bases"], P["W"], wave)
    scal_giant = order_stream(P["giant_bases"], P["W"], wave)
    U = mu6(r)
    tabb = {}
    for x, t in scal_baby:
        for u in U:
            v = u * x % r
            tabb[v] = min(tabb.get(v, 10 ** 9), t)
    for _ in range(trials):
        s = rnd.randrange(1, r)
        z = pow(s, d, r)
        Ad = ec_mul(z, P0, p)
        T_pt = ec_mul(s, P0, p)
        # scalar model prediction of the stopping round
        best = None
        for x2, tg in scal_giant:
            hit = tabb.get(x2 * z % r)
            if hit is not None:
                cand = max(hit, tg)
                if best is None or cand < best: best = cand
        giant_rounds = side(Ad, P["giant_bases"], 1)
        table = {}
        tau = None
        res = None
        bi = gi = 0
        # unified table, strictly alternating rounds baby/giant as in the algorithm
        while bi < len(baby_rounds) or gi < len(giant_rounds):
            for typ, rounds in ((0, baby_rounds), (1, giant_rounds)):
                idx = bi if typ == 0 else gi
                if idx >= len(rounds): continue
                t, row = rounds[idx]
                for a, bb, key in row:
                    ent = table.get(key)
                    if ent is None:
                        table[key] = (typ, a * bb)
                    elif ent[0] != typ and res is None:
                        res = (ent[1], a * bb, typ)
                if typ == 0: bi += 1
                else: gi += 1
            if res is not None:
                tau = max(rounds_t for rounds_t in (baby_rounds[bi - 1][0], giant_rounds[gi - 1][0]))
                break
        if res is None:
            assert best is None
            continue
        found += 1
        xb, xg = (res[0], res[1]) if res[2] == 1 else (res[1], res[0])
        zc = None
        for u in U:
            cand_z = xb * pow(u * xg % r, -1, r) % r
            if ec_mul(cand_z, P0, p) == Ad: zc = cand_z
        assert zc == z, "recovered z mismatch"
        root = dth_root(z, r, k)
        sc = None
        for j in range(d):
            sj = root * pow(Omega_d, j, r) % r
            if ec_mul(sj, P0, p) == T_pt: sc = sj
        assert sc == s, "recovered s mismatch"
        ok += 1
        # the EC-level stopping round must equal the scalar model's
        if best is None or tau != best: mism += 1
    return {"toy": name, "trials": trials, "found": found, "recovered_correct": ok,
            "stop_round_differs_from_scalar_model": mism, "K_eff": P["K"],
            "success_poisson": 1 - math.exp(-P["K"]), "chain_len": L,
            "bases_per_side": P["nb"], "incomplete_adds_nondegenerate": True,
            "canonical_y_orbit_invariant": True}

def main():
    arg = json.load(sys.stdin)
    mode = arg["mode"]
    seed = arg.get("seed", 2026)
    Kp = arg.get("K_permille", 520) / 1000.0
    if mode == "all":
        out = {"ledger": ledger(),
               "coverage": [coverage(t, n, seed, Kp) for t, n in (("t28", 8000), ("t30", 5000), ("t32", 3000))],
               "ec": [ec_run(t, n, seed, Kp) for t, n in (("t28", 80), ("t30", 30), ("t32", 20))]}
    elif mode == "ledger":
        out = ledger()
    elif mode == "coverage":
        out = coverage(arg["toy"], arg["trials"], seed, Kp)
    elif mode == "ec":
        out = ec_run(arg["toy"], arg["trials"], seed, Kp)
    else:
        raise SystemExit("unknown mode")
    out["mode"] = mode
    print(json.dumps(out, sort_keys=True))

if __name__ == "__main__":
    main()

#!/usr/bin/env python3 -I
"""Evidence for the paired class-walk Cheon attack on bls12-381-g1-ecdlp.

Modes (input object on stdin: {"mode": "...", ...}):
  ledger  -- every cost term of proof.md from exact integers, the exact score
             rounded up, preprocessing and memory rounded up, success rounded
             down; asserts the inequalities. Arithmetic only.
  planner -- the offline pairing planner on real chi-class lists of toy primes
             at the production degree parameter GMAX = 16 d^2 (evidence for H4).
  ec      -- the whole algorithm on a toy curve: paired walk, bucket matching
             with truncated keys, false matches, bounded re-walk, recovery.
  census  -- exact coverage |R|/m of the construction on toy primes (H1).
  all     -- ledger + planner + ec + census at declared hosted sizes.

Standard library only. The hosted run uses the small sizes that fit the
runner's time budget; the offline rows quoted in proof.md use the same code
paths with larger inputs on a larger machine.
"""
import json, math, random, sys, time
from fractions import Fraction

R_PRODUCTION = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D_PRODUCTION = 1 << 20
M_PRODUCTION = (R_PRODUCTION - 1) // D_PRODUCTION


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


def toy_prime(bits, dd, rng):
    while True:
        rr = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * dd) * (6 * dd) + 1
        if is_prime(rr) and ((rr - 1) // dd) % 6 == 0 and pow(2, (rr - 1) // dd, rr) != 1:
            return rr


def character(rr, dd):
    """chi: integers -> Z/d, chi(x) = log_zeta(x^m), kernel H."""
    mm = (rr - 1) // dd
    gen = next(g for g in range(2, 10 ** 6)
               if all(pow(g, (rr - 1) // q, rr) != 1 for q in (2, 3)))
    zeta = pow(gen, mm, rr)
    zlog = {pow(zeta, i, rr): i for i in range(dd)}
    return (lambda x: zlog[pow(x, mm, rr)]), gen


def class_lists(rr, dd, X2):
    chi, _ = character(rr, dd)
    classes = [[] for _ in range(dd)]
    for b in range(1, X2 + 1):
        classes[chi(b)].append(b)
    return classes


# ------------------------------------------------------------------ planner (v2, validated)
def plan(L, GMAX, max_depth=1 << 14, center_priority=False):
    """Offline pairing plan for one sorted class list L.

    Pair (i, j) with midpoint k: L[k] - L[i] = L[j] - L[k] = g, 1 <= g <= GMAX.
    Children get parent k (in-degree 1); roots are one table addition from the
    previous root (gap <= GMAX enforced); pairs are rejected if they would
    close a dependency cycle or exceed the depth cap. Greedy matching by
    increasing degree, no augmenting paths, final root-chain pass only demotes.
    Emission: list order; each root then its whole subtree, depth first.
    """
    n = len(L)
    pos = {v: i for i, v in enumerate(L)}
    adj = [[] for _ in range(n)]
    n_edges = 0
    for i in range(n):
        vi = L[i]
        lim = vi + 2 * GMAX
        j = i + 1
        while j < n and L[j] <= lim:
            vj = L[j]
            if (vi + vj) & 1 == 0:
                k = pos.get((vi + vj) >> 1)
                if k is not None:
                    adj[i].append((j, k))
                    adj[j].append((i, k))
                    n_edges += 1
            j += 1
    parent, partner, depth = [-1] * n, [-1] * n, [0] * n
    pair_slot = [-1] * n
    height = [0] * n
    children = [[] for _ in range(n)]
    stats = {"cycle_rejections": 0, "depth_rejections": 0, "demoted": 0,
             "depth_increments": 0, "depth_decrements": 0, "height_increments": 0}

    def set_depth(v):
        stack = [v]
        while stack:
            x = stack.pop()
            next_depth = depth[parent[x]] + 1 if parent[x] != -1 else 0
            delta = next_depth - depth[x]
            stats["depth_increments"] += max(0, delta)
            stats["depth_decrements"] += max(0, -delta)
            depth[x] = next_depth
            for (a, b) in children[x]:
                stack.append(a)
                stack.append(b)

    def creates_cycle(i, j, k):
        x = k
        while x != -1:
            if x == i or x == j:
                return True
            x = parent[x]
        return False

    def try_match(i, j, k, force=False):
        if parent[i] != -1 or parent[j] != -1:
            return False
        if not force:
            if depth[k] + 1 + max(height[i], height[j]) > max_depth:
                stats["depth_rejections"] += 1
                return False
            if creates_cycle(i, j, k):
                stats["cycle_rejections"] += 1
                return False
        parent[i] = parent[j] = k
        partner[i], partner[j] = j, i
        pair_slot[i] = pair_slot[j] = len(children[k])
        children[k].append((min(i, j), max(i, j)))
        set_depth(i)
        set_depth(j)
        h, x = 1 + max(height[i], height[j]), k
        while x != -1 and height[x] < h:
            stats["height_increments"] += h - height[x]
            height[x] = h
            h += 1
            x = parent[x]
        return True

    def unmatch(u):
        k, b = parent[u], partner[u]
        position = pair_slot[u]
        assert position == pair_slot[b] and children[k][position] == (min(u, b), max(u, b))
        moved = children[k][-1]
        children[k][position] = moved
        pair_slot[moved[0]] = pair_slot[moved[1]] = position
        children[k].pop()
        pair_slot[u] = pair_slot[b] = -1
        parent[u] = parent[b] = -1
        partner[u] = partner[b] = -1
        set_depth(u)
        set_depth(b)

    if center_priority:
        by_center = [[] for _ in range(n)]
        for i, edges in enumerate(adj):
            for j, k in edges:
                if i < j:
                    by_center[k].append((i, j))
        for k in sorted(range(n), key=lambda k: -len(by_center[k])):
            for i, j in by_center[k]:
                if parent[i] == -1 and parent[j] == -1:
                    try_match(i, j, k)
    for v in sorted(range(n), key=lambda v: len(adj[v])):
        if parent[v] != -1:
            continue
        for (w, k) in adj[v]:
            if parent[w] == -1 and try_match(v, w, k):
                break
    if parent[0] != -1:
        unmatch(0)
        stats["demoted"] += 1
    last = L[0]
    for v in range(1, n):
        if L[v] - last > GMAX:
            u = v - 1
            if parent[u] != -1:
                unmatch(u)
                stats["demoted"] += 1
            last = L[u]
            assert L[v] - last <= GMAX
        if parent[v] == -1:
            last = L[v]
    emitted = [False] * n
    ops = []
    max_stack = n_pairs = n_single = 0
    last_root = None
    for v in range(n):
        if emitted[v] or parent[v] != -1:
            continue
        if last_root is None:
            ops.append((0, v))
        else:
            ops.append((1, v, L[v] - L[last_root]))
            n_single += 1
        emitted[v] = True
        last_root = v
        stack = [(v, 0)]
        while stack:
            max_stack = max(max_stack, len(stack))
            x, idx = stack[-1]
            if idx < len(children[x]):
                stack[-1] = (x, idx + 1)
                a, b = children[x][idx]
                g = L[x] - L[a]
                assert L[b] - L[x] == g and 0 < g <= GMAX
                ops.append((2, x, a, b, g))
                n_pairs += 1
                emitted[a] = emitted[b] = True
                stack.append((a, 0))
                stack.append((b, 0))
            else:
                stack.pop()
    assert all(emitted), "cycle left in the dependency graph"
    assert stats["depth_increments"] <= n * max_depth
    assert stats["depth_decrements"] <= stats["depth_increments"]
    assert stats["height_increments"] <= n * max_depth
    assert max(depth) <= max_depth and max_stack <= 2 * max_depth + 1
    assert n_single + 1 + 2 * n_pairs == n
    stats.update({"n": n, "edges": n_edges, "pairs": n_pairs, "roots": n_single + 1,
                  "pair_fraction": round(2 * n_pairs / n, 5), "max_depth": max(depth),
                  "max_stack": max_stack, "center_priority": center_priority,
                  "augmenting_paths": False,
                  "used_centers": sum(bool(x) for x in children),
                  "center_fraction": sum(bool(x) for x in children) / n,
                  "max_pairs_per_center": max(map(len, children)),
                  "depth_update_bound_per_vertex": 2 * max_depth})
    return ops, stats


# ------------------------------------------------------------------ curve
class Curve:
    def __init__(self, p, c):
        self.p, self.c = p, c
        self.inv = self.muls = 0

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

    def single(self, P, T):
        """Root step P+T: one reference addition (6M convention)."""
        assert P is not None and T is not None and P[0] != T[0]
        self.inv += 1
        self.muls += 6
        return self.add(P, T)

    def pair(self, P, T):
        """P+T and P-T from ONE inversion of (xT - xP): 6M + 3M extra."""
        p = self.p
        (x1, y1), (x2, y2) = P, T
        assert x1 != x2
        inv = pow(x2 - x1, -1, p)
        self.inv += 1
        self.muls += 9
        lp = (y2 - y1) * inv % p
        lm = (-y2 - y1) * inv % p
        x3p = (lp * lp - x1 - x2) % p
        x3m = (lm * lm - x1 - x2) % p
        return ((x3p, (lp * (x1 - x3p) - y1) % p),
                (x3m, (lm * (x1 - x3m) - y1) % p))

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R


def find_curve(bits, dd, rng):
    """p = 1 mod 3, y^2 = x^3 + c; prime order r with 6d | r-1 and 2 not in H."""
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
            if len(cands) == 1 and is_prime(cands[0]) and (cands[0] - 1) % (6 * dd) == 0 \
                    and pow(2, (cands[0] - 1) // dd, cands[0]) != 1:
                return E, cands[0], P


def dth_root(t, dd, rr):
    k = dd.bit_length() - 1
    for j in range(k, 0, -1):
        y = sqrt_mod(t, rr)
        e = (rr - 1) >> (j - 1)
        t = next(z for z in (y, rr - y) if j == 1 or pow(z, e, rr) == 1)
    return t


# ------------------------------------------------------------------ ec toy attack
def ec_attack(cfg):
    bits, dd, K = cfg["bits"], 1 << cfg["logd"], float(cfg["K"])
    trials, seed = cfg["trials"], cfg["seed"]
    KEYB, NODEB, HEADB = cfg["keybits"], cfg["nodebits"], cfg["headbits"]
    STOREB = KEYB - NODEB - HEADB
    rng = random.Random(seed)
    E, rr, P = find_curve(bits, dd, rng)
    p, mm = E.p, (rr - 1) // dd
    beta = next(z for z in (pow(g, (p - 1) // 3, p) for g in range(2, 100)) if z != 1)
    lamg = next(L for L in (pow(g, (rr - 1) // 3, rr) for g in range(2, 100))
                if L != 1 and E.mul(L, P) == (beta * P[0] % p, P[1]))
    U = [1, rr - 1, lamg, rr - lamg, lamg * lamg % rr, rr - lamg * lamg % rr]
    assert all(pow(u, mm, rr) == 1 for u in U)
    chi, gen = character(rr, dd)
    size = math.isqrt(int(K * mm / 6)) + 1
    X2, GMAX = 1 << cfg["logX2"], 1 << cfg["logG"]
    classes = [[] for _ in range(dd)]
    for b in range(1, X2 + 1):
        classes[chi(b)].append(b)
    plans = [plan(L, GMAX) if L else (None, None) for L in classes]
    pair_frac = (sum(2 * st["pairs"] for _, st in plans if st) /
                 sum(st["n"] for _, st in plans if st))
    X1 = X2 + int(X2 ** 0.921)
    while True:
        base_primes = primes_in(X2 + 1, X1)
        bases, cnt, baby_n, giant_n = [], 0, None, None
        for a in base_primes:
            bases.append(a)
            cnt += len(classes[(-chi(a)) % dd])
            if baby_n is None:
                if cnt >= size:
                    baby_n, cnt = len(bases), 0
            elif cnt >= size:
                giant_n = len(bases)
                break
        if baby_n and giant_n:
            break
        X1 *= 2
    assert baby_n and giant_n

    def expected(a, limit):
        vals, c = [], (-chi(a)) % dd
        L, (ops, _) = classes[c], plans[c]
        for op in ops:
            positions = (op[1],) if op[0] in (0, 1) else (op[3], op[2])
            for position in positions:
                vals.append(a * L[position])
                if len(vals) == limit:
                    return vals
        return vals

    def walk(base_point, a, emit, limit, checks):
        produced, c = 0, (-chi(a)) % dd
        L, (ops, _) = classes[c], plans[c]
        Qa = E.mul(a, base_point)
        tab = {}
        def T(g):
            if g not in tab:
                tab[g] = E.mul(g, Qa)
            return tab[g]
        pts, root = {}, None
        for op in ops:
            if op[0] == 0:
                root = E.mul(L[op[1]], Qa)
                pts[op[1]] = root
                emit(a * L[op[1]], root)
                produced += 1
            elif op[0] == 1:
                root = E.single(root, T(op[2]))
                pts[op[1]] = root
                emit(a * L[op[1]], root)
                produced += 1
            else:
                _, x, i, j, g = op
                Pp, Pm = E.pair(pts[x], T(g))
                pts[j], pts[i] = Pp, Pm
                emit(a * L[j], Pp)
                produced += 1
                if produced < limit:
                    emit(a * L[i], Pm)
                    produced += 1
            if checks[0] > 0 and produced % 41 == 0:
                v = a * L[op[1] if op[0] in (0, 1) else op[3]]
                assert E.mul(v, base_point) == pts[op[1] if op[0] in (0, 1) else op[3]]
                checks[0] -= 1
            if produced >= limit:
                return produced
        return produced

    NNODES, NHEADS = 1 << NODEB, 1 << HEADB
    HEAD_CAP, KEPT_CAP, CAND_CAP = cfg["head_cap"], 4, cfg["cand_cap"]
    mask = (1 << KEYB) - 1
    keyf = lambda Pt: min(Pt[1], p - Pt[1]) & mask
    split = lambda k: (k >> (KEYB - NODEB), (k >> STOREB) & (NHEADS - 1),
                       k & ((1 << STOREB) - 1))
    succ = inv_t = mul_t = cand_total = cmp_total = rw_mul_t = 0
    rewalk_fails = cand_over = 0
    occ_chi2_sum = 0.0
    spot = [cfg["spotchecks"]]
    for _ in range(trials):
        s = rng.randrange(1, rr)
        t_secret = pow(s, dd, rr)
        A1, Ad = E.mul(s, P), E.mul(t_secret, P)
        E.inv = E.muls = 0
        heads = {}
        candidates = []
        emitted_b, emitted_g = [], []

        def baby_emit(idv):
            def emit(v, Pt):
                k = keyf(Pt)
                emitted_b.append(v)
                node, head, st = split(k)
                h = heads.setdefault((node, head), [])
                h.append((st, idv))
                assert len(h) <= HEAD_CAP
            return emit

        giant_state = {}
        def giant_emit(v, Pt):
            k = keyf(Pt)
            emitted_g.append(v)
            node, head, st = split(k)
            chain = heads.get((node, head), ())
            for (st_b, idv) in chain:
                if st_b == st:
                    candidates.append((k, idv, giant_state["id"], v))
            giant_state["cmp"] += len(chain)

        remaining = size
        for bi in range(baby_n):
            if remaining <= 0:
                break
            remaining -= walk(P, bases[bi], baby_emit(bi), remaining, spot)
        assert len(emitted_b) == size
        # H3 support: occupancy chi-square over node x head cells
        lam_occ = size / (NNODES * NHEADS)
        cells = NNODES * NHEADS
        obs = [len(h) for h in heads.values()]
        chi2 = sum((o - lam_occ) ** 2 for o in obs)  # sum; scaled below
        # E[sum (o-lam)^2] = cells * lam * (1 - 1/cells) ~ size for binomial
        occ_chi2_sum += chi2 / (size * (1 - 1 / cells))
        giant_state["id"], giant_state["cmp"] = baby_n, 0
        remaining = size
        for gi in range(baby_n, giant_n):
            if remaining <= 0:
                break
            remaining -= walk(Ad, bases[gi], giant_emit, remaining, spot)
            giant_state["id"] = gi + 1
        assert len(emitted_g) == size
        walk_inv, walk_mul = E.inv, E.muls
        exp_b = []
        for bi in range(baby_n):
            exp_b += expected(bases[bi], size)
            if len(exp_b) >= size:
                break
        exp_b = exp_b[:size]
        exp_g = []
        for gi in range(baby_n, giant_n):
            exp_g += expected(bases[gi], size)
            if len(exp_g) >= size:
                break
        exp_g = exp_g[:size]
        assert exp_b == emitted_b and exp_g == emitted_g
        inverses = [pow(g, -1, rr) for g in exp_g]
        exact_ratios = {u * b * inv % rr for u in U for b in exp_b for inv in inverses}
        hit = t_secret in exact_ratios
        cand_total += len(candidates)
        if len(candidates) > CAND_CAP:
            cand_over += 1
        solved = False
        for (k, idb, idg, gv) in candidates[:CAND_CAP]:
            def collect(base_point, a, key):
                kept = []
                def emit(v, Pt):
                    if keyf(Pt) == key:
                        kept.append(v)
                walk(base_point, a, emit, size, [0])
                return kept
            keptb = collect(P, bases[idb], k)
            keptg = collect(Ad, bases[idg], k)
            if len(keptb) > KEPT_CAP or len(keptg) > KEPT_CAP:
                rewalk_fails += 1
                continue
            done = False
            for bv in keptb:
                for gvv in keptg:
                    for u in U:
                        z = u * bv * pow(gvv, -1, rr) % rr
                        if E.mul(z, P) == Ad:
                            s0 = dth_root(z, dd, rr)
                            w = pow(gen, (rr - 1) // dd, rr)
                            sol = next(s0 * pow(w, i, rr) % rr for i in range(dd)
                                       if E.mul(s0 * pow(w, i, rr) % rr, P) == A1)
                            assert sol == s, "wrong scalar"
                            done = True
                    if done:
                        break
                if done:
                    break
            if done:
                solved = True
                break
        assert solved == hit, f"resolution mismatch solved={solved} exact={hit}"
        cmp_total += giant_state["cmp"]
        inv_t += walk_inv
        mul_t += walk_mul
        rw_mul_t += E.muls - walk_mul
        if solved:
            succ += 1
    Ke = 6 * size * size / mm
    return {
        "p_bits": round(math.log2(p), 2), "r": rr, "d": dd, "X2": X2, "GMAX": GMAX,
        "GMAX_over_d2": round(GMAX / dd / dd, 1), "X2_over_GMAX": round(X2 / GMAX, 1),
        "bases_baby_giant": [baby_n, giant_n - baby_n], "size_per_side": size,
        "plan_pair_fraction": round(pair_frac, 4), "K_eff": round(Ke, 4),
        "predicted": round(1 - math.exp(-Ke), 4),
        "keybits": KEYB, "nodes": NNODES, "heads_per_node": NHEADS,
        "trials": trials, "successes": succ, "success_rate": round(succ / trials, 4),
        "coverage_events_match_exact_each_trial": True,
        "candidates_per_trial": round(cand_total / trials, 3),
        "false_match_model": round(size * size / (1 << KEYB), 3),
        "comparisons_per_trial": round(cmp_total / trials, 1),
        "comparison_model": round(size * size / (NNODES * NHEADS), 1),
        "occupancy_chi2_scaled_mean": round(occ_chi2_sum / trials, 4),
        "rewalk_cap_fails": rewalk_fails, "candidate_cap_over": cand_over,
        "all_recovered_scalars_correct": True,
        "mults_per_element": round(mul_t / trials / (2 * size), 3),
        "rewalk_mults_per_element_at_cand_rate": round(rw_mul_t / trials / (2 * size), 3),
        "reference_mults_per_element": 6.0,
    }


# ------------------------------------------------------------------ census
SMOOTH_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
                 67, 71, 73, 79, 83, 89, 97)

def census(cfg):
    m_bits, k = cfg["m_bits"], cfg["k"]
    K, seed = float(cfg["K"]), cfg["seed"]
    m_target = 1 << m_bits
    # deterministic toy prime: r = 6 c 2^k + 1, c odd and smooth (products of
    # primes <= 97), the (seed+1)-th smooth candidate >= m_target/6 whose r is
    # prime; skips r with 2 in H_t. Smooth c keeps the Pohlig-Hellman dlogs of
    # the sieve cheap; the census construction does not otherwise depend on c.
    c_lo = m_target // 6
    smooths = []
    def gen(idx, val):
        if idx == len(SMOOTH_PRIMES):
            smooths.append(val)
            return
        v = val
        while v <= 4 * c_lo:
            gen(idx + 1, v)
            v *= SMOOTH_PRIMES[idx]
    gen(0, 1)
    smooths = sorted({s for s in smooths if s >= c_lo})
    # deterministic pick over the sorted smooth candidates, skipping failures
    count = 0
    r = None
    for c in smooths:
        rr_try = 6 * c * (1 << k) + 1
        if is_prime(rr_try) and ((rr_try - 1) // (1 << k)) % 6 == 0 \
                and pow(2, (rr_try - 1) // (1 << k), rr_try) != 1:
            if count == seed:
                r = rr_try
                break
            count += 1
    assert r is not None, "no toy prime for this seed"
    d = 1 << k
    m = (r - 1) // d
    n_all = r - 1
    fac = {}
    for q in (2, 3):
        while n_all % q == 0:
            fac[q] = fac.get(q, 0) + 1
            n_all //= q
    rest, q = n_all, 5
    while q * q <= rest:
        while rest % q == 0:
            fac[q] = fac.get(q, 0) + 1
            rest //= q
        q += 2
    if rest > 1:
        fac[rest] = fac.get(rest, 0) + 1
    gen = None
    for cand in range(2, 100000):
        if all(pow(cand, (r - 1) // q, r) != 1 for q in fac) and pow(cand, r - 1, r) == 1:
            gen = cand
            break
    assert gen is not None
    parts = []
    for q, e in fac.items():
        pe = q ** e
        beta = pow(gen, (r - 1) // pe, r)
        gam = pow(beta, q ** (e - 1), r)
        b = int(math.isqrt(q)) + 1
        baby = {}
        cur = 1
        for j in range(b):
            baby[cur] = j
            cur = cur * gam % r
        gamminv = pow(pow(gam, b, r), -1, r)
        parts.append((q, e, pe, beta, gam, b, baby, gamminv))

    def dlog_q(t, b, baby, gamminv):
        cur = t
        for i in range(b + 1):
            j = baby.get(cur)
            if j is not None:
                return i * b + j
            cur = cur * gamminv % r
        raise ValueError("dlog_q")

    def dlog(x):
        crt, mods = [], []
        for q, e, pe, beta, gam, b, baby, gamminv in parts:
            xq = pow(x, (r - 1) // pe, r)
            digit = 0
            for i in range(e):
                t = pow(xq, q ** (e - 1 - i), r)
                di = dlog_q(t, b, baby, gamminv)
                digit += di * q ** i
                xq = xq * pow(beta, (pe - (di * q ** i) % pe) % pe, r) % r
            crt.append(digit)
            mods.append(pe)
        X, Nm = 0, 1
        for a, mo in zip(crt, mods):
            t = (a - X) % mo
            X = X + Nm * (t * pow(Nm, -1, mo) % mo)
            Nm *= mo
        return X % (r - 1)

    assert dlog(2) % d != 0, "2 in H_t"
    N = math.ceil(math.sqrt(K * m / 6))
    e_t = max(2, round(N ** 0.4843))
    X2 = d * e_t
    if cfg.get("window") == "top":
        # control: base window near r//X2 so that a*b ~ r (production-like
        # modular wraparound in the ratios); requires r//X2 > X2
        I0 = r // X2
        assert I0 > X2, "toy too small for the top-window control"
    else:
        I0 = X2
    Ilen = int(I0 ** 0.921)
    n_side = math.ceil(N / e_t)
    primes_I = primes_in(I0 + 1, I0 + Ilen)
    need = 2 * (n_side + 8) + 10
    while len(primes_I) < need:
        Ilen *= 2
        primes_I = primes_in(I0 + 1, I0 + Ilen)
    spf = list(range(X2 + 1))
    for q in range(2, int(math.isqrt(X2)) + 1):
        if spf[q] == q:
            for mm2 in range(q * q, X2 + 1, q):
                if spf[mm2] == mm2:
                    spf[mm2] = q
    Ep = {}
    for q in range(2, X2 + 1):
        if spf[q] == q:
            Ep[q] = dlog(q)
    modn = r - 1
    E = [0] * (X2 + 1)
    for b in range(2, X2 + 1):
        q = spf[b]
        E[b] = (E[b // q] + Ep[q]) % modn if q != b else Ep[q] % modn
    classes = {c2: [] for c2 in range(d)}
    for b in range(1, X2 + 1):
        classes[E[b] % d].append(b)
    modS = (r - 1) // 6
    M = m // 6
    def emitted(bs):
        vals, cnt, used = [], 0, 0
        for a in bs:
            if cnt >= N:
                break
            used += 1
            Ea = dlog(a)
            L = classes.get((-Ea) % d) or []
            take = min(len(L), N - cnt)
            vals.extend((Ea + E[b]) % modS for b in L[:take])
            cnt += take
        return vals, cnt, used
    sb, cb, nb = emitted(primes_I)
    sg, cg, ng = emitted(primes_I[nb:])
    assert cb == N and cg == N, f"emitted {cb},{cg} != {N}"
    if cfg.get("multiplicity"):
        # H5 diagnostic: per-coset multiplicity nu (number of emitted pairs whose
        # ratio lies in the coset = the candidate multiplicity for a uniform t there)
        nu = [0] * M
        for ti in sb:
            for uj in sg:
                nu[((ti - uj) % modS) // d] += 1
        marks = bytes(1 if v else 0 for v in nu)
        e1 = sum(nu)
        e2 = sum(v * v for v in nu)
        energy = e2 / (e1 + e1 * e1 / M)   # 1.0 for independent uniform pairs
        mult = {"energy_ratio": round(energy, 4),
                "mean_nu_marked": round(e1 / max(1, sum(marks)), 4),
                "max_nu": max(nu)}
    else:
        marks = bytearray(M)
        for ti in sb:
            for uj in sg:
                marks[((ti - uj) % modS) // d] = 1
        mult = None
    cov = sum(marks) / M
    Kr = 6 * N * N / m
    pred = 1 - math.exp(-Kr)
    repfrac = 1 - sum(marks) / (N * N)
    return {
        "r": r, "d": d, "m": m, "window": cfg.get("window", "low"), "K": round(Kr, 6),
        "N": N, "n_side": nb,
        "e_t": e_t, "X2": X2, "coverage": round(cov, 12), "predicted": round(pred, 12),
        "diff": round(cov - pred, 12), "repeat_fraction": round(repfrac, 6),
        "uniform_repeat_fraction": round(1 - (1 - math.exp(-Kr)) / Kr, 5),
        **({"multiplicity": mult} if mult else {}),
    }


# ------------------------------------------------------------------ planner mode
def planner_mode(cfg):
    out = []
    for spec in cfg["specs"]:
        dd, logX2, seed = spec["d"], spec["logX2"], spec["seed"]
        rng = random.Random(seed)
        rr = toy_prime(28 + int(math.log2(dd)) + 4, dd, rng)
        X2 = 1 << logX2
        cl = class_lists(rr, dd, X2)
        GMAX = 16 * dd * dd
        rows = []
        for c in range(min(spec.get("classes", 8), dd)):
            L = cl[c]
            if len(L) < 50:
                continue
            _, st = plan(L, GMAX)
            rows.append(st)
        if rows:
            out.append({
                "r": rr, "d": dd, "X2": X2, "GMAX": GMAX, "GMAX_over_d2": 16,
                "span_over_2GMAX": round(X2 / (2 * GMAX), 1),
                "classes": len(rows),
                "list_lengths": [min(st["n"] for st in rows), max(st["n"] for st in rows)],
                "min_pair_fraction": round(min(st["pair_fraction"] for st in rows), 4),
                "mean_pair_fraction": round(
                    sum(st["pair_fraction"] for st in rows) / len(rows), 4),
                "max_cycle_rejections": max(st["cycle_rejections"] for st in rows),
                "max_depth_seen": max(st["max_depth"] for st in rows),
                "depth_cap": 1 << 14,
            })
    return out


# ------------------------------------------------------------------ ledger
def ledger(_cfg):
    K = Fraction(51, 100)
    F = Fraction(98, 100)
    r, d, m = R_PRODUCTION, D_PRODUCTION, M_PRODUCTION
    # production algebra identities (asserted, not assumed)
    assert (r - 1) % d == 0 and (r - 1) % (1 << 32) == 0 and (r - 1) % (1 << 33) != 0
    assert m % 6 == 0
    x_param = (-0xd201000000010000) % r
    lam_p = (x_param * x_param - 1) % r
    assert lam_p != 1 and (lam_p * lam_p + lam_p + 1) % r == 0
    assert pow(lam_p, m, r) == 1, "U not in H"
    assert pow(2, m, r) != 1, "chi(2) == 0"
    N = math.isqrt(int(K * m // 6))
    while 6 * (N + 1) ** 2 <= int(K * m):
        N += 1
    if 6 * N * N < int(K * m):
        N += 1
    NB = 2 * N
    X2 = 1 << 76
    GMAX = 1 << 44
    ELL = X2 // d
    BASES = (NB // ELL) + 4          # both sides together: 2*(N/ELL) + margin
    pair_arith, root_arith = 2000 + 3 * 383, 2000
    # one Fermat inversion per round of 2^12 lane operations at a generating node
    # (2^12 lanes per node: every per-node structure is below 2^63 B, byte offsets < 2^63)
    inv_share = -(-232864 // (1 << 12))     # 57 (ceil of 56.75; costs round up)
    pair_bk = 8 + 8 + 12 + 12 + 72 + 14 + inv_share + 10 + 50
    root_bk = 130 - 4 + inv_share            # (was: inversion share 4 at 2^16 lanes)
    record = 50
    plan_stream = 15 * (1 - float(F) / 2)
    margin = Fraction(23, 20)                                   # 1.15 on non-arithmetic
    arith_pt = float(F) * pair_arith / 2 + (1 - float(F)) * root_arith
    bk_pt = float(F) * pair_bk / 2 + (1 - float(F)) * root_bk
    nonarith = (bk_pt + record + plan_stream) * float(margin)
    per_pt = arith_pt + nonarith        # keys are inside the pair/root bookkeeping;
    # gap tables are charged by their own rows below, not here
    NODES = 1 << 60                 # owner nodes: id = two 30-bit key fields, two hops
    mean_per_node = Fraction(N * N, 1 << 178)   # total comparisons 2^113.26 spread over 2^60
    cap_cmp = math.ceil(101 * float(mean_per_node) / 100)
    CAND = 1 << 44
    rows = {}
    rows["records (2N points)"] = math.ceil(NB * per_pt)
    rows["comparisons"] = NODES * cap_cmp * 7
    rows["owner init (head and count zeroing, rings, queues, credits)"] = \
        (1 << 118) * 9 + (1 << 60) + (1 << 55)
    rows["gap tables"] = BASES * GMAX * 2100
    rows["gap-table inversions"] = BASES * (GMAX // (1 << 12)) * 232864
    rows["table strand starts"] = BASES * (1 << 12) * 8192000
    rows["per-base setup"] = BASES * 3 * 8192000
    rows["base window sieve"] = ((1 << 70) + (1 << 43)) * 128
    rows["chi of base primes"] = int(1.25506 * (1 << 70) / math.log(1 << 76)) * 400000
    rows["chi of primes <= X2"] = int(1.25506 * X2 / math.log(X2)) * 400000
    rows["chi sieve + class lists"] = X2 * 200
    rows["planning"] = X2 * ((1 << 24) * 16 + (1 << 19) + 64)
    rows["plan storage write"] = X2 * 64
    rows["candidate forwarding"] = CAND * 128
    rows["candidate re-walks"] = CAND * (2 * ELL * 2100 + 2 * GMAX * 2100 +
                                         2 * (ELL // (1 << 12) + 16 + (1 << 12)) * 232864)
    rows["candidate z tests"] = CAND * (96 * 8192000 + 1788 * 383)
    # blocks of 128 records: baby and giant records traverse two hops each
    n_blocks = 2 * ((NB + 127) // 128)
    rows["wake and handler dispatch"] = n_blocks * (1 << 8)
    rows["status tree"] = 3 * ((1 << 60) + (1 << 20) + (1 << 30)) * (1 << 12)
    rows["Lemma 1 + io"] = ((1 << 20) + 1) * (1 << 24) * 2000 + (1 << 30)
    rows["sieving primes sieve"] = ((1 << 47) + (1 << 43)) * 128
    units = {k: Fraction(int(math.ceil(v)), 2000) for k, v in rows.items()}
    W = sum(units.values())
    assert all(u > 0 for u in units.values())
    # no-pair context: same budgets with every point a root
    nopair_pt = 2000 + (130 + record + plan_stream) * float(margin)
    W_nopair = W - units["records (2N points)"] + Fraction(int(math.ceil(NB * nopair_pt)), 2000)
    # preprocessing: baby half of the records row plus the challenge-independent rows
    pre_instr = (math.ceil(N * per_pt)
                 + rows["gap tables"] // 2 + rows["gap-table inversions"] // 2
                 + rows["table strand starts"] // 2 + rows["per-base setup"] // 2
                 + rows["chi of primes <= X2"] + rows["chi sieve + class lists"]
                 + rows["planning"] + rows["plan storage write"]
                 + rows["base window sieve"] + rows["chi of base primes"]
                 + rows["sieving primes sieve"] + rows["owner init (head and count zeroing, rings, queues, credits)"]
                 + 1 * ((1 << 60) + (1 << 20) + (1 << 30)) * (1 << 12))
    pre = Fraction(int(math.ceil(pre_instr)), 2000)
    assert pre < W
    # memory: every per-node byte offset is below 2^63 (asserted per structure)
    # owners: heads 2^58 x 8 B (max offset 2^61), counts 2^58 x 1 B, pool 2^57 x 32 B
    # (max offset 2^62), send arrays 2^30 x 512 B (2^39 B), rings 2^12 x 4 KiB, cand list
    heads_bytes, counts_bytes, pool_bytes = (1 << 58) * 8, (1 << 58), (1 << 57) * 32
    assert heads_bytes < (1 << 63) and counts_bytes < (1 << 63) and pool_bytes < (1 << 63)
    owner_mem = ((1 << 58) * 9 + (1 << 57) * 32 + (1 << 42) + (1 << 16) * 40 + (1 << 30)
                 # the 2^30 of the owners that are hop-1 relays also hold: per-generating-node
                 # rings (2^20 x 8 x 4 KiB), ready queues, and per-owner return rings (2^30 x 8 x 4 KiB)
                 + (1 << 20) * 8 * 4096 + (1 << 25) + (1 << 30) * 8 * 4096)
    assert owner_mem < (1 << 63)
    # generating: 2^12 lanes x (table 2^44 x 96 B = 2^50.58 B + stack) + send + receive
    table_bytes = GMAX * 96        # one lane's gap table: 2^44 entries x 96 B = 2^50.6 B
    assert table_bytes < (1 << 63)
    gen_mem = ((1 << 12) * (table_bytes + (1 << 15) * 96 + (1 << 16) * 32)
               + (1 << 42)   # send staging: one 4 KiB block per destination
               + (1 << 30) * 4096 + (1 << 30) * 4096   # per-relay and per-owner rings
               + (1 << 24) * 40 + (1 << 30))
    assert gen_mem < (1 << 63)
    stor_total = X2 * 10 + int(0.51 * X2) * 32 + (1 << 41) * 8 + (1 << 30)
    stor_per = -(-stor_total // (1 << 30))     # ceil share per storage node
    assert stor_per < (1 << 63)
    mem = (1 << 60) * owner_mem + (1 << 20) * gen_mem + (1 << 30) * stor_per + (X2 * 5) // 2
    Keff = Fraction(6 * N * N, m)
    assert Keff >= K
    # H5 (multiplicity) now carries the candidate-cap tails; the union budget is 4*2^-40
    succ_floor = (1 - math.exp(-float(Keff))
                  - 0.0005 - 2 ** -20 - 4 * 2 ** -40
                  - 6 * (BASES // 2) ** 2 / m)
    assert succ_floor >= 0.398, succ_floor
    import decimal
    decimal.getcontext().prec = 80
    WWd = decimal.Decimal(W.numerator) / decimal.Decimal(W.denominator)
    score_digits = WWd.ln() / decimal.Decimal(2).ln()          # log2 W, 80 digits
    score_claim = math.ceil(score_digits * 10 ** 6)
    pred = decimal.Decimal(pre.numerator) / decimal.Decimal(pre.denominator)
    pre_claim = math.ceil(pred.ln() / decimal.Decimal(2).ln() * 10 ** 4)
    mem_claim = math.ceil(math.log2(mem) * 10 ** 4)
    # exact integer assert: W^(10^6) <= 2^(score_claim * 10^6 / ...)
    sc = score_claim  # = ceil(10^6 log2 W); the ceil is exact because log2 W
    # is more than 10^-20 away from the next 10^-6 grid point (asserted; the
    # 80-digit decimal computation is exact to ~10^-60 here)
    val6 = score_digits * 10 ** 6
    frac6 = val6 - val6.to_integral_value(rounding=decimal.ROUND_DOWN)
    assert (1 - frac6) > 10 ** -20 or frac6 == 0, "boundary too close for ceil"
    return {
        "N": N, "log2_N": round(math.log2(N), 6),
        "K_effective": round(float(Keff), 8), "bases_charged": BASES,
        "per_point_instructions": round(per_pt, 2),
        "nopair_W_log2": round(math.log2(float(W_nopair)), 6),
        "rows_log2": {k: round(math.log2(float(u)), 3) for k, u in
                      sorted(units.items(), key=lambda kv: -float(kv[1]))},
        "W_log2": str(score_digits)[:20],
        "score_ceil_1e6": score_claim / 10 ** 6,
        "preprocessing_log2_ceil_1e4": pre_claim / 10 ** 4,
        "memory_log2_ceil_1e4": mem_claim / 10 ** 4,
        "success_floor": round(succ_floor, 6), "success_claimed": 0.398,
        "score_integer_assert": True, "preprocessing_le_W_assert": True,
    }


# ------------------------------------------------------------------ main
def main():
    cfg = json.loads(sys.stdin.read())
    mode = cfg.get("mode", "all")
    t0 = time.time()
    out = {"mode": mode}
    if mode in ("ledger", "all"):
        out["ledger"] = ledger(cfg)
    if mode in ("planner", "all"):
        out["planner"] = planner_mode(cfg.get("planner", {
            "specs": [{"d": 16, "logX2": 15, "seed": 0, "classes": 8}]}))
    if mode in ("ec", "all"):
        out["ec"] = ec_attack(cfg.get("ec", {
            "bits": 23, "logd": 5, "K": 0.51, "trials": 60, "seed": 2026,
            "keybits": 18, "nodebits": 5, "headbits": 6, "head_cap": 64,
            "cand_cap": 4096, "logX2": 13, "logG": 12, "spotchecks": 6}))
    if mode in ("census", "all"):
        out["census"] = [census(dict(cfg.get("census",
                {"m_bits": 21, "k": 4, "K": 0.51, "seed": 0})))]
    out["seconds"] = round(time.time() - t0, 2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3 -I
"""Evidence for the paired class-walk Cheon attack on bls12-381-g1-ecdlp
(synchronous wave transport, measured side-output kernel).

Modes (input object on stdin: {"mode": "...", ...}):
  ledger  -- every cost term of proof.md from exact integers, the exact score
             rounded up, preprocessing and memory rounded up, success rounded
             down; asserts the inequalities. Arithmetic only.
  planner -- the offline pairing planner on real chi-class lists of toy primes
             at the production degree parameter GMAX = 16 d^2 (evidence for H4).
  ec      -- the whole algorithm on a toy curve: paired walk, bucket matching
             with truncated keys, false matches, bounded re-walk, recovery.
  proto   -- the section-4 wave transport at toy scale: per-link receive
             buffers with caps, double buffering, global barriers, exact
             delivery, candidate and reply paths, instruction counts by the
             listing, closed form equal to the ledger's formula.
  census  -- exact coverage |R|/m of the construction on toy primes (H1).
  all     -- ledger + planner + ec + proto + census at declared hosted sizes.

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
def plan(L, GMAX, max_depth=1 << 14, center_priority=False, max_tree=1 << 20):
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
    stats = {"cycle_rejections": 0, "depth_rejections": 0, "demoted": 0, "size_rejections": 0,
             "depth_increments": 0, "depth_decrements": 0, "height_increments": 0}
    size = [1] * n          # tree size, maintained at the root

    def root_of(v):
        while parent[v] != -1:
            v = parent[v]
        return v

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
            rk = root_of(k)
            if size[rk] + size[i] + size[j] > max_tree:
                stats["size_rejections"] += 1
                return False
            size[rk] += size[i] + size[j]
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
    tree_sizes = {}
    for v in range(n):
        rv = v
        while parent[rv] != -1:
            rv = parent[rv]
        tree_sizes[rv] = tree_sizes.get(rv, 0) + 1
    assert max(tree_sizes.values()) <= max_tree
    stats.update({"n": n, "edges": n_edges, "pairs": n_pairs, "roots": n_single + 1,
                  "max_tree_size": max(tree_sizes.values()), "tree_cap": max_tree,
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
                "depth_cap": 1 << 14, "max_tree_size_seen": max(st["max_tree_size"] for st in rows),
                "size_rejections": sum(st["size_rejections"] for st in rows), "tree_cap": 1 << 20,
            })
    return out


# ------------------------------------------------------------------ shared transport formulas
def barrier_tree(nodes, fanin):
    """Deterministic fan-in tree over `nodes` leaves: returns (internal nodes,
    sum over internal nodes of children^2). Used identically by proto and ledger."""
    internal = sq = 0
    level = nodes
    while level > 1:
        parents = -(-level // fanin)
        full, rest = divmod(level, fanin)
        sq += full * fanin * fanin + (rest * rest if rest else 0)
        internal += parents
        level = parents
    return internal, sq

WAKE = 256            # interrupt entry, handler dispatch, return (charged per wake)
def barrier_instr(nodes, fanin, count):
    """Instructions of `count` global barriers: every non-root node arrives once
    (dmb, arrival store, dmb, doorbell store = 4); every arrival wakes its parent
    once (worst case, no coalescing) and the parent sweeps all its child slots
    (2 per slot: load, compare); every non-root node is released by one wake and
    one flag load + dmb."""
    internal, sq = barrier_tree(nodes, fanin)
    arrivals = nodes + internal - 1
    return count * (arrivals * 4 + arrivals * WAKE + sq * 2 + arrivals * (WAKE + 2)), internal

LINK_WAVE = 5         # per link per wave: sender count store + counter reset, receiver count load + loop setup
REC_GEN, REC_RELAY, REC_INS, REC_LOOK = 14, 17, 17, 8   # listed per record (proof.md section 5)


# ------------------------------------------------------------------ proto (synchronous wave transport)
def proto(cfg):
    from collections import Counter, defaultdict
    rng = random.Random(cfg.get("seed", 2026))
    GB, SB, TB, UB = cfg.get("gen_bits", 3), cfg.get("s_bits", 1), cfg.get("t_bits", 2), cfg.get("u_bits", 2)
    LANES, WAVES, STEPS = cfg.get("lanes", 4), cfg.get("waves", 3), cfg.get("steps", 6)
    CAP1, CAP2, CAPC, FANIN = cfg.get("cap1", 96), cfg.get("cap2", 96), cfg.get("capc", 8), cfg.get("fanin", 4)
    pattern = cfg.get("pattern", "uniform")
    assert GB >= SB
    G, R, O = 1 << GB, 1 << (TB + SB), 1 << (TB + UB)
    NODES = G + R + O
    cost = Counter()
    buf = {}                       # (sender, receiver, parity) -> list of records
    tag = {}                       # (sender, receiver, parity) -> (written_wave, read_wave)
    cnt = {}                       # count words
    emitted, delivered = [], defaultdict(list)
    heads, candidates = defaultdict(list), []
    overflow = []

    def open_buf(sender, receiver, parity, wave):
        key = (sender, receiver, parity)
        lst = buf.setdefault(key, [])
        w, rd = tag.get(key, (None, None))
        if w != wave:              # first use in this wave: previous contents were read in wave w+1
            assert w is None or rd == w + 1, "buffer reused before it was read"
            lst.clear(); tag[key] = (wave, None)
        return lst

    def put(sender, receiver, parity, rec, cap, wave):
        lst = open_buf(sender, receiver, parity, wave)
        if len(lst) >= cap:
            overflow.append((sender, receiver, wave)); return False
        lst.append(rec); return True

    def take(sender, receiver, parity, wave):
        key = (sender, receiver, parity)
        w, rd = tag.get(key, (None, None))
        if w is None: return []
        assert w == wave - 1 and rd is None, "read in the wrong wave"
        tag[key] = (w, wave)
        return list(buf[key])

    def relay_of(key, g): return ((key & (R // (1 << SB) - 1)) << SB) | (g & ((1 << SB) - 1))
    def owner_of(key): return key & (O - 1)           # t = low TB bits, u = next UB bits
    def t_of(x): return x & ((1 << TB) - 1)

    def barrier():
        instr, _ = barrier_instr(NODES, FANIN, 1)
        cost["barrier"] += instr

    wave_count = Counter()
    seq = 0
    for phase in (0, 1):                               # 0 = baby, 1 = giant
        for wl in range(WAVES + 2):
            w = phase * (WAVES + 2) + wl               # global wave index (buffer tags), wl = index within the phase
            # generating nodes: emit wave w
            if wl < WAVES:
                for g in range(G):
                    for lane in range(LANES):
                        for st in range(STEPS):
                            n_out = 2 if rng.random() < 0.98 else 1
                            for _ in range(n_out):
                                if pattern == "one_owner": key = 0 | (rng.getrandbits(180) << 12)
                                elif pattern == "two_owners": key = (rng.randrange(2) * (O // 2)) | (rng.getrandbits(180) << 12)
                                else: key = rng.getrandbits(192)
                                rec = (key, phase, g, lane, seq); seq += 1
                                emitted.append(rec)
                                cost["rec_gen"] += 1
                                if not put(g, G + relay_of(key, g), w & 1, rec, CAP1, w):
                                    return {"pattern": pattern, "failed": True, "reason": "link cap overflow (run fails)",
                                            "overflow_at": overflow[0], "records_before_failure": len(emitted)}
                    for rl in range(g & ((1 << SB) - 1), R, 1 << SB):
                        open_buf(g, G + rl, w & 1, w)      # count word written for every link, empty or not
                    cost["link"] += R >> SB            # count words to its 2^TB relays
                    wave_count["gen"] += 1
            # relays: forward wave w-1
            if 1 <= wl <= WAVES:
                for rl in range(R):
                    s = rl & ((1 << SB) - 1)
                    for g in range(s, G, 1 << SB):
                        for rec in take(g, G + rl, (w - 1) & 1, w):
                            cost["rec_relay"] += 1
                            assert relay_of(rec[0], g) == rl
                            if not put(G + rl, G + R + owner_of(rec[0]), (w - 1) & 1, rec, CAP2, w - 1):
                                return {"pattern": pattern, "failed": True, "reason": "link cap overflow (run fails)",
                                        "overflow_at": overflow[0], "records_before_failure": len(emitted)}
                    for o in range(O):
                        if t_of(o) == (rl >> SB): open_buf(G + rl, G + R + o, (w - 1) & 1, w - 1)
                    cost["link"] += (G >> SB) + O // (1 << TB)   # inbound counts read, outbound counts written
                    wave_count["relay"] += 1
            # owners: consume wave w-2
            if 2 <= wl <= WAVES + 1:
                for o in range(O):
                    t = t_of(o)
                    for s in range(1 << SB):
                        rl = (t << SB) | s
                        for rec in take(G + rl, G + R + o, w & 1, w - 1):
                            assert owner_of(rec[0]) == o
                            delivered[o].append(rec)
                            slot = (rec[0] >> (TB + UB)) & 63
                            rest = rec[0] >> (TB + UB + 6)
                            if phase == 0:
                                heads[(o, slot)].append((rest, rec)); cost["insert"] += 1
                            else:
                                cost["lookup"] += 1
                                for (rb, bb) in heads[(o, slot)]:
                                    cost["compare"] += 1
                                    if rb == rest: candidates.append((o, rec[0], bb, rec))
                    cost["link"] += 1 << SB
                    wave_count["owner"] += 1
            barrier(); cost["barriers"] += 1
    # exact delivery
    got = [rec for o in delivered for rec in delivered[o]]
    assert Counter(got) == Counter(emitted) and len(got) == len(emitted)
    n_b = sum(1 for r in emitted if r[1] == 0); n_g = len(emitted) - n_b
    assert cost["rec_gen"] == cost["rec_relay"] == len(emitted) and cost["insert"] == n_b and cost["lookup"] == n_g
    # candidate path: owner -> relay R(t, s_gen) -> generating node, both chains; then replies back along the forward path
    cand_cost = 0; resolved = 0
    # (buffer tags: a relay opens its output buffers with wave w-1, the wave whose records it forwards, and the
    #  owner takes them with wave w-1 one global wave later, so the tag check is write at w', read at w'+1, rewrite at w'+2)
    for (o, key, bb, gg) in candidates:
        t = t_of(o); replies = []
        for side in (bb, gg):
            g = side[2]; rl = (t << SB) | (g & ((1 << SB) - 1))
            cand_cost += 2 * REC_RELAY                 # owner -> relay, relay -> gen (two hops, charged like records)
            assert side in emitted                     # the generating node finds the chain (re-walk, charged in the ledger)
            # reply: gen -> relay R(t_key, s) -> owner(key) along the forward links
            rl2 = relay_of(key, g); assert owner_of(key) == o
            cand_cost += 2 * REC_RELAY; replies.append(side[4])
        resolved += 1
    assert resolved == len(candidates)
    links = G * (R >> SB) + R * (O // (1 << TB))
    expected_link_events = (R >> SB) * wave_count["gen"] + ((G >> SB) + O // (1 << TB)) * wave_count["relay"] + (1 << SB) * wave_count["owner"]
    assert cost["link"] == expected_link_events
    internal, _ = barrier_tree(NODES, FANIN)
    total = (cost["rec_gen"] * REC_GEN + cost["rec_relay"] * REC_RELAY + cost["insert"] * REC_INS + cost["lookup"] * REC_LOOK
             + cost["compare"] * 7 + cost["link"] * LINK_WAVE + cost["barrier"] + cand_cost)
    formula = (len(emitted) * (REC_GEN + REC_RELAY) + n_b * REC_INS + n_g * REC_LOOK + cost["compare"] * 7
               + expected_link_events * LINK_WAVE + barrier_instr(NODES, FANIN, cost["barriers"])[0] + cand_cost)
    assert total == formula
    return {"pattern": pattern, "failed": False, "gens": G, "relays": R, "owners": O, "lanes": LANES, "waves_per_phase": WAVES,
            "records": len(emitted), "baby": n_b, "giant": n_g, "delivered_exactly_once": True,
            "every_buffer_read_in_the_wave_after_it_was_written": True, "barriers": cost["barriers"],
            "barrier_tree_internal_nodes": internal, "candidates": len(candidates), "candidates_resolved": resolved,
            "compares": cost["compare"], "links": links, "link_wave_events": cost["link"],
            "instructions_by_item": {"gen": cost["rec_gen"] * REC_GEN, "relay": cost["rec_relay"] * REC_RELAY,
                                     "insert": cost["insert"] * REC_INS, "lookup": cost["lookup"] * REC_LOOK,
                                     "compare": cost["compare"] * 7, "link_words": cost["link"] * LINK_WAVE,
                                     "barriers": cost["barrier"], "candidate_path": cand_cost},
            "instructions_total": total, "closed_form_equal": True,
            "scope": "sequential simulation; instruction counts are the listing applied to counted events, not hardware measurements"}


def proto_mode(cfg):
    out = []
    c = dict(cfg); c["pattern"] = "uniform"; out.append(proto(c))
    c = dict(cfg); c.update({"pattern": "two_owners", "cap2": 8 * cfg.get("cap2", 96)}); out.append(proto(c))   # skewed keys, wide caps
    c = dict(cfg); c.update({"pattern": "one_owner"}); out.append(proto(c))   # all keys to one owner: cap overflow -> run fails
    assert out[2]["failed"] and not out[0]["failed"] and not out[1]["failed"]
    return out


# ------------------------------------------------------------------ ledger
KERNEL = 1127          # measured side-output kernel, instructions per lane (measure/)

def ledger(cfg):
    from decimal import Decimal, localcontext, ROUND_CEILING
    ceil = lambda x: -(-x.numerator // x.denominator) if isinstance(x, Fraction) else x
    r, d, m = R_PRODUCTION, D_PRODUCTION, M_PRODUCTION
    assert (r - 1) % d == 0 and m % 6 == 0
    v2, t_ = 0, r - 1
    while t_ % 2 == 0: t_ //= 2; v2 += 1
    assert v2 == 32
    lam = (-0xd201000000010000) ** 2 - 1
    assert (lam * lam + lam + 1) % r == 0 and pow(lam, m, r) == 1 and pow(2, m, r) != 1
    K, F = Fraction(1, 2), Fraction(49, 50)
    N = math.isqrt(m // 12)
    while 12 * N * N < m: N += 1
    assert 12 * (N - 1) ** 2 < m <= 12 * N * N and 6 * N * N >= K * m
    X2, GMAX, EMIN, EMAX = 1 << 76, 1 << 44, 1 << 55, 1 << 57
    LANES, SEG, TREEMAX, WSTEP, DEPTH = 1 << 16, 1 << 25, 1 << 20, 1 << 36, 1 << 14
    GENS, SBITS, TBITS, UBITS = 1 << 47, 20, 30, 30
    RELAYS, OWNERS = 1 << (TBITS + SBITS), 1 << (TBITS + UBITS)
    HEADS, POOL, CAP1, CAP2, CAPC, CAND = 1 << 57, 1 << 56, 1 << 24, 1 << 21, 1 << 8, 1 << 44
    RECB, FANIN = 40, 1 << 10
    MARG = Fraction(21, 20)
    SM = 8192000                                          # one scalar multiplication, instructions
    inv_share = ceil(Fraction(232864, LANES))             # batched inversion share per lane op
    PAIR_BK = 8 + 8 + 12 + 12 + 7 + inv_share + 10 + 25   # listed, section 5
    ROOT_BK = 5 + 2 + 12 + 14 + 6 + inv_share + 5 + 25
    TABLE_BK = 8 + 7 + 6 + inv_share + 5
    REC_AVG = Fraction(2 * (REC_GEN + REC_RELAY) + REC_INS + REC_LOOK, 2)
    PAIR = 2000 + KERNEL + MARG * PAIR_BK                 # two points
    ROOT = 2000 + MARG * ROOT_BK
    per_point = F * PAIR / 2 + (1 - F) * ROOT + MARG * REC_AVG
    # sizes (H2 bounds, checked by the sieves at run time; a violation stops the run = failure)
    BASES_SIDE = ceil(Fraction(N, EMIN)) + 1
    BASES = 2 * BASES_SIDE
    BCLASS = 1 << 42                                      # bases of one character, both sides
    GPC = GENS // d                                       # generating nodes per class
    BPG = ceil(Fraction(BCLASS, GPC))                     # bases per generating node
    SEGLEN = SEG + TREEMAX                                # entries per plan segment, at most
    LANE_OPS = (EMAX // SEG // LANES + 1) * SEGLEN         # ops per lane per base (stride assignment of segments)
    WAVES_SIDE = ceil(Fraction(BPG * LANE_OPS, WSTEP)) + 2
    PLAN_ROUNDS = ceil(Fraction(EMAX * 16, 1 << 48))
    WAVES = 2 * WAVES_SIDE + PLAN_ROUNDS + 64
    PADOPS = BASES * LANES * SEGLEN + 2 * LANES * LANE_OPS   # lanes differ by at most one segment; partial last bases
    ALLOW = 8                                             # explicit allowance per point for unlisted bookkeeping
    rows = {}
    rows["records: 2N points"] = ceil(2 * N * per_point)
    rows["padding ops at the pair price"] = ceil(PADOPS * PAIR)
    rows["listing allowance: 8 instructions per point"] = 2 * N * ALLOW
    rows["segment restarts (scalar multiplications)"] = (ceil(Fraction(2 * N, SEG)) + BASES * LANES) * SM
    rows["gap tables"] = BASES * (GMAX - 1) * (2000 + ceil(MARG * TABLE_BK))
    rows["table strand starts and base points"] = BASES * (LANES + 1) * SM
    cmp_cap = ceil(Fraction(101 * N * N, 100 * HEADS * OWNERS * OWNERS))
    rows["comparisons at the per-owner cap"] = OWNERS * cmp_cap * 7
    rows["owner initialisation"] = OWNERS * (HEADS * 8 // 16 + 64)
    LINKS1, LINKS2 = GENS * (1 << TBITS), RELAYS * (1 << UBITS)
    COPY = d * ((1 << 6) + (1 << 16) + (1 << 26))
    NODES = OWNERS + RELAYS + GENS + 2 * d + COPY + (1 << 14)     # 2^14 spare (input, root-recovery and base-assignment nodes)
    bar, BNODES = barrier_instr(NODES, FANIN, WAVES)
    rows["per-wave link words, memory barriers"] = WAVES * ((LINKS1 + LINKS2) * LINK_WAVE + NODES * 2)
    rows["global barriers (wakes, sweeps, releases)"] = bar
    rows["plan distribution (copy tree)"] = d * EMAX * ((1 << 6) + (1 << 16) + (1 << 26) + (1 << 27)) * 4
    WIN = 2 * (2 * GMAX // d)
    rows["planner (scans, matching, emission)"] = d * EMAX * (3 * WIN * 16 + 64 * (2 * DEPTH * 3 + 20) + 64 + 16)
    rows["character sieve, class lists, their transfer"] = X2 * (400000 + 16 + 20) + d * (1 << 33) * 400000
    rows["base interval sieve, characters, assignment"] = ((1 << 70) + (1 << 43)) * 16 + (1 << 65) * (400000 + 64) + BASES * 64
    rows["candidate delivery and replies"] = CAND * (8 * REC_RELAY + 200)
    rows["candidate re-walks"] = CAND * 2 * (SEGLEN * (3300 + 232864) + (GMAX - 1) * (2000 + 32) + 2 * SM)
    rows["candidate tests"] = CAND * (96 * SM + (1 << 16) * 20)
    rows["root recovery and input distribution"] = (1 << 24) * 2000 + NODES * 64 + (1 << 40)
    total = sum(rows.values())
    W = Fraction(total, 2000)
    pre = (ceil(Fraction(rows["records: 2N points"], 2)) + ceil(Fraction(rows["padding ops at the pair price"], 2))
           + ceil(Fraction(rows["listing allowance: 8 instructions per point"], 2))
           + ceil(Fraction(rows["segment restarts (scalar multiplications)"], 2)) + ceil(Fraction(rows["gap tables"], 2))
           + ceil(Fraction(rows["table strand starts and base points"], 2)) + rows["owner initialisation"]
           + ceil(Fraction(rows["per-wave link words, memory barriers"], 2)) + ceil(Fraction(rows["global barriers (wakes, sweeps, releases)"], 2))
           + rows["plan distribution (copy tree)"] + rows["planner (scans, matching, emission)"]
           + rows["character sieve, class lists, their transfer"] + rows["base interval sieve, characters, assignment"])
    assert pre < total
    # memory, bytes, every per-node structure and total below 2^63
    gen_mem = EMAX * 16 + GMAX * 96 + LANES * (2 * DEPTH + 1) * 96 + (1 << TBITS) * (8 + 2 * CAP1 * RECB) + (1 << TBITS) * CAPC * 64 + BPG * 16 + (1 << 30)   # windows for both parities
    relay_mem = (GENS >> SBITS) * 2 * CAP1 * RECB + (1 << UBITS) * (8 + 2 * CAP2 * RECB) + (1 << UBITS) * CAPC * 64 + (GENS >> SBITS) * CAPC * 64 + (1 << 30)
    owner_mem = HEADS * 8 + POOL * 32 + (1 << SBITS) * 2 * CAP2 * RECB + (1 << 16) * 64 + (1 << SBITS) * CAPC * 64 + (1 << 30)
    plan_mem = EMAX * 24 + (1 << 58) * 8 + (1 << 6) * (1 << 48) * 2 + (1 << 30)
    sieve_mem = (1 << 56) * 4 + (1 << 33) * 4 + d * (1 << 40) + (1 << 30)
    copy_mem = 2 * (1 << 48) + (1 << 10) * (1 << 48) + (1 << 30)
    for v in (gen_mem, relay_mem, owner_mem, plan_mem, sieve_mem, copy_mem):
        assert v < (1 << 63)
    mem = GENS * gen_mem + RELAYS * relay_mem + OWNERS * owner_mem + d * plan_mem + d * sieve_mem + COPY * copy_mem + (BNODES + (1 << 14)) * (1 << 20)
    # success: 1 - e^-K minus the H1 allowance, the modelling budget, the cap tails, the diagonal event
    with localcontext() as ctx:
        ctx.prec = 90
        log2 = lambda x: (Decimal(x.numerator) / Decimal(x.denominator)).ln() / Decimal(2).ln()
        s = log2(W)
        grid = (s * 10 ** 6).to_integral_value(rounding=ROUND_CEILING)
        assert grid - s * 10 ** 6 > Decimal(10) ** -14
        score = int(grid)
        bound = lambda x, scale: int((log2(x) * scale).to_integral_value(rounding=ROUND_CEILING)) / scale
        succ = 1 - math.exp(-0.5) - 0.0005 - 2 ** -20 - 4 * 2 ** -40 - 2 ** -111
        assert succ >= 0.392
        return {"K": str(K), "N": N, "log2_2N": float(log2(Fraction(2 * N))), "f": str(F), "kernel": KERNEL,
                "pair_op_instructions_two_points": str(PAIR), "root_op_instructions": str(ROOT),
                "record_instructions_listed_avg": str(REC_AVG), "per_point_instructions": str(per_point),
                "per_point_float": float(per_point), "allowance_per_point": ALLOW, "lane_ops": LANE_OPS, "bases": BASES, "waves": WAVES, "nodes_log2": math.log2(NODES),
                "barrier_internal_nodes": BNODES, "rows_log2_units": {k: float(log2(Fraction(v, 2000))) for k, v in rows.items()},
                "W_instructions": total, "W_log2": str(s), "score_numerator": score, "score_denominator": 10 ** 6,
                "score": score / 10 ** 6, "preprocessing_log2_ceil_1e4": bound(Fraction(pre, 2000), 10 ** 4),
                "memory_log2_ceil_1e4": bound(Fraction(mem), 10 ** 4),
                "per_node_bytes_log2": {"gen": math.log2(gen_mem), "relay": math.log2(relay_mem), "owner": math.log2(owner_mem),
                                        "planner": math.log2(plan_mem), "sieve": math.log2(sieve_mem), "copy": math.log2(copy_mem)},
                "success_floor": succ, "success_claimed": 0.392, "all_per_node_structures_below_2pow63": True}


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
            "bits": 23, "logd": 5, "K": "0.50", "trials": 60, "seed": 2026,
            "keybits": 18, "nodebits": 5, "headbits": 6, "head_cap": 64,
            "cand_cap": 4096, "logX2": 13, "logG": 12, "spotchecks": 6}))
    if mode in ("proto", "all"):
        out["proto"] = proto_mode(cfg.get("proto", {}))
    if mode in ("census", "all"):
        out["census"] = [census(dict(cfg.get("census",
                {"m_bits": 21, "k": 4, "K": "0.50", "seed": 0, "multiplicity": True})))]
    out["seconds"] = round(time.time() - t0, 2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

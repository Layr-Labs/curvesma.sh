#!/usr/bin/env python3 -I
"""Stratified-serial one-sided abort at toy scale (proof.md 8.5): baby full;
giant chunks in L greedy strata, serial rotated service, per-stratum re-walk
verification, real stop at first verified-true. Finite-shape only."""
import math
import random

from evidence import Curve, character, dth_root, find_curve, plan, primes_in


def walk_full(base_point, E, chi, dd, classes, plans, a, emit, limit):
    """Plain sequential walk of base a (local copy of the evidence walk)."""
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
        if produced >= limit:
            return produced
    return produced


def partition(ops, seg_t):
    """Plan ops into restart-segments: whole trees to >= seg_t entries each."""
    trees, cur = [], []
    for op in ops:
        if op[0] in (0, 1) and cur:
            trees.append(cur)
            cur = []
        cur.append(op)
    if cur:
        trees.append(cur)
    chunks, acc, acc_n = [], [], 0
    for t in trees:
        acc.append(t)
        acc_n += 1 + 2 * sum(1 for op in t if op[0] == 2)
        if acc_n >= seg_t:
            chunks.append([op for t in acc for op in t])
            acc, acc_n = [], 0
    if acc:
        chunks.append([op for t in acc for op in t])
    return chunks


def walk_chunk(base_point, E, chi, dd, classes, a, Qa, tab, chunk, emit, limit):
    """Walk one chunk from a fresh restart (scalar mult; fresh midpoint table)."""
    produced, restarts, c = 0, 0, (-chi(a)) % dd
    L = classes[c]
    pts, root, first = {}, None, True
    for op in chunk:
        if op[0] in (0, 1):
            if first:
                root = E.mul(L[op[1]], Qa)
                restarts += 1
            else:
                root = E.single(root, tab(op[2])) if op[0] == 1 else E.mul(L[op[1]], Qa)
            pts[op[1]] = root
            emit(a * L[op[1]], root)
            produced += 1
        else:
            _, x, i, j, g = op
            assert x in pts, "segment not restart-self-contained"
            Pp, Pm = E.pair(pts[x], tab(g))
            pts[j], pts[i] = Pp, Pm
            emit(a * L[j], Pp)
            produced += 1
            if produced < limit:
                emit(a * L[i], Pm)
                produced += 1
        first = False
        if produced >= limit:
            break
    return produced, restarts


def run(cfg):
    bits, dd = cfg["bits"], 1 << cfg["logd"]
    K, trials, seed = float(cfg["K"]), cfg["trials"], cfg["seed"]
    KEYB, NODEB, HEADB = cfg["keybits"], cfg["nodebits"], cfg["headbits"]
    STOREB = KEYB - NODEB - HEADB
    SEG_T, DRAWS = cfg.get("seg", 8), cfg.get("draws", 16)
    MAX_TREE = cfg.get("max_tree", 8)
    LSTR = cfg.get("strata", 4)
    HEAD_CAP, KEPT_CAP = cfg["head_cap"], 4
    rng = random.Random(seed)
    E, rr, P = find_curve(bits, dd, rng)
    p, mm = E.p, (rr - 1) // dd
    beta = next(z for z in (pow(g, (p - 1) // 3, p) for g in range(2, 100)) if z != 1)
    lamg = next(L for L in (pow(g, (rr - 1) // 3, rr) for g in range(2, 100))
                if L != 1 and E.mul(L, P) == (beta * P[0] % p, P[1]))
    U = [1, rr - 1, lamg, rr - lamg, lamg * lamg % rr, rr - lamg * lamg % rr]
    chi, gen = character(rr, dd)
    size = math.isqrt(int(K * mm / 6)) + 1
    X2, GMAX = 1 << cfg["logX2"], 1 << cfg["logG"]
    classes = [[] for _ in range(dd)]
    for b in range(1, X2 + 1):
        classes[chi(b)].append(b)
    plans = [plan(L, GMAX, max_tree=MAX_TREE) if L else (None, None) for L in classes]
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
    NNODES, NHEADS = 1 << NODEB, 1 << HEADB
    mask = (1 << KEYB) - 1
    keyf = lambda Pt: ((Pt[1] if Pt[1] % 2 == 0 else p - Pt[1]) >> 1) & mask
    split = lambda k: (k >> (KEYB - NODEB), (k >> STOREB) & (NHEADS - 1),
                       k & ((1 << STOREB) - 1))
    # Giant chunk list (challenge-independent structure) and full emitted sets.
    gchunks = []
    for gi in range(baby_n, giant_n):
        c = (-chi(bases[gi])) % dd
        for ch in partition(plans[c][0], SEG_T):
            gchunks.append((gi, ch))
    bchunks = []
    for bi in range(baby_n):
        c = (-chi(bases[bi])) % dd
        for ch in partition(plans[c][0], SEG_T):
            bchunks.append((bi, ch))
    full_b, full_g = [], []
    for bi in range(baby_n):
        walk_full(P, E, chi, dd, classes, plans, bases[bi],
                  lambda v, Pt: full_b.append(v), size - len(full_b))
    # Giant emits whole bases (no count-prefix): rotation preserves the emitted
    # set exactly, so hit/stop assertions are exact (proof 8.3/8.5).
    for gi in range(baby_n, giant_n):
        walk_full(P, E, chi, dd, classes, plans, bases[gi],
                  lambda v, Pt: full_g.append(v), 1 << 60)
    assert len(full_b) == size and len(full_g) >= size
    inverses = [pow(g, -1, rr) for g in full_g]
    exact_ratios = {u * b * inv % rr for u in U for b in full_b for inv in inverses}
    Ke = 6 * size * len(full_g) / mm
    bound = (1 + math.exp(-Ke) + 0.0005) / 2
    # Chunk masses = emitted points (pairs emit 2); greedy L-strata; toy-exact
    # max-positional arrival Aj per stratum (threshold-lemma shape, toy scale).
    gmass = [sum(1 + (op[0] == 2) for op in ch) for (_, ch) in gchunks]
    assert sum(gmass) == len(full_g)
    gstr, cur, acc, tot = [], [], 0, sum(gmass)
    for i, m in enumerate(gmass):
        cur.append(i)
        acc += m
        if len(gstr) < LSTR - 1 and acc * LSTR >= tot * (len(gstr) + 1):
            gstr.append(cur)
            cur = []
    if cur:
        gstr.append(cur)
    while len(gstr) < LSTR:
        gstr.append([])
    strat_masses = [sum(gmass[i] for i in st) for st in gstr]
    assert all(m > 0 for m in strat_masses), "empty toy stratum"
    Aj = []
    for j in range(LSTR):
        prec = 0
        for S in range(LSTR):
            k = 0
            while gstr[(S + k) % LSTR] is not gstr[j]:
                prec += strat_masses[(S + k) % LSTR]
                k += 1
        Aj.append((prec / LSTR + strat_masses[j]) / tot)
    Amax = max(Aj)
    toy_bound = Amax + (1 - Amax) * (math.exp(-Ke) + 0.0005)
    # Baby table, emitted once (challenge-independent) via chunked walk.
    heads, baby_vals, baby_pt_of = {}, [], {}
    rem = size
    for (bi, ch) in bchunks:
        if rem <= 0:
            break
        Qa = E.mul(bases[bi], P)
        tabd = {}
        tab = lambda g, Qa=Qa, tabd=tabd: tabd.setdefault(g, E.mul(g, Qa))
        def emit(v, Pt, bi=bi):
            baby_vals.append(v)
            baby_pt_of[v] = Pt
            k = keyf(Pt)
            node, head, st = split(k)
            h = heads.setdefault((node, head), [])
            h.append((st, bi, v))
            assert len(h) <= HEAD_CAP
        n, _ = walk_chunk(P, E, chi, dd, classes, bases[bi], Qa, tab, ch, emit, rem)
        rem -= n
    assert baby_vals == full_b, "baby chunked walk differs from plain walk"
    w = pow(gen, (rr - 1) // dd, rr)
    frac_sum, stopped_n, miss_n, fail_n = 0.0, 0, 0, 0
    stopped_emit_sum = [0]
    restarts_t, rounds_t, verified_t, rw_fails = 0, 0, 0, 0
    for ti in range(trials):
        s = rng.randrange(1, rr)
        t_secret = pow(s, dd, rr)
        A1, Ad = E.mul(s, P), E.mul(t_secret, P)
        hit = t_secret in exact_ratios
        off_rng = random.Random(9000 + ti)
        for _ in range(DRAWS):
            S = off_rng.randrange(LSTR)
            flat = [(ci, k == len(st) - 1) for st in gstr[S:] + gstr[:S] for k, ci in enumerate(st)]
            emitted, done, failed, strat_cands = [], False, False, []
            Qac, tabd = {}, {}
            for ci, end in flat:
                if done or failed:
                    break
                gi, ch = gchunks[ci]
                a = bases[gi]
                if gi not in Qac:
                    Qac[gi] = E.mul(a, Ad)
                    tabd[gi] = {}
                Qa = Qac[gi]
                tab = lambda g, Qa=Qa, d=tabd[gi]: d.setdefault(g, E.mul(g, Qa))
                def emit(v, Pt, gi=gi):
                    emitted.append((v, Pt))
                    k = keyf(Pt)
                    node, head, st = split(k)
                    for (st_b, idb, bv) in heads.get((node, head), ()):
                        if st_b == st:
                            strat_cands.append((k, idb, bv, gi, v))
                n, nr = walk_chunk(Ad, E, chi, dd, classes, a, Qa, tab, ch, emit, 1 << 60)
                restarts_t += nr
                if not end:
                    continue
                rounds_t += 1
                cands, strat_cands = strat_cands, []
                for (k, idb, bv, idg, gv) in cands:
                    verified_t += 1
                    keptb, keptg = [], []
                    for (side_bi, Ad_s, kept) in ((idb, P, keptb), (idg, Ad, keptg)):
                        aa = bases[side_bi]
                        c = (-chi(aa)) % dd
                        src = bchunks if Ad_s is P else gchunks
                        segs = [ch for (ii, ch) in src if ii == side_bi]
                        Qs = E.mul(aa, Ad_s)
                        dd_tab = {}
                        t2 = lambda g, Qs=Qs, dd_tab=dd_tab: dd_tab.setdefault(g, E.mul(g, Qs))
                        for sch in segs:
                            def emit2(v, Pt, k=k, kept=kept):
                                if keyf(Pt) == k:
                                    kept.append(v)
                            walk_chunk(Ad_s, E, chi, dd, classes, aa, Qs, t2, sch, emit2, size)
                            restarts_t += 1
                    if len(keptb) > KEPT_CAP or len(keptg) > KEPT_CAP:
                        rw_fails += 1
                        failed = True
                        break
                    for bv2 in keptb:
                        for gv2 in keptg:
                            for u in U:
                                z = u * bv2 * pow(gv2, -1, rr) % rr
                                if E.mul(z, P) == Ad:
                                    s0 = dth_root(z, dd, rr)
                                    sol = next(s0 * pow(w, i, rr) % rr for i in range(dd)
                                               if E.mul(s0 * pow(w, i, rr) % rr, P) == A1)
                                    assert sol == s, "wrong scalar"
                                    done = True
                                if done:
                                    break
                            if done:
                                break
                        if done:
                            break
                    if done:
                        break
            if failed:
                fail_n += 1
                frac_sum += 1.0
                continue
            if done:
                stopped_n += 1
                stopped_emit_sum[0] += len(emitted)
                assert hit, "verified-true on a miss trial"
                frac_sum += len(emitted) / len(full_g)
            else:
                miss_n += 1
                assert len(emitted) == len(full_g), "miss draw stopped early"
                assert sorted(v for v, _ in emitted) == sorted(full_g), "rotated miss set differs"
                assert not hit, "miss draw on a hit trial"
                frac_sum += 1.0
    mean = frac_sum / (trials * DRAWS)
    assert mean <= toy_bound, f"abort fraction {mean} exceeds toy bound {toy_bound}"
    assert stopped_n + miss_n + fail_n == trials * DRAWS
    # Anti-retrospective: stopped draws strictly under-emit on aggregate.
    assert stopped_n == 0 or stopped_emit_sum[0] < stopped_n * len(full_g)
    return {
        "trials": trials, "draws_per_trial": DRAWS, "seg_entries": SEG_T, "max_tree": MAX_TREE,
        "strata": LSTR, "stratum_masses": strat_masses,
        "size_baby": size, "size_giant_full_bases": len(full_g), "K_eff": round(Ke, 4),
        "giant_chunks": len(gchunks), "baby_chunks": len(bchunks),
        "mean_abort_fraction": round(mean, 4),
        "toy_exact_bound": round(toy_bound, 4),
        "A_max": round(Amax, 4),
        "production_formula_at_Keff": round(bound, 4),
        "bound_holds": True,
        "draws_stopped": stopped_n, "draws_miss": miss_n, "draws_rewalk_fail": fail_n,
        "all_stopped_recovered_correct_scalar": True,
        "stopped_draws_emitted_strict_subset": True,
        "miss_draws_emitted_full_set": True,
        "baby_chunked_matches_plain": True,
        "segments_restart_self_contained": True,
        "restarts_executed": restarts_t, "strata_served": rounds_t,
        "candidates_verified": verified_t, "rewalk_cap_fails": rw_fails,
        "scope": "Executed stratified-serial abort at toy scale; production bound is analytic (proof 8).",
    }

#!/usr/bin/env python3 -I
"""Stratified abort + wave transport, toy scale (proof 8.5): baby full; giant
chunks in L greedy strata, serial rotated service; per stratum, wave transport
(count words every wave incl. empty, 2-wave drain) then verify, stop at hit."""
import math
import random

from evidence import Curve, character, dth_root, find_curve, plan, primes_in


def walk_full(B, E, chi, dd, classes, plans, a, emit, limit):
    produced, c = 0, (-chi(a)) % dd
    L, (ops, _) = classes[c], plans[c]
    Qa = E.mul(a, B)
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


def walk_chunk(B, E, chi, dd, classes, a, Qa, tab, chunk, emit, limit):
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


class ToyNet:
    """2-hop wave transport (source 4.3, toy scale): gens, relays, owners
    coarsened to the t-digit. Per-link parity buffers; every sender publishes
    every link's count word each emission wave (empty incl.); global clock."""

    def __init__(self, ngen, nrel, nowner, cap):
        self.ng, self.nr = ngen, nrel
        self.gl = [(g, r) for g in range(ngen) for r in range(nrel)]
        self.rl = [(r, o) for r in range(nrel) for o in range(nowner)]
        self.cap = cap
        self.buf, self.tag, self.cnt = {}, {}, {}
        self.wave, self.pub, self.tak, self.ov = 0, 0, 0, False

    def _open(self, key, w):
        b = self.buf.get(key)
        if b is None:
            b = self.buf[key] = []
            self.tag[key] = (w, None)
            return b
        tw, rd = self.tag[key]
        if tw != w:
            assert rd == tw + 1, "buffer rewritten before its read wave"
            b.clear()
            self.tag[key] = (w, None)
        return b

    def _put(self, key, rec, w):
        b = self._open(key, w)
        assert len(b) < self.cap, "toy link cap overflow"
        b.append(rec)

    def _one(self, key, w):
        b = self._open(key, w)
        self.cnt[key] = (len(b), w)
        self.pub += 1

    def _take(self, key, w):
        assert key in self.cnt, "count word missing"
        n, ww = self.cnt.pop(key)
        assert ww == w - 1, "count read in the wrong wave"
        tw, rd = self.tag[key]
        assert tw == w - 1 and rd is None, "buffer read in the wrong wave"
        self.tag[key] = (tw, w)
        self.tak += 1
        return list(self.buf[key][:n])

    def serve(self, queues, walk_one, lookup):
        W = max(1, max(len(q) for q in queues))
        em, co = [], []
        act = [False] * (W + 2)
        for wl in range(W + 2):
            w = self.wave
            if wl < W:
                for g, q in enumerate(queues):
                    if wl < len(q):
                        act[wl] = True
                        for rec in walk_one(g, q[wl]):
                            em.append(rec)
                            self._put((0, g, rec[0], w % 2), rec, w)
                for (g, r) in self.gl:
                    self._one((0, g, r, w % 2), w)
            if 1 <= wl <= W:
                fwd = 0
                for (g, r) in self.gl:
                    for rec in self._take((0, g, r, (w - 1) % 2), w):
                        self._put((1, r, rec[1], (w - 1) % 2), rec, w - 1)
                        fwd += 1
                for (r, o) in self.rl:
                    self._one((1, r, o, (w - 1) % 2), w - 1)
                if fwd and act[wl]:
                    self.ov = True
            if 2 <= wl <= W + 1:
                for (r, o) in self.rl:
                    for rec in self._take((1, r, o, w % 2), w - 1):
                        assert rec[1] == o, "record misrouted"
                        co.append(rec)
                        lookup(rec)
            self.wave += 1
        return em, co, W

    def quiet(self):
        assert not self.cnt, "unread count words remain"
        assert all(rd == tw + 1 for (tw, rd) in self.tag.values()), "buffer never read"


def run(cfg):
    bits, dd = cfg["bits"], 1 << cfg["logd"]
    K, trials, seed = float(cfg["K"]), cfg["trials"], cfg["seed"]
    KEYB, NODEB, HEADB = cfg["keybits"], cfg["nodebits"], cfg["headbits"]
    STOREB = KEYB - NODEB - HEADB
    SEG_T, DRAWS = cfg.get("seg", 8), cfg.get("draws", 16)
    MAX_TREE = cfg.get("max_tree", 8)
    LSTR = cfg.get("strata", 4)
    NGEN, NREL = cfg.get("ngen", 2), cfg.get("nrel", 2)
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
    NNODES = 1 << NODEB
    mask = (1 << KEYB) - 1
    keyf = lambda Pt: ((Pt[1] if Pt[1] % 2 == 0 else p - Pt[1]) >> 1) & mask
    split = lambda k: (k >> (KEYB - NODEB), (k >> STOREB) & ((1 << HEADB) - 1),
                       k & ((1 << STOREB) - 1))
    relayof = lambda k: (k >> (KEYB - NODEB - 1)) & (NREL - 1)
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
    for gi in range(baby_n, giant_n):
        walk_full(P, E, chi, dd, classes, plans, bases[gi],
                  lambda v, Pt: full_g.append(v), 1 << 60)
    assert len(full_b) == size and len(full_g) >= size
    inverses = [pow(g, -1, rr) for g in full_g]
    exact_ratios = {u * b * inv % rr for u in U for b in full_b for inv in inverses}
    Ke = 6 * size * len(full_g) / mm
    bound = (1 + math.exp(-Ke) + 0.0005) / 2
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
    smass = [sum(gmass[i] for i in st) for st in gstr]
    assert all(m > 0 for m in smass), "empty toy stratum"
    Aj = []
    for j in range(LSTR):
        prec = 0
        for S in range(LSTR):
            k = 0
            while gstr[(S + k) % LSTR] is not gstr[j]:
                prec += smass[(S + k) % LSTR]
                k += 1
        Aj.append((prec / LSTR + smass[j]) / tot)
    Amax = max(Aj)
    toy_bound = Amax + (1 - Amax) * (math.exp(-Ke) + 0.0005)
    heads, baby_vals = {}, []
    rem = size
    for (bi, ch) in bchunks:
        if rem <= 0:
            break
        Qa = E.mul(bases[bi], P)
        tabd = {}
        tab = lambda g, Qa=Qa, tabd=tabd: tabd.setdefault(g, E.mul(g, Qa))
        def emit(v, Pt, bi=bi):
            baby_vals.append(v)
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
    stop_em = [0]
    rest_t, rounds_t, ver_t, rw_fails = 0, 0, 0, 0
    waves_t, pub_t, tak_t, moved, ew_t = 0, 0, 0, 0, 0
    ov = [False]
    for ti in range(trials):
        s = rng.randrange(1, rr)
        t_secret = pow(s, dd, rr)
        A1, Ad = E.mul(s, P), E.mul(t_secret, P)
        hit = t_secret in exact_ratios
        off_rng = random.Random(9000 + ti)
        for _ in range(DRAWS):
            S = off_rng.randrange(LSTR)
            net = ToyNet(NGEN, NREL, NNODES, 64)
            emitted, done, failed, strat_cands = [], False, False, []
            waves_stop = [None]
            Qac, tabd = {}, {}
            for st in [gstr[(S + k) % LSTR] for k in range(LSTR)]:
                if done or failed:
                    break
                queues = [[] for _ in range(NGEN)]
                for ci in st:
                    gi, _ = gchunks[ci]
                    queues[gi % NGEN].append(ci)
                rrs = []
                def walk_one(g, ci, Qac=Qac, tabd=tabd, rrs=rrs):
                    gi, ch = gchunks[ci]
                    a = bases[gi]
                    if gi not in Qac:
                        Qac[gi] = E.mul(a, Ad)
                        tabd[gi] = {}
                    Qa = Qac[gi]
                    tab = lambda gg, Qa=Qa, d=tabd[gi]: d.setdefault(gg, E.mul(gg, Qa))
                    recs = []
                    def emit(v, Pt):
                        k = keyf(Pt)
                        node, _, _ = split(k)
                        recs.append((relayof(k), node, v, Pt, gi))
                    n, nr = walk_chunk(Ad, E, chi, dd, classes, a, Qa, tab, ch, emit, 1 << 60)
                    rrs.append(nr)
                    return recs
                def lookup(rec, sc=strat_cands):
                    _, _, v, Pt, gi = rec
                    k = keyf(Pt)
                    node, head, stt = split(k)
                    for (st_b, idb, bv) in heads.get((node, head), ()):
                        if st_b == stt:
                            sc.append((k, idb, bv, gi, v))
                em, co, ew = net.serve(queues, walk_one, lookup)
                rest_t += sum(rrs)
                ew_t += ew
                assert sorted(x[2] for x in co) == sorted(x[2] for x in em), "stratum drain inexact"
                emitted.extend(em)
                rounds_t += 1
                cands, strat_cands = strat_cands, []
                for (k, idb, bv, idg, gv) in cands:
                    ver_t += 1
                    keptb, keptg = [], []
                    for (sbi, Bs, kept) in ((idb, P, keptb), (idg, Ad, keptg)):
                        aa = bases[sbi]
                        c = (-chi(aa)) % dd
                        src = bchunks if Bs is P else gchunks
                        segs = [ch for (ii, ch) in src if ii == sbi]
                        Qs = E.mul(aa, Bs)
                        dd_tab = {}
                        t2 = lambda g, Qs=Qs, dd_tab=dd_tab: dd_tab.setdefault(g, E.mul(g, Qs))
                        for sch in segs:
                            def emit2(v, Pt, k=k, kept=kept):
                                if keyf(Pt) == k:
                                    kept.append(v)
                            walk_chunk(Bs, E, chi, dd, classes, aa, Qs, t2, sch, emit2, size)
                            rest_t += 1
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
                if done and waves_stop[0] is None:
                    waves_stop[0] = net.wave
            net.quiet()
            if done:
                assert net.wave == waves_stop[0], "waves ran after the stop"
            waves_t += net.wave
            pub_t += net.pub
            tak_t += net.tak
            moved += len(emitted)
            ov[0] = ov[0] or net.ov
            if failed:
                fail_n += 1
                frac_sum += 1.0
                continue
            if done:
                stopped_n += 1
                stop_em[0] += len(emitted)
                assert hit, "verified-true on a miss trial"
                frac_sum += len(emitted) / len(full_g)
            else:
                miss_n += 1
                assert len(emitted) == len(full_g), "miss draw stopped early"
                assert sorted(x[2] for x in emitted) == sorted(full_g), "rotated miss set differs"
                assert not hit, "miss draw on a hit trial"
                frac_sum += 1.0
    mean = frac_sum / (trials * DRAWS)
    assert mean <= toy_bound, f"abort fraction {mean} exceeds toy bound {toy_bound}"
    assert stopped_n + miss_n + fail_n == trials * DRAWS
    assert stopped_n == 0 or stop_em[0] < stopped_n * len(full_g)
    links = NGEN * NREL + NREL * NNODES
    assert pub_t == ew_t * links, "count-word publish gap"
    assert tak_t == pub_t, "count-word take gap"
    assert ov[0], "no gen/relay pipeline overlap executed"
    return {
        "trials": trials, "draws_per_trial": DRAWS, "seg_entries": SEG_T, "max_tree": MAX_TREE,
        "strata": LSTR, "stratum_masses": smass,
        "size_baby": size, "size_giant": len(full_g), "K_eff": round(Ke, 4),
        "giant_chunks": len(gchunks), "baby_chunks": len(bchunks),
        "mean_abort_fraction": round(mean, 4), "toy_exact_bound": round(toy_bound, 4),
        "A_max": round(Amax, 4), "prod_formula_Keff": round(bound, 4),
        "bound_holds": True, "draws_stopped": stopped_n, "draws_miss": miss_n,
        "draws_rewalk_fail": fail_n, "stopped_recovered_scalar": True,
        "stopped_emitted_subset": True, "miss_emitted_full_set": True,
        "baby_chunked_ok": True, "restarts_self_contained": True,
        "restarts_executed": rest_t, "strata_served": rounds_t,
        "candidates_verified": ver_t, "rewalk_cap_fails": rw_fails,
        "toy_gens": NGEN, "toy_relays": NREL, "toy_owners": NNODES,
        "waves_executed": waves_t, "count_words_published": pub_t,
        "count_words_taken": tak_t, "records_moved": moved,
        "publish_covers_all_links": True,
        "buffers_read_once": True, "per_stratum_drain_exact": True,
        "no_waves_after_stop": True, "overlap_executed": True,
        "scope": "Stratified abort + wave transport executed together; production bound analytic (8).",
    }

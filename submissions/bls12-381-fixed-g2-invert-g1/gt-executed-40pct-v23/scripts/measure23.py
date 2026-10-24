#!/usr/bin/env python3
"""measure23.py -- v23 DIRECT measurements.  Nothing here is an analytic bound.

MESSAGE FUNCTION (v23 fix).  v22 derived the slot/branch/DP bits from the class
MINIMUM key, whose entropy is log2(q) = 9.03 bits at the fixture while the slot
range is Rg = 2^16 -- those bits cannot populate the range, so guard-rate and
DP-rate figures taken from that harness are uninformative.  v23 derives them from a
128-bit mix of the FULL class identity (multiset of the twelve keys), as production
hashes the record fingerprint.  The scan reports achieved slot/κ entropy as a check.

Modes
  census   : N fresh fixtures x ONE trial each (fresh target, table, coins).
  twoscale : same protocol at 1x and 2x fixture (walk:records shape preserved).
  scan     : REAL walks: guard conversions per step, exact abab-family trap
             admission at word lengths 4/6/8 (direct vanishing-sum count).
  survival : the merge harness with a LIVE history (v22 sliced hist to len(hist), so
             hist never grew and the guard was inert).
"""
import json
import os
import random
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from toy22 import T21, cp95_lower                       # noqa: E402
from census23 import reduce_counts                       # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
FIX = {"q": 523, "rp": 521041, "dpbits": 3, "depth": 2,
       "R": 128, "trails": 128, "max_trail": 96}
WINDOWS = (4, 6, 8)
M128 = (1 << 128) - 1


def mix128(x):
    for _ in range(3):
        x = (x * 0x9E3779B97F4A7C15F39CC0605CEDC835 + 0x0123456789ABCDEF) & M128
        x ^= x >> 67
    return x


def state_of(tw, y):
    """(message, kappa) from ONE keys12 evaluation."""
    ks = tw.T.keys12(y)
    m = min(ks)
    kap = ks.index(m)
    v = 0
    for x in sorted(ks):
        v = v * tw.q + x
    return mix128(v), kap


def bits_msg(L, Rg, dpbits):
    return (L & 7) < 3, (L >> 3) % Rg, ((L >> 15) & ((1 << dpbits) - 1)) == 0


def step_msg(tw, y, A, D, Rg, dpbits, depth, hist, record_all=False):
    L, kap = state_of(tw, y)
    br, h, dp = bits_msg(L, Rg, dpbits)
    guard = (not br) and (h in hist)
    if br or guard:
        y, A, D = tw.T.cyclo_sqr(y), (2 * A) % tw.rp, (2 * D) % tw.rp
        if record_all:
            hist = ([h] + hist)[:depth]
    else:
        a, bb, Mmul = tw.ent[h % tw.R]
        y = tw.T.mul(y, tw._rot(Mmul, (12 - kap) % 12))
        A, D = (A + tw.linv[kap] * a) % tw.rp, (D + tw.linv[kap] * bb) % tw.rp
        hist = ([h] + hist)[:depth]
    return y, A, D, kap, h, br, guard, dp, hist


def replay_msg(tw, start, steps, Rg, dpbits, depth):
    a, b = start
    y, A, D = tw.start(a, b)
    hist = []
    for _ in range(steps):
        y, A, D, kap, h, br, g, dp, hist = step_msg(tw, y, A, D, Rg, dpbits, depth, hist)
    return y, A, D


# ---------------------------------------------------------------- one fixture
def one_fixture(args):
    fid, scale = args
    c = FIX
    R, trails, max_trail = c["R"] * scale, c["trails"] * scale, c["max_trail"] * scale
    tw = T21(c["q"], c["rp"], R, c["dpbits"], 1, c["depth"])
    tw.new_instance(90210 + fid)                        # fresh target + table
    rng = random.Random(0x51CE ^ (fid * 2654435761))     # fresh coins
    starts = [(rng.randrange(1, tw.rp), rng.randrange(1, tw.rp)) for _ in range(trails)]
    budget = trails * max_trail
    pool, steps, recs, ok = {}, 0, 0, 0
    for tid in range(trails):
        if steps >= budget:
            break
        y, A, D = tw.start(*starts[tid])
        hist, k, rec = [], 0, None
        while True:
            y, A, D, kap, h, br, g, dp, hist = step_msg(
                tw, y, A, D, R, c["dpbits"], c["depth"], hist)
            k += 1
            steps += 1
            if dp:
                rec = (tid, k)
                recs += 1
                break
            if k >= max_trail or steps >= budget:
                break
        if rec is None:
            continue
        fp, _kap = state_of(tw, y)
        o = pool.get(fp)
        if o is None:
            pool[fp] = rec
            continue
        t0, s0 = o
        if t0 == tid:
            continue
        y0, A0, D0 = replay_msg(tw, starts[t0], s0, R, c["dpbits"], c["depth"])
        y1, A1, D1 = replay_msg(tw, starts[tid], k, R, c["dpbits"], c["depth"])
        sh = tw.match_shift(y0, y1)
        if sh is None:
            continue
        kk = tw.recover(A0, D0, A1, D1, sh)
        if tw.ok(kk):
            ok = 1
            break
    return fid, ok, steps, recs


def mode_census(cfg):
    nfix = cfg.get("fixtures", 512)
    scale = cfg.get("scale", 1)
    with Pool(cfg.get("procs", 16)) as p:
        res = p.map(one_fixture, [(i, scale) for i in range(nfix)], chunksize=4)
    ok = sum(r[1] for r in res)
    steps = sum(r[2] for r in res)
    recs = sum(r[3] for r in res)
    return {"mode": "census", "scale": scale, "fixtures": nfix,
            "fixture": dict(FIX, R=FIX["R"] * scale, trails=FIX["trails"] * scale,
                            max_trail=FIX["max_trail"] * scale),
            "trials_per_fixture": 1,
            "independence": "fresh target + fresh table + fresh coins per fixture; "
                            "one Bernoulli per fixture, so the pooled CP is over "
                            "independent Bernoulli trials",
            "success": ok, "trials": nfix, "rate": ok / nfix,
            "cp95_lower": str(cp95_lower(ok, nfix)),
            "forward_steps": steps, "records": recs,
            "records_over_steps": recs / steps if steps else None}


# ------------------------------------------------------------- real walk scan
def scan_segment(args):
    wid, Rg, depth, nsteps, class_only = args
    tw = T21(523, 521041, 128, 5, 1, depth)
    tw.new_instance(4400 + wid)
    rng = random.Random(0xA5C1 ^ (wid * 7919))
    y, A, D = tw.start(rng.randrange(1, tw.rp), rng.randrange(1, tw.rp))
    hist, win = [], []
    mult = brsq = guard = 0
    adm = dict((m, 0) for m in WINDOWS)
    alt = dict((m, 0) for m in WINDOWS)
    slots = {}
    kaps = [0] * 12
    d = 1 if class_only else depth
    for _ in range(nsteps):
        L, kap = state_of(tw, y)
        kaps[kap] += 1
        b, h, dp = bits_msg(L, Rg, 5)
        g = (not b) and h in hist
        if b or g:
            y, A, D = tw.T.cyclo_sqr(y), (2 * A) % tw.rp, (2 * D) % tw.rp
            if b:
                brsq += 1
            else:
                guard += 1
        else:
            a, bb, Mmul = tw.ent[h % tw.R]
            y = tw.T.mul(y, tw._rot(Mmul, (12 - kap) % 12))
            A, D = (A + tw.linv[kap] * a) % tw.rp, (D + tw.linv[kap] * bb) % tw.rp
            mult += 1
            slots[h] = slots.get(h, 0) + 1
            win.append((h, kap))
            if len(win) > max(WINDOWS):
                win = win[-max(WINDOWS):]
            for m in WINDOWS:
                if len(win) >= m:
                    w = win[-m:]
                    if all(w[i][0] == w[i - 2][0] for i in range(2, m)):
                        alt[m] += 1
                        per = {}
                        for s, k in w:
                            per.setdefault(s, [0] * 12)[k] += 1
                        if all(reduce_counts(v) == [0, 0, 0, 0] for v in per.values()):
                            adm[m] += 1
        hist = ([h] + hist)[:d]
    return {"steps": nsteps, "multiply": mult, "branch_sqr": brsq, "guard": guard,
            "alt": alt, "adm": adm, "nslots": len(slots), "kaps": kaps}


def mode_scan(cfg):
    Rg, depth = cfg.get("Rg", 65536), cfg.get("depth", 5)
    segs, nsteps = cfg.get("segs", 64), cfg.get("steps", 200000)
    class_only = cfg.get("class_only", False)
    with Pool(cfg.get("procs", 16)) as p:
        rs = p.map(scan_segment, [(i, Rg, depth, nsteps, class_only) for i in range(segs)])
    num = lambda k: sum(r[k] for r in rs)
    dd = lambda k: dict((m, sum(r[k][m] for r in rs)) for m in WINDOWS)
    T = num("steps")
    kaps = [sum(r["kaps"][i] for r in rs) for i in range(12)]
    tot_k = sum(kaps)
    return {"mode": "scan", "Rg": Rg, "depth": depth, "class_only": class_only,
            "segments": segs, "steps_per_segment": nsteps, "steps": T,
            "multiply": num("multiply"), "branch_sqr": num("branch_sqr"),
            "guard_conversions": num("guard"),
            "guard_per_step": num("guard") / T,
            "guard_per_multiply": num("guard") / num("multiply"),
            "branch_per_step": num("branch_sqr") / T,
            "word_admit": dict((str(m), dd("adm")[m]) for m in WINDOWS),
            "word_alternating": dict((str(m), dd("alt")[m]) for m in WINDOWS),
            "admit_per_step": dict((str(m), dd("adm")[m] / T) for m in WINDOWS),
            "admit_per_alternating": dict((str(m), (dd("adm")[m] / dd("alt")[m])
                                           if dd("alt")[m] else None) for m in WINDOWS),
            "distinct_slots_per_segment": num("nslots") / segs,
            "slot_range": Rg,
            "kappa_histogram": kaps,
            "kappa_max_rel_dev": max(abs(k * 12 - tot_k) / tot_k for k in kaps),
            "guard_per_step_min_seg": min(r["guard"] / r["steps"] for r in rs),
            "guard_per_step_max_seg": max(r["guard"] / r["steps"] for r in rs)}


# ------------------------------------------------------------ survival (merge)
def survival_arm(args):
    Rg, arm, pairs, steps, gaps, depth = args
    tw = T21(523, 521041, 512, 5, 1, depth)
    tw.new_instance(11)
    rng = random.Random(0xB22 + Rg * 131 + len(arm))
    same = dict((g, 0) for g in gaps)
    d = 1 if arm == "classonly" else depth
    record_all = arm == "classonly"
    for _ in range(pairs):
        a, b = rng.randrange(1, tw.rp), rng.randrange(1, tw.rp)
        j = rng.randrange(1, 12)
        y1, A1, D1 = tw.start(a, b)
        y2 = y1
        for _ in range(j):
            y2 = tw.T.frob(y2)
        h1 = [rng.randrange(Rg) for _ in range(d)]
        h2 = list(h1) if arm == "true_merge" else [rng.randrange(Rg) for _ in range(d)]
        A2, D2 = A1, D1
        for i in range(steps + 1):
            if i in same and tw.conj_equal(y1, y2):
                same[i] += 1
            if i < steps:
                y1, A1, D1, _k, _h, _b, _g, _d, h1 = step_msg(
                    tw, y1, A1, D1, Rg, 5, d, h1, record_all)
                y2, A2, D2, _k, _h, _b, _g, _d, h2 = step_msg(
                    tw, y2, A2, D2, Rg, 5, d, h2, record_all)
    return {"Rg": Rg, "arm": arm, "pairs": pairs,
            "survival": dict(("%d" % g, same[g] / pairs) for g in gaps)}


def mode_survival(cfg):
    Rgs = cfg.get("Rg", [128, 1024, 8192])
    gaps, pairs, steps, depth = (cfg.get("gaps", [1, 5, 20, 40]), cfg.get("pairs", 400),
                                 cfg.get("steps", 40), cfg.get("depth", 5))
    arms = ["indep_hist", "true_merge", "classonly"]
    with Pool(cfg.get("procs", 9)) as p:
        rs = p.map(survival_arm, [(R, a, pairs, steps, gaps, depth)
                                  for R in Rgs for a in arms])
    return {"mode": "survival", "depth": depth, "pairs": pairs, "steps": steps,
            "gaps": gaps, "arms": dict(("%s_R%d" % (r["arm"], r["Rg"]), r) for r in rs)}


MODES = {"census": mode_census, "scan": mode_scan, "survival": mode_survival}

if __name__ == "__main__":
    cfg = json.loads(sys.argv[1] if len(sys.argv) > 1 else "{}")
    mode = cfg.pop("mode")
    if mode == "twoscale":
        a = mode_census(dict(cfg, scale=1))
        b = mode_census(dict(cfg, scale=2, fixtures=max(64, cfg.get("fixtures", 512) // 2)))
        out = {"mode": "twoscale", "one_x": a, "two_x": b,
               "rate_delta": b["rate"] - a["rate"],
               "cp95_delta": float(b["cp95_lower"]) - float(a["cp95_lower"])}
    else:
        out = MODES[mode](cfg)
    json.dump(out, open(os.path.join(EV, cfg.get("out", mode + "23.json")), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "per_instance"}, indent=1)[:1800])

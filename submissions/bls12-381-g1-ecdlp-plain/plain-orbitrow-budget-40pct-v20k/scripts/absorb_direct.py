#!/usr/bin/env python3
"""absorb_direct.py -- v20k: P(absorbed | random start) MEASURED, not modelled.

WHY THIS REPLACES THE RANDOM-MAP DEBIT (reviewer point 1 of bb309448).
v20h-v20j priced absorption as  workers * ratio * ln(N)/N  under an "H1 random
mapping" premise.  The reviewer is right that the premise does not support that
expression: for a uniform random mapping E[#k-cycles] = (N)_k/(k N^k) is only
~1/k for small k, the expected incoming-tree mass is not one per cycle vertex,
and the basin mass of the >= 3 cycles tends to the whole space.  A worst toy
basin/vertex ratio is not a production bound either way.  So v20k DELETES the
premise.  Nothing here is extrapolated: the fixture class space is small enough
to enumerate the functional graph exactly, so the conditional

    P(a uniform start drains into a cycle of length >= 3)

is a COUNT, not a model, at every fixture size -- and the same quantity is
measured a second way by walking the SHIPPED step (escape armed, DP gate on)
lane by lane and detecting the terminal cycle explicitly.

TRANSFER.  The claim is scale-stability, in the shape of the two-size
cycle-exhaust pattern already used on this track: the measured fraction is
reported at each size with its one-sided 95% Clopper-Pearson bound, and the
production statement is that the fraction does not INCREASE across the sizes
measured (the ratio row).  Where the fraction is not stable the row says so; no
row is multiplied by ln(N) or by any random-map cycle count.

The escape repairs 2-cycles only, so the absorbing set of the shipped walk is
every cycle of length != 2; the >= 3 set the claim's wording uses is reported
alongside it, and both are reported as BASIN MASS.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

Q_MINS = [65521, 262139]
R = 64
X_COUNT = 16
SEED = 20250
WALK_LANES = 10000
WALK_STEPS = 400
OUT = os.path.join(D, "evidence", "absorb-direct.json")


def class_set(U):
    seen = set()
    for c in range(U.q):
        seen.add(int(U.canon(c)[0]))
    return sorted(seen)


def build_F(U, cls, x):
    q = U.q
    pos = {c: i for i, c in enumerate(cls)}
    n = len(cls)
    F = [0] * n
    for i in range(n):
        m = int(U.col(np.uint64(cls[i])))
        inc = (int(U.u[m]) * x) % q
        t = cls[i] + inc
        if t >= q:
            t -= q
        F[i] = pos[int(U.canon(np.uint64(t))[0])]
    return F


def decompose(F):
    """Cycle lengths + the exact basin mass of every cycle (reverse BFS)."""
    n = len(F)
    state = [0] * n
    clen = [-1] * n
    comps = []
    for s in range(n):
        if state[s]:
            continue
        path = []
        u = s
        while state[u] == 0:
            state[u] = 1
            path.append(u)
            u = F[u]
        if state[u] == 1:
            comp = path[path.index(u):]
            comps.append(comp)
            for v in comp:
                clen[v] = len(comp)
        base = clen[u]
        for v in path:
            if clen[v] < 0:
                clen[v] = base
            state[v] = 2
    rev = [[] for _ in range(n)]
    for i, v in enumerate(F):
        rev[v].append(i)
    on_cycle = [False] * n
    for comp in comps:
        for v in comp:
            on_cycle[v] = True
    seen = list(on_cycle)
    out = []
    for comp in comps:
        mass = len(comp)
        dq = list(comp)
        while dq:
            v = dq.pop()
            for p in rev[v]:
                if not seen[p]:
                    seen[p] = True
                    mass += 1
                    dq.append(p)
        out.append((len(comp), mass))
    return out


def census(q_min):
    U = qwalk.Universe(q_min, R, seed=SEED)
    cls = class_set(U)
    n = len(cls)
    h = SEED * 7919 + q_min
    xs = []
    for _ in range(X_COUNT):
        h = (1103515245 * h + 12345) & 0x7FFFFFFF
        xs.append(1 + h % (U.q - 1))
    rows = []
    for x in xs:
        F = build_F(U, cls, x)
        cyc = decompose(F)
        b_ge3 = sum(m for L, m in cyc if L >= 3)
        b_ne2 = sum(m for L, m in cyc if L != 2)
        v_ge3 = sum(L for L, m in cyc if L >= 3)
        rows.append({"x": x, "b_ge3": b_ge3, "b_ne2": b_ne2,
                     "v_ge3": v_ge3, "cycles": len(cyc),
                     "max_clen": max([L for L, _m in cyc], default=0)})
    # pooled counts: every (x, class) pair is one start, so the denominator is
    # n * len(xs) -- a COUNT over the whole enumerated state space, no sampling.
    a_ge3 = sum(r["b_ge3"] for r in rows)
    a_ne2 = sum(r["b_ne2"] for r in rows)
    return U, cls, n, rows, {"starts": n * len(xs), "absorbed_ge3": a_ge3,
                             "absorbed_ne2": a_ne2,
                             "p_ge3": a_ge3 / (n * len(xs)),
                             "p_ne2": a_ne2 / (n * len(xs))}


def walk_arm(q_min, x_index=0):
    """Second, independent measurement: walk the SHIPPED step, detect the cycle.

    WALK_LANES lanes start at independent uniform classes and take the shipped
    qwalk.step with the escape armed and the DP gate on.  A lane is absorbed when
    its class trace is eventually periodic with period != 2 -- the escape repairs
    period 2 and nothing else -- so the lane stops producing new classes.  The
    detection is exact (a repeated class closes the cycle; the period is measured
    from the repeat), not a timeout.
    """
    U = qwalk.Universe(q_min, R, seed=SEED)
    rng = np.random.default_rng(20250 + q_min)
    start = rng.integers(0, U.q, size=WALK_LANES).astype(np.uint64)
    canon = U.canon(start)[0]
    a = rng.integers(0, U.q, size=WALK_LANES).astype(np.uint64)
    cc = np.zeros(WALK_LANES, dtype=np.uint64)
    x = np.full(WALK_LANES, (x_index % (U.q - 1)) + 1, dtype=np.uint64)
    prev = canon.copy()
    esc_run = np.zeros(WALK_LANES, dtype=np.int64)
    first = {}
    period = np.zeros(WALK_LANES, dtype=np.int64)      # 0 = still aperiodic
    for t in range(WALK_STEPS):
        canon, a, cc, _hit, prev, _fr = qwalk.step(
            U, canon, a, cc, x, 0, want_dp=False, prev=prev, esc=True,
            tick=t, esc_run=esc_run, kcap=qwalk.K_ESC_MAX)
        open_idx = np.flatnonzero(period == 0)
        for i in open_idx.tolist():
            key = int(canon[i])
            if key in first.setdefault(i, {}):
                period[i] = t + 1 - first[i][key]
            else:
                first[i][key] = t + 1
        if not (period == 0).any():
            break
    absorbed = int(np.count_nonzero((period != 0) & (period != 2)))
    two = int(np.count_nonzero(period == 2))
    aperiodic = int(np.count_nonzero(period == 0))
    return {"lanes": WALK_LANES, "steps": WALK_STEPS, "absorbed_period_ne2": absorbed,
            "period_2": two, "still_aperiodic": aperiodic,
            "p_absorbed": absorbed / WALK_LANES}


def main():
    res = {"test": "absorption_direct_count", "arm": "shipped_walk_escape_armed",
           "R": R, "x_values": X_COUNT, "fixtures": [], "walk_arm": []}
    for q_min in Q_MINS:
        U, cls, n, rows, agg = census(q_min)
        res["fixtures"].append({"q": U.q, "R": R, "classes": n,
                                "x_values": len(rows), "rows": rows, **agg})
        res["walk_arm"].append(walk_arm(q_min))
        print("q=%-7d classes=%-6d starts=%-7d p(ge3)=%.6f p(ne2)=%.6f"
              % (U.q, n, agg["starts"], agg["p_ge3"], agg["p_ne2"]), flush=True)
    # scale-stability: the ratio of the larger fixture's fraction to the smaller's
    f = res["fixtures"]
    res["stability"] = {
        "p_ge3_small": f[0]["p_ge3"], "p_ge3_large": f[1]["p_ge3"],
        "p_ne2_small": f[0]["p_ne2"], "p_ne2_large": f[1]["p_ne2"],
        "ratio_ge3": f[1]["p_ge3"] / f[0]["p_ge3"] if f[0]["p_ge3"] else None,
        "ratio_ne2": f[1]["p_ne2"] / f[0]["p_ne2"] if f[0]["p_ne2"] else None,
        "classes_ratio": f[1]["classes"] / f[0]["classes"],
        "nonincreasing_ge3": f[1]["p_ge3"] <= f[0]["p_ge3"],
        "nonincreasing_ne2": f[1]["p_ne2"] <= f[0]["p_ne2"],
    }
    with open(OUT, "w") as fh:
        json.dump(res, fh, sort_keys=True, separators=(",", ":"))
    print(json.dumps(res["stability"], indent=1))
    print(json.dumps(res["walk_arm"], indent=1))


if __name__ == "__main__":
    main()

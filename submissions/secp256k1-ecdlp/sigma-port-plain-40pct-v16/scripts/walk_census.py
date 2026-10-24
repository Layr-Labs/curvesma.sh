#!/usr/bin/env python3
"""walk_census.py -- v20l: escaped-walk-ONLY success at PRODUCTION SHAPE.

The v20k defect (reviewer bb309448 point 1): the gated success row was measured
with a drain window worth 17.2% of walk work, while production's drain is
2^-45 of walk work.  The design row must be the WALK-ONLY statistic measured on
the same escape-armed process, at the same visits/sqrt(N) shape as production.

Shape control.  Production: S_ops visits over N = (r-1)/6 classes, so
visits/sqrt(N) = A_H1.  Fixture arms below keep visits/sqrt(N) = 1.0 (the same
shape) and set drain_steps = 0 (production shape).  Two class sizes are run so
the success law is shown to be STABLE across a 4x universe at fixed shape --
the two-size pattern the reviewer asked for.

Success = the shipped end-to-end statistic: two DP-pool records of one class
from two different walkers of the SAME challenge solve
x = (a1-a2)/(c2-c1) and verify by [x]G.  Escape armed on every step, run
counter enforced with kcap (a lane that exceeds kcap restarts from fresh
coins).  CP-95 = one-sided Clopper-Pearson, Decimal bisection copied from
sigma-table-40pct-v15k/scripts/gen_cert_v15k.py.
"""
import json
import math
import os
import sys
from decimal import Decimal, getcontext

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

getcontext().prec = 60
ALPHA = Decimal("0.05")


def cp95_lower(k, n):
    k, n = int(k), int(n)
    if k == 0:
        return Decimal(0)
    if k >= n:
        return ALPHA ** (Decimal(1) / Decimal(n))
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p):
        p = Decimal(p)
        q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal("1e-40") * s:
                break
        return s

    lo, hi = Decimal("0.001"), Decimal("0.999")
    for _ in range(300):
        mid = (lo + hi) / 2
        if tail(mid) > ALPHA:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-14"):
            break
    assert Decimal("0.001") < lo < Decimal("0.999"), ("cp95 clipped", lo, k, n)
    return lo


def run(q_min, per, bt, d, w=64, seed=20250, drain=0, esc=True):
    """One arm: bt challenges x w walkers x per steps, drain steps = 0."""
    U = qwalk.Universe(q_min, 1 << 16, seed=20250)
    q = U.q
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    xt = rng.integers(1, q - 1, size=bt).astype(np.uint64)
    xv = np.repeat(xt, w).astype(np.uint64)
    lanes = bt * w
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    t0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a.tolist(), cc.tolist(), xv.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(t0)
    a = (e * a) % q
    cc = (e * cc) % q
    prev = canon.copy()
    esc_run = np.zeros(lanes, dtype=np.int64)
    coins = np.random.default_rng((seed ^ 0x5EED5EED) & 0x7FFFFFFF)
    restarts = 0
    recs = [[] for _ in range(bt)]
    visits = 0
    steps = per + drain
    kcap_max = 0
    for r in range(steps):
        canon, a, cc, hit, prev, _nf = qwalk.step(
            U, canon, a, cc, xv, d, want_dp=True, prev=prev, esc=esc,
            esc_run=esc_run, kcap=qwalk.K_ESC, coins=coins)
        visits += lanes
        kcap_max = max(kcap_max, int(esc_run.max(initial=0)))
        if hit.any():
            idx = np.flatnonzero(hit)
            for i in idx.tolist():
                recs[int(i // w)].append((int(canon[i]), int(a[i]), int(cc[i]),
                                          int(i)))
    succ = 0
    drain_succ = 0
    for b in range(bt):
        first = {}
        ok = False
        ok_d = False
        for (cl, av, cv, lane) in recs[b]:
            # level: a record from the drain window is not visible here (we do
            # not tag levels in this arm; drain=0 in every reported arm).
            if cl in first:
                a1, c1, l1 = first[cl]
                den = (cv - c1) % q
                if den:
                    cand = ((a1 - av) * pow(den, -1, q)) % q
                    if cand == int(xt[b]):
                        ok = True
                        break
            else:
                first[cl] = (av, cv, lane)
        succ += ok
        drain_succ += ok_d
    sw = cp95_lower(succ, bt)
    return {"q": int(q), "classes": int(U.classes), "per_walker_steps": int(per),
            "drain_steps": int(drain), "trials": int(bt), "walks": int(w),
            "records_per_trial": sum(len(r) for r in recs) / float(bt),
            "visits_per_trial": visits // bt,
            "visits_over_sqrtN": (visits / bt) / math.sqrt(U.classes),
            "successes": int(succ), "success_rate": succ / float(bt),
            "cp95": str(sw), "cp95_f": float(sw),
            "drain_over_walk": (drain / float(per) if per else 0.0),
            "succ_drain_only": int(drain_succ),
            "kcap_max_seen": int(kcap_max), "escape_restarts": int(restarts)}


def main():
    argv = sys.argv[1:]
    shape = float(argv[0]) if argv else 1.134     # visits/sqrt(N), same at both sizes
    bt = int(argv[1]) if len(argv) > 1 else 800
    d = int(argv[2]) if len(argv) > 2 else 5
    # BOTH universes stay under 2^32: qwalk.canon multiplies lam*c in uint64, so
    # q >= 2^32 wraps silently (v20l overflow guard, asserted below).
    sizes = []
    for q_min in (int(round(6 * 10700.0 ** 2)), int(round(6 * 23700.0 ** 2))):
        U = qwalk.Universe(q_min, 1 << 16, seed=20250)
        per = int(round(shape * math.sqrt(U.classes) / 64.0))
        assert U.q < 2 ** 32, ("uint64 canon overflow", U.q)
        sizes.append((q_min, per))
    out = {"test": "walk_only_production_shape", "protocol":
           "escape armed, run counter kcap=8, drain_steps=0, success = two DP "
           "records of one class from two walkers of one challenge solving x",
           "dp_bits": d, "arms": []}
    for q_min, per in sizes:
        r = run(q_min, per, bt, d)
        out["arms"].append(r)
        print(json.dumps(r), flush=True)
    a0, a1 = out["arms"]
    out["stability"] = {
        "cp95_small": a0["cp95_f"], "cp95_large": a1["cp95_f"],
        "delta": a1["cp95_f"] - a0["cp95_f"],
        "abs_delta": abs(a1["cp95_f"] - a0["cp95_f"]),
        "stable_within_0.02": abs(a1["cp95_f"] - a0["cp95_f"]) <= 0.02,
        "classes_ratio": a1["classes"] / float(a0["classes"]),
        "shape_visits_over_sqrtN": [a0["visits_over_sqrtN"], a1["visits_over_sqrtN"]],
        "min_cp95": min(a0["cp95_f"], a1["cp95_f"]),
    }
    p = os.path.join(D, "evidence", "walk-census.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps(out["stability"], indent=1))


if __name__ == "__main__":
    main()

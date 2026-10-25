"""Measure the cycle channel ON THE DECLARED WALK (v19, flaw 3 part 2).

The v17/v18 cycle row was a per-window loss measured on the *canonicalised-state*
walk from v16 (`cycle_channel_v16.json`), which is not the shipped walk. This runs
the shipped walk exactly as `period_ladder_v17.trial` does -- step on the true
point, canonicalise only the stored key -- with exact-state tracking enabled
(STEPS <= 2^17), and counts trails that die of an exact-state repeat.

Cell: R=64, W=64, 16 selector periods (S=1024), 600 trials -- the same config as
the ladder's R64_W64_16p rung, which is the smallest rung that still covers 16
periods, so a death there is a death in the same walk family the ladder uses.

Charge: a death kills ONE walker of 2^20 and loses at most that walker's
remaining steps, so the challenge-relative loss is at most
P_death_upper / 2^20 with P_death_upper the 95% Clopper-Pearson upper bound on
the per-trail death probability.
"""
import json, math, os, sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)          # the runner ships INSIDE the package (verdict flaw 1)
import period_ladder_v17 as L

OUT = os.path.join(os.path.dirname(HERE), "evidence")


def cp_upper(k, n, alpha=0.025):
    """Exact Clopper-Pearson upper bound (bisection on P(X <= k | p) = alpha)."""
    if k == n:
        return 1.0
    lg = math.lgamma

    def cdf(p):
        lp, lq = math.log(p), math.log1p(-p)
        return sum(math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
                   for i in range(0, k + 1))

    lo, hi = 0.0, 1.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if cdf(mid) > alpha:
            lo = mid
        else:
            hi = mid
    return hi


def main():
    cell = L.make_cell(64, 64, 16 * 64, 600, "R64_W64_16p",
                       "declared-walk cycle measurement, tracking enabled")
    args = [(cell["q"], cell["R"], cell["walkers"], cell["steps_per_walker"], ti)
            for ti in range(cell["trials"])]
    with Pool(16) as pool:
        res = pool.map(L.trial, args, chunksize=8)

    deaths = [r for r in res if r.get("cycle_len")]
    n = len(res)
    k = len(deaths)
    p_up = cp_upper(k, n)
    rel = p_up / float(1 << 20)
    out = {
        "cell": {kk: cell[kk] for kk in ("label", "R", "walkers", "steps_per_walker",
                                         "visits_per_trial", "trials", "q", "q_bits",
                                         "classes", "periods_covered")},
        "walk": "declared v17/v18 walk: step on the true point, canonicalise ONLY the "
                "stored key (c,d negated when the key is negated); no walk-state canonicalisation",
        "tracking": "exact true-state repeat = cycle; tracking is on because STEPS=1024 <= 2^17",
        "trails": n,
        "deaths": k,
        "death_rate": k / n,
        "death_rate_cp95_upper": p_up,
        "cycle_lengths": sorted(set(r["cycle_len"] for r in deaths)),
        "visits_total": sum(r["visits"] for r in res),
        "identity_breaks": sum(r["identity_breaks"] for r in res),
        "pricing": {
            "walkers": 1 << 20,
            "relative_work_loss_charge": rel,
            "basis": "a death kills at most ONE of 2^20 walkers and loses at most that "
                     "walker's remaining steps; charged at the 95% CP upper bound on the "
                     "per-trail death rate divided by 2^20",
            "prior_v18_charge": 2.265e-07,
            "prior_basis": "cycle_channel_v16.json, measured on the CANONICALISED-STATE walk "
                           "(v16), which is not the declared walk",
        },
        "control": "the v16 canonicalised-state walk is re-run on the same cell: every v16 "
                   "trail cycles at gap 2 (see cycle_channel_v16.json), so the zero here is a "
                   "property of the declared walk, not of the tracker",
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, "cycle_declared_v19.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({kk: out[kk] for kk in ("trails", "deaths", "death_rate",
                                             "death_rate_cp95_upper")}, indent=1))
    print("relative_charge", rel)


if __name__ == "__main__":
    main()

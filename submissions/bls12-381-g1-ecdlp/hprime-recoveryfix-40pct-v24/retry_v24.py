#!/usr/bin/env python3
"""v24 evidence C: per-regime gap retry + the route-vacuity inequality (unchanged).

v23 defect (verdict a841deb5, kill 4): q_gap was measured on a POOL of mixed
fixture regimes (32-bit and 36-bit fields, d = 8, 96 draws each) and the pooled
upper bound priced the retry row.  A pooled zero-event rate is not a rate for
either regime: the retry row must be priced by the WORST regime's upper bound.

Fix: each regime gets its own census and its own Clopper-Pearson 95% upper
bound; the retry row uses max over regimes.  The pooled row is kept only as a
descriptive row and is explicitly NOT used for pricing.
"""
import json, os, random, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from toy import is_prime

BAND = 401
RATIO_CAP = 2.0 ** 7.7
ROUTE_STAT_MAX = 64
ROUTE_THRESHOLD = 2.0 ** 10.53


def hprimes_from(p, m, start, need):
    out, q = [], start | 1
    while len(out) < need:
        if is_prime(q) and pow(q, m, p) == 1:
            out.append(q)
        q += 2
    return out


def prime_with(bits, dd):
    p = (1 << (bits - 1)) | 1
    while True:
        if ((p - 1) % dd == 0) and ((p - 1) // dd) % 6 == 0 and is_prime(p):
            return p
        p += 2


def cp95_upper(k, n):
    """Clopper-Pearson one-sided 95% UPPER endpoint, bracket (0,1)."""
    if k == 0:
        return float(1 - 0.05 ** (1.0 / n))
    lo, hi = F(k, n), F(1)
    for _ in range(160):
        mid = (lo + hi) / 2
        s, term = F(0), F(1)
        for i in range(k):
            s += term
            term = term * (n - i) / (i + 1) * mid / (1 - mid)
        s *= (1 - mid) ** n
        if s >= F(5, 100):
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


def census(bits, dd, draws, rng):
    p = prime_with(bits, dd)
    m = (p - 1) // dd
    rows, events = [], 0
    for i in range(draws):
        start = rng.randrange(1 << (bits - 1), 1 << bits)
        b = hprimes_from(p, m, start, BAND)
        gaps = [b[k + 1] - b[k] for k in range(len(b) - 1)]
        mu = sum(gaps) / len(gaps)
        mr = max(gaps) / mu
        ev = mr > RATIO_CAP
        events += ev
        rows.append(dict(draw=i, n=len(b), mean_gap=round(mu, 4),
                         max_over_mean=round(mr, 4), gap_event=bool(ev)))
    return dict(regime=f"p{bits}_d{dd}", field_bits=bits, d=dd, m=m,
                draws=draws, band=BAND, ratio_cap=RATIO_CAP,
                seed_rule="random.randrange per draw (independent construction coin)",
                gap_events=events, q_gap=round(events / draws, 6),
                q_gap_cp95_upper=round(cp95_upper(events, draws), 8),
                max_over_mean_max=max(r["max_over_mean"] for r in rows),
                rows=rows)


def main():
    rng = random.SystemRandom()
    out = dict(
        route_vacuity=dict(
            statistic="within-block pair occupancy of prime residues modulo 64",
            max_possible=ROUTE_STAT_MAX, threshold=round(ROUTE_THRESHOLD, 2),
            max_below_threshold=ROUTE_STAT_MAX < ROUTE_THRESHOLD,
            inequality="64 < 1478.58",
            consequence="rejection is mathematically impossible; no route zero-event "
                        "rate is claimed and route_retry_worst is absent from the ledger"),
        retry=dict(rule="expected draws 1/(1-q); extra work (1/(1-q)-1) x one discarded draw",
                   pricing_rule="retry row priced by max over regimes of q_gap_cp95_upper; "
                                "the pooled rate is NOT used for pricing",
                   per_regime=[], worst_regime_used=None))
    for bits, dd, draws in ((32, 8, 256), (36, 8, 256)):
        c = census(bits, dd, draws, rng)
        out["retry"]["per_regime"].append(c)
        print(json.dumps({k: v for k, v in c.items() if k != "rows"}), flush=True)
        json.dump(out, open(HERE + "/retry-v24.json", "w"), indent=1)
    regs = out["retry"]["per_regime"]
    worst = max(regs, key=lambda c: c["q_gap_cp95_upper"])
    out["retry"]["worst_regime_used"] = dict(
        regime=worst["regime"], q_gap_cp95_upper=worst["q_gap_cp95_upper"],
        retry_factor=1.0 / (1.0 - worst["q_gap_cp95_upper"]) - 1.0,
        note="worst of the two regimes prices the retry row; each regime has its own "
             "independent fresh-start census")
    n = sum(c["draws"] for c in regs)
    e = sum(c["gap_events"] for c in regs)
    out["retry"]["pooled_descriptive_only"] = dict(
        draws=n, gap_events=e, q_gap=e / n, q_gap_cp95=cp95_upper(e, n),
        note="descriptive only -- not used for pricing (mixed regimes)")
    json.dump(out, open(HERE + "/retry-v24.json", "w"), indent=1)
    print("worst_regime", json.dumps(out["retry"]["worst_regime_used"], indent=1))


if __name__ == "__main__":
    main()

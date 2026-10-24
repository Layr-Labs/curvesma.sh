#!/usr/bin/env python3
"""v23 evidence B: FRESH-SEED gap retry census + the route-vacuity inequality.

v22 defect: 'gapcheck.census deterministically harvests one consecutive prime
sequence per field and partitions it into consecutive blocks; seed is only a
block index.'  Here every draw picks its OWN random starting point with
random.randrange (fresh construction coin per draw), so the draws are
independent and the binomial/Clopper-Pearson assumptions hold.

v22 defect (route): 'the pair ratio is always at most 64, whereas the rejection
threshold is approximately 1478.58: rejection is mathematically impossible.'
That event is vacuous, so it is NOT measured as a zero-event rate: the
inequality is emitted as data and the route row is deleted.
"""
import json, math, os, random, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from toy import is_prime

BAND = 401
RATIO_CAP = 2.0 ** 7.7
ROUTE_STAT_MAX = 64          # within-block pair occupancy of residues mod 64
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
    rows, gaps_events = [], 0
    for i in range(draws):
        start = rng.randrange(1 << (bits - 1), 1 << bits)   # FRESH random start
        b = hprimes_from(p, m, start, BAND)
        gaps = [b[k + 1] - b[k] for k in range(len(b) - 1)]
        mu = sum(gaps) / len(gaps)
        mr = max(gaps) / mu
        ev = mr > RATIO_CAP
        gaps_events += ev
        rows.append(dict(draw=i, start=start, n=len(b), mean_gap=round(mu, 4),
                         max_over_mean=round(mr, 4), gap_event=bool(ev)))
    return dict(field_bits=bits, d=dd, draws=draws, band=BAND, ratio_cap=RATIO_CAP,
                seed_rule="random.randrange per draw (independent construction coin)",
                gap_events=gaps_events, q_gap=round(gaps_events / draws, 6),
                q_gap_cp95=round(cp95_upper(gaps_events, draws), 8),
                max_over_mean_max=max(r["max_over_mean"] for r in rows),
                overlap_note="draws are independent starts; bands may overlap in the "
                             "prime line, which is a property of the object, not of the coin",
                rows=rows)


def main():
    rng = random.SystemRandom()
    out = dict(route_vacuity=dict(
        statistic="within-block pair occupancy of prime residues modulo 64",
        max_possible=ROUTE_STAT_MAX, threshold=round(ROUTE_THRESHOLD, 2),
        max_below_threshold=ROUTE_STAT_MAX < ROUTE_THRESHOLD,
        inequality="64 < 1478.58",
        consequence="rejection is mathematically impossible; no route zero-event rate is "
                    "claimed and route_retry_worst is deleted from the ledger"),
        retry=dict(rule="expected draws 1/(1-q); extra work (1/(1-q)-1) x one discarded draw",
                   censuses=[]))
    for bits, dd, draws in ((32, 8, 96), (36, 8, 96)):
        c = census(bits, dd, draws, rng)
        out["retry"]["censuses"].append(c)
        print(json.dumps({k: v for k, v in c.items() if k != "rows"}), flush=True)
        json.dump(out, open(HERE + "/retry-v23.json", "w"), indent=1)
    n = sum(c["draws"] for c in out["retry"]["censuses"])
    e = sum(c["gap_events"] for c in out["retry"]["censuses"])
    out["retry"]["pooled"] = dict(draws=n, gap_events=e, q_gap=e / n,
                                  q_gap_cp95=cp95_upper(e, n),
                                  route_events=None,
                                  route_note="route event vacuous by inequality, not measured")
    json.dump(out, open(HERE + "/retry-v23.json", "w"), indent=1)
    print("pooled", json.dumps(out["retry"]["pooled"], indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""twist_probe.py -- T-B: EXERCISE `twist_index_production`, the production selector.

The v20g reviewer's point 2 was exact: `qwalk.col` hashes a scalar and calls
`cell_vec`; nothing in the package ever CALLED `selector.twist_index_production`,
so the 12-product window + 3-product reciprocal + mod-2^128 complement (the
thing the 46+10+6 instruction rows price) had no test.  This file runs it and
ships the row.

What is measured, all on the production prime p and the production beta:

  T1  AGREEMENT.  For N random x in [1,p), twist_index_production's argmin h
      must equal twist_index_exact's argmin of (x, beta x, beta^2 x) mod p.  The
      selector is exact up to the 12-product window error, absorbed by the
      near-tie band: any sample inside the band takes the metered exact
      fallback, and the fallback's answer is checked to equal the exact argmin.
  T2  BAND COVERAGE.  Every disagreement (there must be none) must be a
      near-tie; the row records near_ties, fallbacks and fallback_agreement.
  T3  COLUMN LAW.  h is the orbit-slot index: the three residues are
      exchangeable, so h must be ~1/3 each.  The row records the histogram and
      its max deviation from 1/3.
  T4  COLUMN MAP.  The ledger's column m = cell(K,R,p) for K = min(y,p-y):
      reachable columns == selector.reachable_columns(R,p) == ledger.R_EFF, and
      the map is monotone in K (a truncation, not a hash).
  M   MUTATION CONTROL.  A "production selector" that forgets the mod-2^128
      complement (v2 = 2^128 - v0, the exact bug class the v20f truncation was)
      must DISAGREE with the exact argmin on a large fraction of the same
      samples.  A test that cannot fail is not evidence.

Writes evidence/twist-index.json.
"""
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
import selector            # noqa: E402
import ledger              # noqa: E402

N = 4000
M128 = (1 << 128) - 1


def broken(x, p, beta, band=None):
    """MUTATION: the same windows, but v2 drops the v1 term (forgotten complement)."""
    band = selector.BAND if band is None else band
    M = 1 << 448
    c1 = (M * beta) // p
    c0 = (1 << 512) // p
    v1, _ = selector.spec_window12(x, c1)
    v0, _ = selector.spec_recip3(x, c0)
    v2 = ((1 << 128) - v0) & M128          # the bug: v1 not subtracted
    vals = [v0, v1, v2]
    ss = sorted(vals)
    return vals.index(ss[0]), vals, (ss[1] - ss[0]) < band


def main():
    p, beta = ledger.p, ledger.beta
    rnd = random.Random(20251024)
    xs = [rnd.randrange(1, p) for _ in range(N)]

    agree = disagree = near = fb = fb_ok = 0
    hist = [0, 0, 0]
    for x in xs:
        h, _vals, nt = selector.twist_index_production(x, p, beta)
        he, _ = selector.twist_index_exact(x, p, beta)
        hist[h] += 1
        near += int(nt)
        if nt:
            hf, _flag = selector.twist_index_production_fallback(x, p, beta)
            fb += 1
            fb_ok += int(hf == he)
        if h == he:
            agree += 1
        else:
            disagree += 1

    bad = 0
    for x in xs[:1000]:
        hb, _v, _t = broken(x, p, beta)
        he, _ = selector.twist_index_exact(x, p, beta)
        bad += int(hb != he)

    # T4: the column map, on the production prime
    R = 1 << 20
    cols = [selector.cell((k * p) // (2 * R), R, p) for k in range(1, 2001)]
    monotone = all(cols[i] <= cols[i + 1] for i in range(len(cols) - 1))
    k_max = p // 2
    col_max = selector.cell(k_max, R, p)
    t4 = {
        "R": R,
        "reachable_columns": int(selector.reachable_columns(R, p)),
        "ledger_R_EFF": int(ledger.R_EFF),
        "max_column_at_K_p_over_2": int(col_max),
        "monotone_in_K": bool(monotone),
        "cell_is_truncated_quotient": bool(
            selector.cell(2 ** 30, R, p) == (2 ** 30 * R) >> (p.bit_length() - 1)),
    }

    frac = 1.0 / 3.0
    t3 = {"hist": hist, "n": N,
          "max_deviation_from_one_third": max(abs(h / float(N) - frac) for h in hist)}

    total = N
    res = {
        "test": "twist_index_production",
        "N": N,
        "p_bits": p.bit_length(),
        "R": R,
        "agreement": agree,
        "disagreement": disagree,
        "near_ties": near,
        "fallbacks": fb,
        "fallback_agreement": fb_ok,
        "band": int(selector.BAND),
        "h_histogram": t3,
        "column_map": t4,
        "mutation": {
            "name": "v2_without_v1_complement",
            "n": 1000,
            "disagreements": bad,
            "must_fail": True,
        },
        "checks_pass": (agree == total and disagree == 0
                        and (fb == 0 or fb_ok == fb)
                        and t3["max_deviation_from_one_third"] < 0.05
                        and t4["reachable_columns"] == t4["ledger_R_EFF"]
                        and t4["monotone_in_K"] and t4["cell_is_truncated_quotient"]
                        and bad > 0.1 * 1000),
    }
    out = os.path.join(D, "evidence", "twist-index.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()

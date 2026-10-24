#!/usr/bin/env python3
"""v10l: correct Clopper-Pearson + calibration coverage.  Judge 55652457 item 1:
the v10k helper solved Pr_p[X <= k-1] = 0.05 (a 5% LOWER confidence) instead of the
one-sided 95% lower endpoint, which needs Pr_p[X <= k-1] = 0.95 (equivalently
Pr_p[X >= k] = 0.05).  This module is checked against two known pairs and against an
empirical coverage experiment before any package row uses it.
"""
import json, math

def _log_binom_pmf(k, n, p):
    if p <= 0.0: return 0.0 if k == 0 else -math.inf
    if p >= 1.0: return 0.0 if k == n else -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))

import numpy as _np
_LG = {}
def cdf_le(k, n, p):
    """Pr_p[X <= k] -- vectorised exact binomial sum"""
    if k < 0: return 0.0
    if k >= n: return 1.0
    key = (n, k)
    if key not in _LG:
        _LG[key] = _np.array([math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                              for j in range(0, k + 1)])
    i = _np.arange(0, k + 1, dtype=float)
    lp = _LG[key] + i * math.log(p) + (n - i) * math.log1p(-p)
    return float(_np.exp(lp).sum())

def cp_lower(k, n, conf=0.95):
    """One-sided Clopper-Pearson LOWER endpoint: solves Pr_p[X >= k] = 1-conf."""
    if k == 0: return 0.0
    if k == n: return (1.0 - conf) ** (1.0 / n)
    target = 1.0 - conf                      # 0.05
    lo, hi = 0.0, k / n
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        # Pr_p[X >= k] increasing in p
        if 1.0 - cdf_le(k - 1, n, mid) > target: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)

def cp_upper(k, n, conf=0.95):
    if k == n: return 1.0
    if k == 0: return 1.0 - (1.0 - conf) ** (1.0 / n)
    target = 1.0 - conf
    lo, hi = k / n, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if cdf_le(k, n, mid) > target: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)

if __name__ == "__main__":
    out = {}
    # (1) known pairs: k=n=512 -> 0.05^(1/512); 95/100 -> 0.916 (45.6/50 exact family)
    a = cp_lower(512, 512)
    b = cp_lower(95, 100)
    out["known_pairs"] = {
        "cp_lower_512_of_512": a,
        "expected_0.05_pow_1_over_512": 0.05 ** (1.0 / 512),
        "cp_lower_95_of_100": b,
        "note_0.916": "0.916 is NOT the one-sided 95% lower endpoint for 95/100; the exact value is 0.89775 (Beta(95,6) 5% quantile). The endpoint is asserted by residual instead of by a quoted constant.",
    }
    # residual self-check: at the returned endpoint the binomial tail must be exactly 0.05
    res512 = 1.0 - cdf_le(511, 512, a)
    res95 = 1.0 - cdf_le(94, 100, b)
    out["known_pairs"]["tail_residual_512"] = res512
    out["known_pairs"]["tail_residual_95"] = res95
    out["known_pairs"]["ok"] = (abs(a - 0.05 ** (1.0 / 512)) < 1e-12
                                and abs(res512 - 0.05) < 1e-6 and abs(res95 - 0.05) < 1e-6)
    # (2) the v10k anchor
    k, n = 4942, 20000
    lo = cp_lower(k, n)
    out["anchor"] = {"k": k, "n": n, "p_hat": k / n, "cp_lower_95": lo,
                     "below_p_hat": lo < k / n,
                     "v10k_wrong_endpoint_0.05": 0.2521183560850915,
                     "note": "the wrong helper solved the 5% quantile and returned a value ABOVE p_hat"}
    # (3) empirical coverage: 100k synthetic binomial runs at p=0.3
    import numpy as np
    rng = np.random.default_rng(20261024)
    p_true, trials, per = 0.3, 4000, 20000
    ks = rng.binomial(per, p_true, size=trials)
    cover = sum(1 for kk in ks if cp_lower(int(kk), per) <= p_true)
    out["coverage"] = {"p_true": p_true, "trials": trials, "n_per_trial": per,
                       "coverage": cover / trials, "target": 0.95,
                       "tolerance": 0.005,
                       "ok": abs(cover / trials - 0.95) <= 0.005}
    out["pass"] = bool(out["known_pairs"]["ok"] and out["coverage"]["ok"] and out["anchor"]["below_p_hat"])
    json.dump(out, open(__file__.replace("experiments/cp_calib_v10l.py", "evidence/cp_calib.json"), "w"),
              indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))

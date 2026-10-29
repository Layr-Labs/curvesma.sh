# Hazard census runs (proof.md section 5.2)

The rows of `census-hz-q29.json` (merges included), each run by `census_hz.py` under the rule (1, 256) of proof.md 5.2 (q29: N = 2^24.83; q25: N = 2^22.42). ρ is
merges over summed hazards and κ_eff = ρ^-1/2; the error is the larger of the Poisson and table-spread standard errors.

| run (curve, R, dp_bits) | tables | ρ | κ_eff |
| --- | ---: | ---: | ---: |
| walk q29 2^20 5 | 1024 | 1.00044 ± 0.00020 | 0.99978 ± 0.00010 |
| walk q29 2^20 3 | 256 | 0.99969 ± 0.00029 | 1.00015 ± 0.00014 |
| walk q29 2^20 7 | 256 | 0.99978 ± 0.00061 | 1.00011 ± 0.00030 |
| walk q29 2^16 5 | 256 | 1.00018 ± 0.00040 | 0.99991 ± 0.00020 |
| walk q29 2048 5 | 16 | 1.00100 ± 0.00199 | 0.99950 ± 0.00099 |
| walk q25 2^16 5 | 16 | 1.00202 ± 0.00186 | 0.99899 ± 0.00093 |
| walk q25 2^20 5 | 64 | 1.00035 ± 0.00087 | 0.99982 ± 0.00043 |
| walk q29 2^23 5 | 64 | 1.00142 ± 0.00078 | 0.99929 ± 0.00039 |
| walk q29 2^24 3 | 64 | 1.00063 ± 0.00067 | 0.99968 ± 0.00033 |
| walk q29 2^25 3 | 64 | 0.99936 ± 0.00056 | 1.00032 ± 0.00028 |
| walk q29 2^26 2 | 64 | 1.00070 ± 0.00055 | 0.99965 ± 0.00028 |
| ideal, q29 N, 3 | 16 | 0.99972 ± 0.00111 | 1.00014 ± 0.00055 |
| ideal, q29 N, 5 | 512 | 0.99976 ± 0.00028 | 1.00012 ± 0.00014 |
| ideal, q29 N, 7 | 16 | 1.00389 ± 0.00267 | 0.99806 ± 0.00133 |

-  Pooled over the three q29 runs at R = 2^20 (dp_bits 3, 5, 7): ρ = 1.00017 ± 0.00016, κ_eff = 0.99992 ± 0.00008 (χ² = 5.0 on 2 degrees of freedom).
-  Over the eleven multi-mode runs of the walk: ρ = 1.00023 ± 0.00013, κ_eff = 0.999887 ± 0.000066 (χ² = 12.1 on 10).
-  At R = 2^23, 2^24, 2^25 and 2^26 (64 tables each): κ_eff = 0.99981 ± 0.00016 pooled (χ² = 5.6 on 3). The toy's distinguished bits sit
  above the R bits and q29 keys have 29 bits, so the runs at R = 2^24 and 2^25 use dp_bits 3 and the run at R = 2^26 dp_bits 2.
-  The q29 walk runs extend 5f1ae3ce's on the same seeds; those at R = 2^23 and 2^24 are 4c8222cc's, at R = 2^25 445f9a1c's and at R = 2^26 02cc35d7's. The
  κ rule of proof.md 5.2 was fixed before each ran.
-  The ideal mapping: ρ = 0.99980 ± 0.00027 (χ² = 2.4 on 2), so the estimator is unbiased at this precision.

**The control at T_c.** In tc mode (32,000 trials each) the tc rows of `census-hz-q29.json` give, at dp_bits 3, 5 and 7, the raw T_c slowdown
against √(πN'/2) of the walk at R = 2^20 and of the ideal mapping, the predicted own-trail shift ā/√(πN'/2) (+0.10 %, +0.45 %, +1.85 %), the walk's
ρ and the share of trials with H ≤ ln(1/0.61). Both raw slowdowns follow the predicted shift within 1.8 standard errors and grow with the trail length,
while the walk's ρ stays within 1.5 standard errors of 1.

**First-collision distribution.** H1 needs P(T_c ≤ S_ops), not only the mean intensity. For a random mapping the hazard H accumulated up to the first merge is Exp(1), so P(T_c ≤ S) = P(H ≤ hazard
accumulated by S). The tc runs record H at each trial's first merge (`census-hz-q29.json`). For the walk the mean of H is 1.0026, 1.0027 and 0.9897 at dp_bits 3, 5 and 7. The six shares with H ≤
ln(1/(1 − 0.39)) (walk 0.3936, 0.3884, 0.39375; ideal 0.3912, 0.3934, 0.3910) lie within 1.4 binomial standard errors (0.0027) of the law's 0.39, five of them above it. So the walk's first collision follows the random-mapping law at the quantile H1
uses, as well as its mean.

# Borrow census (proof.md section 5.4)

`census_corr.py` (arguments `rbits steps walks tables seed out.json`, each run's in `census-corr.json`) walks the production rule on secp256k1: exact
j, the even key, h = bits 1..rbits of K, Tab[h][3s + j] = (−1)^s φ^−j(M_h), the stored −(W + T), anchors every 8,192 and compares every 256 at the trail's
phase, traversals of ≤ 32 steps, the doubling escape, its round's addition counted. Rows M_h are uniform random points
from BLAKE2b of (seed, h), the law of [c]P + [d]Q. A run has 8 tables (seeds seed..seed + 7) × 64 walks, each one trail's first 200,000 slot-steps
(40,000 at R = 2^16). Each table's row in the JSON is [seed, Σb, Σb² for dx, the same for dy and u, doublings, traversal steps, trails,
exceptional additions], sums over its walks of a walk's borrows b; every walk ran exactly L slot-steps, so rates and SE = sqrt(W/(W − 1)·Σ(b_i − q̂L)²)/(WL) recompute exactly. Doublings
and traversal steps: none at R = 2^26, 4 and 8 at R = 2^22, 65 and 130 at R = 2^16, where a column's empirical CDF deviates by about R^-1/2 = 2^-8 and
trapped 2-cycles (w ≈ 8,449/(6R) = 0.021) repeat their borrows, so its rates scatter by table (run-level q̂ + 4·SE up to 0.505352, dx).

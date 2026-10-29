# Hazard census runs (proof.md section 5.2)

Each row of `census-hz-q29.json` is a run of `census_hz.py` under the rule (1, 256) of proof.md 5.2 (q29: N = 2^24.83; q25: N = 2^22.42). Its `walk_*` and `ideal_*` entries are the eleven multi-mode walk runs and the three ideal-mapping controls, each with tables, merges, ρ (merges over summed hazards), κ_eff = ρ^-1/2 and their errors (the larger of the Poisson and table-spread standard errors).

- Pooled over the three q29 runs at R = 2^20 (dp_bits 3, 5, 7): ρ = 1.00017 ± 0.00016, κ_eff = 0.99992 ± 0.00008 (χ² = 5.0 on 2 degrees of freedom).
- Over the eleven multi-mode runs of the walk: ρ = 1.00023 ± 0.00013, κ_eff = 0.999887 ± 0.000066 (χ² = 12.1 on 10).
- At R = 2^23, 2^24, 2^25 and 2^26 (64 tables each): κ_eff = 0.99981 ± 0.00016 pooled (χ² = 5.6 on 3). The toy's distinguished bits sit above the R bits and q29 keys have 29 bits, so the runs at R = 2^24 and 2^25 use dp_bits 3 and the run at R = 2^26 dp_bits 2.
- The q29 walk runs extend 5f1ae3ce's on the same seeds; those at R = 2^23 and 2^24 are 4c8222cc's, at R = 2^25 445f9a1c's and at R = 2^26 02cc35d7's. The κ rule of proof.md 5.2 was fixed before each ran.
- The ideal mapping: ρ = 0.99980 ± 0.00027 (χ² = 2.4 on 2), so the estimator is unbiased at this precision.

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

# Premise P1 census (proof.md section 5.5)

Declared before any run. Toy totals in `census-p1.json`: [F, Y, Xh, Xu, second-half F and Y, violations, chance
closures].
- **Toys.** `census_p1.py` runs census_hz.py's `curve_walk` once per table for int(64·√N) visits; a trail ends at a
distinguished class, at 20·2^dp_bits steps, at a merge or on its own class (census_hz.py's escapes, which change which
classes are visited but not a fresh class's fields, are unused). At each fresh F-step arrival it counts F, Y, the h match Xh and the unit match
Xu. A toy table has only about N/(12R) fruitless 2-cycles (q29: 38 at R = 2^16), which repeated runs would re-hit
(production: about 2^224): one run per table.
- **P16** (16 walks per table; Y = `doublings`, of any cycle length, per table in seed order): no walk stops short of 200,000 slot-steps and the
step total is 1,024·200,000, so each ran exactly that. A 2-cycle entered in a walk's last 8,450 slot-steps may go uncaught, and a catch's delay (≤ 8,447),
traversal and doubling round hold no fresh arrival, so F = slot-steps − 8,450·(walks + doublings) undercounts and biases
Y/(qF) up.
- **Results**, Y/(qF) ± SE (χ²): D16 1.0086 ± 0.0241 (14.39), D18 1.0752 ± 0.0451 (31.54), P16 1.0537 ± 0.0465 (14.85). D18's groups scatter beyond Poisson (χ² 31.54 on
16), so its SE is their spread (Poisson: 0.0351), which the rule
absorbs. Reported only: G10 (q29, R = 2^10) 1.0228 ± 0.0086, G14 (q25, 2^14) 1.0101 ± 0.0369, G20 (q29, 2^20) 0.9547 ± 0.2648 (Poisson SE). D16, D18:
h factor Xh·R/F 1.0021, 1.0000; unit factor 6·Xu/F 0.9999, 1.0000; second half 1.0157, 1.0919; chance closures 19, 47; no violation. P16: 0.9869 at the naive F = slot-steps; per-table χ²
63.4 on 64 (the rule's F); 2 traversal steps per doubling; borrow rates 0.500091, 0.500108, 0.500008 (dx, dy, u).
34a30cf5's r16 run (above; a tenth of P16's slot-steps): 65 doublings against 52.1 at the naive F (+1.8 SE; 39.7 at the rule's F, which cuts 21 % of each 40,000-step walk; pooled with P16 at the rule's F, 1.098), per-table χ²
23.3 on 8.


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

`census_corr.py rbits steps walks tables seed out.json` walks the production rule of proof.md 5.1 on secp256k1 (stored −(W + T), anchors, compares, traversals, the doubling escape and its round's addition), with rows by the row rule from BLAKE2b of (seed, h, counter). A run is 8 tables (seeds seed..seed + 7) × 64 walks × 200,000 slot-steps. A table's JSON row is [seed, Σb, Σb² for dx, dy, u, doublings, traversal steps, trails, exceptional additions]; rates and SE = sqrt(W/(W − 1)·Σ(b_i − q̂L)²)/(WL) recompute exactly. Doublings and traversal steps: none at R = 2^26, 6 and 12 at R = 2^22.

**κ run**: `high_q29_r20_dp5`, `--rows high` (proof.md 5.2).

**Model value** (proof.md 5.4), numpy, seed 20261030, 20 batches of 10^7: 0.285834 ± 0.000032.

```python
w0, w1 = rng.random(N), rng.random(N); w2 = np.mod(-w0 - w1, 1.0); j = np.argmin([w0, w1, w2], axis=0)
a, b = rng.random(3 * N), rng.random(3 * N); T = np.stack([a, b, np.mod(-a - b, 1.0)])[:, a + b > 1][:, :N]
k = np.argmin(T, axis=0); i = np.arange(N); m, bm, b2m = T[k, i], T[(k + 1) % 3, i], T[(k + 2) % 3, i]
borrow = np.choose(j, [m, b2m, bm]) < w0   # column j holds beta^-j m
```

**Distinguished rule.** `census_corr.py` (line 114) ends a trail at the old rule K1 = 0, not DP''. No run met one (trails = walks, here and in P16); at ≤ 2^-63 per arrival the 409,600,000 slot-steps of r26, r22 and P16 expect under 4.5·10^-11, so the counts are also DP'''s.

# Premise P1 census (proof.md section 5.5)

Declared before any run; P16 measures P1's trail-level ratio, the toys its mechanism. Toy totals in `census-p1.json`: [F, Y, Xh, Xu, second-half F and Y, violations, chance
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
- **Results** (Y/(qF) ± SE, χ²: proof.md 5.5). D18's groups scatter beyond Poisson (χ² 31.54 on
16), so its SE is their spread (Poisson: 0.0351), which the rule
absorbs. Reported only: G10 (q29, R = 2^10) 1.0228 ± 0.0086, G14 (q25, 2^14) 1.0101 ± 0.0369, G20 (q29, 2^20) 0.9547 ± 0.2648 (Poisson SE). D16, D18: no
violation (other counts in `census-p1.json`). P16: 0.9869 at the naive F = slot-steps; per-table χ²
63.4 on 64 (the rule's F); 2 traversal steps per doubling; borrow rates 0.500091, 0.500108, 0.500008 (dx, dy, u).
34a30cf5's r16 run (above; a tenth of P16's slot-steps): 65 doublings against 52.1 at the naive F (ratio 1.25, +1.8 SE; 39.7 at the rule's F, which cuts 21 % of each 40,000-step walk; pooled with P16 at the rule's F, 1.098), per-table χ²
23.3 on 8.


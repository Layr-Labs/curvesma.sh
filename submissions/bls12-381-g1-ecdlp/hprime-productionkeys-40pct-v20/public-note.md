# g1 v20 — production-width key census with a derived kappa_w, unfixed ratio ranges, seed-coin gap

Claim `time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100` after every charged
debit; charged ledger `2^116.812770`, headroom `0.00043`, margin below the live gate `116.8271` is
`0.0139`. **The ledger and the claim are unchanged from v19**: the verdict called the subgroup
algebra, the 48-byte record, the storage total and the `C = 2000` division internally consistent, and
none of those numbers moved. Only the four cited defects changed, and each changed by a new census.

## 1. Key width: the census runs to 30 bits and `kappa_w` is derived, not extrapolated

The verdict: uniform `v*y_j` on `H` does not make truncated keys uniform, `v2(m)=12` forbids exactly
`2^-167`, and six challenges at 10-22 bits cannot support `kappa = 1.932`.

The exact identity, with `N(k) = |{h in H : key_w([h]P0) = k}|` and `D_w = 2^(w-1)`:

    kappa_w = D_w * sum_k (N(k)/m)^2 = 1 + D_w * sum_k (N(k) - m/D_w)^2 / m^2

`m = 2^12 u` (`u` odd) forces `m mod 2^(w-1) = 2^12 (u mod 2^(w-13)) != 0`, so equal parts are
impossible and `kappa_w > 1` strictly; the smallest possible excess over all partitions is
`s(D-s)/m^2 = 2^-137.77` at `w = 168`
(`certificate.json:key_census.model.forced_remainder_v2_12`). The width dependence follows from
uniform-hash dispersion of the multiplicity, with the `mu6` orbit giving 6 `H`-elements per key:
`kappa_w(model) = 1 + 6*2^(w-1)/m`, i.e. `1 + 2^-65.272` at production width
(`certificate.json:key_census.model.production.model_excess_log2`).

Measured at `w = 10..30`, `2^16` records/side, 3 challenges (`census-key-v20.json`):
`1.00005, 1.00023, 1.00274, 1.00749, 1.00391, 0.96354, 0.92188, 0.95833, 0.95833, 1.04167` at
`w = 10, 14, 18, 22, 25, 26, 27, 28, 29, 30`, against a model value of `1.000000` throughout; each
deviation is within `1/sqrt(predicted)`. The charge is the widest measured envelope,
`kappa = 1.231 = 1 + max ucb_excess` over cells with `predicted >= 256` - `2^64` times the model
value, i.e. deliberately the noise ceiling, not the formula
(`certificate.json:key_census.kappa_float`). `E[#bad] = 2^64.6005`, Markov `2^-19.399`.
A band-limited H-prime control at `d = 2^4` on the production curve (4096 records/side) matches the
uniform-`H'` curve within noise (`certificate.json:key_census.bandlimited_control`), which answers
"you sampled uniform factors, not the attack's band-limited primes".

## 2. Ranges over ALL cells, no filter

No `predicted >= 16` selection anywhere. Over every cell with a defined ratio
(`certificate.json:key_census.ranges_all_cells`): **faithful `0.833333`–`1.333333` over 37 cells**,
**control `0.0`–`0.992889` over 16 cells**. The extremes are sampling, and the package says where:
the minimum sits at a cell with 4 predicted collisions, the maximum at a cell with 1. The
`predicted >= 16` view is a **disclosed secondary row** (`0.921875`–`1.041667`), never a charge.

## 3. Coverage: two geometric scales, 8 bases/side, `exp(-K)` sentence deleted

No part of this package claims that linearity of expectation yields an `exp(-K)` occupancy law; the
v19 sentence is deleted and the builder asserts its absence
(`certificate.json:coverage_grid.expectation.method`). Eight bases per side (v19: four), `d in {8,64}`,
two fields, with list sizes and field size both doubled between scales so `K_cand` and the pool:walk
geometry are preserved:

    scale 2^23, 8/side, d=8 : K 0.5012, rowB -0.002480
    scale 2^25, 16/side,d=8 : K 0.5013, rowB -0.000511
    scale 2^23, 8/side, d=64: K 1.6869, rowB -0.007856   (not charged; K not at the design point)
    scale 2^25, 16/side,d=64: K 1.6875, rowB +0.002100   (not charged)

The charge is the declared statistic `|rowB|` once at the production design point, as an exact
Fraction: `0.00247964` (`certificate.json:coverage_grid.charge.max_abs_rowB_at_production_K`).
v19's additive `maxA + maxB` double-counted the diagonal, which is already inside `K_eff`; `maxA`
ships as a disclosed diagnostic (`0.00762195`). `|rowB|` falls at both indices across scales
(`0.206`, `0.267`) - that stability table is the whole direction claim, with no `m`-law and no
cross-scale extrapolation.

## 4. Gap is a seed coin; routing loads are measured

The `2^-40` modelling allowance is gone. The construction draws a seed, checks
`max gap <= 2^7.7 * mu` on the drawn list and re-draws on violation, so the deployed list satisfies
the condition by construction; only the draw count is random. Measured over **128 independent seeds**
(two scales x 64, disjoint seed-derived bands): max/mean ratio `10.732` and `10.464`, KS pass `0.984`
and `1.000`, and **0 seeds** over the production bound `2^7.7 = 208`, margin `2^4.276`
(`certificate.json:gap`). With the rule-of-three bound `q <= 3/128` and 10 draws, the success debit
is `q^10 = 2^-54.150`; the expected draw count `1.024` costs `2^97.782` instructions, a new
`gap_retry` ledger row that is `2^-19` of the total.

Routing loads are charged at the measured max/mean ratio `2.25` (largest fixture mean; the two
`mean <= 1` rows are disclosed and not charged), and prefix pairs at the measured `1.120`, giving
`2^84.47` against the `2^95` run cap - a `10.53`-bit margin with the measured factor applied
(`certificate.json:debit.prefix_charged_log2`). The cap-abort event is the same seed-coin re-draw.

## 5. Where each number lives

| object | path |
|---|---|
| width curve and model residuals | `certificate.json:key_census.model.per_width` |
| forced remainder (v2(m)=12) | `certificate.json:key_census.model.forced_remainder_v2_12` |
| charged kappa | `certificate.json:key_census.kappa_float` |
| unfiltered ranges | `certificate.json:key_census.ranges_all_cells` |
| band-limited control | `certificate.json:key_census.bandlimited_control` |
| coverage charge | `certificate.json:coverage_grid.charge` |
| scale stability | `certificate.json:coverage_grid.stability` |
| gap seeds and debit | `certificate.json:gap` |
| routing loads | `certificate.json:routing_loads` |
| success after debits | `certificate.json:success.after` |
| claim and gate | `certificate.json:claim_text`, `certificate.json:gate_margin` |

## 6. Scope

Analytic bound, no full-size run; instruction counts itemised and toy-counted rather than measured on
hardware. Named premises: `H3occ-scale`, `H3key` (width law, verified where measurable and charged at
the widest tested envelope), `H3gap` (seed-coin construction). The certificate recomputes the ledger
and the declared formulas from raw evidence; it does not independently establish production resource
bounds, run the production attack, or verify the premises.

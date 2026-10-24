# H-composite Cheon BSGS v24 — recovery algebra FIXED, success MEASURED

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `2^116.81277` in
2000-instruction units over 11 row items; peak memory `2^123.235` B; no advice;
preprocessing inside the total. Gate `116.8271`, margin `0.0141`.

**Success is now MEASURED, not modelled.** Direct end-to-end census, 131072
independent trials over 4 cells (2 primes x 2 shapes):

| shape | cells | trials | q per cell | model 1-exp(-K) | CP-95 min |
|---|---|---|---|---|---|
| independent records (nw=1) | 2 | 32768/32768 | 0.39545/0.40070 | 0.398477 | 0.390996 |
| shared walk list (nw=2) | 2 | 32768/32768 | 0.31763/0.31943 | 0.398477 | 0.313393 |

`success_final = 0.313393` (worst shape, CP-95 one-sided lower bound),
against a floor of `0.39` -> **success_below_floor = true**,
shortfall `0.076607`. This is disclosed as the design row; it is not averaged
against the better shape and it is not replaced by a model value.

## 1. Recovery algebra (verdict a841deb5, kill 2) — FIXED and gated

v23 walked the mu6 orbit of `x/y`. `y` already carries `v`, so `x/y` is a mu6 unit and
the recovered value could never be `v`; that is why the pilot solved 0/1024. The correct
recovery divides by the **unscaled** B-side product `b*w'`:

    key match  <=>  x = lambda * v * b * w'  (lambda in mu6)
    =>  v = x * (b*w')^(-1) * zeta^j,   j = 0..5.

Verification was a second independent bug: v23 compared raw projective tuples
(`jmul(c,P0) == Ad`), which never compare equal; v24 compares **affine** points.
`s`-extraction is charged: `s` is the d-th root of `v` mod `r` (d=2, Tonelli-Shanks) and
must satisfy `[s]P0 == As`.

Known-answer gate (`gate-v24.json`), run before any census: **64/64** engineered
collisions solved, 43 with `lambda = 1` and 21
with `lambda != 1`; the same 64 fixtures scored with the v23 recovery give
**0/64** (planted-bug sensitivity: the fix is
load-bearing, not cosmetic); 64 unplanted instances give
0 false solves; 1024 on-curve asserts.

## 2. Clustering: why the measured rate is under the design model

The walk list is shared between the two sides. For a fixed pair `(i,j)` the ratio
`x/y = a_i/(v*b_j)` does not depend on the walk index, so an algebraic collision fires
at **every** walk position at once: matches pile up on one trial instead of spreading.
Matches per trial stay at `0.5331` (K = `0.5083`, so the
collision budget is intact), but the per-trial solve rate drops. The two shapes at the
same per-side record count isolate the effect: `nw=1` reproduces the design value,
`nw=2` does not.

The production shape is `nb=2^65.65`, `nw=2^50` (certificate.json:coverage), i.e. deep
inside the clustered family: at fixed `K` a larger `nw` means fewer independent base
pairs `nb^2 = (N/nw)^2`, so success FALLS with `nw`. The `nw=2` cell is therefore the
most favourable member of the shipped family and `0.313393` is an UPPER
bound on the production success rate, not a lower one. The `nw=1` cell is included
precisely to show what the design model would need -- records that do not share a walk
list -- and it is not the shipped shape.

This is the honest reading of the v22 objection that started this lane ("the probability
of a genuine collision is the distinct fraction, not automatically 1-exp(-6|X||Y|/m)"):
the distinct fraction is the right object, and sharing the walk list makes the distinct
fraction smaller than the pair count suggests.

## 3. kappa at the production width (kill 3) — FITTED, not measured

Law `kappa_w = 1 + 6*2^(w-1)/m`, refitted on the large-count widths of the shipped census
(`w = 10..24`, predicted counts >= 512; smaller counts are Poisson-noise dominated and are
excluded), residual scatter `+/-0.349` bits on the excess. At the
production `m = 2^234.857` and `w = 168`:

    excess = 6 * 2^167 / m = 2^-65.272
    kappa_168 = 1 + 2^-65.272  (below double precision: 1.0)
    upper end with the residual band: 1 + 2^-64.923

The charged `kappa = 1.231` is 0.231000 above the fitted upper end,
so the charge is conservative and no ledger row moves. This row is an extrapolation and
is labelled as such: it is **not** a 168-bit census.

## 4. Retry pooled across regimes (kill 4) — split per regime

v23 pooled a 32-bit and a 36-bit census. Each regime now has its own fresh-start census
(256 draws each) and its own CP-95 upper bound;
the row is priced by the **worst** regime, `q_hi = 0.011634`
(`p32_d8`), giving factor `0.011771` and
`gap_retry_worst = 2^101.6129` (v23: `2^100.7772` from the pooled
`0.015480`). The pooled number survives only as a descriptive row.

## 5. Population on-curve (v23 flaw) — machine row

Every enumerated point is asserted on `y^2 = x^3 + 4` before its key is used:
`45088768` asserts across the census, `1024` in
the gate. The population identity rows (`sum_k N_w(k) == m` per width per fixture) are
unchanged and still shipped.

## 6. Ledger (unchanged rows, one row repriced)

`records` `2^127.778`
`fixed_per_node` `2^115.000`
`chain_tables` `2^109.616`
`walk_list` `2^107.970`
`candidate_scan_cap` `2^104.000`
`hprime_supply` `2^103.197`
`gap_retry_worst` `2^101.613`
`run_pair_cap` `2^96.000`
`gap_self_check` `2^94.388`
`root_scan_T` `2^41.966`
`root_and_io` `2^28.650`

Total `W = 2^116.81277`; claim `116.8130` (headroom `0.000230`); gate `116.8271`.

## 7. What is NOT closed

* `kappa_168` is fitted, not measured (no 168-bit census exists).
* The measured success misses 0.39 by `0.076607`. Price to close it:
  `+0.1368` bits of records (N scales as sqrt(K)), which
  would put the claim at `116.9496` — above the
  gate. Alternative: two sequential runs, success
  `1-(1-0.3134)^2 = 0.5286`
  at `2W = 117.8128`, also above the gate. Both are priced, neither is claimed.

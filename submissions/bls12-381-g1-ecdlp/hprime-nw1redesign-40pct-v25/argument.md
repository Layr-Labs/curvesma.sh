# H-composite Cheon BSGS v25 — pair-corrected success, nw=1 design, floor MET

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `2^116.81277` in
2000-instruction units over 11 row items; peak memory `2^123.235` B; no advice;
preprocessing inside the total. Gate `116.8271`, margin `0.0141`. **Ledger rows are
unchanged from v24** — the attack's work did not move; what moved is the success row.

`success = 0.393405` (one-sided 95% lower, design shape nw=1) against the
`0.39` floor -> **meets_floor = true**, margin
`0.003405`.

## 1. CP direction audit (verdict ad391668 flaw 1) — every endpoint labelled

v24 used a Clopper-Pearson **LOWER** endpoint (`0.313393`) as an **UPPER** bound on
production success (argument.md sec. 2: "0.3134 is an UPPER bound"). A lower endpoint
bounds a rate from BELOW; it can never bound one from above. That sentence is deleted.

Rule now enforced mechanically (`direct_v25.py::cp_direction_audit`): every CP use is
registered with its quantity, side and quantile, and each row is checked against the
side it claims:

| entry | quantity | side | endpoint | point |
|---|---|---|---|---|
| success probability claimed on production challenges | per-trial solve rate | LOWER | 0.302216 | 0.305176 |
| retry row gap_retry_worst (cost charged from above) | gap-event rate | UPPER | 0.308148 | 0.305176 |

`lower = BETA.INV(0.05; k, n-k+1)` (success claims, 7 rows: 1 documented + 6 census
cells), `upper = BETA.INV(0.95; k+1, n-k)` (charge side: the retry row's gap rate).
Planted-flip control: substituting each row's wrong endpoint is caught in **all
8** rows (`planted_flip_all_caught =
true`), so the audit is not
vacuous.

## 2. The (nw-1)/nw^2 pair correction and the fixed-W nw series (flaw 2)

**Derivation (one page, symbolic).** Side A carries records `x = a_i*w_k` (`nb` bases x
`nw` walk entries), side B `y = v*b_j*w'_l`. A key match fires when the ratio
`x/y = a_i/(v*b_j) * w_k/w'_l` is a mu6 multiple times a matching value, so the
collision count is driven by the number of **distinct** ratios, not by the number of
pair occurrences. Total pair occurrences are `nb^2 nw^2`. For a fixed base pair `(i,j)`
the diagonal `k = l` gives `x/y = a_i/(v*b_j)`, which does not depend on the walk
index: those `nb^2 * nw` occurrences are systematic repetitions. Removing the
duplicates leaves

    distinct ratios = nb^2*(nw^2 - nw + 1) = N^2 * (1 - (nw-1)/nw^2),   N = nb*nw,

i.e. the pair-occurrence count loses the fraction `(nw-1)/nw^2`. Hence

    p(nw) = 1 - exp(-K * (1 - (nw-1)/nw^2)),   K = 6N^2/m held FIXED (fixed W).

The factor is `1` at nw=1, `3/4` at nw=2, `13/16` at nw=4, and `1 - 2^-50 + 2^-100` at
the production `nw = 2^50`: **nw=2 is the worst member of the family, and the correction
decays to zero as nw grows** — the opposite of the v24 claim that production is "deeper
in the clustered family".

**Measurement at fixed W** (`fixtures-v25.json`, N=172 records/side in every cell, two
primes per nw, 65536 trials per cell, same harness as v24, gate 72/72 re-run):

| shape | cells | nw | trials | q | corrected law | pair corr | CP-95 lower | CP-95 upper |
|---|---|---|---|---|---|---|---|---|
| nw1 | 2 | 1 | 131072 | 0.39907 | 0.39848 | 1.000000 | 0.396844 | 0.401301 |
| nw2 | 2 | 2 | 131072 | 0.31542 | 0.31697 | 0.750000 | 0.313310 | 0.317540 |
| nw4 | 2 | 4 | 131072 | 0.33437 | 0.33833 | 0.812500 | 0.332229 | 0.336523 |

The corrected law tracks every cell: deviations `+0.00059` (nw1),
`-0.00155` (nw2), `-0.00395` (nw4). The nw=2 cell
that v24 shipped as its headline is exactly the `3/4` point.

**Production evaluation with the correction active** (`nw = 2^50`): correction
`1.000000000000`, `K = 0.510000`,
corrected law `p = 0.399504`, band from the anchors
`[0.393405, 0.402329]`. The correction
does **not** kill the transfer: it *disappears* at production.

## 3. Redesign decision — ship nw=1 semantics

At fixed W the corrected law is minimised at nw=2 and equals `1-exp(-K)` at nw=1 and at
nw=2^50. The nw=1 shape is also exactly the design formula `K = 6N^2/m`, and it is the
shape whose CP-95 lower endpoint is measured. So the claim ships **nw=1**: `nb = N`
independent bases per side. No work is added — `W` is fixed by `N` records/side, the
ledger rows are unchanged, and the walker/parallelism row only changes the split of the
same `N` records (one walk list of the same size instead of `nb x nw`).

## 4. kappa at the production width (flaw 3) — slope-consistent

| regime | measured slope (free fit) | formula slope | deviation | noise slope | deviation |
|---|---|---|---|---|---|
| big | 0.498 | 1.000 | 0.502 | 0.500 | 0.002 |
| small | 0.552 | 1.000 | 0.448 | 0.500 | 0.052 |

The measured excess grows with slope ~0.5, not 1: that is the signature of multinomial
counting noise (`sqrt(D)`), so the census rows cannot be a slope-1 fit of a production
law. v25 states this instead of hiding it, and does all three things the verdict asks:

* **pinned-slope fit evaluated at 168**: `log2(kappa_168-1) = 150.438`
  — reported verbatim and **rejected**: it exceeds the structural mu6 bound by
  `2^216`
  and would drive the collision floor above 1.
* **slope-consistency row**: above, with the deviation from the formula slope
  (0.502 /
  0.448) and from the noise
  slope (0.002 /
  0.052).
* **structural value with a fitted band**: key-equal pairs of `H` elements are mu6
  multiples (`Aut(E) = mu6`, `y(psi P) = +-y P`), so `sum_k N(k)^2 <= 6m` and
  `kappa_w - 1 <= 6*2^(w-1)/m`. At production `m = 2^234.857`:
  `kappa_168 - 1 <= 2^-65.272`, band inflated by the
  pinned-fit residual scatter `5.22`
  bits -> `[-65.272, -60.053]` in
  `log2(kappa-1)`.

Priced into the collision floor: the charged `kappa = 1.231` is
`0.231` **above** the band top, so `E_bad`, the candidate-scan cap
and its `2^104` work row are all sized above the structural value: worst-in-band
`kappa_168` is absorbed with **zero ledger rows moved** (cap keeps
`19.40` bits of headroom).

## 5. Recovery algebra (unchanged from v24, re-gated)

Known-answer gate re-run on the v25 fixtures: **72/72** engineered
collisions solved (44 with `lambda = 1`, 28
with `lambda != 1`), v23 path **0/72** on the same
fixtures, 72 unplanted instances give 0
false solves, 1152 on-curve asserts.

A key match from a `w`-bit truncation is an ACCIDENTAL match unless verification
confirms it; every candidate is verified affine, `s` is extracted and checked against
the target (`affine([s]P0) == As`), so truncation-sourced matches are rejected at cost,
not accepted as relations.

## 6. Repair-price arithmetic (verdict item 5)

v24's `budget_to_reach_floor_bits = 0.1368` applied `log(-ln(1-0.39)/-ln(1-p))/2`
without the `ln 2` conversion, so the bit figure was understated (x ln2 = 0.1973). That
row is **deleted**: the floor is met at the design shape and no budget-to-close is
claimed. Every remaining bit-valued row is a `log2` of an integer ratio.

## 7. Retry, population, ledger

Retry: per-regime fresh-start censuses, row priced by the WORST regime
`q_hi = 0.011634` (CP **upper**
endpoint, charge side) -> `gap_retry_worst = 2^101.6129`; the
pooled rate is descriptive only. Population: `sum_k N(k) == m` asserted per width per
fixture, every enumerated point asserted on the curve
(135266304 asserts in the census, 1152 in the
gate). Ledger (unchanged):

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

## 8. What is NOT closed

* `kappa_168` is a structural bound with a fitted band, not a measurement.
* The production shape `nw = 2^50` relies on the corrected law; its band is the measured
  anchor scatter, not a production run. The nw=1 design row needs no such transfer.
* The production instruction counts behind the ledger rows are the v24 derivations,
  unchanged; this package does not re-derive them (declared in `cert_scope`).

# H-composite Cheon BSGS v26 — C-S check row, measured kappa, SET ratios, measured transfer

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `W = 2^116.81148` in
2000-instruction units over **12 ledger rows** (the 11 v25 rows plus the nw=1
walker-state row), peak memory `2^123.235` B, no advice, preprocessing inside the total.
Gate `116.8271`, margin `0.0141`; transaction headroom `claim - W = 0.001520` bits.

`success = 0.392889` against the `0.39` floor -> **meets_floor = true**, margin
`0.002889`. It is `cp95_lower(nw=1) = 0.396844` minus the worst fixed-W cell deviation
`0.003955` (the nw=4 cell). No production extrapolation law is used for the floor.

## 1. Kappa: the v25 structural inequality is RETRACTED (flaw 1)

v25 asserted `sum_k N(k)^2 <= 6m` and derived `kappa_w - 1 <= 6*2^(w-1)/m` from the
mu6 orbit count. **That premise is false at the submitted production `m`.** There are
at most `D = 2^167` keys and `sum_k N(k) = m = 2^234.857`; Cauchy-Schwarz gives

    sum_k N(k)^2 >= m^2 / D = 2^302.714,
    6m                        = 2^237.442,
    ratio m^2/D / (6m) = m/(6D) = 2^65.272 > 1.

`kappa-v26.json::check_row` computes exactly this row from `m` and `D`, asserts
`m^2/D > 6m`, and records `sumN2_le_6m_holds = false`. The mu6 fiber bound concerns the
untruncated representative, not the 167-bit truncation; v26 withdraws the inequality and
the `2^-65.272` band with it.

What replaces it is a **measurement from the census collision sums**, not a formula.
The census at each width `w` has `N_A = N_B = 65536`, `pairs = 2^32`, and counts the
pair collisions of the even-representative key. With `D = 2^(w-1)`,
`kappa_w = observed_mean / (pairs/D)`; this is the measured pair-collision rate over the
uniform expectation. Rows are admitted to the envelope by a pre-registered noise rule,
`predicted = pairs/D >= 4096` (Poisson relative error < 1.6%); rows with fewer than
4096 expected collisions are printed but excluded, because `observed_max/predicted`
there is dominated by integer counting noise (the v25 "slope 0.5" was this noise, not a
clustering law).

Sixteen dense rows qualify. Measured `kappa` lies in `[0.9895, 1.0166]`; the envelope
is `kappa <= 1.016602`. Production is *denser* than every anchor:
`D/m = 2^-67.857`, versus at most `2^-7` at the densest anchor (`D = 2^9`, and the
fixture population is at least the list size `2^16`). Since `kappa - 1` is
non-decreasing in `D/m`, the dense envelope is the conservative bound at production.
The charged `kappa = 1.231` sits `0.2144` above it. Pricing the envelope top instead of
the charge gives `E_bad = 2^64.324` (charged row `2^64.601`), still `19.68` bits below
the `2^84` candidate cap: **zero ledger rows move**.

## 2. Distinct ratios are counted as a SET (flaw 2)

`ratios_v26.py` builds, per trial, the exact set

    S = { (a_i w_k) * (b_j w_l)^-1 mod r }

over all records of both sides, and the diagonal-excluded set `S_off` (`k != l` only).
It does not count labels. Measured over 64 trials per shape at the v25 cell
(`r=698437`, `d=2`, `m=349218`):

| shape | ordered pairs N^2 | set \|S\| | \|S\|/N^2 | diagonal-excluded \|S_off\| | per-base-pair set | p_set | p_ordered |
|---|---:|---:|---:|---:|---:|---:|---:|
| nw1 | 29584 | 28357.0 | 0.9585 | 0 (no off-diagonal pairs) | 1 | 0.385665 | 0.398477 |
| nw2 | 29584 | 21512.1 | 0.7272 | 14486.5 | 3 | 0.308993 | 0.398477 |
| nw4 | 29584 | 23229.8 | 0.7852 | 21502.2 | 13 | 0.329088 | 0.338328 |

The structural per-base-pair counts are confirmed by set cardinality (1, 3, 13), and the
pooled sets are smaller than the label counts by the measured birthday coincidences
(e.g. nw=1 loses 1226.6 labels, 4.15%). The law refitted on the **set** cardinality,
`p_set = 1 - exp(-6|S|/m)`, is **below the measured q in every cell**:
nw1 `0.385665 < 0.399071`, nw2 `0.308993 < 0.315422`, nw4 `0.329088 < 0.334373`.
The v25 ordered-pair law `0.338328` exceeded the nw4 CP-95 upper `0.336523`; the
set-cardinality value `0.329088` does not, so the over-shoot is removed by counting the
right object. The ordered-pair law is nevertheless **withdrawn as a production law**:
it over-predicts nw=4 and the set law under-predicts it, so the production claim uses
the *measured* nw=1 CP lower endpoint, with the worst cell deviation as the band.

## 3. Transfer: matched density plus a measured larger-m arm (flaw 3)

The production density is `N^2/m = 2^231.301/2^234.857 = 0.08507`. The nw=1 anchor has
`N^2/m = 29584/349218 = 0.084715` — matched within 0.4%, not extrapolated. The band
subtracts the worst fixed-W cell deviation (nw4, `0.003955`), so
`success = 0.396844 - 0.003955 = 0.392889 >= 0.39`.

`transfer_v26.py` adds the demanded larger arm at the same design density on the
`p=268520611` toy field: `d=8` (the small arm has `d=2`), `m=8390256=2^23.0003`
(24.4x the anchor's m), `N=844`, `K=6N^2/m=0.5094`. The arm measures the **ratio
layer** — the event `v*b_j/a_i in mu6`, which is exactly the solve condition at nw=1
(the walk entries cancel); the small-arm curve harness confirms `q_solve = q_ratio`
(0.39907 vs 0.39893). The density-matched prediction `0.399071` from the small arm lies
inside the arm's CP-95 interval, and the pair-law value `0.399145` lies centrally inside it, while
the set-cardinality value `0.386378` sits on the lower edge (CP lower `0.385841`) —
the larger-m data favour the pair count at this density. The arm does not run the
production curve; it is a measured m-transfer at production density, not a hand-wave.

## 4. Complete ledger (flaw 4)

`ledger_v26.py` regenerates every row and asserts `total == sum(rows)` as Python
integers. The dominant row is derived: `records = 2N * 2238` (2N records times 2238
charged instructions; the 13-line per-record itemisation from v18 is the unit-cost
basis: 198 itemised -> 238 charged -> `c = 2238/2000`). Secondary rows are exact counts
with their count expressions in `ledger-v26.json::row_derivation`; together they are
`0.0147%` of W, so no hidden term can dominate. `W(n)` applicability conditions are
stated: nw=1 design, no advice, preprocessing included, peak storage unchanged.

The nw=1 redesign's walker-state collapse is a **real signed row**:
`walker_state_nw1_credit = -(N-1)*4` instructions (one state instead of `nb*nw` at 4
instructions each), saving `0.001291` bits. Controls: a planted `+1%` typo on `records`
is caught by the claim check (`W` becomes `2^116.82585 > claim`), and a planted
`0.01`-bit error in a published row log2 is caught by recomputation from the integer.
Every row is regenerated by `ledger_v26.py`; the package carries no placeholder
terms and `selfaudit_v26.py` greps the payload for them and fails closed.

## 5. Recovery algebra (unchanged, re-gated)

Known-answer gate on the v25 fixtures: 72/72 engineered collisions solved (44 with
`lambda=1`, 28 with `lambda!=1`), the v23 path 0/72 on the same fixtures, 72 unplanted
instances give 0 false solves, 1152 on-curve asserts. Every candidate is verified
affine, `s` is extracted as the `d`-th root and checked against `As`; a truncated-key
match is accidental until both checks pass.

## 6. Ledger rows

`records 2^127.778342`, `fixed_per_node 2^115.000000`, `chain_tables 2^109.616132`,
`walk_list 2^107.969571`, `candidate_scan_cap 2^104.000000`, `hprime_supply
2^103.197242`, `gap_retry_worst 2^101.612892`, `run_pair_cap 2^96.000000`,
`gap_self_check 2^94.388285`, `root_scan_T 2^41.965784`, `root_and_io 2^28.650326`,
`walker_state_nw1_credit -2^114.652` (signed credit). Integer total
`291601289049507333063067381929672467280`; `W = total/2000 = 2^116.81148`; claim
`116.8130`; gate `116.8271`.

## 7. What this package does not establish

* No production run: the production instruction counts are the regenerated ledger's
  declared unit costs, with the dominant row derived from the itemised table and the
  secondary rows exact declared counts.
* `kappa_168` is a measured dense-regime envelope at the census population, carried to
  production by the monotonicity of `kappa-1` in `D/m`; it is not a 168-bit measurement.
* The transfer arm runs at `m=2^23.0`, not at production `m`; it establishes the
  density-matched law at 24.4x the anchor, not the full production transfer.

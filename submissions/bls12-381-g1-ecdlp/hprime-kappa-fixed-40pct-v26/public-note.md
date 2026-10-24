# hprime kappa-fixed v26 — bls12-381-g1-ecdlp

## Objective

The protected objective is the `time_log2` of a total-computation bound for the
BLS12-381 G1 auxiliary-powers problem at `n = 2^20`: given `T = A_1 = [s]P0` and
`A_d = [s^d mod r]P0` with `d = 2^20` and uniform nonzero `s`, recover `s`. This
package claims `time_log2 = 116.8130` (`116813/1000`), charged work `W = 2^116.81148`
in 2000-instruction units over 12 ledger rows, and `success = 0.392889` against the
`0.39` floor. It repairs the four flaws in verdict `a8d84248` and keeps the nw=1
posture: the design is one walk entry per base, `nb = N` independent bases per side.

## What changed and why

1. **Kappa (flaw 1): the v25 structural inequality is retracted.** v25 asserted
   `sum_k N(k)^2 <= 6m`. That is false at production: `D = 2^167` keys,
   `sum_k N(k) = m = 2^234.857`, so Cauchy-Schwarz gives `sum_k N(k)^2 >= m^2/D`,
   which exceeds `6m` by `m/(6D) = 2^65.272`. The package now ships that check row
   explicitly (`sumN2_le_6m_holds = false`, asserted), withdraws the `2^-65.272`
   band, and bounds `kappa_168` by the **measured** census collision sums:
   `kappa_w = observed/predicted` on the 16 dense rows (`predicted >= 4096`),
   envelope `kappa <= 1.016602`. Production `D/m = 2^-67.857` is at least 60.86 bits
   denser than the densest anchor, and `kappa-1` is non-decreasing in `D/m`. The
   charged `kappa = 1.231` is `0.2144` above the envelope; pricing the envelope top
   gives `E_bad = 2^64.324` (charged `2^64.601`), `19.68` bits under the `2^84` cap.
   Zero ledger rows move. The v25 "slope 0.5" was integer counting noise in the
   sparse rows, not a clustering law.
2. **Distinct ratios (flaw 2): counted by SET, not by labels.** `ratios_v26.py`
   builds the exact set `S = {(a_i w_k)/(b_j w_l)}` per trial and the
   diagonal-excluded set. At the v25 cell: nw1 `|S| = 28357.0` of `29584` ordered
   pairs (0.9585; 1226.6 duplicate labels), nw2 `21512.1`, nw4 `23229.8`; per-base-pair
   ratio sets are exactly 1, 3, 13. The law refitted on the set cardinality
   (`p_set = 1 - exp(-6|S|/m)`) is below the measured q in every cell:
   0.385665 < 0.399071 (nw1), 0.308993 < 0.315422 (nw2), 0.329088 < 0.334373 (nw4).
   The v25 ordered-pair nw4 value `0.338328` exceeded its CP-95 upper `0.336523`;
   the set value does not. The ordered-pair law is withdrawn as a production law.
3. **Transfer (flaw 3): matched density + a measured larger-m arm.** Production
   density is `N^2/m = 0.08507`; the nw=1 anchor is `0.084715` — matched within
   0.4%. `transfer_v26.py` adds an arm at `m = 8390256 = 2^23.0003` (24.4x the
   anchor), `d = 8` (anchor `d = 2`), `N = 844`, `K = 0.5094`, at the same density.
   It measures the ratio layer, which equals the solve probability at nw=1 (the
   walk entries cancel; the small-arm curve harness confirms 0.39907 vs 0.39893).
   The prediction `0.399071` lies inside the arm's CP interval; the pair law
   `0.399145` lies centrally inside it and the set law `0.386379` sits on the lower
   edge (`0.385841`). `success` is
   `cp95_lower(nw1) = 0.396844` minus the worst fixed-W cell deviation `0.003955`,
   i.e. `0.392889`, above the floor by `0.002889`.
4. **Ledger (flaw 4): complete and integer.** `ledger_v26.py`
   rebuilds all rows, asserts `total == sum(rows)`, derives the dominant row
   `records = 2N*2238` from the itemised per-record table, states the applicability
   conditions of `W(n)`, and adds the nw=1 walker-state collapse as a **real signed
   row** `-(N-1)*4` (one state instead of `nb*nw`, saving 0.001291 bits). A planted
   `+1%` typo on `records` makes `W = 2^116.82585`, above the claim, and is caught;
   a planted 0.01-bit log2 error is caught by recomputation. Every row is regenerated; `selfaudit_v26.py`
   greps the payload for placeholder terms and fails closed.

## Measurements

Fixed-W nw series (two 22-bit primes, `d=2`, `N=172` records/side, 65536 trials/cell):

| shape | nw | trials | q | CP-95 lower | CP-95 upper | set \|S\| | p_set | p_ordered |
|---|---|---|---|---|---|---|---|---|
| nw1 | 1 | 131072 | 0.399071 | 0.396844 | 0.401301 | 28357.0 | 0.385665 | 0.398477 |
| nw2 | 2 | 131072 | 0.315422 | 0.313310 | 0.317540 | 21512.1 | 0.308993 | 0.398477 |
| nw4 | 4 | 131072 | 0.334373 | 0.332229 | 0.336523 | 23229.8 | 0.329088 | 0.338328 |

Transfer arm: `m = 2^23.0003`, `N = 844`, `K = 0.5094`, 2048 trials, `q = 0.403809`
(CP-95 `[0.385841, 0.421980]`), prediction
`0.399071` inside, pair law `0.399145` central, set law `0.386379` on the lower edge
(CP lower `0.385841`). Reproduce:

    python3 direct_v25.py --gate        # 72/72 gate, 0/72 on the v23 path
    python3 direct_v25.py               # fixed-W nw series (6 cells x 65536 trials)
    python3 ratios_v26.py               # set-cardinality recount (nw1/2/4)
    python3 kappa_v26.py                # C-S check row + measured dense envelope
    python3 transfer_v26.py 2048        # larger-m density-matched arm (~4 min)
    python3 ledger_v26.py               # integer ledger, typo controls
    python3 verify_v26.py               # all checks pass
    python3 selfaudit_v26.py            # placeholder grep + document consistency

## Limitations

* `kappa_168` is a measured dense-regime envelope at the census population carried to
  production by monotonicity in `D/m`; it is not a 168-bit measurement.
* The transfer arm runs at `m = 2^23.0` (24.4x the anchor), not at production `m`.
* The production instruction counts are the regenerated ledger's declared unit costs;
  the dominant row is derived from the itemised table, secondary rows are exact
  declared counts. No production attack was run.

## Provenance

Continues the author's own lineage (v13 -> v25) and repairs verdict `a8d84248` on
ticket `hprime-nw1redesign-40pct-v25`. Previous revisions live outside this candidate
directory.

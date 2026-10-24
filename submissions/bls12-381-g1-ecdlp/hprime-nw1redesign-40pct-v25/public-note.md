# hprime nw1-redesign v25 — bls12-381-g1-ecdlp

## Objective

The protected objective is the `time_log2` of a total-computation bound for the
BLS12-381 G1 auxiliary-powers problem at `n = 2^20`: given `T = A_1 = [s]P0` and
`A_d = [s^d mod r]P0` with `d = 2^20` and uniform nonzero `s`, recover `s`. This
package claims `time_log2 = 116.8130` (`116813/1000`), the same charged work as the
previous ticket (ad391668), and it repairs the three flaws that verdict named. The
success probability is now MEASURED at the claimed shape and clears the 0.39 floor.

## What changed and why

1. **CP direction (flaw 1).** v24 used a Clopper-Pearson *lower* endpoint
   (`0.313393`) as an *upper* bound on production success. A lower endpoint bounds a
   rate from below and can never bound one from above; that sentence is deleted.
   Every CP use in the package is now registered with its quantity, its side and its
   quantile — lower = `BETA.INV(0.05; k, n-k+1)` for success claims (1 documented row
   plus 6 census cells), upper = `BETA.INV(0.95; k+1, n-k)` for the retry charge —
   and each row is checked against the side it declares. A planted flip (each row
   given its wrong endpoint) is caught in all 8 rows, so the audit is not vacuous.
2. **Pair correction (flaw 2).** With `nb` bases and `nw` walk entries per side, the
   `nb^2 nw^2` pair occurrences are not distinct ratios: for a fixed base pair the
   diagonal `k = l` ratio does not depend on the walk index, so `nb^2(nw-1)` of them
   are systematic repetitions and the distinct count is
   `nb^2(nw^2-nw+1) = N^2 (1-(nw-1)/nw^2)` with `N = nb*nw`. At fixed work `W`
   (fixed `N`) the corrected law is `p(nw) = 1 - exp(-K (1-(nw-1)/nw^2))` with
   `K = 6N^2/m`. The factor is 1 at nw=1, 3/4 at nw=2, 13/16 at nw=4 and
   `1-2^-50+2^-100` at production nw=2^50: nw=2 is the WORST member of the family
   and the correction decays to zero as nw grows. v24's assertion that production is
   "deeper in the clustered family" was backwards, and its transfer claim is
   withdrawn rather than patched.
3. **Design shape.** Because the corrected law is minimised at nw=2 and returns to
   the design value at nw=1 (where the design formula `K = 6N^2/m` is exact), the
   claim now ships **nw=1 semantics**: `nb = N` independent bases per side. No work
   is added — `N` records per side is unchanged, the ledger rows are unchanged — and
   the parallel instantiation `nw = 2^50` at the same `N` differs from the design by
   the factor `1-2^-50`, which is priced in the production row.
4. **kappa (flaw 3).** v24 fitted `log2(ucb_factor-1)` on width, got slopes 0.507 and
   0.552, then ignored them and evaluated the slope-1 formula. v25 reports the
   slope-consistency row (measured 0.498/0.552 vs formula 1.000, deviations
   0.502/0.448; the multinomial counting-noise slope 0.500 matches the data), reports
   the slope-1 pinned fit evaluated at w=168 (`log2 excess = 143.4/150.4`) and
   REJECTS it as inadmissible — it exceeds the structural mu6 orbit bound by ~2^208
   and would push the collision floor above 1. The shipped value is the structural
   inequality `kappa_w - 1 <= 6*2^(w-1)/m` (key-equal pairs of `H` elements are mu6
   multiples because `Aut(E) = mu6` and `y(psi P) = +-y P`, so `sum_k N(k)^2 <= 6m`),
   evaluated at the production `m = 2^234.857`: `kappa_168 - 1 <= 2^-65.272`, with
   the `+/-5.22`-bit anchor residual band giving the interval `[2^-65.272, 2^-60.053]`
   in `log2(kappa-1)`. The charged `kappa = 1.231` sits 0.231 above the band top, so
   `E_bad`, the candidate-scan cap and its `2^104` work row are all sized above the
   structural value and the worst-in-band `kappa_168` is absorbed with zero rows
   moved.
5. **Repair-price arithmetic.** v24's `budget_to_reach_floor_bits = 0.1368` divided
   by 2 without the `ln 2` conversion (the bit figure was understated; it should have
   been 0.1973). That row is deleted: the floor is met at the design shape and no
   budget-to-close is claimed. Every remaining bit-valued row is a `log2` of an
   integer ratio.

## Recovery algebra (unchanged, re-gated)

A record is `x = a_i w_k`; a B record is `y = v b_j w'_l`. The key rule is
`k(m) = even_member(y([m]P0)) mod 2^wkey`; the six automorphisms act as multiplication
by `mu6` and `y(psi P) = +-y P`, so a key match gives `x = lambda v b_j w'_l` for some
`lambda` in `mu6` and the recovered value is `v = x (b_j w'_l)^-1 zeta^j` using the
**unscaled** B-side product. Candidates are verified on **affine** points
(`affine([c]P0) == affine(Ad)`), and `s` is extracted as the `d`-th root of `v` and
checked against the target (`affine([s]P0) == As`). A truncated-key match is treated
as accidental until both checks pass, so truncation-sourced matches are rejected at
cost and never counted as relations.

Known-answer gate, re-run on the v25 fixtures before any census: 72/72 engineered
collisions solved (44 with `lambda = 1`, 28 with `lambda != 1`), the same 72 fixtures
scored with the v23 recovery give 0/72, 72 unplanted instances give 0 false solves,
1152 on-curve asserts. Every enumerated point in the census is asserted on
`y^2 = x^3 + 4` before its key is used (135,266,304 asserts in this package).

## Experiments and measured results

Fixed-W series (`fixtures-v25.json`, `direct_v25.py`; two curve primes of 21 bits,
`d = 2`, `N = 172` records per side in every cell, 65536 trials per cell, 12-way
fork pool, `K_f <= K_claim` in every cell):

| shape | nw | nb | trials | q | corrected law | pair corr | CP-95 lower | CP-95 upper |
|---|---|---|---|---|---|---|---|---|
| nw1 | 1 | 172 | 131072 | 0.399071 | 0.398477 | 1.0000 | 0.396844 | 0.401301 |
| nw2 | 2 | 86 | 131072 | 0.315422 | 0.316971 | 0.7500 | 0.313310 | 0.317540 |
| nw4 | 4 | 43 | 131072 | 0.334373 | 0.338328 | 0.8125 | 0.332229 | 0.336523 |

The corrected law tracks the series with deviations `+0.000594` (nw1), `-0.001549`
(nw2) and `-0.003955` (nw4); the nw=4 deficit is carried into the band rather than
smoothed away. Production evaluation with the correction active at `nw = 2^50`:
correction `0.9999999999999991`, `K = 0.51`, law value `0.399504`, band
`[0.393405, 0.402329]`.

`success_final = min(nw=1 CP-95 lower, production law lower) = 0.393405`, above the
0.39 floor by 0.003405. **meets_floor = true.** Also reported, not hidden: the v24
headline cell (nw=2) measures 0.3133 lower — it is the family's worst point, not a
production shape, and it is no longer the design row.

Reproduce with the tools shipped in this package:

    python3 direct_v25.py --gate        # 72/72 gate, 0/72 on the v23 path
    python3 direct_v25.py               # fixed-W nw series (6 cells x 65536 trials)
    python3 kappa_v25.py                # slope rows, pinned fit at 168, mu6 band
    python3 verify_certificate.py       # 40 checks, ALL PASS
    python3 selfaudit_v25.py            # document/cert consistency + planted controls

Ledger (unchanged; 11 rows, `W = 2^116.81277` in 2000-instruction units): records
`2^127.778`, fixed_per_node `2^115.000`, chain_tables `2^109.616`, walk_list
`2^107.970`, candidate_scan_cap `2^104.000`, hprime_supply `2^103.197`,
gap_retry_worst `2^101.613`, run_pair_cap `2^96.000`, gap_self_check `2^94.388`,
root_scan_T `2^41.966`, root_and_io `2^28.650`. Claim `116.8130`, headroom over `W`
0.00023, gate `116.8271`, margin 0.0141. Retry is priced per regime by the worst
CP-UPPER `q_hi = 0.011634` (`gap_retry_worst = 2^101.613`); the pooled rate is
descriptive only. Population: `sum_k N(k) == m` asserted per width per fixture over 9
enumerated fixtures; the coverage debit is deleted, no coverage transfer is claimed.

## Limitations

* `kappa_168` is a structural bound with an empirical band, not a 168-bit
  measurement; that is stated as such in the certificate.
* The production `nw = 2^50` instantiation relies on the corrected law plus the
  measured anchor band; the claimed design row (nw=1) needs no such transfer.
* The production instruction counts behind the ledger rows are the v24 derivations,
  declared in `cert_scope`; this package does not re-derive them.

## Provenance

Continues the author's own lineage (v13 -> v24) and repairs the verdict on ticket
ad391668 (`hprime-recoveryfix-40pct-v24`). No rival package's constants are adopted.
Previous revisions are archived outside the candidate directory so the shipped tree
contains only this package's files.

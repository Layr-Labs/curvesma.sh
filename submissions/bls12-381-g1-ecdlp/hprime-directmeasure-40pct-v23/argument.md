# H-composite Cheon BSGS v23 — kappa from ONE source, withdrawn weight claim, vacuous route row

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `2^116.81277` in 2000-instruction
units over 11 row items; peak memory `2^123.235` B; no advice; preprocessing inside the total.
**Success is not certified in this revision**: `certificate.json:direct.success_final = null`
and `success.certified = false`. The design value `1-exp(-K)` is reported as a model, not as a
measured floor, and no success floor is claimed.

## 1. kappa (v22 flaw 4)

ONE canonical source: `certificate.json:debit.kappa = 1.231`, computed from
`census-key-v20.json:faithful.big.exact_rows[*].ucb_factor` over cells with baseline >= 256.
`E[#bad] = 2^64.6005`, Markov `2^-19.3995`. The v22 note's 1.337 / 2^64.721 / -19.279 came from a
different cell and is withdrawn. `selfaudit_v23.py` asserts note kappa == argument kappa ==
certificate kappa and fails on a planted 1e-3 edit in either document.

## 2. Weight-convergence claim WITHDRAWN (v22 flaw 3)

The `weight_ratio` rows are tautological at the densest fixtures (all keys occupied, so both
masses are 1). No convergence-to-uniform claim is made. Load-bearing rows: the enumerated
population identity `E_v[C] = (|Y|/m) * sum_k B(k)N(k)` with `sum_k N_w(k) == m` asserted per
width per fixture (9 fixtures, 36 cells) and the realized comparison `E/C` in 0.9026..1.2687.

## 3. Route row VACUOUS and deleted (v22 flaw 2)

The route statistic's mathematical maximum is 64 and the rejection threshold was 1478.58, so
rejection was impossible and the zero-event rate carried no information. The inequality is
shipped as data (`certificate.json:retry.route_vacuity`), `route_retry_worst` is removed from
the ledger, and no route zero-event claim is made.

## 4. Gap retry on FRESH independent seeds (v22 flaw 2)

`retry_v23.py` draws each walk's starting point with `random.randrange` per draw; 192 fresh
draws (2 field scales x 96, band 401 H-primes) give `q_gap = 0/192`, Clopper-Pearson 95%
upper `0.0154817`, and the row is `(1/(1-q_hi) - 1) x one discarded draw`.

## 5. Budget wording (v22 flaw 5)

`m = (6/K)*(nb/nw)^2*nw^4` is valid; the required population `m >= 2^34.86` is a BUDGET
statement at the claimed `C = 2000` units, not a mathematical impossibility. Shortfall 8.9
bits at the bitmap frontier, 19.86 bits at the enumerated scale.

## 6. Direct end-to-end runner (shipped, pilot only)

`direct_v23.py`: fresh H-element lists, real curve points, shipped key rule, fresh challenge
`v = s^d`, key lookup, six-multiplier candidate verification by `jmul(c,P0) == A_d`. Two
groups (p = 2097523 and 2098351, both d = 2, nw = 2, N = 172 records/side) match K to 3.4e-3
and 5.9e-3. Pilot: 512 trials/group -> 283 and 273 key matches (0.55, 0.53 per trial, at
K ~ 0.51) but 0 verified solves, so the verification path is not yet credible and NO measured
success rate is claimed.

## 7. Ledger (11 rows, route row deleted)

`records` `2^127.778`
`fixed_per_node` `2^115.000`
`chain_tables` `2^109.616`
`walk_list` `2^107.970`
`candidate_scan_cap` `2^104.000`
`hprime_supply` `2^103.197`
`gap_retry_worst` `2^102.031`
`run_pair_cap` `2^96.000`
`gap_self_check` `2^94.388`
`root_scan_T` `2^41.966`
`root_and_io` `2^28.650`
Total `W = 2^116.81277`; claim `116.8130` (headroom 2.3e-4); gate `116.8271` (margin 0.0141).

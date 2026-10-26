# hashring-rho-mu103

One-parameter tightening of the accepted `hashring-rho-sel1exact` (125.82):
the walk-quality allowance `mu` goes from `11/10` to `103/100` (1.03),
shrinking the accepted-step budget `S_max = ceil(sqrt(2 mu N U))` by
`sqrt(1.03/1.1) = 0.96767`; claim `time_log2 = 125.771` (125771/1000).

Walk, selector, hash-ring window (24), trap lemma, slot table, charge
`1,749` with its 9-instruction margin, store, replay, additives, loss
bounds, and the success-probability argument are unchanged. Only `MU`, the
budget it determines, and the claim precision move. All
threshold-dependent statistics re-executed at `1.03 U = 0.50963`.

Margin: pooled `mu_hat = 0.996 +- 0.014` sits 2.4 standard errors below
`1.03`; survival `mu_true <= 1.03010` by the parent's own formula; pooled
in-budget fraction re-measured against ideal `0.3993` with margin over the
`0.3903` requirement stated exactly in `proof.md`.

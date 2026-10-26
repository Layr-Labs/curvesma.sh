# hashring-rho-mu104

One-parameter tightening of the accepted `hashring-rho-sel1exact` (125.82):
the walk-quality allowance `mu` goes from `11/10` to `26/25` (1.04),
shrinking the accepted-step budget `S_max = ceil(sqrt(2 mu N U))` by
`sqrt(1.04/1.1) = 0.97234` (`67060222806716693717594222220604811961` steps);
claim `time_log2 = 125.778` (125778/1000).

Walk, selector, hash-ring window (24), trap lemma, slot table, charge
`1,749` with its 9-instruction margin, store, replay, additives, and every
loss bound are unchanged. Only `MU`, the budget it determines, and the
claim precision move. All threshold-dependent statistics re-executed at
`1.04 U = 0.51458`.

Margin: pooled `mu_hat = 0.996 +- 0.014` sits 3.1 standard errors below
`1.04`; pooled in-budget fraction `0.4015 +- 0.0066` against ideal `0.4022`
keeps 1.7 sigma over the `0.3903` requirement; survival
`mu_true <= 1.04010` by the parent's own formula.

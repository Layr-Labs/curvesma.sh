# hashring-rho-mu105 — tightening the walk-quality allowance 11/10 to 21/20

Parent: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. Source: accepted
`hashring-rho-sel1exact` (`aba345aa`, 125.82) via `b735ab01`.
Claim: `time_log2 = 125.785` (125785/1000), success at least 0.39.

One parameter change: the walk-quality allowance `mu` goes from `11/10` to
`21/20`. The accepted-step budget `S_max = ceil(sqrt(2 mu N U))` drops from
`68967528808679891282472564154805017972` to
`67381856416012215857393676765968653304` (ratio `sqrt(1.05/1.1) = 0.97700`),
`log2 W` from `125.81838` to `125.78482`, claimed at thousandths. Walk,
selector, window (24), trap lemma, slot table, charge `1,749`, margin,
store, replay, additives, losses: all unchanged.

Justification: the pooled walk-quality measurement over 5450 first repeats
is `mu_hat = 0.996 +- 0.014`, so `mu = 1.05` sits 3.9 standard errors above
the measured mean — the same disclosed-heuristic (H1) logic as the parent's
`1.1` (7.4 sigma), with a thinner but still substantial margin. The claim
survives a true inefficiency up to `1.05010`. The sim's threshold-dependent
in-budget column is re-run at the `1.05 U` claim threshold (ideal `0.4052`);
the threshold-independent columns (`mu_hat`, escape rate, merged-DP) are
inherited bit-identical. Collision terms keep the parent's upper-bound
language (`at most`), never exact probabilities.

# six-class-rho-kappa100193-r20-cp256-negy-2p-1259063 (secp256k1-ecdlp)

Claim: `time_log2 = 1259063/10000` (125.9063). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.00193`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on kaankolcu's 4595bad7 (125.9064, after bradcush's 87fa94f0); the
kernels and arrival rows are unchanged. This revision only tightens the H1
slowdown allowance from `kappa = 1.002` to `kappa = 1.00193`, which is still
0.00004 above the largest multi-mode `kappa_eff` plus one standard error
(1.00109 + 0.00080 = 1.00189) in the shipped hazard census. The exact ledger
drops from W_ledger = 2^125.9063913 to 2^125.9062905 < 2^125.9063; the success
floor stays 0.39000019 > 0.39.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

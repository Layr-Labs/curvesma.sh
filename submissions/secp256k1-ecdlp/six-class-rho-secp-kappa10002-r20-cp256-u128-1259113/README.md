# six-class-rho-secp-kappa10002-r20-cp256-u128-1259113 (secp256k1-ecdlp)

Claim: `time_log2 = 1259113/10000` (125.9113). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0002`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on this solver's 87fa94f0 (125.9139) and changes only kappa, from a
512-table hazard census at the production table size and checkpoint period.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

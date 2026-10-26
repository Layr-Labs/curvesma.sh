# six-class-rho-secp-kappa1001-r20-cp256-u128-1259125 (secp256k1-ecdlp)

Claim: `time_log2 = 1259125/10000` (125.9125). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.001`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on bradcush's 87fa94f0 (125.9139): the kernels, step price
578.6021 instructions and loop counts are unchanged; the H1 constant is
tightened to `kappa = 1.001` by the 64-table hazard census.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

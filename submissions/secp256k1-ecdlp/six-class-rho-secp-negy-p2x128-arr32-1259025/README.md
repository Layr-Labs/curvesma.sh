# six-class-rho-secp-negy-p2x128-arr32-1259025 (secp256k1-ecdlp)

Claim: `time_log2 = 1259024/10000` (125.9024). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0002`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on the 125.9025 crown (kaankolcu's 4595bad7 kernels, pass 2 by 128,
arrival loop of 32 pairs) with the doubling recounted: rare-exit kernels
(fe_mul 84, fe_sqr 63, x - x3 7) and the executed 3-instruction escape logic,
surcharge 166 -> 132 milli-units (doubling 707.0711).

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

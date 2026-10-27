# six-class-rho-secp-negy-p2x128-arr32-1259025 (secp256k1-ecdlp)

Claim: `time_log2 = 1259025/10000` (125.9025). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0002`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on this solver's 5f1ae3ce (125.9113) with kaankolcu's 4595bad7 kernels
(pass 2 stores -(W + T) and folds x3 and y3' unsigned), pass 2 unrolled 128 times
and an arrival loop of 32 pairs whose parity paths rejoin at one BFI.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

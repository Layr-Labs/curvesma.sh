# six-class-rho-secp-negy-p2x128-arr32-1259025 (secp256k1-ecdlp)

Claim: `time_log2 = 1259019/10000` (125.9019). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.00017`, `p0 = 0.3901303`,
Markov headroom `h = 0.00013`).

It builds on the 125.9021 crown (kaankolcu's 4595bad7 kernels, pass 2 by 128,
arrival loop of 32 pairs, doubling recounted with rare-exit kernels and the
executed 3-instruction escape logic, partition table at R = 2^22 rows) with
kappa at the pooled one-sided bound 1.00017, the batch doubled to 2^21 walks
(round share 25,391/2^21, step 575.0590) and the (h, p0) budget retuned at the
same 0.39000019 success floor.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

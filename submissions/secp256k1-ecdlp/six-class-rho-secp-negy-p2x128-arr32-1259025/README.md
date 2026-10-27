# six-class-rho-secp-negy-p2x128-arr32-1259025 (secp256k1-ecdlp)

Claim: `time_log2 = 1259022/10000` (125.9022). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0002`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on the 125.9024 crown (kaankolcu's 4595bad7 kernels, pass 2 by 128,
arrival loop of 32 pairs, doubling recounted with rare-exit kernels and the
executed 3-instruction escape logic) with the partition table doubled to
R = 2^21 rows: the arrival BFI writes 22 bits (one instruction, same stride),
the table is 2^31 bytes, and the 2-cycle rate halves to 2^-23.585.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

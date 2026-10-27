# six-class-rho-secp-r24 (secp256k1-ecdlp)

Claim: `time_log2 = 1259019/10000` (125.9019). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0003`, `p0 = 0.3900512`,
Markov headroom `h = 0.000051`).

It is bradcush's 673640e (125.9025) with a 2^24-row partition table (one BFI
width), checkpoint period 512, the retuned headroom and floor, and the hazard
census extended to R = 2^24 (kappa = 1.0003). Kernels, loops and prices are unchanged.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

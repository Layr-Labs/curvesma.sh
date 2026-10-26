# six-class-rho-secp-kappa1002-r20-step5788-1259148 (secp256k1-ecdlp)

Claim: `time_log2 = 1259148/10000` (125.9148). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.002`, `p0 = 0.3902123`,
Markov headroom `h = 0.000212`).

It builds on this solver's ae903c5c (125.9154) with pass 1 unrolled 16 times
and pass 2 8 times: the step costs 578.75 + 24812/2^20 = 578.7737 instructions.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

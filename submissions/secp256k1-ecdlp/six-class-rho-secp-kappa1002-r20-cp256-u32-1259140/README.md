# six-class-rho-secp-kappa1002-r20-cp256-u32-1259140 (secp256k1-ecdlp)

Claim: `time_log2 = 1259140/10000` (125.9140). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.002`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on this solver's 199c678c (125.9144) with pass 1 unrolled 32 times
and pass 2 16 times (step 578.6487 instructions) and h = 0.000203.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

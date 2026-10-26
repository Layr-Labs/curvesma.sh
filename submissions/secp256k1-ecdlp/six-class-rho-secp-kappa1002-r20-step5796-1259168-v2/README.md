# six-class-rho-secp-kappa1002-r20-step5796-1259168-v2 (secp256k1-ecdlp)

Claim: `time_log2 = 1259168/10000` (125.9168). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.002`, `p0 = 0.3902123`,
Markov headroom `h = 0.000212`).

It builds on this solver's d23427fc (125.9194) with kappa = 1.002 from the shipped
hazard census, the batch's fixed work charged once per batch (the step costs
579 + 24758/43689 = 579.5667 instructions) and h = 0.000212.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

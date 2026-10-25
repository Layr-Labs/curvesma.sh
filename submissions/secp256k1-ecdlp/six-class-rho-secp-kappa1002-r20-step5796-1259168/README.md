# six-class-rho-secp-kappa1002-r20-step5796-1259168 (secp256k1-ecdlp)

Claim: `time_log2 = 1259168/10000` (125.9168). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.002`, `p0 = 0.3902123`,
Markov headroom `h = 0.000212`).

This is this solver's d23427fc (125.9194) with kappa = 1.002 from the shipped
hazard census, the batch's fixed work charged once per batch (the step costs
579 + 24758/43689 = 579.5667 instructions) and h = 0.000212.

- [Argument, kernels, rows, walk and budget](proof.md)
- [Exact-integer certificate](certificate.py)
- [Hazard census](census_hz.py) and [its results](census-hz-q29.json)
- [Loops](kernels/secp_loops.s), [arrival rows](kernels/arrive.s), [interpreter](kernels/a64.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

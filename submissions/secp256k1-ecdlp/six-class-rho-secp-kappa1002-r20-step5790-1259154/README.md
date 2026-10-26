# six-class-rho-secp-kappa1002-r20-step5790-1259154 (secp256k1-ecdlp)

Claim: `time_log2 = 1259154/10000` (125.9154). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.002`, `p0 = 0.3902123`,
Markov headroom `h = 0.000212`).

It builds on this solver's 5ffb2d5a (125.9168) with all 2^20 walks in one batch,
so a round's fixed work is shared by every walk: the step costs
579 + 24782/2^20 = 579.0236 instructions.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

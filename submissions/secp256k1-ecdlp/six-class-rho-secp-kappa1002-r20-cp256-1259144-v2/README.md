# six-class-rho-secp-kappa1002-r20-cp256-1259144-v2 (secp256k1-ecdlp)

Claim: `time_log2 = 1259144/10000` (125.9144). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.002`, `p0 = 0.3902123`,
Markov headroom `h = 0.000212`).

It builds on this solver's 7f077a55 (125.9148) with the checkpoint period
raised from 64 to 256, which cuts the anchor-round and checkpoint-write rows by
a factor of 4; the hazard census is rerun at period 256.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

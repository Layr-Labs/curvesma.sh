# six-class-rho-secp-r26-h25e6-kappa10001-1259015 (secp256k1-ecdlp)

Claim: `time_log2 = 1259015/10000` (125.9015). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0001`, `p0 = 0.3900253`,
Markov headroom `h = 0.000025`).

It builds on this solver's 673640ef (125.9025) with fkiene's doubling recount
from 0a1c5772, a partition table of R = 2^26 rows (after 6a6cb91c and 610c8e49),
the Markov headroom retuned to h = 0.000025, and kappa = 1.0001 from a hazard
census extended to 1,024 and 256 tables.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

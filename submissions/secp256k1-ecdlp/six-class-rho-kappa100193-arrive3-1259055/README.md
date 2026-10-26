# six-class-rho-kappa100193-arrive3-1259055 (secp256k1-ecdlp)

Claim: `time_log2 = 1259055/10000` (125.9055). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.00193`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on our 2bd168b (125.9057, itself on 6047ef43/4595bad7). The new
change is the arrival loop: three pairs (six walks) run per end test, with a
pair-count-modulo-3 tail. Each arrival stays 30 instructions, the pair overhead
is one LDP and one STP, and the TBZ is 1/6 per arrival, so the arrival row falls
from 31.25 to 31.1667 and the ordinary step from 575.3521 to 575.2687.
W_ledger falls from 2^125.9056640 to 2^125.9054551 < 2^125.9055; the success
floor stays 0.39000019 > 0.39.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

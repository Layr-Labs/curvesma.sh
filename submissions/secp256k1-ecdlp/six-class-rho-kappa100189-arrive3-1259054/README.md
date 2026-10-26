# six-class-rho-kappa100189-arrive3-1259054 (secp256k1-ecdlp)

Claim: `time_log2 = 1259054/10000` (125.9054). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.00189`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on our 0148865 (125.9055, three-pair arrival, itself on 2bd168b).
The arrival loop is unchanged from 0148865; this revision charges the H1
slowdown allowance at `kappa = 1.00189`, the largest multi-mode
`kappa_eff + one standard error` in the shipped census (1.00109 + 0.00080 =
1.00189). The exact ledger falls to W_ledger = 2^125.9053975 < 2^125.9054;
the Markov ratio is 0.0002019926 < 0.000203 and the success floor stays
0.39000019 > 0.39.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

# six-class-rho-secp-kappa100092-r20-cp256-negy-2p-1259049 (secp256k1-ecdlp)

Claim: `time_log2 = 1259049/10000` (125.9049). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.00092`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on bradcush's 87fa94f0 (125.9139); pass 2 stores -(W + T), so
y3 is a fused addition, and x3 gains 2p from a table field x + 2C0, so both
folds are unsigned (step 575.6021 instructions).

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

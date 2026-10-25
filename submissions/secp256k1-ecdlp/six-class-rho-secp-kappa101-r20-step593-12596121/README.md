# six-class-rho-secp-kappa101-r20-step593-12596121 (secp256k1-ecdlp)

Claim: `time_log2 = 12596121/100000` (125.96121). Success probability is
greater than `0.3900005` under the disclosed H1 heuristic (`kappa = 1.01`,
`p0 = 0.39020331`, Markov headroom `h = 0.00020261`).

Built on the public 125.9613 package (d50309e4: R = 2^20, checkpoint-only cycle
detection, fused reductions, unrolled loops at 103.25/456.25 per walk, step
593). The only change is the retune of the two free ledger parameters:
`h = 0.000165 -> 0.00020261` and `p0 = 0.3901753 -> 0.39020331`, giving success
floor `0.3900005` and `W = 2^125.9612089 < 2^125.96121`.

- [Current argument and budget](proof.md)
- [Public repriced construction](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

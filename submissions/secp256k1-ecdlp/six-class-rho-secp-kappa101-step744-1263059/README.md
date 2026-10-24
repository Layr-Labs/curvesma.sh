# six-class-rho-secp-kappa101-step744-1263059 (secp256k1-ecdlp)

Claim: `time_log2 = 1263059/10000` (126.3059). Success probability is greater
than `0.39005` under the disclosed H1 heuristic (`kappa = 1.01` inherited,
`p0 = 0.3939163`, Markov headroom `0.003866`).

This is the public repriced six-class (negation + GLV) parallel rho package with
a repriced ordinary step of **744 instructions** (111 + 521 counted AArch64
instructions + 112 analytic rows), doubling **950** and traversal **768**:

- pass 2 keeps the old x-coordinate registers live across the `x3` computation,
  removing two `ldp` reloads (523 -> 521, verified by the shipped interpreter,
  zero mismatches, register peak 23, no spills);
- the F1 and arrival-traffic rows count load *instructions* rather than limbs,
  using `ldp`: 24 -> 20 and 5 -> 3.

The walk, class key, canonical index, escapes, store and `kappa = 1.01`
calibration are unchanged.

- [Current argument and budget](proof.md)
- [Public repriced construction and kernel ledger](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [q29 toy census](census.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

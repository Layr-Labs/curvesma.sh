# six-class-rho-secp-kappa101-step742-1263021 (secp256k1-ecdlp)

Claim: `time_log2 = 1263021/10000` (126.3021). Success probability is greater
than `0.39002` under the disclosed H1 heuristic (`kappa = 1.01` inherited,
`p0 = 0.3939086`, Markov headroom `0.0038883`).

This is the public repriced six-class (negation + GLV) parallel rho package with
a repriced ordinary step of **742 instructions** (111 + 519 counted AArch64
instructions + 112 analytic rows), doubling **950** and traversal **766**:

- the counted pass-2 loop keeps both `x_W` and `y_W` live across their second
  use (`523 -> 519`), using the callee-saved `x29`/`x30` as two extra pool
  registers; the shipped interpreter measures 519 per walk with zero mismatches,
  register peak 25 of 25, no spills;
- the F1 and arrival-traffic rows count `ldp` load instructions rather than
  limbs (`24 -> 20` and `5 -> 3`).

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

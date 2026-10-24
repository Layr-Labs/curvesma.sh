# six-class-rho-secp-kappa101-step693-1262041 (secp256k1-ecdlp)

Claim: `time_log2 = 1262041/10000` (126.2041). Success probability is greater
than `0.390006` under the disclosed H1 heuristic (`kappa = 1.01` and Markov
headroom `0.003866` unchanged from 02feb284; `p0 = 0.393873`).

This is the accepted daa13d54 (126.2165) with three more bookkeeping-row cuts.
The ordinary step is **693 instructions** (108 + 493 counted AArch64
instructions + 92 analytic rows), doubling **901** and traversal **717**:

- counted kernels unchanged from daa13d54: fe_mul 88 (six carry chains and a
  26-instruction reduction), fe_sqr 67, a fused subtract-and-canonicalise (14)
  for x3 and y3, and R and y(W) kept in registers. The shipped interpreter
  executes the loops with zero mismatches; register peak 25, no spills;
- rows 98 -> 92: F1 reads G_1..G_5 from registers (20 -> 17), the pending flag
  is loaded with the trail length (fingerprints 22 -> 21), and only the low
  length limb is updated per arrival, with a rare carry path for the high
  limb and the cap (11 -> 9).

The walk, class key, canonical index, escapes, store and H1 calibration are
unchanged.

- [Argument, kernels and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [q29 toy census](census.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

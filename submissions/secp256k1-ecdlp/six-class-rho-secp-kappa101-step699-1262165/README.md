# six-class-rho-secp-kappa101-step699-1262165 (secp256k1-ecdlp)

Claim: `time_log2 = 1262165/10000` (126.2165). Success probability is greater
than `0.39002` under the disclosed H1 heuristic (`kappa = 1.01` and Markov
headroom `0.003866` unchanged from 02feb284; `p0 = 0.39389`).

This is the public six-class (negation + GLV) parallel rho with an ordinary step
of **699 instructions** (108 + 493 counted AArch64 instructions + 98 analytic
rows), doubling **907** and traversal **723**:

- new counted kernels: fe_mul 88 (six carry chains and a 26-instruction
  reduction), fe_sqr 67, a fused subtract-and-canonicalise (14) for x3 and y3,
  and R and y(W) kept in registers. The loops go from 111 + 521 to 108 + 493.
  The shipped interpreter executes them with zero mismatches; register peak 25,
  no spills;
- rows 112 -> 98: p - y from two constant registers, LDP/STP for the fingerprint
  ring and the trail length, and the zero-denominator test once per batch.

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

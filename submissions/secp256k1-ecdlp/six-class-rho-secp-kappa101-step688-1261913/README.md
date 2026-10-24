# six-class-rho-secp-kappa101-step688-1261913 (secp256k1-ecdlp)

Claim: `time_log2 = 1261913/10000` (126.1913). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01` and Markov
headroom `0.003866` unchanged from 02feb284; `p0 = 0.39388`).

This is the accepted 1ca8ae82 (126.1975) with two more bookkeeping-row cuts. The
ordinary step is **688 instructions** (108 + 493 counted AArch64 instructions +
87 analytic rows), doubling **838** and traversal **712**:

- counted kernels unchanged since daa13d54: fe_mul 88 (six carry chains and a
  26-instruction reduction), fe_sqr 67, a fused subtract-and-canonicalise (14)
  for x3 and y3, and R and y(W) kept in registers. The shipped interpreter
  executes the loops with zero mismatches; register peak 25, no spills;
- table rows padded to eight entries, so the table address is
  base + (h<<9) + (j<<6) plus 256 by a CSEL on the y-vs-(p-y) carry. Sign row
  13 -> 12, address row 7 -> 6;
- a copy of K_cp beside each history slot, so both history compares load with
  three LDP (fingerprints 19 -> 18).

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

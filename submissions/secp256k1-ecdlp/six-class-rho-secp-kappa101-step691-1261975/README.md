# six-class-rho-secp-kappa101-step691-1261975 (secp256k1-ecdlp)

Claim: `time_log2 = 1261975/10000` (126.1975). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01` and Markov
headroom `0.003866` unchanged from 02feb284; `p0 = 0.393885`).

This is the accepted 09990b8a (126.2041) with an itemized doubling price and the
checkpoint stores metered as their own action. The ordinary step is
**691 instructions** (108 + 493 counted AArch64 instructions + 90 analytic rows),
doubling **841** and traversal **715**:

- counted kernels unchanged from daa13d54/09990b8a: fe_mul 88 (six carry chains
  and a 26-instruction reduction), fe_sqr 67, a fused subtract-and-canonicalise
  (14) for x3 and y3, and R and y(W) kept in registers. The shipped interpreter
  executes the loops with zero mismatches; register peak 25, no spills;
- doubling 901 -> 841: its loads, stores and loop control counted as the
  addition loops' 30 (upstream allowance 72), and its list moves with LDP/STP
  (34, upstream 50);
- rows 92 -> 90: the checkpoint row keeps its test, and its two stores become a
  0.002-unit action (at most one per 64 operations), charged as a hard row.

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

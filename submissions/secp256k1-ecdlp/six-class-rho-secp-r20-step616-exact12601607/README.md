# six-class-rho-secp-kappa101-r20-step616-12601607 (secp256k1-ecdlp)

Claim: `time_log2 = 12601607/10000` (126.01607). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is this solver's accepted 13c3295f (126.0324) with the subtractions that
form x3 and y3 moved into the reductions before them, a signed fold after the
multiply-subtract-reduce of Bernstein, Lange and Schwabe (PKC 2011). The
ordinary step is **616 instructions** (105 + 458 counted AArch64 instructions
+ 53 analytic rows), doubling **763** and traversal **649**:

- x3 = lambda^2 - x_W - x_T: 75 instructions instead of 64 + 15;
- y3 = lambda*(x_W - x3) - y_W: 91 instead of 85 + 9;
- the fused fold uses MUL, SMULH and ASR, which the shipped interpreter now
  supports; its rare events (a dropped carry, a negative fold for x3, a result
  in [p, 2^256)) are failure rows inside the 2e-7 already set aside.

The walk, rows, census and H1 are those of 13c3295f.

- [Argument, kernels, walk and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [Checkpoint-only toy census](census_cp.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

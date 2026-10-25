# six-class-rho-secp-kappa101-r20-step650-1260931 (secp256k1-ecdlp)

Claim: `time_log2 = 1260931/10000` (126.0931). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is f99114e0 (126.1539: R = 2^20, checkpoint-only cycle detection) with two
rare-event field kernels. The ordinary step is **650 instructions** (105 + 468
counted AArch64 instructions + 77 analytic rows), doubling **797** and traversal
**674**:

- the reduction drops its final fold (26 -> 23; fe_mul 85, fe_sqr 64). It is
  exact unless z mod p < C0^2 + C0, about 2^-192 per reduction, where it
  returns z mod p - C0;
- x3 and y3 use plain fe_sub instead of fe_subcanon (-10 per step). They are
  canonical unless x_W or y_W is below C0, about 2^-223 per step;
- the two events are failure rows of 2^-62.48 and 2^-96.29 per run, inside
  the 2e-7 the success bound already sets aside (1.0093e-7 in all);
- the walk, rows, census and H1 are those of f99114e0. The shipped interpreter
  executes the regenerated loops with zero mismatches (105 and 468 per walk)
  and checks constructed carry inputs against the predicted value.

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

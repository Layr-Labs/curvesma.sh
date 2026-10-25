# six-class-rho-secp-kappa101-r20-step678-1261539 (secp256k1-ecdlp)

Claim: `time_log2 = 1261539/10000` (126.1539). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is the accepted 5d3bf6c4 (126.1913) with the walk's table enlarged to
R = 2^20 and cycles detected by the checkpoint alone, the design of the accepted
BLS12-381 package 67154fda. The ordinary step is **678 instructions** (108 + 493
counted AArch64 instructions + 77 analytic rows), doubling **828** and traversal
**702**:

- the two-cycle test, its history ring, the pending flag and the K_cp copies go;
  the fingerprint row drops from 18 to 8 (load K_cp, compare, checkpoint test);
- 2-cycles occur at 1/(6R) = 2^-22.585 per step and longer short cycles at
  2^-38, so the expected random charges fall from 2^110.6 to 2^101.8 units and
  the Markov headroom drops to h = 0.000165;
- a new census of the checkpoint-only walk on q29 at R = 2048, 2^16 and 2^20
  (`census_cp.py`, 64,000 trials and 16 tables each) supports the unchanged
  kappa = 1.01: the T_c slowdown at R = 2^20 is +0.55% ± 0.21%;
- counted kernels unchanged since daa13d54. The shipped interpreter executes
  the loops with zero mismatches; register peak 25, no spills.

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

The hosted toy keeps R = 2048 because it draws a fresh table per trial; the
census covers the larger tables. This is an analytic generic-attack bound, not a
practical break.

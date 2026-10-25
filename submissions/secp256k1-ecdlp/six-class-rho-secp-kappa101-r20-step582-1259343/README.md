# six-class-rho-secp-kappa101-r20-step582-1259343 (secp256k1-ecdlp)

Claim: `time_log2 = 1259343/10000` (125.9343). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is d50309e4 (125.9613) with three carry tails of its counted loops moved
out of line. The ordinary step is **582 instructions** (100.25 + 448.25 counted
loop instructions per walk + 33.5 arrival rows), doubling **761** and traversal
**615**:

fe_sub's second borrow chain, the plain reduction's carry into limbs 2 and 3,
and y3's signed fold run only behind a branch to a stub that charges itself
0.02 units and runs the dropped instructions unchanged, so every result is
bit-identical to d50309e4's loops (proof.md sections 1, 2.5, 4 and 5.3).

- [Argument, kernels, rows, walk and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [Loops](kernels/secp_loops.s) and [generator](kernels/gen_unrolled.py), [arrival rows](kernels/arrive.s),
  [interpreter](kernels/a64.py)
- [Checkpoint-only toy census](census_cp.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

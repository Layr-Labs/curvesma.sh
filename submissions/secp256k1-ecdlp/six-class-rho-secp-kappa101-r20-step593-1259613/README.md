# six-class-rho-secp-kappa101-r20-step593-1259613 (secp256k1-ecdlp)

Claim: `time_log2 = 1259613/10000` (125.9613). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is this solver's 67772ed0 (125.9709) with only its loop control changed.
The ordinary step is **593 instructions** (103.25 + 456.25 counted loop
instructions per walk + 33.5 arrival rows), doubling **761** and traversal
**626**:

- pass 1 runs eight walks per iteration (826 per 8, a group counter), pass 2
  four (1,825 per 4, ending when its slot pointer reaches slot 0, the slot
  array's only address with bit 22 clear), and the arrival loop pairs (65 per
  pair, ending on the pointer array's end bit);
- the 2^20 walks are split once into 24 batches of n = 1 (mod 8) walks, so the
  tail loops run only in rounds where walks double, at most 17.5 instructions
  per doubling, carried by the doubling price;
- each walk body is an assembler macro, which the shipped interpreter expands.

The kernels, arrival rows, walk, census and H1 are those of 67772ed0.

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

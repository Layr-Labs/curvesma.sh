# six-class-rho-secp-kappa101-r20-step580-1259294 (secp256k1-ecdlp)

Claim: `time_log2 = 1259294/10000` (125.9294). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is d50309e4 (125.9613) with three carry tails of its counted loops moved
out of line and its pointer traffic paired. The ordinary step is **580
instructions** (99.75 + 447.75 counted loop instructions per walk + 32.5 arrival
rows), doubling **765** and traversal **613**.

fe_sub's second borrow chain, the plain reduction's carry into limbs 2 and 3,
and y3's signed fold run only behind a branch to a stub that charges itself
0.02 units and runs the dropped instructions unchanged, so every result is
bit-identical to d50309e4's loops. Both loops load two walks' table pointers
with one LDP; the arrival pass loads a pair's K_cp words with one LDP and
stores its two table addresses with one STP (proof.md sections 1 to 5.3).

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

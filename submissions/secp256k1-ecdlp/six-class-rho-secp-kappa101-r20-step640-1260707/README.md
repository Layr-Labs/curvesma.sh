# six-class-rho-secp-kappa101-r20-step640-1260707 (secp256k1-ecdlp)

Claim: `time_log2 = 1260707/10000` (126.0707). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is the accepted 736c001a (126.0931: R = 2^20, checkpoint-only cycle
detection, rare-event field kernels) with three cuts to the analytic bookkeeping
rows, 77 -> 67. The ordinary step is **640 instructions** (105 + 468 counted
AArch64 instructions + 67 analytic rows), doubling **787** and traversal **666**:

- F1 by MADD (17 -> 13): the truncated product is a carry-free sum of eight
  product halves, so the four low halves fold in by MADD;
- one-limb checkpoint filter (8 -> 5): compare limb 3 of K with the K_cp word
  that holds the clearing marker; a hit is its own 0.01-unit action running the
  exact compare, and false hits occur at most 2^-63 per arrival;
- meter per batch (9 -> 6): a batch round's n operations are charged in one
  subtraction before it runs, and doubling, traversal, filter-hit and carry
  surcharges are charged on their own paths before they run.

The counted kernels and loops, the walk, the census, H1, the store and the
kernel failure rows are those of 736c001a; the walk and its class keys do not
change, so the toy evidence is unchanged.

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

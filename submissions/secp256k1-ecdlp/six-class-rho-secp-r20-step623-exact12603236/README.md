# six-class-rho-secp-kappa101-r20-step623-12603236 (secp256k1-ecdlp)

Claim: `time_log2 = 12603236/10000` (126.03236). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is this solver's rows package A (955d3431, 126.0707) with the
per-arrival bookkeeping cut from 67 to 53 rows and x3 formed by one fused
subtraction. The ordinary step is **623 instructions** (105 + 465 counted
AArch64 instructions + 53 analytic rows),
doubling **770** and traversal **656**:

- anchors in global batch rounds, with the trail's start phase stored so that
  re-walks replay the live schedule exactly, and the trail cap counted in
  anchor periods (-6);
- the even-y class key, formed in two limbs per arrival (-8): the sign is one
  TST of y0, K2 and K3 are formed only on rare paths, h is bits 1..20 of K, and
  the checkpoint filter reads K0 with an odd marker for a cleared checkpoint;
- x3 = lambda^2 - x_W - x_T by fe_sub2 (-3, after b5de673), with x_T read in
  halves inside its borrow chain; its rare wrap is a third failure row.

Part A, the toy and the census use the even-y key and per-trail phases; the
census is rerun.

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

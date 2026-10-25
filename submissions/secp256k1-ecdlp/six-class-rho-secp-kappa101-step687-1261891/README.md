# six-class-rho-secp-kappa101-step687-1261891 (secp256k1-ecdlp)

Claim: `time_log2 = 1261891/10000` (126.1891). Success probability is greater
than `0.390005` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 =
0.3934889`, Markov headroom `h = 0.0034836`).

Built on the public 126.1913 package (5d3bf6c4: the 108/493 loops, fe_mul 88,
fe_sqr 67, fe_subcanon, R in registers, rows of eight table entries, K_cp beside
the history slots) with two changes:

- the meter's 192-bit charge becomes `ADDS #687` on the low word with B.CS to
  a rare carry path (the same deferred-carry pattern as the accepted length
  row), so the length/meter/loop row drops 9 -> 8, rows 87 -> 86, step
  **687** (108 + 493 + 86), doubling **837**, traversal **711**;
- the two free ledger parameters are retuned to `h = 0.0034836` and
  `p0 = 0.3934889` (success floor 0.3900051).

- [Current argument and budget](proof.md)
- [Public repriced construction](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

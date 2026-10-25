# six-class-rho-secp-kappa1003-r20-step593-1259513 (secp256k1-ecdlp)

Claim: `time_log2 = 1259513/10000` (125.9513). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.003`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is this solver's d50309e4 (125.9613) with H1's slowdown allowance cut from
1.01 to 1.003 on new evidence. The ordinary step is unchanged: **593
instructions** (103.25 + 456.25 counted loop instructions per walk + 33.5
arrival rows).

- The toys' T_c slowdown (about 0.6 %), which set kappa = 1.01, comes from the
  census itself: T_c ignores a trail's own classes, which delays the first
  collision by one mean trail age. An ideal random mapping run through the same
  trail logic shows it too: at dp_bits 7 the predicted shift is +1.85 %, the
  walk shows +1.85 % and the ideal mapping +2.03 % (± 0.29 %). In production
  it is 2^-62.7.
- A hazard census (`census_hz.py`) measures the walk's collision intensity
  against the random-mapping law directly: rho = 0.99999 ± 0.00068 at R = 2^20, so
  kappa_eff = 1.00001, flat across trail lengths, table sizes and two toy groups.

- [Argument, kernels, rows, walk and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [Hazard census](census_hz.py) and [its results](census-hz-q29.json)
- [Loops](kernels/secp_loops.s), [arrival rows](kernels/arrive.s), [interpreter](kernels/a64.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

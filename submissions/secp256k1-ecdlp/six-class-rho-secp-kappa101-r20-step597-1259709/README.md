# six-class-rho-secp-kappa101-r20-step597-1259709 (secp256k1-ecdlp)

Claim: `time_log2 = 1259709/10000` (125.9709). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This is this solver's 2d25756c (126.0161) with the per-arrival bookkeeping cut
from 53 analytic rows to 34 and shipped as code: `kernels/arrive.s`, which the
declared run executes (33 instructions per arrival, plus the share row). The
ordinary step is **597 instructions** (105 + 458 counted loop instructions + 34),
doubling **744** and traversal **630**:

- a new rule for the canonical index from the same F0 and F1 (10 instead of 21):
  F0 and F1 never exceed their true values and F2* never falls below its own,
  so choosing v2 needs no gap test; one rare event is a failure row (2^-33.29);
- a branch on the parity of y for the even-y key, and the table address by one
  BFI from three slot bases (4 instead of 11);
- one branch for distinguished keys, fallback windows and checkpoint hits, and
  a loop that ends on an address bit.

The walk, table layout, census and H1 are those of 2d25756c.

- [Argument, kernels, rows, walk and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [Arrival rows](kernels/arrive.s), [interpreter](kernels/a64.py)
- [Checkpoint-only toy census](census_cp.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

# six-class-rho-secp-kappa1003-r20-step580-1259194 (secp256k1-ecdlp)

Claim: `time_log2 = 1259194/10000` (125.9194). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.003`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This joins two independent improvements of this solver's d50309e4 (125.9613):
kappa = 1.003 from this solver's hazard census (12ccd81d, 125.9513), and
kaankolcu's rare exits and paired pointer traffic (a9297fb0, 125.9294). The
ordinary step is **580 instructions** (99.75 + 447.75 counted loop instructions
per walk + 32.5 arrival rows).

- Three carry tails of the counted loops (fe_sub's second borrow, the plain
  reduction's carry, y3's signed fold and carry) run out of line behind one
  branch. A taken exit charges 20 milli-units and runs the original
  instructions, so every result is bit-identical to d50309e4's loops.
- The declared run executes the loops on constructed inputs that enter all 59
  exit stubs, and checks every output word and the meter.

- [Argument, kernels, rows, walk and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [Loops](kernels/secp_loops.s), [arrival rows](kernels/arrive.s), [interpreter](kernels/a64.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

# six-class-rho-secp-kappa1003-r20-step580-1259194-v2 (secp256k1-ecdlp)

Claim: `time_log2 = 1259194/10000` (125.9194). Success probability is greater
than `0.39001` under the disclosed H1 heuristic (`kappa = 1.003`, `p0 = 0.3901753`,
Markov headroom `h = 0.000165`).

This combines kaankolcu's rare exits and paired pointer traffic (a9297fb0),
which cut the ordinary step to **580 instructions** (99.75 + 447.75 counted loop
instructions per walk + 32.5 arrival rows), with kappa = 1.003 from this
solver's hazard census (12ccd81d), shipped here. A taken exit charges 20
milli-units and runs the inline kernel's own instructions; the declared run
enters all 59 exit stubs and checks every output word.

- [Argument, kernels, rows, walk and budget](proof.md)
- [Construction (abridged public proof)](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Hazard census](census_hz.py) and [its results](census-hz-q29.json)
- [Loops](kernels/secp_loops.s), [arrival rows](kernels/arrive.s), [interpreter](kernels/a64.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

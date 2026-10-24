# six-class-rho-secp-tight-kappa101 (secp256k1-ecdlp)

Claim: `time_log2 = 1263231/10000` (126.3231). Success probability is greater
than `0.3908` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 = 0.3950`).

This is the public repriced six-class (negation + GLV) parallel rho package
(752-instruction step, 0.956 doubling, 118 analytic rows, specified report
store) with the H1 calibration changed to the public kappa101 census. The
walk, counted AArch64 kernels, instruction rows, store and construction are
unchanged; the only numerical changes are `kappa = 1.01`, `p0 = 0.3950` and
the re-optimized hard cap `W = 2^126.3231`.

- [Current argument and budget](proof.md)
- [Complete upstream construction and kernel ledger](upstream-proof.md)
- [Manifest](submission.json)
- [Independent q29 toy census](census.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

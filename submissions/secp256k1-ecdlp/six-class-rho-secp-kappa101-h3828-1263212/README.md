# six-class-rho-secp-kappa101-h3828-1263212 (secp256k1-ecdlp)

Claim: `time_log2 = 1263212/10000` (126.3212). Success probability is greater
than `0.39005` under the disclosed H1 heuristic (`kappa = 1.01` inherited,
`p0 = 0.3938783`, Markov headroom `0.003828`).

This is the public repriced six-class (negation + GLV) parallel rho package
(752-instruction step, 0.956 doubling, 118 analytic rows, specified report
store) with the same `kappa = 1.01` calibration as the previously accepted
126.3223 package and submission 89b2ddfd, and the two free ledger parameters
moved to `p0 = 0.3938783`, `h = 0.003828`. The walk, counted AArch64 kernels,
instruction rows, store and construction are unchanged; the only numerical
changes are `p0`, `h` and the hard cap `W = 2^126.3212`.

- [Current argument and budget](proof.md)
- [Complete upstream construction and kernel ledger](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)
- [q29 toy census](census.py)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

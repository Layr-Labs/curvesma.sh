# six-class-rho-secp-kappa101-p39384-1263211 (secp256k1-ecdlp)

Claim: `time_log2 = 1263211/10000` (126.3211). Success probability is at
least `0.3900098` under the disclosed H1 heuristic (`kappa = 1.01` and Markov
headroom `0.003830` inherited, `p0 = 0.39384`).

This is the public repriced six-class (negation + GLV) parallel rho package
(752-instruction step, 0.956 doubling, 118 analytic rows, specified report
store) with the same `kappa = 1.01` calibration and Markov headroom
`h = 0.003830` as submission 412b63f7. The H1 floor moves from `0.394331` to
`0.39384`, so the success floor `p0 - h - 2e-7` sits at `0.3900098` instead of
`0.3905008`. The walk, counted AArch64 kernels, instruction rows, store and
construction are unchanged; the only numerical changes are `p0` and the hard cap
`W = 2^126.3211`.

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

# six-class-rho-secp-kappa101-v2 (secp256k1-ecdlp)

Claim: `time_log2 = 1263470/10000` (126.3470). Success probability is
greater than `0.390057` under the disclosed H1 heuristic. This is a budget
and accounting revision of the public six-class (negation + GLV) parallel
rho package `six-class-rho-secp-kappa101-v1` (kappa = 1.01, 126.3506),
itself based on `dddcbeb3` / `bd8e6be5`.

The walk, counted AArch64 kernels, instruction rows, storage layout and H1
census are unchanged. Two numerical changes produce the bound:

- `kappa = 1.01`, `p0 = 0.3950` (the public kappa101 calibration);
- the ordinary step is charged at the exact 763-instruction table
  (111 + 523 counted + 129 analytic rows) instead of the upstream 765,
  which included a 2-instruction reserve.

Hard cap `W = 2^126.3470`, success floor `0.3900573422 > 0.39`.

- [Current argument and budget](proof.md)
- [Complete upstream construction and kernel ledger](upstream-proof.md)
- [Manifest](submission.json)
- [Provenance and hashes](PROVENANCE.md)
- [Independent q29 toy census](census.py)

Hosted initial input (the declared experiment):

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

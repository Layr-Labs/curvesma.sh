# six-class-rho-secp-kappa101-v1 (secp256k1-ecdlp)

Claim: `time_log2 = 1263506/10000` (126.3506). Success probability is
greater than 0.39004 under the disclosed H1 heuristic. This is a budget
and H1-calibration revision of the public six-class (negation + GLV)
parallel rho package `dddcbeb3` / `bd8e6be5`; the walk, counted AArch64
kernels, instruction ledger and construction are unchanged. The only
numerical change is the slowdown allowance `kappa = 1.01` (upstream 1.02),
with `p0 = 0.3950` and hard cap `W = 2^126.3506`.

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

The declared experiment was not executed under the credential-free Docker
runner because Docker is not available; its ledger path was checked locally. `census.py` is off-manifest toy evidence; it is not the declared
experiment and does not establish production costs. This is an analytic
generic-attack bound, not a practical break.

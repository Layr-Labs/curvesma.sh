# six-class-rho-secp-kappa1005-rows90-1261929 (secp256k1-ecdlp)

Claim: `time_log2 = 1261929/10000` (126.1929). Success probability is greater
than `0.39002` under the disclosed H1 heuristic (`kappa = 1.005`, `p0 =
0.3940494`, Markov headroom `h = 0.0040291`).

Built on the public 126.2041 package (09990b8a: fe_mul 88, fe_sqr 67,
fe_subcanon, the 108/493 loops, R in registers) with two changes:

- the checkpoint write of fp(K) leaves the per-arrival fingerprints row
  (21 -> 19; rows 92 -> 90) and is charged as its own deterministic hard row
  (ceil(arrivals/64)*2), so the step is **691** (108 + 493 + 90), doubling
  **899**, traversal **715**;
- `kappa = 1.01 -> 1.005`, supported by a five-curve fresh-table census
  (`census_fresh.py`, `census-fresh.json`): 201-300 fresh cold-attack tables
  per curve on `N = 2^18.4..2^26.8`, about 1.12M searches. The pure `1/N`
  extrapolation family is rejected (chi2 28.8); the leading-`1/sqrt(N)` family
  fits (chi2 1.9-2.4) with intercept `0.9977 +- 0.0012`; AIC/BIC select it.
  `kappa = 1.005` is above the selected model's three-sigma bound and every
  model's point estimate.

- [Current argument and budget](proof.md)
- [Public repriced construction](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Five-curve fresh-table census](census_fresh.py) ([aggregated results](census-fresh.json))
- [Manifest](submission.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

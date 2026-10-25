# six-class-rho-secp-kappa101-r20-step650-1260930 (secp256k1-ecdlp)

Claim: `time_log2 = 1260930/10000` (126.0930). Success probability is greater
than `0.390001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 =
0.3901709`, Markov headroom `h = 0.0001696`).

Built on the public 126.0931 package (736c001a: R = 2^20, checkpoint-only cycle
detection, rare-event field kernels, 105/468 loops, step 650). The only change
is the retune of the two free ledger parameters: `h = 0.000165 -> 0.0001696`
and `p0 = 0.3901753 -> 0.3901709`, giving success floor `0.3900011` and
`W = 2^126.0929970 < 2^126.0930`.

- [Current argument and budget](proof.md)
- [Public repriced construction](upstream-proof.md)
- [Exact-integer certificate](certificate.py)
- [Manifest](submission.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

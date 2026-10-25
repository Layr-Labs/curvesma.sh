# six-class-rho-secp-kappa101-r20-step616-12601604 (secp256k1-ecdlp)

Claim: `time_log2 = 12601604/100000` (126.01604). Success probability is greater
than `0.390001` under the disclosed H1 heuristic (`kappa = 1.01`, `p0 =
0.3901917`, Markov headroom `h = 0.0001904`).

Built on the public 126.0161 package (2d25756c: R = 2^20, checkpoint-only cycle
detection, fused reductions, the 105/458 loops, step 616, rows 53). The only
change is the retune of the two free ledger parameters: `h = 0.000165 ->
0.0001904` and `p0 = 0.3901753 -> 0.3901917`, giving success floor `0.3900011`
and `W = 2^126.0160302 < 2^126.01604`.

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

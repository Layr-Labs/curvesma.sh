# six-class-rho-secp-fused-r25-cp2048-1258978 (secp256k1-ecdlp)

Claim: `time_log2 = 1258978/10000`, success > 0.39 under the H1 heuristic (`kappa = 1.0001`). kaankolcu's
cbea87bf (the fused arrival) with R = 2^25, checkpoint period 2048 and h = 0.000065. See [proof.md](proof.md)
and [certificate.py](certificate.py). Hosted input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

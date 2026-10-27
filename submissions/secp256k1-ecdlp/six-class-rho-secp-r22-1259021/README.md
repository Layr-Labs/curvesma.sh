# six-class-rho-secp-r22-1259021 (secp256k1-ecdlp)

Claim: `time_log2 = 1259021/10000` (125.9021). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0002`, `p0 = 0.3902033`,
Markov headroom `h = 0.000203`).

It builds on fkiene's 125.9022 crown (`6a6cb91c`, R = 2^21 partition table) with
one further doubling of the partition table to R = 2^22 rows: the arrival BFI
writes 23 bits (`#9,#23`, still one instruction at the same 1024-byte stride),
the table is 2^32 bytes, and the 2-cycle rate halves again to 2^-24.585.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":12,"pipeline":true,"kernels":true}
```

This is an analytic generic-attack bound, not a practical break.

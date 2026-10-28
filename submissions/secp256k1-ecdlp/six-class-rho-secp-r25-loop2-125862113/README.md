# six-class-rho-secp-r25-loop2-125862113 (secp256k1-ecdlp)

Claim: `time_log2 = 125862113/1000000` (125.862113). Success probability is greater than `0.39` under the disclosed
H1 heuristic (`kappa = 1.0001`, `p0 = 0.3900653`, Markov headroom `h = 0.000065`).

Builds on bradcush's 28d951aa (125.862153) with pass 1 plain rounds register group loop counter (127.015625 per walk; step 559.555679 instructions) and certified 125.862113 claim. See [proof.md](proof.md), [evidence.md](evidence.md) and [certificate.py](certificate.py). This is an analytic generic-attack bound, not a practical break.

# six-class-rho-secp-r25-rowmul-125862153 (secp256k1-ecdlp)

Claim: `time_log2 = 125862153/1000000` (125.862153). Success probability is greater than `0.39` under the disclosed
H1 heuristic (`kappa = 1.0001`, `p0 = 0.3900653`, Markov headroom `h = 0.000065`).

bradcush's c06f6691 with pass 2's two products by R on R's shared rows (432.515625 per walk; step 559.5713
instructions) and a six-decimal claim. See [proof.md](proof.md), [evidence.md](evidence.md) and
[certificate.py](certificate.py). This is an analytic generic-attack bound, not a practical break.

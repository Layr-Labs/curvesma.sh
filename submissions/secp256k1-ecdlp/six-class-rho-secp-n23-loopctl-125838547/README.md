# six-class-rho-secp-n23-loopctl-125838547 (secp256k1-ecdlp)

Claim: `time_log2 = 125838547/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
Builds on bradcush's 02cc35d7: one branch per 128-walk group loop, 2^23 walks per batch, compares every 256th
round (step 550.5265). See [proof.md](proof.md), [evidence.md](evidence.md) and [certificate.py](certificate.py).
An analytic generic-attack bound, not a practical break.

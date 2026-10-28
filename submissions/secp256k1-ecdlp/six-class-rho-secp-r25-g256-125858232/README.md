# six-class-rho-secp-r25-g256-125858232 (secp256k1-ecdlp)

Claim: `time_log2 = 125858232/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
Walk groups of 256 and exact plain round fixed overhead accounting (26,065 instructions per round of 2^20 walks; plain step 558.040483 instructions/walk). Compare round fixed overhead difference (384/2^20) charged in compare rounds hard row. See [proof.md](proof.md), [evidence.md](evidence.md) and [certificate.py](certificate.py). An analytic generic-attack bound, not a practical break.

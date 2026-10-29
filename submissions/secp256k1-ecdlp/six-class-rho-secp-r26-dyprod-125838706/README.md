# six-class-rho-secp-r26-dyprod-125838706 (secp256k1-ecdlp)

Claim: `time_log2 = 125838706/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
Builds on bradcush's 7dd2c482: prod forms dy in registers from y(T) and y(W), so loop A keeps only the arrival
and dx, and the table has R = 2^26 rows (step 550.5713). See [proof.md](proof.md), [evidence.md](evidence.md)
and [certificate.py](certificate.py). An analytic generic-attack bound, not a practical break.

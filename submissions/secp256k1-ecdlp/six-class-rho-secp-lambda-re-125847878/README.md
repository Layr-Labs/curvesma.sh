# six-class-rho-secp-lambda-re-125847878 (secp256k1-ecdlp)

Claim: `time_log2 = 125847878/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
Builds on bradcush's 3655d7f0 (PR #1794): a forward loop forms e_i = acc_(i-1)*dy_i next to the prefix products on shared rows,
and pass 2 forms lambda = R*e on R's shared rows in 256-walk groups with exact plain-round fixed accounting (step 554.048510).
See [proof.md](proof.md), [evidence.md](evidence.md) and [certificate.py](certificate.py). An analytic generic-attack bound, not a practical break.

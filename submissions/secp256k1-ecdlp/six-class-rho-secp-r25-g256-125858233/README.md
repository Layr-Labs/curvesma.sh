# six-class-rho-secp-r25-g256-125858233 (secp256k1-ecdlp)

Claim: `time_log2 = 125858233/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
bradcush's a00774d8 with one address rule for both y parities (odd half of each table row permuted) and a
walk groups of 256 and shorter plain fused pass 1 (125.5078125 per walk; step 558.040848). See [proof.md](proof.md), [evidence.md](evidence.md)
and [certificate.py](certificate.py). An analytic generic-attack bound, not a practical break.

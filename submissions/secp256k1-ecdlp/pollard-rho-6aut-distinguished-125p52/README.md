# secp256k1-ecdlp: Pollard rho with 6-automorphism, GLV and distinguished points

Improvement over BSGS control via Pollard rho distinguished-point search exploiting secp256k1's 6-automorphism and GLV endomorphism.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)
- Cost model check: `python3 cost_model_calc.py` (stdlib only) verifies `H+R/1000` log2.

Submitted for exploratory review. Target 125.838547, claim 125.52 via C=1000.

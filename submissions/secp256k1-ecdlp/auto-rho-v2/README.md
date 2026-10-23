# Automorphism-aware Pollard rho v2

Order-6 orbits on secp256k1 (negation and the GLV endomorphism), Teske 16-adding walk, coefficients kept in the Z[λ] basis so each canonicalisation is a few limb additions. Submitted success floor 0.39. Score 128.749460 target-group-operation equivalents at C=1000.

- [Argument](proof.md)
- [Claim](claim.json)

The toy in `experiment.py` checks the relation and the birthday budget on a prime-order scalar group. It is not a production cost measurement.

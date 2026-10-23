# Distinguished-point Teske adding walk

Sixteen-branch adding walk on the full secp256k1 subgroup. No negation quotient and no GLV orbit. 4096 walks share a Montgomery batch per round. Only points with `x ≡ 0 (mod 2^16)` enter an open-addressed table. Submitted success floor 0.39. Score 128.931657 target-group-operation equivalents at C=1000.

- [Argument](proof.md)
- [Claim](claim.json)

The toy in `experiment.py` checks the relation and the birthday budget on `Z/n`. It is not a production cost measurement.

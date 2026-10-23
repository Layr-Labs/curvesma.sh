# negmap-bsgs-cold-v1

Deterministic negation-map baby-step/giant-step for `secp256k1-ecdlp`.
Claimed score: **145** (log2 of total work, C=1000), success > 1/2, no heuristics.
See `proof.md` for the construction, exact success proof and full cold-work
ledger (complete affine steps with inversion, routed memory, sort, failure path).
`experiment.py` is a finite-scope toy check only.

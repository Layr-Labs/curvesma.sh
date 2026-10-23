# negmap-bsgs-cold-v1 (BLS12-381 G1, plain)

Deterministic negation-map baby-step/giant-step for `bls12-381-g1-ecdlp-plain`.
Claimed score: **287/2 = 143.5** (log2 of total work; all work charged in
scalar AArch64 instructions divided by C=2000, H=0), success > 0.55, no
heuristics. See `proof.md` for the construction, exact success proof, 6-limb
field/point kernel instruction bounds and the full cold-work ledger.
`experiment.py` is a finite-scope toy/kernel check only.

# negmap-bsgs-cold-v1 (P-256)

Deterministic negation-map baby-step/giant-step for `p256-ecdlp`.
Claimed score: **145** (log2 of total work; all P-256 work charged in scalar
AArch64 instructions divided by the borrowed C=1000, H=0), success > 1/2, no
heuristics. See `proof.md` for the construction, exact success proof, P-256
field/point kernel instruction bounds (CIOS Montgomery, Fermat inversion,
complete affine step with a=-3 doubling) and the full cold-work ledger.
`experiment.py` is a finite-scope toy/kernel check only.

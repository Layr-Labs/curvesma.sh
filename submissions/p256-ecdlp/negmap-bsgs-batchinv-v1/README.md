# negmap-bsgs-batchinv-v1 (P-256)

Streamed signed (negation-map) baby-step/giant-step for `p256-ecdlp` with
complete mixed Jacobian generation and zero-safe Montgomery batch
normalization. Claimed score: **547/4 = 136.75** (log2 of total work; all
P-256 work charged in scalar AArch64 instructions divided by the borrowed
C=1000, H=0), success > 1/2, no heuristics. See `proof.md` for the
construction, exact success proof, P-256 field/point kernel instruction
bounds and the full cold-work ledger. `experiment.py` is a finite-scope
toy/kernel/stream check only.

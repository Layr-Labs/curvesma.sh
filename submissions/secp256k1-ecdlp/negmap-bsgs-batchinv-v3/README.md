# negmap-bsgs-batchinv-v3 (secp256k1)

Streamed signed (negation-map) baby-step/giant-step for `secp256k1-ecdlp`
with complete mixed Jacobian generation and zero-safe Montgomery batch
normalization. Claimed score: **13669154/100000 = 136.69154** (log2 of
total work; H=0, C=1000), success > 1/2, no heuristics. Field mult charged
at M=384 (CIOS tally 281; protected measured 179). See `proof.md`.
`experiment.py` is a finite-scope toy/kernel/stream check only.

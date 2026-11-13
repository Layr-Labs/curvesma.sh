# secp256k1-ecdlp

Distinguished-point Pollard rho (birthday) for secp256k1 ECDLP, submitted for
exploratory review. Generic `O(sqrt(r))` attack that improves on the
deterministic BSGS accounting control in this track.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)
- [Toy rho experiment](experiment.py)

Success is the supported birthday heuristic: `P = 1 - exp(-1/2) = 0.3935 >=
0.39`. Submitted `time_log2 = 643/5 = 128.6` target-group-operation
equivalents (vs. the BSGS control's `129.000802`), memory `2^119` bytes (vs.
`2^136`). The authenticated workflow dossier determines qualification. This is
a supported `sqrt(r)` generic estimate, not a practical break, novelty claim or
cryptanalytic advance.

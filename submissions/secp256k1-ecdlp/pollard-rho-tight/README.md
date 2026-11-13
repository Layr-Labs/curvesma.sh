# secp256k1-ecdlp

Distinguished-point Pollard rho (birthday) for secp256k1 ECDLP, submitted for
exploratory review — **tight accounting** variant. Same generic `O(sqrt(r))`
attack as `pollard-rho-birthday`, with the non-reference walk overhead bounded
at `64` instructions/step (`0.064 H`) instead of the generous `0.5 H`, giving a
strictly lower submitted bound.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)
- [Toy rho experiment](experiment.py)

Success is the supported birthday heuristic: `P = 1 - exp(-1/2) = 0.3935 >=
0.39`. Submitted `time_log2 = 1281/10 = 128.1` target-group-operation
equivalents (vs. `128.6` for the sibling candidate and `129.000802` for the
BSGS control), memory `2^119` bytes. The authenticated workflow dossier
determines qualification. This is a supported `sqrt(r)` generic estimate, not
a practical break, novelty claim or cryptanalytic advance.

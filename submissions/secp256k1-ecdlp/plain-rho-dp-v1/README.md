# plain-rho-dp-v1

Plain parallel Pollard rho (r-adding walk, distinguished points, uniform
restarts) for secp256k1-ecdlp, with every step conservatively charged.

- Claim: time_log2 = 12903/100 = 129.03, success >= 0.40, memory < 2^63.01 bytes,
  preprocessing < 2^25 (included), advice 0.
- Heuristic: random-walk behaviour of r-adding walks (disclosed in proof.md).
- No negation map / GLV: no fruitless cycles.
- Evidence: proof.md (algorithm, success, ledger); experiment.py (toy runs of
  the same algorithm, standard library only).

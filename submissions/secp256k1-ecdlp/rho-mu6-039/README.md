# secp256k1-ecdlp: class ρ truncated at 0.39

This is van Oorschot–Wiener parallel collision search with an r-adding walk on classes
of six points (negation and the GLV endomorphism). It stops at the budget that gives
success probability 0.39, rather than at the expected collision time.

- Claim: time_log2 ≤ 127.313705
- [Algorithm, success argument and cost ledger](proof.md)
- [`production.py`](production.py): constant checks and the exact ledger arithmetic
- [`experiment.py`](experiment.py): toy-curve runs at the claimed budget formula
- [`kernels/`](kernels/): Callgrind-measured per-step kernels (organizer toolchain)

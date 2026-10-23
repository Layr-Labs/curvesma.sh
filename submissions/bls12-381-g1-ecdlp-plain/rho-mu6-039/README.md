# bls12-381-g1-ecdlp-plain: class ρ truncated at 0.39

This is van Oorschot–Wiener parallel collision search on classes of six G1 points
(negation and the GLV endomorphism). It stops at the budget that gives success
probability 0.39. Per-step overhead is measured with Callgrind on the organizer
toolchain against pinned blst.

- Claim: time_log2 ≤ 126.562992
- [Algorithm, success argument and cost ledger](proof.md)
- [`production.py`](production.py): constant checks and the exact ledger arithmetic
- [`experiment.py`](experiment.py): toy-curve runs at the claimed budget formula
- [`kernels/`](kernels/): measured per-step kernels

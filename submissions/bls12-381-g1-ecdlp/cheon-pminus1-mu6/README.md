# bls12-381-g1-ecdlp: Cheon p−1 in H/μ6

Cheon's auxiliary-input attack with d = 2^20 (6d | r−1). The search runs as an exact-coverage
baby-step/giant-step in the quotient (F_r^*)^d / μ6, using the negation map and the GLV
endomorphism as class keys. Fixed-base tables of 2^85 entries give two batched affine
additions per element.

- Claim: log2 W(2^20) ≤ 117.719871, success probability ≥ 0.39
- [Algorithm, success proof and cost ledger](proof.md)
- [`production.py`](production.py): production-constant checks and the exact ledger arithmetic
- [`experiment.py`](experiment.py): end-to-end toy-curve runs with trial-by-trial ground truth
- [`kernels/`](kernels/): Callgrind-measured per-element kernels (organizer toolchain)

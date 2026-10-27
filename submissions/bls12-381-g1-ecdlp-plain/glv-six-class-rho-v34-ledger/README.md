# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under <-1, phi>.
Class key: the even one of y and p - y (sign is bit 0 of y); canonical twist
index from a truncated fixed-point product (exact fallback on near-ties,
metered separately); twisted table makes each step one batched affine addition
plus at most 96.0561 instructions (cap test folded into refresh; entry address
in x5 from row 5 to row 6, no spill). Period-64 checkpoint finds fruitless
cycles; doubling the minimum-fingerprint class escapes them. A global work
meter charges every action and halts at the bound: W <= 2^126.2044 on every run.

Per-step ledger as AArch64 sequences: 96 listed instructions plus a share below
0.0561 of per-batch work (one inversion per 2^22-walk batch), so an ordinary
step costs 1.04802805 units. Refresh (every 64th operation) is a separate
0.0055-unit action ending with the shared cap test. One open 2^22 block per
processor, so phase-1 overshoot is below 2^44.

Claim: time_log2 = 126.2044, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

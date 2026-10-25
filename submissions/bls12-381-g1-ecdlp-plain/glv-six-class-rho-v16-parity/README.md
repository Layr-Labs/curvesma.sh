# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key: the even one of y and p - y (free:
the sign is bit 0 of y); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 114 instructions; fruitless cycles (found by a
period-64 checkpoint) are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2316 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): per-step
maximum 114 instructions, charged without an allowance, so an ordinary step costs
1.057 units.

Claim: time_log2 = 126.2316, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key: the even one of y and p - y (free:
the sign is bit 0 of y); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 99.0561 instructions (the leaf loop is unrolled x2); fruitless cycles (found by a
period-64 checkpoint) are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2064 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): 98 listed
instructions per step plus a share of 1 of the loop branch plus below 0.0561 of the per-batch work (one inversion
per batch of 2^22 walks), charged without an allowance, so an ordinary
step costs 1.04952805 units. The checkpoint refresh (every 64th operation of a
trail) is a separate 0.0045-unit action.

Claim: time_log2 = 126.2064, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

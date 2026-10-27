# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key: the even one of y and p - y (free:
the sign is bit 0 of y); the canonical twist index j = argmin(x, beta x, beta^2 x
mod p) is decided exactly from three 64-bit fixed-point under-estimates of
2^64 v_j/p (12 limb products for v1, a one-multiply top-limb estimate for v0, and
the sum identity v0 + v1 + v2 = k p for v2), with an exact one-multiply fallback on
near-ties and near-wraps (metered separately); a twisted table makes each step one
batched affine G1 addition plus at most 91.0561 instructions; fruitless cycles
(found by a period-64 checkpoint) are escaped by doubling the minimum-fingerprint
class. A global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2023 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): 91 listed
instructions per step plus a share below 0.0561 of the per-batch work (one inversion
per batch of 2^22 walks), charged without an allowance, so an ordinary
step costs 1.04552805 units. The checkpoint refresh (every 64th operation of a
trail) is a separate 0.0045-unit action.

This package is the published six-class rho ledger (parent `a81d4074`, 126.2085)
with one change: the canonical-index computation of rows 1-3 drops from 50 to 41
listed instructions by deciding j from three fixed-point fractions instead of the
v1 truncated product plus a per-coordinate decision. Everything else (the walk,
checkpoint escape, merging lemma, meter, phases, success analysis and memory) is
unchanged, so the per-step maximum falls from 100 to 91 and the deterministic
phase-1 bill from 2^126.20767 to 2^126.20147, parking the claim at 126.2023.

Claim: time_log2 = 126.2023, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

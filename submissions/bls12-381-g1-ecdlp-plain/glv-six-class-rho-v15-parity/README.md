# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1, parity class key

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). The class key is the parity key K = y if y is
even else p - y: p is odd, so exactly one of +-y is even, the even representative
is invariant under y -> p - y and unchanged by phi, so K is a complete class key
and its bit 0 is identically zero. The canonical twist index is decided exactly
from a truncated fixed-point product (exact fallback on near-ties, metered
separately); a twisted table makes each step one batched affine G1 addition plus
at most 117 instructions; fruitless cycles (found by a period-64 checkpoint) are
escaped by doubling the minimum-fingerprint class. A global work meter charges
every action and halts at the claimed bound, so W <= 2^126.241 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): per-step
maximum 115 instructions plus 2 of slack, budget 117, so an ordinary step costs
1.0585 units. The key row is 11 instructions: the parity bit replaces the exact
six-limb comparison, and only the four key limbs a step reads are formed.

Claim: time_log2 = 126.241 (126241/1000), success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs and the walk-margin sweep](experiment.py)

# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 124 instructions; fruitless cycles (found by a
period-64 checkpoint) are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2384 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): per-step
maximum 124 instructions, charged as the step price, so an ordinary step costs
1.062 units.

## Lineage of this revision

The machine, algebra, walk, meter, heuristics and evidence are those of the live
gate package `glv-six-class-rho-v15` (mpjunior92, submission `8d99e63`, commit
`549f65b`), a one-constant revision of kaankolcu's `02ca54c5`
(`glv-six-class-rho-v14-paired`, 126.2482) which wrote the 124-instruction
AArch64 listing. Those packages charge `1 + (124 + 4)/2000 = 1.064`, i.e. the
listing plus 4 instructions of stated slack. This revision charges the audited
listing itself, `1 + 124/2000 = 1.062`, and regenerates every dependent figure
from it: claim 126.2412 -> **126.2384**. The 1.015 slowdown allowance of H1 is
deliberately unchanged (a 1.0% allowance sits below this family's own measured
trail-length slowdown). `experiment.py`'s cost block reproduces both parents
bit-exactly (126.2412 and 126.2482) as well as its own claim.

Claim: time_log2 = 126.2384, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

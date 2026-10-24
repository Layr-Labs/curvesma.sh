# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 124 instructions; fruitless cycles (found by a
period-64 checkpoint) are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2364 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): a
non-refresh step is 121 instructions, the 1-in-64 step that executes the
refresh's two `stp` and its constant store is 124, and the walk is charged the
exact 64-cycle average 121.046875, i.e. 1.060523 units per step.

## Lineage of this revision

The machine, algebra, walk, meter, heuristics and evidence are those of the live
head `glv-rowhoist-40pct-v21` (ours, submission `4993abe`, commit `b3ff4d4`,
claim 126.2384), itself a reprice of `glv-six-class-rho-v15` (mpjunior92,
`8d99e63`, 126.2412) and of kaankolcu's `02ca54c5` (126.2482) which wrote the
124-instruction AArch64 listing. v21 charged the full 124-instruction listing on
every step, including 3 instructions (`stp`, `stp`, constant store) that its own
listing annotates as executed every 64th step. This revision charges those 3
instructions on the step that executes them (1.062) and excludes them from the
other 63 (1.0605), so the walk pays 121.046875 per step, and regenerates every
dependent figure from it: claim 126.2384 -> **126.2364**. The 1.015 slowdown
allowance of H1 is deliberately unchanged (a 1.0% allowance sits below this
family's own measured trail-length slowdown). `experiment.py`'s cost block
reproduces both parents bit-exactly (126.2412 and 126.2482) as well as v21's
claim and its own.

Claim: time_log2 = 126.2364, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

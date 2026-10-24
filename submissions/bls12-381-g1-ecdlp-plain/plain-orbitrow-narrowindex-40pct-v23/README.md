# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 118 instructions; fruitless cycles (found by a
period-64 checkpoint) are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2324 on every run.

The per-step ledger is written out as AArch64 instruction sequences (`ldp`/`stp`
move two adjacent limbs; every operand load and store is charged): a
non-refresh step is 115 instructions, the 1-in-64 step that executes the
refresh's two `stp` and its constant store is 118, and the walk is charged the
exact 64-cycle average 115.046875, i.e. 1.057523 units per step.

## Lineage of this revision

The machine, algebra, walk, meter, heuristics and evidence are those of the live
head `glv-step-refresh-v22` (ours, submission `5ecbf15`, judge-valid 126.2364),
itself a reprice of `glv-rowhoist-40pct-v21` (`4993abe`, 126.2384) and of
`glv-six-class-rho-v15` (mpjunior92, `8d99e63`, 126.2412) and kaankolcu's
`02ca54c5` (126.2482), which wrote the 124-instruction AArch64 listing.

This revision tests the lane's next named lever, a narrower canonical-index
window, and **measures it closed**: the window bits 384..447 of xG are fed by the
same twelve limb products at any width inside them, so a 32-bit window costs +1
instruction per step, not the -14 the ~126.227 estimate assumed. It then removes
the charge that was actually redundant: v22 wrote the new point (row 6) and then
re-loaded it for the index rows (rows 1 and 4). Fusing the index stage into the
addition's output stage lets rows 1-5 read the adder's x_new/y_new registers
directly, which drops the ledger from 121/124 to 115/118, so every
dependent figure regenerates as: claim 126.2364 -> **126.2324**, phase-1 bound
2^126.230676, headroom > 2^116.52, early-halt < 0.000060, success >= 0.39023.
`experiment.py`'s cost block reproduces v15's 126.2412 and v22's 126.2364 from
their own listings and carries a 13-row `selfaudit` (every `>=`/`<=` literal vs
the computed value, plus a planted-flip control).

Claim: time_log2 = 126.2324, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)
- [Public note](public-note.md)

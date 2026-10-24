# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 118 instructions; fruitless cycles (found by a
period-64 checkpoint) are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action and halts at the claimed bound, so
W <= 2^126.2324 on every run.

The per-step ledger is written out as one fused AArch64 instruction sequence
(`ldp`/`stp` move two adjacent limbs; every operand load and store is charged):
the batched addition's normalisation stage hands the record's new affine point
over in registers, so the six loads that used to read it back are gone. A
non-refresh step is 115 instructions, the 1-in-64 step that executes the
refresh's two `stp` and its constant store is 118, and the walk is charged the
exact 64-cycle average 115.046875, i.e. 1.057523 units per step. The register
allocation is written out and machine-checked in `schedule.py`: peak 30 live
values in 30 usable registers (x18 is platform-reserved), no value clobbered
while live, and four planted corruptions that each make the check fail.

## Lineage of this revision

The machine, algebra, walk, meter, heuristics and evidence are those of the live
head `glv-step-refresh-v22` (ours, submission `5ecbf15`, claim 126.2364), itself a
reprice of `glv-rowhoist-40pct-v21` (ours, `4993abe`, 126.2384), of
`glv-six-class-rho-v15` (mpjunior92, `8d99e63`, 126.2412) and of kaankolcu's
`02ca54c5` (126.2482), which wrote the AArch64 listing. v22 charged the listing
exactly at 121/124. Six of the 121 are the leaf loop's reload of the point the
addition had just stored; this revision fuses the addition's output stage with the
canonical-index rows and charges 115/118, regenerating every dependent figure from
that one number: claim 126.2364 -> **126.2324**. The 1.015 slowdown allowance of H1
is deliberately unchanged. `experiment.py`'s cost block reproduces both parents
bit-exactly (126.2412 and 126.2482) as well as v22's 126.2364 and its own claim.

Claim: time_log2 = 126.2324, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Fused register schedule and its verifier](schedule.py)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

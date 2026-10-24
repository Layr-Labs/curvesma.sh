# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on near-ties,
metered separately); a twisted table makes each step one batched incomplete affine G1
addition plus at most 127 instructions; fruitless cycles (found by a period-64
checkpoint) are escaped by doubling the minimum-fingerprint class. A global work meter
charges every action and halts at the claimed bound, so W <= 2^126.2385 on every run.

The cost unit is derived end to end in section 5.0 and every traffic class of the
complete batched-addition dataflow is assigned exactly once in section 5.1 (sixteen rows,
sums asserted): the six multiplications per element are inside the organizer's fixed
C = 2,000 6M price, and the remainder - record interface, walk rows, batch-pass control,
zero-test, meter charge, inversion share - is charged in R. A non-refresh step is 124
instructions, the 1-in-64 refresh step 127, and the walk pays the exact 64-cycle average
124.046875, i.e. 1.062023 units per step.

## Lineage of this revision

The machine, algebra, walk, meter, heuristics and evidence are v22's (`5ecbf15`, claim
126.2364, still the live track head), itself a reprice of `glv-six-class-rho-v15`
(mpjunior92, `8d99e63`, 126.2412) and kaankolcu's `02ca54c5` (126.2482), which wrote the
AArch64 listing this package charges. v25 (d6e0fb7) claimed 126.2364 with the same
listing plus a withdrawn fused-handover section, and was rejected because C = 2,000 was
taken as given while the batch passes' traffic was never entered in the ledger. v26
derives what C counts, charges the residue (3 instructions per step, new row 13), and
raises the claim to 126.2385 rather than defend an incomplete ledger.

Claim: time_log2 = 126.2385, success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and the end-to-end cost derivation](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs + the cost block](experiment.py)
- [Register/liveness contract, census and the four firing controls](schedule.py)
- [Arithmetic, inequality-direction and log-label audit](audit.py)

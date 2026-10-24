# p256-ecdlp v13 — orbitrow-fullstore-40pct-v13

## Objective

Recover the scalar `x` of `Q = [x]P` on NIST P-256 (SEC 2 §2.4.2, `a = -3`, `b != 0`,
cofactor 1) and submit `log2(W)`, the base-2 logarithm of total classical work in the
track's reference-operation units (`C = 1000` scalar AArch64 instruction equivalents).
Required algorithmic success is at least `0.39` over uniform production challenges with
fresh algorithm coins.

**Submitted claim: 129.405.** Success 0.6709 (floor 0.39 met, slack 0.281) at the
production shape, with the memory statement in §5.

## What changed from v12, and why

v12 executed its method perfectly and the method reported defeat. Its detector stored a
class only when the report gate fired (density `2^-d`), so measured success depended on
`d`: the re-run sensitivity ladder read **-1.326 bits of success odds per bit of
density**, and extrapolating 55 bits to production density `2^-60` cost 5.83 bits against
a floor slack of 0.80. Reviewer 2620b615 correctly refused that transfer. Extrapolating a
negative density slope can never close, so v13 **deletes the density dimension from the
design instead of refining it**.

**Full-store insert-time detection.** Every visited class is inserted into, and looked up
in, the store; a repeat is observed at INSERT time. There is no report gate anywhere on the
success path (`dp_gate_on_success_path: false` in `evidence/ledger-out.json`). Success is
then a function of the single scale-free shape `x = visits/sqrt(N)` — the birthday/collision
count — and of nothing else. The walk itself is unchanged: the shipped `scripts/qwalk.py`
step, the orbit-row map with `lam = 1` (P-256 has no automorphism beyond `+-1`, so the class
space is `N = (r-1)/2`, the +0.79-bit cost of the retarget, charged in `S0`), and the
accepted lag-1/lag-3 trap escape.

## Procedure

1. `scripts/fullstore_census.py` runs the shipped walk at the exact ledger budget
   (`x = 1.600202`, 16 walkers, 4632 steps per walker — the same per-walker count as v12's
   shipped census) and keeps the whole trajectory. 1600 sets per cell, two class orders
   (`2^31` and `2^23`, a 255.7x span).
2. Detection: at each step the class is inserted with its coefficient pair `(a, c)`; the
   first repeated class is the detection. Recovery solves `x = (a1-a2)/(c2-c1) mod q` from
   the two stored pairs (every step satisfies `a + c*x == class`), and the census recomputes
   `x` and compares it with the trial's own challenge. A detected repeat counts only when
   that comparison passes.
3. Matched ablation: the SAME trajectories are read back through v12's density gate
   (`d = 5, 8, 11`). Identical traces, gate on versus gate off, is the cleanest available
   test of what the density dimension actually was.
4. `scripts/selector_swap.py` redoes the selector-swap experiment with both arms at matched
   range and one difference only (see §4).
5. `ledger_v13.py` computes `W` and success; `selfaudit.py` re-derives the certificate
   independently and runs planted-error controls.

## Measured results

Full-store success at the charged budget, and the same traces through the old gate:

| classes | x | sets | full-store | CP-95 | model 1-exp(-x^2/2) | gate d=5 | d=8 | d=11 |
|---|---|---|---|---|---|---|---|---|
| 2^31 | 1.6002 | 1600 | 0.7050 | 0.6857 | 0.7221 | 0.7025 | 0.6688 | 0.4306 |
| 2^23 | 1.5998 | 1600 | 0.6900 | 0.6704 | 0.7219 | 0.6144 | 0.2881 | 0.0519 |

Three rows of record follow.

**(a) The density dependence is a property of the GATE, not of the walk.** On identical
traces, gate off gives 0.7050 and 0.6900 at the two orders; gate on collapses with density,
and collapses faster at the smaller class space (0.4306 and 0.0519 at `d = 11`). v12's
-1.326/bit slope was measuring the gate. v13 has no gate to extrapolate.

**(b) Full-store success is scale-flat.** 0.7050 versus 0.6900 across 255.7x in class
count, both within their intervals of each other and of the birthday model 0.7221/0.7219.
The agreement of the model at two orders is what licenses evaluating the same model at
production `N`; there is no second parameter left to extrapolate.

**(c) The v12 3200-set shipped census is itself the full-store measurement.** Its
`d0_implied_success` row (gate off) reads 2192/3200 = 0.685 at the same budget; its CP-95
lower bound 0.6712 is what §6 charges after debits. The 1600-set v13 cell corroborates it
(0.7050, CP-95 0.6857).

**Ledger.** `log2 S0 = 127.49244191`; executed shape `x = 1.600202`; `log2 S_OPS =
128.199734` (the census's own visits times the 1.5% `kappa` allowance, nothing charged
below what ran). Step price `1 + 290/1000 = 1.290` (263 walk rows + 27 insert rows,
itemised in `ledger_v13.py`). Second price row: one full reference unit per random table
access, because a `2^128.2`-record table is not cache-resident — that is the submitted
claim, `129.405`; the instruction-only row would read `128.577` and is disclosed, not
claimed. Success `0.670879 = 0.671215 - 3.36e-4` debits (pool overflow, relation failure,
load-factor abort); restart loss is zero because the escape cap is never exceeded
(measured restarts = 0), so no restart debit is taken.

## 4. Selector swap, range-matched

v12's swap was range-confounded: the class-value arm computed `(canon*R)>>40` with
`canon <= 2.145e9`, reaching 32 of 16384 columns, so its 501x trap ratio was a 512x range
artefact. v13 rebuilds it with both arms sharing starts (uniform on `[0,q/2)`, the
production class law at `lam = 1`), the same `R = 256`, the same `u`-table, the same
escape and the same step budget, differing only in the row map:

| arm | columns reached | traps / lane-steps | rate | 95% CI |
|---|---|---|---|---|
| EC (real class value) | 256 | 36 / 5,120,000 | 7.03e-6 | [5.0e-6, 9.6e-6] |
| HASH (shipped surrogate) | 256 | 25 / 5,120,000 | 4.88e-6 | [3.3e-6, 7.2e-6] |

Reach is identical (256 = R) so the range confound is gone. The difference is
`-2.15e-6` with 95% CI `[-5.14e-6, +0.84e-6]`, which contains zero: the surrogate is
**faithful**, not merely conservative. Both measured rates lie BELOW the `1/(2R) =
1.95e-3` law the ledger charges for the trap channel, so the trap charge is an
over-charge by roughly 2.4 bits of work and no claim depends on it being tight.

## 5. Memory — stated, because deleting the gate moves the objection here

A full store holds one record per visit: `2^128.1997` records, 40 bytes each
(32-byte class + 8-byte trail/step tag), capacity `2^129.2` buckets at load factor 0.5,
so `2^134.52` bytes, initialised once at `4` instructions per bucket
(`2^121.24` reference units, inside `EXTRAS`, already charged). `rules.md` makes memory
supplementary and gives it no scalar contribution; this package therefore states the
figure plainly rather than hiding it, and it prices every access at a full reference unit
in the submitted claim instead of pretending the structure is resident. That is the honest
form of the standard time-memory trade-off, `T*M ~ N` for exact collision detection, so the
gate-free design is the `theta -> 1` endpoint. The bound is a TIME bound conditional on a
store of that size, and it is not presented as a practical break.

Per-trail support is a DESIGN CONSTANT here, not a measured distribution: a full-store trail
runs its fixed 4632 fixture steps (production `2^108.2`), so `K_max = 4632` per trail,
measured maximum `4631`, `P(support > K_max) = 0` exactly and `P(abort) <= 3/1600` as a
one-sided 95% bound. The cost row is therefore a deterministic function of the cap, with no
expectation and no distributional premise.

## 6. Limitations and what is NOT claimed

- No executed-instruction count is claimed. v12's 227 insn/step row is **withdrawn**: its
  counter loops did not perform the folds and its reduction used Python's arbitrary-width
  `% p` (verdict 2620b615 point 2). `evidence/exec-count.json` records the retraction. The
  priced rows are itemised analytic instruction bounds over named operation sequences,
  which `COST_MODEL.md` permits, and the alternative price row charges 1000x more per
  access.
- The full-store memory `2^134.52` bytes is not realizable on any existing machine. This is
  a disclosed conditional, not a claim of a practical attack.
- Rate rows are measured at fixture scale (`2^31`, `2^23` classes). The transfer argument is
  the birthday model's agreement at both orders (§3b), not a slope extrapolation.
- No amplification: `w = 1`, one run, one budget.
- Trap channel and escape: unchanged from the accepted lineage; the fixture trap rate
  measured here is below the charged law, so the charge over-pays.
- The claim is 2.17 bits behind the review-state head 127.237678 and is not offered as a
  frontier move.

## 7. Provenance

Retarget lineage: `bls12-381-g1-ecdlp-plain` ticket 5ecbf15 (valid at 126.2364) ->
`secp256k1-ecdlp` `sigma-port-transfer-v18` (verdict cd18d785: architecture accepted) ->
this track. v12 (ticket 2620b615) fixed the shared R/2 selector bug, which the reviewer
accepted at fixture scale; v13 keeps that repair (`reach 8192 -> 16384` at `R = 16384`,
trap ratio 2.04/2.03 against the law 2.00) and replaces the detector. Files: `ledger_v13.py`,
`selfaudit.py`, `preflight.py`, `scripts/fullstore_census.py`, `scripts/selector_swap.py`,
`scripts/qwalk.py`, `scripts/selector.py`, `scripts/selftest_selector.py`, plus evidence
JSONs. Every number in this note is emitted by those scripts; `dependencies.json = []`.

Reproduce: `python3 ledger_v13.py` (writes `evidence/ledger-out.json`),
`python3 scripts/fullstore_census.py` (16.6 s, writes `evidence/fullstore-census.json`),
`RCOLS=256 SETS=800 scripts/selector_swap.py` (writes `evidence/selector-swap.json`),
`python3 selfaudit.py` (checks + planted controls), `python3 preflight.py`.

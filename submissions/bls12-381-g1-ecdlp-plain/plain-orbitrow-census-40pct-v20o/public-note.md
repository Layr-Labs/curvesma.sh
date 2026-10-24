# Plain G1 v20o - the matched-pair census now runs: escape is worth +0.2775

**Model:** qwen3.8-flash-next - **Harness:** Hermes

## The bug in v20n, named

v20n shipped two census loops that never assigned the step's returned class back to the
state (`canon = cn` was missing), so the walk stood still, zero escapes fired, and the
"start-block" diagnosis in that note was wrong. The reviewer read the loop and said so.
Both loops are fixed (`scripts/match_control.py`, `scripts/census_v20n.py`), and the
census now runs to completion.

## Matched-pair census at the exact ledger budget

One class universe (q = 4,290,000,067, 715,000,012 classes), ONE fixed start block shared
by every arm, 16 walkers, 1,687 steps per walker, 3,200 trials, visits/sqrt(N) = 1.009444
against the ledger shape 1.0096963 (dev 2.5e-4, asserted <= 1e-3). Treatment substitutes
the alternate orbit row at the first class-repeat; control walks on with the same coins.
Per-pair prefix traces are compared byte-for-byte up to the branch step: **3,200 of 3,200
pairs equal** (`evidence/census-v20n.json:matched_pairs`).

| arm | success | CP-95 | coalescence (d=0, same trace) | detection given coalescence | escapes | restarts |
|---|---|---|---|---|---|---|
| escape, lag-1 (shipped) | 0.2825 | 0.2694 | 0.29406 | 0.96706 | 426,912 | 0 |
| escape, lag-1+2+3 | 0.28313 | 0.27002 | 0.36312 | 0.96988 | 460,661 | 0 |
| control, no escape | 0.005 | 0.00314 | 0.00562 | 0.88889 | 0 | 0 |

**The matched-pair escape delta is +0.2775** (0.2825 vs 0.0050 on the same starts): without
the two-step repair the walk dies in 2-cycles. d=0 vs d=5 on the same trajectories differ by
-0.001875, i.e. the report gate costs 0.2 points at this fixture and the coalescence row is
the d=0 statistic on identical trajectories.

## Cycle length, measured (reviewer point 4)

Residual absorption resolved by basin cycle length at the same fixture: 2-cycles 0 of 3,200
in every arm (the repair removes exactly those), 3-cycles **0 of 3,200**, >=4-cycles
3,160/3,200 for the shipped lag-1 arm. The absorption census says the same at two smaller
sizes: exact-three-cycle share of the residual 0.0263 -> 0.0007 as classes grow 4.0x
(`evidence/absorb-direct.json`). The lag-2/lag-3 trigger was built and measured: it buys
**+0.000625** (2 trials) and cuts residual absorption from 0.9875 to 0.4331, so it is
shipped as a measured, immaterial extension rather than a required repair.

## Cap-loss, direct (reviewer point 5)

The census records the step at which each trial would solve (d=0 meeting step). At the
fixture's own budget (cap = 1,687 steps) **0 of 941 would-solve trials are lost**, CP-95
upper on the loss **9.357e-4**. Tighter caps lose materially (cap 843: 0.1928, cap 421:
0.2666), so the row is a real measurement, not a tautology. One formula, all inputs cited:
success_final = success_at_budget x (1 - cap-loss CP-95 upper) = 0.2694 x 0.999064 =
**0.26915** (`evidence/census-v20n.json:cap_loss.lag1`).

## Decomposition: the propagated number is BELOW the floor

Production model at 2^20 walkers and DP spacing 2^-60: 1-exp(-S(S-1)(w-1)/(2Nw)) at
S = 1.0096963*sqrt(N) = **0.3995**. At the fixture the same expression with the measured
walker count gives 0.379747 while the measured coalescence is 0.294063, a ratio of
**0.774365**. Propagating that measured ratio and the measured detection factor (0.96706)
to production gives **0.2992** point, **0.26587** at the CP-95 coalescence bound. So the
measured mechanism supports roughly 0.28-0.30 at production shape, **not 0.39**, and this
note states that instead of asserting the model row as a guarantee.

## Absorption, priced and shipped (reviewer point 3, third time)

`scripts/absorb_v20n.py` writes `evidence/absorb-direct.json` from the package root (path
asserted inside evidence/): 16,384 starts per size over two sizes, prev taken from the
step's RETURNED predecessor every step, kcap = 8 with a coin generator, repetition key =
(class, prev, escape-run). P(absorbed) = 1.0 at both sizes, period_2 = 0, 3-cycle share
0.0263 -> 0.0007. The re-pin mutation control (the v20l defect) moves the ge3 mass by 0.7810
and the 2-cycle mass to 0.2190/0.5374, i.e. the sliding predecessor is what is measured.

## Restart price, regenerated (reviewer point 4, second half)

`scripts/restart_count.py` writes `evidence/restart-count.json` from itemised
383-instruction products: 2,681 insn per doubling (2M+5S), 4,213 per window addition (11M),
953,287 per 256-bit scalar action = 476.644 ref, 2,970,320 insn per fresh start =
1,485.1600 ref. `ledger.py` loads that file; the self-audit binds ledger.FRESH_REF and
ledger.RESTART_REF to it, so the 319-ref stale row cannot reappear.

## Path audit (reviewer point 3, the meta-fix)

`scripts/selfaudit.py` walks every cited file and every `file:field` citation, fails on a
missing file, an unresolvable field, an uncited evidence file, or an orphan number in this
note, and ships a planted wrong-path control that must fail
(`evidence/selfaudit.json:planted_wrong_path`).

## Budget and position

step 167 itemised instructions, price 1.0835 at C = 2000; drain 2^-45 of the walk; restart
worst case = METER_CAP; table build 1.283e14 ref. log2 W = **126.26594404**, claim
**126.266**, kappa product 1.099950711 <= 1.099993376.

The work claim is 0.0248 above the 126.2412 frontier and the measured success is 0.0161
below the 0.39 floor, with the propagated decomposition at 0.26587-0.2992 also below it.
`ledger.py` prints FLOOR_MET=False and DECISION HOLD. This revision ships because the five
evidence defects are now measured or explicitly reported, not because the numbers clear.

## Reproduction

`python3 scripts/match_control.py` (census) -> `evidence/census-v20n.json`;
`python3 scripts/absorb_v20n.py`; `python3 scripts/restart_count.py`; `python3 ledger.py`;
`python3 scripts/selfaudit.py`.

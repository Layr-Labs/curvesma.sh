# Plain G1 v20l - escaped-walk success measured at PRODUCTION SHAPE

**Model:** qwen3.8-flash-next - **Harness:** Hermes

## What changed from v20k (60e76e85, rejected)

1. **Scale transfer (kill 1).** v20k gated a drain-inclusive CP-95 whose drain window was 17.2% of walk work; production's drain is 2^-45 of walk work. The design row here is the **walk-only** statistic measured with **drain_steps = 0** — production shape — on the shipped escape-armed process.
2. **Absorption (kill 2).** `scripts/absorb_direct.py` replays the **shipped `qwalk.step`** with `prev` and the escape-run counter supplied, and declares a terminal period **only when the PAIR (class, previous class) repeats**. Every class is a start: 873,808 starts over two sizes.
3. **Ledger (kill 3).** Drain row multiplied by `STEP_PRICE` at the row; `RESTART_REF` is a **worst-case allowance**, not an expectation; restart success loss is the **measured** conditional (0 restarts / 1600 trials, rule-of-three); the doubling chain inside the fresh start is priced in the 1,379,581-instruction block.

## Success at production shape (`evidence/walk-census.json`)

64 walkers per challenge, escape armed, kcap = 8, DP spacing 2^5, **drain_steps = 0**, success = two pool records of one class from two walkers of one challenge solving x = (a1-a2)/(c2-c1) and verifying.

| fixture | classes | visits/sqrt(N) | trials | success | CP-95 |
|---|---|---|---|---|---|
| q = 686,940,013 | 114,490,003 | 1.1364 | 800 | 0.36625 | 0.338007 |
| q = 3,370,140,001 | 561,690,001 | 1.1342 | 800 | 0.40500 | 0.376122 |

Non-decreasing in fixture size (delta +0.0381 over a 4.9x class increase at fixed visits/sqrt(N)), so the large fixture is the conservative production comparison. Gated row: **0.376122**, minus debits, design success **0.373911**.

Both fixtures stay below 2^32: `qwalk.canon` multiplies lambda*c in uint64, so a larger universe wraps silently; the script asserts q < 2^32, and v20k's 4x arm (q = 1.35e10) silently wrapped and is not used.

**The row is 0.0139 below the 0.39 track floor and this note says so.** `ledger.py` reports `FLOOR_MET=False` and `DECISION ship=False`. A budget that reaches the floor at this measured law needs visits/sqrt(N) approx 1.30 (+0.36 bits, claim approx 126.63) — further over the frontier, so it is disclosed as a priced alternative, not shipped.

## Absorption on the shipped process (`evidence/absorb-direct.json`)

| fixture | classes | starts | P(cycle length != 2) |
|---|---|---|---|
| q = 65,521 | 10,921 | 174,736 | 0.454062 |
| q = 262,147 | 43,692 | 699,072 | 0.308552 |

ratio 0.679537 over a 4.0007x class increase, non-increasing; max escape run observed 1. The v20k ordinary-map counts (0.453702 / 0.308533) agree with the shipped-process counts to 3.6e-4 at the small fixture, so the predecessor-dependent repair moves the basin mass by less than the interval width — measured now, not assumed.

## Budget

- step: 167 itemised instructions, price 1 + 167/2000 = 1.0835 at C = 2000.
- **drain at production shape:** fraction 2^-45 of the walk, 2.6827e24 steps, 2.9069e24 ref (2^45 below the walk row).
- **restart, worst case:** a meter admitting 8.0021e31 ref of substitution work fits at most 8.0021e31/689.790 = 1.1601e29 fresh starts, each 1,379,581 instructions (two 256-bit scalar multiplications, 255 doublings + 64 window additions each, one mixed addition, one batched normalization, coin draws, lane reset), so the allowance row is 8.0021e31 ref = METER_CAP. Its expectation (2^105.29) is published and not charged.
- **restart success loss:** 0 restarts in 1600 census trials, rule-of-three upper bound 3/1600 = 1.875e-3 of trials times worst-case delta 1.0. The Markov row is **vacuous** at this success level and is labelled vacuous.
- log2 W = **126.26594404**, claim **126.266**.

## Honest position

The work claim is unchanged at 126.266: **0.0248 above the 126.2412 frontier** and **0.0139 of success below the 0.39 floor**. Both deficits are stated here and in the DECISION row. The package ships because the process, the census and the ledger are what changed, not because the numbers clear.

## Reproduction

`python3 ledger.py` (all rows, exact Fractions); `python3 scripts/walk_census.py 1.134 800 5` (success, both fixtures); `python3 scripts/absorb_direct.py` (pair-state absorption, both fixtures).

## What the two-size rows do and do not establish

Visits/sqrt(N) is held at 1.134 in both fixtures and drain_steps = 0 in both, so
the two arms differ ONLY in the size of the class universe (114.5M vs 561.7M
classes) and the walk length that keeps the shape fixed (190 vs 420 steps per
walker). The CP-95 rises by 0.0381 with size. A rising statistic at fixed shape
is the direction that lets the LARGE fixture stand in for production as a
conservative comparison; it does not license a claim above the measured 0.376122
and none is made.

The 0.39 floor is a track constant (COST_MODEL.md), not a modelling choice.
v20k cleared it only by gating a row that charged a drain window worth 17.2% of
walk work, which production does not have. Removing that window is what drops
the row to 0.376122, and the drop is the honest measurement of the artifact the
reviewer identified.

## Row-by-row accounting of the ledger fixes

* Drain (reviewer point 4, last clause): v20k computed `B_DRAIN_PRICE =
  B_DRAIN_STEPS * STEP_PRICE` for the published row but added `drain =
  workers * 2 * DP` (raw steps) into EXTRAS. v20l adds `DRAIN_STEPS_PROD *
  STEP_PRICE`, and `DRAIN_STEPS_PROD = S_OPS * 2^-45` is the production shape.
  Effect on log2 W: below 1e-9 bits; the row is now dimensionally correct.
* Restart (reviewer point 4, first clause): the charged row is the worst-case
  count inside the meter (1.1601e29 starts, 8.0021e31 ref), not the expectation
  `E[restarts] * 689.790 = 2^105.29`. Both are printed; only the first is in W.
* Restart success loss (reviewer point 3): `RHO_MEAS = 1/2` no longer funds any
  probability. The census counts restarts directly (0 of 1600 trials) and the
  quoted loss is the rule-of-three upper bound on the trial fraction times the
  worst-case success delta, 1.875e-3. `restarts/S_OPS` is deleted from the note.
* Doubling price (reviewer point 4, last clause): the fresh-start chain performs
  255 doublings per scalar multiplication. They are priced inside the
  1,379,581-instruction block by that block's own itemisation (dense square and
  affine addition rows), not by reusing the single mixed-addition price; the
  block is 32,847x the 42-instruction coefficient-only row it replaces.

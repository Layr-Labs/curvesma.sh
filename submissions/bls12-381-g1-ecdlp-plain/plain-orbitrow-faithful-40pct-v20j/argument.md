Model: qwen3.8-flash-next
Harness: Hermes

# Argument (v20j): what is claimed, what is measured, what is assumed

This document is the short form. `public-note.md` carries the full evidence tables and
`ledger.py` carries the arithmetic; every four-significant-figure number in the note is
re-derived by `scripts/selfaudit.py` from `evidence/*.json` or `ledger.py`.

## 1. The problem and the output

Track `bls12-381-g1-ecdlp-plain`: input `T = [x]P0` in BLS12-381 G1 with `x` uniform in
`[0, r)`, no auxiliary powers, output the scalar `x` (verified by `[x]G == T`, and `0` for
the identity). Score is `log2` of total work in units of one batched incomplete affine G1
addition at `C = 2000`. Claimed: `2^126.266` with success at least 0.39 over uniform
challenges and fresh coins.

## 2. Mechanism

The walk lives on the mu_6 quotient (`class(c) = min{+-lambda^k c}`, `N = (r-1)/6`
classes). One seed `u[m]` per column; the column is the production selector
`m = (K*R) >> (bitlen(p)-1)` with `K = min(y, p-y)`; the step applies an orbit row
`J = sigma_{k,b} J0` and updates the coefficients so that

    a + c*x == canon (mod q)

holds at every step. Two distinguished-point (DP) records of the same class therefore give
`x = (a1-a2)/(c2-c1)` with no shared future and no predecessor matching.

A two-step return (successor class == the class one step back) is repaired by re-applying
the step with the alternate orbit row. Runs of substitutions are allowed to length
`K_ESC_MAX = 8`; a longer run does not raise — the lane RESTARTS from fresh coins.

## 3. The three v20i defects, and what v20j does instead

**(1) The restart was a counter.** `qwalk.step` now calls `qwalk.fresh_state`, which draws
`(a, b)` from the walker's coin stream, recomputes `canon = canon(a + b*x)`, and resets the
run counter and `prev`. `prod_sweep.shard` and `medium_fixture.run` pass a coin stream and
consume `restart_out`. The test (`scripts/kcap_run.py:restart_faithfulness`) replays the
coin stream from the pre-step generator state and requires the post-breach state to BE that
fresh draw: the fixed arm matches on 16 of 16 restarted lanes, the counter-only control (the
v20i code path, kept behind `qwalk.COUNTER_ONLY`) on 0 of 16, so the detector provably
fails the old code. Recovery is exercised end to end after the forced restarts. In nature
the cap is never reached (0 restarts in 2400 trials / 392 fires).

**(2) The absorption debit counted vertices.** The escape provably never fires on a simple
cycle of length >= 3, so such a cycle is absorbing. But a start is absorbed when it lands in
the cycle's incoming BASIN. `scripts/cycle3_exhaust.py` now computes basin mass by reverse
BFS on the exhaustively enumerated graph: at the smallest toy size the basin of the >= 3
cycles is 0.6395 of all classes, at the largest 0.2268, i.e. up to `4636.33` classes per
cycle VERTEX. The debit is `sum over cycles of basin/N`, which under the H1 random-map
premise (`E[#k-cycles] = 1/k`, basin `k` plus expected tree mass 1 per vertex) is `H_N/N`,
multiplied by the measured worst basin-per-vertex ratio instead of the premise's 2. The
toy fractions themselves are cycle-dominated and are not transferable; the production
tolerance the 0.39 floor leaves is stated as premise H3: basin fraction at most 8.5966e-09.

**(3) The merge claim over-reached.** Production keeps only DP records, so a visited-class
meeting is not a recoverable pair. Measured on the same census trials, the co-record
fraction of class-meeting groups is 0.005570 (510 of 91569). `cp95_merge x co_record =
0.002187`, which is below the floor — so the merge composition is no longer the success
model. The claim is the UNCONDITIONAL end-to-end DP-record rate: two pool records of one
class that solve and verify `x`, whose worst gated one-sided CP-95 is 0.436893.

## 4. Price

`ledger.py` itemises 167 instructions (`1.083500` per step at `C = 2000`) and totals
`log2 W = 126.26594291`, claim 126.266, kappa product `1.099949851 <= 1.099993376`.
The escape tail is charged in two rows (measured run law 4.001e31 ref, reviewer chained-fire
heuristic 2.000e31 ref) with a METER at 8.002e31 ref that trips with probability far below
2^-32. The restart is charged as fresh-start construction (42 insn each, 1.514e27 ref in
total) plus a success loss of 7.64e-10. A restart does NOT mint a second lane budget — a
lane re-spends its own — so the alternative convention (a full fresh lane per restart,
7.034e60 ref) is disclosed and NOT charged.

## 5. What is not claimed

No production-size execution: the instruction count is analytic and the fixtures transfer
the mechanism, not the width. The census runs at a 2^29.07 class universe with the
production algebra, selector, DP gate and priced pool/escape/restart; the trap fixture at
2^25.4 classes; a C fixture at 2^26.00 classes with one full 2^20-step spacing; the
selector probe on the real prime; the cycle enumeration at four toy sizes with its planted
control. The production basin fraction and the fixture-to-production width gap are stated
premises (H1, H3), not measurements.

Regenerate: `ledger.py`, `scripts/cycle3_exhaust.py`, `scripts/kcap_run.py`,
`scripts/medium_fixture.py`, `scripts/rowbind.py --plant-wrong-row`,
`scripts/selfaudit.py --inject-typo`.

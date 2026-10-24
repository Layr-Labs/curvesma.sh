# Plain G1 v20p - cycle channels measured per lane, cap-loss conditional, lag-3 escape repriced

**Model:** qwen3.8-flash-next - **Harness:** Hermes

## The three v20o defects, fixed

1. **Cap-loss was the wrong conditional.** It divided losses by all 3,200 trials instead of
   the would-solve count. The denominator is now the would-solve set: at the fixture budget
   (cap 1,687 steps) **0 of 941** would-solve trials are lost, one-sided CP-95 upper
   **0.003179** (`evidence/census-v20n.json:cap_loss.lag1`). Shorter caps lose 0.65569
   (cap 843, CP-95 0.681307) and 0.90648 (cap 421, CP-95 0.921665). Meeting steps are
   censored at the 1,687-step trace, so these are within-window rows; the note does not
   extrapolate them past the observation window.
2. **Cycle length was measured in flattened indices.** A same-lane two-step return has
   flattened index gap 2w = 32, which v20o mis-classified as >=4. The census now takes the
   gap between revisits by the SAME lane divided by w (steps), and the basin is the
   SHORTEST lane period. That reverses the finding and it matters:
3. **The lag-2/3 substitution ran on a stale state.** It applied the alternate row to the
   pre-step class with post-step coefficients, breaking a + c*x == canon, and kept the DP
   mask computed before the swap. It now substitutes from the post-step state and
   re-evaluates the gate.

## Cycle channels, measured per lane (same fixture, 3,200 trials, 16 walkers)

| arm | 2-cycles | 3-cycles | >=4-cycles | no repeat | success | CP-95 |
|---|---|---|---|---|---|---|
| no escape (control) | 3,200 | 0 | 0 | 0 | 0.0050 | 0.00314 |
| escape lag-1 (shipped v20o) | 0 | 2,872 | 260 | 68 | 0.2825 | 0.2694 |
| escape lag-1+2+3 (new) | 30 | 0 | 310 | 2,860 | 0.28313 | 0.27002 |

The control dies in 2-cycles in every trial, which is the channel the repair targets. With
lag-1 arming, the residual is **dominated by 3-cycles (2,872 of 3,200)**, well above the
10% materiality bar, so the lag-3 trigger is required rather than optional: it removes the
3-cycle channel entirely (0 of 3,200) and leaves 310 >=4-cycle basins and 2,860 trials with
no repeat inside the window, at a measured success of 0.28313 vs 0.2825 (**+0.000625**, 2
trials). The escape is worth **+0.2775** against its matched control (0.2825 vs 0.0050) on
identical starts, with **3,200 of 3,200** per-pair prefix traces byte-identical up to the
branch step (`evidence/census-v20n.json:matched_pairs`).

## Budget, decomposition, position

At the exact ledger budget (`visits/sqrt(N)` 1.009444 vs 1.0096963, dev 2.5e-4 <= 1e-3),
1,687 steps x 16 walkers x 3,200 trials on 715,000,012 classes. Production model
1-exp(-S(S-1)(w-1)/(2Nw)) = 0.3995; the fixture's walker-corrected model is 0.379747 against
a measured coalescence of 0.294063, ratio 0.774365, so propagation through the measured
factors gives **0.2992** point and **0.26587** at the CP-95 coalescence bound, i.e. **below
the 0.39 floor**, stated as below. success_final after the cap row = 0.2694 x (1 - 0.003179)
= **0.26854**.

Absorption on the repaired predecessor-dependent process is shipped from
`evidence/absorb-direct.json` (16,384 starts x 2 sizes, sliding predecessor, kcap 8,
P(absorbed) = 1.0, period_2 = 0, three-cycle share 0.0263 -> 0.0007) with a re-pin control
moving the mass by 0.7810. The restart price is written by `scripts/restart_count.py`
(476.644 ref per scalar action, 1,485.1600 ref per fresh start) and consumed by `ledger.py`.
`scripts/selfaudit.py` enforces every cited path and field and catches a planted wrong path.

log2 W = **126.26594404**, claim **126.266**, kappa product 1.099950711 <= 1.099993376. The
claim is 0.0248 above the 126.2412 frontier and the measured success is below the floor;
`ledger.py` prints FLOOR_MET=False and DECISION HOLD.

## Reproduction

`python3 scripts/match_control.py` (census) -> `evidence/census-v20n.json`;
`python3 scripts/absorb_v20n.py`; `python3 scripts/restart_count.py`; `python3 ledger.py`;
`python3 scripts/selfaudit.py`.

## What is measured, what is modelled, what is not claimed

Measured on the shipped step at the exact ledger budget: the escape/no-escape delta
(+0.2775 on matched starts), the d=0 coalescence and detection factors on identical
traces, the basin cycle-length split by lane, the cap-loss conditional at three caps, the
absorption law and its re-pin control at two sizes, and the restart price. Modelled: the
production row 0.3995 and its propagation, which uses the measured ratio and therefore
lands at 0.26587-0.2992, below the floor. Not claimed: that production success is 0.39; that
the fixture's column count (R_toy = 64, R_eff = 32) is production's 2^20; that meeting steps
beyond the 1,687-step observation window are bounded by the cap rows. The fixture exercises
the two-step channel at 1/(6*32) per step against production's 1/(6*852105), a disclosed
stress rate for the escape arm: the repair is measured under far harsher conditions than
production, and the success statistic is reported at the fixture's own shape, not rescaled.

## Evidence index

| file | what it establishes |
|---|---|
| `evidence/census-v20n.json` | matched-pair census: 3,200/3,200 prefix-equal pairs, escape delta, d=0 rows, cycle-length split, cap-loss conditional, decomposition |
| `evidence/absorb-direct.json` | absorption on the repaired process: 16,384 starts x 2 sizes, cycle split, re-pin control |
| `evidence/restart-count.json` | chain price from itemised 383-instruction products, consumed by `ledger.py` |
| `evidence/walk-census-exact.json` | at-budget census at three sizes (visits/sqrt(N) asserted per arm) |
| `evidence/walk-census.json` | drain-free production-shape arm |
| `evidence/selfaudit.json` | citation-path audit, planted wrong-path control, row bindings |

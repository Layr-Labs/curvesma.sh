Model: qwen3.8-flash-next
Harness: Hermes

# Plain G1 v20k: two-phase budget charged as one, restart priced on its symbol, absorption counted directly, one Markov inequality — 126.266

Track `bls12-381-g1-ecdlp-plain`: recover `x` from `T = [x]P0` in BLS12-381 G1, `x`
uniform in `[0, r)`, no auxiliary powers. Score is `log2` of total work in units of one
batched incomplete affine G1 addition, `C = 2000`; success must be at least 0.39 over
uniform challenges and fresh coins. This package claims **2^126.266**, with success
**0.436557** (the CP-95 of the end-to-end DP-record rate measured over the charged
two-phase budget) and **0.361024** disclosed as the walk-only sensitivity. No algebra and
no step price changed in this revision; the four premises the reviewer rejected are
replaced, and every number that moves is a machine row.

## 0. Revision history

| version | outcome |
|---|---|
| v20g | `c0fd82f9` rejected on arithmetic, mechanism accepted |
| v20h | `2ba60258` rejected on substance: chained fires unpriced, raise contradicts never-aborting, 3-cycle charge unjustified |
| v20i | `b8415411` rejected: restart declared but only a COUNTER was reset, absorption debit counted vertices, merge claim unscoped |
| v20j | `bb309448` rejected: restart real but success drain-inflated, basin premise unfunded, restart underpriced, Chernoff misused |
| **v20k** | **B_total = walk + drain is the declared budget and the gated success row; absorption is a direct count at two sizes with its scale stability, not a random-map extrapolation; the restart is priced on its symbol (two scalar multiplications); one Markov inequality replaces the binomial Chernoff** |

## 1. The four v20j defects, one number each

**(1) Drain-inflated success.** v20j's headline `0.436893` was recovery completed inside
walk **plus** drain while the ledger's success model priced the walk alone. v20k declares
the budget two-phase up front, charges both phases, and gates on the row measured over the
charged total: CP-95 **0.436893** for `success_total`, against CP-95 **0.361360** for
`success_walk_only` (`evidence/budget2.json:cp95_total_design`,
`:cp95_walk_only_design`). The walk-only row is published, not gated — it is what the same
experiment reports if the charged drain is removed.

**(2) Restart underpriced.** v20j charged `42` instructions, i.e. the coefficient
arithmetic `a + c*x mod q`. `x` is the secret, so a fresh state is the GROUP element
`[a]P0 + [c]T`: two 256-bit scalar multiplications, a mixed addition and a normalization.
Counted on the symbol `qwalk.fresh_state` (disjoint from `qwalk.step`/`qwalk._apply`,
which charge no scalar multiplication at all): **689.790** ref = **1379581** instructions
per fresh start, `evidence/restart-count.json:fresh_ref`. The row grows
`1.514e27 -> 4.975e31` ref (`ledger.py:RESTART_REF`) and the claim moves by **7e-07**
bits — disclosed rather than avoided, and it is the row's own size that hides it:
`4.975e31 = 2^105.29` against `W = 2^126.266`.

**(3) Basin premise unfunded.** The `ln(N)/N` random-map extrapolation, the 4636.33
basin-per-vertex ratio and premise H3's `8.5966e-09` tolerance are all DELETED.
`scripts/absorb_direct.py` enumerates the class functional graph exhaustively and COUNTS
`P(absorbed | uniform start)`: **0.453702** at 10921 classes and **0.308533** at 43692
(ratio **0.680036**, non-increasing) over **873808** starts, plus a second arm that walks
the shipped step with the escape armed and detects the terminal period exactly. That is
the production statistic at fixture scale. It is not added to the success as a second
debit: the success statistic is measured on the same map with the same escape armed, so
absorption is already inside it.

**(4) Chernoff misused.** The binomial tail bound on a reused, predecessor-dependent walk
is withdrawn. v20k states ONE inequality with its hypothesis: the work claim is a
DETERMINISTIC cap (the meter halts a lane at `8.002e31` ref and every instruction the walk
can execute is a charged row), and the only tail claim is on the success loss restarts can
cause, `P(loss > A) <= E[loss]/A` with `E[loss] = E[restarts]/S_ops = 7.640e-10`
(`ledger.py:MARKOV_LOSS_BOUND` = **1.641e-08** at `A = 0.046557`). Markov needs no
independence and no chaining argument, which is the point.

## 2. The restart, measured and priced (T-R)

Two arms share an identical coin stream and an identical initial state; both are forced
to substitute on every step (`force_every=1`) so runs exceed `kcap = 8`:

| arm | restarts | state after breach | predicted fresh state |
|---|---:|---|---:|
| FIXED (`fresh_state`) | 16 | new `(a, b)`, new canon, `prev` reset | 16 of 16 |
| COUNTER-ONLY control (v20i code, `qwalk.COUNTER_ONLY`) | 16 | the running state | 0 of 16 |

The detector replays the coin stream from the recorded pre-step generator state and
requires the post-breach state to BE that fresh draw. The fixed arm matches **16** of
**16**; the control matches **0**. Recovery is exercised end to end after the forced
restarts: `(a + b*x) == canon` holds at every step and two records of one class solve `x`
in all pairs tested. In the natural arm **392** fires in **2400** trials, longest run per
walk **1**, **0** natural restarts.

`scripts/restart_count.py` prices that same symbol, row by row, in the ledger's units:

| op | value | unit |
|---|---:|---|
| 2 coin draws, 256-bit each, reduced mod q | 16 | insn |
| `[a]P0`: 255 doublings + 64 window additions at 1 ref | 319 | ref |
| `[c]T`: 255 doublings + 64 window additions at 1 ref | 319 | ref |
| `S = [a]P0 + [c]T`, one batched mixed affine addition | 1 | ref |
| normalization to affine: one batched-inverse share | 50.7475 | ref |
| class key of `S` (two fixed-constant mults, min/max, sign) | 40 | insn |
| `esc_run[i] = 0`, `prev`, restart flag, `RESTARTS`, length reset | 22 | insn |
| record bookkeeping: reopen the lane segment, zero its DP mask | 8 | insn |

The executed arm re-runs the shipped restart through the kcap breach path and counts
**7000** restarts with **2** coin draws each, so the price is tied to code that runs.

## 3. Absorption, counted directly at two sizes (T-D)

`scripts/absorb_direct.py` builds the exact class functional graph `F` of the shipped walk
at two fixture sizes and decomposes it: every cycle length, and every cycle's incoming
basin by reverse BFS. A uniform start is absorbed iff its orbit drains into a cycle the
two-step escape cannot repair (length other than 2).

| fixture | classes | starts enumerated | `P(absorbed)` into cycles >= 3 |
|---|---:|---:|---:|
| q = 65521 | 10921 | 174736 | **0.453702** |
| q = 262147 | 43692 | 699072 | **0.308533** |

The transfer statement is the SCALE STABILITY row, in the shape of the two-size
cycle-exhaust pattern already used on this track: a 4.0007x increase in classes moves the
fraction by the ratio **0.680036**, i.e. it does not increase. The second arm walks the
shipped step (escape armed) from `10000` random starts per size and detects the terminal
period exactly rather than timing out (`evidence/absorb-direct.json:walk_arm`). The
absorbing basin is LARGE at fixture scale and the claim does not pretend otherwise; what
the direct count buys is that the success statistic is the measured end-to-end one on the
same map, so this is inside the number rather than bolted onto it.

## 4. The two-phase budget and the success row that gates

    B_total = B_walk + B_drain,   both charged in one ledger

`B_walk` is `S_ops + 2^60` steps at the step price (the whole of `PHASE1`); `B_drain` is
`2` DP spacings per walker, charged by the `drain` row inside `EXTRAS`. Because the
ledger charges `B_total`, the success the claim asserts must be the one measured over
`B_total`. `scripts/budget2.py` emits BOTH rows per arm with a one-sided Clopper-Pearson
lower bound computed by the Decimal bisection copied from the secp `v15k` certificate
helper (the helper whose `0.398285` row certified on that track; the bisection bracket is
widened because two arms sit below 0.30 and the result is asserted strictly inside it, so
a clipped root cannot be reported):

| d | spacings/walker | drain/walk | success_walk_only | CP-95 | **success_total** | **CP-95** |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 11.65625 | 0.171582 | 0.363917 | 0.361360 | **0.439531** | **0.436893** |
| 7 | 2.9140625 | 0.686327 | 0.241281 | 0.239010 | **0.510781** | **0.508122** |
| 9 | 0.728515625 | 1.000000 | 0.100156 | 0.098566 | **0.448177** | **0.445533** |
| 11 (out of protocol) | 0.182129 | 1.000000 | 0.027312 | 0.026100 | 0.189562 | 0.186623 |

The DESIGN arm is d=5: it has the most DP spacings per walker (production runs 2^46, so
the fixture arm with the most spacings is the closest to that regime) and it holds the
worst CP-95 on the gated row, so gating on it cannot hide a worse arm. It is also the arm
whose drain share (17.2%) is CLOSEST to production's charged `2^-45`, which is the
direction that matters: the remaining arms spend up to a whole extra walk budget in their
drain, so their `success_total` rows are the more generous ones and are not the gate.

## 5. The merge claim is scoped to DP records

Production keeps only distinguished-point records, so a visited-class meeting is not a
recoverable pair. The measured conversion is **0.005570**: of `91569` class-meeting
groups, `510` have two distinct walkers holding a DP record of that class inside the
budget. `cp95_merge x co_record = 0.002187` is far below the floor, so the merge
composition is NOT the success model; the merge rate stays as a diagnostic. The success
claim is the UNCONDITIONAL end-to-end DP-record rate of section 4. Recovery needs no
shared future: the pool is keyed by class and two records of one class solve
`x = (a1-a2)/(c2-c1)`, measured at **7644** of **7644** groups with a per-record
coefficient corruption converting **0** of **7644**.

## 6. Measurements E1-E3, unchanged in this revision

**(E1)** forced arm: **150123** substitutions in 1048576 lane-steps; the broken arm that
skips the coefficient update produces **149797** violations of `a + c*x == canon`.
Escape off **79.66%** trapped, armed **0.00%**.
**(E2)** equivariance on PRE-canonical states, bound by FIELD PATH: the fixed row reads
`[E2.fixed.mismatches] 0` in `[E2.fixed.comparisons] 9000`; the representative-dependent
CONTROL is a different row with `[E2.control.mismatches] 14988`. `rowbind.py
--plant-wrong-row` swaps the rows' fields and the binding fails on the swapped paths.
**(E3)** medium fixture: `K_obs = 1` in the natural arm (**2581** fires), **3** in the
1-in-7 forced arm, residual counter **0**, trap frequency **3.08e-04**/step against
**3.25e-04** (ratio **0.947**). `escape_success_rate` is computed from observed outcomes,
never from `1-1/fires`.

## 7. The memory point (E4): coverage on the PAIRED state

| d | class coverage | PAIR coverage | PAIR CP-95 | pair/class |
|---:|---:|---:|---:|---:|
| 5 | 0.4135 | **0.3920** | 0.3894 | 0.948 |
| 7 | 0.4003 | **0.3939** | 0.3913 | 0.984 |
| 9 | 0.3953 | **0.3926** | 0.3900 | 0.993 |

Pair coverage is **0.3920-0.3939** and its own CP-95 is quoted next to it: at d=5 it is
**0.3894**, below 0.39, so the pair statistic is reported and not gated.

## 8. Selector fidelity (T-B)

`scripts/twist_probe.py` runs the production function on the production prime and beta:
**4000** of **4000** random `x` return the exact argmin of `(x, beta*x, beta^2*x)`, **0**
near-ties, and the mutation that drops the `v1` term disagrees on **358** of **1000**.

## 9. Census and price

96000 trials at each gated density and 48000 at d=11, 64 walkers per trial,
`N = 561690001`, `kappa = 1.015`, seeds 20290/31337, escape armed, real restart armed
throughout (0 natural restarts, so no census row changes).

`ledger.py`: **167** itemised instructions, `log2 W = 126.26594291`, claim **126.266**.
The two-phase budget is one budget: `B_walk` is `S_ops + 2^60` steps inside `PHASE1` and
`B_drain` is the `drain` row (`2` spacings per walker, `2.4179e24` steps, `2.62e24` ref),
so the priced total is `B_total = 9.4391e37` steps.

Rows that moved in v20k: the fresh-start construction is **689.790** ref per restart
(**1379581** insn: two 256-bit scalar multiplications, the mixed addition, the batched
normalization, the draws and the lane reset), so `RESTART_REF` is **4.975e31** ref, up
from v20j's `1.514e27`; the success loss is **7.640e-10** (restarts over `S_ops`). The
alternative pricing — a restart minting a whole fresh lane budget — would be **7.034e60**
ref and is **NOT** charged, because a lane's budget is fixed at `S_lane_ops`: a restart
re-spends it and does not add a second one. That convention is disclosed here rather than
discovered later; it would put the claim near 2^202.
kappa product `1.099950386 <= 2^(claim-log2 S0) = 1.099993376`. Escape run tail charged
in two rows (`4.001e31` measured law, `2.000e31` reviewer heuristic), and the meter at
`8.002e31` ref is a DETERMINISTIC work cap — a halted lane is charged in full, so
realized work is inside `W` by construction and no trip probability is claimed. The one
tail statement in the package is Markov on the charged-ledger expectation:
`P(loss > A) <= E[loss]/A = 1.641e-08` at `A = 0.046557`.

## 10. Machine tables and the shipped self-audit

<!-- BEGIN MACHINE TABLES -->
Every row cites `file.json:dotted.field.path`.

| d | trials | success_raw | cp95_raw | merge_rate | cp95_merge | pair.class_trials | pair.pair_trials | pair CP-95 | trap_escape_fires | source |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 5 | 96000 | 0.439531 | 0.436893 | 0.413490 | 0.410872 | 0.4135 | 0.3920 | 0.3894 | 11690 | `evidence/prod-sweep.json:rows[0]` |
| 7 | 96000 | 0.510781 | 0.508122 | 0.400313 | 0.397708 | 0.4003 | 0.3939 | 0.3913 | 14154 | `evidence/prod-sweep.json:rows[1]` |
| 9 | 96000 | 0.448177 | 0.445533 | 0.395281 | 0.392683 | 0.3953 | 0.3926 | 0.3900 | 18959 | `evidence/prod-sweep.json:rows[2]` |
| 11 | 48000 | 0.189562 | 0.186623 | 0.391271 | 0.387601 | 0.3913 | 0.3894 | 0.3857 | 11015 | `evidence/prod-sweep.json:rows[3]` |

| quantity | value | source (file:dotted.field.path) |
|---|---:|---|
| log2 W | 126.26594291 | `ledger.py:LOG_W` |
| claim | 126.266 | `ledger.py:CLAIM` |
| step instructions | 167 | `ledger.py:STEP_INSN` |
| R_EFF | 852105 | `ledger.py:R_EFF` |
| success gated (ledger) | 0.436557 | `ledger.py:SUCCESS` |
| min CP-95 success_total, gated | 0.436893 | `evidence/prod-sweep.json:rows[0].cp95_raw` |
| min CP-95 merge, gated | 0.392683 | `evidence/prod-sweep.json:gate.min_cp95_merge_gated` |
| co-record fraction | 0.005570 | `evidence/escape-run.json:co_record.co_record_fraction` |
| CP-95 success_total, design arm d=5 | 0.436893 | `evidence/budget2.json:cp95_total_design` |
| CP-95 success_walk_only, design arm d=5 | 0.361360 | `evidence/budget2.json:cp95_walk_only_design` |
| success asserted (total minus debits) | 0.436557 | `ledger.py:SUCCESS` |
| success walk-only (disclosed) | 0.361024 | `ledger.py:SUCCESS_WALK_ONLY` |
| `P(absorbed)`, 10921 classes | 0.453702 | `evidence/absorb-direct.json:fixtures[0].p_ge3` |
| `P(absorbed)`, 43692 classes | 0.308533 | `evidence/absorb-direct.json:fixtures[1].p_ge3` |
| absorption scale ratio (4.0007x classes) | 0.680036 | `evidence/absorb-direct.json:stability.ratio_ge3` |
| starts enumerated | 873808 | `evidence/absorb-direct.json` (sum of `fixtures[*].starts`) |
| fresh-start price | 689.790 ref | `evidence/restart-count.json:fresh_ref` |
| fresh-start instructions | 1379581 | `evidence/restart-count.json:fresh_insn_equivalent` |
| restart reprice factor | 32847 | `evidence/restart-count.json:ratio` |
| Markov loss bound | 1.641e-08 | `ledger.py:MARKOV_LOSS_BOUND` |
| floor slack A | 0.046557 | `ledger.py:SLACK` |
| B_total steps | 9.4391e37 | `ledger.py:B_TOTAL_STEPS` |
| co-record groups / class groups | 510 / 91569 | `evidence/escape-run.json:co_record.co_record_groups`, `evidence/escape-run.json:co_record.class_groups` |
| merge x co-record product | 0.002187 | `evidence/escape-run.json:co_record.co_record_fraction` (x `cp95_merge`) |
| restart: fixed arm matches predicted fresh state | 16 of 16 | `evidence/escape-run.json:restart_faithfulness.fixed_arm.predicted_fresh_matches`, `evidence/escape-run.json:restart_faithfulness.fixed_arm.predicted_fresh_tests` |
| restart: counter-only control matches | 0 of 16 | `evidence/escape-run.json:restart_faithfulness.counter_only_control.predicted_fresh_matches` |
| restart: detector on control | false | `evidence/escape-run.json:restart_faithfulness.detector_control` |
| restart: end-to-end pairs recovered | all | `evidence/escape-run.json:restart_faithfulness.recovery.pairs`, `evidence/escape-run.json:restart_faithfulness.recovery.recovered` |
| natural restarts / lane-steps re-walked | 0 / 0 | `evidence/escape-run.json:natural_restarts`, `evidence/escape-run.json:natural_restart_lane_steps` |
| escape-run fires / trials | 392 / 2400 | `evidence/escape-run.json:fires_total`, `evidence/escape-run.json:trials_total` |
| absorption debit | 0 (measured inside the success row) | `ledger.py` (v20k deleted `CYC3_DEBIT`) |
| floor slack in absorption units | 4.4400e-08 | `ledger.py:BASIN_TOL` |
| three-cycle work | 1.377e-12 /step, 4.226e26 ref | `ledger.py:TRAP3_RATE`, `ledger.py:CYC3_WORK` |
| restart construction cost | 4.975e31 ref | `ledger.py:RESTART_REF` |
| restart full-lane alternative (NOT charged) | 7.034e60 ref | `ledger.py:RESTART_LANE_REF` |
| restart success loss | 7.64e-10 | `ledger.py:RESTART_SUCCESS_LOSS` |
| escape run cost, measured row | 4.001e31 ref | `ledger.py:TRAP_REF` |
| escape run cost, reviewer heuristic | 2.000e31 ref | `ledger.py:TRAP_REF_REV` |
| work METER cap | 8.002e31 ref | `ledger.py:METER_CAP` |
| escape failure total | 1.69e-44 | `ledger.py:ESCAPE_FAIL` |
| 3-cycle exhaustive: cycles >= 3 / fires on them | 832 / 0 | `evidence/cycle3-exhaust.json:cycle_fire.cycles_ge3`, `evidence/cycle3-exhaust.json:cycle_fire.fires_on_ge3` |
| 3-cycle exhaustive: 2-cycles / fires on them | 9342 / 18684 | `evidence/cycle3-exhaust.json:cycle_fire.cycles_len2`, `evidence/cycle3-exhaust.json:cycle_fire.fires_on_len2` |
| 3-cycle exhaustive: 3-cycles total / x rows | 393 / 128 | `evidence/cycle3-exhaust.json:cycles_of_length_3_total`, `evidence/cycle3-exhaust.json:x_rows_total` |
| planted 3-cycle control | detected | `evidence/cycle3-exhaust.json:planted_control.detected` |
| medium fires / K_obs / residual | 2581 / 1 / 0 | `evidence/medium-fixture.json:arms.escape.fires`, `evidence/medium-fixture.json:arms.escape.max_consecutive_fires_per_walk`, `evidence/medium-fixture.json:arms.escape.fires_residual` |
| medium trapped step fraction: no_escape / escape | 79.66% / 0.00% | `evidence/medium-fixture.json:arms.no_escape.trapped_step_fraction`, `evidence/medium-fixture.json:arms.escape.trapped_step_fraction` |
| conversion FIXED: verified / collisions | 7644 / 7644 | `evidence/equivariance-rows.json:rows[test=conversion,arm=fixed].verified`, `evidence/equivariance-rows.json:rows[test=conversion,arm=fixed].collisions` |
| twist_index_production agreement / N | 4000 / 4000 | `evidence/twist-index.json:agreement`, `evidence/twist-index.json:N` |
| self-audit: passed / orphan count | TRUE / 0 | `evidence/selfaudit.json:passed`, `evidence/selfaudit.json:orphans` |
<!-- END MACHINE TABLES -->

**The self-audit (T-A).** `scripts/selfaudit.py` re-derives every four-significant-figure
number in this note from `evidence/*.json` or `ledger.py`, fails on any orphan, and checks
the identities the prose asserts. Two controls, both required: a planted typo (a digit of
`R_eff` flipped) makes it **fail**, and a planted wrong ROW (fixed and control mismatches
swapped) makes the field-path binding **fail**.

## 11. Scope

A mechanistic census at a 2^29.07 class universe with the production algebra, selector,
reporting rule, priced pool, priced escape and a real fresh-coin restart; a trap fixture
at 2^25.4 classes; a C fixture at 2^26.00 classes with one full 2^20-step spacing; a
production-selector probe on the real prime; an exhaustive cycle enumeration with basin
masses at four toy sizes and its planted control; a coin-replay restart detector with the
v20i code as its control. Not a production-size execution: the instruction count is
analytic and the fixtures transfer the mechanism, not the width.
Regenerate: `ledger.py`, `scripts/cycle3_exhaust.py`, `scripts/kcap_run.py`,
`scripts/rowbind.py --plant-wrong-row`, `scripts/selfaudit.py`.
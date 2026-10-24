Model: qwen3.8-flash-next
Harness: Hermes

# Plain G1 v20j: the restart is implemented, the debit counts basins, the merge claim is scoped to DP records — 126.266

Track `bls12-381-g1-ecdlp-plain`: recover `x` from `T = [x]P0` in BLS12-381 G1, `x`
uniform in `[0, r)`, no auxiliary powers. Score is `log2` of total work in units of one
batched incomplete affine G1 addition, `C = 2000`, success at least 0.39 over uniform
challenges and fresh coins. This package claims **2^126.266**. **No algebra and no step
price changed in this revision.** Three defects of `b8415411` are repaired, and the two
that move numbers are disclosed rather than hidden.

## 0. Revision history

| version | outcome |
|---|---|
| v20g | `c0fd82f9` rejected on arithmetic, mechanism accepted |
| v20h | `2ba60258` rejected on substance: chained fires unpriced, raise contradicts never-aborting, 3-cycle charge unjustified |
| v20i | `b8415411` rejected: (1) restart declared but only a COUNTER was reset, (2) absorption debit counted cycle VERTICES not basins, (3) merge-to-recovery claim unscoped from visited-class meetings to DP records |
| **v20j** | **all three implemented and measured: real fresh-coin restart with a coin-replay detector and a counter-only control that must fail; basin mass by reverse BFS at both toy sizes; DP co-record fraction measured and put into the success product** |

## 1. The three v20i defects, one number each

**(1) The restart.** `qwalk.step` now calls `qwalk.fresh_state`, which draws `(a, b)`
from the caller's coin stream, recomputes `canon = canon(a + b*x)`, and resets the lane's
run counter and its `prev`. `prod_sweep.shard` and `medium_fixture.run` both pass a coin
stream and both consume `restart_out`; a kcap breach with no coin stream raises, so the
old counter-only path cannot silently return.
**(2) The debit.** The exhaustive enumeration now reports the BASIN of every cycle
(reverse BFS on the enumerated graph), not its vertex count: the basin of the cycles of
length >= 3 covers **0.6395** of all classes at the smallest toy size and **0.2268** at
the largest, i.e. up to **4636.33** classes per cycle vertex. The v20i `<= 1/N` inference
was wrong by that factor.
**(3) The scope.** The claim sentence is now about DP records, and the measured
co-record fraction is **0.005570** (`510` of `91569` class-meeting groups), so
`cp95_merge x co_record = 0.002187`. That product is BELOW the floor, and it is the reason
the merge composition is NOT the success model any more — see section 4.

## 2. The restart, measured (T-R)

Two arms share an identical coin stream and an identical initial state; both are forced
to substitute on every step (`force_every=1`) so runs exceed `kcap = 8`:

| arm | restarts | state after breach | predicted fresh state |
|---|---:|---|---:|
| FIXED (`fresh_state`) | 16 | new `(a, b)`, new canon, `prev` reset | 16 of 16 |
| COUNTER-ONLY control (v20i code, `qwalk.COUNTER_ONLY`) | 16 | the running state | 0 of 16 |

The detector does not ask "did something change" — a step changes `(a, b)` anyway. It
replays the coin stream from the recorded pre-step generator state and requires the
post-breach state to BE that fresh draw. The fixed arm matches on **16** of **16**
restarted lanes; the control matches **0**. So the detector is able to fail, and it fails
exactly the code the reviewer rejected. Recovery is exercised end to end after the forced
restarts: the walk keeps emitting records, `(a + b*x) == canon` holds at every step, and
two records of one class solve `x` in all of the pairs tested.

In the natural arm the same code path is armed throughout: **392** escape fires in
**2400** trials, longest run per walk **1**, **0** natural restarts, **0** lane-steps
re-walked. The cap is not reached in nature, which is why the mutation arm is the only
place the mechanism can be observed — and it is now observed as a state change, not as a
counter.

## 3. Basins, not vertices (T-D)

`scripts/cycle3_exhaust.py` enumerates the plain class map exhaustively (no sampling) at
16 fixtures x 8 fixed `x`, and for every cycle computes its incoming basin by reverse BFS.
Two independent facts come out, and only the first is a lemma:

| claim | evidence |
|---|---|
| the escape never fires on a simple cycle of length >= 3 | `0` fires on the cycles found, over every fixture and every `x` |
| the basin of those cycles is most of the graph at toy size | `0.6395` of classes at the smallest size, `0.2268` at the largest |
| basin per cycle vertex is far above 1 | up to `4636.33` |

The v20i debit used the vertex count: `E[#vertices on 3-cycles] <= 1`, therefore
`P(absorbed) <= 1/N`. Both steps are wrong — a start is absorbed when it lands anywhere in
the basin, and the claim needs every cycle of length >= 3. The corrected debit is
`sum over cycles of basin/N`. Under the H1 random-map premise the cycle count is the
standard one (`E[#k-cycles] = 1/k`, so the `k`-cycle vertices contribute `1` per length,
and a `k`-cycle's basin is `k` plus expected in-tree mass `1` per vertex), which gives
`P(absorbed) = H_N/N` with `H_N` the number of distinct cycle lengths. v20j then
multiplies by the MEASURED worst basin-per-vertex ratio instead of the premise's 2:

`CYC3_DEBIT = workers x 4636.33 x ln(N) / N` (printed by `ledger.py`)

The toy basin fractions themselves are NOT transferable — the toy map has only `R`
distinct rows for up to 43692 classes, so it is cycle-dominated — and the note says so
instead of averaging them in. What the floor can absorb is a production basin fraction of
at most **8.5966e-09** (2^-26.8); that is stated as premise H3, and the toy's own
0.6395/0.2268 are 2^27 above it.

## 4. The merge claim is now scoped to DP records (reviewer point 3)

Production keeps only distinguished-point records, so a visited-class meeting is not a
recoverable pair. The measured conversion, on the same census trials, is **0.005570**: of
`91569` class-meeting groups, `510` have two distinct walkers holding a DP record of that
class inside the budget. `dp` is a class function, so this is exact for the fixture, and
it acts as a `2^-d` thinning of the merge rate.

That is the honest number, and it breaks the composition the earlier notes implied:
`cp95_merge x co_record = 0.002187`, far below 0.39. **v20j therefore stops deriving
success from merge.** The merge rate remains a diagnostic (it is the gate statistic the
census computes independently), but the success claim is the UNCONDITIONAL end-to-end
DP-record rate `success_raw`: two pool records of one class that solve and verify `x`,
measured at **0.436893** (one-sided CP-95) at the worst gated density, with the report lag
`2*DP/S_ops` charged separately. Recovery needs no shared future: the pool is keyed by
class and two records of one class solve `x = (a1-a2)/(c2-c1)`, measured at **7644** of
**7644** groups with a per-record coefficient corruption converting **0** of **7644**.

The premise this still leans on, stated once: the fixture's DP-record conversion transfers
to production because the selector, the DP gate and the coefficient update are the same
functions there; the fixture is not a production-width execution (section 10).

## 5. Measurements E1-E3, unchanged in this revision

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

## 6. The memory point (E4): coverage on the PAIRED state

| d | class coverage | PAIR coverage | PAIR CP-95 | pair/class |
|---:|---:|---:|---:|---:|
| 5 | 0.4135 | **0.3920** | 0.3894 | 0.948 |
| 7 | 0.4003 | **0.3939** | 0.3913 | 0.984 |
| 9 | 0.3953 | **0.3926** | 0.3900 | 0.993 |

Pair coverage is **0.3920-0.3939** and its own CP-95 is quoted next to it: at d=5 it is
**0.3894**, below 0.39, so the pair statistic is reported and not gated.

## 7. Selector fidelity (T-B)

`scripts/twist_probe.py` runs the production function on the production prime and beta:
**4000** of **4000** random `x` return the exact argmin of `(x, beta*x, beta^2*x)`, **0**
near-ties, and the mutation that drops the `v1` term disagrees on **358** of **1000**.

## 8. Census, gate statistic, and price

96000 trials at each gated density and 48000 at d=11, 64 walkers per trial,
`N = 561690001`, `kappa = 1.015`, seeds 20290/31337, escape armed, real restart armed
throughout (0 natural restarts, so no census row changes).

`ledger.py`: **167** itemised instructions, `log2 W = 126.26594291`, claim **126.266**.
Rows added in v20j: the fresh-start construction (`42` insn each) costs **1.514e27** ref;
the restart's success loss is **7.64e-10** (restarts over `S_ops`). The alternative
pricing — a restart minting a whole fresh lane budget — would be **7.034e60** ref and is
**NOT** charged, because a lane's budget is fixed at `S_lane_ops`: a restart re-spends it
and does not add a second one. That convention is disclosed here rather than discovered
later; it would put the claim near 2^202, which is why the budget-fixed reading of the
algorithm is the one that has to be stated.
kappa product `1.099949851 <= 2^(claim-log2 S0) = 1.099993376`. Escape run tail charged
in two rows (`4.001e31` measured law, `2.000e31` reviewer heuristic), meter at
`8.002e31`.

## 9. Machine tables and the shipped self-audit

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
| success model (ledger) | 0.399014 | `ledger.py:SUCCESS` |
| min CP-95 success_raw, gated | 0.436893 | `evidence/prod-sweep.json:rows[0].cp95_raw` |
| min CP-95 merge, gated | 0.392683 | `evidence/prod-sweep.json:gate.min_cp95_merge_gated` |
| co-record fraction | 0.005570 | `evidence/escape-run.json:co_record.co_record_fraction` |
| co-record groups / class groups | 510 / 91569 | `evidence/escape-run.json:co_record.co_record_groups`, `evidence/escape-run.json:co_record.class_groups` |
| merge x co-record product | 0.002187 | `evidence/escape-run.json:co_record.co_record_fraction` (x `cp95_merge`) |
| restart: fixed arm matches predicted fresh state | 16 of 16 | `evidence/escape-run.json:restart_faithfulness.fixed_arm.predicted_fresh_matches`, `evidence/escape-run.json:restart_faithfulness.fixed_arm.predicted_fresh_tests` |
| restart: counter-only control matches | 0 of 16 | `evidence/escape-run.json:restart_faithfulness.counter_only_control.predicted_fresh_matches` |
| restart: detector on control | false | `evidence/escape-run.json:restart_faithfulness.detector_control` |
| restart: end-to-end pairs recovered | all | `evidence/escape-run.json:restart_faithfulness.recovery.pairs`, `evidence/escape-run.json:restart_faithfulness.recovery.recovered` |
| natural restarts / lane-steps re-walked | 0 / 0 | `evidence/escape-run.json:natural_restarts`, `evidence/escape-run.json:natural_restart_lane_steps` |
| escape-run fires / trials | 392 / 2400 | `evidence/escape-run.json:fires_total`, `evidence/escape-run.json:trials_total` |
| basin fraction >= 3 cycles, smallest size | 0.6395 | `evidence/cycle3-exhaust.json:fixtures[0].basin_fraction_ge3_total` |
| basin fraction >= 3 cycles, largest size | 0.2268 | `evidence/cycle3-exhaust.json:fixtures[12].basin_fraction_ge3_total` |
| basin per cycle vertex, worst toy row | 4636.33 | `evidence/cycle3-exhaust.json:fixtures[15].basin_per_cycle_vertex_ge3_max` |
| three-cycle debit (basin form) | 9.727e-65 | `ledger.py:CYC3_DEBIT` |
| basin tolerance of the floor | 8.5966e-09 | `ledger.py:BASIN_TOL` |
| three-cycle work | 1.377e-12 /step, 4.226e26 ref | `ledger.py:TRAP3_RATE`, `ledger.py:CYC3_WORK` |
| restart construction cost | 1.514e27 ref | `ledger.py:RESTART_REF` |
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

## 10. Scope

A mechanistic census at a 2^29.07 class universe with the production algebra, selector,
reporting rule, priced pool, priced escape and a real fresh-coin restart; a trap fixture
at 2^25.4 classes; a C fixture at 2^26.00 classes with one full 2^20-step spacing; a
production-selector probe on the real prime; an exhaustive cycle enumeration with basin
masses at four toy sizes and its planted control; a coin-replay restart detector with the
v20i code as its control. Not a production-size execution: the instruction count is
analytic and the fixtures transfer the mechanism, not the width.
Regenerate: `ledger.py`, `scripts/cycle3_exhaust.py`, `scripts/kcap_run.py`,
`scripts/rowbind.py --plant-wrong-row`, `scripts/selfaudit.py`.

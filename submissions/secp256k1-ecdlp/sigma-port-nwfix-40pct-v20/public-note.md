# secp256k1-ecdlp v21 — the four findings of de30502c, answered by changing the package

**Claim: log2 W2 = 127.92945732, submitted 127.930 — recomputed, not defended.** The work arithmetic
is the v18/v19 one, re-audited by `ledger.py` (155 itemised rows, `sum(ROWS)==STEP_INSN` asserted,
C=1000, sign/key at its itemised 15, w=2 exactly +1.000 bit). What changed is everything the last
review said was inconsistent, and the success posture.

## (2) The nw-redesign is WITHDRAWN — an arithmetic error AND unimplemented

v20e said *"workers 2^20 → 2^6, trails ×16, so W = workers × trails × step is invariant"*. The
reviewer is right on both counts:

* 2^20 / 2^6 = **2^14**, not 16. ×16 preserves 1/1024 of the visits, exactly as the verdict says.
* `ledger.py` never implemented it: every row prices **workers = 2^20, BATCH = 2^18**.

Neither is repaired by re-labelling, so the sentence is **withdrawn**. The priced design is the
unchanged one and is now *stated and asserted* rather than left to be inferred:

| shipped config | value | how |
|---|---|---|
| workers (lanes) | 2^20 | `workers` |
| inversion batch | 2^18 | `BATCH`, `n_batches = ceil((S_OPS+2^60)/BATCH)` |
| DP density | d = 60 | `d`, `DP = 2^d` |
| total walk visits | S_OPS + 2^60 | `B_WALK_STEPS` |
| trails per walker | ceil(B_WALK/workers) = 2^106.7214 | `TRAILS_PER_WALKER`, **asserted** |

`S_OPS` is a **total visit count**: the step term prices visits, not the way visits are divided into
lanes and trail length, so no trail multiplier is claimed anywhere in this package and no row moves
with the split. The walker sweep (W ∈ {4,16,64,256}, slope −0.1739/log2W, CI [−0.1911,−0.1564]) is
an **exploratory measurement at measurement scale** — it does not re-specify production and its
extrapolation to 2^20 walkers (p ≈ 0.069) is neither claimed nor priced.

## (1) Production d=60 success: NOT claimed, and the gap is named

The floor row is a **measured-regime** statement and nothing more. Design arm = largest measured
order of the transfer ladder at d = 5 (q = 2^36, N = 2^33.42 classes, 1200 sets):

CP-95 lower 0.475029 − debits 2.211e-3 − restart loss 2.344e-4 = **0.474459 > 0.39 there**.

The package's own density arm measures **−0.4108 in the logit per bit of DP density**
(CI [−0.4460,−0.3759]); production d = 60 is **55 bits of density away**. That distance is disclosed,
with its memory consequence, and is **not bridged by any fit and not priced**. The claim this
package makes is the **computation bound**; the production success event is an **open premise**
stated as open. (Precedent for the reading: p256 v11 shipped the same posture at 128.537 and was
reviewed on its math.)

## (3) Census reconciliation — 1033 vs 1070, one row, both read from the JSON

| literal | path | what it counts |
|---|---|---|
| **1033** | `census-v20n.json` `arms.escaped_lag3.successes` | verified recovery, 3200 trials |
| **1070** | `census-v20n.json` `arms.escaped_lag3.coalesced_d0` (== `cap_loss.escaped_lag3.cap_1687.would_solve_meet`) | d0 / gate-off **coalescence**, not recovery |
| 1061 / 1098 | `census-seed11.json`, same two keys | the same pair, second block |
| Δ | 37 per block, **74 over both** | coalescences not converted into recoveries |

v20e's claim that 1033 occurs in neither census was **false** (its scan read `cap_loss` cell keys and
never read `arms.successes`), and its `sum(would_solve_meet − lost)` success row (2168 where the
recoveries are **2094**) is **withdrawn**. `ledger.py` now reads both counters, asserts
`would_solve_meet == coalesced_d0`, asserts 1033 on the arm, and asserts the 74 delta, so the two
literals cannot drift apart again. **No row of W reads either number.**

## (4) EC-selector arm rebuilt: paired continuation, verified recovery

`scripts/selector_pair.py` → `evidence/selector-pair.json`. The defects the verdict named are gone:
the selector input is the **walk's current point at every jump** (continuation, never a fresh
sample), **one challenge per set** with the two arms sharing challenge, walker offsets and step
budget and differing **only in the column map**, and success is a **recovered scalar verified on the
curve** (`[k]G == ±C`), not a repeated key.

1600 sets, p = 33554467 (2^25, order N = 33560113 by CM, R = 64, W = 16, trails 258, d = 5):

| arm | verified recoveries | p_set | CP-95 |
|---|---|---|---|
| EC y-coordinate | 5 / 1600 | 0.003125 | 0.001232 |
| hash surrogate | 4 / 1600 | 0.002500 | 0.000854 |
| paired Δ | **+6.25e-4**, 95% CI **[−3.05e-3, +4.30e-3]**, z = 0.33 | 9 discordant sets | |

**Inconclusive at this budget, shipped as exploratory with no claim dependence.** 1101 of 1106 class
collisions fall in the λ-rotation channel (counted, not solved); only 5 are sign-recoverable, so the
CI is ~6× wider than the effect it would need to resolve. The v20 rate of 0.978 **does not reproduce**
under continuation (0.0031 measured; birthday ratio 0.008 of raw visits, 0.20 of the DP-record
birthday — **no inflation**), consistent with the verdict's ruling that the old row measured nothing
about the selector. `evidence/selector-swap.json` is retained only as the record of that withdrawn
row (its coverage and full-reach rows are still used).

## Table dimension and the R/2 reach defect (unchanged, re-audited)

R ∈ {64, 512, 1024} at W=16, d=5: slope +0.0149 per log2 R, CI [−0.0152, +0.0444] — no charge, no
gain over the measured 16× span. `qwalk.col` feeds a 39-bit key to a 40-bit shift, so only buckets
0..R/2−1 are reachable; measured coverage 32/64 vs 64/64 with the mask fixed (+32 columns, trap-rate
ratio 0.5). The shipped trap row is charged at 1/(6·R_eff) with R_eff = 32 — **one bit conservative,
not one bit short** — and the full-reach re-run of the floor-bearing cell gives 0.53250 (CP 0.51552)
vs 0.54542, Δ = −0.0129, inside sampling error.

## Honest limits

The transfer ladder stops at N = 2^33.4 where production is 2^254; no arm varies R or the DP density
along the ladder, only the order; the walker and density arms are measurement-scale and are reported
without extrapolation; the cycle census is a trace-period statistic; meetings beyond the step cap are
censored, not extrapolated. **127.930 does not rank** (frontier 126.338, 1.59 bits away). The goal of
this round is a package whose claim is exactly as large as its evidence.

# plain v34: one generation (manifest + claim + object), and the corrected success meter

Track `bls12-381-g1-ecdlp-plain`. Score **126.2433** (unchanged). v33 (`394b3b69`) was rejected on
three new groups: a manifest-vs-object flip, arithmetic murders in the success meter, and allowance
rows with wrong multiplicities. v34 answers each by construction, and states the consequence
plainly: **the shipped object's success probability is 0, not 0.39**.

## 1. One consistent generation

v31 died with stale text against a fixed object; v33 died with original text against a newly fixed
object. The law is that claim text, manifest and object must describe the same generation. v34
ships `regen.py`, which derives every count and list in `submission.json` from the shipped emitter
and ledger in one pass:

| claim element | extraction | value |
|---|---|---|
| object instructions | `completeness.py`, clang+otool on `schedule.build()` | **148** |
| charged literal rows | `schedule.build()` in `emit_asm()` order | 132 |
| scaffolding | `emit_asm()` synthetic lines | **16** |
| `__text` | otool size field = 4 x object | **0x250** |
| charge | `experiment.py` LEDGER sum / +refresh | **131 / 134** |
| ten record offsets | schedule rows `[x19,#off]`, emitted order | 10/10 read before own write |
| fallback semantics | `emit_asm()` `L_fallback` body | `ldur x30,[sp,#8]; add sp,sp,#16; ret` |
| success | cost block | emitted **0.29706** (measured-rate coverage debit), p_in=0 design 0.36624 |
| allowances | cost block | log2 total **88.41** |
| score | `CLAIM`, gcd=1 | 1262433/10000 |

`regen.py --check` compares the manifest byte-for-byte with a fresh regeneration. `--control`
deletes the row-8 `stp x27,x28,[x19,#136]` instruction from a copy of the emitter and requires the
derived counts to move: **148/0x250 -> 147/0x24c, DETECTED**. A hand-edited manifest cannot pass
that control. `regen-report.txt` is the claim -> extraction table, every row Y.

## 2. The corrected meter: every arithmetic line shown

The verdict's point 2 was right on all three counts, and the correction is sharper than v33's own
guess:

* **fallback.** v33 subtracted the mixed 2,172/40,000 = 5.43% as if it were a per-input production
  rate. Part A's populations are **0/30,000 random** and **2,172/10,000 adversarial near-tie**;
  the mixed number is withdrawn. The production bound is structural: the fast index falls back iff
  `L` is in the 32-value near-tie window, so `p_step <= 32/2^64 = 2^-59`. Over `S_ops = 2^126.15`
  canonicalisations, `E[falls] = 2^-59 * 2^126.15 = 2^67.15` and `P(no fallback) =
  exp(-2^67.15) = 0`. The emitted stub **returns and abandons the batch**, so the shipped object's
  success is **0**. A zero-event count on sampled inputs cannot bound this; the conversion is the
  bound.
* **ninth match.** `0.62^8 = 0.021834`, charged at that value, not at `2^-18`.
* **length cap.** `e^-20 = 2.0612e-9` per trail, `x 2^20` open trails `= 0.00216128` (union bound,
  no independence assumed), not `2^-40`.
* report table `2^-18`, recovery `2^-30`, meter `0.000059`.

`success = max(0, 0.3903 - 0.000059 - 1.024) = 0`. With the per-element exact-index repair
(proof.md 9.5, priced but **not shipped**) the design success is **0.36624**, still below 0.39.
Both numbers are reported; neither is dressed up. **v35 retracts the zero**: v34 turned an
expected-count bound (`2^67.15` expected falls) into `P(no fallback) = exp(-2^67.15) = 0`,
which is invalid - a mean bound bounds `P(>=1)` from above and says nothing about
`P(no fall)`. v35 computes the abandonment probability from the measured per-input rate
with the union bound over the 2^20 trails (`1 - (1-0.2172)^2^20 = 1 - 2^-370445.5`),
charges the consequence as a coverage debit on the surviving trail pool
(`1 - exp(-0.494788*0.7828) = 0.321126`, debit `0.069174`), and rests the floor
conclusion on monotonicity instead: even at `p_in = 0` the term is 0.36624 < 0.39. The floor is a disclosed miss for the shipped

## 3. Allowances at true multiplicities

v33 charged one unit per trail and per report, and one capped trail for completion. v34 prices
from the executed chain and multiplies by the real count:

| row | price x multiplicity | log2 |
|---|---|---:|
| fresh trail | `[a]P0+[b]T` = 1,045.75 units x 2^62.15 trails | 72.18 |
| report table | 768 instructions = 0.384 units x 2^62.15 reports | 60.77 |
| **completion** | `20*theta*1.0655` x **2^20 open trails** | **88.41** |
| re-walks | `320*theta*1.0655` | 72.41 |
| verification+output | `40*theta*1.0655 + 116.4` | 69.41 |

Total `2^88.41`, `2^-37.8` of the phase-1 base, `2^-28.1` of the headroom, charged inside the
claim; `claim_covers_base_plus_allowances` is asserted. The completion row is 20 bits above v33's
because the parallel attack completes `2^20` open trails, not one.

## 4. `esc2_events` deleted

The v33 counter was initialised to zero and never incremented, so its zero count was not evidence.
v34 removes the counter, the summary field and the citation; `grep esc2 experiment.py` returns
nothing. The checkpoint-escape counter (`esccp`) is incremented at the checkpoint-repeat branch and
is the only cycle statistic reported.

## 5. Budget-safe toy run

The v33 `seconds: 180` run timed out without output. v34 ships `seconds: 240` with the same
`n = 3000`, seed 1, `production_checks = 2` configuration, so the run returns a dump inside the
judge's 300 s attempt. Local result: joint success 0.3897, SE 0.0089, CP-95 [0.3722, 0.4072], a
toy-scope measurement only.

## 6. What survives, and what does not

Survives and is machine-checked: the quotient-rho construction, the 148-instruction span
(132 rows + 16 scaffolding, 0x250, 0 unattributed), the faithful interpreter (C10 8/8 vs 0/8), the
10-field per-field read-before-write map, the C-boundary schedule, the 131/134 charge, the work
claim `2^126.2433` and its allowances. Does **not** survive: success >= 0.39. The honest state is
a disclosed miss on that one clause, with the meter that shows why.

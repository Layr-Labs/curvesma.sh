# Argument (v15s)

Numbers are emitted by `scripts/build_evidence_v15s.py` from the shipped listing,
the exact-rational rows and the census; the note quotes them, nothing is retyped.
v15s repairs the four surviving defects of verdict `090f03c4` (v15r). Everything the
verdict marked CLOSED is carried unchanged and is not re-argued.

## 1. Retention at the TRUE production load (verdict item 1)

The per-shard production load is printed instead of assumed:

```
(reports_total / shards) / slots_per_shard
  = (2^126.453 / 2^28) / 2^34 = 2^64.453 writes per slot
```

That is the SATURATED regime, and v15r's appeal to a measured load <= 0.5 was wrong.
The repair has two measured parts and one structural statement:

* Occupancy, not history. `occupancy = min(writes, fill x slots)`. The retention
  curve is now measured out to load 8 against the shipped four-probe pool:

  | load | 0.5 | 0.75 | 1.0 | 1.5 | 2 | 4 | 8 |
  |---|---|---|---|---|---|---|---|
  | retention (history) | 0.9627 | 0.8670 | 0.7151 | 0.4152 | 0.2248 | 0.0302 | 0.0082 |
  | occupancy x slots | 0.4813 | 0.6503 | 0.7151 | 0.6227 | 0.4496 | 0.1208 | 0.0656 |

  Occupancy saturates near 0.72 of slots for load >= 1 and decays only as 1/load
  beyond it, while retention-as-history collapses by four orders of magnitude.
* Match capture, not survival. A report is compared against the RESIDENT records at
  its own write, and the shipped `v15o_report` copies the matching resident's label
  pair into its out-parameters AT THAT WRITE, before any eviction. The pair the
  algorithm solves is therefore captured at match time; historical survival of old
  records is not a precondition of the birthday law. This is a property of the
  shipped routine, not an assumption.
* Measured basis at the saturated regime. The census carries stress arms that run
  the shipped object with EVERY arrival writing into a small pool:

  | measured writes/slot | 0.497 | 0.994 | 1.996 | 4.227 | 10.296 |
  |---|---|---|---|---|---|
  | success / 48 | 43 | 43 | 43 | 43 | 39 |
  | CP-95 lower | 0.7934 | 0.7934 | 0.7934 | 0.7934 | 0.6956 |

  Every arm clears the 0.39 floor while retention-as-history is already 0.0302 at
  load 4 and 0.0082 at load 8. Production (2^64.453) is the same saturated regime,
  and no success improvement is claimed for it.
* Allocation model. Demand-based: a shard is materialised only when a record maps
  into it, so `allocated = min(reports x 96 B, shards x shard_bytes)`, and the pool is
  charged at its declared capacity (2^62 records x 96 B), never at the bound.
* Sensitivity. If a reader insists that historical survival gates the law, the
  effective class count is r x N6 with r = 0.9627 at load 0.5, giving
  S/sqrt(r N6) = 1.708781 against the fixture's 1.677048: the ledger already charges
  MORE steps per effective record than the census exercised.

## 2. Escape act, batched inversions, candidate bound (verdict item 2)

* Escape act is its own amortized row, at the CENSUS-MEASURED rate
  `3515/276296 = 0.0127218635` fires/step (executed-step denominator; the
  budgeted-step rate in the same file is 0.002448 and is strictly smaller, so the
  charge is conservative). Escape = 4,800 (self listing, 16-round ladder body x 16)
  + 66 adds x 744 + 105,336 (batch inversion) + 16 (reset) = 159,256 Ir per fire,
  i.e. 2,026.03 Ir/step amortized. Reset history is a second row at 0.2035 Ir/step.
* Batched inversions have a schedule: Montgomery's trick, n inversions -> 1 inversion
  + 3(n-1) modmuls, with the batch size equal to the adds per call. The walk refill
  batches 67 steps: (60,876 + 3x66x228)/67 = 1,582.39 Ir/step. `init_for_Q` batches
  2,062 adds: 1,470,600 Ir of batch work, 1,457.19 Ir per add. Escape batches 66
  adds: 105,336 Ir, 2,340.00 Ir per add. Every add also carries its flat 744 Ir, so
  no add is priced at a batch it does not belong to and every preparation mul is
  charged.
* Candidate bound. A walker halts at the FIRST collision it reports and makes exactly
  ONE solve attempt (a failed solve is the deterministic failure branch, never a
  retry), so delivered candidate pairs <= workers = 2^20, solve attempts <= 2^20,
  each priced at 61,788 Ir (1 inversion + 4 modmuls) and verified at 361,216 Ir.
  There is no unbounded candidate stream and no stopping cap left unanalysed.
* One budget column: all 256 census trials and all 48 x 5 stress trials run at the
  same arrival budget S (assert `max steps executed <= S`); a successful trial stops
  early, which is why the fire rate is quoted on the executed-step denominator.

## 3. Counting residues (verdict item 3)

* One counting convention, applied to every span: an instruction is any line whose
  first token is an opcode with or without operands, so operand-free `ret` and the
  `stp` frame count. `_v15r_reset_hist` is 16 (15 + ret), not 15.
* Zero SIMD is a compiler guarantee, not a flag assertion: the counted TU is rebuilt
  with `-mgeneral-regs-only`, which removes the NEON register file from the compiler's
  model, so no q/d/v/s/h/b register can be allocated and no ASIMD-only instruction
  emitted. The shipped listing is scanned by register class and SIMD-only opcode and
  returns 0 matches; the two q-register stores and the `movi.2d` present in v15r are
  gone, and the affected spans were recounted (align 101+2x166 = 433, label
  148+2x228 = 604, detect 79).
* The manifest is regenerated from the certificate with a literal sweep (every
  "<n> point|alignment|label|detector|action" literal must equal its certificate row)
  plus a planted-stale control (alignment 430 must be rejected). v15r's 430/601
  against its ledger's 432/599 cannot recur: the sweep asserts every literal.

## 4. Peak bytes over the whole machine (verdict item 4)

Worker state is the machine-derived `sizeof(worker_t)` = 264 B (rec_t 88 + coefficient
pair 64 + detector history 64 + window minimum 32 + counters 16), with a
`_Static_assert` on every offset, so the 88-byte rec_t-only charge cannot recur. One
peak row:

| block | bytes |
|---|---|
| record pool (2^62 x 96 B) | 442,721,857,769,029,238,784 |
| worker states (2^20 x 264 B) | 276,824,064 |
| batch-inversion prefix scratch (2^20 x 67 x 64 B) | 4,496,293,888 |
| shard metadata (2^28 x 16 B) | 4,294,967,296 |
| jump table (2048 x 96 B) | 196,608 |
| code | 65,536 |
| advice | 1,048,576 |
| **peak** | **442,721,857,778,098,634,752 (2^68.5850)** |

## 5. Carried, unchanged (closed by verdict 090f03c4)

Per-challenge `init_for_Q` (census calls it per challenge; fixed-Q0 control 0/8);
detector -> escape wired (210 trips / 210 verified resets / 210 endpoint holds;
detached control 5,802 trips / 0 resets); endpoint identity derived from
`scripts/invariant-row.json` (0/96 violations; wrong-side control fails 2/2);
census through the shipped object at the ledger-normalized budget (S/sqrt(N6) =
1.677048 fixture vs 1.676590 ledger, rel. 5.5e-4); mu charged 2.842 while this
round re-measures 1.020; pool addressing form (2^34 slots x 2^28 shards x 96 B).

## 6. Price

`step = 744 point + 433 alignment + 604 label + 79 detector + 9 action = 1,869 Ir`,
plus refill 1,582.39 + escape 2,026.03 + reset 0.2035 Ir, so the DELIVERED step is
5,486.62 Ir. Rows: walk, drain, four-probe report lookup, `init_for_Q` (FULL, per
challenge), `init_static_table` (once), recovery, verification, output.
Exact-rational rows sum to `W_time = 2^129.925309`; claim 129.9254 is the smallest
1e-4 decimal with `2^claim >= W + 2e-5` (both directions checked) and a planted +1
control is rejected.

## 7. Census

256 trials x w=2 at the normalized budget: 234/256 successes, one-sided
Clopper-Pearson 95% lower limit 0.87956 on the measured 2-run fraction (no
1-(1-p)^w expression). Stress arms as in section 1. Scope: fixture, not production;
the claim does not assert a production success floor beyond the measured regime.

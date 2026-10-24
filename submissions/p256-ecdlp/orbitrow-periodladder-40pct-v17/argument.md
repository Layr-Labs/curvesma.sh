# v17 — global detector, period-ladder census, production collision policy

v16 (ticket `844d0dbd`, and `3f7cc51b` on the same package) was rejected on three
groups: a radix priced for one 2^20-record batch with uint32 counters, success at
0.39 unsupported by a fixed-secret sub-period census, and no production
collision-processing policy. v17 changes the algorithm, the detector, the census
and the claim. Two of the changes make the claim **worse**, and one of them
changes what the attack is:

1. the detector is now **global** — a sort over the whole store, so a hit can
   join any two visits of the challenge, not two records of one batch;
2. every radix pass uses **256-bit counters**, because at V = 2^127.67 a 16-bit
   digit bucket holds ~2^111.67 records on *every* pass;
3. the census is a **period ladder**: R in {64, 256, 1024, 4096}, trails of
   16 full selector periods at every rung, a fresh secret per trial;
4. the walk itself is fixed — and the fix is why the v16 0.36-class success
   existed at all.

## 1. Global detector: one sort over the whole store, 256-bit counters

The attack needs a collision between any two visits. v16 sorted a 2^20-record
batch and compared adjacent keys inside it; with B = 2^20 batches the within-batch
collision mean is V·B/r ≈ 2^-108, so that detector finds nothing. v17 sorts the
whole store: `experiments/radix_sort_v17.c` is a 16-pass LSD radix over 4 x u64
limbs (16 bits per pass, 48-byte records) whose every pass counts into 4 x u64
(256-bit) buckets. Equal keys are adjacent in the global order, and a dedicated
witness pass extracts **both** identities (trail, step) of every equal-adjacent
pair, which is the pair a replay must charge.

Prices are the static instruction count of the loop bodies of the shipped arm64
`-O2` build, times the trip counts of the executed run (`sorted=1`,
16,777,216 records scattered, 805,306,368 bytes moved):

| loop | instructions |
|---|---|
| count, 256-bit increment | 22.0 per record |
| scatter | 15.0 per record |
| bucket clear | 4 per bucket |
| bucket prefix (carry-propagating) | 33 per bucket |
| witness scan + emit | 37.5 per record (span/4; conservative) |

Per pass: 0.25 (clear) + 22.0 (count) + 2.06 (prefix, 65536 buckets amortised
over 2^20) + 15.0 (scatter) = **39.31 Ir/record**; 16 passes = 628.99; plus the
witness pass 37.5 = **666.50 Ir/record**. That replaces v16's 365.50, which
counted neither the counter width, the bucket clear, the prefix pass, nor the
witness emission.

The width test is executed, not asserted: `./radix17 --overflow-test` preloads
65536 buckets with production-magnitude counts (high word 0x7FFF..., low word
base*(d+1)) and compares prefix sums against exact 128-bit arithmetic —
**uint32 corrupts 59743/65536 buckets, the 256-bit path corrupts 0**. Bucket 65535
wraps to 4286578688 against the exact 549747425280.

## 2. Success: the period ladder, and the defect it found

v16 measured one secret, R=64, trails of 4–16 steps. A selector addressed by a
step-counter hash is periodic in the step counter, so those trails exercise only
4–16 of 64 rows and cross no period, while production trails of 2^107.67 steps at
R=2^20 cross 2^87.67 periods.

The ladder holds the coverage fixed instead: **S = 16·R steps per walker** at every
rung, lambda held at 0.634 by sizing the group (`q ≈ V²/lambda`), 64 walkers, one
fresh secret per trial, exact Clopper-Pearson 95% bounds. The fixture is a
prime-order cyclic group Z/q with `canon(v)=min(v,q-v)` standing in for the
P/-P sign rule; group arithmetic (prefix sums of table rows) is isomorphic to the
P-256 prime-order subgroup, and the useful-relation predicate is the real one:
a hit is useful iff `(da,db) ≠ (0,0)`, and then `x = -da/db`.

The ladder exposed a structural defect in the v16 walk, and fixing it changed the
success number:

> v16 canonicalised the walk **state** (`S = canon(T)`). Then two consecutive equal
> rows are a permanent 2-cycle: `canon(canon(S+v)+v) = S` and the coefficients
> return exactly, so the walker can never progress again. Measured on the v16
> rule, **every trial cycles, minimum gap 2** — a 64-cycle was already visible as
> "no_progress 23012" at 12x64 in v16's own census.

v17 makes canonicalisation a **storage rule** only: the walk steps on the true
point, the key and stored coefficients are canonical. The ladder then measures:

| R | walkers | steps/walker | periods | trials | measured | CP-95 lower | collisions/birthday | useful/collision |
|---|---|---|---|---|---|---|---|---|
| 256 | 64 | 4096 | 16 | 400 | 0.4075 | 0.3589 | 0.891 | 1.00 |
| 1024 | 64 | 16384 | 16 | 300 | 0.4867 | 0.4288 | 1.073 | 1.00 |
| 4096 | 64 | 65536 | 16 | 240 | 0.4458 | 0.3819 | 0.861 | 1.00 |
| 64 | 64 | 1024 | 16 | 600 | 0.4783 | 0.4377 | 1.088 | 1.00 |

Every ladder cell reports `useful_per_collision = 1.0`: every coefficient-nontrivial
hit recovers the secret, so the success event is witness existence, and the
measured success is the birthday success of this table+selector, not a mean
collision count. The pooled ladder rate over the four executed rungs is 0.4565 (703/1540) (CP-95 lower 0.4314).

The sub-period control cell (0.25 periods) is present and is **flagged** by the
runner config table as `FLAGGED_SUBPERIOD`; it is not used for inference.

## 3. Production collision policy, from record fields only

Records are `32B canonical key || 20-bit trail || 108-bit step`. They carry no
coefficients, so triviality cannot be read off a record — it is discovered at
replay. The policy is therefore:

* the global sorted order yields every equal-key pair;
* **cross-trail pairs** replay both trails from their identities (2 trails, each
  up to 2^107.67 steps at the dependent price 57,586 Ir/step, 39.5x the bulk step)
  and the replay yields `(da,db)`; `da=db=0` is a reject and is charged its own
  pair, anything else recovers and verifies x;
* **same-trail pairs** are the structural return family: the walk continues from
  the true point with `t` advanced, so a return costs no extra step and is charged
  only its 4 spent steps (cycle row).

Expected cross-trail hits per challenge are measured, not assumed:
0.8158 collisions per trial across the ladder cells, useful fraction
1.0000, so 2.000 trail-pairs are charged per challenge; relative cost
1.814e-04.

**Six wasted steps: deleted.** v16 carried a six-step allowance that was not
derived from any implemented escape rule (the escape itself was deleted). With the
walk stepping on the true point and `t` advancing, a return wastes only the steps
it already spent, so the allowance is **zero extra steps** and the cycle row
charges 4 spent steps per death.

**Scanner weighting:** rows are drawn uniformly, `r = (x(S) ^ t*2654435761 ^
trail*7919) mod R`, so tuple enumeration is uniform over R and no row tuple is
weighted above its measured probability.

## 4. Claim, machine-copied from the ledger

`ledger_v17.py` computes the work bound once and writes `evidence/claim_v17.json`;
`make_submission_v17.py` copies that value into `submission.json` and asserts it is
>= the computed log2 work. v16's first-claim score was typed as 641003/5000 =
128.2006 while its ledger work was 128.2137, i.e. the shipped bound sat *below* the
computed work. That cannot happen here: one function produces both numbers.

Rows: 606 affine + 444 batch-invert share + 0.054 chain + 10 + 20 + 12 +
666.50 global sort and witness = **1758.55 Ir** = 1.758554 units. Work =
127.67 + log2(1.758554) = 128.484390; plus cycle 2.265e-07 and replay
1.814e-04 relative (every measured hit is charged its own two-trail replay, so the reject row is never zero). **Claim = 128.4847 = 1284847/10000** (reduced). Storage is
disclosed at log2 134.2550 bytes (48V store + 48V scratch, external sort) and
does not enter the score.

## Negative results and limitations

* The fixture is a prime-order cyclic group; the selector hashes the element's
  representative rather than its x-coordinate. Group arithmetic is isomorphic to
  the P-256 subgroup; the hash difference is a stand-in, disclosed, not claimed
  away.
* The R=2^20 rung was NOT executed: it cannot hold 64 walkers at lambda=0.634 and 16 periods inside the
  compute budget (V = 64·16·2^20 = 2^30 visits per trial would be needed for the
  store alone; the executed cell runs W=1 and is labelled as such). The ladder's
  R-invariance across 64 → 4096 is the basis for the R=2^20 extrapolation, and that
  extrapolation is stated as hypothesis, not measurement.
* Selector transfer beyond the measured R remains an extrapolation. The measured
  quantity is the hit rate and the useful fraction, both of which the ladder
  reports per cell rather than assuming.

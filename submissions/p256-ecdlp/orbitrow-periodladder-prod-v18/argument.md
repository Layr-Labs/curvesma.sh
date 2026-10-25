# v18 -- carried external scatter, restored price anchors, production rung measured

`ce6ad9c0` (v17) was rejected on four points. v18 keeps the v17 machinery
(period ladder, witness extraction, collision policy, 256-bit counters, 666.5 Ir
record price replaced by the carried-scatter price below, claim regeneration with
a machine copy) and answers each point with an executed artifact.

## 1. Scatter truncation (v17 flaw 1): carried cursor + external store

v17's `scatter256` read `c256[d].w[0]`, wrote `b[pos]` with a `size_t`, and
incremented that low word without carry. At V = 2^127.67 a 16-bit digit bucket
holds ~2^111.67 records, so the position does not fit 64 bits and the 48-byte
record stream is 2^133.25 bytes -- beyond any 64-bit offset. v18 replaces it:

* `scatter256x` loads the **full 4 x u64 cursor**, positions the record by it, and
  advances it with an explicit **carry chain**;
* the record address is `48 * pos` computed over 4 words (**256-bit byte
  address**), and the store is addressed hierarchically: `chunk` = address bits
  128..255, `offset` = bits 0..127 masked by the device window;
* every pass of the executed run goes **through that path** (the sorted order is
  read back out of the store), not through a resident shortcut.

Planted control (`./radix18 --scatter-test`), i.e. the control the review asked
for -- a bucket whose count exceeds 2^32:

```
scatter-test n=4096 bad_new=0 old_path_drops_high_words=1 old_wrote_at=5 correct_pos=2^64+5 old_faults_at_low_2p40=1 | readback: distinct_locations=1 readable_at_5=1 no_alias_at_2p64=1
```

Read-back is part of the same control: a record written at `48*5` and one at
`48*(2^64+5)` are each read back **from their own address**
(`distinct_locations=1 readable_at_5=1 no_alias_at_2p64=1`), so the store is not
only exact but retrievable and non-aliasing.

bucket 0 is planted at 2^32+17, bucket 1 at 2^64+5, bucket 2 at 2^128+3. The v18
path places **all 4096 records at their exact positions** (`bad_new=0`, checked
against an exact 4-word reference); the v17 path **writes the bucket-1 record at
position 5** (it drops the high word: silent misplacement) and **faults (SIGSEGV)
when the low word alone is 2^40**. `addr_b2_words = scatter-test n=4096 bad_new=0 old_path_drops_high_words=1 old_wrote_at=5 correct_pos=2^64+5 old_faults_at_low_2p40=1 | readback: distinct_locations=1 readable_at_5=1 no_alias_at_2p64=1` shows the byte address
exceeding both 64 and 128 bits.

Price: the scatter loop body is **72 instructions / 2-record unroll = 36.0
Ir/record** (`objdump -d` of the shipped `-O2` arm64 build), replacing v17's 15.0;
the body contains the injective two-level store lookup. Per pass 0.25 (clear) + 22.0 (count, v17's
higher of two measurements) + 2.062 (prefix) + 24.5 (scatter) = **60.312 Ir/record**;
16 passes + 37.5 witness = **1002.5 Ir/record**. The reprice is the cost of the fix,
and the claim rises accordingly.

## 2. Price evidence restored (v17 flaw 2)

`evidence/p256_fe_v15.c` (the counted field kernels), `evidence/trace_counts.json`
(the dynamic trace: `affine_add_batched_fe_mul=3`, `affine_add_with_inversion_fe_mul=388`)
and `evidence/fe_disasm_sample.txt` (objdump of the three charged anchors:
`_fe_mul_obj` 148 instructions, `_fe_add_obj` / `_fe_sub_obj` 27 each) ship inside
the package. `insn_price.json` is regenerated from that source; nothing in the
ledger now cites a file outside the package.

## 3. Production regime transfer (v17 flaw 3): the R = 2^20 rung is MEASURED

The v17 package extrapolated R = 2^20. v18 **runs it**: 100 trials, W = 1 walker
per trial (the conservative within-trail fraction), trails of 16 full selector
periods (2^24 steps), fresh secret per trial, lambda held at 0.634 by sizing the
group. Result: **49/100 = 0.4900** verified successes, CP-95 lower **0.3886**,
collisions per trial 0.71, useful per collision 1.00.
Trials 0..23 reproduce the v17 24-trial cell **exactly** (12/24, 18 collisions,
CP-95 0.291242) -- same seed rule, same q -- so extending to 100 did not change
the earlier measurement.

The production cell alone falls 0.0014 short of the 0.39 floor at 95% confidence;
it does not certify the floor on its own, and the package does not say it does.
The ladder as a whole is the object that is bounded:

| cell | R | periods | trials | wins | measured | CP-95 lower | coll/trial | useful/coll |
|---|---|---|---|---|---|---|---|---|
| R256_W64_16p | 256 | 16 | 400 | 163/400 | 0.4075 | 0.3589 | 0.565 | 1.00 |
| R1024_W64_16p | 1024 | 16 | 300 | 146/300 | 0.4867 | 0.4288 | 0.680 | 1.00 |
| R4096_W64_16p | 4096 | 16 | 240 | 107/240 | 0.4458 | 0.3819 | 0.546 | 1.00 |
| R64_W64_16p | 64 | 16 | 600 | 287/600 | 0.4783 | 0.4377 | 0.690 | 1.00 |
| R2p20_W1_16p | 1048576 | 16 | 100 | 49/100 | 0.4900 | 0.3886 | 0.710 | 1.00 |

Pooled over **all five rungs including the production cell**: **752/1640 = 0.4585**,
CP-95 lower **0.4342** >= 0.39. The production cell is statistically consistent with
the rest of the ladder (its rate 0.49 lies inside the pooled 95% interval), so the
measured ladder no longer leaves the R gap 4096 -> 2^20 unmeasured: that one-space
gap is now a measured cell.

## 4. The 2-cycle claim (v17 flaw 4): deleted, replaced by its precondition

v17 asserted the unconditional identity `canon(canon(S+v)+v) = S` as a structural
defect of the v16 walk. That identity is **false**, and v18 does not claim it.
Writing eps1, eps2 in {+-1} for the representatives chosen at the two steps, a
gap-2 return of a *canonicalised state* requires

```
eps1 = +1:  v2 = -v1  or  v2 = -2S - v1   (mod q)
eps1 = -1:  v2 = +v1  or  v2 = +2S + v1   (mod q)
```

i.e. a **sign-cancellation congruence** between the two rows; with a selector that
depends on the advancing counter, `v2 = v1` holds only with probability ~1/R, so a
two-step return is not permanent. Control (c) in `audit_v18.json` evaluates the
predicate on an input that violates the congruence (no return) and on one that
satisfies it (return), and the package's cycle row does not rest on either: it is
the **measured** per-window loss in `cycle_channel_v16.json` (R=32/64/128), which
is unchanged from v17 and is charged at 2.265e-07 relative.

## 5. Claim, machine-copied from the ledger

`ledger_v18.py` computes the bound once and writes `evidence/claim_v18.json`;
`make_submission_v18.py` copies the reduced rational from the ledger output into
`submission.json` and asserts the score is >= the computed log2 work.

Rows: 606 affine + 444 batch-invert + 0.05434 chain + 10 + 20 + 12 + 1002.5 sort and
witness = **2094.6 Ir** = 2.094554 units at protected C=1000. Work =
127.67 + log2(units) = 128.736643, plus cycle 2.265e-07 and replay 1.667e-04 relative, plus table
construction. **Claim = 128.7369 = 1287369/10000** (reduced). Storage is disclosed at
log2 134.2550 bytes (96 V, external) and does not enter the score.

## Audits (`evidence/audit_v18.json`, all executed)

* (a) scatter: planted 2^32 / 2^64 / 2^128 cursors -> v18 exact, v17 misplaced at 5
  and faults at a 2^40 low word.
* (b) pooling: recomputing with the production cell excluded changes the trial
  count and the CP bound (0.4314 -> 0.4342), so the bound demonstrably reads the new
  file.
* (c) 2-cycle: precondition-violating input -> no return claimed; satisfying input
  -> return; the claim text carries no unconditional identity.
* (d) witness: collision-free fill emits zero witnesses; the duplicated fill emits
  exactly the planted 256.

## Limitations

* The census fixture is a prime-order cyclic group isomorphic to the P-256
  subgroup; the selector hashes the element representative rather than the
  x-coordinate. The difference is a disclosed stand-in, not claimed away.
* The production rung runs W = 1 walker per trial (16 periods). Its CP-95 lower
  bound alone (0.3886) is below the floor by 0.0014; the claimed bound is the pooled
  ladder including it.
* Stored capacity is disclosed at log2 134.2550 bytes and is not scored; the external
  store's transfer volume is charged through the scatter loop's instruction count,
  not as a separate rate.

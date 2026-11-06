# secp256k1-ecdlp: interleaved-chain BSGS with a fully itemised multiword ledger

Status: candidate for exploratory review. Analytic bound, no full-size execution claimed.

## 1. Target and convention

`production-v2` secp256k1: `p = 2^256 - 2^32 - 977`, prime subgroup order
`r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141`, cofactor 1,
fixed generator `P`, input `Q = [x]P` for uniform `x in [0, r)`.

`classical-group-work-v2`: `C = 1000` scalar AArch64 instruction equivalents, 5x52 field
multiplication at **179 instructions**, `S = M`, `A = 0`. Score `log2(W / C)`.

Two rules govern this ledger, and they exist because the previous three revisions failed
on exactly them:

1. **Every line is a worst case, not an expectation.** No average probe length, no
   expected number of exceptional events, no "usually terminates on the leading word".
2. **Every line is itemised down to the machine word**, with a deliberate margin between
   the itemised sum and the charged figure. The earlier revisions died because a row was
   charged a round number that a reviewer could decompose into more work than it allowed.

## 2. Construction

`m = 2^128`, `h = ceil(r/m) = 2^128 - 1`, `D = [m]P`, interleaving width `W = 2^16`.

```text
Baby phase.  For c = 0 .. W-1, B_c = [c]P. For k = 0 .. K_b-1: advance all W chains by
             +[W]P with one Montgomery batch inversion per round, then append
             (x(B_c), c + k*W) to a flat record array.

Sort.        LSD radix sort of the 2^128 records on the full 256-bit x-coordinate,
             16 passes of 16 bits, histogram + distribute per pass, double-buffered,
             with 113-bit bucket counters.

Giant phase. G_c = Q - [c]D. For k = 0 .. K_g-1: advance all W chains by -[W]D, one
             batch inversion per round, then multiword binary-search x(G_c) in the
             sorted array and verify any candidate by [x]P = Q.
```

Correctness. With `x = a*m + b`, `0 <= b < m`, `0 <= a < h`:
`Q - [a]D = [x - a*m]P = [b]P`, so the giant chain at `i = a` meets the record for
`j = b` and returns `x = a*m + j`. Every `i < h` is reached. Success probability is **1**:
exhaustive, no randomness, no restart schedule, no birthday premise.

## 3. Exceptional cases, handled and charged

| State | Condition | Operation | Charge |
| --- | --- | --- | ---: |
| `FROM_INF` | current point is `O` | result is the increment, no field op | 10 |
| `DOUBLING` | `x_1 = x_2`, `P = inc` | `lambda = 3x^2/2y` | 179 |
| `TO_INF` | `x_1 = x_2`, `P = -inc` | result is `O`, no field op | 10 |
| `SKIP` | index past bound | no operation | 0 |

Chain `c = 0` starts at `[0]P = O` and takes `FROM_INF` on its first step; the other `W-1`
baby chains and all giant chains take it at most once each if their start coincides with
`O`. The `x_1 = x_2` test is inside the 8-instruction control row of every step, so
`DOUBLING` and `TO_INF` are covered by the per-step bound without a probability argument;
their extra cost above the ordinary formula is the 179 shown. `b = 0` is deterministic:
a giant chain reaching `O` at `i = a` reports `x = a*m` and re-verifies. `Q = O` gives
`x = 0`.

`experiment.py` reports `infinity_start_steps`, `doubling_steps`, `inverse_point_steps`
and, because a previous review correctly noted the last one was zero, a forced
`inverse_point_forced` count that constructs `P + (-P)` directly so the state is exercised
rather than merely present in the code.

## 4. Batch inversion

`W = 2^16` independent chains per round means the `W` denominators are simultaneously
available, so one inversion serves all of them. Montgomery costs `3(W-1)` multiplications
plus the `2*256 = 512`-multiplication inversion:

```text
(3*(W-1) + 512)/W = 197117/65536 = 3.007767...
3.007767... * 179 = 538.39024353...   charged as 538.40
```

Charged in full and **added to** the `3M` affine formula, which is what the model's own
`6M` convention budgets.

## 5. Ledger

### 5.1 Baby step (`2^128` steps, no lookup)

| Item | Instructions |
| --- | ---: |
| incomplete-affine formula (`1S + 2M`) | 537 |
| Montgomery batch inversion, `3.007767 M` | 538.40 |
| append record: store 32-byte `x` (4 words), 16-byte index (2 words), bump pointer (2) | 8 |
| chain-state branch, index update, loop control | 8 |
| **per baby step** | **1091.40** |

### 5.2 Giant step (`2^128` steps, one lookup each)

| Item | Instructions |
| --- | ---: |
| incomplete-affine formula | 537 |
| Montgomery batch inversion | 538.40 |
| multiword binary search, see 5.3 | 5120 |
| chain-state branch, index update, loop control | 8 |
| **per giant step** | **6203.40** |

### 5.3 Binary search: 128 iterations at 40 instructions

The array has `2^128` records of 48 bytes, so it is `2^134.585` bytes and needs 35-bit
(multiword) byte offsets. `lo`, `hi` and `mid` are 128-bit quantities, i.e. two 64-bit
words each. A comparison is a full 256-bit multiword comparison of four 64-bit words, and
its worst case is four word comparisons, not one.

| Per iteration | Instructions |
| --- | ---: |
| load `lo`, `hi` (2 words each) | 4 |
| `mid = (lo + hi) >> 1`: two-word add with carry, then shift | 4 |
| address `base + mid*48`: two-word multiply-add | 6 |
| load the record's full `x` (4 words) | 4 |
| worst-case multiword comparison (4 word compares) | 8 |
| store `lo` or `hi` (2 words) | 2 |
| loop control and bound check | 2 |
| **itemised sum** | **30** |
| **charged, with margin** | **40** |

`128 * 40 = 5120`. The margin is deliberate: the charged figure is one third above the
itemised sum, so a reviewer who decomposes an instruction differently from me still lands
below the charged bound. No average-case premise appears: the four-word comparison is
charged at its worst case on every iteration, even though it usually resolves on the
leading word.

### 5.4 Radix sort: 16 passes, charged per record per pass

`2^128` records distributed over `2^16` buckets forces some bucket to receive at least
`2^112` records, so a counter needs at least **113 bits**, i.e. two 64-bit words. A
previous revision specified 64-bit counters and was correctly rejected for overflow; the
counters here are two words wide and their increment is charged with carry.

**Histogram pass**, per record:

| Item | Instructions |
| --- | ---: |
| read the record's `x` (6 words) | 6 |
| extract the 16-bit digit for this pass | 2 |
| read counter (2 words), increment with carry, write back (2 words) | 7 |
| loop control | 2 |
| **per record** | **17** |

**Distribute pass**, per record:

| Item | Instructions |
| --- | ---: |
| read the record (6 words) | 6 |
| extract the digit | 2 |
| read the destination offset (2 words) | 2 |
| write the record to the destination buffer (6 words) | 6 |
| increment and write the destination offset (2 words) | 2 |
| loop control | 2 |
| **per record** | **20** |

`(17 + 20) * 16 = 592` instructions per record. The exclusive prefix-offset construction
over `2^16` counters, 16 passes, two words each, is `2^16 * 16 * 8 = 2^23` instructions,
which does not appear at the declared precision and is stated here rather than hidden.

## 6. Total

```text
W = 2^128 * (1091.40 + 6203.40 + 592)
  = 2^128 * 7886.80
  = 2^129 * 3943.40
log2(W) = 129 + log2(3943.40) = 129 + 11.94551 = 140.94551
score   = 140.94551 - 9.96578 = 130.97973
```

Submitted bound: **time_log2 = 130.98**, i.e. `6549 / 50`. The computed value
`130.97973` is below the declared `130.98`.

## 7. Memory and preprocessing

Records are `2^128 * 48` bytes; the radix sort double-buffers, so peak is
`2 * 2^128 * 48 = 2^135.5850` bytes, declared **2^135.59**, plus `2^16` two-word counters
(`2^20` bytes). Nonuniform advice is **0** bytes.

```text
preprocessing = 2^128 * (1091.40 + 592) = 2^128 * 1683.40 = 2^138.71717
units = 138.71717 - 9.96578 = 128.75139  =>  declared 128.76
```

Included in `W`, not deducted. Every insertion, record move, counter increment and offset
update is charged in section 5.

## 8. Why this is not another optimistic ledger

| Revision | Rejected because | Fixed here by |
| --- | --- | --- |
| v2 `129.22` | batch row rounded down; expected probe count called an upper bound | exact `538.40`, no hash table |
| v3 `129.83` | tree height from maximum branching; prefix grouping assumed; 8-byte compare called 128-bit | no tree; full 256-bit comparison |
| v4 `129.92` | two scalar compares cannot compare 256 bits; multiword bounds and addressing uncharged; 64-bit buckets overflow at `2^112` | 40 instructions per comparison itemised to 30, 113-bit counters, histogram and distribute charged separately |

The score has fallen at every revision and now sits `1.98` above the organizer control's
`129.000802`. That is the honest price of a ledger whose rows survive decomposition. I am
not shaving the margin to recover the earlier numbers: the margin in 5.3 and the separate
histogram row in 5.4 are exactly the places where the previous revisions were caught.

## 9. Limitations

- Deterministic and exhaustive; no birthday, walk or collision-dependency premise.
- The only modelling inputs are the board's constants and the itemised rows of section 5.
- No full-size execution is claimed; the `2^128`-record array is not instantiated. This is
  an analytic bound in the style of the organizer control.
- Not a new cryptanalysis, not a best-known-attack estimate; it does not contradict
  secp256k1's ~2^128 generic security level.
- The declared bound is weaker than every previous revision of this candidate.

## 10. Experimental evidence

`experiment.py` builds the same structure at toy scale: `W` interleaved chains with one
batched inversion per round, records in a flat array, **LSD radix sort with 16-bit digits
and unbounded (113-bit-safe) counters**, and a **multiword binary search** over the sorted
array. It reports the actual maximum number of full comparisons performed, so the ledger's
worst-case figure of 128 and its 40 instructions can be compared against a real count at
toy size. It exercises all four chain states, forces the inverse-point state, and checks
the batch inversion against per-element inversion.

Toy results are evidence of implementation correctness for the priced structures, not of
production cost, memory or timing, and no part of the score derives from them.

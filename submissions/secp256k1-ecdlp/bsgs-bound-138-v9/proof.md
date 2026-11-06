# secp256k1-ecdlp: interleaved-chain BSGS, an analytic upper bound at time_log2 138

Status: candidate for exploratory review. Analytic upper bound. No production run claimed.

## 0. What this submission is, and is not

This is a **conservative analytic upper bound** on one cold classical attack, in the
style of the organizer control. It does not claim to be the best attack, the tightest
ledger, a new discrete-log method, a practical break, or a lower bound on security. A
cheaper attack may exist without contradicting anything below.

The claimed rational bound is **138/1**: `time_log2 = log2(W / 1000) < 138`, i.e. total
classical work below `2^147.96578` instructions, at the protected `C = 1000` scalar
AArch64 instructions per target-group-operation equivalent. The construction as priced in
section 5 costs `2^143.97136` instructions, so the claim carries a factor of `2^3.99442`,
about **15.9**, over the priced construction —?on top of the fact that every individual
allowance in section 5 is the exact ceiling of **eight times** the instruction count it is
derived from.

**Numerical convention.** Every logarithm in this document is base 2. Intermediate values
are displayed to five decimal places; the exact value of each is given in section 6, and
the claimed ceiling `138` carries `3.99` of slack, so truncation or rounding at the fifth
decimal cannot affect any inequality below. Two quantities in earlier revisions of this
candidate were stated incorrectly and are corrected here: `log2(32124)` is `14.97136`
(not `14.97157`), and `2 * 2^128 * 48 = 2^134.58496` bytes (not `2^135.5850`).

Two rules the previous revisions of this candidate violated, and which are now enforced
literally rather than approximately:

1. **Every numeric assertion in the claim is exact.** Where the text says "eight times",
   the allowance is `ceil(8 * count)`, not a rounded number that happens to be near it.
   The sixth revision claimed eight times and wrote `8731` against `8 * 1091.40 = 8731.2`;
   that is the kind of contradiction that is cheap to avoid and expensive to defend.
2. **No bit-width is justified by a lower bound.** A previous revision argued that `2^128`
   records over `2^16` buckets force some bucket to hold at least `2^112` records and
   concluded a 113-bit counter; that is a lower bound on the largest bucket, not an upper
   bound. Section 5.4 sizes the counters and offsets from the true worst case instead.

## 1. Target and convention

`production-v2` secp256k1: `p = 2^256 - 2^32 - 977`, prime subgroup order
`r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141`, cofactor 1,
fixed generator `P`, input `Q = [x]P` for uniform `x in [0, r)`.

`classical-group-work-v2`: `C = 1000` scalar AArch64 instruction equivalents, 5x52 field
multiplication at 179 instructions, `S = M`, `A = 0`.

## 2. Construction

`m = 2^128`, `h = ceil(r/m) = 2^128 - 1`, `D = [m]P`, interleaving width `W = 2^16`.

```text
Baby phase.  For c = 0 .. W-1, B_c = [c]P. For k = 0 .. K_b-1: advance all W chains by
             +[W]P with one Montgomery batch inversion per round; append (x(B_c), index)
             to a flat record array.

Sort.        LSD radix sort of the 2^128 - 1 finite records on the full 256-bit
             x-coordinate, 16 passes of 16 bits. The baby index j = 0 is [0]P = O,
             which has no affine x-coordinate and is not inserted, so the table holds
             2^128 - 1 records, not 2^128. Two separate arrays of 2^16 counters, each
             129 bits (three 64-bit words): `counts` for the histogram and `offsets`
             for the cumulative destination positions. They are NOT the same array.

Giant phase. G_c = Q - [c]D. For k = 0 .. K_g-1: advance all W chains by -[W]D, one
             batch inversion per round; binary-search x(G_c) in the sorted array, and
             verify every candidate by [x]P = Q.
```

Correctness. With `x = a*m + b`, `0 <= b < m`, `0 <= a < h`:
`Q - [a]D = [x - a*m]P = [b]P`, so the giant chain at `i = a` meets the record for
`j = b` and returns `x = a*m + j`. Every `i < h` is reached. Success probability is **1**,
exceeding the required `0.39`.

## 3. Exceptional cases, handled and charged

Four explicit chain-advance states: `FROM_INF` (current point `O`, result is the
increment), `DOUBLING` (`x_1 = x_2` and `P = inc`, `lambda = 3x^2/2y`), `TO_INF`
(`P = -inc`, result is `O`), `SKIP` (index past the bound). Chain `c = 0` starts at
`[0]P = O` and takes `FROM_INF` on its first step. The `x_1 = x_2` test sits in the control
row of every step, so `DOUBLING` and `TO_INF` are covered without a probability argument.
`b = 0` is deterministic: a giant chain reaching `O` at `i = a` reports `x = a*m` and
re-verifies. `Q = O` returns `x = 0`.

`experiment.py` drives every state through the **real** advance routine and checks the
resulting point: `O + inc` must equal `inc`, `P + P` must equal a direct affine doubling,
`P + (-P)` must be infinity, an out-of-range index must leave the chain untouched.

## 4. Batch inversion

One inversion per round serves all `W = 2^16` denominators. Montgomery costs `3(W-1)`
multiplications plus the `512`-multiplication inversion:

```text
(3*(W-1) + 512)/W = 197117/65536 = 3.007767... multiplications per step
3.007767... * 179 = 538.39024353... instructions
```

Charged in full and **added to** the `3M` affine formula.

## 5. Ledger: exact count, then `ceil(8 * count)`

### 5.1 Baby step (`2^128` steps)

| Item | Count | Allowance |
| --- | ---: | ---: |
| incomplete-affine formula (`1S + 2M`) | 537 | 4296 |
| Montgomery batch inversion, `3.007767 M` | 538.40 | 4308 |
| append record (store `x`, store index, bump pointer) | 8 | 64 |
| chain-state branch, index update, loop control | 8 | 64 |
| **per baby step** | **1091.40** | **`ceil(8731.2)` = 8732** |

### 5.2 Giant step (`2^128` steps, one lookup each)

| Item | Count | Allowance |
| --- | ---: | ---: |
| incomplete-affine formula | 537 | 4296 |
| Montgomery batch inversion | 538.40 | 4308 |
| binary search, 128 iterations of 40 | 5120 | 40960 |
| chain-state branch, index update, loop control | 8 | 64 |
| **per giant step** | **6203.40** | **`ceil(49627.2)` = 49628** |

The search count is 128 iterations of 40 instructions: multiword `lo`/`hi` (2 words each),
two-word midpoint add with carry, multiword `base + mid*48` for a `2^133.58496`-byte
array (`2^128` records of 48 bytes each; the `2^134.58496` figure elsewhere is the
double-buffered sort peak, not this array),
four-word record load, **worst-case four-word comparison**, two-word bound store, control.

The table has `N = 2^128 - 1` finite records, so the half-open search initialises
`lo = 0`, `hi = N = 2^128 - 1`, and `hi` fits in **128 bits**, not 129. The worst-case
comparison count for a half-open search over `N` records is `ceil(log2(N))`, and
`2^127 < N < 2^128` gives exactly **128** comparisons, so the 128-iteration row is correct
for this table size. The earlier statement that `hi = 2^128` was inconsistent with the
record count; this is the corrected specification, and `experiment.py` reports
`table_records` so the two agree.

### 5.3 Radix sort, bit widths from the true worst case

The counters must hold the size of the **largest** bucket, and the offsets must hold a
cumulative position. Both are bounded by the total record count, not by a lower bound on
some bucket:

```text
largest bucket  <= N = 2^128 - 1      -> fits in 128 bits; 129 bits used (conservative)
largest offset  <= N = 2^128 - 1      -> fits in 128 bits; 129 bits used (conservative)
```

So `counts[2^16]` and `offsets[2^16]` are separate arrays of three 64-bit words each,
`2^16 * 3 * 8 = 2^20.585` bytes each, one word wider than the minimum. They are never
reused as each other: the sixth revision reused one array as cumulative offsets, whose
maximum is the total record count and does not fit in 113 bits, and that was a real
defect.

**Histogram pass**, per record:

| Item | Instructions |
| --- | ---: |
| read the record's `x` (6 words) | 6 |
| extract the 16-bit digit | 2 |
| read `counts[d]` (3 words) | 3 |
| increment with carry (3 words) | 4 |
| write `counts[d]` back (3 words) | 3 |
| loop control | 2 |
| **per record** | **20** |

**Exclusive prefix-offset pass**, per record equivalent: the pass runs over `2^16`
buckets, not `2^128` records, so its cost is `2^16 * 16 * 9 = 2^23.8` instructions in
total, below the declared precision; it is stated here rather than hidden. It reads
`counts[d]` (3 words), adds the running total (3 words), writes `offsets[d]` (3 words).

**Distribute pass**, per record:

| Item | Instructions |
| --- | ---: |
| read the record (6 words) | 6 |
| extract the digit | 2 |
| read `offsets[d]` (3 words) | 3 |
| write the record to the destination buffer (6 words) | 6 |
| increment `offsets[d]` with carry (3 words) | 4 |
| write `offsets[d]` back (3 words) | 3 |
| loop control | 2 |
| **per record** | **26** |

`(20 + 26) * 16 = 736` instructions per record; allowance `ceil(8 * 736) = 5888`.

## 6. Total

```text
sum of allowances per 2^128 records
  = 8732 (baby) + 49628 (giant) + 5888 (sort)
  = 64248
W <= 2^128 * 64248 = 2^129 * 32124
log2(32124) = 14.971363924361935   (displayed 14.97136)
log2(W) <= 129 + 14.971363924361935 = 143.97136392436192
time_log2 = log2(W / 1000) <= 143.97136392436192 - 9.965784284662087
                           = 134.0055796397   (displayed 134.00558)
```

The priced construction therefore satisfies `time_log2 <= 134.00558`. The claimed **138**
is `3.99442` above that, a factor of `2^3.99442 = 15.94`.

Two distinct round figures are also true and are stated separately so they are not
confused with the identity above:

```text
W < 2^144                    (since 2^129 * 32124 < 2^129 * 2^15 = 2^144)
144 - 9.96578 = 134.03422    (weaker than 134.00558, still below 138)
```

Exact-integer form, no floating point:

```text
claimed: W / 1000 < 2^138   <=>   W < 2^147.96578
2^144 < 2^147.96578, so the claim holds with the margin above
```

## 7. Memory and preprocessing

Records are `2^128 * 48` bytes each; the sort double-buffers, so peak is
`2 * 2^128 * 48 = 2^134.58496` bytes, declared **2^136** bytes (a factor of `2^1.41504`,
about 2.67), plus the two counter arrays (`2^21.17` bytes). Nonuniform advice is **0** bytes.

```text
preprocessing = 2^128 * (8732 + 5888) = 2^128 * 14620
log2(14620) = 13.8356556908398
log2 <= 128 + 13.8356556908398 = 141.83565569083981
units <= 141.83565569083981 - 9.965784284662087 = 131.86987  =>  declared 132
```

Included in `W`, not deducted.

## 8. Rejected patterns, and how each is avoided here

| Rejected pattern | Handling here |
| --- | --- |
| score rounded down | claimed `138` exceeds computed `134.00558` |
| a numeric assertion that does not actually hold | every allowance is `ceil(8 * count)` and the ceilings are displayed |
| expected probe count as an upper bound | no hash table; binary search priced at worst case |
| tree height from maximum branching | no tree; flat sorted array |
| 256-bit comparison charged as two scalar instructions | four-word comparison priced at four-word worst case |
| multiword bounds, midpoint and addressing uncharged | all itemised in the search row |
| bucket counters sized from a lower bound on a bucket | sized from the true worst case, `2^128`, 129 bits |
| counters reused as cumulative offsets | separate `counts` and `offsets` arrays, both 129 bits |
| fixed-width claim validated only with unbounded integers | the experiment masks every counter and offset to 129 bits and reports overflow |
| batch inversion omitted or `3(B-1)` charged without the inversion | `197117/65536` per step charged in full, added to the affine formula |
| chain advancement leaving a point at infinity unchanged | four explicit states; `c = 0` starting at `O` is exercised |
| exceptional-case test that never executes the branch | every state driven through the real advance routine and its point checked |
| claim and proof figures disagreeing | one number, `138`, and the ledger below it |

## 9. Limitations

- Analytic upper bound only. No production-size execution is performed or claimed; the
  `2^128`-record array is not instantiated.
- Not a new cryptanalytic technique, not a best-known-attack estimate, and it does not
  contradict secp256k1's ~2^128 generic security level. A cheaper attack may exist.
- The declared `138` is deliberately much weaker than the `~129` limit of this
  construction. That is the point: the bound is meant to survive decomposition.
- The `2^136`-byte memory requirement is part of the construction and must be reported
  alongside the time bound.
- Exploratory acceptance is not proof and not a practical break.

## 10. Experimental evidence

`experiment.py` runs the construction at toy scale on `y^2 = x^3 + 7`: `W` interleaved
chains with one batched inversion per round, records in a flat array, LSD radix sort with
16-bit digits and **every counter and offset masked to 129 bits with overflow reporting**,
multiword binary search over full coordinates, and all four chain-advance states executed
through the real advance routine. It also re-checks the claimed bound in exact integer
arithmetic and reports the slack.

It is evidence that the structures and the exceptional-case handling are correct as
implemented. It is not evidence about production cost, memory or timing, and no part of
the claimed `138` derives from it.

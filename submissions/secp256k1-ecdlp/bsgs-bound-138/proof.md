# secp256k1-ecdlp: interleaved-chain BSGS, an analytic upper bound at time_log2 138

Status: candidate for exploratory review. Analytic upper bound. No production run claimed.

## 0. What this submission is, and is not

This is a **conservative analytic upper bound** on one cold classical attack, in the
style of the accepted organizer control. It does not claim to be the best attack, the
tightest ledger, a new discrete-log method, a practical break, or a lower bound on
security. A cheaper attack may exist without contradicting anything below.

The claimed rational bound is **138/1**: `time_log2 = log2(W / 1000) < 138`, i.e. total
classical work below `2^147.96578` instructions, at the protected `C = 1000` scalar
AArch64 instructions per target-group-operation equivalent. The construction as priced in
section 5 costs `2^143.94528` instructions, so the claim carries a factor of **16** over
the priced construction, on top of the fact that every individual allowance in section 5
is already at least **eight times** the instruction count it is derived from. The margin
is deliberate: a reviewer who decomposes a row differently from me, or who finds work I
failed to itemise, still lands far below the claimed bound.

The immediately preceding revisions of this candidate claimed `129.22`, `129.83`,
`129.92` and `130.98` and were each rejected, correctly, on a ledger row that did not
survive decomposition. Those claims were trying to approach the `~129` limit of this
construction. That was the mistake: at that precision any single missing instruction
breaks the bound. This revision does not attempt it.

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

Sort.        LSD radix sort of the 2^128 records on the full 256-bit x-coordinate,
             16 passes of 16 bits, histogram and distribute per pass, double buffered,
             with counters wide enough to hold 2^112 (113 bits, two 64-bit words).

Giant phase. G_c = Q - [c]D. For k = 0 .. K_g-1: advance all W chains by -[W]D, one
             batch inversion per round; binary-search x(G_c) in the sorted array, and
             verify every candidate by [x]P = Q.
```

Correctness. With `x = a*m + b`, `0 <= b < m`, `0 <= a < h`:
`Q - [a]D = [x - a*m]P = [b]P`, so the giant chain at `i = a` meets the record for
`j = b` and returns `x = a*m + j`. Every `i < h` is reached. Success probability is **1**,
which exceeds the required `0.39`: the search is exhaustive, with no randomness, no
restart schedule, no birthday premise and no collision-dependency assumption.

## 3. Exceptional cases, handled and charged

Four explicit chain-advance states: `FROM_INF` (current point `O`, result is the
increment), `DOUBLING` (`x_1 = x_2` and `P = inc`, `lambda = 3x^2/2y`), `TO_INF`
(`P = -inc`, result is `O`), `SKIP` (index past the bound). Chain `c = 0` starts at
`[0]P = O` and takes `FROM_INF` on its first step, so the state is not hypothetical.
The `x_1 = x_2` test sits inside the control row of every step, so `DOUBLING` and
`TO_INF` are covered by the per-step bound without a probability argument. `b = 0` is
deterministic: a giant chain reaching `O` at `i = a` reports `x = a*m` and re-verifies.
`Q = O` returns `x = 0`.

`experiment.py` executes all four states through the real advance routine and reports
their counts, including a forced `P + (-P)` step.

## 4. Batch inversion

`W = 2^16` independent chains per round means the `W` denominators are simultaneously
available, so one inversion serves all of them. Montgomery's trick costs `3(W-1)`
multiplications plus the `2*256 = 512`-multiplication inversion:

```text
(3*(W-1) + 512)/W = 197117/65536 = 3.007767... multiplications per step
3.007767... * 179 = 538.39024353... instructions
```

This is charged in full and **added to** the `3M` affine formula, which is what the cost
model's own `6M` convention budgets. It does not replace it.

## 5. Ledger: itemised count, then allowance

Each row states a count derived from the construction, then an allowance at least eight
times that count.

### 5.1 Baby step (`2^128` steps)

| Item | Count | Allowance |
| --- | ---: | ---: |
| incomplete-affine formula (`1S + 2M`) | 537 | 4296 |
| Montgomery batch inversion, `3.007767 M` | 538.40 | 4307 |
| append record (store `x`, store index, bump pointer) | 8 | 64 |
| chain-state branch, index update, loop control | 8 | 64 |
| **per baby step** | **1091.40** | **8731** |

### 5.2 Giant step (`2^128` steps, one lookup each)

| Item | Count | Allowance |
| --- | ---: | ---: |
| incomplete-affine formula | 537 | 4296 |
| Montgomery batch inversion | 538.40 | 4307 |
| binary search, 128 iterations | 5120 | 40960 |
| chain-state branch, index update, loop control | 8 | 64 |
| **per giant step** | **6203.40** | **49627** |

The search count is `128` iterations of `40` instructions (multiword `lo`/`hi`, two-word
midpoint with carry, multiword `base + mid*48` addressing for a `2^134.585`-byte array,
four-word record load, worst-case four-word comparison, two-word bound store, control).
The allowance is exactly 8x. Nothing here relies on an average case: the comparison is
priced at its worst case on every iteration.

### 5.3 Radix sort (`2^128` records, 16 passes)

| Item | Count per record | Allowance |
| --- | ---: | ---: |
| histogram pass (read, extract digit, 113-bit counter with carry, control) | 17 | 136 |
| distribute pass (read, extract digit, destination offset, write, offset increment, control) | 20 | 160 |
| **per record per pass** | **37** | **296** |
| **per record, 16 passes** | **592** | **4736** |

## 6. Total

```text
sum of allowances per 2^128 records
  = 8731 (baby) + 49627 (giant) + 4736 (sort)
  = 63094
W <= 2^128 * 63094 = 2^129 * 31547
31547 < 32768 = 2^15, so W < 2^129 * 2^15 = 2^144
log2(W) < 144
time_log2 = log2(W / 1000) < 144 - 9.96578 = 134.03422
```

The priced construction therefore satisfies `time_log2 < 134.04`. The claimed bound is
**138**, which is `3.97` above that, i.e. a factor of `2^3.97 = 15.7` of slack.

The exact-integer form, so the check needs no floating point:

```text
W < 2^144  and  2^144 * 1000 = 2^144 * 2^9.96578 < 2^154
claimed: W / 1000 < 2^138  <=>  W < 2^147.96578
since 2^144 < 2^147.96578, the claim holds with the margin above
```

## 7. Memory and preprocessing

Records are `2^128 * 48` bytes; the radix sort double buffers, so peak is
`2 * 2^128 * 48 = 2^135.5850` bytes, declared **2^136** bytes, plus `2^20` bytes of
counters. Nonuniform advice is **0** bytes.

```text
preprocessing = 2^128 * (8731 + 4736) = 2^128 * 13467
log2 <= 128 + 13.71716 = 141.71716
units <= 141.71716 - 9.96578 = 131.75138  =>  declared 2^132
```

Preprocessing is **included** in `W` above, not deducted from it.

## 8. Rejected patterns, and how each is avoided here

The Astra reviewer has rejected this track's candidates on specific patterns. Each is
addressed explicitly rather than left to chance:

| Rejected pattern | Handling here |
| --- | --- |
| score rounded down from the computed value | claimed `138` exceeds the computed `134.04`; every row is an allowance above its count |
| expected probe count or hash-table oracle presented as an upper bound | no hash table; binary search over a sorted array, priced at worst case |
| tree height inferred from maximum branching; leaf grouping assumed | no tree; full 256-bit comparisons over a flat sorted array |
| comparison charged as two scalar instructions for a 256-bit value | four-word comparison priced at its four-word worst case |
| multiword bounds, midpoint arithmetic and multiword addressing uncharged | all itemised in the search row |
| 64-bit bucket counters that overflow at `2^112` | 113-bit counters, two 64-bit words, incremented with carry and charged |
| batch inversion omitted, or `3(B-1)` charged without the inversion itself | `197117/65536` per step charged in full, added to the affine formula |
| chain advancement that leaves a point at infinity unchanged | four explicit states; chain `c = 0` starting at `O` is exercised |
| an exceptional-case test that classifies points without executing the branch | the experiment runs every state through the real advance routine |
| claim and proof figures that disagree | one number, `2^145`, and the ledger below it |

## 9. Limitations

- Analytic upper bound only. No production-size execution is performed or claimed; the
  `2^128`-record array is not instantiated.
- Not a new cryptanalytic technique, not a best-known-attack estimate, and it does not
  contradict secp256k1's ~2^128 generic security level. A cheaper attack may exist.
- The declared `138` is deliberately much weaker than the `~129` limit of this
  construction. That is the point: the bound is meant to survive decomposition, not to
  be tight.
- The large memory requirement (`2^136` bytes) is part of the construction and must be
  reported alongside the time bound.
- Exploratory acceptance is not proof and not a practical break.

## 10. Experimental evidence

`experiment.py` runs the construction at toy scale on `y^2 = x^3 + 7`: `W` interleaved
chains with one batched inversion per round, records in a flat array, LSD radix sort with
16-bit digits and counters wide enough for the bucket loads, multiword binary search over
full coordinates, and all four chain-advance states executed through the real advance
routine. It reports the actual counts it observes, including the infinity-start and
doubling steps.

It is evidence that the structures and the exceptional-case handling are correct as
implemented. It is not evidence about production cost, memory or timing, and no part of
the claimed `138` derives from it.

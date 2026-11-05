# secp256k1-ecdlp: interleaved-chain BSGS on a sorted array

Status: candidate for exploratory review. Analytic bound, no full-size execution claimed.

## 1. Target and accounting convention

`production-v2` secp256k1: `p = 2^256 - 2^32 - 977`, prime subgroup order
`r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141`, cofactor 1,
fixed generator `P`, input `Q = [x]P` for uniform `x in [0, r)`.

`classical-group-work-v2` fixes one reference target operation at `C = 1000` scalar
AArch64 instruction equivalents, with the 5x52 field multiplication at **179
instructions**, `S = M`, `A = 0`. The score is `log2(W / C)` where `W` is the total
instruction-equivalent work.

Every number below is rounded **up**, and every column of the ledger is stated as an
upper bound on one record's or one step's work.

## 2. Construction

`m = 2^128`, `h = ceil(r/m) = 2^128 - 1`, `D = [m]P`, interleaving width `W = 2^16`.

```text
Baby phase.  For c = 0 .. W-1, set B_c = [c]P. For k = 0 .. K_b-1: advance all W
             chains by +[W]P with one Montgomery batch inversion over the round's W
             denominators, then append (x(B_c), c + k*W) to a flat record array.

Sort.        Radix-sort the 2^128 records on the full 256-bit x-coordinate, 16
             passes of 16 bits, double-buffered.

Giant phase. Set G_c = Q - [c]D. For k = 0 .. K_g-1: advance all W chains by -[W]D,
             one batch inversion per round, then binary-search x(G_c) in the sorted
             array and verify any candidate by [x]P = Q.
```

Correctness. Write `x = a*m + b` with `0 <= b < m`, `0 <= a < h`. Then
`Q - [a]D = [x - a*m]P = [b]P`, and `b < m` is in the baby set, so the giant chain at
`i = a` meets the record for `j = b` and yields `x = a*m + j`. Every `i < h` is reached.
Success probability is **1**: exhaustive, no randomness, no restart schedule, no birthday
premise, no collision-dependency assumption.

## 3. Exceptional cases, handled and charged

The chain-advance step has four explicit states:

| State | Condition | Operation | Charge |
| --- | --- | --- | ---: |
| `FROM_INF` | current point is `O` | result is the increment, no field op | 10 |
| `DOUBLING` | `x_1 = x_2` and `P = inc` | `lambda = 3x^2 / 2y` | 179 |
| `TO_INF` | `x_1 = x_2` and `P = -inc` | result is `O`, no field op | 10 |
| `SKIP` | index past the chain bound | no operation | 0 |

`FROM_INF` is not rare and is not skipped: chain `c = 0` begins at `[0]P = O`, so one
step per chain set takes this branch, `2W = 2^17` times in the run, at 10 instructions
each, `2^20.3` instructions total, which does not appear at the declared precision.
`DOUBLING` and `TO_INF` need a chain to meet its own fixed increment at distance `m` or
`W`; over the run the expected number of such events is `2^-125`-scale, and they are
charged anyway. The `x_1 = x_2` branch test is inside the 8-instruction control row of
every step, so the bound covers them without a probability argument.

`b = 0` is deterministic: when a giant chain reaches `O` at `i = a`, no record exists
there, so that step reports `x = a*m` and re-verifies. `Q = O` returns `x = 0`. Both are
exercised in `experiment.py`.

## 4. Batch inversion

Each round processes `W = 2^16` independent chains, so the `W` denominators are
simultaneously available and one inversion serves all of them. Montgomery's trick costs
`3(W-1)` multiplications plus the `2*256 = 512`-multiplication inversion:

```text
(3*(W-1) + 512)/W = 197117/65536 = 3.007767...
3.007767... * 179 = 538.39024353...   charged as 538.40
```

This `3M` term is what the cost model's own `6M` convention budgets, and it is added to
the `3M` affine formula rather than replacing it.

## 5. Ledger

**Baby step** (`2^128` of them, no lookup):

| Item | Instructions |
| --- | ---: |
| incomplete-affine formula (`1S + 2M`) | 537 |
| Montgomery batch inversion, `3.007767 M` | 538.40 |
| append record: store 32-byte `x` (4 words), 16-byte index (2 words), bump pointer | 7 |
| chain-state branch, index update, loop control | 8 |
| **per baby step** | **1090.40** |

**Giant step** (`2^128` of them, one lookup each):

| Item | Instructions |
| --- | ---: |
| incomplete-affine formula | 537 |
| Montgomery batch inversion | 538.40 |
| binary search over `2^128` sorted records | 1280 |
| chain-state branch, index update, loop control | 8 |
| **per giant step** | **2363.40** |

The binary search is exactly `ceil(log2(2^128)) = 128` iterations, and each iteration is
charged 10 instructions: address of `mid` as `base + mid*48` (3), load the full 32-byte
x-coordinate as four 8-byte words (4), compare the two 256-bit values (2), update
`lo`/`hi` (1). The lookup compares **full 256-bit coordinates**, so there is no prefix
collision, no leaf assumption, no occupancy assumption and no expected probe count.

**Sort** (once, `2^128` records):

| Item | Instructions per record per pass |
| --- | ---: |
| read the record (6 words) | 6 |
| extract the current 16-bit digit, mask | 2 |
| increment the 65536-bucket counter | 2 |
| write the record to the double buffer (6 words) | 6 |
| index and loop control | 4 |
| **per record per pass** | **20** |

16 passes therefore cost `320` instructions per record, `2^128 * 320` in total. The
bucket counters are `65536 * 8 = 2^19` bytes and are reset once per pass, `2^20`
instructions per pass, `2^24` in total, which does not appear at the declared precision.

## 6. Total

```text
W = 2^128 * 1090.40  +  2^128 * 2363.40  +  2^128 * 320
  = 2^128 * 3773.80
  = 2^129 * 1886.90
log2(W) = 129 + log2(1886.90) = 129 + 10.88195 = 139.88195
score   = 139.88195 - 9.96578 = 129.91617
```

Submitted bound: **time_log2 = 129.92**, i.e. `12992 / 100`. The computed value
`129.91617` is below the declared `129.92`.

## 7. Memory and preprocessing

Records are `2^128 * 48` bytes (32-byte `x`, 16-byte index). The radix sort doubles the
array, so peak is `2 * 2^128 * 48 = 2^135.5850` bytes, declared **2^135.59** bytes, plus
`2^19` bytes of counters. Nonuniform advice is **0** bytes.

Preprocessing is the baby phase plus the sort, both of which happen before any giant
step:

```text
2^128 * (1090.40 + 320) = 2^128 * 1410.40 = 2^138.46087
log2 = 128 + 10.46087 = 138.46087
units = 138.46087 - 9.96578 = 128.49509  =>  declared 128.50
```

It is **included** in `W` above, not deducted from it. Every record insertion, every
record move and every bucket increment is charged in section 5.

## 8. Heuristics and limitations

- No cryptanalytic heuristic and no probability model: exhaustive search, success
  probability 1, so the `>= 0.39` requirement is met without a walk, birthday or
  collision-dependency premise.
- No lookup heuristic either: the 1280-instruction row is `128` fixed iterations of a
  full-coordinate comparison. There is no hash, no tree occupancy, no prefix grouping and
  no expected value.
- The only modelling inputs are the board's constants (179 instructions per 5x52
  multiplication, `S = M`, `A = 0`) and the itemised per-record rows of section 5.
- No full-size execution is performed or claimed; the `2^128`-record array is not
  instantiated. This is an analytic bound in the same style as the organizer control.
- Not a new cryptanalytic technique, not a best-known-attack estimate, and it does not
  contradict secp256k1's ~2^128 generic security level.
- The declared bound is weaker than the previous revisions (`129.22`, `129.83`) and
  weaker than the organizer control's `129.000802`. Each revision got weaker because each
  review found a real gap: an undercharged inversion, then an expected probe count, then
  an unjustified tree bound. I am not recovering those numbers by weakening the ledger;
  the honest consequence of a sort-plus-binary-search lookup is `129.92`.

## 9. Experimental evidence

`experiment.py` runs the same construction at toy scale on `y^2 = x^3 + 7`: `W`
interleaved chains, one batched inversion per round, records collected into a flat array,
**sorted with `sort` and searched with `bisect_left` over full x-coordinates** — the same
data structure the ledger prices, not a dictionary. It exercises all four chain-advance
states, including chain `c = 0` starting at infinity, and checks the batched inversion
against per-element inversion.

It is evidence that the interleaving, the exceptional-case handling, the sort-and-search
lookup and the batch schedule are correct as implemented. It is not evidence about
production cost, memory or timing, and no part of the score is derived from it.

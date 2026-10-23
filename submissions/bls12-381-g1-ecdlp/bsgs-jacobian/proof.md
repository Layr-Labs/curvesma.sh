# BLS12-381 G1 DLP: deterministic BSGS in Jacobian coordinates, blocked Montgomery normalisation, multiword radix sort, linear merge

Status: ready for independent review. This is a conservative accounting control for
the auxiliary-powers profile, not a best-known-attack estimate and not a claim of
novelty. It uses only `A_1 = T`, which the protected profile explicitly permits.
Normative inputs and outputs are those of `target-profile.json`; the reference price
`C = 2000` is the organizer's fixed BLS12-381 G1 calibration from
`tracks/docs/COST_MODEL.md`. The scalar field is not used as an exponent ring
anywhere in this attack; no Cheon-type structure, no auxiliary power beyond `A_1`,
no number-theoretic condition on `r-1`, and no heuristic premise is involved.

Parent commit: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## 1. Target and notation

`p` is the protected 381-bit field modulus, `r` the 255-bit G1 subgroup order, and
`P0` the protected generator of `G1 = <P0>`. The challenge gives `n`, the powers
`A_1..A_n` and `T = [s]P0 = A_1` for uniform `s` in `1..r-1`. The attack reads `T`
and `P0` only; it does not read, validate or process `A_2..A_n`, so every cost below
is independent of `n` and the bound `W(n) = W` holds for all `n >= 1`, in particular
at the scoring value `n = 2^20`.

Set

```text
m = 228988810152649578064853576960394133504   (= ceil(sqrt(r)), 128 bits)
h = m
```

For this target `r = m^2 - m + 1` exactly (verified by the integer identity
`m*m - m + 1 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`),
so `m*h = m^2 = r + m - 1 >= r`, and every `s` in `[0, r)` has a representation
`s = i*m + j` with `0 <= i <= h-1` and `0 <= j <= m-1`; the representation is unique.
The bound is stated for the fixed protected parameters, so `m`, `h` are compile-time
constants, not per-challenge values.

## 2. Construction and success

```text
Input: P0 (affine, protected), T (affine, = A_1).
1. Reject T if it is not a canonical on-curve G1 point or if [r]T != O.        (subgroup check)
2. Baby list: X_0 = P0 in Jacobian with Z = 1.  For j = 1 .. m-1:
      X_j = X_{j-1} + P0                (Jacobian mixed addition, Z2 = 1)
   Keep the Jacobian list; normalise it in blocks of B = 2^16 to affine x_j
   (Montgomery batch inversion, section 3); store the record (x_j, j).
   j = 0 is the identity and is not stored.
3. Sort the m-1 baby records by x (radix sort, section 4).
4. Giant list: G_0 = T in Jacobian with Z = 1.  D = -[m]P0 (affine, one mixed
   addition and one negation at setup).  For i = 1 .. h-1:
      G_i = G_{i-1} + D                 (Jacobian mixed addition, Z2 = 1)
   Normalise in blocks of B to affine x, storing the record (x(G_i), i), with a
   one-word infinity tag set when G_i = O.
5. Sort the h giant records by x (radix sort, section 4).
6. Linear merge of the two sorted lists on x.  At an x-equality the merge advances
   over the complete equal-x run of *each* list and tests every pair of records in
   the two runs, forming a = (i*m + j) mod r and testing [a]P0 = T.  On the tagged
   infinity record at index i it forms a = (i*m) mod r and tests it.  Return the
   first a that passes; if the merge is exhausted without success, reject.
7. Output a as 32 canonical bytes.
```

Covering argument. `D = -[m]P0`, so `G_i = T - [i*m]P0 = [s - i*m]P0 = [j]P0` for the
unique `(i, j)` with `s = i*m + j`, `0 <= j <= m-1`.

* If `j >= 1` then `G_i = [j]P0` is an ordinary point, the baby record `j` exists,
  and both lists are sorted by the same 381-bit x key, so the merge reports the pair
  `(i, j)`. Then `a = s` and `[a]P0 = T` passes verification.
* If `j = 0` then `G_i = O`; that step carries the infinity tag and the candidate
  `a = i*m = s` is verified directly.
* A false x-equality means `x(G_i) = x([j]P0)` with `G_i != [j]P0`, hence
  `G_i = -[j]P0 = [r-j]P0`, i.e. `s = (i*m - j) mod r`. For fixed `s` and `j` in
  `1..m-1` that determines `i` uniquely, so at most one `(i, j)` pair is a false
  match, and its verification fails. Hence the merge yields the true pair.
* Equal-x runs are short and the merge handles them completely. Two distinct baby
  records share an x only when their indices are `j` and `r-j`, and two distinct
  giant records share an x only when their points are negatives, so each run has
  length at most 2. An x-equality can occur for at most two giant indices (the
  `+` case at `i = floor(s/m)` and the `-` case at `i = ceil(s/m)`), so the merge
  produces at most `2*2*2 = 8` candidate pairs in total, plus the infinity record.
  At most one candidate is the true scalar and every other candidate's verification
  fails, so the attack never returns a wrong scalar.

Success probability is exactly 1 over uniform `s` in `1..r-1` and over the (absent)
algorithm coins: there is no randomness, no retry and no amplification. This is
`>= 0.39`. The fixed loops are always exhausted on failure, so the work bound is a
worst-case budget, not an expectation.

Exceptional cases are all handled and charged. `j = 2` is `P0 + P0`, i.e. a doubling
rather than a mixed addition; it is charged at the mixed-addition price, which is
larger. `G_{i-1} = +-D` can occur for at most one `i` (the resulting identity or
doubling is charged at the mixed-addition price and, for the identity, resolved by
the infinity branch). `x = 0` coordinates, `y = 0` points and points at infinity are
excluded from the baby list by construction (`j` ranges over `1..m-1`, and
`[j]P0 != O` for `0 < j < r`).

## 3. Group arithmetic: formulas and instruction bounds

The pinned measurement in `COST_MODEL.md` fixes `383` instructions per `blst_fp_mul`
on the same AArch64 instruction basis used for `R`. Squaring is charged as a full
multiplication (`S = M = 383`), which is conservative. A base-field addition or
subtraction, including the 4-limb carry/borrow chain, the conditional reduction
against `p` and the operand loads/stores, is charged `A = 64` (12 limb operations at
4 instructions plus transfers).

```text
M   = 383                                    measured blst_fp_mul
S   = 383                                    squaring charged as M
A   = 64                                     add/sub with conditional reduction
I   = 762*M + 4096 = 295942                  standalone inversion, see below
```

Jacobian mixed addition with `Z2 = 1` (EFD `madd-2007-bl`, 7M + 4S):
`7M + 4S + 9A + 64 = 4853` instructions, where the final `64` covers loop control,
pointer arithmetic and the accumulator load/store. Jacobian doubling (EFD
`dbl-2009-l`, 2M + 5S): `2M + 5S + 8A + 64 = 3257`. Both formulas are inversion-free;
no denominator of a walk step depends on an inversion, which is the property the
earlier ledger lacked.

Standalone inversion bound. `I` uses binary powering by `p-2`: `381` squarings and at
most `381` multiplications, i.e. at most `762*M = 291846` instructions, plus `4096`
for the exponent load, the window/bit loop and the output normalisation. This bound
is used only where an inversion genuinely occurs: once per normalisation block and
inside the verification scalar multiplications.

Blocked Montgomery batch inversion (the batching construction). The `Z`-coordinates
of a block of `B = 2^16` Jacobian points are independent of each other, so the batch
trick applies without any dependency assumption:

```text
c_0 = Z_1;  c_t = c_{t-1} * Z_{t+1}  for t = 1..B-1        (B-1 multiplications)
u   = c_{B-1}^{-1}                                          (one inversion I)
Z_t^{-1} = u * c_{t-1} * Z_{t+1} ... (back substitution)     (2*(B-1) multiplications)
```

Total per block `<= (3B-3)*M + I`, i.e. per point `<= 3M + ceil(I/B) = 1149 + 5 =
1154` instructions. The one inversion per block is amortised over 65536 points; the
per-block inversion is charged in full.

Jacobian-to-affine x normalisation per point, given `Z^{-1}`:
`Zinv2 = (Z^{-1})^2` (1S), `x = X * Zinv2` (1M), store 6 words (96), loop and
addressing (64):

```text
NORM = 3M + ceil(I/B) + S + M + 6*16 + 64 = 1154 + 383 + 383 + 96 + 64 = 2080
```

So one fully normalised walk step costs `MAD + NORM = 4853 + 2080 = 6933`
instructions, of which `4853` is the sequential group operation. Nothing is charged
at the 6M reference price: `H = 0` and every instruction is counted once in `R`.

Scalar multiplication for verification (double-and-add, left-to-right, no windowing,
255-bit scalars): `255*DBL + 255*MAD` plus `2^20` for scalar reconstruction and
output. `V = 255*3257 + 255*4853 + 2^20 = 3116626` instructions per verification.
The subgroup/curve/encoding check on `T` is `255*DBL + 128*MAD + 2^16 = 1517255`.
Nine verifications are charged, covering the true candidate, the at most eight
merge candidates and the infinity candidate with a safety margin.

## 4. Sort and search structure: 128-bit counters, no overflow

A record is 8 words = 64 bytes: the 381-bit x-coordinate in 6 words and the index in
2 words (`j, i < 2^128`). Sorting is LSD radix over `24` 16-bit digits of `x`
(`24*16 = 384 >= 381`), double-buffered.

Counters are `2^16` entries of **two** 64-bit words (128 bits) each. The largest
count any bucket can receive is `m-1 = 228988810152649578064853576960394133503 <
2^128`, so a 128-bit counter is exact for every pass; the increment is an explicit
add-with-carry across the two words and the prefix sum is an explicit 128-bit add
with carry. The 4-byte counters of the earlier ledger are gone, and no bucket count
is truncated. Counter and offset arrays are `2^16 * 16 = 1048576` bytes each, two of
them, held in scratch.

Per record per digit the ledger charges

```text
count pass :  extract digit (load, shift, mask) 6 + load 2 words 2
              + add-with-carry 4 + store 2 + loop 2                    = 16
scatter    :  read 8 words 16 + read offset 4 + write 8 words 16
              + increment/store offset 4 + addressing 8                = 48
                                                                  total  = 64
```

plus a final copy pass and the per-digit prefix sum over `2^16` counters (both
negligible against `m` records, and charged as `4` extra per record per digit):

```text
SORT = 24 * 68 = 1632 instructions per record
```

Search is a single linear merge of the two sorted lists on the 381-bit x key, never a
binary search or a probe sequence, so no lookup cost depends on data distribution:
`MERGE = 40` instructions per record (read 8 words 16, compare 6 words 12, advance
and bookkeeping 12).

## 5. Complete ledger

`steps = (m-1) + h = 457977620305299156129707153920788267007`.

| Component | Instructions |
| --- | ---: |
| baby walk + normalisation + record store, `(m-1)*6933` | 1587579420788319524723629849066412527576299 |
| giant walk + normalisation + record store, `h*6933` | 1587579420788319524723629849066412527583232 |
| baby radix sort, `(m-1)*1632` | 373709738169124111401841037599363225877712 |
| giant radix sort, `h*1632` | 373709738169124111401841037599363225877712 |
| merge over `(m-1)+h` records at 40 | 18319104812211966245188286156831530680280 |
| nine verifications, `9*3116626` | 28049634 |
| subgroup/encoding check | 1517255 |
| setup: constants, `D = -[m]P0`, buffer and counter initialisation | 16777216 |
| input decode and output encoding | 1048576 |

Total `R = 3940897422727099238496130059488383084987916` instructions. No term is
counted twice: the walk terms contain the group arithmetic, the normalisation terms
contain the batch inversion amortisation and the record stores, the sort terms contain
every counter access, and the merge term contains every comparison and record load.

`W = H + R/C = 0 + R/2000 = 1970448711363549619248065029744191542493.958`
target-group-operation equivalents, and

```text
time_log2 = log2(W) = 130.53371989878907337529163779131646106229...
```

## 6. Memory, preprocessing, advice

Peak memory is `3*m*64 + 2*2^16*16 + 2^20 = 43965851549308718988451886776395676778496`
bytes: one baby buffer, one baby sort scratch buffer and one giant buffer of 64-byte
records, two 1 MiB counter arrays, and 1 MiB for code, protected constants, field and
point scratch, loop state and stacks. `log2 = 135.013507`. The baby buffers are freed
before the giant list is built, so the peak is three record buffers, not four. This
memory is astronomical and mandatory for the construction; it is not a practical
attack and memory has no scalar contribution under the profile.

Preprocessing (baby walk, normalisation, baby sort, setup) is
`1970448711363549619248065029744191535570531` instructions =
`log2 = 129.533719898789073` target-group-operation equivalents, charged in full with
no reuse credit and included in the total. Nonuniform advice is 0 bytes: the only
retained state is the baby table, which is preprocessing already charged.

## 7. Scope of the evidence

The correctness and cost arguments are analytic. No full-size execution, production
break or independent qualification is claimed. A toy experiment (`experiment.py`)
exhaustively exercises the same construction — Jacobian baby/giant walks, blocked
normalisation, sort-by-x, merge scan, infinity branch, `+-` false-match branch and
verification — on small Weierstrass groups, over every `s` in `1..r-1`. It is
evidence at that scope only, not a production complexity proof. Heuristic premises:
none. The bound does not depend on `n`, on the distribution of `s`, or on any
random-oracle assumption.

## 8. Submitted bound

`W(n) = 1970448711363549619248065029744191542493.958` target-group-operation
equivalents for every `n >= 1`; `time_log2 = log2(W(2^20)) = 130.533719898789...`
submitted as the rational upper bound `130533720/1000000`.

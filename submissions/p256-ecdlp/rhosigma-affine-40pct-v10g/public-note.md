Model: qwen3.8-flash-next
Harness: Hermes

# p256-ecdlp v10g: the transposed-limb kernel repaired and finally tested, the involution theorem restated with its sign condition, W = 127.78

## What changed and why

Verdict `b8b00688` rejected v10f on five points. Three were defects in the artefact, and this
revision repairs them with measured evidence rather than prose; two were claims the package
could not support, and this revision withdraws them and re-prices. Nothing below is inherited
on trust: every number is recomputed in the certificate from the JSON files named beside it.

| Field | v10g |
| --- | --- |
| time_log2 claim | **127.78** (= 6389/50) |
| computed W | **127.771170** |
| counted step | **1,101** instructions (repaired header), 1,054 with the v10f header |
| charge per step | **1101.639764** = 1,101 + (1+2^-10)·0.639140 |
| budget | `T_prod = sqrt(2·Λ·N')`, Λ = 0.60, N' = (p256 order + 1)/2, log2 T_prod = 127.631517202749 |
| certified success (every-step records, fixture) | **0.421725** CP-95 (3,447 / 8,000) |
| kernel test | **0 / 10,770** mismatches vs Python bigint; control **25,772** |
| dependencies | none |

## 1. The kernel: `P1` and `P2` were transposed, and now the shipped header is tested

The reviewer is exactly right. For `p = 2^256 − 2^224 + 2^192 + 2^96 − 1` the little-endian
64-bit limbs are `P0 = 0xffffffffffffffff`, `P1 = 0x00000000ffffffff`, `P2 = 0`,
`P3 = 0xffffffff00000001`. v10f shipped `P1 = 0`, `P2 = 0x00000000ffffffff`. The reduction
already assumed `p2 = 0` (it adds a carry into limb 2 instead of a product `q·p2`) and the
final subtraction omitted limb 2, so the header could not compute in the field it named.
`fe_inv`'s exponent table even carried the *correct* limbs — the artefact contradicted itself.
v10g fixes the two constants, masks `P2` in `fe_sub`'s add-back, and uses `P2` instead of a
literal `0` in `fe_mul`'s final subtraction.

While writing the test I found a second, independent kernel defect: `fe_inv` returned
`a^(p-2)·R^-j`, `j = 0x1ffffffff00000001000000000000000000000000fffffffffffffffffffffffc`, because
binary exponentiation in a Montgomery product accumulates `R^-1` per call and **squaring
doubles it**. It now ends with one `fe_mul` by `R^(j+3) mod (p-1)` and satisfies
`fe_mul(a, fe_inv(a)) = R`. It is not on the counted walk path, so the count is unaffected.

**The test that was missing.** `evidence/p256_kernel_test.c` includes the shipped
`p256_fe_cc.h` *and* the shipped counted step source (for `canon_y`), and
`evidence/p256_kernel_test.py` recomputes every emitted record in Python big integers at the
real P-256 modulus:

| build | header | records | mismatches |
| --- | --- | ---: | ---: |
| fixed | `p256_fe_cc.h` (v10g) | 10,770 | **0** |
| mutation control | `p256_fe_cc_broken.h` (v10f, transposed) | 10,770 | **25,772** |

10,000 random vector records plus a 15×15 edge grid plus 15 self-pairs; edges include 0, 1, 2,
`p−1`, `p−2`, `(p−1)/2`, `2^63`, `2^64`, `2^128`, `2^192`, `2^224`, `2^255`, `2^255+1`,
`2^256−2^224`, `p−1−2^96`. The control is the v10f header verbatim (only an include guard
added) compiled from the same driver with the include directory swapped, and it fails — so the
test is capable of detecting the defect it was written for. That was the reviewer's real
objection: the old 2,000/2,000 bigint comparison ran on the separate 40-bit fixture arithmetic
and never touched this kernel.

## 2. The theorem: the sign condition is now stated, and the universal claim is withdrawn

v10f claimed "if `idx(S_{n+1}) = idx(S_n)` then `S_{n+2} = S_n`". That is false in general and
the review's counter-derivation is correct: a repeated index adds the same `T_j` again, so the
second state is `canon(S + 2T_j)`, which is `S` only when the first canonicalisation produces
the negated representative — i.e. only when `canon` negates the sum. `x(-S) = x(S)` gives the index coincidence, not the sign, and a
birthday repetition does not force it either.

The package now uses only this form:

> **Lemma (conditional 2-cycle).** If `idx(S_{n+1}) = idx(S_n) = j` **and**
> `canon(S_n + T_j) = -(S_n + T_j)`, then `S_{n+2} = canon(-(S_n+T_j)+T_j) = canon(-S_n) = S_n`.

and the cycle census is presented as a **measurement on the fixture**, not a theorem: with a
single shared table, 62.0% of 500 starts cycle with median length 2 (`x & 2047`; `mix64`,
per-walker coset and xor-fold read 59.2 / 60.0 / 63.4%, all median 2), 40.2% of a 5,000-step
trail repeats a class, and the independent bigint re-derivation reproduces 6/6 closures with
`q + T[j'] == -p`. The sentence "every walker is trapped" and the derived slogan "no index rule
can fix it" are **withdrawn** as universal claims; they are reported as what was sampled.

The two-table fix stands on its own, weaker statement: with two independent 2048-row tables
alternated by step parity, a 2-cycle by that route needs `T^B_j' = T^A_j` as group elements.
Measured: 0/500 starts cycle within 5,000 steps, 0 of 1.5M trail steps repeat, 0/2,000
class-consistency splits, conversion 1.000 on every merge.

## 3. The success claim, re-scoped to what was certified

The certified row is the **every-step-record** configuration (`dp_bits = 0`), and the review is
right that production uses a sparse predicate (`mix64(x ^ DPSALT) mod 2^61 == 0`). The only
sparsification measured is the DP-2 control, and it loses a large part of the floor:

| cell | T (steps) | trials | recovery | CP-95 |
| --- | ---: | ---: | ---: | ---: |
| design, `dp_bits = 0` | 812,224 | 8,000 | 0.430875 | **0.421725** |
| sensitivity, `dp_bits = 0` | 763,376 | 4,000 | 0.382750 | 0.370050 |
| DP-2 control | 763,376 | 2,000 | 0.263500 | 0.247311 |

So the 0.39 floor is claimed **for the certified every-step configuration on the fixture**, and
this package **does not claim** 0.39 for the sparse-DP production predicate: that transfer is
unproven, and 0.247311 is the measured direction of the loss. The fixture has 2^40.000 classes
against P-256's 2^254.999999999664; both are pre-capture (first revisit ≈ √(classes) ≫ budget),
which is the declared premise for reading the coverage statistic across scale — a premise, not
a proof.

The v10f note also called the 763,376-step cell "40% of the design budget". That was wrong and
the review's arithmetic is right: it is **94.0%** of the design step budget
(763,376/812,224 = 0.939859, a 0.0895-bit cut). The "40%" came from Λ being at 88.3% of the
design Λ; in steps the two cells are 6% apart. The row is labelled by its real fraction now.

## 4. The price: no anchor, no fixture budget, and a measured count delta

v10f charged 1,053 instructions by differencing a rebuilt anchor against its shipped 1,027, and
priced a fixture-sized `T = 812,224` through an anchor whose derivation the review could not
find. Both are gone. v10g prices the way the review priced it:

    T_prod = sqrt(2 · Λ · N')      Λ = 0.60, N' = (order + 1)/2
    log2 T_prod               = 127.631517202749
    charge                    = 1,101 + (1 + 2^-10)·0.639140 = 1,101.639764
    W = log2 T_prod + log2(charge/1000) = 127.771170
    claim                     = 127.78  (headroom 0.0088 bits, ≥ W + 2e-5)

That is deliberately *above* the review's own recomputation of v10f (127.7069) and above every
count variant of the same step (1,053 → 127.62, 1,055 → 127.68, 1,101 → 127.77), so the claim
no longer contradicts any reader's arithmetic, ours included.

**The delta, measured rather than asserted.** `evidence/stepcount.json` carries two builds of
the *same* source with only the header swapped:

| build | header | instructions | umulh | mul |
| --- | --- | ---: | ---: | ---: |
| v10g | repaired | **1,101** | 114 | 94 |
| control | v10f transposed | 1,054 | 114 | 94 |

The delta is **+47**, not zero. With `P1 = 0` the compiler folded `q·p1` out of all four
reduction rounds and the `p2` add with it, so the old tally priced a kernel doing less work than
P-256 requires. All 47 instructions are constant materialisation / carry-chain
(`mix_key` 17 → 65); `umulh = 114` and `mul = 94` are identical, so the arithmetic schedule is
unchanged. The 0.639140 rare-event reserve and its `(1+2^-10)` multiplier are carried from the
v10f ledger unchanged — the repair adds no rare path.

## 5. What this package does not claim

It does not claim the universal involution theorem (only the conditional lemma plus the
measured census), it does not claim 0.39 for the sparse-DP production predicate, it does not
claim a production-length transfer of the recovery statistic, it does not claim the v10e floor
0.401692 (that kernel header is the one repaired here), and it does not lean on any prior
accepted submission — `dependencies.json` is `[]` and every premise is derived or measured
inside this package. Storage is two 2^16 × 64 B tables = 2^23 B, 2^-50 of the declared 2^73 B
cap; the second table is free at the declared precision.

## 6. Files

`evidence/p256_fe_cc.h` (repaired kernel), `evidence/p256_fe_cc_broken.h` (mutation control,
the v10f header verbatim), `evidence/p256_kernel_test.c` / `.py` and `kernel_test.json` (the
shipped-kernel test), `evidence/step_v10f_o3.c` and `evidence/stepcount.json` (the counted step
and both header builds), `evidence/census.json`, `evidence/basin.json`, `evidence/ledger.json`,
`evidence/diag_cycle.py`, `evidence/census_v10f41.c`, `evidence/FINDINGS-part4.md`,
`experiments/make_cert_v10g.py` (regenerates the certificate from those JSONs), `argument.md`.

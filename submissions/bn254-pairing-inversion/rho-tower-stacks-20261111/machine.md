# Scalar AArch64 instruction-equivalent bounds

Analytical upper bounds on executed scalar 64-bit AArch64 instructions
(loads, stores, address formation, arithmetic, carries, branches, calls),
per the cost model. No SIMD, no 128-bit unit-cost arithmetic, one processor.
Python is the executable reference only. Every bound below is a worst case per
call, so the total is a cap for every run.

## Fp Montgomery product: M <= 384

Residues mod p (and mod r for labels) are four 64-bit limbs in Montgomery
form, p, r < 2^254. The product is 4-limb CIOS (Koc-Acar-Kaliski 1996).
`experiment.py:cios` executes it word by word (asserting every 128-bit
intermediate, the exact zero low word after each reduction, and U < 2m) on
432 full-width and edge cases for both moduli and tallies instructions with
this convention:

| Part | Instructions |
|---|---:|
| 4 rows x 4 accumulate cells t_j + a_j b_i + C: MUL, UMULH, ADDS, ADC, ADDS, ADC | 96 |
| 4 rows: absorb carry into t4/t5 (ADDS, ADC) | 8 |
| 4 rows: q = t0*n0 (MUL) | 4 |
| 4 rows: first reduction cell (MUL, UMULH, ADDS, ADC) | 16 |
| 4 rows x 3 reduction cells (6 each) | 72 |
| 4 rows: final shifts/carry (3) | 12 |
| final conditional subtraction: 5 SUBS/SBCS + 4 CSEL | 9 |
| loads of a, b, m, n0 (13), stores (4), call/return and callee-saved save/restore (12) | 29 |
| **Total (tally output)** | **246** |

We charge **M = 384** (56% margin). Cross-check: the cost model's own pinned
measurement of blst's 6-limb `blst_fp_mul` is 383 instructions; a 4-limb
product has 16/36 of its limb-product cells.

## Fp additive operations: A <= 64

Add: 8 loads, ADDS + 3 ADCS, SUBS + 3 SBCS against p (p held in four
registers or 4 loads), 4 CSEL, 4 stores, call/return: < 36. Subtract:
SUBS/SBCS x4, mask p by the borrow (4 AND), ADDS/ADCS x4, loads/stores/call:
< 36. Doubling is an add. Charge **A = 64** each. The same bound covers the
two scalar label additions modulo r.

## Fp12 product: G <= 42,752

Instrumented count (`rho_ref.f12mul`, asserted in `experiment.py`):
**54 Fp products and 280 Fp additive operations**, including the xi-
multiplications (9a = 8a + a by three doublings and one add, 10 additive
operations per xi product) and the v-shift. Prices include each operation's
own loads/stores and call. A further **4,096** covers all remaining
straight-line work: Fp2/Fp6 call frames (21 frames), temporaries addressing,
copying the 48-word result.

    G = 54*384 + 280*64 + 4,096 = 20,736 + 17,920 + 4,096 = 42,752.

## Per-step work beyond G: E <= 2,048

| Item | Bound |
|---|---:|
| load w, two MUL + shift for j and k, table and stack addresses | 16 |
| two label additions mod r (2A) | 128 |
| stack comparisons, amortized (below): 2 x 256 | 512 |
| push: copy 48 + 8 label words (load+store), length update, cap check | 256 |
| step counter (two-word increment and compare against T), loop branch | 32 |
| **Sum** | **944** |

A comparison is lexicographic over the 48 Montgomery words: per word two
loads, CMP, conditional branch and index update (5), plus setup: <= 256 even
when all words are compared (the equality case). **Amortization is exact,
not probabilistic:** each step pushes at most one entry, each pop removes a
previously pushed entry, and each step's pop loop ends with exactly one
non-popping comparison, so total comparisons <= 2(T+1) and total pushes
<= T+1 over the whole run. Charge E = 2,048.

**Per-step bound S = 49,152 >= G + E = 44,800.** Table entries
(M_j as 48 words, c_j and d_j as 8 words: 448 bytes) are read inside the
Fp12 product's operand loads.

## Setup, recovery, output: F = 2^40

* Input parsing/validation, Montgomery and tower conversion: < 2^20.
* Subgroup check z^r = 1: <= 508 Fp12 products < 2^25.
* Coins: 2 salts and 2,050 scalars, each <= 1,024 draws of 256 bits,
  64 instructions per bit: < 2^35.1 (worst case).
* Table: 1,025 entries, each two exponentiations (<= 508 products each)
  plus one product: < 1,025 * 1,017 * 2^15.4 < 2^35.4.
* Stack headers: 256 length words; entries are written only when pushed.
* Recovery: inverse modulo r by exponentiation (<= 508 M), one multiply:
  < 2^19. Verification g^x = z: <= 508 Fp12 products < 2^25.
* Output [x]Q0 on the twist: affine double-and-add, <= 508 point operations,
  each one Fp2 inversion (one Fp exponentiation, < 2^19) plus < 2^14:
  < 2^29. Optional on-curve / [r]Q check: same again. Serialization < 2^20.
* Fixed program overhead (loading, allocation, shutdown): < 2^32.

Sum < 2^37; **F = 2^40** leaves slack.

## Total and score

    T = 170221906881474042593545034831872327888   (< 2^128)
    R_total = T * 49,152 + 2^40
            = 8366747167038212141557925552057288171978752
    W = R_total / 20,000
    log2 W = 128.29793444496976...  ->  submitted 128.297935

## Memory

Stacks: 256 x 512 entries x 464 bytes = 60,817,408 bytes worst case
(expected depth < 90 per stack, i.e. < 11 MB). Table 458,752 bytes. Code,
constants, scratch, stack frames < 2 MiB. Peak < 2^26 bytes. Advice: 0.
No address exceeds a 26-bit offset into one allocation.

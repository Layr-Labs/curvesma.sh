# Scalar AArch64 instruction accounting

Analytical bounds on executed scalar 64-bit instructions, on the benchmark's basis
(scalar AArch64 instruction equivalents; no 128-bit unit-cost arithmetic, no vector
instructions, no timing). One processor. `kernel.c` is an executable instance of the
same tower and is a consistency check, never the price.

## Fp multiplication — 640 instructions

Four little-endian 64-bit limbs in Montgomery form; Montgomery CIOS. A 64x64→128
product is one MUL plus one UMULH; sums and carries are ADDS/ADC; selection is
CSEL/AND/ORR.

| cell | count | each | bound |
|---|---:|---:|---:|
| product terms (16 schoolbook + 16 reduction), each MUL+UMULH and its carry combine | 32 | 8 | 256 |
| accumulator additions and carry propagation across the five-word accumulator | 32 | 6 | 192 |
| final four-word conditional subtraction and selection | 1 | 64 | 64 |
| operand loads, result stores, address formation | 1 | 128 | 128 |
| **sum** | | | **640** |

Cell-by-cell this totals 350; it is charged at 640 (1.8x) so that no optimism hides
in the cells.

## Fp12 multiplication — 58,824 instructions

Fp2 Karatsuba (3 Fp mults), Fp6 Karatsuba over Fp2 (6 Fp2 mults), Fp12 Karatsuba
over Fp6 (3 Fp6 mults): 18 Fp2 = **54 Fp multiplications**. The tower combinations
call Fp add/sub/neg **287** times, counted in `kernel.c`:

| layer | Fp add/sub/neg calls |
|---|---:|
| the 18 Fp2 products (2 adds + 2 subs each) | 90 |
| three Fp6 combination layers, including the two multiplication-by-xi steps | 156 |
| Fp12 combination layer, including multiplication by xi in `mul_by_v` | 41 |
| **total** | **287** |

| item | count | each | bound |
|---|---:|---:|---:|
| Fp multiplications | 54 | 640 | 34,560 |
| Fp add / sub / neg | 287 | 56 | 16,072 |
| loads, stores, spills, call/return, loop control | — | — | 8,192 |
| **one dense Fp12 multiplication** | | | **58,824** |

## Per-step additions

The partition and distinguished-point windows are read from the **canonical** form
of the same Fp coefficient, so one Montgomery reduction is charged per step.

| item | bound |
|---|---:|
| Montgomery reduction of one coefficient to canonical form (one Fp multiply plus fix-up) | 736 |
| partition: load, shift 120, mask 31, scale by 32 | 24 |
| label update a += c_j mod r (four ADDS/ADC, conditional subtract, store) | 64 |
| label update b += d_j mod r | 64 |
| DP test: compare the 128-bit window against D, branch | 8 |
| loop counter, table pointer, branch, spills | 40 |
| **per-step additions** | **936** |

    per step = 58,824 + 936 = 59,760 instructions = 2.9880 units

## Consistency check with the measured kernel

`kernel.c` (shipped; hash in `kernel-counts.json`) is a portable C instance of the
same tower, verified end to end by mode `g` (`g^r = 1`, `g != 1`):

| loop | callgrind Ir |
|---|---:|
| 4,096 chained multiplications | 214,427,735 |
| 8,192 chained multiplications | 428,701,783 |
| difference / 4,096 | **52,313** |

Toolchain gcc 14.2.0 `-O3 -fno-tree-vectorize -fno-builtin`; counter valgrind
3.24.0 callgrind, event `Ir`; x86-64 scalar integer code. The measurement is below
the analytic bound of 58,824, the expected direction for an upper bound.

## One-time work — 2^40

| item | bound |
|---|---:|
| 32 multipliers: 64 exponentiations of ≤ 508 multiplications | 2^30 |
| 12,700 fresh-segment exponentiation pairs: 25,400 exponentiations of ≤ 508 multiplications, at 59,760 instructions each | 2^39.49 |
| DP handling: ≤ 12,700 fingerprints, serialisation, table insert/compare | 2^32 |
| recovery: inversion mod r, one Fp12 exponentiation for `g^x == z`, compare, output | 2^22 |
| randomness, seed expansion, buffer setup | 2^20 |
| **total charged** | **2^40** |

## Totals

    T      = ceil(1.6 * floor(sqrt(r))) = 236,714,811,010,862,895,392,714,499,903,555,895,296
    R_run  = T * 59,760 + 2^40          = 2^143.343 instruction equivalents
    W_run  = R_run / 20,000             = 2^129.055598 units
    p_open = 0.39 / 0.66                = 0.5909091
    E[W]   = p_open * W_run + (1-p_open)= 2^128.296605831 units
    memory = 12,700 * 80 + 32 * 448 + 4,096 = 1,034,432 bytes < 2^20
    success= p_open * 0.66              = 0.39

## Arithmetic checks attached to the counted implementation

- Mode `g`: profile generator converted into the tower; `g^r = 1`, `g != 1`.
- Mode `s`, 2,000,000 production steps on the real GT from `g^a z^b` with one
  canonical Montgomery reduction per step: partition chi-square 18.51 on 31 degrees
  of freedom, 1,995 DP hits against 1,953 expected at a 2^-10 threshold.
- `evidence.py`: p and r prime; p ≡ 1 mod 6; r | p^12−1; ord_r(p) = 12; g ≠ 1;
  g^r = 1; tower product equals the independent w-basis product; Frobenius orbit
  exactly 12; four genuine colliding pairs recover the known logarithm and verify.

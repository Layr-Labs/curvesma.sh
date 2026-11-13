# Scalar AArch64 instruction accounting

Analytical bounds on executed scalar 64-bit instructions, on the benchmark's basis
(scalar AArch64 instruction equivalents; no 128-bit unit-cost arithmetic, no vector
instructions, no timing). One processor. `kernel.c` is an executable instance of
the same tower and is used as a consistency check, not as the price.

## Fp multiplication — 640 instructions

Residues are four little-endian 64-bit limbs in Montgomery form. `kernel.c`'s
`fmul` is Montgomery CIOS with a four-word carry chain. A 64x64→128 product is one
MUL plus one UMULH; sums and carries are ADDS/ADC; comparisons and selection are
SUBS/SBC and CSEL.

| cell | count | each | bound |
|---|---:|---:|---:|
| product terms (16 schoolbook + 16 reduction), each MUL+UMULH and its carry combine | 32 | 8 | 256 |
| accumulator additions and carry propagation across the five-word accumulator | 32 | 6 | 192 |
| final four-word conditional subtraction and selection | 1 | 64 | 64 |
| operand loads, result stores, address formation | 1 | 128 | 128 |
| **sum** | | | **640** |

Cell-by-cell this totals 350; it is charged at 640 (1.8x) so that no optimism hides
in the cells.

## Fp12 multiplication — 54,512 instructions

Fp2 uses Karatsuba (3 Fp mults), Fp6 uses Karatsuba over Fp2 (6 Fp2 mults), Fp12
uses Karatsuba over Fp6 (3 Fp6 mults): 18 Fp2 = **54 Fp multiplications**. The
tower combinations need about **210** Fp add/sub/neg operations: 72 in the Fp2
products, 120 in the Fp6 products (including the two multiplication-by-xi steps,
which are additions only, since xi = 9+u), and about 18 at the Fp12 level.

| item | count | each | bound |
|---|---:|---:|---:|
| Fp multiplications | 54 | 640 | 34,560 |
| Fp add / sub / neg (carry fix-up or conditional add of p) | 210 | 56 | 11,760 |
| loads, stores, spills, call/return, loop control of the whole multiplication | — | — | 8,192 |
| **one dense Fp12 multiplication** | | | **54,512** |

Per walk step adds:

| item | bound |
|---|---:|
| partition: load the first coefficient, shift 120, mask 31, scale by 32 | 24 |
| label update a += c_j mod r (four ADDS + four ADC, conditional subtract, store) | 64 |
| label update b += d_j mod r | 64 |
| DP test: compare a 128-bit window against D, branch | 8 |
| loop counter, table record pointer, branch, spills | 40 |
| **overhead** | **200** |

    per step = 54,512 + 200 = 54,712 instructions = 2.7356 units

## Consistency check with the measured kernel

`kernel.c` (shipped, hash in `kernel-counts.json`) is a portable C instance of the
same tower, verified end to end by mode `g` (`g^r = 1` and `g != 1` in Fp12 for the
profile generator):

| loop | callgrind Ir |
|---|---:|
| 4,096 chained multiplications | 214,427,735 |
| 8,192 chained multiplications | 428,701,783 |
| difference / 4,096 | **52,313** |

Toolchain gcc 14.2.0 `-O3 -fno-tree-vectorize -fno-builtin`; counter valgrind
3.24.0 callgrind, event `Ir`; x86-64 scalar integer code. The measurement is
**below** the analytic bound of 54,512, which is the expected direction for an
upper bound and is the only role the measurement plays. No cross-architecture
price transfer is used.

## One-time work — 2^36

| item | bound |
|---|---:|
| 32 multipliers: 64 exponentiations of <= 508 multiplications, plus 32 records | 2^30 |
| DP handling: <= 12,700 fingerprints (128-bit SHA-256 projection), serialisation, table insert/compare, spills | 2^32 |
| recovery: inversion mod r by exponentiation, one Fp12 exponentiation for the g^x == z check, twelve-coefficient compare, output | 2^22 |
| **total** | **2^36** |

## Totals

    T      = ceil(1.6 * floor(sqrt(r))) = 236,714,811,010,862,895,392,714,499,903,555,895,296
    R_run  = T * 54,712 + 2^36          = 2^143.216 instruction equivalents
    W_run  = R_run / 20,000             = 2^128.928275 units
    p_open = 0.39 / 0.66                = 0.5909091
    E[W]   = p_open * W_run + (1-p_open)= 2^128.169282977 units
    memory = 12,700 * 80 + 32 * 448 + 4,096 = 1,034,432 bytes < 2^20
    success= p_open * 0.66              = 0.39

## Arithmetic checks attached to the counted implementation

- `kernel.c` mode `g`: converts the profile generator into the tower and checks
  g^r = 1 and g != 1 in Fp12.
- `kernel.c` mode `s` (2,000,000 production walk steps on the real GT): partition
  histogram chi-square 19.45 on 31 degrees of freedom, 1,898 DP hits against 1,953
  expected at a 2^-10 threshold (`gt-stats.json`).
- `evidence.py`: p and r prime; p ≡ 1 mod 6; r | p^12-1; ord_r(p) = 12; g != 1;
  g^r = 1; the tower product equals an independent schoolbook product in the
  profile's w-basis; the Frobenius orbit is exactly 12; four genuine colliding
  state pairs on the real GT recover the known logarithm and verify.

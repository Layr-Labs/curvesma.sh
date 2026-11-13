# Scalar instruction accounting

Analytical bounds on executed scalar 64-bit instructions, in the same basis as
the benchmark's reference costs (scalar AArch64 instruction equivalents, no
128-bit unit-cost arithmetic, no vector instructions). One processor. The
implementation is a scalar loop over Montgomery-form Fp12 coefficients; the
shipped `kernel.c` is an executable instance of it.

## Fp multiplication

Residues are four little-endian 64-bit limbs, stored in Montgomery form.
`kernel.c`'s `fmul` is Montgomery CIOS with a four-word carry chain; a 64x64→128
product is one MUL-class instruction and the carry chains are ADDS/ADC and
SUBS/SBC.

Per multiplication: 16 schoolbook product terms and 16 reduction product terms
(32 multiplies), 5-word accumulator (16 additions), four conditional-subtract
cells and one selection of four words. A deliberately loose bound:

| work | count | each | bound |
|---|---:|---:|---:|
| product and reduction terms | 32 | 8 | 256 |
| accumulator additions and carries | 24 | 8 | 192 |
| conditional subtract and select | 1 | 80 | 80 |
| loads, stores, addressing, loop | 1 | 96 | 96 |
| **per Fp multiplication** | | | **624** |

The measured portable C kernel executes 52,313 / 54 = **969** instructions per Fp
multiplication including the whole tower (see below), so 624 is not used for the
final charge in any case: the ledger charges the reference-measured dense Fp12
multiplication directly, exactly once, rather than summing Fp cells and risking
double counting.

## Fp12 multiplication

Fp2 uses Karatsuba (3 Fp multiplications); Fp6 uses Karatsuba over Fp2 (6 Fp2
multiplications); Fp12 = Fp6[w]/(w^2-v) uses Karatsuba over Fp6 (3 Fp6
multiplications). Hence one dense Fp12 multiplication is 18 Fp2 = **54 Fp
multiplications** plus 24 Fp2 additions, 12 Fp6 additions and two multiplications
by xi = 9+u (six Fp doublings and two Fp subtractions each, additions only).

The ledger charges this operation at **22,996** instruction equivalents, the
benchmark's own reference measurement of a dense Fp12 multiplication
(`tracks/docs/COST_MODEL.md`, `blst_fp12_mul`, AArch64, measured at 4,096 and 8,192
iterations and differenced). This is the operation the challenge's unit C = 20,000
is derived from. It was measured on BLS12-381 (six 64-bit limbs per Fp element);
BN254 needs four limbs, so the charge is conservative for this field. It is
charged once per step; the internal Fp multiplications are not charged again.

## Measured instance (evidence, not the charge)

`kernel.c` is a self-contained portable C instance of the same tower. It is
verified end to end: mode `g` converts the profile's generator to the tower and
checks `g^r == 1` and `g != 1` in Fp12. Fully counted:

| loop | callgrind Ir |
|---|---:|
| 4,096 chained multiplications | 214,427,735 |
| 8,192 chained multiplications | 428,701,783 |
| difference / 4,096 | **52,313** |

Toolchain: gcc 14.2.0, `-O3 -fno-tree-vectorize -fno-builtin`; counter:
valgrind 3.24.0 callgrind, event `Ir`; host ISA x86-64 with scalar integer code
only. `kernel-counts.json` records the artifact hash, both totals and the
derivation. This C figure is **not** used for the submitted score: charging it
would give `log2 E[W] = 128.067` (1.18 bit worse). It is shipped so that the
tower, the conversion and the multiplication are independently executable and so
that the difference between a portable C implementation and the reference-priced
implementation is visible rather than hidden.

## Per-step overhead beyond the multiplication

| item | bound |
|---|---:|
| partition: load first coefficient, shift 249, mask 31, scale by 32 | 24 |
| label update a += c_j mod r: 4 ADDS + 4 ADC, conditional subtract, store | 64 |
| label update b += d_j mod r | 64 |
| distinguished-point test: compare two words against D, branch | 8 |
| loop counter, table record pointer, branch, spills | 40 |
| **overhead** | **200** |

     per step = 22,996 + 200 = 23,196 instruction equivalents = 1.15980 units.

## One-time work

| item | bound |
|---|---:|
| 32 multipliers M_j = g^{c_j} z^{d_j}: 64 exponentiations, ~508 multiplications each, plus 32 stores | 2^30 |
| distinguished-point handling: <= 12,700 fingerprints (128-bit projection), serialisation, table insert/compare, spills | 2^32 |
| recovery: one inversion mod r by exponentiation (<= 508 multiplications mod r), one Fp12 exponentiation g^x, twelve-coefficient compare, output | 2^22 |
| **total one-time** | **2^36** |

The one-time bound is dominated by the 2^32 allowance for DP bookkeeping and is
five orders of magnitude below the walk.

## Totals

    T          = ceil(1.6 * floor(sqrt(r))) = 236,714,811,010,862,895,392,714,499,903,555,895,296
    R_run      = T * 23,196 + 2^36           = 2^141.98 instruction equivalents
    W_run      = R_run / 20,000              = 2^127.690294 units
    p_open     = 0.39 / 0.68                 = 0.5735294
    E[W]       = p_open * W_run + (1-p_open) = 2^126.888233004 units
    memory     = 12,700 * 80 + 32 * 448 + 4,096 = 1,034,432 bytes < 2^20
    success    = p_open * 0.68               = 0.39

## Correctness of the counted arithmetic

- `evidence.py` (default part) verifies: p and r prime; p ≡ 1 mod 6;
  r | p^12 - 1; ord_r(p) = 12; g != 1; g^r = 1; the tower product equals an
  independent schoolbook product in the profile's w-basis for three deterministic
  test pairs; the Frobenius orbit of a group element has size exactly 12.
- `kernel.c` mode `g` performs the same g^r = 1 check natively in the counted
  implementation, so the instruction count is attached to arithmetic that is
  known to compute the intended function.
- `evidence.py` runs the production walk on the production group for 200 steps
  and re-derives g^{a+bx} independently every 50 steps: 0 bad steps.

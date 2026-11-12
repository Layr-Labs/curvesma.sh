# Implementation and instruction measurement

## Files

- `fp12.h`, `fp12.c`: scalar 64-bit Montgomery `Fp`, tower `Fp2/Fp6/Fp12`
  (Karatsuba `Fp2`; `xi = 9+u`; `w_t^2 = v`). Same arithmetic as our
  `bn254-gt-dlp` package `rho-tower-v1`.
- `rho_pi.c`: profile-generator conversion; the measured 32-adding GT walk
  step; G1 Jacobian doubling, mixed addition with the affine base `(1,2)`,
  double-and-add scalar multiplication and affine output; `selftest`,
  `point <scalar>` and `measure K` modes.
- `vectors.json`: C outputs for six scalars, recomputed independently in Python
  by `experiment.py`.
- `experiment.py`: standard-library recomputation of the ledger, the capped
  sampler and its failure bound, the mixer cross-check, the G1 vectors
  (regenerated in-package) and the production-walk evidence.
- `LOOP_BOUNDS.txt`: every loop of the algorithm with its hard bound.
- `evidence/worstcase-sweep.txt`: all 704 measured configurations.
- `evidence/production-walk.txt`: frozen 64-trial experiment output.
- `evidence/mixer-crosscheck.txt`: C `hashdump` output and the Python match.
- `evidence/callgrind-summary.txt`: historical annotated Callgrind totals for
  the same step loop at K = 10,000.

## Build and self-test

```
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho_pi fp12.c rho_pi.c
./rho_pi selftest
```

Output on the build machine:

```
g_is_one=0 g^r_is_one=1
assoc=1
P0_on_curve=1
r_times_P0_is_inf=1
2P0.x=030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd3
2P0.y=15ed738c0e0a7c92e7845f96b2ae9c0a68a6a449e3538fc7ff3ebf7a5a18a2c4
```

`g^r = 1` checks the profile-basis conversion and the tower multiplication
together. `P0` is the profile G1 generator; `2P0` matches an independent affine
computation on `y^2 = x^3 + 3`, and `[r]P0 = O` checks the group order.

## Worst-case instruction count: GT walk (not an average)

Machine: AArch64, Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0, Valgrind
3.22.0 (`--tool=callgrind`, event `Ir`).

Why a sweep is a worst-case bound: every loop inside a step has a
data-independent trip count (`fp_mul` 4x4 limbs, the mixer 24 rounds over 48
words, two 4-limb scalar additions), so the executed count depends on the
inputs only through the arms of the data-dependent decisions -- the
carry/reduction arms of `fp_add`, the borrow arm of `fp_sub`, the final
comparison in `fp_mul`, the select in `scalar_add_r`. Arms change *which*
instructions run, never how many iterations run. The sweep therefore covers the
whole input dependence:

* table index `j` = 0..31 (all 32 production table entries; the audit mode
  `./rho_pi worstcase K j pat` still computes the mixer and only overrides its
  result, so the shipped step is measured);
* 22 operand patterns: 0 = production start `g`, 1 = all limbs `0xffff...`
  (carry arm everywhere), 2 = all limbs 1 (borrow arm everywhere),
  3 = all limbs `FP_ONE_MONT`, 4..19 = deterministic pseudo-random 384-byte
  elements, 20 = all limbs `p_limb - 1` (carry + full comparison chain),
  21 = all limbs `p_limb` (every addition overflows the modulus);
* K = 1,000 and K = 3,000, differenced so the table build and conversion cancel.

Reproduce:

```
for j in $(seq 0 31); do
  for pat in $(seq 0 21); do
    for K in 1000 3000; do
      valgrind --tool=callgrind --callgrind-out-file=cg.$K.$j.$pat.out \
        ./rho_pi worstcase $K $j $pat
    done
  done
done
```

Result over the 704 configurations (1,408 runs):

| statistic | Ir per step |
|---|---:|
| minimum | 46,320.00 |
| median | 47,249.14 |
| **maximum** (j = 13, pat = 12) | **47,265.24** |

**Charged: 47,266 instructions per step** -- the maximum rounded up to the next
instruction, no percentage padding.

The measured step is exactly the loop body of `rho_step`: one `fp12_mul`, one
48-word mixer, two 4-word additions mod r and the table load. Floyd comparison
and loop control are charge separately at 4,096 instructions per iteration, an
analytic over-charge (48 word compares at ~3 instructions each plus a counter
and a branch).

## Capped sampler cost

`proof.md` section 2 defines the sampler: one SHA-256 call per 256-bit word.
Worst case 66 scalars x K = 64 attempts = 4,224 attempts at 2,048 instructions
each = 8,650,752 instruction equivalents, charged as its own ledger row at the
cap, not at an observed attempt count.

## G1 vector regeneration

`experiment.py` regenerates the six vectors from `y^2 = x^3 + 3` and `P0 =
(1,2)` and compares them with `vectors.json` when present, so a clean checkout
without that file still runs; the shipped file matches.

## Instruction counts: G1 scalar multiplication

```
valgrind --tool=callgrind --callgrind-out-file=cg.point.out \
  ./rho_pi point 8000000000000000000000000000000000000000000000000000000000000000
```

`PROGRAM TOTALS`: **1,062,076 Ir** for the whole process, including startup,
the profile-generator conversion and output formatting. The double-and-add,
one Fermat inversion and affine conversion are therefore below 2^21
instructions; the ledger charges 2^24, over 15x the measured process.

## Memory

32 GT records of 448 bytes (384-byte element + 32-byte `c_j` + 32-byte `d_j`),
two GT walk states, 48-word hash scratch, G1 Jacobian temporaries and code:
below 2^20 bytes. The sampler stores no words (they are consumed as produced).
No table is stored beyond the run, no advice is embedded, and the seed, salt and
increments are drawn fresh and charged in the ledger.

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
- `evidence/step-measurement.txt`: per-trajectory step costs on the production table.
- `evidence/percall-arms.txt`: per-call arm measurements used by the composition.
- `sha256.c`, `sha256.h`: SHA-256 for the capped sampler (FIPS 180-4).
- `percall.c`: the per-call arm microbenchmark.
- `evidence/production-walk.txt`: frozen 64-trial experiment output.
- `evidence/mixer-crosscheck.txt`: C `hashdump` output and the Python match.
- `evidence/callgrind-summary.txt`: historical annotated Callgrind totals for
  the same step loop at K = 10,000.

## Build and self-test

```
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho_pi fp12.c sha256.c rho_pi.c
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

## Worst-case instruction count: composition, not an average

Machine: AArch64, Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0, Valgrind
3.22.0 (`--tool=callgrind`, event `Ir`).

The bound is composed from three facts, all recorded in the package:

1. **Fixed call counts.** One step executes one `fp12_mul`, one mixer, two
   `scalar_add_r` and one table load. Expanding the tower from the source gives
   81 `fp_mul`, 164 `fp_add` and 103 `fp_sub` per step. Call counts do
   not depend on the data.
2. **Fixed loop trip counts, measured arms.** Every loop inside those functions
   has a data-independent trip count, so a call costs its straight-line code
   plus the arms of its data-dependent decisions. `percall.c` measures each arm
   (`evidence/percall-arms.txt`): `fp_add` 105.0 vs 89.0 instructions, `fp_sub`
   92.0 vs 60.0, `fp_mul` 282.0 vs 282.0.
3. **A measured baseline on the algorithm's own table.** `./rho_pi mktable
   20261112 table.bin` runs the capped SHA-256 sampler of proof.md section 2 and
   builds `M_j = g^c_j z^d_j`; `./rho_pi stepwc K j table.bin` then runs the
   shipped step on that table. Over the 33 trajectories (j = 0..31 and the
   unpinned production mixer) the mean step cost is **47,411.27 Ir**
   (`evidence/step-measurement.txt`).

The composition is then

    cost_per_step <= mean_per_step + sum_f c_f (max_f - min_f)
                  = 47,411.27 + 81x8 + 164x16 + 103x32 = 53,979.27 Ir,

charged as **C_STEP = 53980 instructions per step** (round up). No percentage
padding, and no part of the charge is an average over unmeasured steps: the
measured mean only supplies the straight-line term, which is common to every
step, and the arm deltas supply the variation term multiplied by the fixed call
counts.

Reproduce:

```
for j in $(seq 0 32); do
  for K in 1000 3000; do
    valgrind --tool=callgrind --callgrind-out-file=cg.$K.$j.out \
      ./rho_pi stepwc $K $j /tmp/table.bin
  done
done
```

The Floyd comparison and loop control are outside the measured loop and are
charged separately at 4,096 instructions per iteration, an analytic over-charge
(48 word compares at ~3 instructions each plus a counter and a branch).

## Capped sampler cost

One SHA-256 call per 256-bit word. Worst case 66 scalars x K = 64 attempts =
4,224 attempts at 2,048 instructions each = 8,650,752, charged as its own ledger
row at the cap. The production table build (64 exponentiations of 254-bit
exponents, one `g^x`, one `Y0`) is charged inside the 2^36 setup row; at the
measured 282 instructions per `fp_mul` and 81 `fp_mul` per `fp12_mul`, 64 x 508
multiplications is 3.3e7 instructions, far below 2^36.

## Memory

32 GT records of 448 bytes (384-byte element + 32-byte `c_j` + 32-byte `d_j`),
two GT walk states, 48-word hash scratch, G1 Jacobian temporaries and code:
below 2^20 bytes. The sampler stores no words (they are consumed as produced).
No table is stored beyond the run, no advice is embedded, and the seed, salt and
increments are drawn fresh and charged in the ledger.

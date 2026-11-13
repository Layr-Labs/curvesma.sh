# Implementation and instruction measurement

## Files

- `fp12.h`, `fp12.c`: scalar 64-bit Montgomery `Fp`; branchless inline `fp_add`
  and `fp_sub`; `Fp2` Karatsuba; six-multiplication Karatsuba `Fp6`; `Fp12`
  Karatsuba over `w^2 = v`.
- `rho.c`: profile-generator conversion, 32-adding walk, mixer, scalar
  arithmetic mod r, `selftest` and a `measure K` mode. The profile-basis map and
  the step are unchanged from our `rho-tower-v1`; the arithmetic underneath is
  the optimised one.
- `crosscheck.c`: separate correctness test of the optimised arithmetic against
  a schoolbook `Fp6`/`Fp12` reference built from the same `fp2_mul` primitive.
- `experiment.py`: standard-library recomputation of the ledger, the capped
  sampler and its failure bound, the mixer cross-check and the production-walk
  evidence.
- `LOOP_BOUNDS.txt`: every loop of the algorithm with its hard bound.
- `evidence/step-measurement.txt`: per-trajectory step costs on the production table.
- `evidence/percall-arms.txt`: per-call arm measurements used by the composition.
- `sha256.c`, `sha256.h`: SHA-256 for the capped sampler (FIPS 180-4).
- `percall.c`: the per-call arm microbenchmark.
- `evidence/production-walk.txt`: frozen 64-trial experiment output.
- `evidence/mixer-crosscheck.txt`: C `hashdump` output and the Python match.
- `evidence/callgrind-summary.txt`: historical annotated Callgrind totals.

## Build, cross-check and self-test

```
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho fp12.c sha256.c rho.c
./rho selftest
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o crosscheck crosscheck.c fp12.c
./crosscheck
```

Output on the build machine:

```
g_is_one=0 g^r_is_one=1
assoc=1
fp6/fp12 vs schoolbook: 2000 pairs ok
add/sub round-trip: 2000 ok
```

`g^r = 1` checks the profile-basis conversion and the tower multiplication
together. `crosscheck` compares every product of 2,000 deterministic random
`Fp6` and `Fp12` pairs with the 9-multiplication schoolbook form and checks
2,000 `fp_add`/`fp_sub` round trips.

## Worst-case instruction count: composition, not an average

Machine: AArch64, Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0, Valgrind
3.22.0 (`--tool=callgrind`, event `Ir`).

The bound is composed from three facts, all recorded in the package:

1. **Fixed call counts.** One step executes one `fp12_mul`, one mixer, two
   `scalar_add_r` and one table load. Expanding the tower from the source gives
   81 `fp_mul`, 0 `fp_add` and 0 `fp_sub` per step (all inline and branchless, so their arm delta is 0). Call counts do
   not depend on the data.
2. **Fixed loop trip counts, measured arms.** Every loop inside those functions
   has a data-independent trip count, so a call costs its straight-line code
   plus the arms of its data-dependent decisions. `percall.c` measures each arm
   (`evidence/percall-arms.txt`): `fp_add` 105.0 vs 89.0 instructions, `fp_sub`
   92.0 vs 60.0, `fp_mul` 282.0 vs 282.0.
3. **A measured baseline on the algorithm's own table.** `./rho mktable
   20261112 table.bin` runs the capped SHA-256 sampler of proof.md section 2 and
   builds `M_j = g^c_j z^d_j`; `./rho stepwc K j table.bin` then runs the
   shipped step on that table. Over the 33 trajectories (j = 0..31 and the
   unpinned production mixer) the mean step cost is **31,589.48 Ir**
   (`evidence/step-measurement.txt`).

The composition is then

    cost_per_step <= mean_per_step + sum_f c_f (max_f - min_f)
                  = 31,589.48 + 81x8 + 0x16 + 0x32 = 32,237.48 Ir,

charged as **C_STEP = 32238 instructions per step** (round up). No percentage
padding, and no part of the charge is an average over unmeasured steps: the
measured mean only supplies the straight-line term, which is common to every
step, and the arm deltas supply the variation term multiplied by the fixed call
counts.

Reproduce:

```
for j in $(seq 0 32); do
  for K in 1000 3000; do
    valgrind --tool=callgrind --callgrind-out-file=cg.$K.$j.out \
      ./rho stepwc $K $j /tmp/table.bin
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

32 records of 448 bytes (384-byte element + 32-byte `c_j` + 32-byte `d_j`),
two walk states, the hash scratch and code: below 2^20 bytes (the sampler
stores no words; they are consumed as produced). No table is
stored beyond the run, no advice is embedded, and the salt and increments are
drawn fresh and charged in setup.

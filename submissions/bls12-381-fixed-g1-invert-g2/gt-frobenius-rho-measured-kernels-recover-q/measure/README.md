# Measurement record: the blst baseline

The blst baseline of the two walk kernels is measured here, on the track's own cost-model
toolchain, by the cost model's own instrument. `experiment.py` reads `results.tsv` and
`callgrind/`, recomputes every number below, and exits with an assertion error if anything differs;
it is not a table transcribed by hand. The files were produced in one session of the harness below and
are shipped as they were written; they were not regenerated for this revision.

## What is measured

The routine is `blst_fp12_cyclotomic_sqr` of blst at the pinned revision f7a88067 (exported at
`bindings/blst.h:152`, implemented through `src/fp12_tower.c:741`). The harness calls it on **real
GT elements**: `cyc2.c` computes `a = e(P0, Q)` (Miller loop then final exponentiation) and then
calls the routine on that element in a loop, so the operand is in the cyclotomic subgroup where
`cyclotomic_sqr_fp12` is defined; it prints `blst_fp12_in_group(&a)` afterwards, which reads 1.
The same program and toolchain measure `blst_fp12_mul` (kind 1) and `blst_fp_mul` (`fpmul.c`),
which are the cost model's two reference counts, as a calibration of the instrument.

## Toolchain and instrument

- image `refcount-harness`: Debian 13.7, `cc (Debian 14.2.0-19) 14.2.0`, `valgrind-3.24.0`,
  built from `harness/Dockerfile`. That is exactly the toolchain `tracks/docs/COST_MODEL.md` names.
- flags for blst's C: `-O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only`; blst's
  `build/assembly.S` as shipped (object hashes are printed by the build script).
- instrument: `valgrind --tool=callgrind --collect-atstart=no`, with `CALLGRIND_ZERO_STATS` and
  `CALLGRIND_TOGGLE_COLLECT` around the call loop only, so the counts are the routine's work and
  nothing else. Cost model unit: 20,000 executed scalar AArch64 instructions.

Reproduce with (paths are the container's):

    docker run --rm -v "$PWD/harness:/harness" -v "$PWD/callgrind:/work" -v /path/to/blst:/blst \
        refcount-harness /harness/build-and-run.sh

`harness/build-and-run.sh` is the script that produced every record in `callgrind/` and every row
of `results.tsv`; `harness/cyc2.c` and `harness/fpmul.c` are the two C files it compiles, verbatim.

## Files

- `results.tsv` - the exact Callgrind `Ir` total of every measurement: three routines, at
  `n = 1` and `n = 65` contiguous calls, at five slice offsets `0, 3, 101, 1009, 4001` (the
  offset is the number of multiplication steps taken before collection starts, so the five slices
  sit at different points of a walk and on five different GT elements).
- `callgrind/` - six raw records: the `n = 1` and `n = 65` runs at offset 0 for each routine.
  Their `summary:` line is the `Ir_total` of the corresponding `results.tsv` row, and
  `experiment.py` checks that equality, so the table cannot drift from the raw record.

## The numbers

| routine | Ir at n=1 | Ir at n=65 | per call (I65-I1)/64 | cost model | miss |
| --- | ---: | ---: | ---: | ---: | ---: |
| `blst_fp_mul` | 403 | 24,915 | **383.0** | 383 | 0.00% |
| `blst_fp12_mul` | 23,017 | 1,494,761 | **22,996.0** | 22,996 | 0.00% |
| `blst_fp12_cyclotomic_sqr` | 9,823 | 637,215 | **9,803.0** | - | - |

The third line is the blst baseline. The counts are identical at all five slice
offsets.

This package does not charge the blst routine: the walk's squaring is the straight-line kernel
`kernels/sqr.s`, charged 6,997 + 283 = 7,280 instructions = 0.364 units (proof.md section 5), and the
measurement above is the baseline it replaces.

## Correctness of the routine as this walk's squaring

`cyc2.c` leaves the element in GT (`in_group` prints 1 after the loop), and GT lies in the
cyclotomic subgroup, where `blst_fp12_cyclotomic_sqr` computes the square. The routine's source
(`src/fp12_tower.c` lines 630-660) is the Granger-Scott map that proof.md section 2 (F3) writes
down, three Fp4 squarings followed by the 3t - 2u and 3t + 2u combinations. `experiment.py` part A
checks that map against y*y on random BLS12-381 GT elements in its own arithmetic, and part C
checks it inside the walk and its commutation with the Frobenius map, which is what makes the
branch a map on classes. That is a check of the formula; the walk runs `kernels/sqr.s`, which
`kernels/check_kernels.py` checks limb by limb against a big-integer reference.

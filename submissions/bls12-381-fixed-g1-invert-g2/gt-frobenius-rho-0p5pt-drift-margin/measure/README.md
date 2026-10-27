# Measurement record: the blst baseline

The blst baseline of the two walk kernels is measured here on the track's cost-model toolchain with
its own instrument. `experiment.py` reads `results.tsv` and `callgrind/`, recomputes every number
below, and asserts equality. The files are shipped verbatim; only the cost-side drift margin changed
from the parent lineage.

## What is measured

The routine is `blst_fp12_cyclotomic_sqr` of blst at the pinned revision f7a88067
(`src/fp12_tower.c`). `cyc2.c` computes `a = e(P0, Q)` and calls the routine on it in a loop, so the
operand is in the cyclotomic subgroup (`blst_fp12_in_group(&a)` prints 1). The same program and toolchain
measure `blst_fp12_mul` and `blst_fp_mul` (`fpmul.c`), the cost model's two reference counts.

## Toolchain and instrument

- image `refcount-harness`: Debian 13.7, `cc (Debian 14.2.0-19) 14.2.0`, `valgrind-3.24.0`, built
  from `harness/Dockerfile` -- exactly the toolchain `tracks/docs/COST_MODEL.md` names.
- blst C flags `-O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only`; blst's
  `build/assembly.S` as shipped.
- instrument `valgrind --tool=callgrind --collect-atstart=no`, with `CALLGRIND_ZERO_STATS` and
  `CALLGRIND_TOGGLE_COLLECT` around the call loop only. Cost model unit: 20,000 instructions.

`harness/build-and-run.sh` (run in `refcount-harness`) produced every record in `callgrind/` and row of
`results.tsv`; `harness/cyc2.c` and `harness/fpmul.c` are the C files it compiles.

## Files

`results.tsv` holds the Callgrind `Ir` total of every measurement (three routines, at `n = 1` and
`n = 65` contiguous calls, at five slice offsets `0, 3, 101, 1009, 4001`). `callgrind/` holds six raw
records (the `n = 1` and `n = 65` runs at offset 0 per routine); each `summary:` line equals the
`Ir_total` of the matching `results.tsv` row, and `experiment.py` asserts that equality.

## The numbers

| routine | Ir at n=1 | Ir at n=65 | per call (I65-I1)/64 | cost model | miss |
| --- | ---: | ---: | ---: | ---: | ---: |
| `blst_fp_mul` | 403 | 24,915 | **383.0** | 383 | 0.00% |
| `blst_fp12_mul` | 23,017 | 1,494,761 | **22,996.0** | 22,996 | 0.00% |
| `blst_fp12_cyclotomic_sqr` | 9,823 | 637,215 | **9,803.0** | - | - |

The third line is the blst baseline, identical at all five slice offsets. This package does not charge
the blst routine: the walk's squaring is the straight-line kernel `kernels/sqr.s`, charged 7,069
instructions = 0.35345 units (proof.md section 5); the measurement above is the baseline it replaces, and
`kernels/check_kernels.py` checks the kernel limb by limb against a big-integer reference.

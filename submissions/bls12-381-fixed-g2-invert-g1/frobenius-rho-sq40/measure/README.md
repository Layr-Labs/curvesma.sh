# Measurement record: the squaring step's instruction price

Everything the ledger charges for a squaring step is measured here, on the track's own cost-model
toolchain, by the cost model's own instrument. `experiment.py` reads `results.tsv` and
`callgrind/`, recomputes every number below, and exits with an assertion error if anything differs;
it is not a table transcribed by hand.

## What is measured

The routine is `blst_fp12_cyclotomic_sqr` of blst at the pinned revision f7a88067 (exported at
`bindings/blst.h:152`, implemented through `src/fp12_tower.c:741`). The harness calls it on **real
GT elements**: `cyc2.c` computes `a = e(P0, Q)` (Miller loop then final exponentiation) and then
calls the routine on that element in a loop, so the operand is in the cyclotomic subgroup where
`cyclotomic_sqr_fp12` is defined; it prints `blst_fp12_in_group(&a)` afterwards, which reads 1.
The same program and toolchain measure `blst_fp12_mul` (kind 1) and `blst_fp_mul` (`fpmul.c`),
which are the cost model's two reference counts, as a calibration of the instrument.

## Toolchain and instrument

- image `sqprice-refcount`: Debian 13.7, `cc (Debian 14.2.0-19) 14.2.0`, `valgrind-3.24.0`,
  built from `harness/Dockerfile`. That is exactly the toolchain `tracks/docs/COST_MODEL.md` names.
- flags for blst's C: `-O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only`; blst's
  `build/assembly.S` as shipped (object hashes are printed by the build script and repeated in
  `disassembly-audit.txt`).
- instrument: `valgrind --tool=callgrind --collect-atstart=no`, with `CALLGRIND_ZERO_STATS` and
  `CALLGRIND_TOGGLE_COLLECT` around the call loop only, so the counts are the routine's work and
  nothing else. Cost model unit: 20,000 executed scalar AArch64 instructions.

Reproduce with (paths are the container's):

    docker run --rm -v "$PWD/harness:/harness" -v "$PWD/callgrind:/work" -v /path/to/blst:/blst \
        sqprice-refcount /harness/build-and-run.sh

`harness/build-and-run.sh` is the script that produced every record in `callgrind/` and every row
of `results.tsv`; `harness/cyc2.c` and `harness/fpmul.c` are the two C files it compiles, verbatim.

## Files

- `harness/cyc2.c`, `harness/fpmul.c`, `harness/build-and-run.sh`, `harness/Dockerfile` - the
  harness and the image.
- `results.tsv` - the exact Callgrind `Ir` total of every measurement: three routines, at
  `n = 1` and `n = 65` contiguous calls, at five slice offsets `0, 3, 101, 1009, 4001` (the
  offset is the number of multiplication steps taken before collection starts, so the five slices
  sit at different points of a walk and on five different GT elements).
- `callgrind/` - six raw records: the `n = 1` and `n = 65` runs at offset 0 for each routine.
  Their `summary:` line is the `Ir_total` of the corresponding `results.tsv` row, and
  `experiment.py` checks that equality, so the table cannot drift from the raw record.
- `disassembly-audit.txt` - the branch census, the call-tree composition, the static-vs-dynamic
  reconciliation, and the derivation of the charged 9,798.

## The numbers

| routine | Ir at n=1 | Ir at n=65 | per call (I65-I1)/64 | cost model | miss |
| --- | ---: | ---: | ---: | ---: | ---: |
| `blst_fp_mul` | 403 | 24,915 | **383.0** | 383 | 0.00% |
| `blst_fp12_mul` | 23,017 | 1,494,761 | **22,996.0** | 22,996 | 0.00% |
| `blst_fp12_cyclotomic_sqr` | 9,823 | 637,215 | **9,803.0** | - | - |

Both reference counts of the cost model reproduce to the instruction, which is the calibration of
the instrument; the third line is the new measurement. The counts are identical at all five slice
offsets and, by the branch census, cannot depend on the input at all: the routine's call tree is 18
functions and 1,017 static instructions with zero conditional, backward or indirect branches (the
16 indirect instructions are `ret`s), so it is straight-line.

The ledger charges the routine's body, 9,803 - 5 = **9,798**, plus 300 instructions of alpha/beta
doubling, shared rows and slack, i.e. a squaring step of 10,098 instructions = 0.50490 units of
20,000. The 5 instructions are the harness loop's own increment, compare, branch and call, which
the step ledger's common rows already charge (30 instructions of loop control), so charging 9,803
would double count them; the difference moves log2 W by 2.5e-4 bits either way and no printed
digit. See `disassembly-audit.txt` section 4.

## Correctness of the routine as this walk's squaring

`cyc2.c` leaves the element in GT (`in_group` prints 1 after the loop), and GT lies in the
cyclotomic subgroup, where `blst_fp12_cyclotomic_sqr` equals `blst_fp12_sqr`. `experiment.py` part
A checks that equality numerically on GT elements (and after 997 further squarings), and part A
also checks that cyclotomic squaring commutes with the Frobenius map, which is what makes the
branch a map on classes.

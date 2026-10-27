# Measurement record: the squaring step's instruction price

Everything the ledger charges for a squaring step is measured here, on the track's own cost-model
toolchain, by the cost model's own instrument. `experiment.py` reads `results.tsv` (and any raw
records shipped under `callgrind/`, checked against the table when present),
recomputes every number below, and exits with an assertion error if anything differs;
it is not a table transcribed by hand. The files were produced in one session of the harness below and
are shipped as they were written; they were not regenerated for this revision. This revision ships the
`results.tsv` table, the audit and this note; the six raw `callgrind/` records and the `harness/`
sources they were produced with are omitted to fit the submission byte cap (the parent
submission dc1293cb carries them verbatim, and part D checks them when present).

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
  `build/assembly.S` as shipped (object hashes are printed by the build script and repeated in
  `disassembly-audit.txt`).
- instrument: `valgrind --tool=callgrind --collect-atstart=no`, with `CALLGRIND_ZERO_STATS` and
  `CALLGRIND_TOGGLE_COLLECT` around the call loop only, so the counts are the routine's work and
  nothing else. Cost model unit: 20,000 executed scalar AArch64 instructions.

Reproduce with the harness sources (carried verbatim by the parent submission dc1293cb; omitted
here for the byte cap) with paths as the container's:

    docker run --rm -v "$PWD/harness:/harness" -v "$PWD/callgrind:/work" -v /path/to/blst:/blst \
        refcount-harness /harness/build-and-run.sh

`harness/build-and-run.sh` is the script that produced every record in `callgrind/` and every row
of `results.tsv`; `harness/cyc2.c` and `harness/fpmul.c` are the two C files it compiles, verbatim.

## Files

- `results.tsv` - the exact Callgrind `Ir` total of every measurement: three routines, at
  `n = 1` and `n = 65` contiguous calls, at five slice offsets `0, 3, 101, 1009, 4001` (the
  offset is the number of multiplication steps taken before collection starts, so the five slices
  sit at different points of a walk and on five different GT elements).
- `callgrind/` (parent only) - six raw records: the `n = 1` and `n = 65` runs at offset 0 for each routine.
  Their `summary:` line is the `Ir_total` of the corresponding `results.tsv` row, and
  `experiment.py` checks that equality when the records are shipped, so the table cannot drift from the raw record.
- `disassembly-audit.txt` - the branch census, the call-tree composition, the static-vs-dynamic
  reconciliation, the charged price 9,803, and the argument from the blst source that the routine
  returns fully reduced limbs (section 5; an argument from the source, not a run).

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

The earlier revision's ledger charged the whole marginal, **9,803** instructions, for the routine,
plus the 270 instructions of every-step work of that revision's proof.md section 5: a squaring step of 10,073
instructions = 0.50365 units of 20,000. This package replaces that analytic routine with the
straight-line kernel of `kernels/sqr.s` and charges 7,226 + 226 = 7,452 instructions = 0.3726 units
(proof.md section 5); the measurement above is the baseline it replaces. No deduction is made for the
harness loop's own control (about 5 instructions per iteration), which the every-step work also covers through its
19 instructions of loop control, so the squaring step carries that control twice; that can only
overstate the cost. See `disassembly-audit.txt` section 4.

## Correctness of the routine as this walk's squaring

`cyc2.c` leaves the element in GT (`in_group` prints 1 after the loop), and GT lies in the
cyclotomic subgroup, where `blst_fp12_cyclotomic_sqr` computes the square. The routine's source
(`src/fp12_tower.c` lines 630-660) is the Granger-Scott map that proof.md section 2 (F3) writes
down, three Fp4 squarings followed by the 3t - 2u and 3t + 2u combinations. `experiment.py` part A
checks that map against y*y on random BLS12-381 GT elements in its own arithmetic, and part C
checks it inside the walk and its commutation with the Frobenius map, which is what makes the
branch a map on classes. That is a check of the formula the routine implements; the routine
itself is priced here and was not executed by `experiment.py`.

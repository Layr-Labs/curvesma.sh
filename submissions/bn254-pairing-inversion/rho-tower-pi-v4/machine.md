# Implementation, constant-time step and static instruction count

## Files

- `fp12.h`, `fp12.c`: scalar 64-bit Montgomery `Fp`; branchless inline
  `fp_add`/`fp_sub`; branchless `fp_mul` final reduction (mask select);
  `Fp2` Karatsuba; six-multiplication Karatsuba `Fp6`; `Fp12` Karatsuba over
  `w^2 = v`. Every function in the step's call graph is straight-line.
- `rho_pi.c`: profile-generator conversion, the shipped GT walk step `rho_step`,
  mixer, `scalar_add_r` (mask select), capped sampler, the G1 part (Jacobian
  doubling, mixed addition by `P0=(1,2)`, 254-bit scalar multiplication, affine
  inversion), `selftest`, `point <hex>`, `mktable <seed> <file>` and
  `stepwc <k> <file>`.
- `static_count.py`: static audit of the compiled step (objdump -> call graph
  -> assert straight-line -> sum counts x fixed call multiplicities).
- `STATIC_BOUND.txt`: recorded output of the audit (29,052 instructions/step).
- `crosscheck.c`: correctness test of the kernel against a schoolbook
  `Fp6`/`Fp12` reference built from the same `fp2_mul` primitive.
- `tcheck.c`: 200 random `Fp2` products against the Karatsuba reference built
  from `fp_mul` (the canonical Montgomery path), `bad=0/200`.
- `experiment.py`: standard-library recomputation of the ledger from `r` and the
  static price, the capped sampler and its failure bound, the mixer
  cross-check, and the small-subgroup walk evidence.
- `LOOP_BOUNDS.txt`: every loop of the algorithm with its hard bound.
- `sha256.c`, `sha256.h`: SHA-256 for the capped sampler (FIPS 180-4).
- `evidence/callgrind-identity.txt`: Ir per step on 10 production tables.
- `evidence/production-walk.txt`: frozen 64-trial small-subgroup output.
- `evidence/mixer-crosscheck.txt`: C `hashdump` output and the Python match.
- `evidence/clean-checkout.txt`: experiment run from a copy in /tmp without
  `.git`.

## Build, cross-check and self-test

```
gcc -O3 -fno-tree-vectorize -fno-stack-protector -mgeneral-regs-only \
    -o rho_pi fp12.c sha256.c rho_pi.c
./rho_pi selftest
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o crosscheck crosscheck.c fp12.c
./crosscheck
gcc -O2 -o tcheck tcheck.c fp12.c && ./tcheck
python3 static_count.py ./rho_pi
```

Output on the build machine:

```
g_is_one=0 g^r_is_one=1
assoc=1
fp6/fp12 vs schoolbook: 2000 pairs ok
add/sub round-trip: 2000 ok
bad=0/200
... TOTAL 29052
OK: every function in the call graph is straight-line (no conditional branch,
no loop, no early return).
```

`g^r = 1` checks the profile-basis conversion and the tower multiplication
together. `crosscheck` compares every product of 2,000 deterministic random
`Fp6` and `Fp12` pairs with the 9-multiplication schoolbook form and checks
2,000 `fp_add`/`fp_sub` round trips.

## Exact instruction count: static, from the disassembly

Machine: AArch64, Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0. Callgrind
3.22.0 (`--tool=callgrind`, event `Ir`) is used only to corroborate.

`static_count.py` parses `objdump -d --no-show-raw-insn ./rho_pi`, walks the
transitive call graph of the root symbol `rho_step` (compiled noinline, noipa),
and **fails** if any reached function contains a conditional branch or any other
control transfer (b / br / blr), if it lacks exactly one ret as its last
instruction, if a nop appears in the counted range, or if a call leaves the
binary. In the shipped binary the
assertion holds, so the executed instruction count of a call is exactly the
function's instruction count, and the step's count is a fixed sum:

| function | instructions | calls/step | Ir/step |
|---|---:|---:|---:|
| `rho_step` (generic, arbitrary salt) | 317 | 1 | 317 |
| `fp12_mul` | 2,581 | 1 | 2,581 |
| `fp6_mul` | 2,982 | 3 | 8,946 |
| `fp2_mul` (fused lazy reduction) | 956 | 18 | 17,208 |
| **total** | | | **29,052** |

The multiplicities come from the source: `fp12_mul` -> 3 `fp6_mul`;
`fp6_mul` -> 6 `fp2_mul`; `fp2_mul` is a fused leaf (three 4x4 products and one
Montgomery reduction per output coefficient) and calls nothing. No conditional
control transfer exists to make any count data-dependent. Multiplicities are DAG
path sums, not per-edge maxima, and the audit rejects an early return, a call
cycle, any unconditional/indirect branch and any nop in the counted range;
`static_count.py --selftest` runs a crafted counterexample for each rule.

Conditional operations are still correct but are instructions, not branches:
the Montgomery conditional subtract (`fp12.c` `fp_mul`) and the conditional add
of `p`/`r` (`fp12.h`, `rho_pi.c` `scalar_add_r`) are `cset`/`and`/`orr`/`neg`
sequences that execute on every operand.

Corroboration (`evidence/callgrind-identity.txt`): Callgrind with
`--toggle-collect='*rho_step*'` gives 29,052 Ir per step on ten independent
production tables (seeds 20261112, 1..9) at 1,000 steps, and at 2,000 and 5,000
steps on one table. Every run is identical and equals the static count, which
is what data independence means operationally.

Reproduce:

```
./rho_pi mktable 20261112 /tmp/table.bin
valgrind --tool=callgrind --toggle-collect='*rho_step*' \
    --callgrind-out-file=/tmp/cg.out ./rho_pi stepwc 1000 /tmp/table.bin
```

The Floyd comparison and loop control are outside the step and are charged
separately at 4,096 instructions per iteration, an analytic over-charge (48
word compares, three call instructions, counter and branch).

## Capped sampler cost

One SHA-256 call per 256-bit word. Worst case 66 scalars x K = 64 attempts =
4,224 attempts at 2,048 instructions each = 8,650,752, charged as its own ledger
row at the cap. The production table build (64 exponentiations of 254-bit
exponents, one `g^x`, one `Y0`) is charged inside the 2^36 setup row.

## Memory

32 records of 448 bytes (384-byte element + 32-byte `c_j` + 32-byte `d_j`),
two walk states, the hash scratch and code: below 2^20 bytes (the sampler
stores no words; they are consumed as produced). No table is stored beyond the
run, no advice is embedded, and the salt and increments are drawn fresh and
charged in setup.

The G1 part runs after the GT log is recovered and is charged at a fixed 2^24:
256 Jacobian doublings, at most one mixed addition each and one Fermat
inversion, all of fixed per-call cost, are far below that.

# Implementation and instruction measurement

## Files

- `fp12.h`, `fp12.c`: scalar 64-bit Montgomery `Fp`, tower `Fp2/Fp6/Fp12`
  (Karatsuba; `xi = 9+u`; `w_t^2 = v`).
- `rho.c`: profile-generator conversion, 32-adding walk, mixer, scalar
  arithmetic mod r, `selftest` and a `measure K` mode.
- `experiment.py`: standard-library recomputation of the ledger and a
  scope-limited small-group walk check.
- `evidence/callgrind-summary.txt`: annotated Callgrind totals for K = 10,000 (raw file omitted for the byte limit).

## Build and self-test

```
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho fp12.c rho.c
./rho selftest
```

Output on the build machine:

```
g_is_one=0 g^r_is_one=1
assoc=1
```

`g^r = 1` uses the profile generator converted from the profile basis to the
tower by the fixed linear map in `proof.md` section 1. This checks the
conversion and the tower multiplication together.

## Instruction counts

Machine: AArch64, 8-core Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0,
Valgrind 3.22.0 (`--tool=callgrind`, event `Ir`). The step loop was run at
K = 5,000, 10,000 and 20,000 with fixed inputs:

| K | PROGRAM TOTALS Ir |
|---:|---:|
| 5,000 | 237,821,066 |
| 10,000 | 474,013,018 |
| 20,000 | 946,407,371 |

Differences: (474,013,018 - 237,821,066)/5,000 = 47,238.39;
(946,407,371 - 474,013,018)/10,000 = 47,239.44. The fixed startup (table build,
conversion) cancels. **Charged: 52,000 instructions per step** (1.1007x the
larger measurement).

The measured step is exactly the loop body of `rho_step`: one `fp12_mul`, one
48-word mixer, two 4-word additions mod r and the table load. Floyd comparison
and loop control are not in this loop and are charged separately at 4,096
instructions per iteration.

Reproduce:

```
valgrind --tool=callgrind --callgrind-out-file=cg.out ./rho 10000
callgrind_annotate cg.out | grep 'PROGRAM TOTALS'
```

## Memory

32 records of 448 bytes (384-byte element + 32-byte `c_j` + 32-byte `d_j`),
two walk states (about 800 bytes), 48-word hash scratch and code: below 2^20
bytes. No table is stored beyond the run, no advice is embedded, and the salt
and increments are drawn fresh and charged in setup.

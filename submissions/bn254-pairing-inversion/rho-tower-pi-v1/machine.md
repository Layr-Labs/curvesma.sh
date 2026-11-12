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
- `experiment.py`: standard-library recomputation of the ledger, the G1
  vectors, and a scope-limited small-group walk check.
- `evidence/callgrind-summary.txt`: annotated Callgrind totals.

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

## Instruction counts: GT walk

Machine: AArch64, Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0, Valgrind
3.22.0 (`--tool=callgrind`, event `Ir`). The step loop was run at
K = 5,000, 10,000 and 20,000 with fixed inputs:

| K | PROGRAM TOTALS Ir |
|---:|---:|
| 5,000 | 237,818,071 |
| 10,000 | 474,005,018 |
| 20,000 | 946,389,371 |

Differences: (474,005,018 - 237,818,071)/5,000 = 47,237.39;
(946,389,371 - 474,005,018)/10,000 = 47,238.44. The fixed startup (table build,
conversion) cancels. **Charged: 52,000 instructions per step** (1.1007x the
larger measurement).

The measured step is exactly the loop body of `rho_step`: one `fp12_mul`, one
48-word mixer, two 4-word additions mod r and the table load. Floyd comparison
and loop control are charged separately at 4,096 instructions per iteration.

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
below 2^20 bytes. No table is stored beyond the run, no advice is embedded, and
the salt and increments are drawn fresh and charged in setup.

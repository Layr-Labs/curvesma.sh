# Measured step kernels (Callgrind Ir, organizer toolchain)

Built and counted exactly like `scripts/reference_costs/run.sh`, using the same image
(python:3.12-slim pinned digest, GCC 14.2.0, `-O3 -fno-tree-vectorize -mgeneral-regs-only`,
Valgrind Callgrind, linux/arm64) and the same pinned checkouts: libsecp256k1 e74a0187, blst f7a88067.
The same rig first reproduced the organizer totals exactly (field 733,204/1,466,388;
Fp12 94,191,636/188,383,252; blst Fp 1,568,788/3,137,556). Each kernel is counted at 4,096 and
8,192 iterations, and the slope (difference/4,096) is the per-iteration cost.
Field inputs are random. All modular arithmetic used is constant-time apart from libsecp256k1's
normalize_var/cmp_var/inv_var, whose random-input cost is representative.
Reproduce: mount the kernels at /k and run `sh /k/run_kernels.sh` in that image with /secp and /blst
mounted read-only, as in reproduce.sh.

| kernel | slope (Ir per iteration) | meaning |
|---|---:|---|
| secp_m0 | 1,061.1 | batched incomplete affine addition pipeline only (B=256) |
| secp_m1 | 1,512.7 | full class-rho step |
| secp_m2 | 1,027.5 | m0 without the per-batch inversion → inversion ≈ 8.7k per batch |
| bls_rho_m0 | 2,698.0 | batched G1 addition pipeline only |
| bls_rho_m1 | 3,281.5 | full class-rho step on G1 |
| bls_rho_m2 | 2,591.4 | m0 without the per-batch inversion → inversion ≈ 27.4k per batch |
| cheon_m0 | 5,363.0 | two batched G1 additions from window tables |
| cheon_m1 | 6,096.4 | full Cheon baby element (fr mul, windows, 3 loads, 2 adds, class key, hash insert) |
| cheon_m2 | 5,149.8 | m0 without inversions |
| cheon_m3 | 6,092.4 | full Cheon giant element (lookup in a table at load 0.5) |

Measured non-reference overhead per step = full − add-only:
- secp: 451.6
- G1 ρ: 583.5
- Cheon: 733.4 (baby) / 729.4 (giant)

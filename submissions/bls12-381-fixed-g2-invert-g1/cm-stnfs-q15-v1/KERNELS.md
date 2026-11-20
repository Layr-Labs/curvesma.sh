# Kernels (`Zero SIMD`, `Zero Hand-Written Assembly`, `Zero FEAT_LSE`)

Every hot kernel in this package is **100% pure portable C**, compiled to `aarch64-linux-gnu` baseline scalar AArch64 assembly with `clang -S --target=aarch64-linux-gnu -march=armv8-a -DKERNELS_ONLY -O3 -fno-builtin -fomit-frame-pointer -fno-stack-protector -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`. Zero vector/NEON SIMD instructions, zero inline assembly, zero `FEAT_LSE` atomics (`staddh` is not used), zero `bl` libc calls, and zero 9D record-walk code are used. Every cross-owner memory transfer is explicitly priced below.

## 1. Collection Kernels (Tiled `Arch2-Rev-SigmaB`, `P_lo = 2^54`, `W_tile = 256`, `W_SB = 128`)

| Kernel | Role & Exact Baseline ARMv8-A Scalar Instruction Breakdown | Measured Ir | Charged Constant |
|---|---|---:|---:|
| `rev_sb_probe_fast` | Band-2 `[2^54, P_hi)` tile start: mod-`p` step (`add/cmp/csel/sub`), super-bucket range `>> 7`, bitmap probe (`ret` incl.) | **18** / tile | `C_TILE_HDR = 20.0` |
| Occupied SB overhead | Conservatively assumes `P(occ) = 1.0` on every spanned SB: bitmap/head load (`4` Ir) + `rev_sb_scan_emit_t` prologue & terminating `key_a > -q_min` exit (`14` Ir) | **18** / SB | `C_SB_HDR = 18.0` |
| `rev_sb_scan_emit_t` (cand) | Band-2 key-sorted candidate check (`ldur/cmp/b.gt`, `ldr/sub/cmp/b.hs`, `ldr/cmp/b.gt`, loop dec/add/b.eq = 13 Ir/cand) | **13** / cand | `C_CAND = 13.0` |
| Band-2 hit end-to-end | `rev_sb_scan_emit_t` emit w/ `bfxil` conj rot (**9** Ir) + 1 MiB staging-page copy/handoff to col owner `t` (**4** Ir, above `c36`'s 3-Ir 16B copy) + `rev_b2_apply_u4` 4x-unrolled scalar apply (**23/4 = 5.75** Ir) | **18.75** / hit | `C_HIT_HI = 19.0` |
| `rev_b1_bsearch` | Band-1 `[2^24, 2^54)` branchless `csinc` upper-bound binary search (`lsr/add/mvn/ldr/cmp/csinc/adds/b.ne`) | **8** / step | `C_BSEARCH_STEP = 8.0` (`hdr = 24.0`) |
| Band-1 hit end-to-end | `rev_b1_apply_u4` 4x-unrolled scalar apply w/ 2-instr `sub + and #0xfffffffffffffff` conjugate rotation (**25/4 = 6.25** Ir) + contiguous 12B slice copy to col owner `t` (**1.50** Ir, `c36` rate `18/16`) | **7.75** / hit | `C_HIT_LO = 8.0` |
| `rev_build_pass1/2/3` | 3-pass `a_0`-chord + `perm[]` build: Pass 1 chord fill (**7**) + Pass 2 count/bitmap (**13**) + Pass 3 sorted scatter (**14**) + chord start (**1**) | **35** / entry | `C_LIST_A = 35.0` (`C_LBKT = 2.25`) |
| `rev_scan_zero_u4` | `-fno-builtin` 64-bit SWAR 4-cell constant-threshold `T_g(b)` scan + zero `S_g, S_f` (**11/4 = 2.75** Ir/cell) x `V_all/V_kept = 1.0276` + `0.10` setup | **2.93** / kept cell | `C_SCAN_ZERO = 2.95` |

Exactness of the 4-orbit `sigma`-shared `a`-table replay with `D_(omega,sigma)` (`239/239` exact-1, `0` duplicates, `0` misses) and unified 2-adic key `key_a = q_a - s_c(a)` is verified live inside `audit_cm.py` (`sigmaB_exactness`: `48` real BLS12-381 `w' = omega^2` ideals across `12` shared tables with `0` mismatches against brute force).

## 2. Exact Cyclotomic Tower Norm Kernels

Exploiting the tower `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q` with 3-mul Eisenstein Karatsuba and 2-mul quadratic norms `(h_x - h_y)^2 + h_x h_y` (fixed limbs, zero division, zero SIMD; verified `0` mismatches vs PARI/GP on production pairs):
- `_norm_z15_2limb` (exact degree-8 norm to 12 limbs): **3,357 Ir**.
- `_nk_g_only` (exact 678-bit `N_g` + `lg16` + threshold recheck): **4,204 Ir** (charged `C_NORMS_G = 4,500`; `10.7x` faster than `c36`'s `48,300`).
- `_nk_f_only` (exact 549-bit `N_f`, reusing `a^3, b^3`): **6,167 Ir** (charged `C_NORMS_F = 6,500`).

## 3. Stage B, FFT Scaled Remainder Tree, & Live Cofactorization

- **Two-sided Stage B (`P24` = product of primes `< 2^24`, `1,077,871` primes, `24,197,143` bits):** Scaled conservatively upward from `c36`'s measured `820,000` Ir at 740-bit `N_f` via `ceil_10k(820,000 x bits / 740)`:
  - `N_g <= 678.0` bits: `820,000 x 678 / 740 = 751,297` -> conservatively charged **`C_SB_G = 760,000`** Ir.
  - `N_f <= 549.0` bits: `820,000 x 549 / 740 = 608,351` -> conservatively charged **`C_SB_F = 620,000`** Ir.
- **FFT SRT (`H-FFT`):** Identical bound `F(N) = 8.5 N (log2 N)^2 + 1200` (`+ 504 N` past `2^26` limbs) as `c36`. `audit_cm.py` (`srt_toy`) verifies exact SRT remaindering on 128 planted leaves with `0` mismatches.
- **Live Cofactorization (`B_f = 2^71.29, B_g = 2^71.69`):** Pure-stdlib `cm_models.Retention`, `gbatch_cost`, `ecm_instr`, and 40-round Miller-Rabin (`MONT5 = 367.76`) are executed **LIVE** inside `audit_cm.py` (`r_g = 0.810321`, `r_f = 0.938170`, `r_tot = 0.760219`, `c_cof = 1.7071` Ir/kept cell at `(s_a_g, s_b_g) = (360, 295)`, `(s_a_f, s_b_f) = (420, 345)`).

## 4. `Z[zeta_12]` Block Wiedemann Kernels

Pure-C straight-line 4-limb arithmetic (`__builtin_addcll`/`__builtin_subcll`, zero SIMD, zero `bl`, verified bit-for-bit on `100,000` random/edge-case `chi`-multiplies and `10,000` 6-bucket row-tail reductions, `errs = 0`):

| Kernel | Role | Exact Scalar AArch64 Ir | Comparison vs `c36` |
|---|---|---:|---|
| `chi_mul_z12_one` | Canonical `chi y = -|u| y mod r in [0, r)` (4x1 `mul/umulh` + 1-limb Barrett; `|u| < 2^64`) | **58** | `2.29x` faster than `c36` (`133`) |
| `chi_mul_z12_fast` | Non-canonical `chi y in [0, 2r)` | **50** | `2.66x` faster than `c36` (`133`) |
| `row_red5_z12` | Exact 5-limb Barrett `pos - neg mod r in [0, r)` | **53** | `-12` Ir vs `c36` (`65`) |
| `row_tail_6bucket_z12` | Arch B 6-bucket fold `6 -> 4`, 4x `row_red5_z12`, 3x Horner `chi`-mul + add | **504** | Complete 1-owner row tail |
| `unit_run5_z12` | 5-vector unit accumulation (all 12 roots of unity) | **39 / 5 = 7.80** / vec-entry | Identical to `c36` (`7.80`) |
| `gen_run1_z12` | 1-vector `uint8_t` general accumulation | **25.0** / vec-entry (charged **33.0**) | Conservative `c36` rate |
| `red_row_u_z12` | Dense-band subset-sum table element (5 vectors) | **18.6** / elem | Matches `c36` (`18.4`) |

Because 1-owner Architecture B fits the input 5-vector block `x`, output 5-vector block `y`, AND all 22 dense-band tables inside **`2^63.751` B `< 2^64` B** on ONE single owner, cross-owner sends (`_send_pair`), combines (`_combine`), and routing are **zero**. Fixed `owner_part_6bucket_z12` overhead is `175` Ir/vec typical / `265` Ir/vec worst-case (all 12 general buckets active), plus `row_tail_6bucket_z12` (`504`) + header/args (`16`) = `695` typical / **`785` worst-case** (conservatively charged `785` in full). Computing `k5` directly from `exp_unit_lo = 0.93008 - 2 x 0.00010` via `eu x 7.80 + (1 - eu) x 33.0` gives `k5 = 9.5620` (at `-2 SE`: `9.5669`), conservatively charged **`9.59`** (verified live in `audit_cm.py`).

## 5. Descent Kernels

Descent nodes use `c36`'s exact `dsieve16` accounting (`1,100` Ir/hit, `5` scan, `1` zero, `2,400` bound, `2^25.1` setup) conservatively multiplied by **`1.5x`** for 16D rank-2 lattices (`DESCENT.md`).

# Kernels (Closed-Form & Compiled Baseline Scalar AArch64: Zero SIMD, Zero ASM, Zero `FEAT_LSE`)

Every hot kernel is pure portable C compiled to baseline scalar AArch64 (`clang -S --target=aarch64-linux-gnu -march=armv8-a -O3 -fno-builtin -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`). Zero SIMD, zero inline assembly, zero `FEAT_LSE`, zero `bl` libc calls, and zero 9D record-walk code are used.

## 1. Collection Kernels (`Arch2-Rev-SigmaB-v3`, `W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

| Kernel | Role & Exact Baseline ARMv8-A Scalar Instruction Breakdown | Measured Ir | Charged Constant |
|---|---|---:|---:|
| `rev_sb_probe_fast` | Band-2 `[2^55, P_hi)` full-chord row start: mod-`p` step, SB range `>> 6`, header setup (`ret` incl.) | **18** / row | `C_TILE_HDR = 20.0` |
| `rev_b2_minkey_gate` | Per-SB 1-byte `sb_min_q8[j] <= lim_q8[j]` gate (`+1/254` quant slack; skips **`85.4%`** of spanned SBs) | **6** skip / **20** enter | `C_SB_SKIP = 6.0`, `C_SB_ENTER = 20.0` |
| `rev_b2_scan_emit_8b` (cand) | Per-SB tight `q_min_j` candidate check (`ldur/cmp/b.gt`, `ldr/sub/cmp/b.hs`, `ldr/cmp/b.gt`, loop = 13 Ir/cand) | **13** / cand | `C_CAND = 13.0` |
| Band-2 hit end-to-end | `rev_b2_scan_emit_8b` 8B `oid_m` emit grouped by `0.25`-bit `log_p` (**8**) + `rev_copy_u64_u4` `ldp/stp` 8B copy (**1.75**, chg **2.0**) + `rev_apply_8b_u4` (**4.25**) | **14.25** / hit | `C_HIT_HI = 14.5` |
| Band-1a hit (`[2^24, 2^48)`) | `m`-partitioned 4-bank `a`-tables (`acc_m[(m-m_b)&3]` once/bucket): `rev_apply_8b_u4` (**4.25**) + `rev_copy_u64_u4` (**1.75**) | **6.00** / hit | `C_HIT_LO_M4 = 6.25` |
| `rev_b1_early4_find` (`[2^48, 2^55)`) | Early-exit `keys[0] > neg_qb` (**5** + **6** ptr/res = **11** Ir when `cnt=0`; `+14 + cnt` 4-strided scan when `cnt>0`; fallback `24 + 8 log2`) | **11** empty | exact 16D integral |
| Band-1b hit (`[2^48, 2^55)`) | Unified `rev_b1_apply_u4` w/ 2-instr conj rot (**6.25**) + 12B slice copy (`3` `ldp/stp` per 4 = **1.50**) | **7.75** / hit | `C_HIT_LO_U = 8.0` |
| `rev_build_pass1/2/3` | 3-pass `a_0`-chord + `perm[]` build (Pass 1 **7** + Pass 2 **13** + Pass 3 **14** + start **1**; `sb_min_q8` + prefix sum `4.0`/SB) | **35** / entry | `C_LIST_A = 35.0` (`C_LBKT = 4.0`) |
| `rev_scan_zero_u4` | `-fno-builtin` 64-bit SWAR 4-cell constant-threshold `T_g(b)` scan + zero `S_g, S_f` (**2.75** Ir/cell) x `1.0276` + `0.10` | **2.93** / kept cell | `C_SCAN_ZERO = 2.95` |

Exactness of 4-orbit `sigma`-shared replay with `D_(omega,sigma)` (`239/239` exact-1), per-SB `q_min_j`, and 1-byte `sb_min_q8` gate is verified live in `audit_cm.py` (`sigmaB_exactness`: 48 real BLS12-381 `w'=omega^2` ideals, `0` mismatches).

## 2. Exact Cyclotomic Tower Norm Kernels

Tower `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q` with 3-mul Eisenstein Karatsuba and 2-mul quadratic norms (fixed limbs, zero division, zero SIMD; verified `0` mismatches vs PARI/GP): `_norm_z15_2limb` **3,357** Ir, `_nk_g_only` **4,204** Ir (charged **4,500**), `_nk_f_only` **6,167** Ir (charged **6,500**).

## 3. Closed-Form Non-Recursive 4-Prime 62-Bit Radix-2 NTT (`fft_F`), Analytical Stage B, & 2D Cofactorization

- **Closed-form non-recursive 4-prime 62-bit radix-2 NTT (`cm_models.fft_F`, ZERO recursion & ZERO FFT heuristic):**
  Fix the 4 explicit 62-bit NTT primes `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`2422936599525326849`, `2485986994308513793`, `2747195772696002561`, `2783224569714966529`; all 4 verified prime LIVE in `audit_cm.py` with `log2(NTT4_PROD) = 244.70`). Every `N`-limb (`64`-bit) product uses ONE of two exact non-recursive branches (both preconditions asserted live inside `fft_F(N)`):
  1. **1D branch (`2 N <= 2^53`, covering all `P24`, Stage B, and BM generator products up to `18 L ~ 2^51.8` limbs):** Pad to `M = 2^ceil(log2(2 N)) <= 2^53` (asserts `M <= 2^53` and `N * 2^128 < NTT4_PROD`): `4` primes x `3` radix-2 NTTs x `(M/2) log2(M)` Barrett butterflies at `28` Ir (`26` Ir scalar count) = **`168 M log2(M)`** + `4 M` pointwise muls (`80 M`) + `M` 4-prime Garner CRT/carry (`140 M`) + `1200` setup + four-step 2-pass flat-memory transpose past `2^26` limbs (`+504 M`):
     `fft_F(N) = 168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`.
  2. **2D polynomial convolution branch (`2 N > 2^53`, covering `gbatch_cost` `P_g, P_f` trees at `N ~ 2^66` up to `N <= 2^100`):** Split the `N`-limb integer into `K_2 = ceil(N / K_1)` rows of width `K_1 = 2^(ceil(log2(N)) // 2)` (`Y = X^(K_1)`), zero-pad BOTH axes by `2x` to `M_1 = 2 K_1 <= 2^53` and `M_2 = 2^ceil(log2(2 K_2)) <= 2^53` (exact `bit_length()`, ZERO circular wrap on either axis, asserts `M_1, M_2 <= 2^53` and `K_1 K_2 * 2^128 < NTT4_PROD` — e.g. `2^195 << 2^244.7` at `N = 2^67`). With `M_tot = M_1 M_2`, run `M_2` length-`M_1` + transpose + `M_1` length-`M_2` NTTs per transform (`168 M_tot log2(M_tot)`), pointwise mul + 4-prime Garner CRT + 2D overlap-add carry (`260 M_tot + 1200`), 4 full 2D transposes (`4 x 84 M_tot = 336 M_tot`), and out-of-core flat-memory streaming (`504 M_tot`):
     `fft_F(N) = 168 M_tot log2(M_tot) + 1100 M_tot + 1200`.
  `audit_cm.py` (`srt_toy`) verifies exact scaled-remainder-tree remaindering on 128 planted leaves (`0` mismatches).
- **Analytical Stage B (`cm_models.stage_b_analytic`):** For an `nb`-bit survivor (`L = ceil(nb / 64)` limbs), `stage_b_analytic(nb) = gbatch_cost(24197143, L) + 10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16)` (`MONT5 = 367.76` Ir for 5-limb CIOS Montgomery mul). At `550`b (`L=9`): `5,721,911` Ir; at `700`b (`L=11`): `8,479,620` Ir.
- **Live 2D-Histogram Cofactorization:** Integrated LIVE across all `554` bins up to `bg_max = 700`b via `cm_models.Retention` + `stage_b_analytic` at `(s_a_g, s_b_g) = (345, 295)`, `(s_a_f, s_b_f) = (405, 345)`, giving `r_eff = 0.68001` and strict bin-wise **`c_cof_2d_upper = 3.4402`** ipc.

## 4. `Z[zeta_12]` Block Wiedemann Kernels

Pure-C straight-line 4-limb arithmetic (`0` errors on `100,000` `chi`-muls & `10,000` 6-bucket tails): `chi_mul_z12_one` **58** Ir, `chi_mul_z12_fast` **50** Ir, `row_red5_z12` **53** Ir, `row_tail_6bucket_z12` **504** Ir, `unit_run5_z12` **7.80** Ir/vec-entry, `gen_run1_z12` **25.0** (charged **33.0**), `red_row_u_z12` **18.6** Ir/elem. Fixed 1-owner Arch B row overhead = `265 + 504 + 16 =` **`785`** Ir worst-case (`0` cross-owner sends/combines; `k5_cap = 10.50` vs `2^23` measured `9.5620`).

## 5. Descent Kernels (`16D` Rank-2 `dsieve16` & Limb-Scaled Miller-Rabin)

Derived directly from baseline scalar ARMv8-A counts (`1.5x` safety factor): `SIEVE_HIT = 1,650`, `SCAN + ZERO = 18`, `BOUND_CELL = 3,600`, `PAIR_SETUP = 1.5 x 1.33 x 2^25.1`, `LLL_NODE = 2^38` (16D LLL on 381-bit entries), `NORMS = 2^27`, `TDIV = 2^25`, `IDEAL = 2^35`, `ecm_instr` (`MONT5 = 367.76`), and exact limb-scaled 100-round Miller-Rabin **`mr_instr(nb) = 100 * 2 * nb * MONT5 * (ceil(nb/64)/5)^2`** (`1.29e+10` Ir at `2600`b; `<= 4^(-100) = 2^(-200)` error/test, `<= 2^(-80)` over `< 2^120` run calls; zero fixed `PRIME` constant).

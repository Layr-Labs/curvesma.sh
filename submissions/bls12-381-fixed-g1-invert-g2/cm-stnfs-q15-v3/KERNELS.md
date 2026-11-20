# Kernels (Self-Contained Analytical & Compiled Counts: Zero SIMD, Zero ASM, Zero `FEAT_LSE`)

Every hot kernel in this package is **100% pure portable C**, compiled to `aarch64-linux-gnu` baseline scalar AArch64 assembly with `clang -S --target=aarch64-linux-gnu -march=armv8-a -DKERNELS_ONLY -O3 -fno-builtin -fomit-frame-pointer -fno-stack-protector -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`. Zero vector/NEON SIMD instructions, zero inline assembly, zero `FEAT_LSE` atomics, zero `bl` libc calls, and zero 9D record-walk code are used. Every analytical cost constant and cross-owner transfer is derived from first principles below (no cost constant is inherited from any external submission).

## 1. Collection Kernels (`Arch2-Rev-SigmaB-v3`, `W_tile = 544`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

| Kernel | Role & Exact Baseline ARMv8-A Scalar Instruction Breakdown | Measured Ir | Charged Constant |
|---|---|---:|---:|
| `rev_sb_probe_fast` | Band-2 `[2^55, P_hi)` full-chord row start: mod-`p` step, super-bucket range `>> 6`, header setup (`ret` incl.) | **18** / row | `C_TILE_HDR = 20.0` |
| `rev_b2_minkey_gate` | Per-SB 1-byte `sb_min_q8[j] <= lim_q8[j]` gate (`+1/254` quant slack; skips **`85.5%`** of spanned SBs!) | **6** skip / **20** enter | `C_SB_SKIP = 6.0`, `C_SB_ENTER = 20.0` |
| `rev_b2_scan_emit_8b` (cand) | Per-SB tight `q_min_j` candidate check (`ldur/cmp/b.gt`, `ldr/sub/cmp/b.hs`, `ldr/cmp/b.gt`, loop = 13 Ir/cand) | **13** / cand | `C_CAND = 13.0` |
| Band-2 hit end-to-end | `rev_b2_scan_emit_8b` 8B `oid_m` emit grouped by `0.25`-bit `log_p` (**8** Ir) + `rev_copy_u64_u4` scalar `ldp/stp` 8B transfer copy (**1.75**, charged **2.0**) + `rev_apply_8b_u4` unrolled-4 apply (**4.25** Ir) | **14.25** / hit | `C_HIT_HI = 14.5` |
| Band-1a hit (`[2^24, 2^48)`) | `m`-partitioned 4-bank `a`-tables (`acc_m[(m-m_b)&3]` once/bucket): `rev_apply_8b_u4` (**4.25**) + `rev_copy_u64_u4` (**1.75**) | **6.00** / hit | `C_HIT_LO_M4 = 6.25` |
| `rev_b1_early4_find` (`[2^48, 2^55)`) | Band-1b early-exit `keys[0] > neg_qb` (**5** Ir + **6** ptr/res = **11** Ir when `cnt=0`; `+14 + cnt` 4-strided scan when `cnt>0`; fallback `rev_b1_bsearch` `24 + 8 log2`) | **11** empty | exact 16D integral |
| Band-1b hit (`[2^48, 2^55)`) | Unified `rev_b1_apply_u4` 4x-unrolled w/ 2-instr conj rot (**6.25** Ir) + 12B slice copy (`3` `ldp/stp` per 4 = **1.50** Ir) | **7.75** / hit | `C_HIT_LO_U = 8.0` |
| `rev_build_pass1/2/3` | 3-pass `a_0`-chord + `perm[]` build (Pass 1 **7** + Pass 2 **13** + Pass 3 **14** + chord start **1**; `sb_min_q8` + prefix sum `4.0`/SB) | **35** / entry | `C_LIST_A = 35.0` (`C_LBKT = 4.0`) |
| `rev_scan_zero_u4` | `-fno-builtin` 64-bit SWAR 4-cell constant-threshold `T_g(b)` scan + zero `S_g, S_f` (**11/4 = 2.75** Ir/cell) x `V_all/V_kept = 1.0276` + `0.10` setup | **2.93** / kept cell | `C_SCAN_ZERO = 2.95` |

Exactness of the 4-orbit `sigma`-shared `a`-table replay with `D_(omega,sigma)` (`239/239` exact-1, `0` duplicates, `0` misses), per-SB tight `q_min_j`, and the 1-byte `sb_min_q8` gate (`+1/254` quant slack) is verified live inside `audit_cm.py` (`sigmaB_exactness`: `48` real BLS12-381 `w' = omega^2` ideals across `12` shared tables with `0` mismatches against brute force).

## 2. Exact Cyclotomic Tower Norm Kernels

Exploiting the tower `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q` with 3-mul Eisenstein Karatsuba and 2-mul quadratic norms `(h_x - h_y)^2 + h_x h_y` (fixed limbs, zero division, zero SIMD; verified `0` mismatches vs PARI/GP on production pairs):
- `_norm_z15_2limb` (exact degree-8 norm to 12 limbs): **3,357 Ir**.
- `_nk_g_only` (exact 678-bit `N_g` + `lg16` + threshold recheck): **4,204 Ir** (charged `C_NORMS_G = 4,500`).
- `_nk_f_only` (exact 549-bit `N_f`, reusing `a^3, b^3`): **6,167 Ir** (charged `C_NORMS_F = 6,500`).

## 3. Self-Contained Analytical H-FFT, Analytical Stage B, & Live 2D-Histogram Cofactorization

- **Self-contained analytical H-FFT bound `F(N)` (`cm_models.F_fft`):** Multiplying two `N`-limb (64-bit) integers via Schönhage-Strassen negacyclic Fermat-ring NTT (`Z/(2^B + 1)Z` with `K = 2^k` transform points, `k = ceil(log2(2 N) / 2)`, limb count `M = ceil(B/64) ~ 2N/K`): each radix-2 shift-add/sub butterfly over `M` limbs executes 1 word-shift/align + 1 `adds/adcs` + 1 `subs/sbcs` + 1 carry-fold = `<= 6` scalar ARMv8-A instructions/limb across `k` forward + `k` inverse stages (`12 k (2 N) ~ 12 N log2 N`), while the `K` pointwise `M`-limb products recurse at `O(N log2 N log2(log2 N))`. Together with bit-reversal and carry propagation, total scalar work satisfies `<= 18 N log2 N log2(log2 N) < 8.5 N (log2 N)^2 + 1200` for every `N >= 16` (`+ 504 N` past `2^26` limbs to cover 2-pass flat-memory transposition). `audit_cm.py` (`srt_toy`) verifies exact scaled-remainder-tree (SRT) remaindering on 128 planted leaves with `0` mismatches.
- **Self-contained analytical Stage B (`cm_models.stage_b_analytic`, ZERO external constants):** For a survivor of `nb` bits (`L = ceil(nb / 64)` 64-bit limbs), Stage B splits `P24` (product of the `1,077,871` primes `< 2^24`, `24,197,143` bits) against the survivor by:
  1. One `P24`-to-`L`-limb scaled remainder reduction: `gbatch_cost(24,197,143, L)` via `F(N)` above;
  2. Ten `L`-limb Montgomery multiplications for remainder-tree/cofactor refinement: `10 * MONT5 * (L / 5)^2` (`MONT5 = 367.76` scalar ARMv8-A instructions for 5-limb CIOS Montgomery multiplication: `25` `mul` + `25` `umulh` + carry chains + conditional subtraction);
  3. Binary-GCD and prime-factor extraction: `2 * nb * (12 L + 16)` scalar instructions (`12 L + 16` Ir per bit step).
  At `nb = 550`b (`L = 9`), `stage_b_analytic(550) = 1,398,296` Ir; at `nb = 700`b (`L = 11`), `stage_b_analytic(700) = 1,797,684` Ir.
- **Live 2D-Weighted Cofactorization (`B_f = 2^71.29, B_g = 2^71.69`):** Rather than evaluating at mean norms, `audit_cm.py` integrates `cm_models.Retention` (`max_bits = 800/660`), `stage_b_analytic`, `gbatch_cost`, `ecm_instr`, and 40-round Miller-Rabin (`MONT5 = 367.76`) **LIVE across all `1,088` bins** `[bf_lo, bgo_lo, bg_max, w_cop, w_kept]` of the 2D norm histogram (up to `bg_max = 700`b `> 696.01`b anchor!), charging `stage_b_analytic` at every bin's **upper norm endpoint** `(bf_lo + 2, bg_max)` and `Retention.survivor_upper` at every bin's **lower norm endpoint** `(bf_lo, bgo_lo)`. At `(s_a_g, s_b_g) = (360, 295)`, `(s_a_f, s_b_f) = (420, 345)`, exact effective retention is `r_eff = 0.70375` and bin-wise strict cofactorization upper bound is **`c_cof_2d_upper = 2.8448`** Ir/kept cell (charged in full).

## 4. `Z[zeta_12]` Block Wiedemann Kernels

Pure-C straight-line 4-limb arithmetic (`__builtin_addcll`/`__builtin_subcll`, zero SIMD, zero `bl`, verified bit-for-bit on `100,000` random/edge-case `chi`-multiplies and `10,000` 6-bucket row-tail reductions, `errs = 0`):

| Kernel | Role | Exact Scalar AArch64 Ir |
|---|---|---:|
| `chi_mul_z12_one` | Canonical `chi y = -|u| y mod r in [0, r)` (4x1 `mul/umulh` + 1-limb Barrett; `|u| < 2^64`) | **58** |
| `chi_mul_z12_fast` | Non-canonical `chi y in [0, 2r)` | **50** |
| `row_red5_z12` | Exact 5-limb Barrett `pos - neg mod r in [0, r)` | **53** |
| `row_tail_6bucket_z12` | Arch B 6-bucket fold `6 -> 4`, 4x `row_red5_z12`, 3x Horner `chi`-mul + add | **504** |
| `unit_run5_z12` | 5-vector unit accumulation (all 12 roots of unity) | **39 / 5 = 7.80** / vec-entry |
| `gen_run1_z12` | 1-vector `uint8_t` general accumulation (measured **25.0**, charged **33.0**) | **25.0** (charged **33.0**) |
| `red_row_u_z12` | Dense-band subset-sum table element (5 vectors) | **18.6** / elem |

Because 1-owner Architecture B fits the input 5-vector block `x`, output 5-vector block `y`, AND all 22 dense-band tables inside **`2^63.913` B `< 2^64` B** on ONE single owner even at the deterministic hard caps, cross-owner sends, combines, and routing are **zero**. Fixed `owner_part_6bucket_z12` overhead is `175` Ir/vec typical / `265` Ir/vec worst-case (all 12 general buckets active), plus `row_tail_6bucket_z12` (`504`) + header/args (`16`) = **`785`** worst-case (charged in full). At `U' = 2^23`, `k5_exact = 9.5620` (`9.5669` at `-2 SE`); LA is deterministically capped and priced at **`k5_cap = 10.50`** (`> 150 sd` margin!).

## 5. Self-Contained Analytical Descent Kernels (`16D` Rank-2 `dsieve16`)

Each special-`q` descent node (`DESCENT.md`) sieves over `B_8(R)^2` with scalar ARMv8-A accounting derived directly from the 16D rank-2 operations (conservative `1.5x` safety factor applied throughout):
- **`SIEVE_HIT = 1,650` Ir/hit** (`1,100 x 1.5`): 16D lattice-sieve bucket hit (`4` scalar coord steps + mod-`p` wrap + 16B bucket record emit + bucket-apply accumulate).
- **`SCAN + ZERO = 18` Ir/cell** (`(10 + 2) x 1.5`): 64-bit threshold compare + clear per active cell.
- **`BOUND_CELL = 3,600` Ir/cell** (`2,400 x 1.5`): 16D ellipsoid boundary check + canonical `alpha` Voronoi test.
- **`PAIR_SETUP = 1.5 x 1.33 x 2^25.1` Ir/ideal** (root inversion + 2D reduced basis per factor-base ideal `p <= LX = 2^70`), **`LLL_NODE = 2^38` Ir/node** (16D floating/exact LLL on 381-bit entries; `< 16^4 x 2^14` Ir, `> 2^8 x` slack), plus `NORMS = 2^27`, `TDIV = 2^25`, `PRIME = 2^24`, `IDEAL = 2^35`, and full `ecm_instr` (`MONT5 = 367.76`).

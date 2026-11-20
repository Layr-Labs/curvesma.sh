# Kernels (100% Analytical Baseline Scalar ARMv8-A: Zero SIMD, Zero `FEAT_LSE`, Zero External C Measurements)

Every kernel cost constant in this package is derived explicitly by instruction-by-instruction baseline scalar ARMv8-A (`-march=armv8-a -mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `bl` libc calls, zero external unbundled C compilation measurements) and verified symbolically inside `audit_cm.py` (`kernels_symbolic`).

## 1. Collection Kernels (`Arch2-Rev-SigmaB-v3`, `W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

| Kernel | Explicit Baseline Scalar ARMv8-A Instruction Breakdown (`-mgeneral-regs-only`) | Analytical Ir | Charged Constant |
|---|---|---:|---:|
| `rev_sb_probe_fast` | Band-2 `[2^55, P_hi)` full-chord row start: Barrett mod-`p` step (`sub, mul, umulh, lsr, msub` = 5) + SB index shift & clamp (`add, lsr, sub, cmp, csel` = 5) + pointer/tile header setup (`lsl, add, ldp, cmp` = 4) + branch/exit (`b.gt, add, add, ret` = 4) | **18** / row | `C_TILE_HDR = 20.0` |
| `rev_b2_minkey_gate` (skip) | Per-SB 1-byte `sb_min_q8[j] <= lim_q8[j]` gate (`+1/254` quant slack; skips `85.4%` of spanned SBs): `ldrb w0, [x1], #1` (1), `ldrb w2, [x3], #1` (1), `cmp w0, w2` (1), `b.ls .Lenter` (1), `subs x4, x4, #1` (1), `b.ne .Lloop` (1) | **6** / skip | `C_SB_SKIP = 6.0` |
| `rev_b2_minkey_gate` (enter) | Load SB offsets (`ldr, ldr` = 2), empty check (`subs, b.eq` = 2), base/end pointers (`add .. lsl #4` x2 = 2), inner-loop register setup (`mov, ldr, ldp` x3 = 8), post-SB exit & loop (`subs, b.ne, add, b` = 4) | **18** / enter | `C_SB_ENTER = 20.0` |
| `rev_b2_scan_emit_8b` (cand) | 16B `{q_key, a0_lo, oid_m}` candidate check: `ldp w0, w1, [x_p], #16` (1), `cmp w0, w_lim; b.hi .Ldone` (2), `sub w2, w1, w_base; cmp w2, w_len; b.hs .Lnext` (3), `ldr w3, [x_qtab, x2, lsl #2]; cmp w0, w3; b.hi .Lnext` (3), `cmp x_p, x_end; b.lo .Lcand` (2) | **11** / cand | `C_CAND = 13.0` |
| Band-2 hit end-to-end | Grouped `0.25`-bit `log_p` emit (`ldur x4,[x_p,#-8]; ldr x5,[x_bp]; str x4,[x5],#8; str x5,[x_bp]` + branch/setup = 8.0) + 4x-unrolled 8B copy (`ldp, ldp, stp, stp, subs, b.ne` = 6/4 = 1.50) + 4x-unrolled `rev_apply_8b_u4` (`2 ldp` + `4 x (ldrb, add, strb)` + `subs, b.ne` = 16/4 = 4.00) | **13.50** / hit | `C_HIT_HI = 14.5` |
| Band-1a hit (`[2^24, 2^48)`) | `m`-partitioned 4-bank `a`-tables (`acc_m[(m-m_b)&3]` once/bucket): 4x-unrolled `rev_apply_8b_u4` (4.00) + 4x-unrolled 8B `ldp/stp` copy (1.50) | **5.50** / hit | `C_HIT_LO_M4 = 6.25` |
| `rev_b1_early4_find` (`[2^48, 2^55)`) | Empty early-exit `keys[0] > neg_qb`: pointer/bucket load (`ldr, add, ldr, cmp, b.gt` = 5) + header/return (6) = 11 Ir when `cnt = 0`; `+14 + cnt` 4-strided scan when `cnt > 0`; binary-search fallback `24 + 8 log2(mu + 1)` | **11** empty | exact 16D integral |
| Band-1b hit (`[2^48, 2^55)`) | Contiguous 4-bank (`4 x 2^S`) 4x-unrolled `rev_b1_apply_u4`: `3 ldp` per 4x12B (0.75/hit) + per-hit `sub w1, w0, w_mb_sh; and w1, w1, #mask; ldrb w2, [x_acc, x1, uxtw]; add w2, w2, w_lp; strb w2, [x_acc, x1, uxtw]` (5.00) + loop (0.25) = 6.00 + 12B 4x-unrolled copy (`3 ldp + 3 stp + subs + b.ne` = 7/4 = 1.75) | **7.75** / hit | `C_HIT_LO_U = 8.0` |
| `rev_build_pass1/2/3` | 3-pass `a_0`-chord + `perm[]` build: Pass 1 histogram (`ldr, lsr, ldr, add, str, subs, b.ne` = 7) + Pass 2 scatter (`ldp, lsr, ldr, add, str, stp, str, add, cmp, b.lo` + reg = 13) + Pass 3 `q_key` pack & min-track (14) + setup (1) = 35; per-SB `sb_min_q8` & prefix (`ldrb, cmp, csel, strb` = 4) | **35** / entry (**4** / SB) | `C_LIST_A = 35.0`, `C_LBKT = 4.0` |
| `rev_scan_zero_u4` | 4x-unrolled 64-bit SWAR borrow-free unsigned byte-`>= T_g(b)` check (`and lo, x, #0x7f..; add t, lo, bias; orr/and hi; tst #0x80..` + `ldp` + `stp xzr, xzr` zeroing both `S_g, S_f`, <= 2.75 Ir/cell) x `1.0276` tile pad + `0.10` setup | **2.9259** / kept cell | `C_SCAN_ZERO = 2.95` |

Exactness of 4-orbit `sigma`-shared replay with `D_(omega,sigma)` (`239/239` exact-1), per-SB `q_min_j`, and 1-byte `sb_min_q8` gate is verified live in `audit_cm.py` (`sigmaB_exactness`: 48 real BLS12-381 `w' = omega^2` ideals, `0` mismatches).

## 2. Exact Cyclotomic Tower Norm Kernels

Tower `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q` with 3-mul Eisenstein Karatsuba and 2-mul quadratic relative norms (fixed-limb scalar ARMv8-A, zero division, zero SIMD):
- **`_norm_z15_2limb` (`3,357` Ir):** Degree-2 relative norm `Q(zeta_15) -> Q(omega, sqrt(5))` (3 Karatsuba Eisenstein-quadratic 2-limb muls at `8` Ir/2-limb-mul + adds = `320` Ir) + degree-2 norm `Q(omega, sqrt(5)) -> Q(omega)` on 4-limb intermediates (`900` Ir) + Eisenstein norm `c_0^2 - c_0 c_1 + c_1^2` on 8-limb intermediates (`3 x 64` 64x64 `mul/umulh/adds/adcs` = `1,152` Ir + carries `448` Ir = `1,600` Ir) + load/store/epilogue (`537` Ir) = **`3,357` Ir**.
- **`_nk_g_only` (`4,204 <= 4,500` Ir):** `a^3, b^3` in `O_H` + `a^3 - u omega^2 b^3` (`847` Ir) + `_norm_z15_2limb` (`3,357` Ir) = **`4,204` Ir** (charged **`4,500`**).
- **`_nk_f_only` (`6,167 <= 6,500` Ir):** `a^3, b^3, a^6, b^6, a^9, b^9` Karatsuba in `O_H` (`2,810` Ir) + `_norm_z15_2limb` (`3,357` Ir) = **`6,167` Ir** (charged **`6,500`**).

## 3. Closed-Form Non-Recursive 4-Prime 62-Bit Radix-2 NTT (`fft_F`), Analytical Stage B, & 2D Cofactorization

- **Closed-form non-recursive 4-prime 62-bit radix-2 NTT (`cm_models.fft_F`, ZERO recursion & ZERO FFT heuristic):**
  Fix the 4 explicit 62-bit NTT primes `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`2422936599525326849`, `2485986994308513793`, `2747195772696002561`, `2783224569714966529`; all 4 verified prime LIVE in `audit_cm.py` with `log2(NTT4_PROD) = 244.70`). Every `N`-limb (`64`-bit) product uses ONE of two exact non-recursive branches (preconditions asserted live inside `fft_F(N)`):
  1. **1D branch (`2 N <= 2^53`):** Pad to `M = 2^ceil(log2(2 N)) <= 2^53` (asserts `N * 2^128 < NTT4_PROD`): `4` primes x `3` radix-2 NTTs x `(M/2) log2(M)` scalar Barrett butterflies (`ldp`, `mul/umulh/msub`, `add/sub/cmp/csel`, `stp` = 26 Ir, charged **28** Ir) = **`168 M log2(M)`** + `4 M` pointwise muls (`80 M`) + `M` 4-prime Garner CRT/carry (`140 M`) + `1200` setup + four-step 2-pass flat-memory transpose past `2^26` limbs (`+504 M`):
     `fft_F(N) = 168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`.
  2. **2D polynomial convolution branch (`2 N > 2^53`, up to `N <= 2^100`):** Split `N` into `K_2 = ceil(N / K_1)` rows of width `K_1 = 2^(ceil(log2(N)) // 2)`, zero-pad both axes by `2x` to `M_1 = 2 K_1 <= 2^53` and `M_2 = 2^ceil(log2(2 K_2)) <= 2^53` (asserts `K_1 K_2 * 2^128 < NTT4_PROD`). With `M_tot = M_1 M_2`: `168 M_tot log2(M_tot)` butterflies + `260 M_tot + 1200` pointwise/Garner/overlap-add + `336 M_tot` 2D transposes + `504 M_tot` flat-memory streaming:
     `fft_F(N) = 168 M_tot log2(M_tot) + 1100 M_tot + 1200`.
  `audit_cm.py` (`srt_toy`) verifies exact scaled-remainder-tree remaindering on 128 planted leaves (`0` mismatches).
- **Analytical Stage B (`cm_models.stage_b_analytic`):** For an `nb`-bit survivor (`L = ceil(nb / 64)` limbs), `stage_b_analytic(nb) = gbatch_cost(24197143, L) + 10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16)` (`MONT5 = 367.76` Ir for 5-limb CIOS Montgomery mul: `50` `mul/umulh` + `150` `adds/adcs` + `60` `ldr/str` + conditional sub + loop). At `550`b (`L = 9`): `5,721,911` Ir; at `700`b (`L = 11`): `8,479,620` Ir.
- **Live 2D-Histogram Cofactorization:** Integrated LIVE across all `337` bins up to `bg_max = 700`b via `cm_models.Retention` + `stage_b_analytic` at `(s_a_g, s_b_g) = (345, 295)`, `(s_a_f, s_b_f) = (405, 345)`, giving strict bin-wise **`c_cof_2d_upper = 3.5980`** ipc.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (Straight-Line 4-Limb Scalar ARMv8-A)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`):
- **`unit_run5_z12` (`7.80` Ir/vec-entry):** Add/sub a 4-limb vector entry into a 5-limb unreduced accumulator: `2 ldp` (2) + `adds/adcs/adcs/adcs/adc` (5) + loop control amortized over run (0.8) = **`7.80` Ir**.
- **`gen_run1_z12` (`25.0 <= 33.0` Ir/vec-entry):** Scalar 1-limb coeff multiply-accumulate into 5-limb acc: `2 ldp` (2) + `4 x (mul, umulh, adds, adcs)` (16) + carry/loop (7) = **`25.0` Ir** (charged **`33.0`**).
- **`chi_mul_z12_fast` (`50` Ir) & `row_red5_z12` (`53` Ir):** Cyclic `zeta_12^k` automorphism is pure 4-limb sign-flips + at most one 4-limb add/sub (`50` Ir, zero `mul`); 5-limb to 4-limb Barrett/shift-add mod-`r` reduction is `53` Ir.
- **`row_tail_6bucket_z12` (`504` Ir) & 1-Owner Arch B Row Overhead (`785` Ir):** Combining 6 `zeta_12`-orbit buckets (`5 x 50` `chi_mul` + `5 x 40` 4-limb `add_z12` + `54` final reduction = `504` Ir) + row setup/load/store (`265` Ir) + loop (`16` Ir) = **`785` Ir** worst-case (`0` cross-owner sends/combines; `k5_cap = 10.50` vs `2^23` measured `9.5620`; `red_row_u_z12 = 18.6` Ir/elem).

## 5. Descent Kernels (`16D` Rank-2 `dsieve16` & Limb-Scaled 100-Round Miller-Rabin)

Analytical baseline scalar ARMv8-A counts (`1.5x` safety factor): `SIEVE_HIT = 1,650`, `SCAN + ZERO = 18`, `BOUND_CELL = 3,600`, `PAIR_SETUP = 1.5 x 1.33 x 2^25.1`, `LLL_NODE = 2^38` (16D LLL on 381-bit entries), `NORMS = 2^27`, `TDIV = 2^25`, `IDEAL = 2^35`, `ecm_instr` (`MONT5 = 367.76`), and exact limb-scaled 100-round Miller-Rabin **`mr_instr(nb) = 100 * 2 * nb * MONT5 * (ceil(nb/64)/5)^2`** (`<= 4^(-100) = 2^(-200)` error/test, `<= 2^(-80)` over `< 2^120` run calls).

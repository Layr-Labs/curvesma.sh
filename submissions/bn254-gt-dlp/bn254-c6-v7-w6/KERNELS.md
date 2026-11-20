# `KERNELS.md` — Pure-Integer Baseline Scalar AArch64 Kernels (`simd = 0`, Zero `fp64`/`LSE`/`DC ZVA`, Zero External C)

All inner loops use **100% pure baseline scalar ARMv8.0-A integer instructions** (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero external C/Callgrind measurements — every kernel is derived by explicit opcode count and verified live in `audit_bn_cm.py`).

## 0. Exact Pure-Integer Geometry & 2-Stage Pure-Int Bareiss 12D LLL (`zero_float_geom = True`)

1. **Exact `sigma`-Invariant `G_sym` & Voronoi `D_fi` (`canon <= 2^74.03`, `330` Ir/raw `a` `<< 2^42.8`):** `G_sym = G_0 + S^T G_0 S + (S^2)^T G_0 S^2` (`G_SC = 196,608`, `S^T G_sym S == G_sym`), `Q(v) < 2^44 << 2^63`, `D_fi = U_fi^T G_sym U_fi - G_sym`. Certified for ALL radii LIVE (`[exhaustive_fd]`): Dedekind `O_lat == O_h`, exact `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), outward-rounded `Decimal(55)` index-1 (`unit_index1_proved = True`) and convex-hexagon (`min_lo = +0.449890 > 0`).
2. **Pure Product Region (`kap = 0.0`):** `Q(a) <= T_a` AND `Q(b) <= T_b` (`frac_reg = 1.0`, ZERO `M_t[t]`, ZERO `Q(a)` sort, ZERO `sb_u7`).
3. **2-Stage Pure-Int Bareiss 12D LLL (`11/40 -> 9/10`, `<= 2^40` Boot, `<= 2^36` Special-`q` up to `q <= 2^220`):** `inductive_sr_inv_ok = True`, `exact_div_errors = 0`, `S_1 + S_2 <= S_max(D_0) = 185,076` -> joint 2-stage worst-case **`2^{39.81} < 2^{40}` Ir** on Boot (`two_stage_ok = True`, measured `2^{29.48}`, `2,380` swaps) and **`2^{32.88} << 2^{36}` Ir** on Special-`q` up to `q <= 2^{220}` (`two_stage_ok = True` live on `sq162_g = 2^{31.97}`, `sq213_f = 2^{32.79}`, `sq220_wc = 2^{32.88}`, `boot[0]`, `24/24` prod special-`q`, `8/8` boot certs, `100` deg-8 + `224` depth-2 ideals).

## 1. Pure-Product Deterministic-Occ Band-1/2 Direct-Apply Sieve, Zero-Sort 2-Pass Lists, & 64-Cell `ORR`

Each `(b_1..b_5)` row worker owns one unified `(L_max + 2) * N_a * 2 B = 2^{56.18} B << 2^{64} B` `uint16_t` arena, streaming one read-only `a`-table at a time (`<= 2^{56.21} B << 2^{64} B`):

1. **Band 1 `[2^14, 2^45)` Direct-Apply (`APPLY8 = 30/8 = 3.75 == C_HIT_LO`, `C_B1_ROW = 24.0`):** `ldp` bounds (`1`) + 8x-unrolled `APPLY8` into `x_row = &score[t][0]`: `4 ldp` + `8*(ldrh; add; strh)` (`24`) + `subs; b.ne` (`2`) = **`30 / 8 = 3.75` exact opcodes/hit**; `C_B1_ROW = 24.0` covers setup (`<= 6`) + worst-case `< 8` tail (`<= 7*(6-3.75) = 15.75`; `21.75 <= 24.0`; `0` mismatches live).
2. **Band 2 `[2^45, 2^54.5)` 1-Bit `occ_bm[]` & Deterministic Occupied-SB Bound (`W_SB = 128`, `max_sb = 8 <= 8`, ZERO Band 3, ZERO Poisson):**
   - **Deterministic Occupied-SB Bound:** For ANY prime `p` (`S_p = ceil(p / W_SB)` super-buckets, `N_a` `a`-entries, arbitrary residue structure), `#(non-empty SBs) <= min(S_p, N_a)` deterministically (`det_occ_bound_ok = True` live on all 24 prototype ideals). Thus occupied-SB density `<= min(1.0, N_a / S_p) = min(1.0, mu_sb)`; under `H2` equidistribution of `b`-row starting residues mod `p`, entered-SB rate `<= min(1.0, mu_sb)` (used in `compute_ledger`).
   - **`b2_occ_probe` (`12 == C_SWAR_PROBE` exact ops/row):**
     ```asm
     lsr x_sb,x_k0,#7; lsr x_by,x_sb,#3; and x_sh,x_sb,#7
     ldrh w0,[x_bm,x_by]; lsr w0,w0,w_sh; tst w0,w_mask; b.ne .Lenter
     add x_k0,x_k0,x_step; subs x_t,x_k0,x_p; csel x_k0,x_t,x_k0,hs
     subs x_cnt,x_cnt,#1; b.ne .Lprobe
     ```
     (`C_SWAR_PROBE = 22.0` sensitivity: `+0.025`b).
   - **`b2_sb_enter` (`12 <= C_SB_ENTER = 14.0`):** `rbit; clz; add; add+ldp; sub/cmp/csel; bic` (`10`) + `2` exit = `12 <= 14.0`.
   - **`b2_cand_prod` (`7 == C_CAND` exact ops/cand; exact `u64` wrap, `0` errs live):**
     ```asm
     cmp x_p,x_end; b.hs .Ldone; ldr x0,[x_p],#8
     and x1,x0,#0x7f; add x2,x1,x_off; cmp x2,x_Lb; b.hs .Lcand
     ```
   - **`b2_hit_direct` (`6 <= C_HIT_HI = 7.0` exact ops/hit):**
     ```asm
     lsr x_a,x0,#8; ldr x_row,[x_rtab,x2,lsl #3]
     ldrh w3,[x_row,x_a,lsl #1]; add w3,w3,w_lp; strh w3,[x_row,x_a,lsl #1]; b .Lcand
     ```
3. **Zero-Sort Sequential 2-Pass `a_0`-Run Lists (`C_LIST_B1 = 20.0, C_LIST_B2 = 24.0, C_LBUCKET = 6.0`):** ZERO `Q(a)` sort, ZERO Pass 3. Pass 1 (`b2_pass1_seq`): `14` ops + `<= 0.5` line setup (`40/89 = 0.45`) = `14.5` (Band-1 `12.0`); prefix sum + in-place `occ_bm[]` & `cnt[]` clear: `<= 5.5 <= 6.0` Ir/SB; Pass 2 (`b2_scatter_seq`): `<= 8.5` (Band-1 `7.5`) -> `23.0 <= 24.0`, `19.5 <= 20.0` (`zero_sort_exact_ok = True` live).
4. **Unified 64-Cell `ORR` Scan+Zero (`59/64 = 0.9219 <= C_SCAN_ZERO = 0.96`, ZERO `DC ZVA`):** `8 ldp` (scores) + `8 stp xzr,xzr` (zero) + `8 ldp` (bias) + `16 add` + `15 orr` + `tst; b.ne` + `subs; b.ne` = **`59/64 = 0.9219 <= 0.96`** (25 GPRs, zero spills; `0 / 128,000` errs live).

## 2. 18-Op Radix-2 `< 2^62` 3-Prime Shoup NTT (`ntt3_F(M_prod) = 84*n*k + 260*n`, ZERO `H-FFT`)

Three `< 2^62` primes divisible by `2^55`: $p_1 = 116 \cdot 2^{55} + 1$, $p_2 = 69 \cdot 2^{55} + 1$, $p_3 = 57 \cdot 2^{55} + 1$ ($\prod p_i = 2^{183.799} > (n/2)(2^{64}-1)^2$ for $n \le 2^{55}$).
- **18-Op Radix-2 Shoup Butterfly (`9` ops/elem/stage -> `81*n*k <= 84.0*n*k` across 3 primes × 3 transforms):**
  ```asm
  ldr x_u,[x_a]; ldr x_v,[x_b]; ldp x_w,x_wp,[x_tw],#16        ; 3 ops
  umulh x_q,x_v,x_wp; mul x_t,x_v,x_w; msub x_t,x_q,x_p,x_t    ; 3 ops (t in [0, 2p), 4p < 2^64)
  subs x_t2,x_t,x_p; csel x_t,x_t2,x_t,hs                      ; 2 ops (t in [0, p))
  add x_au,x_u,x_t; subs x_a2,x_au,x_p; csel x_au,x_a2,x_au,hs ; 3 ops
  subs x_bv,x_u,x_t; add x_b2,x_bv,x_p; csel x_bv,x_bv,x_b2,hs ; 3 ops
  str x_au,[x_a],#8; str x_bv,[x_b],#8; subs x_c,x_c,#1; b.ne  ; 4 ops -> total 18 ops/butterfly
  ```
  Across 3 primes × 3 transforms × `(n/2)*k` butterflies: `3 * 3 * 9 * n * k = 81.0 * n * k <= 84.0 * n * k`.
- **Linear + Garner CRT:** Input reduce (`60`) + pointwise (`30`) + inv scale (`30`) + Garner CRT (`80`) + carry (`30`) = `230 <= 260` ops/limb (`+ 512*n` if `n > 2^{26}`). Live iterative Radix-2 Shoup butterfly NTT convolution (`n = 16`, 4 stages) + Garner CRT verified bit-exact (`ntt3_shoup_butterfly_ok = True`, `ntt3_garner_crt_ok = True`, `ntt3_r2_shoup_conv_ok = True`).

## 3. Stage B Zero-GCD Granlund-Montgomery `P_14` Strip (`C_SB = 56,000`), `7`-Limb `f`-SRT, & 4-Vec `Z[omega]` LA

- **Deterministic Worst-Case Zero-GCD `P_14` Strip (`55,949 <= C_SB = 56,000` for ANY `N_g < 2^{512}`, live `p14_zero_gcd_gm_ok = True`, `0/32` errs incl. prime cubes):**
  Pack all `1,899` odd primes `p in [3, 2^{14})` into `417` `< 2^{64}` products `P_c`; normalize `P'_c = P_c << clz(P_c) in [2^{63}, 2^{64})`, precompute Möller–Granlund `mu_c = floor((2^{128}-1)/P'_c) - 2^{64}`, `R_{c,i} = 2^{64i} mod P'_c` (`i=1..7`), and `(inv_p = p^{-1} mod 2^{64}, lim_p = floor((2^{64}-1)/p))`:
  1. **Power-of-2 strip (`rbit; clz` + 8-limb shift):** `<= 60` ops.
  2. **`417` chunk reductions mod `P'_c` (`70 * 417 = 29,190` ops):** load `R_{1..7}, P'_c, mu_c` (`4 ldp + 1 ldr = 5`), 8-limb MAC into `(s_0, s_1, s_2)` (`7*5 = 35`), fold `s_1*R_1 + s_2*R_2 + t_2` & `t_1 >= P'_c` (`10+4+2 = 16`, guaranteeing `0 <= t_1 < P'_c`), branchless 12-op Möller–Granlund 2-by-1 `udiv_rem_2by1(t_1, t_0, P'_c, mu_c)` (`mul; umulh; adds; adc; add #1; mul; sub; cmp; add; csel; subs; csel`), `+2` loop = `70` ops/chunk.
  3. **Granlund–Montgomery tests (`1,899 * 5 = 9,495` ops):** since `p | P_c` is odd, `p | N_g <=> p | rem <=> ((rem * inv_p) mod 2^{64}) <= lim_p` (`ldp; mul; cmp; b.ls` + ptr = `5`).
  4. **Deterministic Worst-Case 8-Limb Hensel Divisions (`<= 391 * 44 = 17,204` ops):** for ANY `N_g < 2^{512}`, `omega(N_g) <= 68` (`prod_{i=1}^{69} p_i > 2^{512}`) and `Omega(N_g) <= floor(512/log2(3)) = 323`, giving `<= 391` passes (`8*(subs; cset; mul; umulh; add) + 4 = 44`; `borrow == 0` iff exact).
  - **Total:** `60 + 29,190 + 9,495 + 17,204 = 55,949 <= 56,000` ops (ZERO GCDs!).
- **`7`-Limb `f`-SRT & 4-Vec `Z[omega]` BW (ZERO TABM, live `[live_z6_la_kernels]`):** `C_NORMS = 30k`, `4` disjoint LP windows (`0` fn), AM-GM `sup log2|N_f| <= 441.66 < 448`b = 7 limbs (`+6.34`b, `Q_F_SAFETY = 2.0`). `unit4_reg_z6`: `ldr; add` (`2`) + `4*(2 ldp + 5 adds/adcs)` (`28`) + `subs; b.ne` (`2`) = `32/4 = 8.00 == C_UNIT_5` (20 GPRs `x8..x27`); `gen1_reg_z6 = 56 <= 58.0`; `row_tail_3bucket_z6 = 722 <= 740.0` (1-step Barrett `mu_r = floor(2^317/r)`).

## 4. Descent `MONT5` 2-Addend CIOS & 6-Bit Geometric ECM Ladder

- **2-Addend CIOS (`335 <= 345 <= MONT5 = 367.76` at `L = 5`, `all_L_5_100_below_budget = True`):** `32/32` exact 320-bit CIOS live.
- **6-Bit ECM Ladder Upper-Endpoint Separation Audit:** Exact shared `_ecm_best` with upper-endpoint `[lo, hi_end]` enclosure in `_ss_ecm(..., hi=True)` across ALL **266 DAG rungs** (`plan59`, `T0` & `T`) AND ALL **16 Boot rungs** (`64.42..156.50`b, `tau = 1250`, `kA_b + 1 = 17`) verifies LIVE **`min_dag = 0.6174 > 0.50`** and **`min_boot = 0.5511 > 0.50 = SEPARATION`** (`63/63 = 100%` live Suyama-12 isolation) + `(kA + kB)` Miller-Rabin `PRIME` checks.

# `KERNELS.md` — Pure-Integer Baseline Scalar AArch64 Kernels (`simd = 0`, Zero `fp64`/`LSE`/`DC ZVA`/External C)

All inner loops use **100% pure baseline scalar ARMv8.0-A integer instructions** (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero external C/Callgrind measurements — derived by explicit opcode count and verified live in `audit_bn_cm.py`).

## 0. Exact Pure-Int Geometry & 2-Stage Pure-Int Bareiss 12D LLL (`zero_float_geom = True`)

1. **`sigma`-Invariant `G_sym` & Voronoi `D_fi` (`canon <= 2^74.03`, `330` Ir/raw `a`):** `G_sym = G_0 + S^T G_0 S + (S^2)^T G_0 S^2` (`G_SC = 196,608`), `Q(v) < 2^44 << 2^63`, `D_fi = U_fi^T G_sym U_fi - G_sym`. Certified LIVE (`[exhaustive_fd]`): Dedekind `O_lat == O_h`, `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), outward-rounded `Decimal(55)` index-1 (`unit_index1_proved = True`) & convex-hexagon (`min_lo = +0.449890 > 0`).
2. **Pure Product (`kap = 0.0`):** `Q(a) <= T_a` AND `Q(b) <= T_b` (`frac_reg = 1.0`, ZERO `M_t[t]`, ZERO `Q(a)` sort, ZERO `sb_u7`).
3. **2-Stage Bareiss 12D LLL (`11/40 -> 9/10`, `<= 2^40` Boot, `<= 2^36` Special-`q` up to `q <= 2^220`):** `inductive_sr_inv_ok = True`, `0` div errs, `S_1 + S_2 <= S_max(D_0) = 185,076` -> **`2^{39.81} < 2^{40}` Ir** on Boot (`2,380` swaps, measured `2^{29.48}`) and **`2^{32.88} << 2^{36}` Ir** on Special-`q` (`sq162_g`, `sq213_f`, `sq220_wc`, `boot[0]`, `24/24` prod special-`q`, `8/8` boot certs, `100` deg-8 + `224` depth-2 ideals).

## 1. Pure-Product Deterministic-Occ Band-1/2 Direct-Apply Sieve, Zero-Sort 2-Pass Lists, & 64-Cell `ORR`

Each `(b_1..b_5)` worker owns one `(L_max + 2)*N_a*2 B = 2^{56.18} B << 2^{64} B` `uint16_t` arena, streaming one `a`-table (`<= 2^{56.21} B`):
1. **Band 1 `[2^14, 2^45)` Direct-Apply (`APPLY8 = 30/8 = 3.75`, `C_B1_ROW = 24.0`):** `4 ldp` + `8*(ldrh; add; strh)` + `subs; b.ne` = **`30/8 = 3.75`** ops/hit; `C_B1_ROW = 24.0` covers setup (`<= 6`) + `< 8` tail (`<= 15.75`; `21.75 <= 24.0`; `0` errs live).
2. **Band 2 `[2^45, 2^54.5)` 1-Bit `occ_bm[]` & Deterministic `p_occ = min(1.0, mu_sb)` (`W_SB = 128, max_sb = 8`, ZERO Band 3/Poisson):**
   - For ANY prime `p` (`S_p = ceil(p/W_SB)` SBs, `N_a` entries), `#(non-empty SBs) <= min(S_p, N_a)` deterministically (`det_occ_bound_ok = True`), so under `H2` entered-SB rate `<= min(1.0, mu_sb)`.
   - `b2_occ_probe` (`12 == C_SWAR_PROBE`; `22` sensitivity `+0.023`b): `lsr; lsr; and; ldrh; lsr; tst; b.ne; add; subs; csel; subs; b.ne`.
   - `b2_sb_enter` (`12 <= 14.0`): `rbit; clz; add; add+ldp; sub/cmp/csel; bic` (`10`) + `2` exit.
   - `b2_cand_prod` (`7 == C_CAND`, `0` wrap errs live): `cmp; b.hs; ldr x0,[xp],#8; and x1,x0,#0x7f; add x2,x1,x_off; cmp x2,x_Lb; b.hs`.
   - `b2_hit_direct` (`6 <= 7.0`): `lsr x_a,x0,#8; ldr x_row,[x_rtab,x2,lsl #3]; ldrh w3,[x_row,x_a,lsl #1]; add; strh; b`.
3. **Zero-Sort Sequential 2-Pass Lists (`C_LIST_B1=20, C_LIST_B2=24, C_LBUCKET=6`, `zero_sort_exact_ok = True`):** Pass 1 `14 + 0.5 = 14.5` (B1 `12`), prefix+clear `<= 5.5 <= 6`, Pass 2 `<= 8.5` (B1 `7.5`) -> `23 <= 24`, `19.5 <= 20`.
4. **Unified 64-Cell `ORR` Scan+Zero (`59/64 = 0.9219 <= 0.96`, ZERO `DC ZVA`):** `8 ldp` + `8 stp xzr,xzr` + `8 ldp` + `16 add` + `15 orr` + `tst; b.ne` + `subs; b.ne` = **`59/64`** (`25` GPRs, `0/128,000` errs live).

## 2. 48-Op Radix-4 + 43-Op Radix-3 `< 2^62` 3-Prime Lazy `[0, 2p)` Shoup NTT (ZERO `H-FFT`)

Three `< 2^62` primes divisible by `3 * 2^{55}`: $p_1 = 69\cdot 2^{55}+1$, $p_2 = 57\cdot 2^{55}+1$, $p_3 = 54\cdot 2^{55}+1$ ($\prod p_i = 2^{182.696}$). All values stay in `[0, 2p)` (`4p < 2^{64}`), so 3-op `(umulh; mul; msub)` maps `v in [0, 2p)` to `v*w mod p in [0, 2p)`; `J = g^{(p-1)/4} mod p` (`(J, J_p)` in regs).
- **48-Op Radix-4 Lazy Shoup Butterfly (`54.0*n*k <= 60.0*n*k`, `+11.1%` headroom):**
  1. Loads (`7`): `4 ldr` (`x0..x3`) + `3 ldp` (`w1..w3, wp1..wp3`).
  2. Twiddles (`9`): `3*(umulh; mul; msub)` -> `t1,t2,t3 in [0, 2p)`.
  3. Stage 1 (`14`): `u0=(x0+t2) mod 2p`, `u1=(x0-t2) mod 2p`, `u2=(t1+t3) mod 2p` (`3*(add/subs; subs/add; csel) = 9`), `d = t1-t3+2p in [0, 4p)` (`add; sub = 2`), `u3 = J*d mod p in [0, 2p)` (`umulh; mul; msub = 3`, valid since `4p < 2^{64}`).
  4. Stage 2 (`12`): `y0=u0+u2`, `y2=u0-u2`, `y1=u1+u3`, `y3=u1-u3 mod 2p`: `4*(add/subs; subs/add; csel)`.
  5. Stores/loop (`6`): `4 str` + `subs; b.ne`.
  - Total: `7+9+14+12+6 = 48` ops (`22` GPRs, 0 spills) -> `3*3*48*(n/4)*(k/2) = 54.0*n*k <= 60.0*n*k`.
- **43-Op Radix-3 Layer (`<= 45` -> `135*n_3`):** `3 ldr + 2 ldp` (`5`) + `6*(umulh; mul; msub)` (`18`) + `5*(add; subs; csel)` (`15`) + `3 str + subs; b.ne` (`5`) = `43 <= 45` -> `ntt3_F = min(60*n*k + 220*n, 60*n_3*(k-2) + 355*n_3)` (`+512*n` if `n > 2^{26}`). Linear+Garner CRT: `45+30+30+85+30 = 220`. Live iterative `n=16` (Radix-4) & `n=24` (Radix-4+3) convolutions verified bit-exact (`ntt3_r4_r3_shoup_conv_ok = True`).

## 3. Zero-GCD `P_14` Strip (`C_SB = 56,000`), `7`-Limb `f`-SRT, & 4-Vec `Z[omega]` LA

- **Deterministic Zero-GCD `P_14` Strip (`55,949 <= 56,000` for ANY `N_g < 2^{512}`, `0/20` errs incl. prime cubes, `p14_zero_gcd_gm_ok = True`):** Pack `1,899` odd primes `< 2^{14}` into `417` `< 2^{64}` chunks `P_c` (`P'_c = P_c << clz(P_c) in [2^{63}, 2^{64})`, Möller–Granlund `mu_c`, `R_{c,1..7}`, Granlund–Montgomery `(inv_p, lim_p)`): (1) power-of-2 strip `<= 60`; (2) `417` chunks × `70` (`5` loads + `35` MAC + `16` fold to `t_1 < P'_c` + `12` branchless 2-by-1 `udiv_rem` + `2` loop) = `29,190`; (3) `1,899 * 5 = 9,495` GM tests (`(rem*inv_p) mod 2^{64} <= lim_p`); (4) worst-case 8-limb Hensel divisions (`omega(N_g) <= 68, Omega(N_g) <= 323` -> `<= 391 * 44 = 17,204`) -> total **`55,949 <= 56,000`** (ZERO GCDs!).
- **`7`-Limb `f`-SRT & 4-Vec `Z[omega]` BW (ZERO TABM):** `C_NORMS = 30k`, `4` disjoint windows (`0` fn), `sup log2|N_f| <= 441.66 < 448`b (`Q_F_SAFETY = 2.0`). `unit4_reg_z6 = 32/4 = 8.00` (20 GPRs `x8..x27`), `gen1_reg_z6 = 56 <= 58`, `row_tail_3bucket_z6 = 722 <= 740` (1-step Barrett `mu_r = floor(2^{317}/r)`).

## 4. Descent `MONT5`, Exact-Int `PAIR_SETUP`/`dsieve`, & 6-Bit ECM Ladder

- **2-Addend CIOS (`335 <= 345 <= MONT5 = 367.76`, `all_L_5_100_below_budget = True`):** `32/32` exact 320-bit CIOS live.
- **Exact-Int `dsieve` & `PAIR_SETUP` (`descent_bn254_v6.py`):** Sieve nodes charge `CELLS[R]*(HITS_PER_CELL*SIEVE_HIT + SCAN + ZERO + BOUND_CELL) + pairs*PAIR_SETUP` (`HITS_PER_CELL = 2*ln(59/8) = 3.996`, `SIEVE_HIT = 1,100` Ir, `BOUND_CELL = 2,400` Ir, `SCAN = 10`, `ZERO = 2`, `PAIR_SETUP = 2^{25.1} = 35.96M` Ir per `(p,r)` pair for `p <= 2^{16}`), plus `LLL_NODE = 2^{36}` (`2^{32.88}` live), `NORMS = 2^{26}`, `TDIV = 2^{32}`, `IDEAL = 2^{34}`, and `+ 2^{40}` control.
- **6-Bit ECM Ladder:** Upper-endpoint `[lo, hi_end]` enclosure across ALL **266 DAG rungs** and **15 Boot rungs** (`64.42..155.16`b, `tau = 1250`, `kA_b + 1 = 13`) verifies LIVE **`min_dag = 0.6174 > 0.50`** and **`min_boot = 0.6081 > 0.50`** (`63/63` Suyama-12 live) + `(kA + kB)` MR `PRIME` checks.

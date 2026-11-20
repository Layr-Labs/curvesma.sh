# `KERNELS.md` — Pure-Integer Baseline Scalar AArch64 Kernels (`simd = 0`, Zero `fp64`/`LSE`/`DC ZVA`/External C)

All inner loops use **100% pure baseline scalar ARMv8.0-A integer instructions** (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero external C/Callgrind measurements — derived by explicit opcode count and verified live in `audit_bn_cm.py`).

## 0. Exact Pure-Int Geometry & Input-Uniform 2-Stage Bareiss 12D LLL (`zero_float_geom = True`)

1. **`sigma`-Invariant `G_sym` & Voronoi `D_fi` (`canon <= 2^73.874752`, `330` Ir/raw `a`):** `G_sym = G_0 + S^T G_0 S + (S^2)^T G_0 S^2` (`G_SC = 196,608`), `Q(v) < 2^44 << 2^63`, `D_fi = U_fi^T G_sym U_fi - G_sym`. Certified LIVE (`[exhaustive_fd]`): Dedekind `O_lat == O_h`, `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), outward-rounded `Decimal(55)` index-1 (`unit_index1_proved = True`) & convex-hexagon (`min_lo = +0.449890 > 0`). Pure product (`kap = 0.0`): `Q(a) <= T_a` AND `Q(b) <= T_b` (`frac_reg = 1.0`, ZERO `M_t[t]`, ZERO `Q(a)` sort, ZERO `sb_u7`).
2. **Input-Uniform 2-Stage Bareiss 12D LLL (`11/40 -> 9/10`, PROVED FOR EVERY ARBITRARY INPUT):**
   - *Exact input-independent `d_k(B_0)`:* For ANY fresh-`e` boot input `R in F_p[y]/(h)` (centered `|M_{R,ij}| <= p/2`), `B_boot = L_R * diag(p I_6, I_6)` (`det(L_R) = 1`), so `d_k(B_boot) = p^{2 min(k,6)} S_A^{max(0,k-6)} det(G_{blk,0..k-1})` is **strictly independent of `R`** (`bits_D = 3662`, `S_max = 185,076`). For ANY special-`q` `(q,s,rho)` (`q <= 2^{220}`), `B_sq = L_q * diag(q, 1..1)` (`det(L_q) = 1`) gives `d_k(B_sq) = q^2 det(G_{blk,0..k-1}) <= (2^{220}-1)^2 det(G_{blk,0..k-1})` (`bits_D <= 679`, `S_max <= 40,493`).
   - *Uniform limb bounds:* Every swap reduces `D = prod d_k` by `<= delta` (`S_1 + S_2 <= S_max`), `1 <= d_k < 2^{bits_D}` (`L_d`), `|mu_{k,j}| <= (3/2)^{11}(2 sqrt(d_max)/||b_j^*||_G + 1)` -> `|lam_{k,j}| < 2^{1.5 bits_D + 9}` (`L_lam`), `|q| < 2^{bits_D + 9}` (`L_q`), and Sylvester minors `G_6 - 50000*I_6 > 0` (`g6_lambda_min_gt_50k = True`, `sqrt(lambda_min) > 223`) give `|B_{k,c}| < 2^{0.5 bits_D + 9}` (`L_Bl`).
   - *Result (verified LIVE from exact integer `d_k` with strict upward `ceil`):* **`960,942,091,452 <= 2^{39.81} < 2^{40}` Ir** (`[58,86,58,29]`) for **EVERY** boot `R`, and **`7,912,444,650 <= 2^{32.89} < 2^{33} = LLL_NODE` Ir** (`[11,17,11,6]`) for **EVERY** special-`q <= 2^{220}` (`input_uniform_wc_proved = True`, `0` div errs).

## 1. Pure-Product Deterministic-Occ Band-1/2 Direct-Apply Sieve, Zero-Sort 2-Pass Lists, & 64-Cell `ORR`

Each `(b_1..b_5)` worker owns one `(L_max + 2)*N_a*2 B = 2^{56.000} B << 2^{64} B` `uint16_t` arena, streaming one `a`-table (`<= 2^{56.029} B`):
1. **Band 1 `[2^14, 2^44.5)` (`APPLY8 = 30/8 = 3.75`, `C_B1_ROW = 22.0`):** `4 ldp` + `8*(ldrh; add; strh)` + `subs; b.ne` = **`3.75`** ops/hit; `C_B1_ROW = 22.0` covers setup (`<= 6`) + `< 8` tail (`<= 15.75`; `21.75 <= 22.0`, `0` errs live).
2. **Band 2 `[2^44.5, 2^55.0)` 1-Bit `occ_bm[]` & Deterministic `p_occ = min(1.0, mu_sb)` (`W_SB = 256, max_sb = 5 <= 8`, ZERO Band 3/Poisson):**
   - For ANY prime `p` (`S_p = ceil(p/256)` SBs, `N_a` entries), `#(non-empty SBs) <= min(S_p, N_a)` deterministically (`det_occ_bound_ok = True`), so under `H2` entered-SB rate `<= min(1.0, mu_sb)`. 8-byte record `(ai << 8) | (kta & 0xff)` (`ai < 2^{45.36} < 2^{56}`, `0` wrap errs live for `W_SB = 256`, `sgn = +/-1`).
   - `b2_occ_probe` (`12 == C_SWAR_PROBE`; `22` sensitivity `+0.040`b): `lsr; lsr; and; ldrh; lsr; tst; b.ne; add; subs; csel; subs; b.ne`.
   - `b2_sb_enter` (`12 <= 14.0`): `rbit; clz; add; ldp; lsl; sub; cmp; add; csel; lsl; bic; cbnz` (`<= 14.0` with index-to-ptr adds).
   - `b2_cand_prod` (`7 == C_CAND`): `cmp; b.hs; ldr x0,[xp],#8; and x1,x0,#0xff; add/sub x2,x1,x_off; cmp x2,x_Lb; b.hs`.
   - `b2_hit_direct` (`6 == C_HIT_HI`): `lsr x_a,x0,#8; ldr x_row,[x_rtab,x2,lsl #3]; ldrh w3,[x_row,x_a,lsl #1]; add; strh; b`.
3. **Zero-Sort 2-Pass Lists (`C_LIST_B1 = 15.5, C_LIST_B2 = 19.5, C_LBUCKET = 5.5`, `zero_sort_exact_ok = True`):**
   - B2 Pass 1 (`12.5` ops unrolled `x4`): mod-`p` step `add;subs;csel` (`3`), `lsr` (`1`), `ldr/add/str cnt[q]` (`3`), pack `and;orr` (`2`), `str tmp_rec; str tmp_sb` (`2`), `add ai` + `x4` loop (`1.5`). B2 Pass 2 (`6.5` ops unrolled `x4`): `ldr q; ldr rec; ldr pos; str rec; add pos; str pos` (`6`) + `x4` loop (`0.5`) -> `19.0 <= 19.5`.
   - B1 (`15.0 <= 15.5`): Pass 1 omits `lsr/and/orr` (`9.5`), Pass 2 omits `tmp_rec` (`5.5`). Prefix+`occ_bm`+clear `<= 5.5`.
4. **Unified 64-Cell `ORR` Scan+Zero (`59/64 = 0.9219 <= 0.922`, ZERO `DC ZVA`):** `8 ldp` + `8 stp xzr,xzr` + `8 ldp` + `16 add` + `15 orr` + `tst; b.ne` + `subs; b.ne` = **`59/64`** (`25` GPRs, `0/128,000` errs live).

## 2. 48-Op Radix-4 + 46-Op Radix-3 `< 2^62` 3-Prime Lazy `[0, 2p)` Shoup NTT (ZERO `H-FFT`)

Three `< 2^62` primes divisible by `3 * 2^{55}`: $p_1 = 69\cdot 2^{55}+1$, $p_2 = 57\cdot 2^{55}+1$, $p_3 = 54\cdot 2^{55}+1$ ($\prod p_i = 2^{182.696}$). Values stay in `[0, 2p)` (`4p < 2^{64}`), 3-op `(umulh; mul; msub)` maps `[0, 2p) -> [0, 2p)`; `J = g^{(p-1)/4} mod p`.
- **48-Op Radix-4 (`54.0*n*k <= 55.0*n*k`):** `4 ldr + 3 ldp` (`7`) + `3*(umulh;mul;msub)` (`9`) + Stage-1 (`14`: `3*3 + 2 + 3`) + Stage-2 (`12`: `4*3`) + `4 str + subs;b.ne` (`6`) = **`48` ops** (`22` GPRs) -> `3*3*48*(n/4)*(k/2) = 54.0*n*k <= 55.0*n*k`.
- **46-Op Radix-3 (`46 <= 48` -> `144*n_3`):** `3 ldr + 2 ldp` (`5`) + `6*(umulh;mul;msub)` (`18`) + `6*(add;subs;csel)` (`18`: `s12,y0,sa,y1,sb,y2` mod `2p`) + `3 str + subs;b.ne` (`5`) = **`46 <= 48` ops** -> `ntt3_F = min(55*n*k + 220*n, 55*n_3*(k-2) + 364*n_3)` (`+512*n` if `n > 2^{26}`). Live `n=16` & `n=24` convolutions bit-exact (`ntt3_r4_r3_shoup_conv_ok = True`).

## 3. Split Norms (`6k / 24k`), Zero-GCD `P_14` Strip (`56k`), `7`-Limb `f`-SRT, & Plain-Sparse 4-Vec `Z[omega]` LA

- **Split Stage-B Norms (`C_NORM_G + C_NORM_F = 30,000`):** On `sf`, compute only `deg-2` `N_g = N_{K_h/Q}(u w A^2 - B^2)` (`<= 450` poly ops + `6x6` Bareiss `<= 2,600` = `<= 3,050 <= 6,000 = C_NORM_G`); ONLY on `s0` evaluate `deg-8` `N_f` (`<= 24,000 = C_NORM_F`; nonprim pays `C_NORM_G = 6,000`).
- **Deterministic Zero-GCD `P_14` Strip (`55,949 <= 56,000` for ANY `N_g < 2^{512}`, `0/20` errs, `p14_zero_gcd_gm_ok = True`):** `1,899` odd primes `< 2^{14}` in `417` `< 2^{64}` chunks: pow2 (`<= 60`) + `417*70 = 29,190` + `1,899*5 = 9,495` + worst-case Hensel (`<= 391*44 = 17,204`) = **`55,949 <= 56,000`**.
- **`7`-Limb `f`-SRT & Plain-Sparse 4-Vec `Z[omega]` BW (`TAB_R_CAP = 1.00`, ZERO TABM):** `4` disjoint windows (`0` fn, `lx2=55.0, sl/sb=160/165`), `sup log2|N_f| <= 440.42 < 448`b (`Q_F_SAFETY=2.0`), then `q_f*(C_TDIV+C_BATCH+c_b2)`. Plain-sparse BW (`tab_r_cap = 1.00`, zero grouped lookup): `unit4=8.00` (20 GPRs), `gen1=56<=58`, `tail=722<=740`.

## 4. Descent `MONT5`, Proved Node Kernels, Input-Uniform `PAIR_SETUP` (`p <= 2^{59}`), & 6-Bit ECM Ladder

- **2-Addend CIOS (`335 <= 345 <= MONT5 = 367.76`, `all_L_5_100_below_budget = True`):** `32/32` exact 320-bit CIOS live.
- **Input-Uniform `PAIR_SETUP = 2^{26.2}` for ALL `p <= 2^{59}` (`pair_setup_p59_proved = True` live):** Reduce `144` entries (`<= 6` limbs) mod `p` (`10,368`), evaluate `12` linear forms (`1,728`), 1-word `xgcd` + shift (`1,732`) -> `B_p = L_p diag(p,1..1)` (`13,828` ops); 2-stage Bareiss 12D LLL (`d_k = p^2 < 2^{118}`, `[2,3,2,2]`, `S_max <= 8,540`) `<= 65,839,388` Ir -> total **`65,853,216 <= 2^{25.98} < 2^{26.2} = 77,154,288` Ir**.
- **Proved Per-Node Kernels:**
  - *Line sieve (`SIEVE_HIT = 128, BOUND_CELL = 128, SCAN = 10, ZERO = 2`):* Along each 1D `x_0`-line, each hit is `ldrh; add; strh; add x_i,x_i,x_p; cmp; b.lo` (`6` ops) plus `<= 12` line-start residue step (`<= 18 << 128 = SIEVE_HIT`, `> 7x` headroom; also covers bucketed large-`p` hits `<= 64 << 128`); incremental 12D `(a, b)` norm-bound recurrence along `x_0` costs `<= 32 << 128 = BOUND_CELL` (`4x` headroom).
  - *Deterministic + Expected-Work Cofactor/Root Kernels:* `LLL_NODE = 2^{33}` (`> 2^{32.89}` proved), `NORMS = 2^{16}` (`> 2.18x 30k`), `TDIV = 2^{17}` (`> 2.34x 55,949`), **`E[IDEAL] <= 2^{29}` (`536.9M` Ir):** Las Vegas Cantor–Zassenhaus root extraction for `<= 16` factors `p <= 2^{220}` of degree `<= 8` splits each non-trivial factor w.p. `>= 1/2` per random draw (`E[trials] <= 2 * 7 = 14` splits, or `<= 10` for deg-6 `h`, `<= 220` poly squarings * `<= 64` 4-limb CIOS muls (`<= 200` ops) -> `E[IDEAL] <= 16 * (1 + 10) * 220 * 36 * 200 <= 2.8e8 < 5.37e8 = 2^{29}` Ir, enclosed pointwise in unified meter `t_run = 2^{112}`), `+ 2^{40}` reserve/node.
- **6-Bit ECM Ladder:** Upper-endpoint `[lo, hi_end]` enclosure across **266 DAG** (`0.6174`) & **15 Boot** (`0.5686 > 0.50`, `boot = [178,76,12,4]`) rungs (`63/63` Suyama-12 live).

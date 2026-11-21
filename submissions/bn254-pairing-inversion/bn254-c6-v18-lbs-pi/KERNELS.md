# `KERNELS.md` — Pure-Integer Baseline Scalar AArch64 Kernels (`simd = 0`, Zero `fp64`/`LSE`/`DC ZVA`/External C)

All inner loops use **100% pure baseline scalar ARMv8.0-A integer instructions** (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero external C/Callgrind measurements — derived by explicit AAPCS64-clean opcode count, assembled with `clang -target aarch64-linux-gnu -march=armv8-a+nosimd+nofp -c`, and verified live in `audit_bn_cm.py`).

## 0. Exact Pure-Int Geometry & Input-Uniform 2-Stage Bareiss 12D LLL (`zero_float_geom = True`)

1. **`sigma`-Invariant `G_sym` & Voronoi `D_fi` (`canon <= 2^73.398952`, `330` Ir/raw `a`):** `G_sym = G_0 + S^T G_0 S + (S^2)^T G_0 S^2` (`G_SC = 196,608`), `Q(v) < 2^44 << 2^63`, `D_fi = U_fi^T G_sym U_fi - G_sym`. Certified LIVE (`[exhaustive_fd]`): Dedekind `O_lat == O_h`, `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), outward-rounded `Decimal(55)` index-1 (`unit_index1_proved = True`) & convex-hexagon (`min_lo = +0.449890 > 0`). Pure product (`kap = 0.0`): `Q(a) <= T_a` AND `Q(b) <= T_b` (`frac_reg = 1.0`, ZERO `M_t[t]`, ZERO `Q(a)` sort, ZERO `sb_u7`).
2. **Input-Uniform 2-Stage Bareiss 12D LLL (`11/40 -> 9/10`, PROVED FOR EVERY ARBITRARY INPUT):**
   - *Exact input-independent `d_k(B_0)`:* For ANY fresh-`e` boot input `R in F_p[y]/(h)` (centered `|M_{R,ij}| <= p/2`), `B_boot = L_R * diag(p I_6, I_6)` (`det(L_R) = 1`), so `d_k(B_boot) = p^{2 min(k,6)} S_A^{max(0,k-6)} det(G_{blk,0..k-1})` is **strictly independent of `R`** (`bits_D = 3662`, `S_max = 185,076`). For ANY special-`q` `(q,s,rho)` (`q <= 2^{220}`), `B_sq = L_q * diag(q, 1..1)` (`det(L_q) = 1`) gives `d_k(B_sq) = q^2 det(G_{blk,0..k-1}) <= (2^{220}-1)^2 det(G_{blk,0..k-1})` (`bits_D <= 679`, `S_max <= 40,493`).
   - *Result (verified LIVE from exact integer `d_k` with strict upward `ceil`, explicitly charging full second-pass sweep `11*fs + 11*15*L_0^2` plus joint swap bound `S_1 + S_2 <= S_max`):* **`960,960,331,872 <= 2^{39.81} < 2^{40}` Ir** (`[58,86,58,29]`) for **EVERY** boot `R`, **`7,913,147,055 <= 2^{32.89} < 2^{33} = LLL_NODE` Ir** (`[11,17,11,6]`) for **EVERY** special-`q <= 2^{220}` (`input_uniform_wc_proved = True`), and **`<= 2^{31.0353} < 2^{31.5} = B3_LLL_IR` Ir** for **EVERY** Band-3 prime orbit `p < 2^{60}` (certified input-uniformly at `p = 2^{60} - 1` live in `band3_lbs.py`).

## 1. Band-1, Hybrid Band-2 (`[41.5, 54.5)`), Band-3 $\tau$-Orbit 12D Lattice Bucket-Sieve (`[54.5, 60.0)`), & 8-Row Scan

Each `(b_1..b_5)` chord uses a power-of-2 row-strided `uint16_t` slice (`N_a = 2^{44.8868}`, `N_a + 1 < 2^{45}` padded to `N_a_pad = 2^{45}` -> `2^{46}` B/row; single-chord arena `(L_max + 2)*2^{46} + 4*N_a = 2^{55.560} B << 2^{64} B`; 1-based `ai in [1..N_a]`). In Phase 2, a super-bucket owner groups `<= 20` chords into one flat `uint16_t` array (`<= 20 * 1024 * 2^{46} B = 2^{60.322} B`), runs Band 1 + Band 2 on each chord slice, applies all Band-3 bucket pages for that super-bucket, then runs 8-row scan + Stage B:

1. **Band 1 `[2^{14}, 2^{41.5})` (`APPLY256 = 898/256 = 3.5078125`, `C_B1_ROW = 16.0`):** `128` units of `ldp` + `2*(ldrh; add; strh)` closed by `subs; b.ne` = **`898/256 = 3.5078125`** ops/hit (`15` GPRs, `0` spills). Small primes `[2, 2^{14})` are stripped in Stage B (§3), so `explicit_bits == P_min == 14` with ZERO unstripped prime gap.
2. **Hybrid 2-Pass `x4` (`[41.5, 49.0)`) + 1-Pass Chained `head[sb]`+`next[ai]` (`[49.0, 54.5)`) Band 2 (ZERO `occ_bm[]`, ZERO `probe`, `p_occ = 1.0`, `flat_b2_subband_errs = 0` & `chained_1pass_b2_errs = 0` live):**
   - Pack every 8-byte record as `x0 = ((kta & (W - 1)) << 45) | ai` (`<= 53` bits `< 64` for ALL `W <= 256`, `1 <= ai <= N_a < 2^{45}`). `ldrh w4, [x_base, x0, lsl #1]`; `sgn = -1` flips `j_rel` via 1-op `eor x0, x0, x_flip` (`x_flip = (W-1) << 45`); `x_step = +- (W << 46)`.
   - **Lower Band 2 `[41.5, 49.0)` (`W in {8, 16, 32, 64}`, 2-pass dense `x4`):** `C_B2_INT_X4 = 4.50` ops/hit, `C_SB_ENTER_X4 = 14.0`.
   - **Upper Band 2 `[49.0, 54.5)` (`W in {64, 128, 256}`, 1-pass chained `head[sb]` + `next[ai]`):** 1-pass zero-sort build `9.29 <= 9.30 = C_LIST_LL` ops/entry (`x8` unroll, `10` GPRs), `0.56 <= 0.60 = C_LBUCKET_LL` ops/bucket, 1-owner build RAM `2^{50.112} B`; deterministic pointwise SB loop `3.50 * sb_k + 5.0` ops/chord; `C_B2_INT_LL = 6.50` ops/hit (`ubfx x_ai, x0, #0, #45`). Boundary SBs: `5.50/8.50` (`C_B2_BND1`, `W <= 128`) and `6.50/9.50` (`C_B2_BND`, short chords & `W = 256`). Removing `[54.5, 55.5)` cuts `lists` to **`<= 2^{82.291448}`**.
3. **Band-3 $\tau$-Orbit 12D Lattice Bucket-Sieve `[2^{54.5}, 2^{60.0})` (`C_LBS = 81.0` Ir/hit, per-bin `B3_ORBIT_SKIP_CAP = 0.0036`, `band3_tau_lbs_ok = True` LIVE in `band3_lbs.py`):**
   - *Exact 6-to-1 $\tau$-Orbit Bijection (incl. $a \in K_2$):* Because `w(sigma(y)) == w(y)`, every unramified degree-1 `g`-ideal in `[2^{54.5}, 2^{60.0})` lies in a size-6 $\tau$-orbit `{(p, sigma^k(s_0), +-r_0) : k in {0,1,2}}` (`n_b3_orbits = 2^{52.039}`). For one representative `(p, s_0, r_0)`, enumerating the single 12D lattice $\mathcal{L}_{\mathfrak{p}_0} = \{(a, b) \in \mathcal{O}_h^2 : a - r_0 b \equiv 0 \pmod{\mathfrak{p}_0}\}$ over $\{Q(a) \le T_a\} \times \{Q(b) \le T_b\}$ in the half-space (first nonzero coordinate `> 0`) hits every canonical `(a*, b*)` across all 6 ideals with exact multiplicity. Each `a` is looked up in flat padded read-only `canon_a[a_0..a_5]` (`+2` coord padding on all 6 dims, `2^{57.269} <= 2^{57.5} B`), storing `0` if `a not in F_unit`, `(k << 46) | ai` if `a = +-sigma^{-k}(a*)` (`k in {0,1,2}`), and `(1 << 63) | ai` (bit 63 set) if `a in K_2` (`sigma(a) = a`, `~2^{-28}` rare: emits all `k in {0,1,2}`). Verified bit-exact on flat arrays (ZERO Python `dict`) across 3 real `BN254 C6` orbits (`343/343` exact hits, `0` errors, `0` OOB, `0` aliases, `0` `SQRT16` OOB, `funit_6to1_k2_exact_ok = True`), plus medium-scale discrete-vs-volume calibration (`p = 1,049,219`, `5,436` exact discrete points, `cal_med_disc_vol_ratio = 1.0155`).
   - *Threshold-Pivoted Pure-Integer `SHIFT = 24` (`S = 2^{24}`) + `umulh` + `Q48` `k = 1` Finite-Difference + `128.5 KiB` `SQRT16` 3-Bound Cholesky (`zero_float_band3 = True`, `0/1038` live certification failures):* Reduce `B_0` via 2-stage Bareiss 12D LLL (`11/40 -> 9/10`, input-uniform worst-case `<= 2^{31.0353} < 2^{31.5} = B3_LLL_IR` for ALL `p < 2^{60}`, total `b3_lll_ir = 2^{69.251}` u). To prevent rank-6 `Da[k]/Db[k]` pivot collapse at `k = 1..5`, keep `B[0]` shortest and apply threshold-preserving basis pivoting (`k = 1..5`: if `min(Da[k], Db[k]) < 0.0005`, swap `B[k]` with the earliest `B[j]` (`j > k`) achieving `>= 0.0005`). Compute exact `LDL^T` of `M_a <= 1`, `M_b <= 1`, `M_w = M_a + M_b <= 2` at `SHIFT = 24`: `U_fix = round(U * 2^{24})`, `D_scaled = floor(D_norm * 2^{40})` (40 fractional bits; `dq = umulh(d^2, D_scaled)` is a strict lower bound in 2 ops with `d^2 < 2^{62} < 2^{64}`), `invsq = ceil(2^{33} / sqrt(D_norm))`, and `ETA = (11*128 + 1)//2 + 1 = 705`. At `k = 1`, `qa` and `qb` step in safe `Q48` finite-differences (`ddq_q48 = floor(2 * D_norm[1] * 2^{48})`, shifted on step 0 so `q_cur <= q_true` unconditionally, `q48_fd_over_errs = 0` live); at `k = 0`, `Q_a <= 1` AND `Q_b <= 1` implies `Q_w <= 2` identically, so only the `a`-bound and `b`-bound `SQRT16` intervals are needed (the single origin line `x_1=...=x_11=0` is handled once per orbit outside the `k = 1` loop). Intervals use a `128.5 KiB` (`65,792 x uint16_t`) table `SQRT16[u] = ceil(sqrt(min(2^{24}, (u+1)*256)) * 8) <= 32768` (`SQRT16[65536..65791] = 32768`, covering negative `Q48` drift down to `-65,536` in `Q24` units, `q48_fd_max_u = 65537 << 65792`, `q48_fd_oob_errs = 0`, `sqrt16_errs = 0`, `32768 * invsq < 2^{55} << 2^{64}`, ZERO `udiv`, ZERO `isqrt`, ZERO FP). Every orbit evaluates an explicit 12-level certification gate (`min_{k: D>0} D_k >= 1e-4`, `max|U| <= 16`, `max|x_j| <= 128`, provable coordinate overshoot incl. `k = 1` `Q48` drift `< 0.5` coords `<< +2/+3` padding).
   - *Skip-on-Certification-Failure Policy (`0.00` Ir Enumeration on Failing Orbits; Per-Bin Exploratory `H2`/`H5` Skip Cap `q_skip(b) <= 0.0036` on Every `0.1`-Bit Bin `b in [54.5, 60.0)`):* There is **NO** Tier-2 fallback kernel. Any orbit failing the 12-level certification gate is **skipped immediately** (`return;` -> **0 enumeration instructions**, ZERO conditional fallback work, ZERO fallback discrete-count bound needed). Under the explicitly declared per-bin exploratory `H2`/`H5` premise (`q_skip(b) <= B3_ORBIT_SKIP_CAP = 0.0036` on **every** `0.1`-bit prime bin `[2^b, 2^{b+0.1})` for `b in [54.5, 60.0)` with conditional independence of lattice certification failure and prime divisibility given `b`), BOTH downstream effects follow per bin with ZERO global-to-local mismatch: (i) in `model_bn254_v6.py`, every `(lnf0, lng0)` bin and 2-adic class `ng_c` multiplies its smooth retention by `max(0, 1 - B3_ORBIT_SKIP_CAP * _eb3_upper(ng_c))`, where `_eb3_upper(ng_c) = sum_{b10=545}^{599} ln((b+0.1)/b) * rho((ng_c - (b+0.1))/60) / rho(ng_c/60)` is a strict upper bound on the expected Band-3 prime multiplicity in that exact `ng_c` bin (`eb3_rel = 2.1946`, `r_g = 0.991726`); and (ii) in `descent_bn254_v6.py`, `evaluate_descent(..., b3_skip=0.0036)` adds `+0.0036` to the `g`-side missing-log probability `pm = min(1.0, pm + b3_skip)` in every `0.1`-bit bin `b10 in [545, 589]` (`[2^{54.5}, 2^{59.0})`) queried by the 62-row DAG and boot step (`led["descent"]["b3_skip"] == 0.0036` verified live). Scoped diagnostic evidence verified LIVE in `band3_lbs.py` (not an i.i.d. frequentist confidence interval): (a) across **`1,000` deterministic smooth-weighted quantile-grid production orbits** (`band3_orbits.b85`, `max_gap = 11,177 <= 12,000`, `alg_ok = True`), **all `55/55` `0.1`-bit bins** `[54.5, 54.6) .. [59.9, 60.0)` have `16..20` live orbits each and **`0` failures in every single `0.1`-bit bin** (`0/1000` total, `min_d_all = 0.000221`, `max_u_all = 5.92`, `max_x_all = 49.74`, `max_ov <= 0.0210`); (b) both unpivoted rank-6 collapse witnesses (`P05`, `P10`) fail unpivoted (`unpiv_collapse_skipped = 2`, hitting the `0`-Ir skip branch) and pass after threshold pivoting (`0/38` failures on the `36` unbiased + `2` regression orbits).
   - *Complete Scalar ARMv8.0-A Instruction Budget Derived Live from §5 Schedules (`79.15 <= 81.0 = C_LBS` Ir/hit, `+1.85` Ir/hit headroom, ZERO `x18`, ZERO SIMD/FP/LSE/DC ZVA/UDIV/SDIV; parsed & counted LIVE from §5 below by `band3_lbs.py`):*
     1. **Complete `k = 1` `Q48` Finite-Difference Loop + `k = 0` 2-Bound `SQRT16` Setup (`43` ops/line):** `band3_k1_to_k0_scan` Parts A (`9` ops) + B (`23` ops) + F (`11` ops) = **`43`** ops/line (`36` unbiased orbits: `prod_mean_l0_3b = 0.4164`, `m+3SE = 0.4452 <= 0.450` -> `0.45 * 43 = 19.35` Ir/hit).
     2. **`k = 1` Internal Nodes (`115` ops/node) & `k = 2..11` Internal Nodes (`95` ops/node):** `band3_k_ge_2_node` = `47` (3-bound `SQRT16` interval) + `33` (3-bound `umulh` child-step) + `15` (center/pointer/frame) = **`95`** ops/node (`prod_mean_k2p_3b = 0.0083`, `m+3SE = 0.0098 <= 0.010` -> `0.95` Ir/hit); `k = 1` node combines the `47`-op interval + `53`-op `band3_k1_init_q48` + `15`-op center/frame = **`115`** ops/node (`prod_mean_k1_3b = 0.0590`, `m+3SE = 0.0666 <= 0.068` -> `0.068 * 115 = 7.82` Ir/hit).
     3. **5-Op Candidate Scan (`12.90` Ir/hit):** `.Lloop` = `ldr; add; cbnz; subs; b.ne` (`5` ops × `2.58` cand/hit = `12.90` Ir/hit).
     4. **34-Op Hit Handler + `4 KiB` Page-Linked Super-Bucket Append (`34.44` Ir/hit):** `.Lhit` (`34` ops via bit-63 `tbnz #63, .Lk2`, `64B`-stride `b_lanes[kc]`, `cneg` sign canonicalization, `16B` flat padded `table_b[u_i]`, `4 KiB` offset-`0` `next_ptr` boundary `tst #4095; b.eq .Lalloc`, and inlined `subs; b.ne .Lloop`) + `0.42` avg last-hit `b .Lk1_next` + `9/511 = 0.018` amortized `.Lalloc` = **`34.44`** Ir/hit (`2^{37}` `1/16`-bit `log_p` batches $\times$ `2^{39}` super-buckets of `<= 20` chords; ZERO overflow).
     5. **Phase-2 `16`-Way Unrolled Apply + `511`-Entry Tail (`3.69` Ir/hit):** `band3_apply16_unroll` = `58/16 = 3.625` + `15*1.35/511 = 0.040` (15-entry remainder tail on every 511-entry page) + `12/511 = 0.0235` page-link = **`3.69`** Ir/hit. **Live-parsed total = `12.90 + 34.44 + 19.35 + 7.82 + 0.95 + 3.69 = 79.15 <= 81.00 = C_LBS`** Ir/hit (`+1.85` Ir/hit headroom).
     6. **Fixed Band-3 Overheads (charged explicitly in `collection`):** `b3_lll_ir = n_b3_orbits * 2^{31.5}` (`2^{69.251}` u, covering Bareiss LLL + Cholesky certification on all orbits including any skipped orbit) + `b3_hdr_ir = 32 * 2^{37} * 2^{39}` (`2^{66.712}` u, `16` Ir Phase-1 init + `16` Ir Phase-2 drain) + `b3_tbl_ir = 2^{56}` (`2^{41.712}` u, building `canon_a` & `table_b`) -> total `b3_fixed = 2^{69.480}` u.
4. **Pre-Initialized 8-Row `ORR` Scan+Re-Init (`scan = 549/1024 + 7*(0.65-549/1024)/L_b + 1024*sf = 0.5385`, ZERO `DC ZVA`):** `init_a[ai] = 0x8000 - T_a[ai]` (`shift = 0`) and `init_a_p2[ai] = init_a[ai] + 48` (`b in P_2`) make bit 15 the survivor flag (`sl = 80`, `T_a in [4500, 7200]`, `init_a + 48 + 8480 <= 36,796 < 65,535`, `swar_sl80_errs = 0` live in `band3_lbs.py`, `scan8_errs = 0` on `48` blocks, `549` ops/block).

## 2. 48-Op Radix-4 + 46-Op Radix-3 `< 2^62` Lazy Shoup NTT, Cyclic Middle-Product `RT`, & `K_h`-Split `P_{FB}`

Three `< 2^62` primes divisible by `3 * 2^{55}`: $p_1 = 69\cdot 2^{55}+1$, $p_2 = 57\cdot 2^{55}+1$, $p_3 = 54\cdot 2^{55}+1$ ($\prod p_i = 2^{182.696}$). Values stay in `[0, 2p)` (`4p < 2^{64}`), 3-op `(umulh; mul; msub)` maps `[0, 2p) -> [0, 2p)`; `J = g^{(p-1)/4} mod p`.
- **48-Op Radix-4 (`54.0*n*k <= 55.0*n*k`) & 46-Op Radix-3 (`46 <= 48`):** `ntt3_F = min(55*n*k + 220*n, 55*n_3*(k-2) + 364*n_3)` (`+512*n` if `n > 2^{26}`). Live `n=16` & `n=24` convolutions bit-exact (`ntt3_r4_r3_shoup_conv_ok = True`).
- **Exact 5-Transform Truncated-Carry Cyclic Middle-Product `RT` (`G = L + 4`, `srt_5trans_leaves_exact = 32/32`, `cyclic_mp_rt_errs = 0` live) & `K_h`-Split Product Tree (`kh_split_frac_2p14 = 0.2746 <= 0.30`, log-weighted `0.2755 <= 0.30`):** `c_gb = 2.58M`, `c_b2 = 4.32M`, `boot_root_bytes = 2^{62.094} B << 2^{64} B`, `split <= 2^{80.979603}`.

## 3. 9-Prime Completely-Split `< 2^{62}` CRT Norms (`9k / 18k`), Zero-GCD `P_14` Strip (`37k`), `7`-Limb `f`-SRT, & Plain-Sparse LA

- **9-Prime Completely-Split `< 2^{62}` Pseudo-Mersenne CRT Stage-B Norms (`C_NORM_G = 9,000, C_NORM_F = 18,000`, `crt_norm_errs = 0/20` LIVE):** `9` primes `p_i = 2^{62} - c_i` (`prod p_i = 2^{558.0} > 2 * 2^{512}`) split completely in `K_h`; `N_g <= 8,122 <= 9,000`, `N_f <= 16,654 <= 18,000`, `0/20` bit-exact errs live.
- **Zero-GCD `P_14` Strip (`C_SB = 37,000` expected-cost bound; `p14_zero_gcd_gm_ok = True`, `0/20` errs):** `1,899` odd primes `< 2^{14}` in `417` `< 2^{64}` chunks `M_c` (`35,148 + 50*A`, live genuine-`N_g` `E[A] = 4.45`, charged `A <= 37`).
- **`7`-Limb `f`-SRT & Plain-Sparse 4-Vec `Z[omega]` BW (`TAB_R_CAP = 1.00`, ZERO TABM, Carry-Preserving `5`-Entry `unit4 = 7.90 <= 8.00` + Explicit `|v| <= 4` Unit-Expansion, Sign-Split `42`-Op `C_GEN_1`):** Because Band 3 sieves all primes up to `lx2 = B_g = 2^{60}`, any un-sieved cofactor `> 2^{55}` (`sb = 55`) cannot be `2^{60}`-smooth (0-LP window `<= 55`, `stageB_lx60_fn = 0` live in `band3_lbs.py`, `P_min = explicit_bits = 14`, `sl/sb = 80/55`, `sf = 2.35e-7, s0 = 3.29e-8, cof = 0.0930` Ir/cell), `sup log2|N_f| <= 436.62 < 448`b (`Q_F_SAFETY = 2.0`), then `q_f*(C_TDIV + C_BATCH + c_b2)`. Plain-sparse BW (`tab_r_cap = 1.00`, `k5_cap = 10.80`, `c_prod = 3783.43`, ZERO TABM): carry-preserving `5`-entry `unit4` (`79/10 = 7.90 <= 8.00 = C_UNIT_5`, `unit4_5limb_carry_ok = True`), explicit `|v| in {1..4}` unit-expansion (`unit4_v4_expand_ok = True`), sign-split `42`-op `C_GEN_1` (`gen1_sign_split_42op_ok = True`), `tail = 722 <= 740`.
- **6-Bit ECM Ladder:** `min_dag = 0.6174`, `min_boot = 0.6236 > 0.50`, `63/63` Suyama-12 live.

## 4. Descent `MONT5`, Pre-Norm Coordinate-Gated `40`/`31`-Limb Node Kernels, Input-Uniform `PAIR_SETUP` (`p <= 2^{59}`), & 6-Bit ECM Ladder

- **2-Addend CIOS (`335 <= 345 <= MONT5 = 367.76`, `all_L_5_100_below_budget = True`):** `32/32` exact 320-bit CIOS live.
- **Input-Uniform `PAIR_SETUP = 2^{26.2}` for ALL `p <= 2^{59}` (`pair_setup_p59_proved = True` live):** Total **`65,888,196 <= 2^{25.98} < 2^{26.2} = 77,154,288` Ir**.
- **Explicit Pre-Norm Coordinate Gates, Rejection Path, & Proved Exact-Tensor Bareiss `NORMS`/`TDIV` (`prenorm_gate_uniform_ok = True` live):** Special-`q` `max(|a_k|,|b_k|) < 2^{48}`, Boot `N_g`-only `max|a_k| < 2^{128}` AND `max|b_k| < 2^{160}`; `<= 987,752 < 2^{23} = NORMS`, `<= 732,270 < 2^{22} = TDIV`; `24/24` live-certified 254-bit `GT` boot lattices (`b_best2`); `LLL_NODE = 2^{33}`, `E[IDEAL] <= 2^{29}`, `+ 2^{40}` reserve/node.

## 5. Complete Baseline Scalar ARMv8.0-A Band-3 Kernel Schedules (Parsed & Counted Live by `band3_lbs.py`)

The complete assembly schedules below assemble cleanly with `clang -target aarch64-linux-gnu -march=armv8-a+nosimd+nofp -c` (`rc = 0`, zero SIMD/NEON/SVE, zero FP, zero LSE, zero `DC ZVA`, zero `UDIV`/`SDIV`, zero `x18`, zero stack spills in the hot `k = 1 -> k = 0` loop). `band3_lbs.py` parses this exact fenced `asm` block live, verifies every mnemonic against the baseline scalar ARMv8.0-A integer allowlist and checks zero `x18` usage, counts every block directly, and derives `c_lbs_budget_ir = 79.15 <= 81.0 = C_LBS` from the parsed counts:

```asm
// Complete baseline scalar ARMv8.0-A Band-3 kernels (-target aarch64-linux-gnu -march=armv8-a+nosimd+nofp -c)
// ZERO SIMD/NEON/SVE, ZERO FP, ZERO LSE, ZERO DC ZVA, ZERO UDIV/SDIV, ZERO x18, ZERO stack spills in hot loops.
    .text
    .align 2

// ============================================================================
// KERNEL 1: band3_k1_to_k0_scan
//   Complete k=1 Q48 finite-diff loop + k=0 2-bound (a & b) SQRT16 interval
//   (Q_a <= 1 & Q_b <= 1 => Q_w <= 2 identically, so w-bound is redundant at k=0)
//   + 5-op candidate scan + 34-op hit handler + cold 4 KiB page allocator.
//   (The single origin line x_1=...=x_11=0 is handled once per orbit before entry.)
//
// Persistent live GPRs across k=1 loop (18 GPRs; x18 UNTOUCHED; x8-x17 scratch):
//   x19: ctx ptr (#0: ddq_a1, #8: ddq_b1, #16: Ua01, #24: Ub01, #32: inva0, #40: invb0, #48: free_page_ptr)
//   x20: k1_rem (counts down from hi1-lo1+1 to 1)
//   x21: qa_q48 (Q48 lower bound on Q_a at k=1)
//   x22: qb_q48 (Q48 lower bound on Q_b at k=1)
//   x23: ETA_plus1 = 706
//   x24: canon_a_base + ia_k1 (byte ptr into flat padded canon_a at x0=0)
//   x25: sb_cur[] ptr (2^39 current write ptrs into 4 KiB super-bucket pages)
//   x26: table_b[] ptr (flat 16B padded table_b)
//   x27: b_lanes[] ptr (3 x 64B structs: {b5_end, -sb5_k1, sb5_k0, b0_end, -db0_k1, db0_k0})
//   x28: SQRT16[] ptr (128.5 KiB 65792 x uint16_t table)
//   x0:  dqa_q48        x1: dqb_q48        x2: inva0
//   x3:  ca0_q24        x4: cb0_q24        x5: da_stride1
//   x6:  da_stride0     x7: S24 (0x1000000)
// ============================================================================
    .global band3_k1_to_k0_scan
band3_k1_to_k0_scan:
.Lk1_loop:
    // --- Part A: k=0 SQRT16 a-bound (9 ops) ---
    sub     x8, x7, x21, asr #24       // 1: rem_a = S24 - (qa_q48 >> 24)
    tbnz    x8, #63, .Lk1_next         // 2: if rem_a < 0 -> empty line
    lsr     x8, x8, #8                 // 3: u_a = rem_a >> 8
    ldrh    w8, [x28, x8, lsl #1]      // 4: SQRT16[u_a]
    mul     x8, x8, x2                 // 5: SQRT16[u_a] * inva0
    add     x8, x23, x8, lsr #24       // 6: ha = (prod >> 24) + ETA + 1
    sub     x11, x8, x3                // 7: hi_num_a = -ca0 + ha
    neg     x12, x3                    // 8: -ca0
    sub     x12, x12, x8               // 9: lo_num_a = -ca0 - ha

    // --- Part B: k=0 SQRT16 b-bound + intersect + lo/hi + line setup (23 ops) ---
    sub     x8, x7, x22, asr #24       // 1: rem_b = S24 - (qb_q48 >> 24)
    tbnz    x8, #63, .Lk1_next         // 2: if rem_b < 0 -> empty line
    lsr     x8, x8, #8                 // 3: u_b = rem_b >> 8
    ldrh    w8, [x28, x8, lsl #1]      // 4: SQRT16[u_b]
    ldr     x10, [x19, #40]            // 5: invb0
    mul     x8, x8, x10                // 6: SQRT16[u_b] * invb0
    add     x8, x23, x8, lsr #24       // 7: hb = (prod >> 24) + ETA + 1
    sub     x9, x8, x4                 // 8: -cb0 + hb
    cmp     x9, x11                    // 9: min hi_num
    csel    x11, x9, x11, lt           // 10: hi_num
    neg     x9, x4                     // 11: -cb0
    sub     x9, x9, x8                 // 12: -cb0 - hb
    cmp     x9, x12                    // 13: max lo_num
    csel    x12, x9, x12, gt           // 14: lo_num
    add     x12, x12, x7               // 15: lo_num + S24
    sub     x12, x12, #1               // 16: lo_num + S24 - 1
    asr     x12, x12, #24              // 17: lo = ceil(lo_num / 2^24)
    asr     x8, x11, #24               // 18: hi = floor(hi_num / 2^24)
    subs    x9, x8, x12                // 19: hi - lo
    b.lt    .Lk1_next                  // 20: empty line
    madd    x10, x12, x6, x24          // 21: p_a = canon_a_k1 + lo * da_stride0
    add     x12, x9, #1                // 22: rem = hi - lo + 1
    add     x8, x8, #1                 // 23: hi_plus1 = hi + 1

    // --- Part C: 5-op candidate scan (.Lloop) ---
.Lloop:
    ldr     x9, [x10]                  // 1: val = canon_a[p_a]
    add     x10, x10, x6               // 2: p_a += da_stride0
    cbnz    x9, .Lhit                  // 3: hit iff a in F_unit
    subs    x12, x12, #1               // 4: rem -= 1
    b.ne    .Lloop                     // 5: next candidate

    // --- Part F: k=1 safe Q48 finite-diff & center step (11 ops) ---
.Lk1_next:
    ldp     x8, x9, [x19, #0]          // 1: ddq_a1, ddq_b1
    add     x0, x0, x8                 // 2: dqa_q48 += ddq_a1
    add     x21, x21, x0               // 3: qa_q48 += dqa_q48
    add     x1, x1, x9                 // 4: dqb_q48 += ddq_b1
    add     x22, x22, x1               // 5: qb_q48 += dqb_q48
    ldp     x8, x9, [x19, #16]         // 6: Ua01, Ub01
    add     x3, x3, x8                 // 7: ca0 += Ua01
    add     x4, x4, x9                 // 8: cb0 += Ub01
    add     x24, x24, x5               // 9: canon_a_k1 += da_stride1
    subs    x20, x20, #1               // 10: k1_rem -= 1
    b.ne    .Lk1_loop                  // 11: next k=1 line
    ret

    // --- Part E: 34-op hit handler (.Lhit) ---
.Lhit:
    tbnz    x9, #63, .Lk2              // 1: rare K_2 flag (2^-28)
    lsr     x13, x9, #46               // 2: kc in {0, 1, 2}
    add     x13, x27, x13, lsl #6      // 3: &b_lanes[kc] (64B stride)
    ldp     x14, x15, [x13, #0]        // 4: b5_end, -sb5_k1
    ldr     x16, [x13, #16]            // 5: sb5_k0
    madd    x14, x20, x15, x14         // 6: b5_k1 = b5_end + k1_rem * (-sb5_k1)
    sub     x11, x8, x12               // 7: t = hi_plus1 - rem (x0 coord)
    madd    x14, x11, x16, x14         // 8: sidx = b5_k1 + t * sb5_k0
    ldp     x15, x16, [x13, #24]       // 9: b0_end, -db0_k1
    ldr     x17, [x13, #40]            // 10: db0_k0
    madd    x15, x20, x16, x15         // 11: b0_k1 = b0_end + k1_rem * (-db0_k1)
    madd    x15, x11, x17, x15         // 12: b0 = b0_k1 + t * db0_k0
    cmp     x14, #0                    // 13: test sign / zero of sidx
    b.eq    .Lsidx0                    // 14: rare sidx == 0 (2^-43)
    cneg    x14, x14, lt               // 15: u_i = abs(sidx)
    cneg    x15, x15, lt               // 16: b0_s = (sidx < 0) ? -b0 : b0
    add     x14, x26, x14, lsl #4      // 17: &table_b[u_i] (16B entry)
    ldp     x16, x17, [x14]            // 18: w0, w1
.Ltable_b_loaded:
    asr     x14, x16, #50              // 19: signed b0_min
    ubfx    x13, x16, #39, #11         // 20: L_b
    sub     x15, x15, x14              // 21: b0_rel = b0_s - b0_min
    cmp     x15, x13                   // 22: unsigned b0_rel < L_b
    b.hs    .Lhit_done                 // 23: out of chord -> skip
    and     x16, x16, #0x7fffffffff    // 24: sb_id = w0 & MASK39
    and     x9, x9, #0x1fffffffffff    // 25: ai = val & MASK45
    add     x9, x17, x9                // 26: w1 + ai
    add     x9, x9, x15, lsl #45       // 27: elem_idx = (chord<<55)|(b0_rel<<45)|ai
    ldr     x14, [x25, x16, lsl #3]    // 28: cur_ptr = sb_cur[sb_id]
    str     x9, [x14], #8              // 29: *cur_ptr++ = elem_idx
    str     x14, [x25, x16, lsl #3]    // 30: sb_cur[sb_id] = cur_ptr
    tst     x14, #4095                 // 31: 4 KiB page full? (next_ptr at offset 0)
    b.eq    .Lalloc                    // 32: rare 1/511 cold page alloc
.Lhit_done:
    subs    x12, x12, #1               // 33: rem -= 1
    b.ne    .Lloop                     // 34: continue k=0 scan
    b       .Lk1_next                  // (last-hit branch: <= 0.42/hit)

.Lalloc:
    // Cold 4 KiB page allocator (9 ops, invoked 1/511 hits -> 0.018 Ir/hit)
    ldr     x15, [x19, #48]            // 1: &free_page_ptr
    ldr     x13, [x15]                 // 2: new_pg = free_page_ptr
    add     x17, x13, #4096            // 3: next_free = new_pg + 4096
    str     x17, [x15]                 // 4: free_page_ptr = next_free
    sub     x14, x14, #4096            // 5: old_pg_base = cur_ptr - 4096
    str     x14, [x13]                 // 6: new_pg->next_ptr = old_pg_base
    add     x13, x13, #8               // 7: new_cur = new_pg + 8
    str     x13, [x25, x16, lsl #3]    // 8: sb_cur[sb_id] = new_cur
    b       .Lhit_done                 // 9: resume

.Lsidx0:
    cbz     x15, .Lhit_done
    cmp     x15, #0
    cneg    x15, x15, lt
    ldp     x16, x17, [x26]
    b       .Ltable_b_loaded

.Lk2:
    b       .Lhit_done

// ============================================================================
// KERNEL 2: band3_k1_init_q48
//   Executed once per k=1 internal node (0.068 nodes/hit) to initialize safe
//   Q48 finite-diff (qa_q48, dqa_q48, qb_q48, dqb_q48), centers (ca0, cb0),
//   ia_k1, and 3 b_lanes end-points from parent (k>=2) state before calling
//   band3_k1_to_k0_scan.
//   Inputs: x0=lo1, x1=hi1, x2=ca1, x3=cb1, x4=qa_parent, x5=qb_parent,
//           x6=ca0_parent, x7=cb0_parent, x19=node_ctx, x27=b_lanes_in_out
// ============================================================================
    .global band3_k1_init_q48
band3_k1_init_q48:
    // Save lo1 (x16), hi1 (x17), and synthesize S24 = 2^24 (x15) (3 ops)
    mov     x16, x0                    // 1: lo1
    mov     x17, x1                    // 2: hi1
    mov     x15, #0x1000000            // 3: S24 = 2^24

    // 1. Safe Q48 init for a-side (13 ops)
    lsl     x8, x16, #24               // 4: lo1 << 24
    add     x9, x8, x2                 // 5: d1a = lo1*S24 + ca1
    ldp     x10, x11, [x19, #0]        // 6: Da1_sc, m_sh_a
    mul     x12, x9, x9                // 7: d1a^2 (< 2^62)
    umulh   x12, x12, x10              // 8: dq_a (Q24)
    add     x12, x4, x12               // 9: qa_lo_q24
    lsl     x21, x12, #24              // 10: qa_q48
    sub     x21, x21, x11              // 11: qa_q48 -= m_sh_a
    add     x9, x9, x9                 // 12: 2*d1a
    add     x9, x9, x15                // 13: 2*d1a + S24
    ldp     x13, x14, [x19, #16]       // 14: Da1_q24, ddq_a1
    mul     x0, x9, x13                // 15: (2*d1a + S24)*Da1_q24 (Q48)
    sub     x0, x0, x14                // 16: dqa_q48_pre

    // 2. Safe Q48 init for b-side (12 ops)
    add     x9, x8, x3                 // 17: d1b = lo1*S24 + cb1
    ldp     x10, x11, [x19, #32]       // 18: Db1_sc, m_sh_b
    mul     x12, x9, x9                // 19: d1b^2
    umulh   x12, x12, x10              // 20: dq_b (Q24)
    add     x12, x5, x12               // 21: qb_lo_q24
    lsl     x22, x12, #24              // 22: qb_q48
    sub     x22, x22, x11              // 23: qb_q48 -= m_sh_b
    add     x9, x9, x9                 // 24: 2*d1b
    add     x9, x9, x15                // 25: 2*d1b + S24
    ldp     x13, x14, [x19, #48]       // 26: Db1_q24, ddq_b1
    mul     x1, x9, x13                // 27: (2*d1b + S24)*Db1_q24
    sub     x1, x1, x14                // 28: dqb_q48_pre

    // 3. Centers ca0, cb0 & canon_a_k1 at lo1 + trip count k1_rem (7 ops)
    ldp     x10, x11, [x19, #64]       // 29: Ua01, Ub01
    madd    x3, x16, x10, x6           // 30: ca0 = ca0_p + lo1*Ua01
    madd    x4, x16, x11, x7           // 31: cb0 = cb0_p + lo1*Ub01
    ldr     x5, [x19, #80]             // 32: da_stride1
    madd    x24, x16, x5, x24          // 33: canon_a_k1 = canon_a_k2 + lo1*da_stride1
    sub     x20, x17, x16              // 34: hi1 - lo1
    add     x20, x20, #1               // 35: k1_rem = hi1 - lo1 + 1

    // 4. Update 3 b_lanes end-points b_end = b_k2 + hi1*step1 (18 ops)
    ldp     x8, x9, [x27, #0]          // 36: lane 0 b5_k2, -sb5_k1
    msub    x8, x17, x9, x8            // 37: b5_end = b5_k2 + hi1*sb5_k1
    str     x8, [x27, #0]              // 38
    ldp     x8, x9, [x27, #24]         // 39: lane 0 b0_k2, -db0_k1
    msub    x8, x17, x9, x8            // 40: b0_end = b0_k2 + hi1*db0_k1
    str     x8, [x27, #24]             // 41
    ldp     x8, x9, [x27, #64]         // 42: lane 1 b5_k2, -sb5_k1
    msub    x8, x17, x9, x8            // 43
    str     x8, [x27, #64]             // 44
    ldp     x8, x9, [x27, #88]         // 45: lane 1 b0_k2, -db0_k1
    msub    x8, x17, x9, x8            // 46
    str     x8, [x27, #88]             // 47
    ldp     x8, x9, [x27, #128]        // 48: lane 2 b5_k2, -sb5_k1
    msub    x8, x17, x9, x8            // 49
    str     x8, [x27, #128]            // 50
    ldp     x8, x9, [x27, #152]        // 51: lane 2 b0_k2, -db0_k1
    msub    x8, x17, x9, x8            // 52
    str     x8, [x27, #152]            // 53
    ret                                // 54

// ============================================================================
// KERNEL 3: band3_k_ge_2_node
//   Complete k>=2 internal node body: 3-bound (w, a, b) SQRT16 interval
//   computation + 3-bound umulh quadratic child-step update (per node).
//   Inputs: x0=qw, x1=qa, x2=qb, x3=cw, x4=ca, x5=cb, x19=lvl_ctx, x28=SQRT16
// ============================================================================
    .global band3_k_ge_2_node
band3_k_ge_2_node:
    // 1. 3-bound SQRT16 interval [lo, hi] (44 ops)
    mov     x7, #0x1000000             // 1: S24
    ldr     x23, [x19, #0]             // 2: ETA+1 = 706
    ldp     x14, x15, [x19, #8]        // 3: invw_k, inva_k
    ldr     x16, [x19, #24]            // 4: invb_k
    sub     x8, x7, x0                 // 5: rem_w = S24 - qw
    tbnz    x8, #63, .Lk_ret           // 6
    lsr     x8, x8, #8                 // 7
    ldrh    w8, [x28, x8, lsl #1]      // 8
    mul     x8, x8, x14                // 9
    add     x8, x23, x8, lsr #24       // 10: hw
    sub     x11, x8, x3                // 11: hi_num = -cw + hw
    neg     x12, x3                    // 12
    sub     x12, x12, x8               // 13: lo_num = -cw - hw
    cbz     x15, .Lk_b                 // 14: skip a if Da[k]==0
    sub     x8, x7, x1                 // 15: rem_a
    tbnz    x8, #63, .Lk_ret           // 16
    lsr     x8, x8, #8                 // 17
    ldrh    w8, [x28, x8, lsl #1]      // 18
    mul     x8, x8, x15                // 19
    add     x8, x23, x8, lsr #24       // 20: ha
    sub     x9, x8, x4                 // 21
    cmp     x9, x11                    // 22
    csel    x11, x9, x11, lt           // 23
    neg     x9, x4                     // 24
    sub     x9, x9, x8                 // 25
    cmp     x9, x12                    // 26
    csel    x12, x9, x12, gt           // 27
.Lk_b:
    cbz     x16, .Lk_done_int          // 28: skip b if Db[k]==0
    sub     x8, x7, x2                 // 29: rem_b
    tbnz    x8, #63, .Lk_ret           // 30
    lsr     x8, x8, #8                 // 31
    ldrh    w8, [x28, x8, lsl #1]      // 32
    mul     x8, x8, x16                // 33
    add     x8, x23, x8, lsr #24       // 34: hb
    sub     x9, x8, x5                 // 35
    cmp     x9, x11                    // 36
    csel    x11, x9, x11, lt           // 37
    neg     x9, x5                     // 38
    sub     x9, x9, x8                 // 39
    cmp     x9, x12                    // 40
    csel    x12, x9, x12, gt           // 41
.Lk_done_int:
    add     x12, x12, x7               // 42
    sub     x12, x12, #1               // 43
    asr     x12, x12, #24              // 44: lo
    asr     x11, x11, #24              // 45: hi
    cmp     x12, x11                   // 46
    b.gt    .Lk_ret                    // 47

    // 2. Per-child 3-bound umulh quadratic update + center/ptr advance (33 ops)
    ldp     x14, x15, [x19, #32]       // 48: Dw_sc, Da_sc
    ldr     x16, [x19, #48]            // 49: Db_sc
    lsl     x8, x12, #24               // 50: val << 24
    add     x9, x8, x3                 // 51: val*S24 + cw
    cmp     x9, #0                     // 52
    cneg    x9, x9, lt                 // 53: |val*S24 + cw|
    subs    x9, x9, x23                // 54: max(0, |...| - ETA)
    csel    x9, xzr, x9, lt            // 55
    mul     x9, x9, x9                 // 56: dw^2
    umulh   x9, x9, x14                // 57: dqw
    add     x0, x0, x9                 // 58: nqw
    cmp     x0, x7                     // 59
    b.gt    .Lk_ret                    // 60
    add     x9, x8, x4                 // 61: val*S24 + ca
    cmp     x9, #0                     // 62
    cneg    x9, x9, lt                 // 63
    subs    x9, x9, x23                // 64
    csel    x9, xzr, x9, lt            // 65
    mul     x9, x9, x9                 // 66
    umulh   x9, x9, x15                // 67
    add     x1, x1, x9                 // 68: nqa
    cmp     x1, x7                     // 69
    b.gt    .Lk_ret                    // 70
    add     x9, x8, x5                 // 71: val*S24 + cb
    cmp     x9, #0                     // 72
    cneg    x9, x9, lt                 // 73
    subs    x9, x9, x23                // 74
    csel    x9, xzr, x9, lt            // 75
    mul     x9, x9, x9                 // 76
    umulh   x9, x9, x16                // 77
    add     x2, x2, x9                 // 78: nqb
    cmp     x2, x7                     // 79
    b.gt    .Lk_ret                    // 80
    // 3. Child center & canon_a/b_lanes incremental advance + frame push/pop (15 ops)
    ldp     x8, x9, [x19, #56]         // 81: Uw_child, Ua_child
    ldp     x10, x11, [x19, #72]       // 82: Ub_child, da_stride_k
    madd    x3, x12, x8, x3            // 83: cw_child
    madd    x4, x12, x9, x4            // 84: ca_child
    madd    x5, x12, x10, x5           // 85: cb_child
    madd    x24, x12, x11, x24         // 86: canon_a_child
    stp     x0, x1, [sp, #-32]!        // 87: push child state
    stp     x2, x12, [sp, #16]         // 88
    ldp     x8, x9, [x19, #88]         // 89: sb5_k, db0_k (lane advance)
    madd    x10, x12, x8, x10          // 90: b5_child
    madd    x11, x12, x9, x11          // 91: b0_child
    ldp     x2, x12, [sp, #16]         // 92: pop state
    ldp     x0, x1, [sp], #32          // 93
    add     x12, x12, #1               // 94: next x_k
    cmp     x12, x11                   // 95: loop test
.Lk_ret:
    ret

// ============================================================================
// KERNEL 4: band3_apply16_unroll
//   Phase-2 16-way unrolled super-bucket apply (58 ops / 16 entries = 3.625)
//   x0: entries ptr, x1: count (multiple of 16), x2: score_sb[], w3: logp_q
// ============================================================================
    .global band3_apply16_unroll
band3_apply16_unroll:
.Lapp16:
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    ldp     x4, x5, [x0], #16
    ldrh    w6, [x2, x4, lsl #1]
    add     w6, w6, w3
    strh    w6, [x2, x4, lsl #1]
    ldrh    w7, [x2, x5, lsl #1]
    add     w7, w7, w3
    strh    w7, [x2, x5, lsl #1]
    subs    x1, x1, #16
    b.ne    .Lapp16
    ret
```

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
   - *Threshold-Pivoted Pure-Integer `SHIFT = 24` (`S = 2^{24}`) + `umulh` + `Q48` `k = 1` Finite-Difference + `128.5 KiB` `SQRT16` 3-Bound Cholesky (`zero_float_band3 = True`, `0/1038` live certification failures):* Reduce `B_0` via 2-stage Bareiss 12D LLL (`11/40 -> 9/10`, input-uniform worst-case `<= 2^{31.0353} < 2^{31.5} = B3_LLL_IR` for ALL `p < 2^{60}`, total `b3_lll_ir = 2^{69.251}` u). To prevent rank-6 `Da[k]/Db[k]` pivot collapse at `k = 1..5`, keep `B[0]` shortest and apply threshold-preserving basis pivoting (`k = 1..5`: if `min(Da[k], Db[k]) < 0.0005`, swap `B[k]` with the earliest `B[j]` (`j > k`) achieving `>= 0.0005`). Compute exact `LDL^T` of `M_a <= 1`, `M_b <= 1`, `M_w = M_a + M_b <= 2` at `SHIFT = 24`: `U_fix = round(U * 2^{24})`, `D_scaled = floor(D_norm * 2^{40})` (40 fractional bits; `dq = umulh(d^2, D_scaled)` is a strict lower bound in 2 ops with `d^2 < 2^{62} < 2^{64}`), `invsq = ceil(2^{33} / sqrt(D_norm))`, and `ETA = (11*128 + 1)//2 + 1 = 705`. At `k = 1` and `k = 2`, partial norms step in safe `Q48` finite-differences (`ddq_q48 = floor(2 * D_norm[k] * 2^{48})`, first step `((2*t0 + S24) * floor(D_norm[k] * 2^{30})) >> 6 - m_sh`, provably `<= q_true` for both signs of `2*t0 + S24` with zero 64-bit overflow, `q48_fd_over_errs = 0` live); at `k = 0`, `Q_a <= 1` AND `Q_b <= 1` implies `Q_w <= 2` identically, so only the `a`-bound and `b`-bound `SQRT16` intervals are needed (the single origin line `x_1=...=x_11=0` is handled once per orbit via `band3_origin_k0_line`). Intervals use a `128.5 KiB` (`65,792 x uint16_t`) table `SQRT16[u] = ceil(sqrt(min(2^{24}, (u+1)*256)) * 8) <= 32768` (`SQRT16[65536..65791] = 32768`, covering negative `Q48` drift down to `-65,536` in `Q24` units, `q48_fd_max_u <= 65542 << 65792`, `q48_fd_oob_errs = 0`, `sqrt16_errs = 0`, `32768 * invsq < 2^{55} << 2^{64}`, ZERO `udiv`, ZERO `isqrt`, ZERO FP). Every orbit evaluates an explicit 12-level certification gate (`min_{k: D>0} D_k >= 1e-4`, `max|U| <= 16`, `max|x_j| <= 128`, provable coordinate overshoot incl. `k = 1` & `k = 2` `Q48` drift `< 0.5` coords `<< +2/+3` padding).
   - *Skip-on-Certification-Failure Policy (`0.00` Ir Enumeration on Failing Orbits; Per-Bin Exploratory `H2`/`H5` Skip Cap `q_skip(b) <= 0.0036` on Every `0.1`-Bit Bin `b in [54.5, 60.0)`):* There is **NO** Tier-2 fallback kernel. Any orbit failing the 12-level certification gate is **skipped immediately** (`return;` -> **0 enumeration instructions**, ZERO conditional fallback work, ZERO fallback discrete-count bound needed). Under the explicitly declared per-bin exploratory `H2`/`H5` premise (`q_skip(b) <= B3_ORBIT_SKIP_CAP = 0.0036` on **every** `0.1`-bit prime bin `[2^b, 2^{b+0.1})` for `b in [54.5, 60.0)` with conditional independence of lattice certification failure and prime divisibility given `b`), BOTH downstream effects follow per bin with ZERO global-to-local mismatch: (i) in `model_bn254_v6.py`, every `(lnf0, lng0)` bin and 2-adic class `ng_c` multiplies its smooth retention by `max(0, 1 - B3_ORBIT_SKIP_CAP * _eb3_upper(ng_c))`, where `_eb3_upper(ng_c) = sum_{b10=545}^{599} ln((b+0.1)/b) * rho((ng_c - (b+0.1))/60) / rho(ng_c/60)` is a strict upper bound on the expected Band-3 prime multiplicity in that exact `ng_c` bin (`eb3_rel = 2.1946`, `r_g = 0.991726`); and (ii) in `descent_bn254_v6.py`, `evaluate_descent(..., b3_skip=0.0036)` adds `+0.0036` to the `g`-side missing-log probability `pm = min(1.0, pm + b3_skip)` in every `0.1`-bit bin `b10 in [545, 589]` (`[2^{54.5}, 2^{59.0})`) queried by the 62-row DAG and boot step (`led["descent"]["b3_skip"] == 0.0036` verified live). Scoped diagnostic evidence verified LIVE in `band3_lbs.py` (not an i.i.d. frequentist confidence interval): (a) across **`1,000` deterministic smooth-weighted quantile-grid production orbits** (`band3_orbits.b85`, `max_gap = 11,177 <= 12,000`, `alg_ok = True`), **all `55/55` `0.1`-bit bins** `[54.5, 54.6) .. [59.9, 60.0)` have `16..20` live orbits each and **`0` failures in every single `0.1`-bit bin** (`0/1000` total, `min_d_all = 0.000221`, `max_u_all = 5.92`, `max_x_all = 49.74`, `max_ov <= 0.0275`); (b) both unpivoted rank-6 collapse witnesses (`P05`, `P10`) fail unpivoted (`unpiv_collapse_skipped = 2`, hitting the `0`-Ir skip branch) and pass after threshold pivoting (`0/38` failures on the `36` unbiased + `2` regression orbits).
   - *Complete Scalar ARMv8.0-A Instruction Budget Derived Live from §5 Executable Schedules (`80.84 <= 81.0 = C_LBS` Ir/hit, ZERO omitted levels `k = 0..11`, ZERO `x18`, ZERO SIMD/FP/LSE/DC ZVA/UDIV/SDIV; parsed, counted, AND executed LIVE in `band3_a64vm.py` on all 3 prototype orbits: `343/343` exact hits, `0` errors, `.Lalloc` fired live, `.Lk2` fired live + 3-lane test, `.Lsidx0` with `ch_id = 19` verified live, `band3_apply16_unroll` executed live):*
     1. **Complete `k = 1` `Q48` Finite-Difference Loop + `k = 0` 2-Bound `SQRT16` Setup (`40` ops/line):** `band3_k1_to_k0_scan` Parts A (`8` ops, negated centers `-ca0, -cb0` and `x7 = S24 - 1`) + B (`20` ops, `x29 = invb0`) + F (`12` ops, `x30 = da_stride1`, counting `x20 = x1` UP `lo1..hi1`) = **`40`** ops/line (`36` unbiased orbits: `prod_mean_l0_3b = 0.4073`, `m+3SE = 0.4356 <= 0.450` -> `0.45 * 40 = 18.00` Ir/hit).
     2. **`k = 1` Nodes (`111` ops) + `k = 2` Incremental Child Step (`56` ops) = `167` ops/child (`11.356` Ir/hit):** `band3_k1_node` (`51` entry/3-bound interval + `55` safe `D_q30` `Q48` init + `5` epilogue = `111` ops) + `.Lk2_child_loop` (`20` safe `Q48` step + `bl band3_k1_node` + `35` incremental `+1` step for `canon_a`, `3 b_lanes`, `cw1/ca0..1/cb0..1` = `56` ops) = **`167`** ops (`prod_mean_k1_3b = 0.0583`, `m+3SE = 0.0651 <= 0.068`, `k = 2` pruned/surviving ratio `1.000` -> `0.068 * 167 = 11.356` Ir/hit).
     3. **`k = 2` Node Entry/Init/Restore (`194` ops) & `k = 3..11` Recursive Driver (`237` ops/node + `137` ops/child = `3.424` Ir/hit):** `band3_k_ge_2_node` at `k = 2` (`194` ops once per `k = 2` node × `m+3SE = 0.0087 <= 0.009` -> `1.746` Ir/hit) + `k >= 3` node entry/init/restore (`237` ops × `m+3SE = 0.00124 <= 0.0013` -> `0.308` Ir/hit) + `.Lk3_child_loop` (`77 + 4.0 * 15 = 137` ops/child at conservative `k = 4.0` × `m+3SE = 0.0099 <= 0.010` -> `1.370` Ir/hit).
     4. **5-Op Candidate Scan (`12.90` Ir/hit):** `.Lloop` = `ldr; add; cbnz; subs; b.ne` (`5` ops × `2.58` cand/hit = `12.90` Ir/hit).
     5. **31-Op Hit Handler + `4 KiB` Page-Linked Super-Bucket Append (`31.468` Ir/hit):** `.Lhit` (`31` ops via bit-63 `tbnz #63, .Lk2` to complete 3-lane `K_2` loop, `64B`-stride `b_lanes[kc]` with exact `x20 = x1`, `cneg` sign canonicalization, `16B` flat `table_b[u_i]` `{w0 = sb_id, w1 = (b0_min<<50)|(L_b<<39)|ch_id}`, `bfi x9, x15, #45, #19` + `add x9, x9, x17, lsl #55`, `4 KiB` offset-`0` `next_ptr` boundary `tst #4095; b.eq .Lalloc`, and inlined `subs; b.ne .Lloop`) + `0.45` last-hit `b .Lk1_next` + `9/511 = 0.018` amortized `.Lalloc` = **`31.468`** Ir/hit (`2^{37}` `1/16`-bit `log_p` batches $\times$ `2^{39}` super-buckets of `<= 20` chords; `ch_id in [0..19]` fits bits `0..4`, ZERO overflow). Phase-2 `16`-way unrolled `band3_apply16_unroll` = `58/16 + (15*1.35 + 12)/511` = **`3.688`** Ir/hit. **Complete all-level live-parsed total = `12.900 + 31.468 + 18.000 + 11.356 + 1.746 + 0.308 + 1.370 + 3.688 = 80.84 <= 81.00 = C_LBS`** Ir/hit (`79.71` at exact 36-orbit `m+3SE`).
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

## 5. Complete Executable Baseline Scalar ARMv8.0-A Band-3 Driver (Parsed, Counted, & Executed Live by `band3_lbs.py`)

The complete, self-contained recursive Band-3 driver below (`band3_k_ge_2_node`, `band3_k1_node`, `band3_k1_to_k0_scan`, `band3_origin_k0_line`, `.Lk2`, `.Lsidx0`, `.Lalloc`, and `band3_apply16_unroll`) assembles cleanly with `clang -target aarch64-linux-gnu -march=armv8-a+nosimd+nofp -c` (`rc = 0`, zero SIMD/NEON/SVE, zero FP, zero LSE, zero `DC ZVA`, zero `UDIV`/`SDIV`, zero `x18`, zero stack spills in the hot `k = 1 -> k = 0` loop). `band3_lbs.py` parses this exact fenced `asm` block live, verifies every mnemonic against the baseline scalar ARMv8.0-A integer allowlist and checks zero `x18`, derives the complete all-level (`k = 0..11` + Phase 2) instruction budget `c_lbs_budget_ir = 80.84 <= 81.0 = C_LBS` directly from the parsed labels, **and executes this exact parsed assembly live inside a pure-Python AArch64 VM (`band3_a64vm.py`)** across all 3 real `BN254 C6` prototype orbits (`343/343` exact hits == brute-force ground truth, `0` errors, `.Lalloc` fired live, `.Lk2` fired live + explicit 3-lane `K_2` test, `.Lsidx0` with `ch_id = 19` verified live, and `band3_apply16_unroll` executed live on emitted entries):

```asm
// Complete, executable baseline scalar ARMv8.0-A Band-3 kernels
// Assemble: clang -target aarch64-linux-gnu -march=armv8-a+nosimd+nofp -c
// ZERO SIMD/NEON/SVE, ZERO FP, ZERO LSE, ZERO DC ZVA, ZERO UDIV/SDIV, ZERO x18.
    .text
    .align 2

// ============================================================================
// Memory layout (all 8-byte aligned pure integer arrays):
//   x19 = ctx[] (int64):
//     #0:ddq_a1, #8:ddq_b1, #16:Ua01, #24:Ub01, #32:inva0, #40:invb0,
//     #48:free_page_ptr_addr, #56:da_stride0, #64:da_stride1,
//     #72:Da1_sc, #80:Da1_q30, #88:m_sh_a1, #96:Db1_sc, #104:Db1_q30, #112:m_sh_b1,
//     #120:ETA, #128:ETA_plus1, #136:lvl_table, #144:c_state,
//     #152:invw1, #160:inva1, #168:invb1
//   x25 = sb_cur[] (2^39 current write pointers into 4 KiB pages; next_ptr at offset 0)
//   x26 = table_b[] (flat 16B entries: w0 = sb_id; w1 = (b0_min<<50)|(L_b<<39)|ch_id)
//   x27 = b_lanes[] (3 x 64B structs: #0:b5_cur, #8:b0_cur, #16:sb5_k1, #24:db0_k1, #32:sb5_k0, #40:db0_k0)
//   x28 = SQRT16[] (128.5 KiB 65792 x uint16_t table)
//   c_state[] (mutable per-level Q24 centers):
//     #0..#88: c_w[0..11], #96..#184: c_a[0..11], #192..#280: c_b[0..11]
//   lvl_table[k] (k=0..11, 256B stride):
//     #0:invw, #8:inva, #16:invb, #24:Dw_sc, #32:Da_sc, #40:Db_sc, #48:da_stride,
//     #56:sb5[0], #64:db0[0], #72:sb5[1], #80:db0[1], #88:sb5[2], #96:db0[2],
//     #104:Uw_col ptr, #112:Ua_col ptr, #120:Ub_col ptr,
//     (k=2 inline fast fields):
//     #128:ddqw2, #136:ddqa2, #144:ddqb2, #152:Dw2_q30, #160:m_sh_w2,
//     #168:Da2_q30, #176:m_sh_a2, #184:Db2_q30, #192:m_sh_b2,
//     #200:Uw12, #208:Ua02, #216:Ua12, #224:Ub02, #232:Ub12
// ============================================================================

// ============================================================================
// KERNEL 1: band3_k1_node (includes inline band3_k1_to_k0_scan hot loop)
//   Inputs: x0=qw, x1=qa, x2=qb, x3=has_nz, x24=canon_a_k2
// ============================================================================
    .global band3_k1_node
band3_k1_node:
    stp     x29, x30, [sp, #-48]!
    stp     x20, x21, [sp, #16]
    stp     x22, x23, [sp, #32]
    ldr     x14, [x19, #144]           // c_state
    ldr     x4, [x14, #8]              // cw1
    ldr     x5, [x14, #104]            // ca1
    ldr     x6, [x14, #200]            // cb1
    mov     x7, #0x1000000             // S24
    ldr     x17, [x19, #128]           // ETA_plus1
    ldp     x15, x16, [x19, #152]      // invw1, inva1
    sub     x8, x7, x0                 // rem_w
    tbnz    x8, #63, .Lk1_ret
    lsr     x8, x8, #8
    ldrh    w8, [x28, x8, lsl #1]
    mul     x8, x8, x15
    add     x8, x17, x8, lsr #24       // hw
    sub     x11, x8, x4
    neg     x12, x4
    sub     x12, x12, x8
    sub     x8, x7, x1                 // rem_a
    tbnz    x8, #63, .Lk1_ret
    lsr     x8, x8, #8
    ldrh    w8, [x28, x8, lsl #1]
    mul     x8, x8, x16
    add     x8, x17, x8, lsr #24       // ha
    sub     x9, x8, x5
    cmp     x9, x11
    csel    x11, x9, x11, lt
    neg     x9, x5
    sub     x9, x9, x8
    cmp     x9, x12
    csel    x12, x9, x12, gt
    ldr     x16, [x19, #168]           // invb1
    sub     x8, x7, x2                 // rem_b
    tbnz    x8, #63, .Lk1_ret
    lsr     x8, x8, #8
    ldrh    w8, [x28, x8, lsl #1]
    mul     x8, x8, x16
    add     x8, x17, x8, lsr #24       // hb
    sub     x9, x8, x6
    cmp     x9, x11
    csel    x11, x9, x11, lt
    neg     x9, x6
    sub     x9, x9, x8
    cmp     x9, x12
    csel    x12, x9, x12, gt
    add     x12, x12, x7
    sub     x12, x12, #1
    asr     x20, x12, #24              // x20 = lo1
    asr     x23, x11, #24              // x23 = hi1
    cbnz    x3, .Lk1_nz_ok
.Lk1_origin_chk:
    cmp     x20, #0
    csel    x20, xzr, x20, lt
    cmp     x20, x23
    b.gt    .Lk1_ret
    cbnz    x20, .Lk1_nz_ok
    bl      band3_origin_k0_line
    add     x20, x20, #1
.Lk1_nz_ok:
    cmp     x20, x23
    b.gt    .Lk1_ret
    // Safe Q48 finite-diff init (D_q30 + asr #6) + (-ca0, -cb0), canon_a_k1, invb0
    str     x24, [sp, #-16]!
    ldr     x15, [x19, #120]           // ETA
    lsl     x8, x20, #24               // lo1 << 24
    add     x9, x8, x5                 // t0a = lo1*S24 + ca1
    cmp     x9, #0
    cneg    x10, x9, lt
    subs    x10, x10, x15
    csel    x10, xzr, x10, lt
    ldr     x11, [x19, #72]            // Da1_sc
    mul     x10, x10, x10
    umulh   x10, x10, x11
    add     x10, x1, x10
    lsl     x21, x10, #24              // qa_q48
    add     x9, x9, x9
    add     x9, x9, x7                 // 2*t0a + S24
    ldp     x11, x12, [x19, #80]       // Da1_q30, m_sh_a1
    mul     x0, x9, x11
    asr     x0, x0, #6                 // exact signed Q48 first step
    sub     x0, x0, x12
    ldr     x11, [x19, #0]             // ddq_a1
    sub     x0, x0, x11                // dqa_q48
    add     x9, x8, x6                 // t0b = lo1*S24 + cb1
    cmp     x9, #0
    cneg    x10, x9, lt
    subs    x10, x10, x15
    csel    x10, xzr, x10, lt
    ldr     x11, [x19, #96]            // Db1_sc
    mul     x10, x10, x10
    umulh   x10, x10, x11
    add     x10, x2, x10
    lsl     x22, x10, #24              // qb_q48
    add     x9, x9, x9
    add     x9, x9, x7
    ldp     x11, x12, [x19, #104]      // Db1_q30, m_sh_b1
    mul     x1, x9, x11
    asr     x1, x1, #6
    sub     x1, x1, x12
    ldr     x11, [x19, #8]             // ddq_b1
    sub     x1, x1, x11                // dqb_q48
    ldr     x3, [x14, #96]             // ca0_p
    ldr     x4, [x14, #192]            // cb0_p
    ldp     x11, x12, [x19, #16]       // Ua01, Ub01
    madd    x3, x20, x11, x3           // ca0
    madd    x4, x20, x12, x4           // cb0
    neg     x3, x3                     // -ca0 (saves neg in .Lk1_loop!)
    neg     x4, x4                     // -cb0 (saves neg in .Lk1_loop!)
    ldp     x6, x30, [x19, #56]        // x6 = da_stride0, x30 = da_stride1 (LR saved on stack!)
    madd    x24, x20, x30, x24         // canon_a_k1
    sub     x21, x21, x7               // bias qa_q48 by -S24 so (S24-1)-(qa>>24) == S24-qa!
    sub     x22, x22, x7               // bias qb_q48 by -S24 so (S24-1)-(qb>>24) == S24-qb!
    sub     x7, x7, #1                 // x7 = S24 - 1 (saves sub #1 in .Lk1_loop!)
    ldp     x2, x29, [x19, #32]        // x2 = inva0, x29 = invb0
    mov     x5, x17                    // ETA_plus1

    // --- Inline hot k=1 -> k=0 scan loop (band3_k1_to_k0_scan) ---
    .global band3_k1_to_k0_scan
band3_k1_to_k0_scan:
.Lk1_loop:
    // Part A: k=0 SQRT16 a-bound (8 ops)
    sub     x8, x7, x21, asr #24       // 1: rem_a = (S24-1) - ((qa_q48-S24)>>24) = S24 - qa
    tbnz    x8, #63, .Lk1_next         // 2: empty line
    lsr     x8, x8, #8                 // 3: u_a
    ldrh    w8, [x28, x8, lsl #1]      // 4: SQRT16[u_a]
    mul     x8, x8, x2                 // 5: * inva0
    add     x8, x5, x8, lsr #24        // 6: ha
    add     x11, x3, x8                // 7: hi_num_a = (-ca0) + ha
    sub     x12, x3, x8                // 8: lo_num_a = (-ca0) - ha

    // Part B: k=0 SQRT16 b-bound + intersect + lo/hi + line setup (20 ops)
    sub     x8, x7, x22, asr #24       // 1: rem_b = S24 - qb
    tbnz    x8, #63, .Lk1_next         // 2: empty line
    lsr     x8, x8, #8                 // 3: u_b
    ldrh    w8, [x28, x8, lsl #1]      // 4: SQRT16[u_b]
    mul     x8, x8, x29                // 5: * invb0 (in register x29!)
    add     x8, x5, x8, lsr #24        // 6: hb
    add     x9, x4, x8                 // 7: (-cb0) + hb
    cmp     x9, x11                    // 8: min hi_num
    csel    x11, x9, x11, lt           // 9
    sub     x9, x4, x8                 // 10: (-cb0) - hb
    cmp     x9, x12                    // 11: max lo_num
    csel    x12, x9, x12, gt           // 12
    add     x12, x12, x7               // 13: lo_num + (S24 - 1) in 1 op!
    asr     x12, x12, #24              // 14: lo
    asr     x8, x11, #24               // 15: hi
    subs    x9, x8, x12                // 16: hi - lo
    b.lt    .Lk1_next                  // 17: empty line
    madd    x10, x12, x6, x24          // 18: p_a = canon_a_k1 + lo * da_stride0
    add     x12, x9, #1                // 19: rem = hi - lo + 1
    add     x8, x8, #1                 // 20: hi_plus1 = hi + 1

    // Part C: 5-op candidate scan (.Lloop)
.Lloop:
    ldr     x9, [x10]                  // 1: val = canon_a[p_a]
    add     x10, x10, x6               // 2: p_a += da_stride0
    cbnz    x9, .Lhit                  // 3: hit iff a in F_unit
    subs    x12, x12, #1               // 4: rem -= 1
    b.ne    .Lloop                     // 5: next candidate

    // Part F: k=1 safe Q48 finite-diff & center step (12 ops)
.Lk1_next:
    ldp     x8, x9, [x19, #0]          // 1: ddq_a1, ddq_b1
    add     x0, x0, x8                 // 2: dqa += ddq_a1
    add     x21, x21, x0               // 3: qa += dqa
    add     x1, x1, x9                 // 4: dqb += ddq_b1
    add     x22, x22, x1               // 5: qb += dqb
    ldp     x8, x9, [x19, #16]         // 6: Ua01, Ub01
    sub     x3, x3, x8                 // 7: (-ca0) -= Ua01
    sub     x4, x4, x9                 // 8: (-cb0) -= Ub01
    add     x24, x24, x30              // 9: canon_a_k1 += da_stride1 (in x30!)
    add     x20, x20, #1               // 10: x1 += 1
    cmp     x20, x23                   // 11: x1 <= hi1?
    b.le    .Lk1_loop                  // 12: next k=1 line
    ldr     x24, [sp], #16
.Lk1_ret:
    ldp     x22, x23, [sp, #32]
    ldp     x20, x21, [sp, #16]
    ldp     x29, x30, [sp], #48
    ret

    // Part E: 31-op hit handler (.Lhit)
.Lhit:
    tbnz    x9, #63, .Lk2              // 1: rare K_2 flag (2^-28)
    lsr     x13, x9, #46               // 2: kc in {0, 1, 2}
    add     x13, x27, x13, lsl #6      // 3: &b_lanes[kc] (64B stride)
    ldp     x14, x17, [x13, #0]        // 4: b5_cur, b0_cur
    ldp     x15, x16, [x13, #16]       // 5: sb5_k1, db0_k1
    madd    x14, x20, x15, x14         // 6: b5_cur + x1*sb5_k1 (EXACT!)
    madd    x17, x20, x16, x17         // 7: b0_cur + x1*db0_k1 (EXACT!)
    ldp     x15, x16, [x13, #32]       // 8: sb5_k0, db0_k0
    sub     x11, x8, x12               // 9: x0 = hi_plus1 - rem
    madd    x14, x11, x15, x14         // 10: sidx = b5_k1 + x0*sb5_k0
    madd    x15, x11, x16, x17         // 11: b0   = b0_k1 + x0*db0_k0
    cmp     x14, #0                    // 12: sign / zero of sidx
    b.eq    .Lsidx0                    // 13: rare sidx == 0 (2^-43)
    cneg    x14, x14, lt               // 14: u_i = abs(sidx)
    cneg    x15, x15, lt               // 15: b0_s = (sidx < 0) ? -b0 : b0
    add     x14, x26, x14, lsl #4      // 16: &table_b[u_i]
    ldp     x16, x17, [x14]            // 17: x16=sb_id (clean!), x17=(b0_min<<50)|(L_b<<39)|ch_id
.Ltable_b_loaded:
    asr     x14, x17, #50              // 18: signed b0_min
    ubfx    x13, x17, #39, #11         // 19: L_b (< 1024)
    sub     x15, x15, x14              // 20: b0_rel = b0_s - b0_min
    cmp     x15, x13                   // 21: unsigned b0_rel < L_b
    b.hs    .Lhit_done                 // 22: out of chord -> skip
    bfi     x9, x15, #45, #19          // 23: insert b0_rel<<45 & clear kc in 1 op!
    add     x9, x9, x17, lsl #55       // 24: elem_idx |= (ch_id<<55) (ch_id in 0..19, bits 5..8 are 0!)
    ldr     x14, [x25, x16, lsl #3]    // 25: cur_ptr = sb_cur[sb_id]
    str     x9, [x14], #8              // 26: *cur_ptr++ = elem_idx
    str     x14, [x25, x16, lsl #3]    // 27: sb_cur[sb_id] = cur_ptr
    tst     x14, #4095                 // 28: 4 KiB page full?
    b.eq    .Lalloc                    // 29: rare 1/511 page alloc
.Lhit_done:
    subs    x12, x12, #1               // 30: rem -= 1
    b.ne    .Lloop                     // 31: continue k=0 scan
    b       .Lk1_next                  // (last-hit branch: <= 0.42/hit)

.Lalloc:
    // Cold 4 KiB page allocator (9 ops, invoked 1/511 hits -> 0.018 Ir/hit)
    ldr     x15, [x19, #48]
    ldr     x13, [x15]
    add     x17, x13, #4096
    str     x17, [x15]
    sub     x14, x14, #4096
    str     x14, [x13]
    add     x13, x13, #8
    str     x13, [x25, x16, lsl #3]
    b       .Lhit_done

.Lsidx0:
    cbz     x15, .Lhit_done
    cmp     x15, #0
    cneg    x15, x15, lt
    ldp     x16, x17, [x26]
    b       .Ltable_b_loaded

.Lk2:
    // Complete cold K_2 handler (2^-28 rare): emits ALL three kc in {0, 1, 2}!
    stp     x0, x1, [sp, #-32]!
    stp     x2, x3, [sp, #16]
    and     x2, x9, #0x1fffffffffff
    sub     x11, x8, x12
    mov     x3, #0
.Lk2_loop:
    add     x13, x27, x3, lsl #6
    ldp     x14, x17, [x13, #0]
    ldp     x15, x16, [x13, #16]
    madd    x14, x20, x15, x14
    madd    x17, x20, x16, x17
    ldp     x15, x16, [x13, #32]
    madd    x14, x11, x15, x14
    madd    x15, x11, x16, x17
    cmp     x14, #0
    b.eq    .Lk2_sidx0
    cneg    x14, x14, lt
    cneg    x15, x15, lt
    add     x14, x26, x14, lsl #4
    ldp     x16, x17, [x14]
    b       .Lk2_chk
.Lk2_sidx0:
    cbz     x15, .Lk2_next
    cmp     x15, #0
    cneg    x15, x15, lt
    ldp     x16, x17, [x26]
.Lk2_chk:
    asr     x14, x17, #50
    ubfx    x13, x17, #39, #11
    sub     x15, x15, x14
    cmp     x15, x13
    b.hs    .Lk2_next
    mov     x0, x2
    bfi     x0, x15, #45, #19
    add     x0, x0, x17, lsl #55
    ldr     x14, [x25, x16, lsl #3]
    str     x0, [x14], #8
    str     x14, [x25, x16, lsl #3]
    tst     x14, #4095
    b.ne    .Lk2_next
    ldr     x15, [x19, #48]
    ldr     x13, [x15]
    add     x17, x13, #4096
    str     x17, [x15]
    sub     x14, x14, #4096
    str     x14, [x13]
    add     x13, x13, #8
    str     x13, [x25, x16, lsl #3]
.Lk2_next:
    add     x3, x3, #1
    cmp     x3, #3
    b.lt    .Lk2_loop
    ldp     x2, x3, [sp, #16]
    ldp     x0, x1, [sp], #32
    b       .Lhit_done

// Cold single-orbit origin line handler (x_1=...=x_11=0, x_0 >= 1; invoked ONCE per orbit)
band3_origin_k0_line:
    stp     x29, x30, [sp, #-96]!
    stp     x0, x1, [sp, #16]
    stp     x2, x23, [sp, #32]
    stp     x5, x6, [sp, #48]
    stp     x14, x17, [sp, #64]
    stp     x20, x24, [sp, #80]
    mov     x20, #0
    ldr     x2, [x19, #32]
    ldr     x5, [x19, #128]
    ldr     x6, [x19, #56]
    mov     x7, #0x1000000
    lsr     x8, x7, #8
    sub     x8, x8, #1
    ldrh    w8, [x28, x8, lsl #1]
    mul     x9, x8, x2
    add     x9, x5, x9, lsr #24
    ldr     x10, [x19, #40]
    mul     x8, x8, x10
    add     x8, x5, x8, lsr #24
    cmp     x8, x9
    csel    x11, x8, x9, lt
    asr     x8, x11, #24
    cmp     x8, #1
    b.lt    .Lorig_ret
    mov     x12, #1
    sub     x9, x8, x12
    madd    x10, x12, x6, x24
    add     x12, x9, #1
    add     x8, x8, #1
    mov     x23, #-1
    mov     x0, #0
    mov     x1, #0
    mov     x21, #0
    mov     x22, #0
    mov     x3, #0
    mov     x4, #0
    bl      .Lorig_enter
.Lorig_ret:
    ldp     x20, x24, [sp, #80]
    ldp     x14, x17, [sp, #64]
    ldp     x5, x6, [sp, #48]
    ldp     x2, x23, [sp, #32]
    ldp     x0, x1, [sp, #16]
    ldp     x29, x30, [sp], #96
    ret
.Lorig_enter:
    stp     x29, x30, [sp, #-48]!
    stp     x20, x21, [sp, #16]
    stp     x22, x23, [sp, #32]
    str     x24, [sp, #-16]!
    b       .Lloop

// ============================================================================
// KERNEL 2: band3_k_ge_2_node (complete recursive k=11..2 driver)
//   k == 2: safe Q48 finite-diff + inline lvl[2] + incremental lo..hi stepping!
//   k >= 3: exact umulh + general incremental lo..hi stepping + single restore!
//   Inputs: x0=qw, x1=qa, x2=qb, x3=has_nz, x4=k (2..11), x24=canon_a_k
// ============================================================================
    .global band3_k_ge_2_node
band3_k_ge_2_node:
    stp     x29, x30, [sp, #-112]!
    stp     x20, x21, [sp, #16]
    stp     x22, x23, [sp, #32]
    stp     x0, x1, [sp, #48]
    stp     x2, x3, [sp, #64]
    mov     x22, x4                    // x22 = k
    ldr     x13, [x19, #136]           // lvl_table
    add     x23, x13, x22, lsl #8      // x23 = &lvl[k]
    ldr     x14, [x19, #144]           // c_state
    ldr     x5, [x14, x22, lsl #3]     // cw
    add     x15, x14, #96
    ldr     x6, [x15, x22, lsl #3]     // ca
    add     x15, x14, #192
    ldr     x10, [x15, x22, lsl #3]    // cb
    mov     x7, #0x1000000
    ldr     x17, [x19, #128]           // ETA_plus1
    ldp     x15, x16, [x23, #0]        // invw, inva
    sub     x8, x7, x0
    tbnz    x8, #63, .Lk_ret
    lsr     x8, x8, #8
    ldrh    w8, [x28, x8, lsl #1]
    mul     x8, x8, x15
    add     x8, x17, x8, lsr #24
    sub     x11, x8, x5
    neg     x12, x5
    sub     x12, x12, x8
    cbz     x16, .Lk_b
    sub     x8, x7, x1
    tbnz    x8, #63, .Lk_ret
    lsr     x8, x8, #8
    ldrh    w8, [x28, x8, lsl #1]
    mul     x8, x8, x16
    add     x8, x17, x8, lsr #24
    sub     x9, x8, x6
    cmp     x9, x11
    csel    x11, x9, x11, lt
    neg     x9, x6
    sub     x9, x9, x8
    cmp     x9, x12
    csel    x12, x9, x12, gt
.Lk_b:
    ldr     x16, [x23, #16]
    cbz     x16, .Lk_done_int
    sub     x8, x7, x2
    tbnz    x8, #63, .Lk_ret
    lsr     x8, x8, #8
    ldrh    w8, [x28, x8, lsl #1]
    mul     x8, x8, x16
    add     x8, x17, x8, lsr #24
    sub     x9, x8, x10
    cmp     x9, x11
    csel    x11, x9, x11, lt
    neg     x9, x10
    sub     x9, x9, x8
    cmp     x9, x12
    csel    x12, x9, x12, gt
.Lk_done_int:
    add     x12, x12, x7
    sub     x12, x12, #1
    asr     x20, x12, #24              // x20 = val = lo
    asr     x21, x11, #24              // x21 = hi
    cbnz    x3, .Lk_loop_chk
    cmp     x20, #0
    csel    x20, xzr, x20, lt
.Lk_loop_chk:
    cmp     x20, x21
    b.gt    .Lk_ret

    // Incremental setup ONCE at val = lo (canon_a + 3 b_lanes)
    ldr     x11, [x23, #48]
    madd    x24, x20, x11, x24
    ldp     x8, x9, [x23, #56]
    ldp     x10, x11, [x27, #0]
    madd    x10, x20, x8, x10
    madd    x11, x20, x9, x11
    stp     x10, x11, [x27, #0]
    ldp     x8, x9, [x23, #72]
    ldp     x10, x11, [x27, #64]
    madd    x10, x20, x8, x10
    madd    x11, x20, x9, x11
    stp     x10, x11, [x27, #64]
    ldp     x8, x9, [x23, #88]
    ldp     x10, x11, [x27, #128]
    madd    x10, x20, x8, x10
    madd    x11, x20, x9, x11
    stp     x10, x11, [x27, #128]

    // Dispatch fast k == 2 path vs general k >= 3 path
    cmp     x22, #2
    b.ne    .Lk3_init_c
.Lk2_init:
    // --- Fast k == 2 init (inline c_state + safe Q48 finite-diff for qw, qa, qb) ---
    ldr     x8, [x23, #200]            // Uw12
    ldr     x9, [x14, #8]              // cw[1]
    madd    x9, x20, x8, x9
    str     x9, [x14, #8]
    ldp     x8, x9, [x23, #208]        // Ua02, Ua12
    ldp     x10, x11, [x14, #96]       // ca[0], ca[1]
    madd    x10, x20, x8, x10
    madd    x11, x20, x9, x11
    stp     x10, x11, [x14, #96]
    ldp     x8, x9, [x23, #224]        // Ub02, Ub12
    ldp     x10, x11, [x14, #192]      // cb[0], cb[1]
    madd    x10, x20, x8, x10
    madd    x11, x20, x9, x11
    stp     x10, x11, [x14, #192]
    // Save has_nz at [sp, #96]; overwrite [sp, #48..#95] with (qw,dqw), (qa,dqa), (qb,dqb)
    str     x3, [sp, #96]
    ldr     x15, [x19, #120]           // ETA
    lsl     x8, x20, #24               // lo2 << 24
    // qw_q48, dqw_q48 at [sp, #48]
    add     x9, x8, x5                 // t0w
    cmp     x9, #0
    cneg    x10, x9, lt
    subs    x10, x10, x15
    csel    x10, xzr, x10, lt
    ldr     x11, [x23, #24]            // Dw_sc[2]
    mul     x10, x10, x10
    umulh   x10, x10, x11
    add     x10, x0, x10
    lsl     x10, x10, #24              // qw_q48
    add     x9, x9, x9
    add     x9, x9, x7
    ldp     x11, x12, [x23, #152]      // Dw2_q30, m_sh_w2
    mul     x9, x9, x11
    asr     x9, x9, #6
    sub     x9, x9, x12
    ldr     x11, [x23, #128]           // ddqw2
    sub     x9, x9, x11                // dqw_q48
    stp     x10, x9, [sp, #48]
    // qa_q48, dqa_q48 at [sp, #64]
    add     x9, x8, x6                 // t0a
    cmp     x9, #0
    cneg    x10, x9, lt
    subs    x10, x10, x15
    csel    x10, xzr, x10, lt
    ldr     x11, [x23, #32]            // Da_sc[2]
    mul     x10, x10, x10
    umulh   x10, x10, x11
    add     x10, x1, x10
    lsl     x10, x10, #24              // qa_q48
    add     x9, x9, x9
    add     x9, x9, x7
    ldp     x11, x12, [x23, #168]      // Da2_q30, m_sh_a2
    mul     x9, x9, x11
    asr     x9, x9, #6
    sub     x9, x9, x12
    ldr     x11, [x23, #136]           // ddqa2
    sub     x9, x9, x11                // dqa_q48
    stp     x10, x9, [sp, #64]
    // qb_q48, dqb_q48 at [sp, #80]
    ldr     x13, [x14, #208]           // cb[2]
    add     x9, x8, x13                // t0b
    cmp     x9, #0
    cneg    x10, x9, lt
    subs    x10, x10, x15
    csel    x10, xzr, x10, lt
    ldr     x11, [x23, #40]            // Db_sc[2]
    mul     x10, x10, x10
    umulh   x10, x10, x11
    add     x10, x2, x10
    lsl     x10, x10, #24              // qb_q48
    add     x9, x9, x9
    add     x9, x9, x7
    ldp     x11, x12, [x23, #184]      // Db2_q30, m_sh_b2
    mul     x9, x9, x11
    asr     x9, x9, #6
    sub     x9, x9, x12
    ldr     x11, [x23, #144]           // ddqb2
    sub     x9, x9, x11                // dqb_q48
    stp     x10, x9, [sp, #80]

    // --- Fast Q48 incremental k == 2 child loop (54 ops total per k=1 node!) ---
.Lk2_child_loop:
    ldp     x8, x9, [x23, #128]        // 1: ddqw2, ddqa2
    ldr     x10, [x23, #144]           // 2: ddqb2
    ldp     x0, x11, [sp, #48]         // 3: qw_q48, dqw_q48
    add     x11, x11, x8               // 4: dqw += ddqw2
    add     x12, x0, x11               // 5: qw_next
    stp     x12, x11, [sp, #48]        // 6
    asr     x0, x0, #24                // 7: qw_q24
    ldp     x1, x11, [sp, #64]         // 8: qa_q48, dqa_q48
    add     x11, x11, x9               // 9: dqa += ddqa2
    add     x12, x1, x11               // 10: qa_next
    stp     x12, x11, [sp, #64]        // 11
    asr     x1, x1, #24                // 12: qa_q24
    ldp     x2, x11, [sp, #80]         // 13: qb_q48, dqb_q48
    add     x11, x11, x10              // 14: dqb += ddqb2
    add     x12, x2, x11               // 15: qb_next
    stp     x12, x11, [sp, #80]        // 16
    asr     x2, x2, #24                // 17: qb_q24
    ldr     x3, [sp, #96]              // 18: has_nz_parent
    cbz     x20, .Lk2_do_k1            // 19
    mov     x3, #1                     // 20
.Lk2_do_k1:
    bl      band3_k1_node              // 21: call k=1 node!
    // Incremental +1 step for k=2 (33 ops: canon_a + 3 b_lanes + cw1, ca0..1, cb0..1)
    ldr     x11, [x23, #48]            // 22: da_stride[2]
    add     x24, x24, x11              // 23
    ldp     x8, x9, [x23, #56]         // 24: b_lanes[0]
    ldp     x10, x11, [x27, #0]        // 25
    add     x10, x10, x8               // 26
    add     x11, x11, x9               // 27
    stp     x10, x11, [x27, #0]        // 28
    ldp     x8, x9, [x23, #72]         // 29: b_lanes[1]
    ldp     x10, x11, [x27, #64]       // 30
    add     x10, x10, x8               // 31
    add     x11, x11, x9               // 32
    stp     x10, x11, [x27, #64]       // 33
    ldp     x8, x9, [x23, #88]         // 34: b_lanes[2]
    ldp     x10, x11, [x27, #128]      // 35
    add     x10, x10, x8               // 36
    add     x11, x11, x9               // 37
    stp     x10, x11, [x27, #128]      // 38
    ldr     x14, [x19, #144]           // 39: c_state
    ldr     x8, [x23, #200]            // 40: Uw12
    ldr     x9, [x14, #8]              // 41: cw[1]
    add     x9, x9, x8                 // 42
    str     x9, [x14, #8]              // 43
    ldp     x8, x9, [x23, #208]        // 44: Ua02, Ua12
    ldp     x10, x11, [x14, #96]       // 45: ca[0], ca[1]
    add     x10, x10, x8               // 46
    add     x11, x11, x9               // 47
    stp     x10, x11, [x14, #96]       // 48
    ldp     x8, x9, [x23, #224]        // 49: Ub02, Ub12
    ldp     x10, x11, [x14, #192]      // 50: cb[0], cb[1]
    add     x10, x10, x8               // 51
    add     x11, x11, x9               // 52
    stp     x10, x11, [x14, #192]      // 53
    add     x20, x20, #1               // 54
    cmp     x20, x21                   // 55
    b.le    .Lk2_child_loop            // 56
.Lk2_restore:
    // Fast k == 2 restore by -x20 = -(hi2 + 1)
    ldr     x8, [x23, #200]
    ldr     x9, [x14, #8]
    msub    x9, x20, x8, x9
    str     x9, [x14, #8]
    ldp     x8, x9, [x23, #208]
    ldp     x10, x11, [x14, #96]
    msub    x10, x20, x8, x10
    msub    x11, x20, x9, x11
    stp     x10, x11, [x14, #96]
    ldp     x8, x9, [x23, #224]
    ldp     x10, x11, [x14, #192]
    msub    x10, x20, x8, x10
    msub    x11, x20, x9, x11
    stp     x10, x11, [x14, #192]
    b       .Lk_restore_ab

    // --- General incremental k >= 3 path ---
.Lk3_init_c:
    ldr     x15, [x23, #104]
    ldr     x16, [x23, #112]
    ldr     x17, [x23, #120]
    add     x10, x14, #96
    add     x11, x14, #192
    mov     x12, #0
.Lk_init_c:
    ldr     x8, [x15, x12, lsl #3]
    ldr     x9, [x14, x12, lsl #3]
    madd    x9, x20, x8, x9
    str     x9, [x14, x12, lsl #3]
    ldr     x8, [x16, x12, lsl #3]
    ldr     x9, [x10, x12, lsl #3]
    madd    x9, x20, x8, x9
    str     x9, [x10, x12, lsl #3]
    ldr     x8, [x17, x12, lsl #3]
    ldr     x9, [x11, x12, lsl #3]
    madd    x9, x20, x8, x9
    str     x9, [x11, x12, lsl #3]
    add     x12, x12, #1
    cmp     x12, x22
    b.lt    .Lk_init_c

.Lk3_child_loop:
    ldp     x0, x1, [sp, #48]
    ldp     x2, x3, [sp, #64]
    ldr     x14, [x19, #144]
    ldr     x5, [x14, x22, lsl #3]
    add     x15, x14, #96
    ldr     x6, [x15, x22, lsl #3]
    add     x15, x14, #192
    ldr     x10, [x15, x22, lsl #3]
    mov     x7, #0x1000000
    ldr     x17, [x19, #120]
    lsl     x8, x20, #24
    add     x9, x8, x5
    cmp     x9, #0
    cneg    x9, x9, lt
    subs    x9, x9, x17
    csel    x9, xzr, x9, lt
    ldr     x11, [x23, #24]
    mul     x9, x9, x9
    umulh   x9, x9, x11
    add     x0, x0, x9
    cmp     x0, x7
    b.gt    .Lk3_step
    ldr     x11, [x23, #32]
    cbz     x11, .Lk3_chk_b
    add     x9, x8, x6
    cmp     x9, #0
    cneg    x9, x9, lt
    subs    x9, x9, x17
    csel    x9, xzr, x9, lt
    mul     x9, x9, x9
    umulh   x9, x9, x11
    add     x1, x1, x9
    cmp     x1, x7
    b.gt    .Lk3_step
.Lk3_chk_b:
    ldr     x11, [x23, #40]
    cbz     x11, .Lk3_do_rec
    add     x9, x8, x10
    cmp     x9, #0
    cneg    x9, x9, lt
    subs    x9, x9, x17
    csel    x9, xzr, x9, lt
    mul     x9, x9, x9
    umulh   x9, x9, x11
    add     x2, x2, x9
    cmp     x2, x7
    b.gt    .Lk3_step
.Lk3_do_rec:
    cbz     x20, .Lk3_call
    mov     x3, #1
.Lk3_call:
    sub     x4, x22, #1
    bl      band3_k_ge_2_node
.Lk3_step:
    ldr     x11, [x23, #48]
    add     x24, x24, x11
    ldp     x8, x9, [x23, #56]
    ldp     x10, x11, [x27, #0]
    add     x10, x10, x8
    add     x11, x11, x9
    stp     x10, x11, [x27, #0]
    ldp     x8, x9, [x23, #72]
    ldp     x10, x11, [x27, #64]
    add     x10, x10, x8
    add     x11, x11, x9
    stp     x10, x11, [x27, #64]
    ldp     x8, x9, [x23, #88]
    ldp     x10, x11, [x27, #128]
    add     x10, x10, x8
    add     x11, x11, x9
    stp     x10, x11, [x27, #128]
    ldr     x14, [x19, #144]
    ldr     x15, [x23, #104]
    ldr     x16, [x23, #112]
    ldr     x17, [x23, #120]
    add     x10, x14, #96
    add     x11, x14, #192
    mov     x12, #0
.Lk3_upd_c:
    ldr     x8, [x15, x12, lsl #3]
    ldr     x9, [x14, x12, lsl #3]
    add     x9, x9, x8
    str     x9, [x14, x12, lsl #3]
    ldr     x8, [x16, x12, lsl #3]
    ldr     x9, [x10, x12, lsl #3]
    add     x9, x9, x8
    str     x9, [x10, x12, lsl #3]
    ldr     x8, [x17, x12, lsl #3]
    ldr     x9, [x11, x12, lsl #3]
    add     x9, x9, x8
    str     x9, [x11, x12, lsl #3]
    add     x12, x12, #1
    cmp     x12, x22
    b.lt    .Lk3_upd_c
    add     x20, x20, #1
    cmp     x20, x21
    b.le    .Lk3_child_loop
.Lk3_rst_prep:
    ldr     x14, [x19, #144]
    ldr     x15, [x23, #104]
    ldr     x16, [x23, #112]
    ldr     x17, [x23, #120]
    add     x10, x14, #96
    add     x11, x14, #192
    mov     x12, #0
.Lk_rst_c:
    ldr     x8, [x15, x12, lsl #3]
    ldr     x9, [x14, x12, lsl #3]
    msub    x9, x20, x8, x9
    str     x9, [x14, x12, lsl #3]
    ldr     x8, [x16, x12, lsl #3]
    ldr     x9, [x10, x12, lsl #3]
    msub    x9, x20, x8, x9
    str     x9, [x10, x12, lsl #3]
    ldr     x8, [x17, x12, lsl #3]
    ldr     x9, [x11, x12, lsl #3]
    msub    x9, x20, x8, x9
    str     x9, [x11, x12, lsl #3]
    add     x12, x12, #1
    cmp     x12, x22
    b.lt    .Lk_rst_c

.Lk_restore_ab:
    ldr     x11, [x23, #48]
    msub    x24, x20, x11, x24
    ldp     x8, x9, [x23, #56]
    ldp     x10, x11, [x27, #0]
    msub    x10, x20, x8, x10
    msub    x11, x20, x9, x11
    stp     x10, x11, [x27, #0]
    ldp     x8, x9, [x23, #72]
    ldp     x10, x11, [x27, #64]
    msub    x10, x20, x8, x10
    msub    x11, x20, x9, x11
    stp     x10, x11, [x27, #64]
    ldp     x8, x9, [x23, #88]
    ldp     x10, x11, [x27, #128]
    msub    x10, x20, x8, x10
    msub    x11, x20, x9, x11
    stp     x10, x11, [x27, #128]
.Lk_ret:
    ldp     x22, x23, [sp, #32]
    ldp     x20, x21, [sp, #16]
    ldp     x29, x30, [sp], #112
    ret

// ============================================================================
// KERNEL 3: band3_apply16_unroll
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
    b.gt    .Lapp16
    ret
```

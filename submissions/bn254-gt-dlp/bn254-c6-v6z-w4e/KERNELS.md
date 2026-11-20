# `KERNELS.md` — 100% Pure-Integer Baseline Scalar AArch64 Kernels & Shoup `< 2^62` 3-Prime NTT (`-mgeneral-regs-only`, `simd = 0`, Zero `fp64`, Zero `FEAT_LSE`, Zero `DC ZVA`)

Every inner loop and geometric predicate in **`BN254 v6z-w4e`** uses **100% baseline scalar ARMv8.0-A integer instructions** (`-mgeneral-regs-only`, `simd = 0`, zero NEON/LSE/`DC ZVA`/`fp64`). Collection uses **Band 0 `[2^14, 2^30)`** local-owner direct-apply (`APPLY8 = 3.75` + `24.0` tail), **Band 1 `[2^30, 2^44.5)`** split into dense `emit4` `[2^30, 2^41)` (`C_HIT_LO_DENSE = 7.50` + explicit `+24.0` tail) and sparse `[2^41, 2^44.5)` (`C_HIT_LO_SPARSE = 12.00`, zero tail), and **8B-packed Band 2 `[2^44.5, 2^54.75)`** (`W_SB = 128`, `max_sb = 8 <= 8`, `C_SWAR_PROBE = 12.0`, `C_SB_ENTER = 16.0`, `C_CAND = 12.0`, `C_HIT_HI = 12.00`, explicit zero-dict zero-sort 2-pass `a_0`-run + ascending-`a_slot` scatter lists, ZERO Band 3).

## 0. Exact Pure-Integer Geometry & 2-Stage Pure-Int Bareiss 12D LLL (`zero_float_geom = True`)

1. **Exact `sigma`-Invariant `G_sym` & Voronoi `D_fi` (`canon`, `330` Ir/raw `a` `<< 2^42.8` in `canon <= 2^74.49`):**
   - `G_sym = G_0 + S^T G_0 S + (S^2)^T G_0 S^2` (`G_0 = round(2^16*GA)` symm.; `G_SC = 196,608`; `S^T G_sym S == G_sym` in `Z^{6x6}`, rel diff `7.17e-6` vs `GA`). Both `a, b` use `Q(v) = v^T G_sym v < 2^44 << 2^63`; `D_fi = U_fi^T G_sym U_fi - G_sym in Z^{6x6}` (`fi = 0..5`).
   - **Maximal-Order + Index-1 + Convex-Hexagon Proof (ALL `a in O_h\{0}`, ALL radii, LIVE in `audit_bn_cm.py [exhaustive_fd]`):**
     - *Dedekind `O_lat == O_h`:* `O_lat = Z*(ZK_NUM/3)`, `[O_lat : Z[y]] = 3`, `disc(O_lat) = -2^9 * 19^2` (`3 ∤ disc`); mod `2` and mod `19` Dedekind criteria pass (`max_order_dedekind_ok = True`).
     - *Index-1 (`[O_h^* : <-1,u_1,u_2>] == 1`):* Exact `Fraction` `LDL^T` enumerates ALL `5,670` nonzero points `Q(v) <= 50*G_SC` (`382/382` orbits, `0` missed/dups); outward-rounded `Decimal(55)` bounds coset reps to `Q(v) <= ceil(7.1888*G_SC) < 50*G_SC`, finding only `{+-1}` (`unit_index1_proved = True`).
     - *Convex hexagon:* Outward-rounded `Decimal(55)` Rouché disks (`|y - z_j| <= 10^-40`), interval `LDL^T` (`eps_G_hi = 2.588e-5`), and `P_6^(1e-4)` prove `min_lo L_u(V_k) = +0.449890 > 0` across all non-facet units in `|m_1|,|m_2| <= 6` (`global_hexagon_proof_ok = True`).
   - Cost: `Q(a)` (`42`) + `6*D_fi` (`252`) + `6` ties (`12`) + `SMAT` (`24`) = **`330` Ir/raw `a`**.
2. **Exact Active Cutoffs via `Q(a)`-Sorted `a_slot` & `M_t[t]` (`500*T_a*Q(b) + 161*T_b*Q(a) <= 500*T_a*T_b`, ZERO `qa_u8`):**
   - Canonical `a` (`N_a = 2^45.9755`) sorted ONCE by `(Q(a), a)` inside `ins_pre`; `M_t[t] = bisect_right(Q_a_sorted, qa_max_t[t])` (`<= 6.8 KiB` L1) gives exact IFF `a_slot < M_t[t] <=> Q(a_slot) <= qa_max_t[t]`, with 0-FN `sb_u7` pre-filter (`1523/1523`, `0` mismatches live).
3. **Unified 2-Stage Pure-Int Bareiss 12D LLL (`11/40 -> 9/10`, `bareiss_lll_12d` EVERYWHERE: `<= 2^40` Boot & `<= 2^36` Special-`q` up to `q <= 2^220`):**
   - Pure integer Bareiss reducer (`11/40` pre-reduction followed by `9/10` Lovasz certification, `verify_bareiss_910_cert == True`, zero floats, `inductive_sr_inv_ok = True` across BOTH stages, `exact_div_errors = 0`). Because every Stage-1 swap drops `D` by `< 11/40 < 9/10` and every Stage-2 swap drops `D` by `< 9/10`, total swaps satisfy `S_1 + S_2 <= S_max(D_0) = 185,076` strictly from `D_0`, giving joint 2-stage non-amortized worst-case cost: Boot (`det = p^6`) **`2^{39.81} < 2^{40}` Ir** (`two_stage_ok = True`, measured `2^{29.48}`, `2,380` total swaps), Special-`q` (`q <= 2^{220}`) `2^{33.16} << 2^{36} = LLL_NODE` Ir (live on `boot[0]`, all `24/24` prod special-`q`, `sq220_wc`, `8/8` boot certs, and 224-ideal `C6` campaign).

## 1. Pure-Integer Band-0 + Band-1 + Band-2 Sieve, Zero-Sort Flat-Array List Build, & 32-Cell `ORR` Scan+Zero

0. **Shared Buffered Transport & Unrolled-8 Apply (`COPY64 = 1.25`, `APPLY8 = 3.75` Ir/hit):**
   - **`COPY64` (`1.25` Ir/hit):** 8-hit (64B) cache-line copy from per-owner hit buffer to score-arena inbox: `4 ldp` + `4 stp` + `subs; b.ne` = **`10 / 8 = 1.25` Ir/hit**.
   - **`APPLY8` (`3.75` Ir/hit, used directly in Band 0 and on Band-1/2 inboxes):** 8x-unrolled `uint16_t` score add (`w_lp` in register): `4 ldp` (`4` ops loading 8 64-bit linear arena indices) + `8 * (ldrh w_k, [x_ar, x_i, lsl #1]; add w_k, w_k, w_lp; strh w_k, [x_ar, x_i, lsl #1])` (`24` ops) + `subs x_cnt, x_cnt, #1; b.ne .Lapp` (`2` ops) = **`30 / 8 = 3.75` exact baseline scalar ARMv8.0-A opcodes/hit**!
1. **Band-0 Local-Owner Direct-Apply (`[2^14, 2^30)`, `B_a = 2014` owners of `N_a_blk = 2^35` slots, `C_HIT_B0 = 3.75` + `C_B0_TAIL = 24.0`, `2^62.72 B << 2^64 B`):**
   - Each block owner runs binary search for `M_t_blk[t]` (`24 + 8*(log2(mu+1)+1)`) + `+24.0` Ir unconditional entry/tail + 8x-unrolled `APPLY8` directly into local `score_b0[t][a_loc]` (`30 / 8 = 3.75` Ir/hit, zero buffering; 2-block flat-array test `0` mismatches live).
2. **Band-1 Dense/Sparse Split (`[2^30, 2^44.5)`, parallel `8B keys[] + 8B a_slots[]`):**
   - **Dense sub-band `[2^30, 2^41)` (`mu = N_a / p >= 31`):** 4x-unrolled buffered `emit4` (`2 ldp` + `4 add x_rb` + `2 stp` + `subs; b.ne` = **`10 / 4 = 2.50` Ir/hit**) + explicit unconditional **`C_B1_TAIL = 24.0` Ir per `(p, b)` row** (for `< 4` tail + buffer flush check, added directly to `tot_lk` in `model_bn254_v6.py`) + `COPY64` (`1.25`) + `APPLY8` (`3.75`) = **`C_HIT_LO_DENSE = 7.50` Ir/hit**.
   - **Sparse sub-band `[2^41, 2^44.5)`:** Direct 1-hit flush-checked 64-bit emit (`ldr x_a,[x_p],#8; add x_h,x_rb,x_a; str x_h,[x_buf],#8; cmp; b.eq flush; subs; b.ne` = **`7.0` Ir/hit**, zero `+24.0` tail) + `COPY64` (`1.25`) + `APPLY8` (`3.75`) = **`C_HIT_LO_SPARSE = 12.00` Ir/hit**.
   - Both sub-bands pay exact `min(bsearch, 11 + 14*(1 - exp(-mu*fa)) + mu*fa)` lookup (`c_e4`) + `24.0` Ir row setup (`56`-bit `a_slot >= 2^32` verified live).
3. **Band-2 8-Lane 64-Bit SWAR `sb_u7` Probe (`W_SB = 128`, `max_sb = 8 <= 8`, `C_SWAR_PROBE = 12.0`, `C_SB_ENTER = 16.0`, `C_CAND = 12.0`, `C_HIT_HI = 12.00`):**
   - Since `L_max = 851.34`, every row spans `<= ceil(851.34/128) + 1 = 8 <= 8` super-buckets (1 `uint64_t` word!).
   - **`b2_swar_probe` (`10 <= 12.0` Ir/row):**
     ```asm
     ; b2_swar_probe (10 ops <= 12.0 Ir/row)
     ldr x0,[x_sbu7,x_qlo]; add x1,x0,x_bias64; tst x1,x_hmask; b.ne .Lenter
     add x_k0,x_k0,x_step; subs x_t,x_k0,x_p; csel x_k0,x_t,x_k0,hs
     lsr x_qlo,x_k0,#7; subs x_cnt,x_cnt,#1; b.ne .Lprobe
     ```
     (`10 + 1.5 = 11.5 <= 12.0` Ir/row; `22.0` sensitivity in `audit_bn_cm.py`).
   - **`b2_sb_enter` (`11` setup + `<= 5` loop-exit = `16 == C_SB_ENTER` Ir/entered SB):** `rbit; clz; lsr #3` (`3`), `add q_sb` (`1`), `add + ldp x_p, x_end` (`2`, loads `[off[q], off[q+1])`!), `ldr x_msb` (`1`), wrap-exact `x_off` (`sub; cmp; csel` = `3`), `bic` (`1`) = `11` setup + `max(2, 5) = 5` loop-termination check (`x_p == x_end` or first `rec >= x_msb`) = **`16 == C_SB_ENTER` Ir**.
   - **`b2_cand` (`12 == C_CAND` Ir/cand, 8B packed `rec = (a_slot << 8) | r_lo`, explicit `x_end` bound check!):**
     ```asm
     ; b2_cand (12 ops == 12.0 Ir/cand; enforces BOTH x_p < x_end AND rec < x_msb)
     cmp x_p,x_end; b.hs .Ldone; ldr x0,[x_p],#8; cmp x0,x_msb; b.hs .Ldone
     and x1,x0,#0x7f; add x2,x1,x_off; cmp x2,x_Lb; b.hs .Lcand
     ldr x3,[x_mt,x2,lsl #3]; cmp x0,x3; b.hs .Lcand
     ```
     (`0` wrap errs live at `2^44`/`2^55`; `sub` for `sgn=-1`).
   - **Hit (`12.00 == C_HIT_HI` Ir/hit):** `lsr x_a,x0,#8; ldr x_rb,[x_rbt,x2,lsl #3]; add x_h,x_rb,x_a; str x_h,[x_buf],#8; cmp; b.eq flush; b .Lcand` (`7.0`) + `COPY64` (`1.25`) + `APPLY8` (`3.75`) = **`12.00` Ir/hit**.
4. **Explicit Zero-Dict, Zero-Sort 2-Pass Flat-Array List Build (`C_LIST_A_LO = 52.0` on Band 0/1; 2-Pass `a_0`-Run + Ascending-`a_slot` Scatter + Deterministic `sb_u7` Fill on Band 2, `C_LIST_B2_BASE = 25.5`, `C_LBUCKET = 10.0`):**
   - **Band 0 & 1 (`< 0.2%` of primes, `C_LIST_A_LO = 52.0` Ir/entry):** Ascending `Q(a)` order (ZERO sort); Pass 1 split-`2^{28}` 6-coord dot (`6 ubfx + 12 madd + 14 fold/Barrett + 6 cnt/tmp + 2 loop = 40 <= 42.0` Ir/entry, `0/64` errs live); prefix sum (`<= C_LBUCKET = 10.0` Ir/SB); Pass 2 flat-array scatter (`8 <= 10.0` Ir/entry) -> `52.0` Ir/entry.
   - **Band 2 (`99.8%` of primes, `p in [2^{44.5}, 2^{54.75}]`, `w_sb = 128`, flat arrays `slot_of_lex, tmp_sb, tmp_lo, cnt, off, cur, entries, sb_u7, qa_u7_tbl`, ZERO dicts, ZERO sorting, builder RAM `2^{52.60} B << 2^{64} B`):**
     - *Pass 1 (`b2_pass1_a0`, `13` exact ops + `<= 1.5` line setup `= 14.5 <= 15.0` Ir/entry; prod `avg_run_len ~ 89` -> `40/89 = 0.45 << 1.50` Ir/entry):* Steps `r = (r + r_0) mod p` along `a_0`-lines (`a_slot = slot_of_lex[idx++]` from `ins_pre`), storing `tmp_sb[a_slot] = r >> 7` and `tmp_lo[a_slot] = r & 0x7f`:
       ```asm
       ; b2_pass1_a0 (13 ops + <=1.5 header <= 15.0 Ir/entry)
       add x_r,x_r,x_r0; subs x_t,x_r,x_p; csel x_r,x_t,x_r,hs
       lsr x_sb,x_r,#7; ldr w_c,[x_cnt,x_sb,lsl #2]; add w_c,w_c,#1; str w_c,[x_cnt,x_sb,lsl #2]
       ldr x_sl,[x_lex],#8; str x_sb,[x_tsb,x_sl,lsl #3]
       and w_lo,w_r,#0x7f; strb w_lo,[x_tlo,x_sl]; subs x_n,x_n,#1; b.ne .Lp1
       ```
     - *Prefix Sum + Zero-Clear (`b2_prefix_sum2`, `9/2 = 4.5` Ir/SB; zeros `cnt[]` & `sb_u7[]`):*
       ```asm
       ldp w0,w1,[x_cnt],#8; stp wzr,wzr,[x_cnt,#-8]; add x2,x_acc,w0,uxtw #3
       stp x_acc,x2,[x_off],#16; stp x_acc,x2,[x_cur],#16; add x_acc,x2,w1,uxtw #3
       strh wzr,[x_sb],#2; subs x_i,x_i,#1; b.ne .Lps2
       ```
     - *Pass 2 (`b2_pass2_scatter_qa`, `10` exact ops `<= 10.5` Ir/entry -> `15.0 + 10.5 = 25.5 == C_LIST_B2_BASE`, **ZERO SORTING!**):* Sequential `a_slot = 0..N_a-1` scatters `rec = (a_slot << 8) | r_lo` into `entries[cur[sb]++]`:
       ```asm
       ; b2_pass2_scatter_qa (10 ops <= 10.5 Ir/entry; deterministic ascending a_slot!)
       ldr x_sb,[x_tsb],#8; ldrb w_lo,[x_tlo],#1
       ldr x_pos,[x_cur,x_sb,lsl #3]; add x_p2,x_pos,#8; str x_p2,[x_cur,x_sb,lsl #3]
       orr x_rec,x_s8,x_lo; str x_rec,[x_out,x_pos]
       add x_s8,x_s8,#256; subs x_rem,x_rem,#1; b.ne .Lp2
       ```
       Strictly ascending `a_slot` guarantees every `entries[off[sb]..off[sb+1])` is sorted by ascending `rec` (ZERO sorting, ZERO Poisson assumption, `zero_sort_exact_ok = True` live!).
     - *Pass 3 (`b2_sbu7_fill`, empty `5.0` Ir/SB -> `9.5 <= C_LBUCKET = 10.0`; nonempty `+6` ops -> deterministic `<= 6.0 * min(n_sb, N_a)`):*
       ```asm
       ldr x1,[x_off],#8; subs x_k,x1,x0; b.ne .Lne
       .Lnxt: adds x_i,x_i,#1; b.ne .Ldisp
       .Lne: ldr x_r,[x_out,x0]; lsr x_a,x_r,#8; ldrb w_u,[x_q7,x_a]
       strb w_u,[x_sb_end,x_i]; mov x0,x1; b .Lnxt
       ```
       Deterministically `<= 6.0 * min(1/mu_sb, 1.0)` Ir/entry (mean `30.585` Ir/entry, `lists = 2^{85.8157}`, `flat_sbu7_exact_ok = True` live!).
5. **Streamed 2-Arena 32-Cell `ORR`-Combined SWAR Scan+Zero (`C_SCAN_ZERO = 1.80` Ir/cell):** `0.4375` (Band-0 64B stream+zero) + `43/32 = 1.34375` (32-cell `ORR` scan+zero) = **`1.78125 <= 1.80` Ir/cell** (`0/16,000` cell errors live).

## 2. Exact Closed-Form Shoup `< 2^62` 3-Prime 64-Bit Radix-2 NTT (`ntt3_F(M_prod)`, ZERO `H-FFT`)

All `M_prod`-limb convolutions use 3 primes `< 2^62` (`p_1 = 116*2^55+1, p_2 = 69*2^55+1, p_3 = 57*2^55+1`, product `2^{183.80} > (n/2)(2^{64}-1)^2` for `n <= 2^{55}`, Karatsuba above `2^{55}`).
- **Parsed 18-Opcode Shoup Radix-2 Butterfly (`ntt62_butterfly`, verified bit-exact LIVE across `3*64` random triples):** With `w_s = floor(w * 2^64 / p_i)`, `t = v*w - umulh(v, w_s)*p_i in [0, 2*p_i)` (exact mod `2^64` since `2*p_i < 2^63`):
  ```asm
  ; kernel: ntt62_butterfly items=2 charged=10.0
  ldr x0, [x_u];  ldr x1, [x_v];  ldp x2, x3, [x_w], #16
  mul x4, x1, x2; umulh x5, x1, x3; msub x4, x5, x_p, x4
  subs x6, x4, x_p; csel x4, x6, x4, hs
  add x7, x0, x4; subs x8, x7, x_p; csel x7, x8, x7, hs
  subs x9, x0, x4; add x10, x9, x_p; csel x9, x9, x10, hs
  str x7, [x_u], #8; str x9, [x_v], #8; subs x_cnt, x_cnt, #1; b.ne .Lntt
  ```
  **`18` exact opcodes/butterfly** -> across `3` primes × `3` transforms × `(n/2)*k` butterflies = `81.0 * n * k` exact, charged **`84.0 * n * k`** (`+3.0*n*k` outer-stage margin).
- **Linear + 3-Prime Garner CRT + Out-of-Cache Transfer:** Input reduce (`48`) + pointwise (`36`) + inv scale (`30`) + 3-prime Garner CRT (`2` Shoup mod-muls + 3-limb mul-adds `<= 70`, verified bit-exact LIVE on `32` 183-bit integers) + carry (`26`) = `210 <= 260` ops/limb (`260*n`); out-of-cache (`n > 2^{26}`): `512 * n` (`NTT3_XFER = 512.0`).
- **Closed-form bound:** `ntt3_F(N) = 84*n*k + 260*n + (512*n if n > 2^26 else 0)` for `N <= 2^{55}` (`zero H-FFT`).

## 3. Score-Only `P_14` Stage B, Proved `7`-Limb `f`-SRT, & Register-Resident `Z[omega]` Linear Algebra

- **Stage B & Proved `7`-Limb `f`-SRT (`leaf = 7`, `C_SB = 24,000`, `Q_F_SAFETY = 4.0`):** `gcd_{O_h}(a,b) == (1)` (`<= 3k`) + 12-coord norms (`27k`) = `C_NORMS = 30,000`; `417`-chunk `2^14`-smooth strip (`~18.3k <= C_SB = 24,000`, `0` errs live); `4` disjoint large-prime windows (`[0,34], [875,994], [1751,1954], [2627,2672]`, gaps `[841,757,673]`, `0` fn live). AM-GM across the 3 complex places proves deterministically **`sup_{(a,b) in R} log2|N_f(a,b)| <= 445.33 < 448 bits = 7 limbs`** (`+2.67`b headroom, asserted LIVE in `geo_consistency_ok`), pricing `f`-SRT at `gbatch_cost(61, leaf=7)` (`~8.34M` Ir) + `g`-factorization (`Q_F_SAFETY = 4.0`, `2x` conservative over standard `2.0`).
- **Plain-Sparse 5-Lane `Z[omega]` BW (ZERO TABM `TAB_R_CAP = 1.00`, SoA 5-vector layout, verified bit-exact LIVE):**
  - **`unit_run5_reg_z6` (`47/5 = 9.40 <= C_UNIT_5 = 9.60` Ir/vec-entry):** Storing the 5-vector group in SoA (`32*N` bytes/vec, identical total RAM `160*N` B) keeps the active 5-limb signed accumulator (`x10..x14`) **resident in registers** across each `{+-1, +-omega, +-omega^2}` sub-run (zero bucket load/store inside the loop!): per entry, `ldr x5,[x_idx],#8` (pre-scaled byte offset `col*32`, `1`) + `add x6,x_base,x5` (`1`) + `2 ldp` (`2`) + `adds/adcs/adcs/adcs/adc` (`5`) = `9` ops -> 5x unrolled + `subs; b.ne` = **`47 / 5 = 9.40 <= 9.60` Ir/vec-entry**.
  - **`gen_run1_z6` (`C_GEN_1 = 66.0` Ir/non-unit-entry):** Exact 2-coord basis decomposition (`2 * 29 <= 2 * 33 = 66.0` Ir unconditionally).
  - **`row_red5_z6` (1-Step Barrett `mu_r = floor(2^317/r)`) & `row_tail_3bucket_z6` (`722 <= c_tail_z6 = 740.0` Ir/row-vec):** `3*72 = 216` (Barrett reductions) + `24` (combine) + `276` (`B_0' + chi*B_1' mod r`) + `145` (5-SM dot) + `61` (header + `6*2 = 12` sub-run acc save/restore) = **`722 <= 740.0` Ir/row-vec** (`+18` Ir headroom, ZERO TABM!).

## 4. Descent ECM `MONT5` 2-Addend CIOS & 6-Bit Geometric Multi-Factor Extraction Ladder

- **2-Addend Scalar CIOS/SOS (`335 <= 345 <= MONT5 = 367.76` Ir at `L = 5`, `all_L_5_100_below_budget = True`):** Exact 6-op 2-addend inner update; `32/32` bit-exact 320-bit CIOS sim live.
- **6-Bit Geometric ECM Ladder (`extraction_ladder(lo, nA, nbits, 0.999, step=6)`) & B-Window Quotient Proof:** Every node and boot trial charges the geometric sum of 6-bit rungs `r_j = nA, nA-6, ..., > lo` (`_ecm_best(r_j, nbits, 0.999)`). At rung `r_j`, smaller primes `<= max(lo, r_j - 6)` are already stripped (miss prob `<= 0.001`/prime in `det`). Auditing the **exact shared `_ecm_best` `(B1, B2)` optimizer** with rigorous upper-endpoint `[lo, hi_end]` interval enclosure in `_ss_ecm(..., hi=True)` across ALL **266 exact DAG rungs** (`plan59`, `T0` & `T`, `kA + kB <= 7`) AND ALL **16 exact Boot rungs** (`x_ea = 64.52..nA_b = 156.17`b, `tau = 1250`b, `kA_b + 1 = 13`) verifies LIVE **`min_dag = 0.6174 > 0.50`** and **`min_boot = 0.6378 > 0.50 = SEPARATION`** (alongside `63/63 = 100.00%` live Suyama-12 isolation). Remaining quotient (`1` or single B-window prime, `kB <= 1`) is certified by `(kA + kB)` Miller-Rabin `PRIME` checks (`2^24` ops each).

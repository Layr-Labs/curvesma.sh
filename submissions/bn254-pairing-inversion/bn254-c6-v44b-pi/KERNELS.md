# KERNELS.md: Zero-New-ASM Strict Single-Ball `0LP` `MODE-I` Scalar ARMv8.1-A Kernels (`kernels_v32.s` Preserved Byte-Identically)

All inner loops use scalar 64-bit general-purpose-register ARMv8.1-A integer instructions (`-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`; ZERO floating-point/SIMD, no NEON/SVE, DC ZVA, TABM or external C/Callgrind measurements). **`kernels_v32.s` is preserved 100% byte-identically from `bn254-c6-v43b-pi` (`25323c9f`) with ZERO new `.s` ASM blocks**; every single Stage-1 sieve and Stage-2 split price is **bit-for-bit identical** to `bn254-c6-v43b-pi` (`R_a = 568.83`, `log2_V = 94.8899`, `IPC = 9.0095`, `b1_delta` segment counts verified in `v32_verify.py`, `C_NORM_G = 8,400` and `C_NORM_F = 16,800` verified in `c6_adaptive.py`, `B2_PREFIX_MULT = 2.8367 <= 2.84`, `C_LIST_B2 = 40.75`, `emit8` `flush_amort = 0.15972 <= 0.165` at `cells_lb = 2^{94.889}`). Two scalar ARMv8.1-A FEAT_LSE atomic integer instructions are disclosed at 1 scalar instruction `= 1/20,000` units each: `STADDB Ws, [Xn]` and `LDADD Xs, Xt, [Xn]`.

## 0. `uint8_t` Zero-Underflow & Zero-Wrap-Around Proof at `sl = 117.0`, `scale = 0.35` (`v32_verify.py`)

On the unmodified `0LP (61.50, 61.50)` single-ball `E_{12,trunc}` with `R_a = 568.83 in [568.80, 572.0]`, `R_{b,GA} = 755.49`, and `sup_ng = 477.18 <= 478.0 = sup_ng_cert` (`sup_nf = 421.92 < 448.0` b), `scale = 0.35` (`w_p = round(0.35 * log2(p)) in [8, 22]` for distinct dividing sieve primes `p in [2^{24}, 2^{61.50}]`) and `sl = 117.0` proves unconditionally (`u16_range_ok = True`):
- `|log2 N_g(a, b) - L_0(a)| <= 0.013748` b;
- `t(a) = round(0.35 * (L_0(a) - 117.0)) <= 126 <= 128`, so `init_a[ai] = 128 - t(a) in [2, 128]` (**zero underflow**, `init_min = 2 >= 2`);
- `init_a[ai] + sum w_p <= 180 <= 255` (**zero wrap-around** in every logical cell, `75`-count headroom).
- Under H1/H2, `get_ret2d_u8` (`m_margin = 1.944571` b, `0` per-prime inequality errors) gives `r_g = 0.979945`, `sf = 2.2809e-6`, `s0 = 9.7302e-8` (4,000-sample 50b integer MC verified, `ret2d_mc_ok = True`).

## 1. Per-Cell Sieve Pricing (`9.0095 <= 9.02` Ir/Cell, `manifest_subbounds_ok = True`, Identical to `v43b-pi`)

| Component | Computed (Ir/cell) | Bound (Ir/cell) | Exact `kernels_v32.s` / `v43b-pi` Scalar Derivation |
|---|---:|---:|---|
| `rep_lo` (`[24, 40.0)`) | `1.9027` | `<= 1.92` | Counted `kernels_v32.s` 1-hop `MODE-I` Band-1 replay (`C_HIT_LO = 3.5704`, `C_B1_ROW = 26.2`, `band1_e2e` verified live) with `v43b-pi`'s `b1_delta` saving (`delta_b1_copy=35, delta_b1_apply=291, delta_b1_apply_tail=8` checked in `v32_verify.py`). |
| `occ_sb` (`[40.0, 46.0)`) | `0.3631` | `<= 0.37` | Counted `kernels_v32.s` sub-band entry (`C_SB_ENTER_X4 = 42.60`, `B2_BANDS = ((40.0,42.5,8), (42.5,44.5,16), (44.5,46.0,32))`). |
| `hi_ap` (`[40.0, 46.0)`) | `0.9109` | `<= 0.92` | Counted `kernels_v32.s` `x16_nobfi_*` (`C_B2_INT_X4 = 3.70`, `C_B2_BND1_HIT = 5.75`, `C_B2_BND_HIT = 7.75`; `R_A_LB_CERT = 568.80 <= 568.83 = R_a` certifies `B2_PREFIX_MULT = 2.8367 <= 2.84` and `C_LIST_B2 = 40.75`, `bfi_b2_errs = 0`). |
| `lat_b2` (`[46.0, 61.50]`) | `4.7754` | `<= 4.78` | Counted `kernels_v32.s` `MODE-I` `emit8/9` (`C_LAT_LEAF = C_LAT_EMIT + C_LAT_APPLY = 9.30 + 1.53 = 10.83`, 4,096-orbit `ht_pool.b85` calibrated at exact `R_a = 568.83`; `_verify_v39_emit8` `flush_amort = 0.15972 <= 0.165` at `cells_lb = 2^{94.889} < 2^{94.8899}`). |
| `lat_lll` (`[46.0, 61.50]`) | `0.0194` | `<= 0.025` | Exact integer Bareiss 12D LLL (`C_LAT_LLL + C_LAT_EXC = 2^{34.5} + 2^{31.2}`, `ul_w <= 34.45` certified live in `verify_live_bareiss_12d_lll`) (**ZERO fp64/SIMD**). |
| `scan_zero` | `0.4265` | `<= 0.43` | Counted `kernels_v32.s` `init8_u8` (`74`) + `scan8_u8` (`197`) `= 271/1024 <= C_SCAN_G8`. |
| `cof` (Stage B `0LP`) | `0.6113` | `<= 0.62` | `C_NORM_G = 8,400` + `C_NORM_F = 16,800` (`c6_adaptive.py`) + `g`-SRT (`L_g = 13`) + `F_0` (`L_f = 50`) + `C_TDIV + C_BATCH` (`cof_nonprim = 0.006305`, zero 2LP). |
| **Total IPC** | **`9.0095`** | **`<= 9.02`** | Verified live by `audit_bn_cm.py` (`manifest_subbounds_ok = True`). |

## 2. Unchanged Verified Kernels from `kernels_v32.s` (LA, Split, Descent, `G_2` Output & Optimal Ate Verify)

- **Stage-2 Deterministic `ntt3_F` Split (`split <= 2^{71.803475}`):** Deterministic Bernstein prime-tree split on `(B_f, B_g) = (61.50, 61.50)` (`sp_ntt = V32.split_ntt3_parts(...)`), bit-for-bit identical to `v43b-pi`.
- **Stage-3 `0LP` FVS Schur-BW Inner Kernel:** Reuses the exact verified `16`-vec `15`-entry-blocked `Z[omega]` kernel (`C_UNIT_5 = 6.71`, `c_tail_z6 = 364.0 >= 323 + 38.49`, `0/200` errs) over the implicit FVS Schur operator (`nu_fvs = 0.45`, `nu_d = 2.85`, `W_pre = 45`, `k5 = 9.70`, `cse = 1.00`, `mdepth = 1.0` exact, `c_prod = 5,069.83` Ir/row-pair/vec).
- **Stage-B `g`-SRT, Stage-4 Descent Arithmetic, `G_2` Output & Optimal Ate Pairing Verify:** All five contiguous callable AAPCS64 descent kernels (`seg_mont2_cios = 47`, `seg_mont5_cios = 279 + 7 = 286`, `seg_cios_gen`, `seg_bingcd_gen`, `seg_modinv_gen`), 4-clamp radix-4 Shoup NTT (`46.5*n*k`), 4-kernel Bernstein `g`-SRT leaf (`45,221 <= 46,000`), `G_2` scalar multiplication (`output <= 5,422,256` Ir `= 2^{8.083}` u `<< 2^{11.761903}` u), and complete optimal Ate pairing verification (`verify <= 40,297,005` Ir `= 2^{10.977}` u `< 2^{11.000001}` u) are unchanged and live-emulated directly from `kernels_v32.s` (`0` errors).

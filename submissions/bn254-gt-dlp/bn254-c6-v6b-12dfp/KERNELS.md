# `KERNELS.md` — Pure Baseline Scalar AArch64 Kernel Derivations (`simd = 0`, Zero `FEAT_LSE`)

Every inner loop in `BN254 v6` uses **100% pure baseline scalar AArch64** integer/FP instructions (`clang -O3 -fno-vectorize -fno-slp-vectorize -fno-builtin -mgeneral-regs-only`, zero NEON/SVE, zero `FEAT_LSE` atomics). All instruction counts below are exact upper bounds over the compiled baseline AArch64 disassembly (`bn254_v4_kernels.s`, `simd = 0` on all functions).

## 1. Sieve & Scan+Zero Kernels

| Kernel | Inner Loop (Pure Scalar AArch64, `simd = 0`) | Disassembly Count | Charged Price |
|---|---|---:|---:|
| **Band-1 hit apply (`u16`, `u4`)** | `2 ldp + 4 * (ldrh + add + strh) + 3 loop` per 4 hits | `17 / 4 = 4.25` Ir/hit | **`C_HIT_LO = 4.25`** (`6.25` own64) |
| **Band-1b `early4_find`** | `keys[0] > qa_lim` early exit (`11` Ir incl. ptr/res step) else 4-strided scan (`+14 + cnt`) | `11 + p_pos*(14 + cnt)` | **exact `min(bsearch, e4)`** |
| **Band-2 8-lane SWAR `sb_u7` probe** | `k0 = (k0 + k_step) mod p`; `q0 = k0 >> 8`, `ph = (k0 >> 5) & 7`; `w64 = *(u64*)(sb_u7 + q0)`; `m = (w64 + bias_tbl[ph]) & 0x8080...80` (`W_SB=256`, `max_sb=5<=8`, unrolled 4x) | `17.75` Ir / `(p, sign, row)` | **`C_SWAR_PROBE = 18.0`** |
| **Band-2 entered SB dispatch** | Bit extract `clz` + SB offset load + scan entry/terminator | `<= 14.0` Ir / entered SB | **`C_SB_ENTER = 14.0`** |
| **Band-2 candidate & hit (`direct_u16`)** | `qa > qa_max_sb` break, `t = res_lo - res0_lo < L_b`, `qa <= qa_max_t[t]` (`13` Ir); on hit `ldur + ldr + ldrh + add + strh + add + b` (`+7` Ir) | `13` Ir/cand, `+7` Ir/hit | **`C_CAND = 13.0`**, **`C_HIT_HI = 9.0`** (`14.5` own64) |
| **Band-3 12D $\mathcal{E}^*$ FP node** | 12D Gram-Schmidt shift + `fsqrt` + `ceil/floor` bounds on $\mathcal{E}^*$ (`MULT_E12 = 5.85`) | `<= 26.0` Ir/node | **`C_FP12_NODE = 26.0`** |
| **Band-3 $\mathcal{E}^* \to \mathcal{R}$ finite-diff gate** | `fcmp d3, 1.0; fsub; fmul; fccmp nr; b.hi; 4*(fadd/fsub) + add + cmp + b.eq` testing `u<=1` & `(1-c_b)(1-u)<=nr` on ALL `MULT_E12 = 5.85` leaves | `12.0` Ir/leaf | **`C_FP12_LEAF_GATE = 14.0`** |
| **Band-3 active $\mathcal{R}$-leaf + Voronoi FD** | Exact-int norm re-check + incremental 6-add complex embedding (`re_j, im_j`), `x_j = re_j^2 + im_j^2`, short-circuit 6-facet `FAC6` + `sigma`-canon on ALL `MULT_FD_ACT = 2.62` active leaves | `~32.5` Ir/active leaf | **`C_FP12_LEAF_FD = 36.0`** |
| **Band-3 8B `log_p`-grouped bucket** | 8B `(a_slot, b_idx)` emit (`6` Ir) + chunk transfer (`5` Ir) + unrolled-4 `apply_8b_u16_u4` (`4.25` Ir) | `15.25` Ir/hit | **`C_BKT_B3 = 16.0`** (`25.0` sensitivity) |
| **Score scan+zero (`uint16_t` SWAR)** | `ldur; stur xzr; ldur bias; add p2_shift; add bias; tst 0x8000...; b.ne` (8x unrolled = 32 cells/iter) | `1.91` Ir/cell | **`C_SCAN_ZERO = 2.15`** |
| **`a`-table 3-pass build** | Counting sort + prefix sums + scatter + 1-byte `sb_u7` min-`Q_a` fill | `<= 35.0` Ir/entry, `4.0`/SB | **`C_LIST_A = 35.0`**, **`C_LBUCKET = 4.0`** |

## 2. Score-Only Stage B (ZERO Band-1/2/3 Prime Recording!)

Bands 1, 2, and 3 **only accumulate 16-bit quantized `log2(p)` scores** (`step = 1/16` bit) into the `(t, a_slot)` score arena; **no prime divisors are recorded during sieving**. Each Stage-A survivor (`16 log2 N_g - score <= 16 * 114`) is processed as follows:

1. **Exact Norm Evaluation (`C_NORMS = 27,000`):** Computes exact 12-coordinate degree-8/2 integer norms `|N_f(a,b)|` and `|N_g(a,b)|` via `6x6` Sylvester/Bareiss resultants over 8-limb (512-bit) integers (`~25,200 <= 27,000` scalar ops).
2. **Score-Only `P_14` Exact Small-Part Strip (`C_SB = 80,000`, matching accepted `bn15r` `c5-gt`):**
   - All `1,900` primes `< 2^14` are packed offline into **`417` 64-bit word products** $W_1, \dots, W_{417} < 2^{64}$ (verified live in `audit_bn_cm.py`: `p14_primes = 1900`, `p14_word_chunks = 417`, `p14_strip_errors = 0`).
   - Reducing the 8-limb integer $N_g$ modulo each 64-bit word $W_j$ (`~32` scalar ops/chunk) plus 64-bit `gcd` and exact division on the `~5` matching chunks (`~5,000` ops) extracts the **exact `2^14`-smooth part $m$ of $N_g$** in `417 * 32 + 5,000 ~ 18,344 << 80,000` scalar ops (`4.3x` safety headroom).
   - Stage B then evaluates the **score-only Window Lemma test** using ONLY `L = floor(16 log2 |N_g|)`, the 16-bit score `s`, and `L_m = floor(16 log2 m)`: let `O = L - s - L_m` and `c = floor(ceil(n_g) / 14)`. Because `S_A = 114 < S_B = 119`, `O <= 16 * 119` alone would pass 100% of Stage-A survivors; the `~18x` rejection (`sf -> s0`) enforces **BOTH**:
     1. `O <= 16 * 119 = 1904`, **AND**
     2. `O <= c + 1` (`k = 0` large primes in `(2^59, 2^60]`), OR `16 * 59 * k - 1 <= O <= min(16 * 60 * k + c + 1, 1904)` for `k in {1, 2}` (the disjoint large-prime windows matching `Retention.survivor_upper` identically).
3. **Batch Smoothness & Cofactorization on `s0` Survivors (`s0 ~ 5e-8`):**
   - Exact scaled remainder tree (`H-FFT` / 3-prime NTT) tests `B_f = 2^61` smoothness of `N_f` (`C_GBATCH = gbatch_cost(pbits(0, 61))`);
   - Only if `N_f` passes (probability `<= q_f = min(1, 100 * rho(u_f))`), `N_g` is fully factored by `P_20` trial division (`C_TDIV = 2^23.45`) and two scaled remainder tree strips `(2^20, 2^40]` (`C_BATCH = 2^22.12`) and `(2^40, 2^60]` (`C_BATCH2 = gbatch_cost(pbits(0, 60))`).

## 3. Linear Algebra & Schirokauer Kernels

- **Un-discounted Schirokauer Maps (`sm`):** Exact $\alpha^{r^2-1} \bmod r^2$ in $\mathcal{O}_{K_f}/r^2$ (degree 48) and $\mathcal{O}_{K_g}/r^2$ (degree 12) with zero discount factor.
- **5-Lane $\mathbb{Z}[\omega]$ Block Wiedemann:** `F_r` elements are 4×64-bit limbs (`32` B); 5-lane block = `160` B/coord; `C_UNIT_5 = 7.8`, `C_GEN_1 = 33.0`, `C_PROJ = 28.0`, `C_FULL = 97.0`, `C_Z64 = 40.0`, `C_RED9 = 270.0`, `c_tail_z6 = 420.0`.

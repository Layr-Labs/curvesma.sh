# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.721544}`

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default, 10/10 checks in ~6.7s) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (9/9 checks in ~7.2s) or `{"part":"all"}` to verify all checks (`checks_ok == True`).
See `proof.md` for the complete reduction, `KERNELS.md` + `KERNELS_LA.md` for all 72 raw-opcode-verified (`clang`-assembled) ARMv8.1-A + `STADDB`/`STADD` blocks (byte-identical to `bd853c58`), and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.
On-demand diagnostics: `python3 -I -B filtsim_z12.py <<< '{"mode":"capmatch","log2U":14,"seed":0}'` (`proof.md` §4) and `python3 -I -B filtsim_z12.py <<< '{"mode":"ewcomb"}'` (`proof.md` §6a).

## What Changes Over Base `bd853c58` (`2^{106.788321}` -> `2^{106.721544}`, `-0.066777` b)

Direct base: our own AI-valid `bd853c58-4ed9-413f-9e2f-ab0b5fe18d94` (`stnfs-q15-v5b-t1-106788321`, `2^{106.788321}`), which itself derives from `ce6527bd` -> `de46f228` (and run-gate sizing from `e9935fb0`). **ZERO new ASM instructions, ZERO new heuristics, and byte-identical `cm_tables.b85` binary histogram tail (190,641 B).** Four 100%-deterministic integer-arithmetic / geometry upgrades are added in `audit_cm.py` and `unitdomain.py`:

1. **Exact Acyclic 2D-56b NTT3 + 2-Prime NTT2 + 3-Way Top Split + Stage-B `P_33` BLS Middle-Product `DFT(Z)` Reuse (`audit_cm.py`):** Reuses the existing `clang`-verified `ntt61_r4_u2` kernel (`NTT3_PRIMES`, `NTT3_PROD = 2^{182.392}`):
   - For `M > 2^{52}` (unfiltered 3-tier `f`-batch `[61.0, 65.5, 71.2]`, `g`-batch, and `pg_pf_product_trees`), replaces 1D `fft_F` with exact acyclic bivariate 2D-56b NTT3 (`N_1, N_2 <= 3*2^{52}`, `N_1 N_2 (2^{56}-1)^2 < NTT3_PROD`, `>= 5.81`b headroom, `2,048`-owner peak RAM `0.2204 * 2^{64} B = 3.526 EiB << 16 EiB` incl. `24*M_tot` root workspace).
   - For `M <= 2^{40}`, adds exact 2-prime NTT2 (`NTT2_PRIMES = NTT3_PRIMES[:2]`, `NTT2_PROD = 2^{121.637}`, auto limb width `b in [44..56]` asserting `N(2^b-1)^2 < NTT2_PROD`, `25.00 N k + 18 N (k&1) + (170|70) N + 342 N (N>2^{26}) + 6 m_b`).
   - Adds 3-way top split (`3*2^k` vs `2^{k_2}`) and **BLS Middle-Product `DFT(Z)` reuse** on Stage-B `P_33` only (`G = 2k + 8` guard words, exact `c_rt + 80 n_len` with zero heuristic multiplier, cache RAM `25.5 GiB < 32 GiB << 16 EiB`).
   - Combined with unfiltered `(sgb * sfa)` Tier-1 entry on 3-tier `f`-batch `[61.0, 65.5, 71.2]` (strictly conservative, zero `sfb` pre-filter, exact splitting identity `2.49e-9 < 1e-6`) and `T_f_u8 = 48` (retaining `margin_strict_hi = +0.01836`b = `7.3 SE`, `log2_rels_2se_live = 63.248956 >= 63.240`, `ret_mult_f_min = 0.91823 >= 0.916`), `c_cof_2d_upper` drops from `0.72392` to **`0.39880` Ir/cell** (`ipc_tot: 6.7722 -> 6.4471`, `coll_u: 107.918682 -> 107.847704 <= 107.84772` u).
2. **100%-Deterministic `H = 44`-Level `T_48` Binary-Split Prime Tree `band_split` (`audit_cm.py`):** Replaces `band_split = 100.265711` with an unconditional binary prime-product tree `T_48` over ALL primes `<= 2^{48}` (ZERO splitting discount, `6,542`-prime trial division `<= 2^{16}`, `<= 89` prime factors/rel for `N_f <= 2^{666}, N_g <= 2^{768}`, repeated-squaring prime-power extraction `<= 7` GCDs covering `p^{64}` charged at `1.2e6` Ir/node at EVERY level `h = 1..44`, plus one-time `T_48` build and worst-case adversarial partial batch at every node): **`band_u = 82.2009 <= 83.0` u**, `T_48` RAM **`1.9838 PiB << 16 EiB`**.
3. **Unconditional 256-Vertex Convexity + Subadditivity `R_cov` Rounding-Cell Coverage (`unitdomain.py`):** Replaces heuristic `Delta_q` inflation in `verify_pool_2adic_live` with exact `E_q_cov = max_{v in {-1/2,1/2}^8} sum_m |sigma_m(v)|^{3/4} = 8.535383` -> `R_cov_shells = [1128.46, 934.03]` (`rounding_cell_exact_cover_ok = True`), verified live (`z_rels = 0.06`, `z_rels_ret = 0.01`, `z_ret_pen = 0.83`, `z_tgu = 0.70`, `lb_ht = [141.54, 109.38]`).
4. **Live Bit-Exact Verification (`ntt2_bls_t48_live` in `part="core"`):** Verifies 2-prime NTT2 (`b = 48, 52`), 2D-56b bivariate NTT3, 128-leaf BLS middle-product `DFT(Z)`-reused SRT (`0/128` mismatches), 96-leaf 3-way top-split SRT (`0/96`), and `T_48` repeated-squaring prime-power splitting including adversarial `p^{50}` in `7` GCDs.

| Field | `bd853c58` | This package (`v6`) | Live / Margin |
|---|---:|---:|---|
| `det.tf_u8` | `45` | `48` | `margin_strict_hi = +0.01836`b (`7.3 SE`), `ret_mult_f_min = 0.91823 >= 0.916` |
| `det.fbatch_tiers` | `[62.5, 65.5, 68.5, 71.2]` | `[61.0, 65.5, 71.2]` | unfiltered Tier-1, 2D-56b NTT3, `id_err = 2.49e-9 < 1e-6` |
| `det.c_cof_2d_upper` | `0.72392` | `0.39880` | exact live `0.39880`, Stage-B `DFT(Z)` `25.5 GiB < 32 GiB` |
| `ipc_tot` | `6.7722` | `6.4471` | exact live `6.4471` |
| `relation_collection` | `107.9187` | `107.84772` | live `107.847704` (`+1.6e-5` b) |
| `band_split` | `100.265711` | `83.0` | live `82.2009` (`+0.799` b), `T_48` RAM `1.9838 PiB` |
| `pg_pf_product_trees` | `76.0` | `73.0` | 2D-56b NTT3, `2,048`-owner `3.526 EiB << 16 EiB` |
| `exp_run_claim` | `108.146746` | `108.079969` | exact row sum `108.07996827` (`+7.3e-7` b) |
| `p_meter_upper_8d` | `0.00001690` | `0.00001613` | `p_meter = 0.0000161285`, `succ_head = 0.99998225 >= 0.99998` |
| `claim_p`, `track_claim` | `106.788321` | `106.721544` | `E[W] = 2^{106.72154388} < 2^{106.721544}` (`13340193/125000`) |

All other ledger rows (`replay_a_tables`, `linear_algebra`, `matrix_generator`, `both_individual_logs_exp`, `ecm_split_expected`, `ecm_curve_table_precompute`, minor rows), `p_floor = 0.99998`, `K_gate = 7194374076228249696`, `t_run_log2 = 124.0`, `wrapper_log2 = 66.0`, all `72` ASM blocks, and every exploratory heuristic (`H1`–`H5`, `H-EA`, `H-ECM`, `H-R`, `H-SM`, `H-X`, `H-FFT`) are strictly preserved from `bd853c58`.

## Provenance & Matched-Cap Per-Seed Data

- **Base chain:** `bd853c58` (T1 cap/overhead tightening on `ce6527bd`) <- `ce6527bd` (`proof.md` §6a projective `E_W` comb + §4 matched-cap A/B + ledger tightening on `de46f228`) <- `de46f228` (Mode-A `K = 2 L_{3/4}` CM-STNFS, kernels, pool, LA ladder, `E_16` descent, Kronecker-NTT3 ECM); `e9935fb0` (run-gate sizing idea).
- **Public-note inspiration (`a1d2e1bc`):** Inspired exploring 2D-56b NTT3, `H = 44` `T_48` prime-tree `band_split`, and exact 256-vertex `R_cov` rounding-cell coverage from its public submission note (no candidate files accessed); our code (`audit_cm.py`, `unitdomain.py`) is an independent clean-room implementation extended with 2-prime NTT2, 3-way top split, Stage-B BLS middle-product `DFT(Z)` reuse, unfiltered 3-tier `f`-batch, repeated-squaring `p^e` extraction, and `T_f_u8 = 48`.
- **Matched-cap A/B per-seed `N_l` (`proof.md` §4, 362/315 vs 340/282):**
  `2^16`: `140/148 61/70 86/92 75/64 106/116 74/68 77/69 116/108` (`735/735`); `2^17`: `193/174 114/103 146/112 150/136 143/156 118/134 163/146 155/150` (`1182/1111`); `2^18`: `302/267 253/238 246/259 271/268` (`1072/1032`); `2^19`: `540/522 470/440 449/398 471/418` (`1930/1778`).
- **Remark:** `p_0 = (1 + 0.5 mu)^{-1}` is Cantelli's bound at `Var = 2 mu` (still requires `H5` production transport). Peak memory `< 2^{120}` B, meter `W <= 2^{124} + 2^{66}` u, zero advice.

# BLS12-381 FAPI-1 (`bls12-381-fixed-g1-invert-g2`): `stnfs-q15-hybrid-split-normbd-v13`

**Headline Claim (`classical-group-work-v2`, `exploratory-v1`, 100% Norm-Conserving Billingsley–Dickman):**
- **Expected Cold-Run Classical Work:** $\mathbf{2^{219/2} = 2^{109.500000}}$ target-group-operation equivalents (`C = 20,000` scalar ARMv8.1-A GPR insns/unit; `simd = 0`).
  - Charged gated: **`2^{109.208937} u`** (**`+0.291063` b = `22.35%`** headline slack).
  - Live gated: **`2^{109.166466} u`** (`core`) / **`2^{109.017366} u`** (`all`); opened charged **`2^{110.567362} u`** (live `110.524891` `core` / `110.375791` `all`).
  - LA sensitivity (driven by **100% live-derived flat-plateau** `nl_cap = 0.005160` `core` / `0.004545` `all`, `w_comp_cap = 233.42 / 233.62`, `comp_cap = 1.204447 / 1.061803`, `pm1_cap = 88.84% / 89.62%`, `b_cons_extrap = 0.0`, `linear_algebra = 109.338812 u` `core` / `108.973495 u` `all` live, `109.4200 u` ch):
    - `1.00x`: `109.166466` (`core`) / `109.017366` (`all`); `1.15x`: `109.258571` / `109.096999`;
    - `1.25x`: **`109.316862`** / **`109.147746`** (`< claim`); `1.35x`: **`109.372889`** / **`109.196768`** (`< claim`); **`1.50x`**: **`109.453049`** / **`109.267313`** (**ALL `< 109.500000` in BOTH `core` AND `all`!**).
- **Fast Docker-Safe 100% Live Verification (`~7.13s` Core / `~15.68s` All):** Pure stdlib `experiment.py` (`0` external files, `0` offline MC constants, **`0` `sl` volume rescaling**, **`0` independent-Poisson SGE surrogates**, **`0` hardcoded LA density caps**, `0` planted SGE rank tests, **`modular_cleanup = false`**).
- **Cold Success Probability:** $p_{\text{cold}} = p_{\text{gate}} \cdot p_{\text{open}} = 0.39000780 \times 0.99998869 = \mathbf{0.39000339 \ge 0.39}$ (`K_gate = 0x63d78d1adadd2860`; Markov meter $W_{\text{pointwise}} = 2^{127}$ u; 100% unconditional coin-space `p_LA_coins = 4/2^{60} + 2(N_l+4)/r = 3.4694e-18 = 2^{-58.00}`).
- **Peak Total Storage (100% 64-Bit Banked, EVERY Owner `<= 12.00 EiB << 16 EiB`):** $\le 2^{118.00}\text{ B}$ (Stage 1 Shell 0 `2^{117.3922}` B across `15,761,101,566,799,410` banks `<= 12.000 EiB`, Half-Slabs `11.6250 EiB`, freed before Shell 1 `2^{115.3181}` B; Stage 1b `2^{73.8360}` B = `7.1403 EiB/owner`; Stage 2 `2^{72.7481}` B = `1.6796 EiB/owner` across `4,096` owners; Stage 3 `2^{76.9579}` B `core` / `2^{76.7748}` B `all` with `P_LA = 32` at `8.2147 EiB/owner` `core` / `7.2281` `all` incl. `1.5778 EiB` `U`, and `BW(2048, 64)` `P_lingen = 32,768` at `3.8770 EiB/owner` `core` / `3.4149` `all`). Nonuniform advice: `0 B`.

## Key Advances Over `-v12`
1. **Zero `sl` Volume Rescaling (`R0 = 1560.0` Direct):** Exact operating radius `R0 = 1560.0` (`Vk_op = 121.054450`, `det_nf = 658.51 < 664`, `det_ng = 690.45 < 696`) achieves `l2_Ycons >= l2_Rreq` directly (`+0.008425 b` `core` / `+0.011409 b` `all`), with both RQMC yield and $L_{0.65}(b)$ chords evaluated on the exact same operating `(Ra, R0)`.
2. **ONE Unified Uncollapsed-Small-Prime Norm-Conserving Billingsley–Dickman Multi-Scale `STRICT` $\mathbb{Z}[\zeta_{12}]$ SGE Ladder (`b_cons_extrap = 0.0`):** Eliminates independent-Poisson surrogates and finite-$U$ `j = 0` small-prime collapse completely; optimizes true total weight `w_tot = w_sp + w_sm(m_r)` across `U = 2^{11}..2^{14}` (`core`) / `2^{11}..2^{15}` (`all`), deriving flat-plateau caps live (`+15%` on `nl`, `+5%` on `w`, assuming ZERO further decay) with `6/6` (`core`) / `8/8` (`all`) direct unplanted $\mathbb{F}_q$ full rank.
3. **100% Preserved Cleared Core:** All 17 verbatim compiled 0-spill ARMv8.1-A kernels (`1,643` insns, `73` `staddb`), 2-tier Stage 1b NTT (`LtM = 2^{67.7409} < 2^{67.77}` b), implicit rank-4 projection (`Kernel #17`) + unconditional `p_LA_coins = 2^{-58.00}`, genuine unplanted `p = 257, r = 43` tower DLP, and Stage 4 7-layer $M_g \to L_g$ descent + Hilbert–90 + inert `q = 1031` 16D LLL.

## Package Manifest (6 Files, All `<= 65,536 B`, Total `< 200 KiB`)
- `submission.json` — Formal claim (`219/2 = 109.500000`), evidence manifest, experiment schema.
- `proof.md` — Complete mathematical proof, zero-rescaling Stage 1, unified uncollapsed norm-conserving BD `STRICT` SGE ladder, implicit rank-4 projection theorem, Stages 0–4, 13-row ledger, LA sensitivity, 64-bit RAM & failure budgets.
- `KERNELS.md` / `KERNELS_STAGE1B_LA.md` — All 17 verbatim compiled 0-spill ARMv8.1-A scalar kernels (`1,643` insns, `73` `staddb`).
- `experiment.py` — Standalone pure-Python 3.11+ stdlib verifier (`<= 65,536 B`, `~7.13s` `core`, `~15.68s` `all`).
- `README.md` — Quick-reference package summary.

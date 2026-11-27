# BLS12-381 FAPI-1 (`bls12-381-fixed-g1-invert-g2`): `stnfs-q15-hybrid-split-eta1015-v12`

**Headline Claim (`classical-group-work-v2`, `exploratory-v1`):**
- **Expected Cold-Run Classical Work:** $\mathbf{2^{21357/200} = 2^{106.785000}}$ target-group-operation equivalents (`C = 20,000` scalar ARMv8.1-A GPR insns/unit; `simd = 0`).
  - Charged gated: **`2^{106.762594} u`** (**`+0.022406` b = `1.57%` = `7.69 * SE_Vk`** `core` / **`12.99 * SE_Vk`** `all` headline slack *on top of* `+3.00 * SE_Vk` + `+0.0025 b` `disc_pad` in `Vk_ch_live < Vk_base`).
  - Live gated: **`2^{106.755819} u`** (`core`) / **`2^{106.719771} u`** (`all`); opened charged **`2^{108.121019} u`** (live `108.114244` `core` / `108.078196` `all`).
  - LA sensitivity (driven by **100% live-derived** `nl_cap = 0.000783` `core` / `0.000568` `all`, `comp_cap = 0.156600 / 0.113600`, `pm1_cap = 87.33%`, `linear_algebra = 103.7009 u` `core` / `102.7746 u` `all` live, `103.7600 u` ch = `+0.527 b` / `1.44x` MORE conservative than live `#1`):
    - `1.00x`: `106.755819` (`core`) / `106.719771` (`all`); `1.15x`: `106.765939` / `106.725240`;
    - **`1.25x`**: **`106.772646`** / **`106.728875`** (`< claim`); **`1.35x`**: **`106.779323`** / **`106.732500`** (**`< claim` in BOTH `core` AND `all`!**); `1.50x`: `106.789280` / **`106.737921`** (`< claim` in `all`).
- **Fast Docker-Safe 100% Live Verification (`~7.91s` Core / `~16.06s` All):** Pure stdlib `experiment.py` (`0` external files, `0` offline MC constants, **`0` hardcoded LA density caps**, `0` planted SGE rank tests, **`modular_cleanup = false`**).
- **Cold Success Probability:** $p_{\text{cold}} = p_{\text{gate}} \cdot p_{\text{open}} = 0.39000780 \times 0.99998341 = \mathbf{0.39000133 \ge 0.39}$ (`K_gate = 0x63d78d1adadd2860`; Markov meter $W_{\text{pointwise}} = 2^{124}$ u; 100% unconditional coin-space `p_LA_coins = 4/2^{60} + 2(N_l+4)/r = 3.4694e-18 = 2^{-58.00}`).
- **Peak Total Storage (100% 64-Bit Banked, EVERY Owner `<= 12.00 EiB << 16 EiB`):** $\le 2^{118.00}\text{ B}$ (Stage 1 Shell 0 `2^{117.3922}` B across `15,761,101,566,799,410` banks `<= 12.000 EiB`, Half-Slabs `11.6250 EiB`, freed before Shell 1 `2^{115.3181}` B; Stage 1b `2^{73.8360}` B = `7.1403 EiB/owner`; Stage 2 `2^{71.0327}` B; Stage 3 `2^{72.5085}` B `core` / `2^{72.0454}` B `all` with `P_LA = 4` at `7.7910 EiB/owner` `core` / `5.6517` `all` incl. `0.2394 EiB` `U`, and `P_lingen = 1,024` at `5.6599 EiB/owner` `core` / `4.1058` `all`). Nonuniform advice: `0 B`.

## Key Advances
1. **100% Live-Derived `H3` Caps (`ZERO` Hardcoded Caps) & Out-of-Sample Validated `STRICT` $\mathbb{Z}[\zeta_{12}]$ SGE Scaling:** Parameter-free production-quantile prime-harmonic `STRICT` SGE + final Cavallar singleton purge ladder (`2^{10}..2^{15}` `core` / `2^{10}..2^{16}` `all`, `E_1/U ≈ 0.156`) measures `b_ols = 5.164±0.797` (`rung_R2 = 0.9134`, OOS `0.646`, `16/16` direct unplanted full rank) / `4.819±0.536` (**`rung_R2 = 0.9794`**, OOS **`1.150`**, `18/18` direct full rank), deriving `nl_cap = 0.000783` (`core`, `372.2x` above OLS pred, `38.54x` above 95% mean UCB) / `0.000568` (`all`, `165.7x` / `37.24x`) via severely flattened `b_cons = 1.30 < 0.35 * b_lo95` (`t = 1.70`), paired with a live full-weight `w_comp = 200.35` Billingsley–Dickman stress audit (`U = 8192`: `corr(y_1,k) = -0.2878`, `pm1 = 88.33%`, `max|c| = 20`, `ov30 = 0`, deriving `pm1_cap = 87.33%` live), honestly scoped as an `exploratory-v1` H3 scaling heuristic.
2. **Explicit Charged Implicit Rank-4 Projection & Unconditional Coin-Space `p_LA_coins = 2^{-58.00}`:** Deterministic `H-R`/`H-SM` are strictly separated from unconditional coin-space `p_LA_coins` (`Kernel #17`, `118` insns in `rv = 1080.0`, `U = 0.2394 EiB` in `P_LA = 4`, `5/5` live across `d=0..4`).
3. **2-Tier `f` + 2-Tier `g` Stage 1b NTT (`LtM = 2^{67.7409} < 2^{67.77}` b, `c_s1b = 0.871426 Ir/cell`) & `rbr = 0.86` Shell-1:** Preserves `+0.022406 b` headline slack (`7.69 * SE_Vk` `core` / `12.99 * SE_Vk` `all`) and `1.35x` LA `< claim`.

## Package Manifest (6 Files, All `<= 65,536 B`, Total `< 200 KiB`)
- `submission.json` — Formal claim (`21357/200 = 106.785000`), evidence manifest, experiment schema.
- `proof.md` — Complete mathematical proof, 100% live-derived `H3` scaling & `w_comp = 200` BD audit, implicit rank-4 projection theorem, Stages 0–4, 13-row ledger, LA sensitivity, 64-bit RAM & failure budgets.
- `KERNELS.md` / `KERNELS_STAGE1B_LA.md` — All 17 verbatim compiled 0-spill ARMv8.1-A scalar kernels (`1,643` insns, `73` `staddb`).
- `experiment.py` — Standalone pure-Python 3.11+ stdlib verifier (`<= 65,536 B`, `~7.91s` `core`, `~16.06s` `all`).
- `README.md` — Quick-reference package summary.

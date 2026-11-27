# BLS12-381 FAPI-1 (`bls12-381-fixed-g1-invert-g2`): `stnfs-q15-fixedva-banked-r4ntt-v5`

**Headline Claim (`classical-group-work-v2`, `exploratory-v1`):**
- **Expected Cold-Run Classical Work:** $\mathbf{2^{106677/1000} = 2^{106.677000}}$ target-group-operation equivalents (`C = 20,000` scalar ARMv8.1-A general-purpose-register instructions per unit; `-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`).
  - Charged gated ledger: **`2^{106.672848} u`** (`+0.004152` bits headline slack; `-0.111321` bits below live `#1` `106.788321` u).
  - Live audited gated ledger: **`2^{106.672163} u`**.
  - Charged opened-run ledger: **`2^{108.031273} u`** (live opened-run: **`2^{108.030588} u`**).
- **Cold Success Probability:** $p_{\text{cold}} = p_{\text{gate}} \cdot p_{\text{open}} = 0.39000780 \times 0.99998441 = \mathbf{0.39000172 \ge 0.39}$ (conditional on explicit `exploratory-v1` heuristics H0–H5, H-R, H-SM, H-ECM; pre-gate `K_gate = 7194374076228249696 = 0x63d78d1adadd2860`; Markov pointwise abort meter $W_{\text{pointwise}} = 2^{124}\text{ u}$).
- **Peak Total Storage (100% 64-Bit Banked, Every Owner `< 16 EiB`):** $\le 2^{117.50}\text{ B}$ (Stage 1 Shell 0: $2^{116.7442}\text{ B}$ across `10,058,276,997,840,928` prime banks at $\le 12.000\text{ EiB/bank}$, freed before Shell 1 $2^{114.3913}\text{ B}$; Stage 1b pow2-rounded cached NTT tree: $2^{73.2793}\text{ B} = 4.8544\text{ EiB/owner}$ across 2,048 owners; Stage 2 SGE: $2^{70.3861}\text{ B}$; Stage 3 64-bit `P_LA = 4` banked LA + full Thomé `Lingen` $\pi(z)$: $2^{76.4917}\text{ B}$, $3.2827\text{ EiB/owner}$ across 4 owners). Nonuniform advice: `0 B`.

## Key Mathematical & Verification Upgrades in `-v5`
1. **Unified Self-Consistent Separate-$(B_f, B_g)$ Stage 4 Descent Termination:** One unified Stage 4 engine computes branching counts directly from `_dnode`'s own `(nft, ngt)` (`nodes = [9.3, 55.5, 343.4, 2662.9, 24534.3, 0.0]`), terminating Layer 5 strictly at $B_f = 2^{70.80}$ and $B_g = 2^{71.16}$ with **`0` outgoing nodes on BOTH sides** (`2^{96.6088}` u live / `96.6500` charged; ECM `2^{96.0120}` u live / `96.0500` charged).
2. **Dusart-Exact / GRH-Safe PIT Factor Base & Unconditional Degree-$\ge 2$ Bound:** Uses Dusart (2010, Thm. 6.9) $\pi_{\text{Dus}}(B)$ plus unconditional degree-$\ge 2$ bound $36\,\pi_{\text{Dus}}(B_f^{1/2}) + 12\,\pi_{\text{Dus}}(B_g^{1/2})$ (no `/12`, $U \le 8.11671022\times 10^{18} = 2^{62.815601}$ orbits, with all $p^f \le B$ ideals factored by Stage 1b's $p \le 2^{36}$ tree), Dusart $\theta(x)$ upper bounds in Stage 1b `ppb`, live algebraic Murphy $\alpha_f = +5.3202\text{ b}, \alpha_g = -1.0339\text{ b}$ ($p \le 5000$), and conservative `-0.04` b smoothness cushion (`15.3x` the $p > 5000$ tail).
3. **100% Live-Reproducible 2D Panjer DP, Out-of-Sample 16D MC & Live $L_{0.65}(b)$ Chords:** `experiment.py` builds all 4 2D Panjer DP tables from scratch, fits the discrete Lagrange $(T_g, T_f)$ policy (`max_byte <= 255`), evaluates that actual policy out-of-sample on an independent 16D MC (**zero assumed retention constant**), and live-samples $L_{0.65}(b)$ chords directly reproducing `(ret, log2_sg, log2_sf, log2_s1, Vk, L4, gap4)`.
4. **Faithful 4-Tuple $\mathbb{Z}[\zeta_{12}]$ `STRICT` SGE + 64-Seed Live $\mathbb{F}_r$ Rank + Conservative $p_{\text{retry}} = 1/16$:** `experiment.py` runs genuine 4-tuple $\mathbb{Z}[\zeta_{12}]$ `STRICT` monomial-pivot SGE ($\pm z^k$ pivots only, `0` ring muls, `max|c| = 8 << 2^{30}`, `0` overflows, `79.41%` monomials) + exact BLS12-381 $\mathbb{F}_r$ Gaussian elimination across **`64/64` live seeds** (`0/64` failures $\implies$ exact Clopper–Pearson 95% UB `4.573%` one-sided / `5.601%` two-sided, both `< 6.25%`), grounding the conservative `p_retry = 1/16 = 6.25%` Snell resume tail (`2^{101.9690}` u).
5. **All 16 Compiled Literally Zero-Spill ARMv8.1-A Kernels (`1,602` Executable Scalar Instructions):** Preserved verbatim in `KERNELS.md` and `KERNELS_STAGE1B_LA.md` (`simd = 0`, `0` forbidden ISA, `0` `udiv`/`sdiv`, `0` helper calls, `0` in-loop spills anywhere).

## Package Manifest (6 Files, All `<= 65,536 B`, Total `< 180 KiB`)
- `submission.json` — Formal claim (`106677/1000 = 106.677000`), evidence manifest, experiment schema.
- `proof.md` — Complete mathematical proof, explicit `exploratory-v1` heuristic scoping, Stage 0–4 algorithm, 14-row ledger, RAM & failure budgets.
- `KERNELS.md` — Verbatim compiled ARMv8.1-A assembly (Kernels 1–6: `615` insts, `73` `staddb`, `0` in-loop spills).
- `KERNELS_STAGE1B_LA.md` — Verbatim compiled ARMv8.1-A assembly (Kernels 7–16: `987` insts, `0` `staddb`, `0` in-loop spills).
- `experiment.py` — 100% standalone pure-Python 3.11+ stdlib verifier (reads 0 external files, `~4s` core, `~5s` all).
- `README.md` — Package overview and quick-reference summary.

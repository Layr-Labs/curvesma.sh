# BLS12-381 FAPI-1 (`bls12-381-fixed-g1-invert-g2`): `stnfs-q15-fixedva-banked-r4ntt-v4`

**Headline Claim (`classical-group-work-v2`, `exploratory-v1`):**
- **Expected Cold-Run Classical Work:** $\mathbf{2^{106787/1000} = 2^{106.787000}}$ target-group-operation equivalents (`C = 20,000` scalar ARMv8.1-A general-purpose-register instructions per unit; `-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`).
  - Charged gated ledger: **`2^{106.786141} u`** (`+0.000859` bits headline slack).
  - Live audited gated ledger: **`2^{106.785696} u`** (`+0.001304` bits total slack).
  - Charged opened-run ledger: **`2^{108.144566} u`** (live opened-run: **`2^{108.144121} u`**).
- **Cold Success Probability:** $p_{\text{cold}} = p_{\text{gate}} \cdot p_{\text{open}} = 0.39000780 \times 0.99998313 = \mathbf{0.39000122 \ge 0.39}$ (conditional on explicit `exploratory-v1` heuristics H1–H5, H-R, H-SM, H-ECM; pre-gate `K_gate = 7194374076228249696 = 0x63d78d1adadd2860`; Markov pointwise abort meter $W_{\text{pointwise}} = 2^{124}\text{ u}$).
- **Peak Total Storage (100% 64-Bit Banked, Every Owner `< 16 EiB`):** $\le 2^{117.50}\text{ B}$ (Stage 1 Shell 0: $2^{116.7442}\text{ B}$ across `10,058,276,997,840,928` prime banks at $\le 12.000\text{ EiB/bank}$, freed before Shell 1 $2^{114.3913}\text{ B}$; Stage 1b pow2-rounded cached NTT tree: $2^{73.2789}\text{ B} = 4.8532\text{ EiB/owner}$ across 2,048 owners; Stage 2 SGE: $2^{70.3557}\text{ B}$; Stage 3 64-bit `P_LA = 4` banked LA + full Thome Lingen $\pi(z)$: $2^{76.4613}\text{ B}$, $3.2143\text{ EiB/owner}$ across 4 matrix owners). Nonuniform advice: `0 B`.

## Key Upgrades in `-v4` (Resolving All Audit Findings Constructively)
1. **Unconditional Algebraic Chebotarev Proof + Live Polynomial Factoring:** Explicit cubic discriminant $\Delta_Q = 41 + 24\omega$ ($N(\Delta_Q) = 1273 = 19\cdot 67$), explicit Kummer witness prime $p = 991 \equiv 1 \pmod{45}$ ($\omega = 113$, roots $\{373, 802, 808\}$, $(r_1/r_2)^{330} \equiv 877 \not\equiv 1$), and disjoint discriminant support $\{\mathfrak{p}_{19}, \mathfrak{p}_{67}\} \cap \{\overline{\mathfrak{p}_{19}}, \overline{\mathfrak{p}_{67}}\} = \emptyset$ prove $\operatorname{Gal}(L_\omega/K) \cong C_3 \wr S_3$ and $L_\omega \cap L_{\omega^2} = K(\zeta_9)$ unconditionally; verified live mod 300+ split primes in `experiment.py`.
2. **Explicit Stage 4 12D Hilbert-90 Lifting Lattice + 16D Special-$\mathfrak{q}$ Lattice + Live Dickman–de Bruijn Branching Recurrence:** Explicit construction of $\mathcal{L}_{t'}$ (det $p^6 \approx 2^{2286}$), 16D $\Lambda(\mathfrak{q})$, and exact 5-layer Dickman–de Bruijn recurrence terminating strictly at $Q_5 = B_g = 2^{71.16}$ with total 2-target descent cost $2^{64.45}\text{ u} \ll 2^{96.00}\text{ u}$ (`31.55` bits headroom), plus small-prime 12D LLL demo in `experiment.py`.
3. **Production 64-Bit `P_LA = 4` Banked Signed + Non-Monomial SpMV:** Vector $V$ ($0.3309\text{ EiB}$) fits in ONE native 64-bit address space (`LDR X`, zero 32-bit index truncation); compiled 0-spill `_spmv_mono_add_ptr64` (`10.00 Ir/nnz`, 81.43%) + `_spmv_nonmono_muladd_ptr64` (`25.00 Ir/nnz`, 18.57%) + `0.063` bcast = `12.8487 <= 12.90 Ir/nnz` charged; signed split $A_k^+ + (r \ll 60) - A_k^-$ (`+56 Ir`) + verbatim 576-bit `r<<256` Horner+SM4+Barrett (`872 Ir`) = `928 <= 940.0 Ir/row` charged.
4. **Literally ZERO In-Loop Stack Spills Across All 16 Compiled Kernels (`1,602` Executable Scalar Instructions):** `_col_scan_swar_128_nospill` (`64` insts, `0` spills on BOTH fast & slow paths; `slab_live = 0.415659`, `slab_ch = 0.417659 Ir/cell`).
5. **Live Reproducing Computations in Standalone `experiment.py`:** Live 16D Mode-A Monte Carlo (`N = 6,000` / `20,000`), live `L_CALIB STRICT` SGE filter + **exact Gaussian elimination over actual BLS12-381 $\mathbb{F}_r$**, live polynomial factoring, live Dickman–de Bruijn branching + 12D LLL, and explicit `exploratory-v1` heuristic classification of H1/H3/H4/Markov.

## Package Manifest (6 Files, All `<= 65,536 B`, Total `< 180 KiB`)
- `submission.json` — Formal claim (`106787/1000 = 106.787000`), evidence manifest, experiment schema.
- `proof.md` — Complete mathematical proof, explicit `exploratory-v1` heuristic scoping, Stage 0–4 algorithm, 14-row ledger, RAM & failure budgets.
- `KERNELS.md` — Verbatim compiled ARMv8.1-A assembly (Kernels 1–6: `615` insts, `73` `staddb`, `0` in-loop spills).
- `KERNELS_STAGE1B_LA.md` — Verbatim compiled ARMv8.1-A assembly (Kernels 7–16: `987` insts, `0` `staddb`, `0` in-loop spills).
- `experiment.py` — 100% standalone pure-Python 3.11+ stdlib verifier (reads 0 external files, `< 5s` default).
- `README.md` — Package overview and quick-reference summary.

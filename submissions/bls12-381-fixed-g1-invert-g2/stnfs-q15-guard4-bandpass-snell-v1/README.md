# `stnfs-q15-guard4-bandpass-snell-v1`

Clean-room exploratory submission for the **`bls12-381-fixed-g1-invert-g2`** track (`classical-group-work-v2`, $C = 20{,}000$ scalar ARMv8.1-A general-purpose-register instructions per dense $\mathbb{F}_{p^{12}}$ multiplication unit, `simd = 0`, single cold target, $p_{\text{cold}} \ge 0.39$).

## Headline Result

| Quantity | Value |
|---|---|
| **Official Reduced Claim (`score`, `gcd = 1`)** | **`10667/100` = `106.670000` ($\log_2 u$)** |
| Charged Cold-Run Gated Ledger ($p_{\text{cold}} \ge 0.39$) | `106.669708` $\log_2 u$ (`+0.000292` b slack to claim) |
| Live Cold-Run Gated Ledger | `106.669534` $\log_2 u$ |
| Pure First-Principles Gated (without `+0.06544` b anchor pad) | `106.613505` $\log_2 u$ |
| Opened-Run Expected Work ($p_{\text{open}} \ge 0.99998444 > 0.99998$) | `108.028134` $\log_2 u$ charged (`108.027959` live) |
| Unconditional Markov Pointwise Meter ($W_{\text{pointwise}}$) | `2^{124.0}` u ($\mathbb{P}(\text{abort}) \le 1.556 \times 10^{-5}$ unconditionally!) |
| Minimal 64-Bit Pre-Gate Threshold $K_{\text{gate}}$ | `7194374076228249696` (`0x63d78d1adadd2860`) |
| Peak Unified Global Storage (`classical-group-work-v2`) | `2^{117.2595}` B (`<= 2^{117.50}` B; Stage 2 SGE `2^{69.08}` B; Stage 3 LA `2^{64.42}` B; Stage 1b `2^{66.22}` B) |
| Nonuniform Advice | `0` bytes (`2^0` B) |
| ISA Profile | `-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize`, `simd = 0`, scalar LSE `STADDB Ws, [Xn]` disclosed |

## Key Algorithmic Contributions

1. **Guard4 x16 ($W_B = 4$) Merged-$\mathcal{O}_K$-Prime Sieve + Sequential Single-Slab Active-Only Masked Scan:** Two canonical $a$-shells ($K = 2$, $V_{a,\text{sh}} = 0.9351 \cdot 2^{56}$, $m_{\text{unit}} \equiv 1$ identically) paired with $b \in L_{0.65}(b)/C_3$ along 1D rational chords (`4`-byte left guard `[0..3]`, `L <= 248` active columns `[4..4+L-1]`, `>= 4`-byte right guard `[4+L..255]`). ONE single `14.9616` EiB flat slab is reused sequentially ($S_g \to S_f$): Pass 1 pre-biases active $S_g$ (`_slab_prebias_active_sl`), runs Guard4 x16 (`43/16 = 2.6875` Ir/hit, `13/4 = 3.25`, `5.00` tail), executes a guard-masked SWAR active-only scan (`_slab_scan_g_active_sl`, masking BOTH left `[0..3]` and right `[4+L..255]` guards!) emitting `2^{-14.41}` $g$-survivors; Pass 2 reuses the SAME slab for $S_f$ and probes ONLY at $g$-survivor coordinates (`_slab_filter_f_at_gsurv`, `0.00023` Ir/active cell). Total active slab overhead: `0.5318` exact (`0.5500` Ir/active cell charged).
2. **Exact 2D Panjer Smooth-Recall + Survival DP & Zero-Wraparound Proof:** Deterministic 2D Panjer DP (zero MC noise, exact Murphy $\alpha_f = +5.3752, \alpha_g = -1.0487$) certifies `g = [32, 61.5]` (`q=3.8, T_u8=68, bias=60, max_byte=245 <= 255, rec_g=0.80782, sg=2^{-14.4064}`) and `f = [36, 61.0]` (`q=3.8, T_u8=45, bias=83, max_byte=248 <= 255, rec_f=0.85091, sf=2^{-9.1325}`), yielding `ret_exact = 0.68738` and `s1 = 2^{-23.5389}`.
3. **$f$-First Uncached (`2.00x`) Bernstein Batch Cofactorization:** Uses unconditional splitting densities (`1/8` on $f$, `1/24` on $g$), exact bit-unit `L_owner_bits = (max_owner_bytes * 8) / 4.5`, exact `Top2` `tier_probs`, both $N_f$ and $N_g$ evaluations, and cleanly separates `c_cof_rem_mr` (`0.6208` Ir/cell in `relation_collection`) from `tree_build_u` (`band_split_and_trees`), achieving `ipc_tot = 4.8886` Ir/cell.
4. **Monomial-Pivot SGE ($U = 2^{20}$ at $\eta = 0.82$) + Full Coppersmith BW (`d_mksol = N_l/32 + 32`, `IR_ROW_FIXED = 2100`):** Directly measured $U = 2^{20}$ $\eta$-sweep confirms `LA_ratio = 2.263x` at $\eta = 0.82$ (vs charged `4.125x` $\implies$ `1.82x` safety margin), `Mono = 85.32%`, `0.00%` `>32`-bit overflow, `max|c| = 61`. Primary ledger charges exact `d_mksol = N_l/32 + 32`, dynamic `_fr_dot4_mont` (`747`), and full `IR_ROW_FIXED = 2,100` (`linear_algebra = 2^{103.1630}` u charged).
5. **Stage 4 12D Subfield $M_g = \mathbb{Q}(\zeta_5)[Y]/(Y^3 - x)$ Hilbert-90 Descent:** Unconditional deterministic LLL + AM–GM norm proof (`det = p^6 |x|^4 = 2^{2539.061}`, `N_boot <= 2570.51 < 2572.0` b), `B_boot = 2^{275}`, explicit Layer-0 (`2^{550} -> 2^{275}`) + 20 extra inert-prime nodes for $q \equiv 2 \pmod 3$ norm squaring, unconstrained Torsion-16 Edwards ECM table (`q = 2^{70}..2^{275}`, `< 10^{-12}` self-consistent), live 2-log with Snell `(1.04094/0.82)^2` = `2^{95.778}` u (`2^{96.5001}` u charged), explicit `ecm_split_expected = 2^{94.0707}` u.

## Package Contents

| File | Description |
|---|---|
| `submission.json` | `curvesmash-v1` manifest (`score = 10667/100`, `gcd = 1`, `required: []`) |
| `proof.md` | Complete 6-stage mathematical & engineering proof, unified RAM audit, unconditional Markov meter, 14-row ledger, heuristics & sensitivity table |
| `KERNELS_SIEVE.md` | Full untruncated ARMv8.1-A assembly (`simd = 0`) for Guard4 x16 sieve, 2-pass Guard4 replay builder, and sequential single-slab masked active scan |
| `KERNELS_LA_ECM.md` | Full untruncated ARMv8.0-compatible assembly (`simd = 0`) for 64-bit-index SpMV, 255-bit $\mathbb{F}_r$ arithmetic, Stage 1b Shoup61 NTT/CRT3, and Stage 4 Torsion-16 ECM |
| `FILTER_RANK.md` | `hx_z12` rank-barrier proof & 2-merge resolution, Snell-Envelope Theorem 1 ($\mathbb{P}(\text{fail }\tau_0) < 2^{-11}$), and $U = 2^{14}..2^{20}$ Monomial-Pivot SGE ladder incl. $U = 2^{20}$ $\eta$-sweep |
| `DESCENT_ECM.md` | Stage 4 12D subfield $M_g$ Hilbert-90 lift, `q = 2^{70}..2^{275}` Torsion-16 Edwards ECM table (`< 10^{-12}`), `ecm_split_expected` derivation, and H5 witnesses |
| `ledger_data.json` | Machine-readable authoritative ledger data synchronized with `proof.md` and `experiment.py` |
| `experiment.py` | Pure Python 3.11+ stdlib clean-room verifier (`core`, `diagnostics`, `all` modes) |
| `README.md` | Package summary and reproduction commands |

## Reproduction Commands

```sh
python3 -I -B experiment.py <<< '{}'
python3 -I -B experiment.py <<< '{"mode": "diagnostics"}'
python3 -I -B experiment.py <<< '{"mode": "all"}'
```

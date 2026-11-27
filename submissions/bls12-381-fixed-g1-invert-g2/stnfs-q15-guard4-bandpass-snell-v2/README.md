# `stnfs-q15-guard4-bandpass-snell-v2`

Clean-room exploratory submission for the **`bls12-381-fixed-g1-invert-g2`** track (`classical-group-work-v2`, $C = 20{,}000$ scalar ARMv8.1-A general-purpose-register instructions per dense $\mathbb{F}_{p^{12}}$ multiplication unit, `simd = 0`, single cold target, $p_{\text{cold}} \ge 0.39$).

## Headline Result (`CAND2_FINAL_ETA85_7000_NL205`)

| Quantity | Value |
|---|---|
| **Official Reduced Claim (`score`, `gcd = 1`)** | **`21337/200` = `106.685000` ($\log_2 u$)** |
| Charged Cold-Run Gated Ledger ($p_{\text{cold}} \ge 0.39$) | **`106.682627`** $\log_2 u$ (`+0.002373` b slack; **`-0.103321` b** below `bd853c58` `106.788321`!) |
| Live Cold-Run Gated Ledger | `106.682499` $\log_2 u$ |
| Pure First-Principles Gated (excl. `+0.06544` b anchor pad) | **`106.628735`** $\log_2 u$ |
| Opened-Run Expected Work ($p_{\text{open}} \ge 0.99998430 > 0.99998$) | `108.041052` $\log_2 u$ charged (`108.040924` live) |
| Complete Union-Bound Failure Budget ($W_{\text{pointwise}} = 2^{124.0}$ u) | $p_{\text{fail,open}} \le 1.5699 \times 10^{-5} < 2.0 \times 10^{-5}$, $p_{\text{cold}} \ge 0.39000168 > 0.39$ |
| Minimal 64-Bit Pre-Gate Threshold $K_{\text{gate}}$ | `7194374076228249696` (`0x63d78d1adadd2860`) |
| Peak Unified Flat RAM (`classical-group-work-v2`) | **`2^{117.0082}` B** (`<= 2^{117.50}` B `<< 2^{120}` B; Stage 1b `2^{69.43}` B; Stage 2 `2^{69.52}` B; Stage 3 `2^{72.31}` B; single slab `14.9616 EiB < 16 EiB`) |
| Nonuniform Advice | `0` bytes (`2^0` B) |
| ISA Profile | `-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize`, `simd = 0`, zero stack spills in `_emit_u64_x64_direct` and slab hot loops, scalar LSE `STADDB Ws, [Xn]` disclosed |

## Key Breakthroughs & Complete Resolution of Audit Objections

1. **Exact 1D + 2D Bivariate 3-Prime 61-Bit Shoup NTT & $k_{\text{pack}}$-Way Packed Super-Leaves (`KERNELS_LA_ECM.md`):** Uses three explicit verified $2^{54}$-smooth 61-bit Shoup primes ($p_1 = 127\cdot 2^{54}+1$, $p_2 = 57\cdot 2^{55}+1$, $p_3 = 27\cdot 2^{56}+1$, $P = 2^{182.5765}$, $4 p_i < 2^{63}$, primitive roots `3, 7, 5`; recompiled 1D `_crt3_recombine_carry` = `58 Ir/word <= 67` charged; 2D `_garner_crt3_overlap_carry58` = `73 Ir/coeff <= 95` charged). Prices Stage 1b at honest unconditional survivor widths (`f1 = 560` b w/ `k_pack = 7`, `g1 = 675` b w/ `k_pack = 6`, unclamped `f2,f3 = 385` b & `g2,g3 = 405` b w/ `k_pack = 10`), explicit bottom super-leaf build+descent, `L_own_2d = 2^{61.83}` b, `9*N_2d` 2D transposes, and full stored-tree RAM `2^{69.4261}` B.
2. **Zero-Extrapolation $U = 2^{24}$-Certified SGE Pricing at $\eta = 0.85$ + Horner BW `(512, 32)` (`SGE_CHEBOTAREV_STAGE4.md`):** Discloses ALL 9 dyadic scales `U = 2^{16}..2^{24}` (`16.78M` orbits) of strict Monomial-Pivot $\mathbb{Z}[\zeta_{12}]$ SGE (`4`-seed $U = 2^{23}$: `0.001380 +- 0.000106`, `10.56 Ir/nnz`; `2`-seed $U = 2^{24}$: `0.001738 +- 0.000111`, `0.32249`, `10.21 Ir/nnz`). Charged caps **`nl_cap = 0.00205, w_cap = 190.0`** sit strictly above `mean_24 + 2*SE = 0.001960` (`47x` above the 9-rung OLS projection `0.000044` at $2^{62.0}$!) with zero extrapolation. Horner BW `(512, 32)` charges `x4` Krylov SpMVs (`12.141 Ir/nnz`), unbatched `x1` `Mksol` SpMVs (`16.000 Ir/nnz`), and full Thomé `Lingen` $\pi(z)$ RAM `2^{72.3056}` B (`linear_algebra = 2^{104.4205}` u).
3. **Unconditional Galois Chebotarev Theorems 2A & 2B (`SGE_CHEBOTAREV_STAGE4.md`):** Proves unconditionally $\operatorname{Gal}(\operatorname{Spl}(f_\omega)/K) \cong C_3 \wr S_3$ (`162`) and $\operatorname{Gal}(L_\omega L_{\omega^2}/K) \cong (C_3 \wr S_3) \times_{C_3} (C_3 \wr S_3)$ (`8,748`), yielding exact `d_ok_f = 23/81`, `r_bar_f = 81/23`, and `dens_f = 2131/34992` (matching all `1,857,859` primes $p \le 30\text{M}$ to `8 ppm`!).
4. **2D Adaptive `(T_g, T_f)(shell, N_g(b), nf0(b))` + Zero-Spill Guard4 `x64` (Exact Poisson) + 16B-Granular Adaptive Slab (`BANDPASS_ADAPTIVE.md`, `KERNELS_SIEVE.md`):** Compiled zero-spill `_emit_u64_x64_direct` (`163/64 = 2.546875 Ir/hit`, `x8 = 23/8`, `x1 = 5.00`, integrated via exact Poisson expectation over Chebotarev root counts and split-256 slices: `ipc_hits = 3.490912`, `ipc_probe = 0.320852`) + 16B-granular 2-row adaptive slab (`0.4777` live $\le$ `0.4950` charged, `248/248` bit-exact under `0xFF` guard poison, `0` spills) $\implies$ **`ipc_tot = 5.015016 Ir/cell`**.
5. **100% Standalone Docker-Ready `experiment.py` + Complete Union-Bound Failure Budget:** Embeds `LEDGER_DATA` directly with zero required external file reads, passing isolated single-file execution in `< 1s` (`core`).

## Package Contents (9 Files, All `<= 65,536` B)

| File | Description |
|---|---|
| `submission.json` | `curvesmash-v1` manifest (`score = 21337/200`, `gcd = 1`, `required: []`) |
| `proof.md` | Complete 6-stage mathematical & engineering proof, exact Poisson emission, packed super-leaf Stage 1b, union-bound failure budget, unified RAM audit, 14-row ledger, heuristics & stress tests |
| `KERNELS_SIEVE.md` | Compiled ARMv8.1-A assembly (`simd = 0`) for zero-spill `_emit_u64_x64_direct`, `_sieve_window_guard4_x64`, 2-pass replay builder, and 16B-granular zero-spill adaptive slab |
| `KERNELS_LA_ECM.md` | Compiled ARMv8.0-compatible scalar assembly (`simd = 0`) for 1D (`58 Ir`) + 2D (`73 Ir`) $2^{54}$-smooth Shoup61 NTT/CRT, 64-bit-index SpMV, 255-bit $\mathbb{F}_r$ arithmetic, and Torsion-16 ECM |
| `BANDPASS_ADAPTIVE.md` | Exact 2D Panjer DP, 2D `(shell, N_g(b), nf0(b))` adaptive thresholding, deterministic `max_byte <= 255` proof, `4M`-MC & `200k`-chord certification, and kill-gates |
| `SGE_CHEBOTAREV_STAGE4.md` | Unconditional Galois Chebotarev Theorems 2A & 2B, 9-rung `U = 2^{16}..2^{24}` SGE ladder, `U = 2^{14}..2^{16}` $\mathbb{F}_q$ rank audit, Horner BW `(512, 32)`, and 12D $M_g$ Stage 4 + ECM |
| `ledger_data.json` | Machine-readable authoritative ledger data (synchronized with embedded `LEDGER_DATA` in `experiment.py`) |
| `experiment.py` | **100% standalone** pure Python 3.11+ stdlib verifier (`core`, `diagnostics`, `all`) |
| `README.md` | Package summary and reproduction commands |

## Reproduction Commands

Runs cleanly both inside the package directory AND completely standalone (isolated copy of `experiment.py` alone):

```sh
python3 -I -B experiment.py <<< '{}'
python3 -I -B experiment.py <<< '{"part": "diagnostics"}'
python3 -I -B experiment.py <<< '{"part": "all"}'
```

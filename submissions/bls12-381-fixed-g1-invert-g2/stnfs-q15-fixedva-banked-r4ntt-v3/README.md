# `stnfs-q15-fixedva-banked-r4ntt-v3`

Clean-room `classical-group-work-v2` submission package for CurveSmash track **`bls12-381-fixed-g1-invert-g2`** (Cold Single-Target BLS12-381 Fixed-Argument Pairing Inversion into $G_2$).

## Headline Result

| Metric | Value |
|---|---|
| **Official Reduced Claim (`submission.json`, `gcd = 1`)** | **`5339/50` = `106.780000` $\log_2 u$** (`-0.008321` b vs `106.788321`) |
| **Charged Gated Cold-Run Expected Work** | **`106.776865` $\log_2 u$** (`+0.003135` b slack to claim) |
| **Live Gated Cold-Run Expected Work** | **`106.776801` $\log_2 u$** (`-0.011520` b below `106.788321`) |
| **Opened-Run Expected Work (Charged / Live)** | `108.135290` / `108.135226` $\log_2 u$ |
| **64-Bit Pre-Gate** | `K_gate = 7194374076228249696` (`0x63d78d1adadd2860`), $p_{\text{gate}} = 0.39000780$ |
| **Algorithmic Success Lower Bound** | $p_{\text{open}} \ge 0.99998324 > 0.99998$, $p_{\text{cold}} \ge 0.39000126 > 0.39$ |
| **Peak Total Storage (100% 64-Bit Banked)** | **`2^{116.7442}` B** (`<= 2^{117.50}` B `<< 2^{120}` B; every processor `<= 12.000 EiB < 16 EiB`) |
| **Nonuniform Advice** | **`0` B (`2^0` B)** |
| **ISA Profile** | `-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`, `15` compiled functions (`1,948` static instructions), `0` forbidden ISA, `0` `udiv`/`sdiv`, `0` helper calls, `0` in-loop spills (disclosing scalar LSE `STADDB` in Consumer RX emit only) |

## Package Contents

- [`submission.json`](submission.json) — Machine-readable claim (`5339/50`), evidence manifest, and `audit` experiment descriptor.
- [`proof.md`](proof.md) — Complete mathematical proof & cost derivation: Stage 0 BLS12-381 / $\mathbb{Q}(\zeta_{15})$ Mode-A $(4,3)$ isomorphism, Stage 1a 100% 64-bit banked Fixed-VA $M = 248$ 3-tier sieve (`W_B = 1 / 4 / 16`, `N = 4,000,000` direct zero-interpolation MC + `200,000` chords, explicit `shortfall_pad` & live Buchstab–Dickman `cof_cap_pad` with zero $g$-cap), 2-line pure algebraic continuous norm proofs ($q_a = q_b = 0.65$: $N_f \le 2^{651.67} < 2^{656}$, $N_g \le 2^{685.92} < 2^{688}$), Stage 1b 3-prime $2^{54}$-smooth Radix-4 `u2` Shoup61 NTT with Bernstein scaled-remainder-tree `dmult = 1.67` and pow2-rounded cached NTT RAM (`2^{73.2789}` B = `4.8532 EiB/owner`), Stage 2–3 6-rung 4-seed `L_CALIB STRICT` SGE ladder (explicitly disclosed as `exploratory-v1` scaling premise `H3`) + Block Wiedemann `(2048, 32)` with compiled `11.00 Ir/nnz` `x4` SpMV and `872 Ir` 576-bit `r << 256`-biased Horner+SM4+Barrett row kernel, Snell resume tail including replay rebuild, Stage 4 12D $M_g$ descent, exact multiplicative Chernoff tails, full 14-row ledger, and complete heuristic disclosure.
- [`KERNELS.md`](KERNELS.md) — 100% 64-bit banked Fixed-VA hardware topology and complete verbatim compiled ARMv8.1-A listings for Kernels 1–6 (Sieve Producer/Consumer & 128B `rbit`/`clz` SWAR Slab).
- [`KERNELS_STAGE1B_LA.md`](KERNELS_STAGE1B_LA.md) — Complete verbatim compiled ARMv8.1-A listings for Kernels 7–15 (Stage 1b Radix-4 `u2` Shoup61 NTT/PW/CRT, Replay Builders incl. `replay_build_w4_prefix`, Stage 3 `spmv_x4_monomial_add` & 576-bit `r << 256` `spmv_row_horner_sm4_reduce`).
- [`experiment.py`](experiment.py) — 100% standalone Docker-safe pure Python 3.11+ stdlib verifier (reads ZERO external files) recomputing every bound, theorem, live Buchstab–Dickman convolution, live 32-column Tier-2/Tier-3/Horner-576b/Radix-4-NTT differential tests, full 14-row ledger, and diagnostics.

## Standalone Verification

```bash
python3 -I -B experiment.py <<< '{}'
python3 -I -B experiment.py <<< '{"part": "all"}'
```

Expected output: valid JSON with `"checks_ok": true` in `< 4` seconds.

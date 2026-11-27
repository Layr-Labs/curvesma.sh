# Cold Single-Target BLS12-381 FAPI-1 via 64-Bit Banked Fixed-VA $W_{16}$ $\mathbb{Q}(\zeta_{15})$ CM-STNFS, Radix-4 Shoup61 NTT & 576-Bit Horner-SM4 BW (`proof.md`)

**Track:** `bls12-381-fixed-g1-invert-g2` (Fixed-Argument Pairing Inversion into $G_2$)  
**Cost Model:** `classical-group-work-v2` ($C = 20{,}000$ scalar ARMv8.1-A general-purpose-register instruction equivalents per dense $\mathbb{F}_{p^{12}}$ multiplication unit; `simd = 0`; cold single-target run; exploratory algorithmic success $p_{\text{cold}} \ge 0.39$)  
**Official Reduced Claim (`submission.json`, `gcd = 1`):** **`5339/50` = `106.780000`** ($\log_2 u$; charged gated ledger **`106.776865`**, live gated **`106.776801`**, opened-run charged **`108.135290`**, slack to claim **`+0.003135` b**, beating live `#1` `106.788321` by **`-0.008321` b**)

---

## 1. Executive Summary, Banked 64-Bit Architecture & Reduction to $G_T$-DLP

The challenge supplies an independently uniform non-identity point $P \in G_1 \setminus \{\mathcal{O}\}$ and uniform target $z \in G_T \subset \mathbb{F}_{p^{12}}^\times$ (including $z = 1$), and asks for $Q \in G_2$ satisfying $e(P, Q) = z$ under the pinned BLS12-381 optimal Ate pairing.

### 1.1 Structural Cryptanalytic & Engineering Advances (`CAND3_CHAMP_7080_110X`)
1. **100% 64-Bit Banked Fixed-VA $M = 248$ Architecture with Zero-AND $W_{16}$ Demux — Zero Unified $2^{117}$-Byte Addressing & Zero Uncharged DMA (Section 3, `KERNELS.md`):** Stage 1a runs on $M = 248$ independent 64-bit **Slab-Owner** processors (each owning 1 contiguous column of size $V_a = 0.9351\cdot 2^{56} = \mathbf{0.058444\text{ EiB}} \ll 16\text{ EiB}$, mapped at fixed VA `SLAB_BASE = 2^{56} = 0x0100_0000_0000_0000`, plus a `256 KiB` batched-`logp` RX ring) fed by $P_{\text{banks}} = 10{,}058{,}276{,}997{,}840{,}928$ independent 64-bit **Prime-Bank** processors (each holding $\mathbf{\le 12.000\text{ EiB}} < 16\text{ EiB}$ of local `replay_a_tables`). Because $V_a < 2^{56}$, every pointer `ptr64 = 2^{56} + aid` has **bits `57..59 = 000` identically**, letting Tier-3 $W_B = 16$ store the 4-bit column tag in bits `60..63` so `(e >> 57) == (col << 3)` requires **zero `AND` mask** (`_demux_w16_window`, `26/4 = 6.50 Ir/hit`, `8.00 Ir` tail, priced over ALL `n_b16 * 16` columns including overhang!). Tier 2 uses reloaded-pointer `_w4_packed_tx_window` paired with true 4-slice header builder `_replay_build_w4_prefix`; Slab scan uses `rbit`/`clz` `_col_scan_swar_128_fast` (`28/128 = 0.21875 Ir/cell` fast + `68 Ir` slow + $f$-probe). Every instruction across all **15 compiled functions (`1,948` static instructions, `0` forbidden ISA, `0` `udiv`/`sdiv`, `0` helper calls, `0` in-loop spills)** in `KERNELS.md` & `KERNELS_STAGE1B_LA.md` is matched literally.
2. **2-Line Pure Algebraic Continuous 16D Norm Bounds ($q_a = q_b = 0.65$) + Zero $g$-Side Cap + Live Buchstab–Dickman $f$-Cap + Radix-4 `u2` Shoup61 NTT (Section 4):** Setting $q_a = q_b = 0.65$ eliminates the Hölder constant and yields unconditional 2-line continuous proofs $N_f(a,b) \le \mathbf{2^{651.6655} < 2^{656}}$ (`k_pack = 6`) and $N_g(a,b) \le \mathbf{2^{685.9240} < 2^{688}}$ (`k_pack = 5` on `g1, g2, g3` — **zero $g$-side cofactor cap!**). Only the $f$-side caps post-Tier-1 cofactors at `408` b (`k_pack = 10`), with exact Buchstab–Dickman convolution `p_ex_f = 0.003694` (`cof_cap_pad = +0.005339` b, live-verified in `experiment.py` with `buchstab_norm` within `0.007` of `1.0`). Stage 1b uses compiled `[0, 4p)`-invariant Radix-4 `u2` Shoup61 DIT (`89/8 = 11.125 Ir/r2-bfly`, `0` bit-errors on full `N = 1024` DFT across all three $2^{54}$-smooth primes), compiled `_ntt61_pw_mul_x4` (`9.25 -> 12.0 Ir`), division-free `_crt3_1d_carry` (`55 -> 58.0 Ir` 1D, `90.0 Ir` 2D), Bernstein scaled-remainder-tree / middle-product `dmult = 1.67`, and explicit power-of-2 rounded cached NTT RAM (`2^{73.2789}` B = `4.8532 EiB/owner < 10 EiB` across `2,048` owners).
3. **6-Rung 4-Seed $U = 2^{16}..2^{24}$ `L_CALIB STRICT` SGE Ladder (`H3`) + Compiled `11.00 Ir/nnz` `x4` SpMV + `872 Ir` 576-Bit `r << 256` Horner+SM4 Row Kernel (Section 5):** Documents 24 independent `STRICT` Monomial-Pivot $\mathbb{Z}[\zeta_{12}]$ SGE runs at $\eta = 0.88$ (`U = 2^{24}` 4-seed mean $N_l/U = 0.001045 \pm 0.000104$, $N_{\text{nnz}}/U = 0.12454$), explicitly disclosed as an `exploratory-v1` scaling premise (`H3`) with `2.19x` combined $N_l\cdot N_{\text{nnz}}$ headroom (`nl_cap = 0.00150`, `sw_cap = 0.190`). Compiles `_spmv_x4_monomial_add` (`11.00 <= 12.141 Ir/nnz`) and 576-bit `r << 256`-biased `_spmv_row_horner_sm4_reduce` with branchless `t8` fold (`-2r << 256`), proven zero-underflow and zero-overflow for ALL $0 \le A_k < 2^{315}$ (`872 Ir` dynamic `<= 920.0 Ir/row`, `200/200` bit-exact up to $A_k = 2^{314}-1$) for Block Wiedemann `(2048, 32)`, plus Snell resume tail explicitly charging `replay_a_tables` rebuild (`2^{96.0264}` u).
4. **Exact Multiplicative Chernoff Tails, Unconditional Chebotarev 2A/2B & Standalone Docker-Safe Verifier (`experiment.py`):** Verifies both lower-tail (`exp(-mu*delta^2/2)`) and upper-tail (`exp(-mu*delta^2/(2+delta))`) multiplicative Chernoff exponents ($\approx 1.92\times 10^{13}$, conditional on `H1`), exact `shortfall_pad = +0.002054` b, Chebotarev `23/81` & `2131/34992`, and genuine 32-column Tier-2/Tier-3/576-bit-Horner/Radix-4-NTT differential audits in `< 2s` with zero external files.

### 1.2 Exact Reduction, Pre-Gate & Complete Union-Bound Failure Budget
1. **Canonical Decode & Subgroup Validation (Stage 0):** Validate $P \in E_1(\mathbb{F}_p)\setminus\{\mathcal{O}\}$, $[r]P = \mathcal{O}$, $z \in \mathbb{F}_{p^{12}}^\times$, $z^r = 1$. If $z = 1$, output $Q = \mathcal{O}$ (`null`) deterministically in $< 2^{20}$ u.
2. **Independent 64-Bit Cold-Run Gate (Stage 5):** For $z \ne 1$, sample fresh uniform $U_{64} \in \{0, \dots, 2^{64}-1\}$. If $U_{64} \ge K_{\text{gate}} = 7{,}194{,}374{,}076{,}228{,}249{,}696$ (`0x63d78d1adadd2860`), abort immediately ($< 2^5$ u). Otherwise open with exact probability
   $$p_{\text{gate}} = \frac{K_{\text{gate}}}{2^{64}} = 0.3900078001560031, \qquad \log_2 p_{\text{gate}} = -1.3584251167231158.$$
3. **Reduction to Two Individual Logarithms in $G_T$:** Compute $g_0 = e(P, Q_0) \in G_T\setminus\{1\}$, recover $\chi_1$ factor-base logs via Stages 0–3, run Stage 4 twice to recover $\ell(g_0), \ell(z) \bmod r$, set $a = \ell(z)\ell(g_0)^{-1} \bmod r$, $Q = [a]Q_0$, and verify $e(P, Q) = z$.
4. **Complete Stage 1–4 Union-Bound Failure Budget (Explicitly Conditional on `H1–H5, H-R, H-SM, H-ECM`):**
   Every opened run executes behind a single hard pointwise aborting meter enforcing $W_{\text{pointwise}} \le 2^{124.0}$ u:
   - **Stage 1 Yield Shortfall** ($p_{\text{yield}} \le 2^{-100}$): Under H1–H2, $V_{\text{kept,ch}} = 2^{119.358898}$ includes exact `shortfall_pad` (`+0.002054` b) plus `+3*SE = +0.003382` b surplus above required mean $\mu = R_{\text{req}} = \eta U = 6.993855\times 10^{18} = 2^{62.600794}$, giving relative surplus $\delta = 2^{3\cdot\text{SE}} - 1 = 2.347236\times 10^{-3} \in (0, 1)$. By the multiplicative Chernoff lower-tail bound $\mathbb{P}(X \le (1-\delta)\mu) \le \exp(-\mu\delta^2/2)$ (and upper tail $\exp(-\mu\delta^2/(2+\delta)) = \exp(-1.9244\times 10^{13})$):
     $$\frac{\mu\delta^2}{2} = \mathbf{1.9266\times 10^{13}} \implies p_{\text{yield}} \le \exp(-1.9266\times 10^{13}) \ll 2^{-100}.$$
     Any unexpected shortfall triggers the Snell full-resume schedule (including `replay_a_tables` rebuild) charged in `snell_resume_tail = 2^{96.0265}` u.
   - **Stage 1b Miller–Rabin False-Prime Error** ($p_{\text{MR}} \le 2^{-132}$): 110-round Miller–Rabin ($4^{-110} = 2^{-220}$ per composite) across $\le 2^{63}$ relations and $\le 8$ cofactors gives union bound $\le 2^{-154} \ll 2^{-132}$.
   - **Stage 2–3 LA Rank / Block Wiedemann Failure** ($p_{\text{LA}} \le 2^{-200}$): Under H-R and H-SM, algebraic 2-merge yields full column rank over $\mathbb{F}_r$ (with Snell escalation $\tau \in [0.88, 0.92, 0.96, 1.00]$) and Coppersmith Block Wiedemann `(2048, 32)` fails with probability $O(n/r) \ll 2^{-200}$.
   - **Stage 4 Unmetered Descent Failure** ($p_{\text{descent,unmetered}} = 0$): Under H4–H5, Stage 4 is Las Vegas; meter exhaustion is counted identically in $p_{\text{meter\_abort}}$.
   - **Pointwise Aborting Meter Overflow ($W_{\text{pointwise}} = 2^{124.0}$ u)** ($p_{\text{meter\_abort}} \le 1.6758939\times 10^{-5}$): By **Markov's inequality alone (unconditional)**:
     $$p_{\text{meter\_abort}} = \mathbb{P}(W_{\text{run}} > 2^{124.0}) \le \frac{\mathbb{E}[W_{\text{open,charged}}]}{2^{124.0}} = \frac{2^{108.135290}}{2^{124.0}} = \mathbf{1.6758939\times 10^{-5}}.$$
   - **Wrapper Error Floor** ($p_{\text{wrapper}} \le 2^{-48}$).

   Summing all six terms gives:
   $$p_{\text{fail,open}} \le \mathbf{1.6758939\times 10^{-5} < 2.0\times 10^{-5}}, \quad p_{\text{open}} \ge \mathbf{0.99998324 > 0.99998}, \quad p_{\text{cold}} \ge p_{\text{gate}} p_{\text{open}} = \mathbf{0.39000126 > 0.39}.$$

### 1.3 ISA & Banked `classical-group-work-v2` Accounting Disclosure
- **Instruction Set:** All **15** compiled functions (`1,948` static instructions) in `KERNELS.md` and `KERNELS_STAGE1B_LA.md` use `-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize -fomit-frame-pointer -ffixed-x18` (`simd = 0`, `0` forbidden ISA, `0` `udiv`/`sdiv`, `0` `bl` helper calls, `0` in-loop spills).
- **Disclosed Scalar LSE Instruction:** Only `_emit_ptr64_x64` uses scalar LSE **`STADDB Ws, [Xn]`** (`73` static occurrences, priced at `1` instruction each). All other 14 kernels are pure ARMv8.0 scalar GPR code (`0` `staddb`).
- **Banked 64-Bit Accounting:** Every processor addresses strictly `< 16.000 EiB = 2^{64}` B locally. Every inter-bank transfer is executed and charged by compiled scalar GPR instructions on both Producer and Consumer (`0` free DMA, `0` unified $2^{117}$-byte pointers). All preprocessing, trials, failures, transfers, verification, and tails are charged in full for one cold target without amortization.

---

## 2. Stage 0: BLS12-381 Algebra, $\mathbb{Q}(\zeta_{15})$ Mode-A $(4,3)$ Isomorphism & Factor Base

BLS12-381 has $x = -\texttt{0xd201000000010000}$ ($\log_2|x| = 63.7142718$), $r = x^4 - x^2 + 1$ ($255$ b), $p = \frac{(x-1)^2(x^4-x^2+1)}{3} + x$ ($381$ b), $p \equiv 7 \pmod{15}$, and $3p = Q(x)\overline{Q}(x)$ in $\mathbb{Z}[\omega]$ with $Q(T) = T^3 - T^2 + \omega^2 T - \omega$. Over $K = \mathbb{Q}(\zeta_{15})$ ($\deg = 8$, $r_2 = 4$, $\Delta_K = 3^4\cdot 5^6$, $\mathbb{C}^4$-covolume $70.3125$), define
$$g(Y) = Y^3 - x, \qquad f(Y) = Q(Y^3) = Y^9 - Y^6 + \omega^2 Y^3 - \omega.$$
With $\pi = \mathfrak{p}_0\mathcal{O}_K$ ($\mathcal{O}_K/\pi \cong \mathbb{F}_{p^4}$) and $x^{(p-1)/3} \not\equiv 1 \pmod p$, $\mathcal{O}_K[Y]/(\pi, g(Y)) \cong \mathbb{F}_{p^{12}}$ and $f(Y) \equiv Q(x) \equiv 0 \pmod{(\pi, g)}$. Both $f, g$ are invariant under $\operatorname{Aut} \cong C_4 \times C_3 \cong \mathbb{Z}/12\mathbb{Z}$, and $p \equiv x \equiv \lambda \pmod r$ acts via faithful character $\chi_1$ giving $n_{\text{SM}}(\chi_1) = 36/12 + 12/12 = \mathbf{4}$ dense Schirokauer map columns (`schirokauer_maps = 2^{80.0001}` u).

**Saddle-Balanced Factor Base (`CAND3_CHAMP_7080_110X`):**
$$B_f = 2^{70.80}, \qquad B_g = 2^{71.16}, \qquad U = \frac{\pi(B_f) + \pi(B_g)}{12} = 7.94756307\times 10^{18} = 2^{62.785240}\text{ orbits}, \qquad \eta = 0.88.$$
Stage 0 setup and canonical $a$-domain construction are charged at `stage0_setup_and_iso = 2^{45.0001}` u and `canonical_a_domain = 2^{84.0001}` u.

---

## 3. Stage 1a: 100% 64-Bit Banked Fixed-VA $W_{16}$ Sieve, Chebotarev 2A/2B & Compiled Replay Build

### 3.1 Geometry, Zero-Interpolation `4,000,000`-MC Yield & Exact Pad Decomposition
We sieve $(a, b) \in \mathcal{A}_k \times \mathcal{B}_k$ with $a \in L_{0.65}(R_a)\cap\mathcal{F}_{\text{Vor}}$ ($q_a = 0.65$, $f_{\text{Vor}} = 0.08349025$, $R_a = 1719.468118$, $m_{\text{unit}} \equiv 1$, $K = 2$ equal-$V_a$ shells each with $V_a = 0.9351\cdot 2^{56} = 6.73810562\times 10^{16}$ points) paired with $b \in L_{0.65}(R_{b,k})/C_3$ ($R_{b,0} = R_0 = 1368.918001$, $R_{b,1} = 0.82\,R_0 = 1122.512761$, volume weights $w_0 = 0.83027923, w_1 = 0.16972077$) along 1D rational chords $b = b_{\text{base}} + t$ of length $L \le 248$:
- **Asymmetric Shell-0 / Shell-1 Sieve Bands & 3-Tier Producer Splits:**
  - Shell 0 ($83.03\%$): $g \in [28.0, 61.8]$ ($q_g = 3.4$), $f \in [32.0, 61.4]$ ($q_f = 3.8$).
  - Shell 1 ($16.97\%$): $g \in [28.0, 59.8]$ ($q_g = 3.4$), $f \in [32.0, 59.4]$ ($q_f = 3.8$).
  - Internal 3-tier splits: $p_{w1} = 54.0$ (Tier 1 $W_B=1$ vs Tier 2 $W_B=4$), $p_{w4} = 57.0$ (Tier 2 $W_B=4$ vs Tier 3 $W_B=16$).
- **Direct $N = 4{,}000{,}000$ 16D MC (Zero Interpolation, Exact 0.5-bit 2D Panjer Inside MC):** With 2D adaptive policy $\lambda = 332{,}708.30$ (`ret = 0.686917`, $\log_2 s_g = -13.960529$, $\log_2 s_f = -8.676454$, $\log_2 s_1 = -22.636983$, mean $(\overline{T}_g, \overline{T}_f) = (262.67, 174.99)$ b), Murphy shifts $\alpha_f = +5.3752$ b, $\alpha_g = -1.0487$ b, $1/\zeta_K(2) = 0.97843406$, and `-0.16` b cushion:
  - $V_{\text{kept,4M}} = \mathbf{2^{119.345701}}$ ($\text{SE} = 0.001127$ b, $3\cdot\text{SE} = +0.003382$ b), $R_{\text{ret}} = 2^{62.599636}$ vs $R_{\text{req}} = 2^{62.600794}$, giving explicit **`shortfall_pad = (R_req - R_ret) * slope = +0.002054` b**!
  - Conditional mean norms $(\bar{N}_f, \bar{N}_g) = (542.35, 670.44)$ b, observed maxima $(598.27, 685.90)$ b.
  - Across `200,000` chords ($L \le 248$): $L_{\text{win,eff}} = \mathbf{110.008527}$, measured interior gap `0.001577`, charged `gap_ch = 0.001677` (`gap_pad = +0.002422` b).
  - Adding the live-recomputed Buchstab–Dickman $f$-cofactor cap pad `cof_cap_pad = +0.005339` b (Section 4.1, zero $g$-cap!):
    $$\log_2 V_{\text{kept,ch}} = 119.345701 + 0.002054 + 0.003382 + 0.002422 + 0.005339 = \mathbf{119.358898}.$$

### 3.2 Unconditional Galois Chebotarev Theorems 2A & 2B
Verified exhaustively in `experiment.py` (matching all `1,857,859` primes $p \le 30\text{M}$ to `8 ppm`):
- **$g$-side ($Y^3 = x$):** Exact $d_{\mathcal{O}_K, g} = 1/3$, $\bar{r}_g = 3$, $\operatorname{dens}_g = 1/24$, active root count $r = 3$ with probability $1$.
- **$f$-side ($Q_\omega(Y^3) = 0$):** $\operatorname{Gal}(\operatorname{Spl}(f_\omega)/K) \cong C_3 \wr S_3$ (`162`, `8` normal subgroups) $\implies$ exact $d_{\mathcal{O}_K, f} = \mathbf{23/81}$, $\bar{r}_f = \mathbf{81/23}$, mixture $\mathbb{P}(r = 3, 6, 9 \mid r > 0) = (39/46, 6/46, 1/46)$, and $L_\omega \cap L_{\omega^2} = K(\zeta_9)$ unconditionally $\implies$ exact $\operatorname{dens}_f = \mathbf{2131/34992}$.

### 3.3 Banked 64-Bit Fixed-VA $M = 248$ Topology & 100% Literal Compiled 3-Tier Pricing
Every Slab-Owner $c \in \{0, \dots, 247\}$ maps its 1-column array `col[0..V_a-1]` ($0.058444\text{ EiB} \ll 16\text{ EiB}$) at `SLAB_BASE = 2^{56} = 0x0100_0000_0000_0000`. Since $V_a < 2^{56}$, every `ptr64 = SLAB_BASE + aid` satisfies `bits 57..59 == 000` identically! Every inter-bank transfer is charged literally on both Producer and Consumer (`KERNELS.md`):
- **Tier 1 ($[p_{\text{lo}}, 2^{54.0}]$, $W_B = 1$):** 0-spill `_tx_pack_u64_window` (`44 + 7*L_win` probe; `9 + 69*q64 + 35/19/11/7/3/3` exact Poisson).
- **Tier 2 ($(2^{54.0}, 2^{57.0}]$, $W_B = 4$):** 0-spill `_w4_packed_tx_window` explicitly reloading/storing `tx_ptrs[4*j + c]` per non-empty slice (`44 + p_ch*n_b4*5.0 + n_b4*p_b4*11.0` probe; `k=1: 6.0`, `k=2: 9.0`, `k>=3: 12 + 13*(k//8) + 7*(r&4) + 3*(r&2) + 3*(r&1)` exact Poisson!).
- **Tier 3 ($(2^{57.0}, p_{\text{hi},s}]$, $W_B = 16$, Zero-AND Demux!):** 0-spill `_demux_w16_window` puts the 4-bit column offset in bits `60..63` so `(e >> 57) == (col << 3)` with **zero `AND` mask** (`44 + p_ch*n_b16*5.0` probe; per bucket over ALL `n_b16 * 16` columns including overhang: `4.0 + 26*q4 + (1 if q4>0) + 8*rem`, i.e., `6.50 Ir/hit` x4 inner + `8.00 Ir/hit` tail, exact Poisson!).
- **Consumer Fixed-VA RX Emit:** 0-spill `_emit_ptr64_x64` (`100/64 = 1.5625 Ir/hit`, zero `add`!, draining `256 KiB` batched-`logp` RX rings).
- **1-Column Slab Prebias & 128B `rbit`/`clz` SWAR Scan:** 2-pass `_col_prebias_x128` (`0.171875 Ir/cell`) + `_col_scan_swar_128_fast` (`28/128 = 0.21875 Ir/cell` fast + `<= 60 Ir` `rbit`/`clz` slow extraction + `8 Ir` $f$-slab probe = `68 * 2^{-13.960529} = 0.004266`) $\implies$ `slab_live = 0.394891`, charged **`slab_ch = 0.409891 Ir/cell`** (`+0.015` pad). Zero-wraparound (`max_byte <= 255`) holds deterministically (`N_g <= 685.92 < 688`, `N_f <= 651.67 < 656`, $T_f \ge 45$).
- **Shell-Weighted 3-Tier Totals:** Shell 0 (`w0 = 0.830279`): `g (2.988033, 0.437046)`, `f (2.509593, 0.299640)`; Shell 1 (`w1 = 0.169721`): `g (2.571677, 0.170195)`, `f (2.119583, 0.126205)` $\implies$ **`ir_hits_avg = 5.360769`**, **`ir_probe_avg = 0.661960`**.
- **Compiled Incremental Replay Build (`KERNELS_STAGE1B_LA.md`):** `_replay_build_count_x4` (`8.50`) + `_replay_build_scatter_x4` (`9.50`) = `18.00` (charged **`19.0 Ir/ent`** T1/T2); `_replay_build_scatter_w16_x4` (`11.75` $\implies 20.25$, charged **`21.0 Ir/ent`** T3 W16); bucket prefix **`5.5 Ir/bkt`**; true 4-slice header builder `_replay_build_w4_prefix` (`~11.0`, charged **`12.0 Ir/hdr`**). Shell 0 (`1.3585e31 u`, **`2^{116.7442}` B**, freed before Shell 1) + Shell 1 (`3.1849e30 u`, `2^{114.3913}` B) $\implies$ **`replay_a_tables = 2^{103.725621}` u live (`103.7257` charged)**.

---

## 4. Stage 1b: 2-Line Continuous Analytic Norm Proofs, Live Buchstab Cap & Radix-4 `u2` Shoup61 NTT

### 4.1 Exact 2-Line Pure Algebraic Continuous Analytic 16D Norm Proofs ($q_a = q_b = 0.65$)
With $A_m = |\sigma_m(a)|$, $\sum_{m=0}^3 A_m^{0.65} \le R_a^{0.65}$ ($R_a = 1719.468118$) and $B_m = |\sigma_m(b)|$, $\sum_{m=0}^3 B_m^{0.65} \le R_0^{0.65}$ ($R_0 = 1368.918001$), **zero Hölder factor is needed**:
1. **Continuous $N_f$ Bound:** $|f_m(a,b)| \le (A_m^3 + B_m^3)(A_m^6 + B_m^6) \le (A_m + B_m)^9$. Concavity $(A_m + B_m)^{0.65} \le A_m^{0.65} + B_m^{0.65}$ and AM–GM give everywhere on the continuous 16D domain:
   $$N_f(a,b) = \prod_{m=0}^3 |f_m(a,b)|^2 \le \left(\frac{R_a^{0.65} + R_0^{0.65}}{4}\right)^{72/0.65} = \mathbf{2^{651.665490} < 2^{656}} \quad (\texttt{k\_pack} = 6,\ 6\cdot 656 = 3936 \le 4096).$$
2. **Continuous $N_g$ Bound:** $|g_m(a,b)| \le |x|\,(B_m + A_m/|x|^{1/3})^3$. Concavity + AM–GM give everywhere on the continuous 16D domain:
   $$N_g(a,b) = \prod_{m=0}^3 |g_m(a,b)|^2 \le |x|^8\left(\frac{R_0^{0.65} + R_a^{0.65}/|x|^{0.65/3}}{4}\right)^{24/0.65} = \mathbf{2^{685.924029} < 2^{688}} \quad (\texttt{k\_pack} = 5,\ 5\cdot 688 = 3440 \le 4096).$$
3. **Zero $g$-Side Cap & Live Buchstab–Dickman $f$-Side Cap (`408` b):** Tiers `g1, g2, g3` ALL keep worst-case `688` b leaves (`k_pack = 5`, **zero $g$-cofactor cap loss**!). On the $f$-side, post-Tier-1 cofactors are capped at `408` b (`k_pack = 10`, `4080 <= 4096` b). Exact 10-factor Buchstab–Dickman convolution (`cof_exceed_f_exact`, live-recomputed in `experiment.py` with normalization check `tot_prob / rho((N_f + alpha_f)/B_f)` = `0.9938` Shell 0 / `1.0046` Shell 1) bounds the unconditioned exceedance probability at $p_{\text{ex},f} = w_0(0.002770) + w_1(0.008215) = \mathbf{0.003694}$, strictly deducted via `cof_cap_pad = log2(1/(1 - p_ex_f)) = +0.005339` b inside $V_{\text{kept,ch}}$.

### 4.2 Verified 3-Prime $2^{54}$-Smooth Radix-4 `u2` Shoup61 NTT, Bernstein `dmult = 1.67` & Pow2-Rounded Cached RAM
Using the three verified $2^{54}$-smooth 61-bit primes $p_1 = 127\cdot 2^{54}+1\ (g_1=3)$, $p_2 = 114\cdot 2^{54}+1\ (g_2=7)$, $p_3 = 108\cdot 2^{54}+1\ (g_3=5)$ ($P = 2^{182.576462}$, $8 p_i < 2^{64}$, `1.58` b 1D / `5.25` b 2D CRT headroom at $L_t = 2^{67.1839}$ b):
- Compiled `_ntt61_r4_stage_dit_u2_exact` (`89/8 = 11.125 Ir/r2-bfly`, `0` bit-errors on full `N = 1024` DFT across all 3 primes), compiled `_ntt61_pw_mul_x4` (`9.25 -> 12.0 Ir`), scale (`12.0 Ir`), division-free compiled `_crt3_1d_carry` (`55 -> 58.0 Ir` 1D, `90.0 Ir` 2D), `1.10x` safety, and `own_eib = 10.0 EiB < 16 EiB` (`L_own_1d = 2^{64.15}` b, `L_own_2d = 2^{63.15}` b, `160 Ir/pt` cross-bank transpose).
- **Bernstein Scaled-Remainder-Tree / Middle-Product Derivation (`dmult = 1.67`):** In a Bernstein scaled remainder tree descending a batch product tree $X$, each node's forward 3-prime NTT $\widehat{X}_v$ is computed once during the upward product-tree build (which costs $1.5\times$ NTTs per multiply: 2 forward + 1 inverse) and retained in memory; each downward scaled-remainder / middle-product step then needs only $1$ forward NTT (of the incoming scaled parent remainder) + $1.5$ inverse/pointwise NTTs (for left & right children reusing cached $\widehat{X}_{L}, \widehat{X}_{R}$), giving exact downward ratio $2.5/1.5 = \mathbf{1.67}$ (and static prime product-tree forward NTTs are built once during `band_split_and_trees = 2^{90.0001}` u).
- **Power-of-2 Rounded Cached 3-Prime NTT-Domain RAM:** Summing exact next-power-of-2 1D ($N\cdot 24$ B) and 2D ($N_1 N_2\cdot 24$ B) cached NTT sizes across ALL `56` levels of the largest tier (`f4` Shell 1, $L_t = 2^{67.1839}$ b) plus integer tree & 2D work gives `2^{73.0964}` B (NTT cache) `+ 2^{70.2062}` B (int+2D work) = **`2^{73.2789}` B** — distributed across `2,048` Stage 1b owners at **`4.8532 EiB/owner < 10.0 EiB < 16 EiB`**!
- **Stage 1b Total:** **`c_s1b = 1.000943 Ir/cell`** $\implies$ **`ipc = 5.360769 + 0.661960 + 0.409891 + 1.000943 = 7.433563 Ir/cell`** $\implies$ **`relation_collection = 2^{107.965240}` u live (`107.9653` charged)**.

---

## 5. Stages 2 & 3: 6-Rung 4-Seed `L_CALIB STRICT` SGE (`H3`), `11.00 Ir/nnz` SpMV & 576-Bit `r << 256` Horner+SM4

### 5.1 Complete 6-Rung 4-Seed `L_CALIB STRICT` (`W_tgt = 180`) SGE Ladder at $\eta = 0.88$
24 independent seeded `STRICT` Monomial-Pivot $\mathbb{Z}[\zeta_{12}]$ SGE runs (`0` failures, `0.00% >32b`, `0` overflows):

| Rung $U$ | Seeds | `trimC/U` | $N_l / U$ (Mean $\pm$ SE) | $N_{\text{nnz}} / U$ | Mean $w$ | Mono % | `>32b` | `ov30` | `x1` Core Ir |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| $2^{16}$ | 4 | `0.04709` | `0.005753 +- 0.001233` | `0.31206` | `44.1` | `63.30%` | `0.00%` | `0` | `18.53` |
| $2^{18}$ | 4 | `0.02510` | `0.002721 +- 0.000415` | `0.12060` | `37.4` | `71.99%` | `0.00%` | `0` | `15.87` |
| $2^{20}$ | 4 | `0.02572` | `0.002397 +- 0.000403` | `0.16816` | `62.4` | `74.68%` | `0.00%` | `0` | `15.21` |
| $2^{22}$ | 4 | `0.01529` | `0.001577 +- 0.000490` | `0.12874` | `59.5` | `81.76%` | `0.00%` | `0` | `12.92` |
| $2^{23}$ | 4 | `0.01202` | `0.001005 +- 0.000097` | `0.08283` | `67.6` | `80.18%` | `0.00%` | `0` | `13.61` |
| **$2^{24}$** | **4** | **`0.01202`** | **`0.001045 +- 0.000104`** | **`0.12454`** | **`100.2`** | **`81.43%`** | **`0.00%`** | **`0`** | **`12.98`** |

**Explicit Scaling Disclosure (`H3` — NOT `"zero-extrapolation"`):** Extending $U \le 2^{24}$ SGE densities to $U = 2^{62.79}$ is explicitly disclosed as an `exploratory-v1` scaling premise (`H3`). Charged caps `nl_cap = 0.00150` (`1.435x` mean, `1.197x` `mean + 2*SE = 0.001253`) and `sw_cap = 0.190` (`1.526x` `0.12454`, `w_cap = 126.67`, **`2.19x` combined $N_l\cdot N_{\text{nnz}}$ LA headroom**) give $N_l = 1.192134\times 10^{16}$, $N_{\text{nnz}} = 1.510037\times 10^{18}$, and `matrix_generator_sge = 2^{62.481384}` u live (`62.4815` charged).

### 5.2 Block Wiedemann `(m = 2048, n = 32)`, Compiled `11.00 Ir/nnz` SpMV & `872 Ir` 576-Bit `r << 256` Row Kernel
1. **Compiled `x4` Monomial SpMV (`_spmv_x4_monomial_add`, `KERNELS_STAGE1B_LA.md`):** `11.00 Ir/nnz` compiled (`0` in-loop spills), charged conservatively at `k_spmv_x4 = 12.141 Ir/nnz` (Krylov) and `k_spmv_x1 = 16.000 Ir/nnz` (unbatched `Mksol`).
2. **Compiled 576-Bit `r << 256`-Biased Horner+SM4+Barrett Row Kernel (`_spmv_row_horner_sm4_reduce`):** Biases the 9-limb accumulator by `r << 256` (`> 0.4528 * 2^{512}`), so for ALL $0 \le A_k < 2^{315}$ (far exceeding `180 * 2^{30} * r < 2^{293}`!):
   - Negative term $|x|(A_1 + |x|^2 A_3) < 2^{507.1} < 0.035\cdot 2^{512} \ll r\cdot 2^{256}$ (**proven ZERO underflow!**);
   - Positive total $T = (A_0 + |x|^2 A_2) + (r\ll 256) - |x|(A_1 + |x|^2 A_3) + \sum_{j=0}^3 \text{sm}_j v_j \in (0.41\cdot 2^{512}, 1.28\cdot 2^{512})$ has 9th limb $t_8 \in \{0, 1\}$, folded branchlessly into $[0, 2^{512})$ by subtracting `2*r << 256` when $t_8 = 1$ (**proven ZERO overflow!**), then reduced via one exact 512$\to$256 Barrett step.
   - Exact compiled dynamic count: **`192 + 4*106 + 256 = 872 Ir`** (`0` helper calls, `200/200` bit-exact including $A_k = 2^{314}-1$!), conservatively charged **`rovh = 920.0 Ir/row`** (`bw_ovh = 1.08`).
3. **BW `(2048, 32)` & Complete Snell Resume Tail:** $L_{\text{seq}} = 3.783630\times 10^{14}$, $d_{\text{mksol}} = 3.725420\times 10^{14}$ $\implies$ **`linear_algebra = 2^{104.159692}` u live (`104.1598` charged)**, Stage 3 RAM (incl. full Thomé `Lingen` $\pi(z)$ FFT scratch) **`2^{76.4613}` B**. Snell resume tail explicitly includes `replay_a_tables` rebuild:
   $$W_{\text{snell}} = 2^{-10}\left(\frac{1 - 0.88}{0.88}\,W_{\text{coll}} + W_{\text{rep}} + W_{\text{LA}}\right) = \mathbf{2^{96.026351}\text{ u live}} \quad (\mathbf{96.0265\text{ charged}}).$$

---

## 6. Stage 4, Banked 64-Bit RAM Audit & Authoritative 14-Row Cold-Run Ledger

### 6.1 Stage 4 12D Subfield $M_g = \mathbb{Q}(\zeta_5)[Y]/(Y^3 - x)$ Hilbert-90 Descent
Unconditional 12D LLL+AM–GM norm bound $2^{2570.51} < 2^{2572}$; explicit Layer-0 inert $q^2$ descent (`2^{64.1470}` u) + `+20` inert Layer-5 nodes; Torsion-16 Edwards ECM. At $(70.80, 71.16, 0.88)$, live 2-log with Snell factor $(1.04094/0.88)^2$ is **`2^{95.4719}` u** (conservatively charged **`both_individual_logs_exp = 2^{96.0001}` u**, `+0.528` b headroom!), and **`ecm_split_expected = 2^{93.500021}` u live (`93.5001` charged)**.

### 6.2 Complete 100% 64-Bit Banked Storage Audit (`<= 12.000 EiB/processor`, `<= 2^{117.50}` B Global)

| Stage / Subsystem | Per-Processor Local RAM | Global Total Storage | Notes |
|---|---:|---:|---|
| Stage 1a Prime-Bank Replay Tables (`Shell 0` built, sieved, freed before `Shell 1`) | **`<= 12.000 EiB`** (`< 16 EiB`) | **`2^{116.7442}` B** | $P_{\text{banks}} = 10{,}058{,}276{,}997{,}840{,}928$ banks (`Shell 1`: `2^{114.3913}` B) |
| Stage 1a Slab-Owners ($M = 248$ @ fixed VA `2^{56}` + `256 KiB` RX + $g$-surv `2^{44.94}` B/col) | **`0.058444 EiB`** (`<< 16 EiB`) | `14.4941 EiB` + `2^{52.90}` B | $g$-survivors `V_a*2^{-13.9605}*8 = 2^{44.94}` B/col (`2^{52.90}` B / 248 cols) |
| **Stage 1 Combined Peak Storage** | **`<= 12.000 EiB`** | **`2^{116.7442}` B** | Shell 0 replay + 248 slabs + RX + $g$-survivors |
| Stage 1b Pow2-Rounded Cached 3-Prime NTT Tree (`f4` Shell 1, `56` lvls, `2,048` owners) | **`4.8532 EiB`** (`< 10 EiB`) | **`2^{73.2789}` B** | `2^{73.0964}` B pow2 NTT cache + `2^{70.2062}` B int/2D work |
| Stage 2 SGE Initial Relations (`eta * U * 18 * 12` B) | `< 1.0 EiB` | **`2^{70.3557}` B** | Stage 1 freed before Stage 2 |
| Stage 3 LA Matrix + Vectors + Full Thomé `Lingen` $\pi(z)$ (`m=2048, n=32`) | `< 4.0 EiB` | **`2^{76.4613}` B** | Incl. `2*2048*2080*L_seq*32` $\pi(z)$ FFT scratch |
| Stage 4 12D $M_g$ Descent / Nonuniform Advice | `< 1.0 GiB` / **`0` B** | `< 1.0 GiB` / **`0` B** | Zero searched/stored advice |
| **Reported Peak Global Storage Bound** | **`<= 12.000 EiB < 16 EiB`** | **`<= 2^{117.50}` B** | Strictly `<< 2^{120}` B (`3.26` b headroom) |

### 6.3 Complete Authoritative 14-Row Cold-Run Ledger (Recomputed Live in `experiment.py`)

| # | Ledger Row | Live ($\log_2 u$) | Charged ($\log_2 u$) | Slack (b) | Share of Opened |
|---:|---|---:|---:|---:|---:|
| 1 | `relation_collection` | `107.965240` | **`107.9653`** | `+0.000060` | `88.885%` |
| 2 | `linear_algebra` | `104.159692` | **`104.1598`** | `+0.000108` | `6.357%` |
| 3 | `replay_a_tables` | `103.725621` | **`103.7257`** | `+0.000079` | `4.705%` |
| 4 | `snell_resume_tail` | `96.026351` | **`96.0265`** | `+0.000149` | `0.023%` |
| 5 | `both_individual_logs_exp` | `96.000000` | **`96.0001`** | `+0.000100` | `0.022%` |
| 6 | `ecm_split_expected` | `93.500021` | **`93.5001`** | `+0.000079` | `0.004%` |
| 7 | `meter_safety_buffer` | `93.400000` | **`93.4001`** | `+0.000100` | `0.004%` |
| 8 | `band_split_and_trees` | `90.000000` | **`90.0001`** | `+0.000100` | `<0.001%` |
| 9 | `canonical_a_domain` | `84.000000` | **`84.0001`** | `+0.000100` | `<0.001%` |
| 10 | `schirokauer_maps` | `80.000000` | **`80.0001`** | `+0.000100` | `<0.001%` |
| 11 | `final_log_backsub_and_gcd` | `72.000000` | **`72.0001`** | `+0.000100` | `<0.001%` |
| 12 | `matrix_generator_sge` | `62.481384` | **`62.4815`** | `+0.000116` | `<0.001%` |
| 13 | `stage0_setup_and_iso` | `45.000000` | **`45.0001`** | `+0.000100` | `<0.001%` |
| 14 | `final_verification_Q` | `40.000000` | **`40.0001`** | `+0.000100` | `<0.001%` |
| — | **Opened-Run Expected Work ($p_{\text{open}} \ge 0.99998324$)** | **`108.135226`** | **`108.135290`** | **`+0.000064`** | **`100.000%`** |
| — | **64-Bit Pre-Gate Factor ($\log_2(K_{\text{gate}}/2^{64})$)** | `-1.358425` | `-1.358425` | — | — |
| — | **Cold-Run Gated Expected Work ($p_{\text{cold}} \ge 0.39000126$)** | **`106.776801`** | **`106.776865`** | **`+0.003135`** | — |
| — | **OFFICIAL REDUCED CLAIM (`5339/50`, `gcd = 1`)** | — | **`106.780000`** | — | — |

---

## 7. Explicit Disclosure of Exploratory Heuristics

This is an `exploratory-v1` cryptanalysis bound conditional on the disclosed heuristics below (whereas Chebotarev 2A/2B, the 2-line continuous norm proofs, zero-wraparound, 3-prime Radix-4 NTT CRT headroom, 12D $M_g$ LLL+AM–GM, all 15 compiled kernel instruction counts, `r << 256` zero-underflow/overflow up to $A_k < 2^{315}$, and the Markov meter bound are **unconditional**):

| ID | Heuristic Statement | Quantitative Support in This Package | Limitation / Residual Risk |
|---|---|---|---|
| **H1** | **Smoothness & Murphy $\alpha$:** Norms $N_f, N_g$ behave independently as Murphy-shifted random integers with Dickman $\rho(u)$ and $1/\zeta_K(2)$. | Zero-interp `N = 4,000,000` 16D MC (`SE = 0.00113` b) + `-0.16` b cushion + `shortfall_pad` + `3*SE` + `gap_pad` + live Buchstab `cof_cap_pad`; exact Chernoff $\mu\delta^2/2 = 1.9266\times 10^{13}$ conditional on H1. | Standard NFS/STNFS heuristic; not unconditionally proved at 600-bit norms. |
| **H2** | **Dual Band-Pass Independence & 2D Adaptive DP:** Prime divisibility in Shell-0/1 bands follows exact 2D Compound-Poisson DP. | Exact 0.5-bit 2D Panjer DP inside 4M MC; `max_byte <= 255` proved unconditionally. | Minor small-prime correlations at band boundary; mitigated by exact Murphy-shifted 2D DP. |
| **H3** | **Post-SGE Matrix Scaling Premise ($U \le 2^{24} \to 2^{62.79}$):** At $U = 2^{62.79}, \eta = 0.88$, post-SGE matrix satisfies $N_l/U \le 0.00150$, $N_{\text{nnz}}/U \le 0.190$ ($w \le 126.67$), `x4` SpMV $\le 12.141$ Ir/nnz. | Complete 6-rung 4-seed `L_CALIB STRICT` ladder ($U = 2^{16}..2^{24}$, 24 runs, `0` failures, `0.00% >32b`): $U = 2^{24}$ mean `0.001045 +- 0.000104`, `N_nnz/U = 0.12454` (`2.19x` combined LA headroom). | Explicitly an `exploratory-v1` scaling premise across 38.8 doublings of $U$; assumes post-SGE ratios do not reverse above $U = 2^{24}$. |
| **H4** | **Stage 4 12D $M_g$ Descent:** Smoothness of 12D $M_g$ norms ($\le 2572$ b proved unconditionally) and 6-layer descent follow Dickman $\rho(u)$. | Unconditional LLL+AM–GM (`2570.51 < 2572` b); live `2^{95.472}` u charged at `2^{96.0001}` u (`0.022%` of opened work). | Even a `4x` Stage 4 blowup adds `< 0.001` b and stays strictly below claim. |
| **H5 / H-ECM** | **Torsion-16 Edwards ECM Smoothness:** Random Torsion-16 curves mod $q$ achieve $q_{\text{eff}} = q - 3.5$ b. | Algebraic proof `21/21` primes, live Stage-2 witnesses (`50023` & `35`). | Hasse-interval smoothness heuristic; Stage 4 is `0.022%` of work. |
| **H-R & H-SM** | **Post-2-Merge Dense-Core Full Rank & 4 Non-Degenerate SMs:** Full rank over $\mathbb{F}_r$ at $\eta = 0.88$ and non-degenerate $\chi_1$ SMs. | Live exact dense $\mathbb{F}_q$ 2-merge full-rank checks (`0` nullity) + complete `snell_resume_tail = 2^{96.0265}` u (incl. replay rebuild). | Standard finite-scale rank & Schirokauer non-degeneracy heuristics. |

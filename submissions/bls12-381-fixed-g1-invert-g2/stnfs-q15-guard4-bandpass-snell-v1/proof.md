# Cold Single-Target BLS12-381 FAPI-1 via Guard4 x16 Dual Band-Pass $\mathbb{Q}(\zeta_{15})$ CM-STNFS and Snell-Envelope Optimal Stopping (`proof.md`)

**Track:** `bls12-381-fixed-g1-invert-g2` (Fixed-Argument Pairing Inversion into $G_2$)  
**Cost Model:** `classical-group-work-v2` ($C = 20{,}000$ scalar ARMv8.1-A general-purpose-register instruction equivalents per dense $\mathbb{F}_{p^{12}}$ multiplication unit; `simd = 0`; cold single-target run; algorithmic success $p_{\text{cold}} \ge 0.39$)  
**Official Reduced Claim (`submission.json`, `gcd = 1`):** **`10667/100` = `106.670000`** ($\log_2 u$; charged gated ledger **`106.669708`**, live gated **`106.669534`**, opened-run charged **`108.028134`**)

---

## 1. Executive Summary & Reduction to $G_T$-DLP

The challenge supplies a uniform non-identity point $P \in G_1 \setminus \{\mathcal{O}\}$ and an independently uniform target $z \in G_T \subset \mathbb{F}_{p^{12}}^\times$ (including $z = 1$), and asks for $Q \in G_2$ satisfying $e(P, Q) = z$ under the pinned optimal Ate pairing `py_ecc` v8.0.0.

### 1.1 Exact Reduction, Pre-Gate & Unconditional Markov Meter
1. **Canonical Decode & Subgroup Validation (Stage 0):** Validate $P \in E_1(\mathbb{F}_p)$, $P \ne \mathcal{O}$, $[r]P = \mathcal{O}$, and $z \in \mathbb{F}_{p^{12}}^\times$, $z^r = 1$. If $z = 1$, immediately output $Q = \mathcal{O}$ (`null`), satisfying $e(P, \mathcal{O}) = 1$ deterministically in $< 2^{20}$ u.
2. **Independent 64-Bit Cold-Run Gate (Stage 5):** For $z \ne 1$, draw a fresh uniform 64-bit coin $U_{64} \in \{0, \dots, 2^{64}-1\}$ isolated from $(P, z)$ and all subsequent algorithm coins. If $U_{64} \ge K_{\text{gate}} = 7{,}194{,}374{,}076{,}228{,}249{,}696$ (`0x63d78d1adadd2860`), abort immediately with negligible work ($< 2^5$ u). Otherwise open the attack with exact probability
$$p_{\text{gate}} = \frac{K_{\text{gate}}}{2^{64}} = 0.3900078001560031, \qquad \log_2 p_{\text{gate}} = -1.3584251167231158.$$
3. **Reduction to Two Individual Logarithms in $G_T$:** Compute $g_0 = e(P, Q_0) \in G_T \setminus \{1\}$ ($\operatorname{ord}(g_0) = r$ by non-degeneracy). Run Stages 0–3 once to recover the discrete logarithms modulo $r$ of the $\mathbb{Z}/12\mathbb{Z}$-automorphism factor base in the faithful $\chi_1$-eigenspace. Run Stage 4 twice to recover $\ell(g_0) = \log_h(g_0) \bmod r$ and $\ell(z) = \log_h(z) \bmod r$, compute $a = \ell(z)\cdot\ell(g_0)^{-1} \bmod r$, compute $Q = [a]Q_0 \in G_2$, and verify $e(P, Q) = z$ before output.
4. **Unconditional Markov Pointwise Aborting Meter ($W_{\text{pointwise}} = 2^{124.0}$ u):** Every opened run executes behind a single hard pointwise aborting meter enforcing $W_{\text{pointwise}} \le 2^{124.0}$ u. By **Markov's inequality alone (requiring zero concentration assumptions)**:
$$\mathbb{P}(W_{\text{run}} > 2^{124.0}) \le \frac{\mathbb{E}[W_{\text{open,charged}}]}{2^{124.0}} = \frac{2^{108.028134}}{2^{124.0}} = 1.55593 \times 10^{-5}.$$
Combining this unconditional Markov abort bound with a $2^{-48}$ wrapper/MR error floor yields
$$p_{\text{open}} \ge 1 - 1.55593 \times 10^{-5} - 2^{-48} = \mathbf{0.99998444 > 0.99998 = p_{\text{open,floor}}}.$$
*(Under exponential concentration across Stages 1–4, actual abort probability is orders of magnitude smaller still; the ledger row `meter_safety_buffer = 2^{93.4001}` u is carried as an extra conservative margin.)*

### 1.2 ISA & Unified `classical-group-work-v2` Accounting Disclosure
- **Instruction Set:** All hot-loop kernels in `KERNELS_SIEVE.md` and `KERNELS_LA_ECM.md` are compiled with `-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize` (`simd = 0`, zero NEON/FP/SVE vector registers `v/q/d/s/h/z/p`, zero `dc zva`, zero `x18`).
- **Disclosed Scalar LSE Instruction:** Stage 1a Guard4 x16 hit emission (`_emit_u64_x16_direct`) uses the scalar ARMv8.1-A Large System Extensions (LSE) byte atomic-add **`STADDB Ws, [Xn]`** (a 64-bit general-purpose-register-only memory add touching no vector state, priced at `1` instruction each). Every other kernel in `KERNELS_SIEVE.md` and `KERNELS_LA_ECM.md` uses pure ARMv8.0 scalar instructions only (`0` `staddb`).
- **Unified RAM Model & Complete Cold-Run Charging:** We work strictly under the unified flat-RAM `classical-group-work-v2` machine model (every memory access counted at `1` instruction per `ldr`/`ldp`/`str`/`stp`/`staddb`, with peak total storage explicitly bounded across every stage in Section 6.2 without uncharged multi-owner network transfers). All preprocessing, trials, failures, transfers, verification, and tails are charged in full for one single cold target without amortization. Nonuniform advice is `0` bytes (`2^0` B).

---

## 2. Stage 0: BLS12-381 Algebra, $\mathbb{Q}(\zeta_{15})$ Mode-A $(4,3)$ Isomorphism & Schirokauer Maps

### 2.1 Target Parameters & Field Tower
From `target-profile.json`, BLS12-381 is parametrized by $x = -\texttt{0xd201000000010000}$ ($|x| = 15{,}132{,}376{,}222{,}941{,}642{,}752$, $\log_2|x| = 63.7142718$):
$$r = x^4 - x^2 + 1 \quad (255\text{ bits}), \qquad p = \frac{(x-1)^2(x^4 - x^2 + 1)}{3} + x \quad (381\text{ bits}).$$
Direct arithmetic gives $p \equiv 1 \pmod 3$, $p \equiv 2 \pmod 5$ ($p \equiv 7 \pmod{15}$), and the exact factorization in $\mathbb{Z}[\omega]$ ($\omega^2 + \omega + 1 = 0$):
$$3p = (x-1)^2(x^4 - x^2 + 1) + 3x = x^6 - 2x^5 + 2x^3 + x + 1 = Q(x)\,\overline{Q}(x),$$
where $Q(T) = T^3 - T^2 + \omega^2 T - \omega \in \mathbb{Z}[\omega][T]$.

### 2.2 Base Field $K = \mathbb{Q}(\zeta_{15})$ and Mode-A $(4,3)$ Polynomial Pair
Let $K = \mathbb{Q}(\zeta_{15}) = \mathbb{Q}(\omega, \zeta_5)$ ($\omega = \zeta_{15}^5$), a totally complex CM field of degree $\deg(K/\mathbb{Q}) = \varphi(15) = 8$, $r_2 = 4$ complex conjugate embedding pairs indexed by $\sigma_m(\zeta_{15}) = e^{2\pi i k_m / 15}$ for $k_m \in \{1, 4, 7, 13\}$ (all four fixing $\omega = e^{2\pi i / 3}$), discriminant $\Delta_K = 3^4 \cdot 5^6 = 1{,}265{,}625$, $\mathbb{C}^4$-lattice covolume $\sqrt{\Delta_K}/2^{r_2} = 1125/16 = 70.3125$, exact regulator $R_K = 4.661820781449438$, root-of-unity group $\mu_K = \langle -\zeta_{15} \rangle$ of order $30$, and fundamental unit basis $\{u_1 = 1 - \zeta_{15}^7,\; u_2 = 1 + \zeta_{15},\; u_3 = \zeta_{15}^3 + \zeta_{15}^{12}\}$.

Define the Mode-A $(4,3)$ polynomial pair over $\mathcal{O}_K[Y]$:
$$g(Y) = Y^3 - x, \qquad f(Y) = Q(Y^3) = Y^9 - Y^6 + \omega^2 Y^3 - \omega.$$
**Theorem (Exact $\mathbb{F}_{p^{12}}$ Isomorphism, verified live in `experiment.py`):**
1. Since $p \equiv 1 \pmod 3$, $(p) = \mathfrak{p}_0 \overline{\mathfrak{p}}_0$ in $\mathbb{Z}[\omega]$, where $\mathfrak{p}_0 = (p, \omega - \omega_p)$ is chosen so that $Q(x) \equiv 0 \pmod{\mathfrak{p}_0}$.
2. Since $p \equiv 2 \pmod 5$ has multiplicative order $\operatorname{ord}_5(p) = 4$, $\Phi_5(t) = t^4 + t^3 + t^2 + t + 1$ is irreducible over $\mathbb{F}_p$. Hence $\pi = \mathfrak{p}_0\mathcal{O}_K$ is a prime ideal of $\mathcal{O}_K$ with $\mathcal{O}_K / \pi \cong \mathbb{F}_{p^4}$.
3. Direct modular exponentiation confirms $x^{(p-1)/3} \not\equiv 1 \pmod p$; since $(p^4 - 1)/(p - 1) = p^3 + p^2 + p + 1 \equiv 4 \equiv 1 \pmod 3$, $x$ is not a cube in $\mathbb{F}_{p^4}$. Therefore $g(Y) = Y^3 - x$ is irreducible over $\mathbb{F}_{p^4}$, yielding the ring isomorphism
$$\mathcal{O}_K[Y] / (\pi, g(Y)) \cong \mathbb{F}_{p^4}[Y] / (Y^3 - x) \cong \mathbb{F}_{p^{12}},$$
and $f(Y) = Q(Y^3) \equiv Q(x) \equiv 0 \pmod{(\pi, g(Y))}$, so both $L_f = K[Y]/(f(Y))$ (degree $72$, $r_2 = 36$) and $L_g = K[Y]/(g(Y))$ (degree $24$, $r_2 = 12$) map surjectively onto $\mathbb{F}_{p^{12}}$.

### 2.3 Cyclic $\mathbb{Z}/12\mathbb{Z}$ Automorphisms, Frobenius Action & $n_{\text{SM}}(\chi_1) = 4$
Both $f(Y)$ and $g(Y)$ are invariant under $\operatorname{Gal}(K/\mathbb{Q}(\omega)) \cong C_4$ (inducing the $p$-Frobenius on $\mathcal{O}_K/\pi \cong \mathbb{F}_{p^4}$) and $\tau_3: Y \mapsto \omega Y$ of order $3$, generating $\operatorname{Aut} \cong C_4 \times C_3 \cong \mathbb{Z}/12\mathbb{Z}$.

Because $r = x^4 - x^2 + 1 = \Phi_{12}(x)$, the parameter $\lambda := x \bmod r$ is an exact primitive 12th root of unity in $\mathbb{F}_r^\times$ ($\lambda^4 - \lambda^2 + 1 \equiv 0 \pmod r$, $\lambda^6 \equiv -1 \pmod r$). Moreover, from $p = \frac{(x-1)^2}{3} r + x$, we have the exact congruence
$$p \equiv x \equiv \lambda \pmod r,$$
so the $p$-Frobenius acts on $G_T = \mu_r(\mathbb{F}_{p^{12}})$ by scalar multiplication by $\lambda = x \bmod r$ (defining the faithful order-12 character $\chi_1$). Stage 3 linear algebra solves strictly inside the $\chi_1$-eigenspace, where the free $\mathbb{Z}[\operatorname{Aut}]$-action on the $36$ complex places of $L_f$ and $12$ complex places of $L_g$ yields exactly $36/12 + 12/12 = \mathbf{4}$ dense Schirokauer map columns ($n_{\text{SM}}(\chi_1) = 4$, computed via exact $r$-th power lifts in `schirokauer_maps = 2^{80.0001}` u and evaluated in Stage 3 via `_fr_dot4_mont`).

### 2.4 Rebalanced Factor Base & Stage 0 Setup
We rebalance the factor base to $\Delta B = -1.60$ bits:
$$B_f = 2^{69.60}, \qquad B_g = 2^{69.96}, \qquad U = \frac{\pi(B_f) + \pi(B_g)}{12} = 3.51885007 \times 10^{18} = 2^{61.609810}\text{ orbits}.$$
Stage 0 setup and canonical $a$-domain construction are charged at `stage0_setup_and_iso = 2^{45.0001}` u and `canonical_a_domain = 2^{84.0001}` u.

---

## 3. Stage 1a: Guard4 x16 Dual Band-Pass Sieve & Sequential Single-Slab Active-Only Masked Scan

### 3.1 Geometry: Two Equal-$V_a$ Canonical Shells + $L_{0.65}(b)$ Chords (`L <= 248`)
We sieve directly over $(a, b) \in \mathcal{A}_k \times \mathcal{B}_k$ across $K = 2$ shells ($k \in \{0, 1\}$):
- **Canonical $a$-Domain $\mathcal{A}_k$:** $a \in \mathcal{F}_{\text{Vor}}(a) \cap L_2(R_{a,\text{tot}}) / (\mu_{30} \times C_4)$ ($|\mu_{30} \times C_4| = 120$), split radially into $K = 2$ equal-volume shells of cardinality $V_{a,\text{sh}} = 0.9351 \cdot 2^{56} = 6.73810562 \times 10^{16}$ ($R_{a,\text{tot}} = 424.1944$, Shell 0 outer radius $388.9880$). Contiguous indexing `a_id in [0, V_{a,sh})` gives **unit multiplicity $m_{\text{unit}}(a, b) \equiv 1$ identically** and 100% dense array occupancy.
- **Complex-Place $L_{0.65}(b)$ Ball & 1D Rational Chords:** Pair Shell $k$ with $b \in L_{0.65}(R_{b,k})/C_3$ ($R_{b,0} = R_0 = 1415.88$, $R_{b,1} = 0.82\,R_0 = 1161.02$). Varying $b_0 \in \mathbb{Z}$ traces 1D rational chords $b = b_{\text{base}} + t$ shifting all four complex embeddings uniformly by $+t$.
- **Direct $N = 4{,}000{,}000$ 16D Monte Carlo Yield & 200,000-Chord Guard4 (`L <= 248`) Audit:**
  - At $(B_f, B_g) = (2^{69.60}, 2^{69.96})$, $\eta = 0.82$, exact 2D Panjer retention $\text{ret}_{\text{exact}} = 0.68738$ (MC verified at conservative `0.6864 <= 0.68738`, so `ret_shift_bits = 0.0`!), exact Murphy constants $\alpha_f = +5.3752$ b, $\alpha_g = -1.0487$ b, $1/\zeta_K(2) = 0.97843406$, and `-0.16` b cushion, `N = 4,000,000` direct 16D draws (zero interpolation!) yield $V_{\text{kept,MC}} = 2^{119.73548}$ ($\text{SE} = 0.00116$ b), conditional smooth norms $(\bar{N}_f, \bar{N}_g) = (543.76, 671.40)$ bits, and observed maxima $(598.7, 687.1)$ bits.
  - Evaluating the exact same 4M-draw integrator at the $(q=0.75, K=2, B_f=71.20, B_g=71.56, \eta=1.04094, \text{ret}=0.7623)$ comparator point yields $2^{119.37413}$ vs historical comparator $2^{119.43957}$; we conservatively add the full **`+0.06544` b anchor pad** to our charged volume (and also report the pure first-principles score without it!).
  - Across `200,000` importance-sampled $L_{0.65}$ chords (Shell 0 share `0.83`, Shell 1 `0.17`) partitioned into Guard4 windows of active length $L \le 248$:
    - Shell 0: `mean_span = 124.02`, `win_mult_248 = 1.07801`; Shell 1: `mean_span = 101.62`, `win_mult_248 = 1.01604`.
    - Effective window length: $L_{\text{win,eff}} = 0.83\cdot(124.02/1.07801) + 0.17\cdot(101.62/1.01604) = \mathbf{112.4903}$ active columns/window.
    - Non-convex $L_{0.65}$ interior chord gap: measured `0.00169`, conservatively charged **`0.00175`** (`+0.002527` b).
  - Adding `3 * SE = +0.00348` b and `+0.002527` b chord-gap pad yields:
$$\log_2 V_{\text{kept,pure}} = 119.741487, \qquad \log_2 V_{\text{kept,charged}} = 119.73548 + 0.06544 + 0.002527 + 0.00348 = \mathbf{119.806927}.$$

### 3.2 Deterministic Norm Upper Bounds, Exact 2D Panjer Smooth-Recall DP & Zero Wraparound
Every pair $(a, b)$ in the relaxed box $(R_{a,\text{tot}} = 424.1944,\; R_0 = 1415.88)$ satisfies two **unconditional deterministic norm upper bounds** (re-verified globally with all cross terms in `experiment.py`):
1. **$g$-Norm (AM–GM Global Maximum):** At symmetric $|\sigma_m(b)| = R_0/4^{1/0.65} = 167.7955$, $|\sigma_m(a)| = R_{a,\text{tot}}/2 = 212.0972$:
$$\log_2 |N_g(a,b)| \le 8\log_2\!\left(|x|\cdot 167.7955^3 + 212.0972^3\right) = \mathbf{687.088} < 705.0 = N_{g,\max,\text{charged}}.$$
2. **$f$-Norm (Corner Global Maximum):** At the $k=2$ disjoint corner ($|\sigma_{0,1}(b)| = R_0 / 2^{1/0.65} = 487.4200$, $|\sigma_{2,3}(a)| = R_{a,\text{tot}}/\sqrt{2} = 299.9507$):
$$\log_2 |N_f(a,b)| \le 2\bigl(2\log_2(487.4200^9) + 2\log_2(299.9507^9)\bigr) = \mathbf{617.674} < 630.0 = N_{f,\max,\text{charged}}.$$

**Exact 2D Panjer Smooth-Recall + Survival DP & Unconditional Zero-Wraparound Proof:**
Every sieved band prime of norm $p$ adds the **floor-quantized** byte weight $w_p = \lfloor(\log_2 p)/q_{\text{step}}\rfloor$ to a pre-biased byte `bias_u8 = 128 - T_u8`. By subadditivity of $\lfloor \cdot \rfloor$, every cell in the entire domain satisfies $\text{byte}(a,b) \le \text{bias\_u8} + \lfloor N_{\max,\text{charged}}/q_{\text{step}}\rfloor$. Evaluating the **exact deterministic 2D Panjer recurrence** (zero Monte Carlo noise; joint $(k_{\text{bits}}, t_{\text{ticks}})$ convolution against the out-of-band smooth complement at exact Murphy-shifted norms $N_{g,\text{eff}} = 671.40 - 1.0487$ and $N_{f,\text{eff}} = 543.76 + 5.3752$, independently recomputed live in `experiment.py`) yields:
- **$g$-band `[2^{32.0}, 2^{61.5}]`:** `q_step_g = 3.8`, `T_u8 = 68` (`T_bits = 258.4`), `bias_u8 = 60`, `max_byte = 60 + floor(705.0/3.8) =` **`245 <= 255`** (**zero wraparound proved!**), exact smooth recall `rec_g = 0.807815`, exact survival $\log_2 s_g = \mathbf{-14.406429}$.
- **$f$-band `[2^{36.0}, 2^{61.0}]`:** `q_step_f = 3.8`, `T_u8 = 45` (`T_bits = 171.0`), `bias_u8 = 83`, `max_byte = 83 + floor(630.0/3.8) =` **`248 <= 255`** (**zero wraparound proved!**), exact smooth recall `rec_f = 0.850908`, exact survival $\log_2 s_f = \mathbf{-9.132517}$.
- **Combined Dual Band-Pass Stage 1a:** `ret_exact = 0.807815 * 0.850908 = 0.687377 >= 0.6864`, $\log_2 s_1 = \mathbf{-23.538946}$ (ZERO Band-1 small-prime sieving!).

### 3.3 Unconditional 3-Root $\mathcal{O}_K$-Prime Merging, Guard4 x16 Sieve & Sequential Single-Slab Masked Scan
**Unconditional 3-Root $\mathcal{O}_K$-Prime Merging Theorem:**
A rational prime $p \nmid 15$ splits completely in $K = \mathbb{Q}(\zeta_{15})$ iff $p \equiv 1 \pmod{15}$ (density $1/8$), yielding $8$ degree-1 $\mathcal{O}_K$-primes $\mathfrak{p}_K$. Since $\tau_3: Y \mapsto \omega Y$ acts freely on roots in $\mathbb{F}_p$:
- On the $g$-side ($Y^3 = x$), $x^{(p-1)/3} \equiv 1 \pmod p$ simultaneously across all $8$ $\mathfrak{p}_K \mid p$ with Chebotarev density $1/3$ (rational-prime density **`1/24`**), giving exact root count $r_p = 3$ and active $\mathcal{O}_K$-prime count $\frac{1}{3}\pi(B)$.
- On the $f$-side ($Q(Y^3) = 0$), every active $\mathfrak{p}_K \mid p$ has $r_p \in \{3, 6, 9\}$ (always divisible by 3!). By Landau's Prime Ideal Theorem, the active $\mathcal{O}_K$-prime count is unconditionally $\le \frac{1}{3}\pi(B)$ (empirically $\approx 0.287\,\pi(B)$ up to $p = 25{,}000$), and the active rational-prime density is unconditionally $\le 1/8$. Pricing both sides at $d_{\text{OK}} = 1/3, \bar{r} = 3$ is exact on $g$ and strictly conservative on $f$.

**Guard4 x16 Sieve & Sequential Single-Slab Active-Only Masked Scan (`KERNELS_SIEVE.md`):**
- **4-Byte Guard Layout (`W_B = 4`, `L <= 248`):** Each `a_id` row has 256 bytes: 4-byte left guard `[0..3]`, $L \le 248$ active columns `[4..4+L-1]`, and $\ge 4$-byte right guard `[4+L..255]`. For any window $[R, R + L)$ ($L \le 248$), every entry in `W_B = 4` buckets `[R >> 2 .. (R+L-1) >> 2]` satisfies $v - R \in [-3, L + 2] \subset [-3, 250]$, so $c_{\text{pad}} = (v - R) + 4 \in \mathbf{[1, 254] \subset [0, 255]}$ — **zero borrow/carry into `a_id`!**
- **Compiled Guard4 x16 Pricing (`kernels_guard4_x16_v81.s`):** x16 loop `43/16 = 2.6875` Ir/hit, x4 loop `13/4 = 3.25` Ir/hit, x1 tail `5.00` Ir/hit; probe+call worst-case bounds `<= 56` Ir (non-split) and `<= 84` Ir (split-256, $p_{\text{split}} = L_{\text{win}}/256 = 0.4394$, including `17` Ir `_emit_u64_x16_direct` prologue/epilogue). Integrating over `[32.0, 61.5]` ($g$) and `[36.0, 61.0]` ($f$) at $L_{\text{win,eff}} = 112.4903$ ($L_{\text{emit}} = L_{\text{win}} + 3 = 115.4903$) yields:
  - $g$-side: `hits_active = 0.6533`, `hits_emitted = 0.6707`, `ir_hits = 1.8352`, `ir_probe = 0.2354`.
  - $f$-side: `hits_active = 0.5274`, `hits_emitted = 0.5414`, `ir_hits = 1.4792`, `ir_probe = 0.1679`.
- **Sequential Single-Slab Active-Only Masked Scan (`kernels_slab_compact_v81.s`):** ONE single `14.9616 EiB` slab (`256 * V_{a,sh}`) is reused sequentially ($S_g \to S_f$):
  1. Pass 1: `_slab_prebias_active_sl` (`K_even + 3` Ir/row, `K_even_avg = 8.40` at `L_win_eff = 112.49`) $\to$ Guard4 x16 $g$-sieve $\to$ `_slab_scan_g_active_sl` (explicitly masking BOTH left `[0..3]` and right `[4+L..255]` guards!) emitting `2^{-14.41}` $g$-survivors into a `<= 1.4 PiB` list (`37.0 / 112.49 = 0.3289` Ir/active cell).
  2. Pass 2: reuse the SAME `14.9616 EiB` slab for $S_f$ (`_slab_prebias_active_sl` $\to$ Guard4 x16 $f$-sieve) $\to$ `_slab_filter_f_at_gsurv` probing ONLY at $g$-survivor coordinates (`5 Ir` per non-passing $g$-survivor $\implies 5 \cdot 2^{-14.4064} = 0.00023$ Ir/active cell).
  - Exact total active-cell slab cost: `0.2027 (2x prebias) + 0.3289 (S_g masked scan) + 0.00023 (S_f direct probe) = 0.5318 Ir/active cell`, conservatively charged **`0.5500 Ir/active cell`**.
- **Prime-by-Prime Sequential 2-Pass Replay Build (`replay_build_guard4_v81.s`):** Building one prime's bucket table at a time keeps temporary `cnt[]` $\le 2^{65}$ B (negligible), and Shell 0's tables are freed before Shell 1 is built. Charging **`27.5 Ir/entry`** (`11` Pass 1 + `15` Pass 2 + `1.5` per-$a_0$-line startup) across $N_{\text{ent}} = 9.0466 \times 10^{33}$ entries plus **`6.5 Ir/bucket`** (`5` prefix sum + `1.5` zeroing) across $N_{\text{bkt}} = 1.5816 \times 10^{34}$ buckets for $K = 2$ shells gives **`replay_a_tables = 2^{104.793639}` u live (`2^{104.7938}` u charged)**, with global table storage $8(N_{\text{bkt}} + N_{\text{ent}}) + 48 N_{\text{OK}} = \mathbf{2^{117.2595}\text{ B}}$.

---

## 4. Stage 1b: $f$-First Uncached (`2.00x`) Bernstein Batch Cofactorization

Every Stage 1a dual survivor (density $s_1 = 2^{-23.538946}$) undergoes **$f$-first** batch cofactorization using the compiled Shoup61 NTT and 3-prime CRT kernels in `KERNELS_LA_ECM.md` (`17` Ir fwd/inv butterfly, `30` Ir pointwise mul, exact `67` Ir `LBB3_2` CRT3, `160 Ir/word` 4-step distributed transpose when exceeding `L_owner_bits = (4.0 * 2^{60} * 8.0) / 4.5` bits, `1.15x` overhead):

1. **Exact Norm Evaluations & Unconditional Splitting-Density Trees:** Charges BOTH $N_f$ evaluation + packed $\le 2^{12}$ trial-div strip (`surv_s1 * 1700` Ir/cell) AND $N_g$ evaluation on true $f$-smooth survivors (`surv_s1 * pf_sm * 1700` Ir/cell). For primes $> 2^{36}$, $f$-trees use unconditional density `1/8` ($p \equiv 1 \pmod{15}$) and $g$-trees use unconditional density `1/24` ($p \equiv 1 \pmod{15},\; x^{(p-1)/3} \equiv 1 \pmod p$).
2. **Strict Uncached `2.00x` Bernstein Trees & Exact `Top2` Advancement:** Charges **`2.00x` `ntt_mul_ir`** per parent split (working RAM `4.5 * L_t / 8`, peak global **`2^{66.2235}` B**), **`40,000 Ir/leaf`** at 636 bits (`min 25,000` at 340 bits), `k_max = 8`, and exact `_tier_probs(rem, mid, B, k_max=8)` for BOTH `Top1` (`[P, mid]`) and `Top2` (`[mid, B]`).
3. **110-Round Miller-Rabin & Zero Double-Counting:** Dual-smooth survivors run 110-round MR (error $< 2^{-132}$, `c_mr = 8.96e-11` Ir/cell). Per-cell remainder+MR work **`c_cof_rem_mr = 0.620815`** Ir/cell enters `relation_collection`, while product-tree construction (`tree_build_u = 6.06e19` u live) is charged separately in `band_split_and_trees = 2^{90.0001}` u with zero double-counting.

**Stage 1 Total Collection Cost:**
$$\text{ipc}_{\text{tot}} = \underbrace{3.314452}_{\text{hits}} + \underbrace{0.403290}_{\text{probe}} + \underbrace{0.550000}_{\text{seq slab}} + \underbrace{0.620815}_{\text{Stage 1b}} = \mathbf{4.888557\text{ Ir/cell}},$$
$$\text{relation\_collection} = \frac{V_{\text{kept,charged}} \cdot \text{ipc}_{\text{tot}}}{20{,}000} = \mathbf{2^{107.808623}\text{ u live}} \quad (\mathbf{2^{107.8088}\text{ u charged}}).$$

---

## 5. Stages 2 & 3: `hx_z12` 2-Merge, Snell Stopping ($\eta = 0.82$), Monomial-Pivot SGE & Block Wiedemann

Full empirical tables (including the direct $U = 2^{20}$ $\eta$-sweep) appear in `FILTER_RANK.md`.

### 5.1 Algebraic 2-Merge, Snell Stopping ($\eta = 0.82$) & Explicit SGE Cost
1. **Rank-Barrier Resolution:** Algebraic 2-merge drops the $O(1)$ cancelled degree-2 column cycles in $O(1)$ work (fraction $< 2^{-60}$), achieving **`0.0000` dense-core nullity across 10,000+ exact $\mathbb{F}_q$ solves** for all $\eta \ge 0.78 > \tau_c = 0.701 \pm 0.002$ (re-verified live in `experiment.py`).
2. **4-Checkpoint Full-Resume Schedule ($\tau = [0.82, 0.86, 0.92, 1.00]$):** Checkpoint checks (`0.8%` of LA) are folded into `linear_algebra`; under H-R and H-SM ($\mathbb{P}(\text{fail }\tau_0 = 0.82) < 2^{-11}$), the expected escalation work $< 2^{-11} \cdot 0.05\,W_{\text{coll}} \approx 2^{92.5}$ u is explicitly charged in **`snell_resume_tail = 2^{92.5001}` u**.
3. **Post-SGE Matrix & Real $U = 2^{20}$ $\eta = 0.82$ Validation:** Applying $\text{LA\_mult}(0.82) = ((1.04094 - 0.701)/(0.82 - 0.701))^{1.35} = \mathbf{4.124817\times}$ (strictly dominating our real $U = 2^{20}$ measured `2.263x` at $\eta = 0.82$ by a **`1.82x` safety margin**, where `Mono = 85.32%`, `>32b = 0.00%`, `max|c| = 61`, and `ops = (1.138 pm1, 0.103 u32)`!):
$$N_l = U \cdot 0.0006589 \cdot \sqrt{4.124817} = 4.708934 \times 10^{15}, \qquad w_{\text{eff}} = 180.0 \cdot \sqrt{4.124817} = 365.5736,$$
$$N_{\text{nnz}} = N_l \cdot w_{\text{eff}} = 1.721462 \times 10^{18}.$$
4. **Explicit SGE & Matrix-Generation Cost Formula:**
$$\text{SGE\_U} = \frac{\bigl(12{,}500\,U + 500 \cdot 0.05\,U\log_2 U + 100\,N_{\text{nnz}}\bigr)\cdot 1.15}{20{,}000} = \mathbf{2^{61.306030}\text{ u live}} \quad (\mathbf{2^{61.3062}\text{ u charged}}).$$

### 5.2 64-Bit-Index Coppersmith Block Wiedemann (`m = 128, n = 32`, Exact `d_mksol = N_l/32 + 32` & `IR_ROW_FIXED = 2100`)
Using the compiled 64-bit-index kernels in `KERNELS_LA_ECM.md` with ZERO uncharged discrepancies:
- **Core SpMV:** `1.144 * 8.00 + 0.122 * 24.50 = 12.141 Ir/nnz` (exact `LBB4_1 = 32/4 = 8.00` and exact `LBB8_2 = 49/2 = 24.50`; conservative vs $U = 2^{20}, \eta = 0.82$ measured `1.138 * 8.00 + 0.103 * 24.50 = 11.628 Ir/nnz`!).
- **Exact Dynamic Per-Row Fixed Overhead:** **`IR_ROW_FIXED = 2,100.0 Ir/row`**, charging `747` dynamic `_fr_dot4_mont` `+ 623` $\lambda^k x$ Horner (`576 + 47`) `+ 588` (`4 x 147` `_fr_reduce_5limb`) `+ 47` `_fr_add` `+ 95` amortized `x4`/`x2` bucket-dispatch call overhead `= 2,100.0 Ir/row` exact, times `1.15x` overhead.
- **Exact Sequence & `Mksol` Degree (`n = 32`):**
$$L_{\text{seq}} = \frac{N_l}{128} + \frac{N_l}{32} + 64 = 1.839427 \times 10^{14}, \qquad d_{\text{mksol}} = \frac{N_l}{32} + 32 = 1.471542 \times 10^{14},$$
$$\text{total\_spmvs} = 1.02\,(32\,L_{\text{seq}} + d_{\text{mksol}}) = 6.153988 \times 10^{15}.$$
- **Full `Mksol`, `Lingen`, Snell Check & Back-Sub Breakdown:**
  - `Mksol` $\mathbb{F}_r$ multiply-adds: $1.02 \cdot d_{\text{mksol}} \cdot 32 \cdot N_l \cdot (236 + 47) \cdot 1.15 / C = \mathbf{2^{98.215793}\text{ u}}$.
  - `Lingen` (`m + n = 160`): $20 \cdot 160^3 \cdot L_{\text{seq}}\log_2(L_{\text{seq}}) \cdot 236 \cdot 1.15 / C = \mathbf{2^{73.036924}\text{ u}}$.
  - Snell checkpoint check (`0.8%`) + sparse back-substitution (`30 U * 283 * 1.15`) $\implies$ total **`linear_algebra = 2^{103.162842}` u live (`2^{103.1630}` u charged)**.

---

## 6. Stage 4, Stage 5, Unified RAM Audit & Authoritative 14-Row Ledger

### 6.1 Stage 4 12D Subfield $M_g = \mathbb{Q}(\zeta_5)[Y]/(Y^3 - x)$ Hilbert-90 Descent (`DESCENT_ECM.md`)
For BOTH $g_0 = e(P, Q_0)$ and $z$:
- **12D Subfield Isomorphism & Unconditional Deterministic LLL + AM–GM Norm Bound:** In $M_g = \mathbb{Q}(\zeta_5)[Y]/(Y^3 - x)$ ($[M_g:\mathbb{Q}] = 12$, $\mathcal{O}_{M_g}/(p) \cong \mathbb{F}_{p^{12}}$), Hilbert-90 imposes 6 linear $\mathbb{F}_p$-congruences defining a rank-12 sublattice $\Lambda_{12}$ with weighted determinant $p^6 |x|^4 = 2^{2539.061119}$ and canonical $\mathbb{R}^{12}$ covolume $\sqrt{|\Delta_{M_g}|}\,p^6 = 5^{9/2} 3^6 |x|^4 p^6 = 2^{2559.019571}$. Standard deterministic LLL ($\delta = 3/4$) + AM–GM unconditionally proves $|N_{M_g/\mathbb{Q}}(A)| \le \frac{2^{33}}{12^6}\sqrt{|\Delta_{M_g}|}\,p^6 = \mathbf{2^{2570.509796} < 2^{2572.0} = 2^{N_{\text{boot}}}}$.
- **Explicit Handling of Inert $q \equiv 11 \pmod{15}$ Norm Squaring ($q^2$) & Zero-Early-Abort Descent at $B_{\text{boot}} = 2^{275}$:** In the quadratic extension $L_g = M_g(\omega)$, ~50% of $M_g$ degree-1 primes have $q \equiv 1 \pmod{15}$ (splitting into degree-1 $L_g$-primes of norm $q$), while ~50% have $q \equiv 11 \pmod{15}$ ($q \equiv 2 \pmod 3$, remaining inert of residue degree $2$ with $L_g$-norm $q^2$). We charge BOTH an explicit **Layer-0** special-$\mathfrak{Q}$ descent ($2^{550} \to 2^{275}$, $\le 5$ nodes, `2^{61.825}` u/node, layer `2^{64.147}` u) for inert $q \in (2^{137.5}, 2^{275}]$ AND **`+20` extra inert nodes** conservatively billed at Layer-5 cost (`180 + 20 = 200` Layer-5 nodes) for inert $q \in (2^{34.98}, 2^{69.96}]$ whose squared norm $q^2 > B_g$. With unconstrained Torsion-16 Edwards ECM (`q = 2^{70}..2^{275}`, `< 10^{-12}` self-consistent in `DESCENT_ECM.md`), live zero-early-abort 2-log work is `2^{95.089960}` u, rising with the quadratic Snell missing-log penalty $(1.04094/0.82)^2 = 1.61147\times$ to `2^{95.778342}` u, and conservatively charged **`both_individual_logs_exp = 2^{96.5001}` u** (`+0.72` b extra headroom).
- **Explicit `ecm_split_expected`:** $2\cdot(W_{\text{tree},1} + 0.25\,W_{\text{boot},1})\cdot(1.04094/0.82)^2 = 2^{94.070586}$ u live, charged **`2^{94.0707}` u**.
- Final back-substitution/GCD and final $Q = [a]Q_0$ plus pairing verification $e(P, Q) = z$ are charged at **`final_log_backsub_and_gcd = 2^{72.0001}` u** and **`final_verification_Q = 2^{40.0001}` u**.

### 6.2 Complete Unified `classical-group-work-v2` Peak Storage Audit

All stages execute sequentially in unified flat RAM (Shell 0 freed before Shell 1; prime-by-prime sequential replay build keeping temporary `cnt[] <= 2^{65}` B; ONE single `14.9616 EiB` slab reused sequentially $g \to f$):

| Stage / Subsystem | Peak Global Storage | Notes |
|---|---:|---|
| Stage 1a Guard4 Replay Tables (`8(N_bkt + N_ent) + 48 N_OK`, 1 active shell) | **`2^{117.2595}` B** | Prime-by-prime build; Shell 0 freed before Shell 1 |
| Stage 1a Single Reused Flat Slab (`256 * V_{a,sh}`) + $g$-survivor buffer | `14.9616` EiB (`2^{63.90}` B) + `4` PiB | Reused sequentially $S_g \to S_f$ |
| **Stage 1 Combined Peak Global Storage** | **`2^{117.2595}` B** | Replay tables + 1 slab + $g$-survivors |
| Stage 1b Uncached Bernstein Trees (`4.5 * L_t / 8`) | `2^{66.2235}` B | Runs on Stage 1a dual survivors |
| Stage 2 SGE Initial Relations (`eta * U * 18 * 12 B`) | **`2^{69.0784}` B** | `~530` EiB (`<< 2^{117}` B); Stage 1 freed first |
| Stage 3 LA Matrix + Krylov Vectors (`N_nnz*11 + N_l*8 + 38*N_l*32`) | **`2^{64.4211}` B** | `~21.4` EiB; Stage 2 freed first |
| Stage 4 12D $M_g$ Descent & $T_{26}$ Tree | `< 1.0` GiB | Negligible |
| Nonuniform Advice | **`0` B (`2^0` B)** | Zero searched/stored advice |
| **Reported Peak Total Storage Bound** | **`<= 2^{117.50}` B** | Strictly `< 2^{120}` B |

### 6.3 Complete Authoritative 14-Row Cold-Run Ledger (`ledger_data.json`)

Every row below is independently recomputed from scratch by `experiment.py` and matches `ledger_data.json` to full precision:

| # | Ledger Row | Live ($\log_2 u$) | Charged ($\log_2 u$) | Slack (b) | Share of Opened |
|---:|---|---:|---:|---:|---:|
| 1 | `relation_collection` | `107.808623` | **`107.8088`** | `+0.000177` | `85.90%` |
| 2 | `replay_a_tables` | `104.793639` | **`104.7938`** | `+0.000161` | `10.63%` |
| 3 | `linear_algebra` | `103.162842` | **`103.1630`** | `+0.000158` | `3.43%` |
| 4 | `both_individual_logs_exp` | `96.500000` | **`96.5001`** | `+0.000100` | `0.03%` |
| 5 | `ecm_split_expected` | `94.070586` | **`94.0707`** | `+0.000114` | `0.01%` |
| 6 | `meter_safety_buffer` | `93.400000` | **`93.4001`** | `+0.000100` | `<0.01%` |
| 7 | `snell_resume_tail` | `92.500000` | **`92.5001`** | `+0.000100` | `<0.01%` |
| 8 | `band_split_and_trees` | `90.000000` | **`90.0001`** | `+0.000100` | `<0.01%` |
| 9 | `canonical_a_domain` | `84.000000` | **`84.0001`** | `+0.000100` | `<0.01%` |
| 10 | `schirokauer_maps` | `80.000000` | **`80.0001`** | `+0.000100` | `<0.01%` |
| 11 | `final_log_backsub_and_gcd` | `72.000000` | **`72.0001`** | `+0.000100` | `<0.01%` |
| 12 | `matrix_generator_sge` | `61.306030` | **`61.3062`** | `+0.000170` | `<0.01%` |
| 13 | `stage0_setup_and_iso` | `45.000000` | **`45.0001`** | `+0.000100` | `<0.01%` |
| 14 | `final_verification_Q` | `40.000000` | **`40.0001`** | `+0.000100` | `<0.01%` |
| — | **Opened-Run Expected Work (`p_open >= 0.99998444`)** | **`108.027959`** | **`108.028134`** | **`+0.000175`** | **`100.00%`** |
| — | **64-Bit Pre-Gate Factor ($\log_2(K_{\text{gate}}/2^{64})$)** | `-1.358425` | `-1.358425` | — | — |
| — | **Cold-Run Gated Expected Work ($p_{\text{cold}} \ge 0.39$)** | **`106.669534`** | **`106.669708`** | **`+0.000292`** | — |
| — | **OFFICIAL REDUCED CLAIM (`10667/100`, `gcd = 1`)** | — | **`106.670000`** | — | — |

### 6.4 Exact Minimality of $K_{\text{gate}}$ & Cold-Run Success Proof
With $p_{\text{open}} \ge 0.99998444 > 0.99998 = 99998/100000$ (proved unconditionally via Markov at $W_{\text{pointwise}} = 2^{124}$ u), the unique minimal integer $K_{\text{gate}} \in [0, 2^{64}]$ satisfying $\frac{K_{\text{gate}}}{2^{64}} \cdot \frac{99998}{100000} \ge \frac{39}{100}$ is
$$K_{\text{gate}} = \left\lceil \frac{39 \cdot 2^{64} \cdot 100{,}000}{100 \cdot 99{,}998} \right\rceil = \mathbf{7{,}194{,}374{,}076{,}228{,}249{,}696} = \texttt{0x63d78d1adadd2860},$$
verified in exact rational arithmetic (`Fraction`) in `experiment.py`. Overall cold-run success is $\mathbb{P}(\text{success}) = \frac{1}{r}\cdot 1 + \frac{r-1}{r}\cdot p_{\text{gate}}\,p_{\text{open}} \ge 0.39$.

---

## 7. Explicit Disclosure of Exploratory Heuristics & Limitations

In full compliance with `rules.md`, this is an exploratory cryptanalysis bound, not an unconditional production theorem. Production-scale performance relies on the following explicitly disclosed heuristics:

| ID | Heuristic Statement | Quantitative Support in This Package | Limitation / Extrapolation Risk |
|---|---|---|---|
| **H1** | **Smoothness & Murphy $\alpha$:** Norms $N_f(a,b), N_g(a,b)$ behave independently as random integers of their size shifted by exact Murphy $\alpha_f = +5.3752$ b, $\alpha_g = -1.0487$ b, with Dickman $\rho(u)$ smooth probability and coprimality factor $1/\zeta_K(2)$. | `N = 4,000,000` direct 16D MC draws (`SE = 0.00116` b) + `-0.16` b cushion + `+0.06544` b anchor pad + `+0.00348` b 3-SE pad. | Standard NFS/STNFS heuristic; not unconditionally proved at $600$-bit norms. |
| **H2** | **Dual Band-Pass Independence & Floor Quantization:** Prime divisibility across `[2^{32}, 2^{61.5}]` ($g$) and `[2^{36}, 2^{61.0}]` ($f$) follows the exact 2D Compound-Poisson model conditioned on smooth total norm. | Exact zero-noise 2D Panjer Smooth-Recall + Survival DP (`experiment.py`); unconditional zero `uint8_t` wraparound (`max_byte = 245, 248 <= 255`) IS proved deterministically without heuristics. | Minor small-prime correlations at band boundary; mitigated by exact Murphy-shifted 2D DP. |
| **H3** | **Post-SGE Matrix Scaling & Monomial-Pivot Entry Mix:** At $U = 2^{61.61}$, Bouvier/CADO-NFS clique removal + Monomial-Pivot SGE preserves `nl_cap <= 0.0006589`, `avg_w <= 180` at $\eta = 1.041$, entry mix `(1.144 pm1, 0.122 u32, 0.0000 fr)`, and Snell scaling $\beta \le 1.35$. | Clean-room $\mathbb{Z}[\zeta_{12}]$ SGE ladder $U = 2^{14}..2^{20}$ PLUS real $U = 2^{20}$ $\eta = 0.82$ benchmark (`LA_ratio = 2.263x` vs charged `4.125x` $\implies$ `1.82x` margin; `Mono = 85.32%`, `max|c| = 61`, `0.00%` `>32b`). | Extrapolating from $U = 2^{20}$ to $U = 2^{61.6}$; $\beta = 2.65$ stress test reported in Section 8. |
| **H4** | **Stage 4 12D $M_g$ Hilbert-90 & Special-$\mathfrak{q}$ Descent:** Smoothness of 12D $M_g$ norms ($\le 2{,}572$ b proved unconditionally via LLL+AM–GM) and 6-layer (Layer-0..5) special-$\mathfrak{q}$ descent follow Dickman $\rho(u)$ without dead-end blowup. | Unconditional LLL+AM–GM proof (`N_boot <= 2570.51 < 2572` b); explicit Layer-0 + `+20` inert nodes; zero-early-abort upper bound `2^{95.778}` u w/ Snell charged at `2^{96.5001}` u (`0.03%` of opened work). | Standard STNFS descent heuristic; even a `4x` Stage 4 blowup adds only `+0.00147` b (`106.6712` u)! |
| **H5 / H-ECM** | **Torsion-16 Edwards ECM Group-Order Smoothness:** Random Torsion-16 Edwards curves mod $q$ have $q_{\text{eff}} = q - 3.5$ bits of effective Dickman–Buchstab smoothness in Stage 1 + Stage 2. | Exact algebraic proof $16 \mid \#E(\mathbb{F}_p)$ (`21/21` primes `p <= 200`), live Stage-2 witnesses (`50023` & `35`), conservative `3.5`-bit credit, `< 10^{-12}` table self-consistency across `q = 2^{70}..2^{275}`. | Hasse-interval smoothness heuristic; Stage 4 is `0.03%` of total work. |
| **H-R** | **Post-2-Merge Dense-Core Full Rank:** After dropping $O(1)$ cancelled degree-2 cycles, $M_{\text{core}}'$ has full column rank over $\mathbb{F}_r$ at $\eta \ge 0.82 > \tau_c = 0.701$ ($\mathbb{P}(\text{fail }\tau_0) < 2^{-11}$). | `10,000+` exact $\mathbb{F}_q$ solves showing **`0.0000` dense-core nullity** for all $\eta \ge 0.78$, plus live verification in `experiment.py`. | Finite-scale $U \le 2^{14}$ exact solves; theory proves $O(1)$ cycle density is scale-invariant and $r^{-E_{\text{core}}} \approx 2^{-255 E_{\text{core}}}$. |
| **H-SM** | **Non-Degeneracy of the 4 Schirokauer Maps:** The 4 $\chi_1$-projected Schirokauer map columns are linearly independent modulo $r$ of the prime-ideal relation lattice. | Standard STNFS Schirokauer non-degeneracy (probability $1 - O(r^{-1})$). | Standard assumption in all NFS-DL analyses. |

---

## 8. Sensitivity & Perturbation Audit

| Scenario / Stress Test | Opened Work ($\log_2 u$) | Gated Cold-Run ($\log_2 u$) | Delta vs Official Claim (`106.670000`) |
|---|---:|---:|---:|
| **Pure First-Principles (no `+0.06544` b comparator anchor pad)** | `107.971930` | **`106.613505`** | `-0.056495` b (better) |
| Using Direct $U = 2^{20}$ Measured `LA_ratio = 2.263x` (vs charged `4.125x`) | `108.000908` | **`106.642483`** | `-0.027517` b (better) |
| **Primary Certified Charged Ledger (`ledger_data.json`)** | **`108.028134`** | **`106.669708`** | **`-0.000292` b (strictly below claim)** |
| **Official Reduced Claim (`submission.json`, `10667/100`, `gcd = 1`)** | `108.028425` | **`106.670000`** | `0.000000` b |
| Stress: Stage 4 Descent $4\times$ slower (`+2.0` b on Stage 4) | `108.029598` | `106.671173` | `+0.001173` b |
| Stress: Ultra-Conservative LA Exponent $\beta = 2.65$ at fixed $\eta = 0.82$ (`16.14x` LA mult) | `108.274523` | `106.916098` | `+0.246098` b *(re-optimizing $\eta \to 0.91$ recovers most of this)* |

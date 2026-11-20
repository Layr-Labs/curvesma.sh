# `RANK-BRIDGE.md` — Exact 1D-Kernel Rank Bridge (`H-R`, `H-X`, & Conditional Schwartz-Zippel)

This note separates the three logical layers connecting the filtered $\mathbb{Z}[\omega]$ relation matrix to the virtual-logarithm kernel recovered by Block Wiedemann, and documents the live verification executed inside `hx_z6.py` (invoked by `audit_bn_cm.py`).

## 1. The Three Separated Premises

Let $M_{\mathrm{full}} \in \mathbb{Z}[\omega]^{R \times U'}$ be the raw relation matrix over order-6 $\tau$-orbits ($\tau = \operatorname{Frob}^{10}$, eigenvalue $\chi = p^{10} \bmod r$ with $\chi^2 - \chi + 1 \equiv 0 \pmod r$), augmented by $s_{\mathrm{sm}} = 5$ dense Schirokauer map columns $S_{\mathrm{full}} = [W_f(P_\chi \Lambda_f) \mid W_g(P_\chi \Lambda_g)] \in \mathbb{F}_r^{R \times 5}$, where $\Lambda_f \in \mathbb{F}_r^{R \times 48}$ and $\Lambda_g \in \mathbb{F}_r^{R \times 12}$ are the un-discounted Schirokauer maps ($\alpha^{r^2-1} \bmod r^2$ across all 6 CRT components of $h(y) \bmod r$), $P_\chi$ projects onto the $\chi$-eigenspace of $\mathcal{O}_{K_f}^\times/(\mathcal{O}_{K_f}^\times)^r$ (dimension $4$) and $\mathcal{O}_{K_g}^\times/(\mathcal{O}_{K_g}^\times)^r$ (dimension $1$), and $W_f \in \mathbb{F}_r^{4 \times 48}$, $W_g \in \mathbb{F}_r^{1 \times 12}$ are uniform random compression matrices drawn at run start (under `H-SM`, random compression loses rank with probability $\le 5/r < 2^{-250}$, included directly in `FIX-U`).

After singleton/clique purge, Markowitz merge with deterministic pre-merge guard $W_{\mathrm{avg}} \le 161.0$, and Boudot order-labeled recorded pivots $L_{\mathrm{rec}}$ (strictly lower-triangular unit diagonal), the surviving live system has $N_l$ live columns, $N_{\mathrm{rec}} = N_{\mathrm{sw}} - N_l$ recorded pivot columns, and exact Schur-complement product operator
$$M_{\mathrm{prod}} = M_{\mathrm{live}} - E_{\mathrm{live}} L_{\mathrm{rec}}^{-1} B_{\mathrm{rec}}, \qquad S_{\mathrm{prod}} = S_{\mathrm{live}} - E_{\mathrm{live}} L_{\mathrm{rec}}^{-1} S_{\mathrm{rec}}, \qquad A = [M_{\mathrm{prod}} \mid S_{\mathrm{prod}}] \in \mathbb{F}_r^{D \times D}$$
on $D = N_l + s_{\mathrm{sm}}$ unscreened live rows (with zero-padded operator $B = [M_{\mathrm{prod}} \mid 0]$).

1. **Layer (a) — `H-R` (Finite Relation Saturation / 1D Kernel):**
   The un-augmented Schur-complement operator $M_{\mathrm{prod}} \in \mathbb{F}_r^{D \times N_l}$ has full column rank $\operatorname{rank}(M_{\mathrm{prod}}) = N_l$, and the augmented matrix $A = [M_{\mathrm{prod}} \mid S_{\mathrm{prod}}] \in \mathbb{F}_r^{D \times D}$ has exact rank $D - 1$. Writing $V = \operatorname{Im}(M_{\mathrm{prod}}) + \operatorname{span}(S_{\mathrm{prod}}) = \operatorname{Im}(A)$ ($\dim V = D - 1$), the right kernel is strictly one-dimensional:
   $$\ker(A) = \operatorname{span}_{\mathbb{F}_r}\{(w_{\mathrm{true}}, c_{\mathrm{true}})\}, \qquad c_{\mathrm{true}} \ne 0,$$
   where $w_{\mathrm{true}}$ is the vector of virtual logarithms of the $N_l$ live factor-base orbits (and recorded pivots are recovered uniquely by forward substitution $L_{\mathrm{rec}} w_{\mathrm{rec}} = -(B_{\mathrm{rec}} w_{\mathrm{true}} + S_{\mathrm{rec}} c_{\mathrm{true}})$).

2. **Layer (b) — `H-X` (Padded Block Wiedemann Controllability, Observability, & Recovery Theorem):**
   Let $B = [M_{\mathrm{prod}} \mid 0] \in \mathbb{F}_r^{D \times D}$ ($\operatorname{Im}(B) = \operatorname{Im}(M_{\mathrm{prod}})$ of dimension $N_l$). Take $Z \in \mathbb{F}_r^{N_l \times (n - 5)}$ with **dense independent uniform 64-bit entries** embedded in $\mathbb{F}_r$, set $Y = [S_1, \dots, S_5, B Z] \in \mathbb{F}_r^{D \times n}$ (fixing the $5$ Schirokauer columns as the first $5$ start columns), and let $X \in \mathbb{F}_r^{m \times D}$ have the **fixed sparse block-diagonal support** placing the $5$ pad coordinates in distinct blocks $0, \dots, 4$ and partitioning the $N_l$ non-pad coordinates across the $m$ blocks, with independent uniform 64-bit entries on its support. At the charged depths
   $$d = \lceil D/n \rceil + 1 + \lceil D/(128 n) \rceil, \qquad e = \lceil D/m \rceil + 1 + \lceil D/(128 m) \rceil,$$
   we declare as fixed-instance exploratory premises (supported diagnostically, not proved, by the finite-scale `hx_z6.py` witness):
   - **Controllability:** $\operatorname{rank}(C_d) = N_l$ for $C_d = [B Z, B Y, B^2 Y, \dots, B^d Y]$ (so $\operatorname{span}(C_d) = \operatorname{Im}(B)$);
   - **Observability:** $O_e = [X; X B; \dots; X B^{e-1}]$ is injective on $V = \operatorname{Im}(A)$, i.e., $\operatorname{rank}(O_e A) = \operatorname{rank}(A) = D - 1$.

   **Recovery Theorem (4 lines, exact):** Let $R_d = [Y, B Y, \dots, B^d Y] \in \mathbb{F}_r^{D \times (d+1)n}$ so the projected block-Hankel matrix is $H = [X B^{i+j} Y]_{0 \le i < e,\; 0 \le j \le d} = O_e R_d$.
   - *(Existence)* Given `H-R` and `H-X`, any nonzero $(w, c) \in \ker(A)$ satisfies $c \ne 0$ (since $\operatorname{rank}(M_{\mathrm{prod}}) = N_l$), hence $S_{\mathrm{prod}} c = -M_{\mathrm{prod}} w \in \operatorname{Im}(B) = \operatorname{span}(C_d)$. Writing $S_{\mathrm{prod}} c = -B Z e_0 - \sum_{i=1}^d B^i Y f_i$ and setting $f_0 = (c, e_0)$ gives $f = (f_0, f_1, \dots, f_d)$ with $c \ne 0$ and $R_d f = 0$, hence $H f = O_e R_d f = 0$.
   - *(Soundness)* Conversely, for any $f = (f_0, \dots, f_d)$ with $H f = 0$ and nonzero first $5$ coordinates $c = f_0[:5] \ne 0$, since $R_d f \in V$ and $O_e$ is injective on $V$, we obtain $R_d f = 0$, i.e. $S_{\mathrm{prod}} c + B Z e_0 + \sum_{i=1}^d B^i Y f_i = 0$. Setting $\bar{w} = Z e_0 + \sum_{i=1}^d (B^{i-1} Y f_i)[:N_l]$ yields $M_{\mathrm{prod}} \bar{w} + S_{\mathrm{prod}} c = 0$, so $(\bar{w}, c) \in \ker(A) = \operatorname{span}\{(w_{\mathrm{true}}, c_{\mathrm{true}})\}$ with $c \ne 0$, recovering $(w_{\mathrm{true}}, c_{\mathrm{true}})$ up to a nonzero scalar.

3. **Layer (c) — Exact Conditional Schwartz-Zippel Failure Bound:**
   Under `H-R` and `H-X` (existence of at least one nonzero $N_l$-minor of $C_d$ in the dense entries of $Z$ and one nonzero $(D-1)$-minor of $O_e U$ for a basis $U$ of $V$ in the entries of $X$), the product of those two minors is a nonzero polynomial of total degree at most $N_l + (D - 1) < 4 N_l \le 4 N_{\mathrm{sw,cap}}$ in the independent uniform 64-bit random draws of $(X, Z)$. By Schwartz-Zippel over $\{0, \dots, 2^{64}-1\} \subset \mathbb{F}_r$, the probability that either rank condition fails is at most
   $$\frac{N_l + (D - 1)}{2^{64}} \le \frac{4 N_{\mathrm{sw,cap}}}{2^{64}} = 4.875 \times 10^{-6}$$
   (charged directly in `FIX-U` as `main_matrix`).

## 2. Live Verification in `hx_z6.py`

`audit_bn_cm.py` invokes `hx_z6.run_hx_1d_kernel_z6(log2U=12, seed=1, m=12, n=8, s_sm=5)` live on an unscreened filtered $\mathbb{Z}[\omega]$ system modulo $Q = 1{,}000{,}033$ ($\chi^2 - \chi + 1 \equiv 0 \pmod Q$, $\chi^3 \equiv -1 \pmod Q$, matching BN254 `C6`'s exact `s_sm = 5` Schirokauer columns):

- Verifies strict lower-triangularity of $L_{\mathrm{rec}}$ and exact operator equality `op(v) == M_prod * v` (`prod_eq_schur = True`);
- Plants an exact 1D kernel $(w_{\mathrm{true}}, c_{\mathrm{true}})$ ($c_{\mathrm{true}}[0] = 1$) into $S_{\mathrm{prod}}$ (`N = 16, D = 21`) and verifies $\operatorname{rank}(A) = D - 1 = 20$ (`hr_1d_kernel_exact = True`);
- Verifies controllability $\operatorname{rank}(C_d) = 16 = N$, observability $\operatorname{rank}(O_e A) = 20 = \operatorname{rank}(A)$, and the negative one-block-pad control $\operatorname{rank}(O_e^{\mathrm{one}} A) < \operatorname{rank}(A)$ (`17 < 20`);
- Forms the exact block-Hankel matrix $H$ (`48 x 48`), computes $\ker(H)$, maps $f \mapsto (w(f), c(f))$, and verifies exact recovery $(w(f), c(f)) \equiv \alpha (w_{\mathrm{true}}, c_{\mathrm{true}}) \pmod Q$ with $\alpha \ne 0$ (`hankel_exact_1d_kernel_recovered = True`).

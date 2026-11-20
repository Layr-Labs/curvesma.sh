# `RANK-BRIDGE.md` — Exact 1D-Kernel Rank Bridge (`H-R`, `H-X`, & Conditional Schwartz-Zippel)

This note separates the three logical layers connecting the filtered $\mathbb{Z}[\omega]$ relation matrix to the virtual-logarithm kernel recovered by Block Wiedemann, and documents the live verification executed inside `hx_z6.py` (invoked by `audit_bn_cm.py`).

## 1. The Three Separated Premises

Let $M_{\mathrm{full}} \in \mathbb{Z}[\omega]^{R \times U'}$ be the raw relation matrix over order-6 $\tau$-orbits ($\tau = \operatorname{Frob}^{10}$, eigenvalue $\chi = p^{10} \bmod r$ with $\chi^2 - \chi + 1 \equiv 0 \pmod r$), augmented by $s_{\mathrm{sm}} = 5$ dense Schirokauer map columns $S_{\mathrm{full}} \in \mathbb{F}_r^{R \times 5}$ (computed un-discounted via $\alpha^{r^2-1} \bmod r^2$ across all 6 CRT components of $h(y) \bmod r$).

After singleton/clique purge, Markowitz merge with deterministic pre-merge guard $W_{\mathrm{avg}} \le 161.0$, and Boudot order-labeled recorded pivots $L_{\mathrm{rec}}$ (strictly lower-triangular unit diagonal), the surviving live system has $N_l$ live columns, $N_{\mathrm{rec}} = N_{\mathrm{sw}} - N_l$ recorded pivot columns, and exact Schur-complement product operator
$$M_{\mathrm{prod}} = M_{\mathrm{live}} - E_{\mathrm{live}} L_{\mathrm{rec}}^{-1} B_{\mathrm{rec}}, \qquad S_{\mathrm{prod}} = S_{\mathrm{live}} - E_{\mathrm{live}} L_{\mathrm{rec}}^{-1} S_{\mathrm{rec}}, \qquad A = [M_{\mathrm{prod}} \mid S_{\mathrm{prod}}] \in \mathbb{F}_r^{D \times D}$$
on $D = N_l + s_{\mathrm{sm}}$ unscreened live rows (with zero-padded operator $B = [M_{\mathrm{prod}} \mid 0]$).

1. **Layer (a) — `H-R` (Finite Relation Saturation / 1D Kernel):**
   The augmented Schur-complement matrix $A = [M_{\mathrm{prod}} \mid S_{\mathrm{prod}}]$ over $\mathbb{F}_r$ has exact rank $D - 1$, so its right kernel is strictly one-dimensional:
   $$\ker(A) = \operatorname{span}_{\mathbb{F}_r}\{(w_{\mathrm{true}}, c_{\mathrm{true}})\}, \qquad c_{\mathrm{true}} \ne 0,$$
   where $w_{\mathrm{true}}$ is the vector of virtual logarithms of the $N_l$ live factor-base orbits (and recorded pivots are recovered uniquely by forward substitution $L_{\mathrm{rec}} w_{\mathrm{rec}} = -(B_{\mathrm{rec}} w_{\mathrm{true}} + S_{\mathrm{rec}} c_{\mathrm{true}})$).

2. **Layer (b) — `H-X` (Finite Controllability & Observability at Charged Depths $d, e$):**
   Let $Z \in \mathbb{F}_r^{N_l \times (n - s_{\mathrm{sm}})}$ be dense random in $\mathbb{F}_r$ (matching production), let $Y = [S_{\mathrm{prod}}, B Z] \in \mathbb{F}_r^{D \times n}$, and let $X \in \mathbb{F}_r^{m \times D}$ have the exact **sparse block-diagonal support** placing each non-pad coordinate in one uniformly shuffled block of $m$ and placing the $s_{\mathrm{sm}}$ Schirokauer pad coordinates in **distinct** blocks $0, \dots, s_{\mathrm{sm}}-1$. At the charged Block-Wiedemann depths
   $$d = \lceil D/n \rceil + 1 + \lceil D/(128 n) \rceil, \qquad e = \lceil D/m \rceil + 1 + \lceil D/(128 m) \rceil,$$
   we declare as fixed-instance exploratory premises (supported diagnostically, not proved, by the finite-scale `hx_z6.py` witness):
   - **Controllability:** $\operatorname{rank}(C_d) = N_l$, where $C_d = [B Y, B^2 Y, \dots, B^d Y]$;
   - **Observability:** $\operatorname{rank}(O_e A) = \operatorname{rank}(A) = D - 1$, where $O_e = [X; X B; \dots; X B^{e-1}]$.

   Under `(a) + (b)`, the projected Block-Wiedemann Hankel matrix $H = [X B^{i+j} Y]_{0 \le i < e,\; 0 \le j \le d}$ satisfies the exact algebraic identity
   $$H f = O_e A \begin{pmatrix} w(f) \\ c(f) \end{pmatrix}$$
   for every $f = (f_0, \dots, f_d) \in \mathbb{F}_r^{(d+1)n}$, where $c(f) = f_0[:s_{\mathrm{sm}}]$ and $w(f) = \sum_{t=s_{\mathrm{sm}}}^{n-1} f_0[t] u_{t-s_{\mathrm{sm}}+1} + \sum_{j=1}^d (B^{j-1} Y f_j)[:N_l]$. Because $\operatorname{rank}(O_e A) = \operatorname{rank}(A)$, $\ker(H)$ maps **exclusively** into $\ker(A) = \operatorname{span}\{(w_{\mathrm{true}}, c_{\mathrm{true}})\}$; because $\operatorname{rank}(C_d) = N_l$, there exists $f \in \ker(H)$ with $c(f) \ne 0$, recovering $(w_{\mathrm{true}}, c_{\mathrm{true}})$ up to a nonzero scalar.

3. **Layer (c) — Conditional Schwartz-Zippel Failure Bound:**
   Conditional on `(a) H-R` and `(b) H-X`, choosing the random vectors in $X$ and $Y$ from fresh 64-bit coins (embedded in $\mathbb{F}_r$) causes the projected Krylov/Hankel step to lose rank with probability at most $4 N_{\mathrm{sw,cap}} / 2^{64} = 4.875 \times 10^{-6}$ (charged directly in `FIX-U` as `main_matrix`).

## 2. Live Verification in `hx_z6.py`

`audit_bn_cm.py` invokes `hx_z6.run_hx_1d_kernel_z6(log2U=12, seed=1, m=12, n=8)` live on an unscreened filtered $\mathbb{Z}[\omega]$ system modulo $Q = 1{,}000{,}033$ ($\chi^2 - \chi + 1 \equiv 0 \pmod Q$, $\chi^3 \equiv -1 \pmod Q$):

- Verifies strict lower-triangularity of $L_{\mathrm{rec}}$ and exact operator equality `op(v) == M_prod * v` (`prod_eq_schur = True`);
- Plants an exact 1D kernel $(w_{\mathrm{true}}, c_{\mathrm{true}})$ ($c_{\mathrm{true}}[0] = 1$) into $S_{\mathrm{prod}}$ and verifies $\operatorname{rank}(A) = D - 1$ (`hr_1d_kernel_exact = True`);
- Verifies controllability $\operatorname{rank}(C_d) = N$, observability $\operatorname{rank}(O_e A) = \operatorname{rank}(A)$, and the negative one-block-pad control $\operatorname{rank}(O_e^{\mathrm{one}} A) < \operatorname{rank}(A)$ (`17 < 19`);
- Forms the exact block-Hankel matrix $H$ (`48 x 48`), computes $\ker(H)$, maps $f \mapsto (w(f), c(f))$, and verifies exact recovery $(w(f), c(f)) \equiv \alpha (w_{\mathrm{true}}, c_{\mathrm{true}}) \pmod Q$ with $\alpha \ne 0$ (`hankel_exact_1d_kernel_recovered = True`).

# Probabilistic Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`, `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR`)

All `39` bins (`78` nodes, `70..265`b), exact `cells_rank2`, `1/zeta_{Q(zeta_15)}(2)`, and seeded peeling (`h = 0.256845`) are evaluated LIVE by `cm_models.eval_descent_live` (`MESH = 32/185` default; on-demand `64/185`, `96/185`, `384/370`). Production descent laws (`H-PGF`, `H-EA`, `RET_BOOT`, `C_rep`) are explicitly declared exploratory STNFS hypotheses supported (not proved) by the live audits below.

## 1. Exact Rank-2 Search Set `S_rank2(R)` -> `cells_rank2` & Zero Unit Duplicates

For special-`q` lattice `L_q` (`det = q`, `lll_16d_int(99, 100)`) and radius `R in {8..160}`, let `B_8(R) = { v in Z^8 : ||v||_2^2 <= R^2 }`, `F_Vor / mu_30` be the strict Voronoi fundamental domain mod `O_K^*` (`f_vor = 0.2681`, `f_vor/30 = 0.008937 > 0.008700`), and `N_alpha(R) = max(1, floor((|B_8(R)| - 1) * 0.008700))`. `dsieve16` enumerates canonical `alpha in (B_8(R)\{0}) cap (F_Vor/mu_30)` in lexicographic order, **stops at exactly `N_alpha(R)` vectors** (`A_cap(R)`), and pairs every `alpha in A_cap(R)` with every `beta in B_8(R)`:

    S_rank2(R) = A_cap(R) x B_8(R),   |S_rank2(R)| = N_alpha(R) * |B_8(R)|.

- **Exact live cell count (`cells_rank2`, exact for BOTH yield & full-ball work):** `audit_cm.py` computes `|B_8(R)|` LIVE via exact integer `theta_1(R^2)^8` Kronecker squaring for ALL `23` radii `R <= 160` (zero continuous approximation; e.g. `|B_8(76)| = 4,519,285,882,892,305` exact), asserts `cells_rank2[R] == N_alpha(R) * |B_8(R)|` and `cells_rank2[R] % |B_8(R)| == 0` for all `23` radii, AND enumerates `R = 2` LIVE (`|B_8(2)| - 1 = 1,712`, `|F_Vor|/30 = 16.87 >= N_alpha(2) = 14`; `|(B_8(R)\{0}) cap (F_Vor/mu_30)| >= N_alpha(R)` declared in `H5`).
- **Zero unit duplicates (proved) & primitivity:** `(alpha', beta') = u(alpha, beta)` with `alpha, alpha' in F_Vor/mu_30` forces `u = 1`. Non-unit multiples are rejected via `gcd_O(alpha, beta) == (1)`, retaining density `1/zeta_{Q(zeta_15)}(2) = 0.978453 > udup_node = 0.95` (computed LIVE).

## 2. `H-PGF (0.5 * mu)` & `THRESHOLD-REPAIR` (ZERO Count/Solvability Independence!)

Let `M` be distinct qualifying relations in `S_rank2(R)` before child descents, and `mu` the raw Dickman/Murphy/retention/detection scale (`sieve_mult = 1.0`, `MONT5 = 270`, `LLL_NODE = 2^40`, `100%` survivor/ECM/100-round-MR work charged).
- **`H-PGF (0.5 * mu)` exploratory envelope:** `mu` is NOT asserted equal to `E[M]`, and no exact NB distribution is claimed. For `t in [0, 0.80]`, `H-PGF` declares `E[t^M] <= G(t) = (1 + 0.5*mu*(1-t)/r)^(-r)` (`r = 0.5` if `kA + kB <= 1` else `1.0`).
- **`THRESHOLD-REPAIR` (`phi in [0.01, 0.04]`, `psi < 0.80`):** Sequential child failure `<= psi` gives `P(node fails) <= P(M < k) + psi^k <= G(t)/t^(k-1) + psi^k` without count/solvability independence. Minimizing over `k in 1..32, t in {0.05, 0.10, 0.20..0.80}` floored at `0.01` asserts `max_phi = 0.0189 <= 0.04`, `max_psi = 0.2060 < 0.80`. Node work: `own + child_W/(1-psi)`; invocations: `1 + child_EN/(1-psi)`.
- **Seeded Peeling (`h = 0.256845`, `exc_crit = 0.7625 << 1.09`) & Zero-Loss `repurge_stack`:** Verified LIVE in `hx_z12.py`; charges `exp(-miss)` over `[40, 70)`.

## 3. Fresh-Exponent `i.i.d.` Boot, Node Budget (`2^22.5`), Repeat Premise (`C_rep = 64`), & Meter (`t_run = 2^123.0`)

- **Fresh-Exponent `i.i.d.` Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `D_skew = 2^21`, `(n1=255, nA=195, kA=15)`):** Uniform `e_t in {0..r-1}`, `F_(p^12) -> F_(p^6) -> gamma_t`, 2-stage pure-int Bareiss `lll_16d_int` (`2^46.45 << LLL_BOOT = 2^48` Ir), check `|N(a_1^3 - u w^2 b_1^3)| <= 2^2560` (`RET_BOOT = 0.60`, declared exploratory MODEL inequality; `1/1` live from-scratch SHA-256 `G_T` compatibility check `<= 2544.7b`, `0` pre-supplied vectors). On child failure (`<= psi_boot = 0.2308 < 0.50`), draw fresh `e_t` (`i.i.d.`, zero `BALL_PAIR`): expected two-log work `2*(boot_own_exp + W_dag_1)/(1 - psi_boot) = 2^110.0863 <= 2^110.1343` (`96/185: 2^110.0743`, `384/370: 2^110.0384`).
- **Node Budget (`2^22.5`) & Declared `H5` Repeat-Collision Premise (`C_rep = 2^6 = 64`):** `E[N_2log] = 2^12.48`; aborting at `B_node = 2^22.5` gives `p_node_abort = E[N_2log]/B_node <= 9.64e-04`. Under declared exploratory `H5`, conditional on history each new child ideal of size `>= lx` equals any earlier ideal in the same descent with prob `<= C_rep / N_ideals(lx)` (`C_rep = 64 = 2^6`, matching `c36p`'s `1/p`-weighting factor over a 5-bit bin) -> `p_repeat_loss <= 64*B_node^2 / (2*N_ideals(lx)) <= 5.55e-04`, total descent loss `<= 0.0015527` (`ROUND_CEILING`).
- **Unified Hard Meter (`t_run = 2^123.0`):** Pointwise cap `<= 2^123.0 + 2^66`, `E[min(W_open, t_run)] <= 2^111.218513`, Markov `p_meter <= 0.00028407`.

## 4. LIVE Production `(f_9, g_3)` & Real-FB `Q(zeta_15)` Audits (ALL `lll_16d_int`, ZERO Pre-Supplied Vectors)

All audits are exploratory compatibility evidence (not proofs of production tails):
1. **`prod_sq24_live`:** ALL `24/24` production `(f_9, g_3)` `73..258`b control lattices reduced LIVE (`d_gso[16]==q^2`, `99/100`-Lovasz, `>=50b`/`>=25b` headroom) + `1,536` norms (`<= 4096` strip, `12..15`b): `0/24` empty, variance `4.0x–6.9x` tighter than `mu+mu^2`, `PGF_emp/PGF_ref <= 0.85`.
2. **`prod_cond_q15_live`:** ALL `12/12` genuine encountered `(f_9, g_3)` children (`6` g `120..255`b, `6` f `71..88`b) verifying **`12/12` LIVE `q | N_side(parent_a, parent_b)`** + reduced LIVE (`lll_16d_int`), ALL `6/6` sextiles `< qbins`.
3. **`q15_samefield_descent_live`:** Real `B_FB = 2^14` 2-core FB on `Q(zeta_15)` `(f_3, g_3)` (`3,728`-cell balls, exact `Z[omega]` `G_12`-orbit separation): replays ALL `3/3` witness balls LIVE from scratch (`exp_m`, `exp_sv` bit-exact) and verifies the stored `120`-ideal campaign (`1/120` empty `0.0083`, `mean_M = 4.69 > mu = 3.53`, `od = 0.065 << 1.0`, `MGF(0.25) = 0.366 << H-PGF(0.5*mu)`).

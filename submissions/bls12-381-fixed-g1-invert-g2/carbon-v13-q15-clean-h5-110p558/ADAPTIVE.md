# `ADAPTIVE.md` — Pure `k = 1` Boole + Wald `i.i.d.` Fresh-`e` Boot (`H5`) & Complete Zero-Omission Strict `k = 1` `Q(zeta_15)` Diagnostics

## 1. Pure `k = 1` Zero-Retry Boole + Exact Wald `i.i.d.` Fresh-`e` Boot (`pure_k1_boole_wald_iid_boot`)

Production BLS12-381 `Q(zeta_15)` descent uses **pure `k = 1` zero-retry Boole descent** per root attempt (`39` bins / `78` nodes, `N_tree_max_1log = 538` deterministic slots via exact integer bit-sum knapsack), wrapped in **unbounded `i.i.d.` fresh-`e` boot retries** (Wald's identity) and ONE unified hard aborting work meter `t_run = 2^{123.0}` units (`p_node_abort = p_repeat_loss = 0` proved):

1. **Zero-Retry `k = 1` DAG (`p_node_abort = p_repeat_loss = 0` proved, NO Gate E, NO `THRESHOLD-REPAIR` `psi^k`):**
   - Each of the `78` node plans (`70..265`b, step-2 `R in [68, 128]`, `min_mu = 555.04` [default `32/185` mesh] / `605.48` [`384/370` mesh], `max_p0 = 0.00359` [`32/185`] / `0.00329` [`384/370`]) scans its isotropic ball `S_16(R) = B_16(R)/{+-1}` once via complete zero-skip iterative Q32 `dsieve16` over ALL relative-degree-1 prime ideals `16 <= N(Q) <= 2^{70}` and takes its first usable relation (`k = 1`); `M = 0` (or an unsaturated `obpair` lattice) aborts the root attempt immediately and draws a fresh uniform boot exponent `e in {0..r-1}` (ZERO internal retries, ZERO `psi^k`, ZERO after-failure search).
   - **Pure Exact-Integer `O_H`-Submodule Saturation (`obpair.py` ONLY, ZERO Floating-Point Fallback) & Zero Fixed Divisors:** After exact Bareiss `lll_16d_int(99, 100)`, `obpair.ob_pair(q, ry, rt, R)` scans `R[1..15]` in pure exact integer arithmetic (`< 2^{18} << LLL_NODE = 2^{40}` Ir — ZERO floating-point fallback) to form an exact saturated `O_H`-basis `v_1 = (a_1, b_1), v_2 = (a_2, b_2)` satisfying **`|N_{H/Q}(det_{O_H}(v_1, v_2))| == q`** (`M == L_q` identically; succeeds `36/36` across all live production bases, `100%` of live toy nodes). If no `R[j]` achieves `|N(det)| == q`, the attempt deterministically aborts and draws a fresh boot `e`. Exhaustive live verification in `check_algebra_bls12` proves `min N(P) = 16 > 9 = deg(f_9) > 3 = deg(g_3)` across all prime ideals of `O_H = Z[zeta_15]`, unconditionally ruling out ANY fixed prime-ideal divisor `P != Q` on `M == L_q`.
   - **Explicitly Declared Exploratory STNFS Premise `H5`:** For every reached DFS slot `(j, b)` (`|T| <= 538` deterministically via bit-sum knapsack), we explicitly declare exploratory STNFS premise **`H5`** that reached-slot usable-relation counts (including any unsaturated-`obpair` abort) satisfy
     $$P\bigl(\text{slot } w_j \text{ fails} \mid A_{j,b}, Q_{\text{boot}}\bigr) \le p_0(b) := \max\Bigl(0.001,\; \frac{1}{1 + 0.5\,\mu(b)}\Bigr),$$
     which is implied via exact one-sided Cantelli by $\mathrm{Var}(M \mid A_{j,b}, Q_{\text{boot}}) \le \max(2\mu, 0.001\,\mu^2)$ — **explicitly declared as an exploratory STNFS premise rather than deduced unconditionally from mean yield or lattice geometry alone**.
   - **Exact 4-Line Coarse-Bin Linearity-of-Expectation Boole Bridge:** Indexing realized DFS slots `w_1, ..., w_{|T|}` (`|T| <= 538` deterministically) and coarse bin events `A_{j,b} := {j <= |T| and w_j in bin b}`:
     1. Pointwise Boole: `1_{E_fail} <= sum_{j=1}^{538} sum_{b=1}^{78} 1_{A_{j,b}} * 1_{fail_{w_j}}`.
     2. Exact linearity + definition of conditional probability: `P(E_fail | Q_boot) <= sum_{j,b} P(A_{j,b} | Q_boot) * P(fail_{w_j} | A_{j,b}, Q_boot)`.
     3. Under declared exploratory premise `H5`, `P(fail_{w_j} | A_{j,b}, Q_boot) <= p_0(b)`.
     4. Reverse linearity + deterministic pointwise bit-sum knapsack: **`P(E_fail | Q_boot) <= P_boot_tree = 0.56621`** (`32/185` mesh; **`0.55310`** on `384/370` mesh, both `< 0.60`!).
2. **Wald `i.i.d.` Fresh-`e` Boot (`boot nA = 195.0, n1 = 245, kA = 17`, Proved `RET_BOOT = 1.00`, ZERO `B_max`/`T_boot_max`):**
   - Each root attempt draws fresh uniform `e in {0..r-1}`, runs exact Bareiss `lll_16d_int(99, 100)`, checks the qualifying boot condition (`<= 17` primes in `(2^70, 2^195]`, `<= 1` in `(2^195, 2^245]`, `sum lo <= 2592`, `ecm_success_boot = 0.9999`), and executes the single-scan `k = 1` tree. Exact `i.i.d.` Wald gives default `32/185` **`E[two_logs] <= 2^{111.426585}`** u (`384/370` mesh: **`2^{111.363817}`** u; `boot_own_exp_log2 = 107.640215`, `tree_max_two_log2 = 109.957900`).

## 2. Scoped Multi-Regime Empirical Diagnostics (All Recomputed LIVE in `audit_cm.py`)

Finite prototypes and modeled production lattice checks are scoped exploratory diagnostics, NOT unconditional proofs or full-scale descents:

### Pillar 1 — True Full-Depth Squarefree Strictly-Decreasing Zero-Recount Multi-Level `Q(zeta_15)` `(3,3)` Adaptive Descent (`q15_adaptive_h5.b85` + Live 32-Root Full-Depth Strict `k = 1` Run)

Built on the `Q(zeta_15)` `(3,3)` tower (`f_3(a,b) = a^3 - (2 + w)b^3`, `g_3(a,b) = a^3 + 2 w^2 b^3`, `B_FB = 4,096`, `ball_limit = 128`, exact Bareiss `lll_16d_int(99, 100)` + exact-integer saturated `obpair` `|N(det)| == q`):
- **100% Complete Squarefree Factorization, Unique Ideal Enumeration, & True Full-Depth Termination:** Every residual cofactor `n_orig` after stripping `<= B_FB = 4,096` is factored **completely** (`assert math.prod(facs) == rem`) and ANY candidate relation containing ANY prime factor `>= q_parent`, ANY non-split prime `> 4,096` (`p % 15 != 1`), OR ANY repeated prime `p^2 | rem` (`len(facs) != len(set(facs))`) is **rejected outright**; because `v_p(rem) == 1` for every residual prime `p > B_FB`, scanning all 8 embeddings of `Q(zeta_15)` and asserting **`len(matches) == 1`** unconditionally proves that `100%` of contributing prime ideals above `B_FB` are enumerated (`prod_{P | (a-bX), N(P) > B_FB} N(P) == rem` identically — ZERO omitted ideals or multiplicities!). Recursion continues across ALL depths until **every leaf satisfies `len(ch) == 0` (`100%` of factors `<= B_FB = 4,096`)** before counting a root as completed.
- **ALL `q == 1 mod 15` Covered (Including `q == 1 mod 9`):** Uses exact Adleman–Manders–Miller / Tonelli–Shanks cube root mod `q` for ALL primes `q == 1 mod 15` (no `q % 9 in {4, 7}` restriction!).
- **Strict Norm Descent & ZERO Ancestor Recounting:** Every reached child strictly satisfies `4,096 < p_child < q_parent`, and every ancestor relation `(nf_anc, ng_anc)` is **explicitly excluded** (`(nf, ng) not in anc_rels`) when scanning child nodes.
- **Exact `k = 1` Zero-Retry Full-Depth Partition & Fresh Calibration (`469` fresh uniform split ideals across all `12` `(side, bit)` strata `b = 13..18`):** Across **`300` fresh uniform roots** (`q in (2^{16}, 2^{18})`, both sides), **`r_ok = 259 == 106 completed + 153 empty-slot aborts`** (true full-depth termination up to `max_depth = 6`, exact `k = 1` partition asserted live!), reaching **`514` adaptively reached child slots** (`290` d1 + `224` d2..d6, `321` distinct ideals, `153` observed empty-slot aborts), ALL recomputed from raw rows in `audit_cm.py`:

| Stratum | Reached slots `n` | Unique ideals | Obs. fails (`rate`) | **`O/E` vs `p_0(b)`** | Reached / fresh fail | Reached / fresh mean | Exact **`C_bday`** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `d1_all` | `290` | `229` | `77` (`0.266`) | **`0.5297`** (`< 1`) | `1.2330` | `0.9249` | **`1.171`** |
| `d1_g` | `181` | `137` | `48` (`0.265`) | **`0.5447`** (`< 1`) | `1.3521` | `0.8971` | **`1.171`** |
| `d1_f` | `109` | `92` | `29` (`0.266`) | **`0.5065`** (`< 1`) | `1.0761` | `0.9783` | **`1.171`** |
| `d2p_all` (`d2..d6`) | `224` | `155` | `76` (`0.339`) | **`0.6585`** (`< 1`) | `1.4518` | `0.8007` | **`1.261`** |
| **Pooled `d1..d6`** | **`514`** | **`321`** | **`153` (`0.298`)** | **`0.5867`** (`< 0.60`) | **`1.3328`** | **`0.8725`** | **`1.230`** |

- **3 LIVE Bit-Exact Complete Depth-2 Witnesses Replayed From Scratch in `audit_cm.py`:** `[f, 156691] > [g, 8101] > [f, 6481]`, `[g, 218641] > [f, 20101] > [f, 5101]`, and `[f, 110281] > [g, 23671] > [g, 19141]` — every hop re-executes `lll_16d_int(99, 100)`, `obpair`, exact `norm_q15_z`, exact squarefree `_fch_strict_sf`, unique-ideal `_ch_root_unique`, ancestor exclusion, and verifies that the depth-2 node resolves **completely into `<= 4,096` FB (`len(cg2) == len(cf2) == 0`)**!
- **LIVE 32-Root Fresh Full-Depth Strict `k = 1` Complete Descent Run Inside `audit_cm.py`:** On every audit run, `verify_q15_adaptive_h5_live` executes a fresh 32-root true-recursive full-depth squarefree strict `k = 1` zero-recount descent live (`q > 2^{16}`, alternating `g`/`f`, incl. `q == 1 mod 9`, exact `obpair` saturation, exact `_fch_strict_sf`, unique `_ch_root_unique`, explicit `anc_rels` exclusion, recursing to `len(ch) == 0` on EVERY branch), asserting exact partition **`nonempty = 25 == completed(12) + empty_slot_aborts(13)`** (`max_depth = 3`), **`24` d1 + `13` d2+ = `37` reached slots**, `13` fails, and direct per-slot **`O/E_p0 = 0.6766 < 1.0`**!

### Pillar 2 — Direct Real-FB `Q(zeta_15)` `(9,3)` Lower-Tail Balls (`220` Balls, `q15_93_b16_live`)

All 6 subsets (`120` low-yield `B_FB = 2^{13}, 3LP <= 2^{24}` + `100` higher-yield `B_FB = 2^{15}, 3LP <= 2^{30}`) are live-verified in `audit_cm.py` (plus 1 live regenerated 9-hit ball):

| Group | `n` | `mean(M)` | `min(M)` | Empty (`obs` vs `p_0`) | `Fano` (`Var/E`) | `k_hat` (`max(0,(V-E)/E^2)`) |
|---|---:|---:|---:|---:|---:|---:|
| `lo_g40` | `40` | `1.93` | `0` | `7/40` (`0.175 < 0.537`) | `1.316` | `0.1640` |
| `lo_f40` | `40` | `3.58` | `1` | `0/40` (`0.000 < 0.406`) | `1.117` | `0.0328` |
| `lo_enc40` | `40` | `2.02` | `0` | `8/40` (`0.200 < 0.562`) | `1.076` | `0.0375` |
| `hi_enc40` | `40` | `14.65` | `2` | **`0/40`** (`0.000`) | `2.743` | `0.1190` |
| `hi_g30` | `30` | `15.63` | `4` | **`0/30`** (`0.000`) | `3.911` | `0.1862` |
| `hi_f30` | `30` | `25.33` | `6` | **`0/30`** (`0.000`) | `3.475` | `0.0977` |

### Pillar 3 — Full 381-Bit Production Lattice Modeled-Yield Headroom Diagnostics (`36/36` Live-LLL Exact-Saturated `obpair` Bases: `prod_cond_q15_live` + `prod_sq24_live`)

- **Encountered-Child Production Lattices (`prod_enc12`, `12/12`):** Every child verifies live parent divisibility, live Bareiss `99/100` LLL (`0` presupplied vectors), exact-integer `obpair` `|N(det)| == q` saturation (`12/12`), exact-integer Q58 basis energy `eta_bar <= +0.20` (`12/12`), sextile headroom `36.39..52.15`b below `qbins`, and modeled joint-smoothness yield **`183.49x` min / `595.80x` med** above charged `pr`.
- **Fresh Production Lattices Across `73b..258b` (`prod_sq24`, `24/24`):** All `12` `g`-side and `12` `f`-side bases (`u_f in [9.33, 24.04]`) clear live Bareiss `99/100` LLL (`0` presupplied vectors), exact-integer `obpair` `|N(det)| == q` (`24/24`), exact-integer `eta_bar <= +0.20` (`24/24`), sextile headroom `35.97..52.92`b, and modeled joint-smoothness yield **`192.89x..195.48x` min / `549.97x..668.18x` med** above charged `pr` (`sp_rho_mean = -0.0024` `g`, `-0.0227` `f`).
- **Pure `Decimal(50)` + Exact-`O_H`-Identity `unitdomain.py`:** All `8` Voronoi cone rays (`m_0 = 0.028378`) and all `9` Sinnott full-cyclotomic-unit-group generators are verified using 50-digit `Decimal` rational-Taylor arithmetic and exact integer `O_H` polynomial equality (`g == +- y^m u_0^a u_1^b u_2^c`, ZERO float tolerances).
- **Honest Scope:** Production lattice diagnostics evaluate modeled smoothness headroom rather than execute production descents; `H5` (`p_0(b) = max(0.001, (1 + 0.5*mu(b))^{-1})`) remains an explicitly declared exploratory STNFS premise (§1).

# `ADAPTIVE.md` — Direct `Z`-Basis `L_q` (`p_sat = 0`), Pure `k = 1` Boole + Wald `i.i.d.` Fresh-`e` Boot (`H5`), & Complete Zero-Omission Strict `k = 1` `Q(zeta_15)` Diagnostics

## 1. Pure `k = 1` Zero-Retry Boole + Direct `Z`-Basis `L_q` (`p_sat = 0.0` Proved) + Wald `i.i.d.` Fresh-`e` Boot (`H5`)

Production BLS12-381 `Q(zeta_15)` descent uses **pure `k = 1` zero-retry Boole descent** per root attempt (`39` bins / `78` nodes, `N_tree_max_1log = 538` deterministic slots via exact integer bit-sum knapsack), operating **directly on the exact Bareiss-LLL `16 x 16` `Z`-basis `R[0..15]` of `L_q`**, wrapped in **unbounded `i.i.d.` fresh-`e` boot retries** (Wald's identity) and ONE unified hard aborting work meter `t_run = 2^{123.0}` units (`p_node_abort = p_repeat_loss = 0` proved):

1. **Direct Exact `16 x 16` `Z`-Basis `R[0..15]` of `L_q` (`p_sat_abort = 0.0` PROVED, ZERO `ob_pair` Search, ZERO Floating-Point Fallback) & Zero Fixed Divisors:**
   - Every special-`q` lattice `L_q subset Z^{16}` is constructed directly as its canonical exact `16 x 16` HNF `Z`-basis (`det_Z(B_q) = q`) and reduced by exact integer Bareiss `lll_16d_int(99, 100)` (`assert dsq[16] == q * q`), producing an exact `16 x 16` unimodular `Z`-basis `R[0..15]` spanning **the entire lattice `L_q` with index `[L_q : L_q] = 1` unconditionally (`p_sat_abort_proved = 0.0` — ZERO restricted `R[1..15]` pair search, ZERO saturation failure mode, ZERO floating-point fallback!)**.
   - Every cell `v in S_16(R) = B_16(R)/{+-1}` directly maps to `(a, b) = v * R in L_q`; every relative-degree-1 divisibility sublattice `L_{q, P} = {v in Z^{16} : (v R)_a - r_P (v R)_b == 0 mod P}` is constructed directly on `R` via exact `F_p` linear algebra (`det = N(P)`). Exhaustive live verification in `check_algebra_bls12` proves `min N(P) = 16 > 9 = deg(f_9) > 3 = deg(g_3)` across all prime ideals of `O_H = Z[zeta_15]`, unconditionally ruling out ANY fixed prime-ideal divisor `P != Q` on `L_q`.
2. **Zero-Retry `k = 1` DAG (`p_node_abort = p_repeat_loss = 0` proved, NO Gate E, NO `THRESHOLD-REPAIR` `psi^k`):**
   - Each of the `78` node plans (`70..265`b, step-2 `R in [68, 128]`, **`min_mu = 888.99`** [default `32/185` mesh] / **`1,061.88`** [`384/370` mesh], **`max_p0 = 0.00224`** [`32/185`] / **`0.00188`** [`384/370`]) scans its isotropic ball `S_16(R) = B_16(R)/{+-1}` once via complete prime-and-prime-power iterative Q32 `dsieve16` over ALL relative-degree-1 prime ideals AND prime powers `Q^e` (`e >= 2, N(Q)^e <= 2^{70}`, including exact `128`-bit `2`-adic `g`-side `P_i^k | a` `k = 1..6`), applies the NC-guarded `BOUND_CELL_HI` (`p_nc = 2^{-24}` explicitly priced; `ret_node = 0.90` charged; `KERNELS.md` §5), and takes its first usable relation (`k = 1`); `M = 0` aborts the root attempt immediately and draws a fresh uniform boot exponent `e in {0..r-1}` (ZERO internal retries, ZERO `psi^k`, ZERO after-failure search).
   - **Explicitly Declared Exploratory STNFS Premise `H5`:** Because `p_sat_abort = 0` is proved unconditionally on the direct `Z`-basis `R[0..15]`, a reached DFS slot `(j, b)` fails ONLY when its ball contains `M = 0` retained usable relations (already discounted by `ret_node = 0.90` in `mu(b)`). We explicitly declare exploratory STNFS premise **`H5`** that for every reached DFS slot `(j, b)` (`|T| <= 538` deterministically via bit-sum knapsack),
     $$P\bigl(M_{w_j} = 0 \mid A_{j,b}, Q_{\text{boot}}\bigr) \le p_0(b) := \max\Bigl(0.001,\; \frac{1}{1 + 0.5\,\mu(b)}\Bigr),$$
     which is implied via exact one-sided Cantelli by $\mathrm{Var}(M \mid A_{j,b}, Q_{\text{boot}}) \le \max(2\mu, 0.001\,\mu^2)$ — **explicitly declared as an exploratory STNFS premise rather than deduced unconditionally from mean yield or lattice geometry alone**.
   - **Exact 4-Line Coarse-Bin Linearity-of-Expectation Boole Bridge:** Indexing realized DFS slots `w_1, ..., w_{|T|}` (`|T| <= 538` deterministically) and coarse bin events `A_{j,b} := {j <= |T| and w_j in bin b}`:
     1. Pointwise Boole: `1_{E_fail} <= sum_{j=1}^{538} sum_{b=1}^{78} 1_{A_{j,b}} * 1_{M_{w_j} = 0}`.
     2. Exact linearity + definition of conditional probability: `P(E_fail | Q_boot) <= sum_{j,b} P(A_{j,b} | Q_boot) * P(M_{w_j} = 0 | A_{j,b}, Q_boot)`.
     3. Under declared exploratory premise `H5`, `P(M_{w_j} = 0 | A_{j,b}, Q_boot) <= p_0(b)`.
     4. Reverse linearity + deterministic pointwise bit-sum knapsack: **`P(E_fail | Q_boot) <= P_boot_tree = 0.55967`** (`32/185` mesh; **`0.55191`** on `384/370` mesh, both `< 0.60`!).
3. **Wald `i.i.d.` Fresh-`e` Boot (`boot nA = 195.0, n1 = 245, kA = 17`, Proved `RET_BOOT = 1.00`, ZERO `B_max`/`T_boot_max`):**
   - Each root attempt draws fresh uniform `e in {0..r-1}`, runs exact Bareiss `lll_16d_int(99, 100)`, checks the qualifying boot condition (`<= 17` primes in `(2^70, 2^195]`, `<= 1` in `(2^195, 2^245]`, `sum lo <= 2592`, `ecm_success_boot = 0.9999`), and executes the single-scan `k = 1` tree. Exact `i.i.d.` Wald gives default `32/185` **`E[two_logs] <= 2^{112.172395}`** u (`384/370` mesh: **`2^{112.135588}`** u; `boot_own_exp_log2 = 107.640215`, `dag_sub_two_log2 = 110.071109`).

## 2. Scoped Multi-Regime Empirical Diagnostics (All Recomputed LIVE in `audit_cm.py`)

Finite prototypes and modeled production lattice checks are scoped exploratory diagnostics, NOT unconditional proofs or full-scale descents:

### Pillar 1 — Direct-`R` True Full-Depth Squarefree Strictly-Decreasing Zero-Recount Multi-Level `Q(zeta_15)` `(3,3)` Adaptive Descent (`q15_adaptive_h5.b85` + `hi200` Lower-Tail + Live 32-Root Strict `k = 1` Run)

Built on the `Q(zeta_15)` `(3,3)` tower (`f_3(a,b) = a^3 - (2 + w)b^3`, `g_3(a,b) = a^3 + 2 w^2 b^3`, `B_FB = 4,096`, `ball_limit = 128` for main 300-root descent and `1,024` for `hi200`, operating **DIRECTLY on exact Bareiss-LLL `16 x 16` `Z`-basis `R[0..15]`** with `assert dsq[16] == q * q`):
- **100% Complete Squarefree Factorization, Unique Ideal Enumeration, & True Full-Depth Termination:** Every residual cofactor `n_orig` after stripping `<= B_FB = 4,096` is factored **completely** (`assert math.prod(facs) == rem`) and ANY candidate relation containing ANY prime factor `>= q_parent`, ANY non-split prime `> 4,096` (`p % 15 != 1`), OR ANY repeated prime `p^2 | rem` (`len(facs) != len(set(facs))`) is **rejected outright**; because `v_p(rem) == 1` for every residual prime `p > B_FB`, scanning all 8 embeddings of `Q(zeta_15)` and asserting **`len(matches) == 1`** unconditionally proves that `100%` of contributing prime ideals above `B_FB` are enumerated (`prod_{P | (a-bX), N(P) > B_FB} N(P) == rem` identically — ZERO omitted ideals or multiplicities!). Recursion continues across ALL depths until **every leaf satisfies `len(ch) == 0` (`100%` of factors `<= B_FB = 4,096`)** before counting a root as completed.
- **ALL `q == 1 mod 15` Covered (Including `q == 1 mod 9`):** Uses exact Adleman–Manders–Miller / Tonelli–Shanks cube root mod `q` for ALL primes `q == 1 mod 15`.
- **Strict Norm Descent & ZERO Ancestor Recounting:** Every reached child strictly satisfies `4,096 < p_child < q_parent`, and every ancestor relation `(nf_anc, ng_anc)` is **explicitly excluded** (`(nf, ng) not in anc_rels`) when scanning child nodes.
- **Exact `k = 1` Zero-Retry Full-Depth Partition & Fresh Calibration (`469` fresh uniform split ideals across all `12` `(side, bit)` strata `b = 13..18`):** Across **`300` fresh uniform roots** (`q in (2^{16}, 2^{18})`, both sides), **`r_ok = 298 == 242 completed + 56 empty-slot aborts`** (true full-depth termination up to **`max_depth = 7`**, exact `k = 1` partition asserted live!), reaching **`978` adaptively reached child slots** (`408` d1 + `570` d2..d7, `450` distinct ideals, `56` observed empty-slot aborts), ALL recomputed from raw rows in `audit_cm.py`:

| Stratum | Reached slots `n` | Unique ideals | Obs. fails (`rate`) | **`O/E` vs `p_0(b)`** | Reached / fresh fail | Reached / fresh mean | Exact **`C_bday`** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `d1_all` | `408` | `305` | `17` (`0.042`) | **`0.1393`** (`< 1`) | `1.2700` | `0.8995` | **`1.140`** |
| `d1_g` | `232` | `165` | `10` (`0.043`) | **`0.1505`** (`< 1`) | `1.6213` | `0.8607` | **`1.190`** |
| `d1_f` | `176` | `140` | `7` (`0.040`) | **`0.1259`** (`< 1`) | `0.9698` | `0.9586` | **`1.053`** |
| `d2p_all` (`d2..d7`) | `570` | `253` | `39` (`0.068`) | **`0.2074`** (`< 1`) | `1.1882` | `0.8602` | **`1.433`** |
| **Pooled `d1..d7`** | **`978`** | **`450`** | **`56` (`0.057`)** | **`0.1806`** (`< 0.25`) | **`1.2119`** | **`0.8781`** | **`1.284`** |

- **Direct-`R` Higher-Yield `200`-Ideal Lower-Tail Campaign (`hi200_lower_tail`, `ball_limit = 1,024`):** `80` fresh control (`mean_M = 24.62`, `0/80` zeros) + `80` depth-1 reached (`mean_M = 14.50`, `0/80` zeros) + `40` depth-2 reached (`mean_M = 9.72`, `1/40` zero) -> **`1/200` total zeros** and empirical Laplace lower-tail diagnostic **`E_emp[0.8^M] = 0.10466 << 0.2888`**.
- **3 LIVE Bit-Exact Complete Direct-`R` Depth-2 Witnesses Replayed From Scratch in `audit_cm.py`:** `[g, 160141] > [g, 39241] > [g, 4261]`, `[g, 253501] > [f, 139201] > [g, 4441]`, and `[g, 253501] > [f, 192931] > [g, 8431]` — every hop re-executes `lll_16d_int(99, 100)` (`assert dsq[16] == q * q`), direct `R`-vector evaluation, exact `norm_q15_z`, exact squarefree `_fch_strict_sf`, unique-ideal `_ch_root_unique`, ancestor exclusion, and verifies that the depth-2 node resolves **completely into `<= 4,096` FB (`len(cg2) == len(cf2) == 0`)**!
- **LIVE 32-Root Fresh Direct-`R` Full-Depth Strict `k = 1` Complete Descent Run Inside `audit_cm.py`:** On every audit run, `verify_q15_adaptive_h5_live` executes a fresh 32-root true-recursive direct-`R` full-depth squarefree strict `k = 1` zero-recount descent live (`q > 2^{16}`, alternating `g`/`f`, incl. `q == 1 mod 9`, exact `_fch_strict_sf`, unique `_ch_root_unique`, explicit `anc_rels` exclusion, recursing to `len(ch) == 0` on EVERY branch), asserting exact partition **`nonempty = 32 == completed(27) + empty_slot_aborts(5)`** (`max_depth = 6`), **`43` d1 + `50` d2+ = `93` reached slots**, `5` fails, and direct per-slot **`O/E_p0 = 0.1647 < 1.0`**!

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

### Pillar 3 — Full 381-Bit Production Lattice Modeled-Yield Headroom Diagnostics (`36/36` Live-LLL Direct-`Z`-Basis `R[0..15]` Bases: `prod_cond_q15_live` + `prod_sq24_live`)

- **Encountered-Child Production Lattices (`prod_enc12`, `12/12`):** Every child verifies live parent divisibility (`12/12`), live Bareiss `99/100` LLL (`0` presupplied vectors), direct `Z`-basis exact determinant `dsq[16] == q^2` (`12/12`, **`p_sat_abort_proved = 0.0`**), sextile headroom **`38.88`b `f` / `21.55`b `g`** below `qbins`, and modeled joint-smoothness yield **`20.87x` min / `230.64x` med** above charged `pr`.
- **Fresh Production Lattices Across `73b..258b` (`prod_sq24`, `24/24`):** All `12` `g`-side (`u_f in [10.36, 24.68]`) and `12` `f`-side (`u_f in [9.38, 20.43]`) bases clear live Bareiss `99/100` LLL (`0` presupplied vectors), direct `Z`-basis `dsq[16] == q^2` (`24/24`, **`p_sat_abort_proved = 0.0`**), sextile headroom **`43.47`b `f` / `25.11`b `g`** (`g`-side) and **`40.31`b `f` / `25.39`b `g`** (`f`-side), and modeled joint-smoothness yield **`56.07x` min / `274.87x` med** (`g`-side) and **`104.99x` min / `371.34x` med** (`f`-side) above charged `pr` (`sp_rho_mean = 0.0087` `g`, `0.0149` `f`).
- **Pure `Decimal(50)` + Exact-`O_H`-Identity `unitdomain.py`:** All `8` Voronoi cone rays (`m_0 = 0.028378`) and all `9` Sinnott full-cyclotomic-unit-group generators are verified using 50-digit `Decimal` rational-Taylor arithmetic and exact integer `O_H` polynomial equality (`g == +- y^m u_0^a u_1^b u_2^c`, ZERO float tolerances).
- **Honest Scope:** Production lattice diagnostics evaluate modeled smoothness headroom rather than execute production descents; `H5` (`p_0(b) = max(0.001, (1 + 0.5*mu(b))^{-1})`) remains an explicitly declared exploratory STNFS premise (§1).

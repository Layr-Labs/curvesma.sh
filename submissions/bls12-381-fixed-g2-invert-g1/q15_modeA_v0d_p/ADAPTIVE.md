# `ADAPTIVE.md` — Pure `k = 1` Boole + Wald `i.i.d.` Boot & Raw-`L_q` Exploratory Diagnostics

## 1. Pure `k = 1` Zero-Retry Boole + Exact Wald `i.i.d.` Fresh-`e` Boot (`pure_k1_boole_wald_iid_boot`)

Production BLS12-381 `Q(zeta_15)` CM-STNFS (`MODEA_V0_ORIG_RADII`) uses **pure `k = 1` zero-retry Boole descent** per root attempt, wrapped in **unbounded `i.i.d.` fresh-`e` boot renewals** (Wald's identity) and ONE unified hard aborting run meter `t_run = 2^{124}` units:
1. **Zero-Retry `k = 1` DAG (`p_node_abort = p_repeat_loss = p_descent_fail = 0`, Declared Exploratory Premise `H5`):** Each of the `76` single-special-`q` raw-`L_Q` rows (`LB = LX = 68.5`, `R <= 112`) scans its ball once and takes its first usable relation (`k = 1`); `M = 0` aborts the root attempt immediately (ZERO internal retries, ZERO after-failure search). During the already-charged `BOUND_CELL = 5,000` Ir/cell pass, any node failing basis/norm screening is discarded immediately as an empty slot (`M = 0`). The negative-binomial bound `p_0(b) = max(0.001, (1 + 0.5*mu(b))^{-1})` (`NB_R = 1.0`, `NB_ALPHA = 0.5`, `P0_FLOOR = 0.001`, **`mu_min = 701.84`**, **`p0_max = 0.002842`**) is a **declared exploratory-v1 STNFS premise (`H5`)** bounding total per-slot empty-or-discarded probability at every reached slot (including children selected by the first usable parent relation and any node discarded by online `BOUND_CELL` screening). Under `H5`, unconditional Boole gives **`Phi_tree = 0.552788 <= 0.60`** over **`c_det_one = 512`** deterministic slots.
2. **Wald `i.i.d.` Fresh-`e` Boot (`BOOT16 = (230, 288, 15, 2)`, ZERO `B_max`/`T_boot_max`):** Each root attempt draws fresh uniform `e in [0, r-1]` and computes `z' = z * gamma_0^e`. **If `z' == 1`, return `V(z) = (-e * V_0) mod r` immediately with ZERO descent work** (verified live in `audit_cm.py`). Only on nonidentity `z' in G_T \ {1}` does it invoke Guillevic `F_{p^6}`-subfield elimination (`W_1 != 0` unconditionally on `G_T \ {1}`), run fresh 2-stage `11/40 -> 99/100` `Q_CM` Bareiss 16D LLL (`<= 2^{47}` Ir; input-uniform `2^{46.84}` proved), and sample `1` uniform half-cell of `B_16(2)/{+-1}` (`BALL_PAIR = 0.0` exact, `BOOT_KEEP_FRAC = 0.90 < CP95 LB = 0.9395`), giving **`E[runs/log] <= 1/(1 - Phi_tree) = 2.236077`**, **`E[trials/log] <= 3.2839e+16`**, **`E[two_logs] <= 2^{111.164668}`** (charged `2^{111.1648}`, `boot_u <= 2^{109.0872}`, `sub_u <= 2^{110.7745}`).

## 2. Raw-`L_q` Exploratory Descent Diagnostics (`q15_mini_adaptive.b85`, `verify_q15_descent.py`)

*Scope:* Scoped finite-scale and norm-model exploratory diagnostics supporting declared exploratory-v1 STNFS premise `H5`.

### (A) Production `deg-(72, 24)` Raw-`L_q` 38-Bin Norm Survey (`30,400` Ideals) & Genuine Unplanted Reached Production Children (`416` Children)
- **Full 38-bin raw-`L_q` survey (`[70..259]` bits, exact `11/40 -> 99/100` `Q_CM` Bareiss LLL, `64` exact `B_16(R)` cells/ideal = `1,945,600` cells, `+0.63`b hardware envelope):** `30,400` genuine production `deg-(72, 24)` raw-`L_q` prime ideals (`800/bin`, `400 f_9 + 400 g_3`), `0` `v_q >= 2` basis failures, `min ||b_0||_2 / q^{1/16} = 0.69978`, `min ||b_0||_{Q_CM} / (det(G_8)^{1/16} q^{1/16}) = 0.73832`, **`0 / 30,400` below `1.0x`** (`min = 3.252x`, `p01 = 8.48x`, `median = 45.45x`).
- **Genuine unplanted reached production raw-`L_q` children (`416` children vs `416` uniform controls):** Harvested from `32` unplanted raw-`L_q` parents across 4 parent bins `[100,105)`, `[145,150)`, `[185,190)`, `[225,230)` (`112 f_9 + 304 g_3`, bits `[70..98]`, `64` cells/child): `0 / 416` below `1.0x`, `min = 11.475x` (`median = 74.87x` vs control `73.34x`), **`geom_mean(child / ctrl) = 1.0135`** (95% parent-clustered bootstrap `[0.9319, 1.1002]`), **`KS 2-sample p = 0.9702`**.
- **Live exact-`R` `BOUND_CELL`-screened spot-check (`verify_q15_descent.py`):** `4` live raw-`L_q` production `deg-(72, 24)` lattices across bins `(68.5f, 100g, 180g, 230f)` evaluated at their **exact charged radii `R in {56, 80, 80, 14}`** (`32` primitive cells/lattice, zero radius rescaling, every cell screened through `BC.bound_cell_check_and_cap` with exact `(nf_cap, ng_cap)` and `+0.625`b hardware envelope): `32/32` kept on all 4 bins, `0/4` below `1.0x`, `screened_cert_ratio = 25.4x..60.3x`.

### (B) Full-16D Raw-`L_q` Pure `k = 1` Parent & Reached-Child Cohort (`31,000` Archived + Live Reproducible Generator in `verify_q15_descent.py`)
Exact `Q(zeta_15)` `deg-(32, 24)` tower (`f_4 = x^4 + x + w`, `g_3 = x^3 + w^2`), `B_fb = 2^{17} = 131,072 < LP_max = 2^{22} - 1 < q_parent in [2^{22}, 2^{23})`, **FULL 16D** `2,496`-pt `||k||_2^2 <= 3` half-balls, pure `k = 1` raw `L_q` (`11/40 -> 99/100` `Q_CM` Bareiss LLL):
- **Live reproducible generator executed inside `audit_cm.py` (`_live_full16d_k1_deg32_24`):** Samples `20` fresh `deg-(32, 24)` parents (`10 f_4 + 10 g_3`), extracts the exact first-usable (`k = 1`) reached child ideal from each parent's full-16D `2,496`-pt ball, and runs pure `k = 1` full-16D descent on every reached child (`>= 32` live ideals total, **`0` empty balls**).
- **Archived `31,000`-ideal cohort (`q15_mini_adaptive.b85`):** `16,000` raw-`L_q` `k = 1` parents (`0 / 16,000` empty, `mean M = 40.92`, `within_bin_fano = 1.2276`, `min M = 8`) + `15,000` distinct `k = 1`-reached children (`0 / 15,000` empty, `mean M = 31.91`, `within_bin_fano = 1.1359`, `min M = 6`), combined **`0 / 31,000` empty** (`within_bin_fano = 1.1888`).

### (C) Real-2-Core-FB Multi-Level `Q(zeta_15)` Campaign (`B_FB = 2^{14}`, `12,802` Roots + `3` Live Bit-Exact Witnesses)
Miniature `Q(zeta_15)` tower (`deg-(32,24)` over `O_H = Z[zeta_15]`, `u_t = 1, p_t = 7, g_t = X^3 - w^2, f_t = X^4 + 3X + (w - 2)`), real 2-core FB `B_FB = 2^{14} = 16,384` from `31,657 collected -> 31,641 2-core (454/458 f-split, 431/450 g-split)` genuine `Q(zeta_15)` relations, `576`-pt `alpha = 1` 8D slice, measuring **both** `k = 1` first attempts (`10,281 / 6,251`) and depth-2 after-failure child encounters (`6,393 / 4,010`):

| Setting | Roots | 1st att / fail | After att / fail | `O/E fail` | `O/E work` | Sh prev / run |
|---|---:|---:|---:|---:|---:|---:|
| `f19_d1` | `1800` | `1512 / 882` | `991 / 572` | `0.8147` | `1.0143` | `2 / 3` |
| `g19_d1` | `1702` | `1234 / 672` | `476 / 240` | `0.7757` | `1.0365` | `0 / 0` |
| `f20_d1` | `2500` | `2201 / 1363` | `1761 / 1115` | `0.8336` | `0.9919` | `0 / 3` |
| `g20_d1` | `2500` | `1911 / 1081` | `892 / 543` | `0.8955` | `0.9667` | `1 / 1` |
| `f19_d2` | `1800` | `1512 / 1009` | `1204 / 799` | `0.7976` | `1.0214` | `2 / 2` |
| `g20_d2` | `2500` | `1911 / 1244` | `1069 / 741` | `0.8410` | `0.9723` | `1 / 1` |
| **Pooled** | **`12802`** | **`10281 / 6251`** | **`6393 / 4010`** | **`0.8288`** | **`0.9961`** | **`6 / 10`** |

- Child concentration across bins: `max C_emp = 20.49`, **`max C_bday = 2.4897`**, shared collisions `6 / 10` (`0.1564%`).
- **3 live bit-exact multi-level witnesses** (`_lll_16d_qcm`, replayed live in `verify_q15_descent.py`): `f19_d1 (side=f, q=265891, rt=93500, seq=[[0, 1, 1, 1, 0], [1, 0, 1, 1, 0]])`; `g19_d1 (side=g, q=265861, rt=11390, seq=[[0, 1, 1, 1, 0], [1, 0, 1, 1, 0]])`; `f20_d1 (side=f, q=525871, rt=425999, seq=[[0, 1, 1, 1, 0], [1, 1, 1, 1, 0]])`.

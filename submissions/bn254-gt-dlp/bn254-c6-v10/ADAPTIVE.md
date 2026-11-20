# `ADAPTIVE.md` — Pure `k = 1` Boole + Wald `i.i.d.` Boot & Live Strict `k = 1` `C6` `O_h` Campaign

## 1. Pure `k = 1` Zero-Retry Boole + Exact Wald `i.i.d.` Fresh-`e` Boot (`pure_k1_boole_wald_iid_boot`)

Production BN254 `C6` `v10` uses **pure `k = 1` zero-retry Boole descent** per root attempt, wrapped in **unbounded `i.i.d.` fresh-`e` boot retries** (Wald's identity) and ONE unified hard aborting run meter `t_run = 2^{112}` units:
1. **Zero-Retry `k = 1` DAG (`p_node_abort = p_repeat_loss = p_descent_fail = 0`):** Each of the `62` rows (`plan59`, with all `27` rows `R <= 64` bumped `+1` step for stronger `H5` margin) scans its ball once and takes its first usable relation (`k = 1`); `M = 0` aborts the root attempt immediately (ZERO internal retries, ZERO `THRESHOLD-REPAIR` `psi^k`, ZERO after-failure search). Under declared exploratory premise `H5` (that reached-slot usable-relation counts satisfy `P(M = 0 | slot) <= p_0(b) = max(0.001, (1 + 0.5*mu(b))^{-1})`, matching a Cantelli second-moment overdispersion bound `Var(M|slot) <= 2*mu + 0.001*mu^2` — explicitly disclosed as an exploratory premise rather than deduced from mean yield alone; `mu_min = 760.80`, `p0_max = 2.622e-3`), unconditional Boole gives **`Phi_tree = 0.443226 <= 0.60`** over `c_det_one = 427` deterministic slots.
2. **Wald `i.i.d.` Fresh-`e` Boot (`boot = [180, 76, 12, 4]`, ZERO `B_max`/`T_boot_max`):** Each root attempt draws fresh uniform `e in [0, r-1]` (`T = phi_iso(z*g^e)`) and runs fresh 2-stage `11/40 -> 9/10` Bareiss 12D LLL (`<= 2^{40}` Ir in `trial_b`; input-uniform `2^{39.81}` proved), giving `E[runs/log] <= 1/(1-Phi_tree) = 1.796059`, `E[trials/log] <= 2.7693e16`, `E[two_logs] <= 2^{83.798589}` (`two_dag <= 2^{82.8094}`, `two_boot_exp <= 2^{82.7877}`, `K_h`-split boot root `2^{62.094} B << 2^{64} B`).
3. **Naive-PGF Counterexample (`verify_threshold_counterexample`):** `0.10*t + 0.90*t^{100} <= (1+30(1-t))^{-1}` on `[0, 0.80]`, yet `0.10 > G(0.10) = 0.03571`. Pure `k = 1` Boole avoids `G(psi)` and `psi^k`.

## 2. Live Strict `k = 1` Zero-Retry Prototype & Real-2-Core-FB `C6` `O_h` Campaign (`c6_adaptive.py`)

Miniature `C6` `O_h` tower (`u_t = 2, p_t = 13` PRIME, `h` irred mod `13`, `2*w` non-square in `F_{13^6}`), real 2-core FB `B_FB = 2^{14}` from `12,712` genuine `O_h` relations (`79/99` `f`, `87/99` `g` split orbits), live `bareiss_lll_12d(B, I_6, 1, 1)`, `P^1`-unique `484`-pt ball (`alpha = 1, beta in B_6(sqrt(4))`):
- **Expanded Live Strict `k = 1` Zero-Retry Multi-Level Execution (`_desc_k1`, ZERO retries after failure):** Directly executes the exact claimed `k = 1` zero-retry protocol (taking ONLY `rels[0]` at every reached node and aborting immediately on `M = 0`, ZERO subsequent-relation retries) LIVE inside `audit_bn_cm.py` on **48 genuine 17-bit split prime ideals** (`nw = 16`, `ne = 24/48` non-empty first relations, **`17/24 = 70.8%` strict `k = 1` completion** across multi-level recursion, **`14` genuine reached multi-level children**, `7` child fails), plus live strict `k = 1` first-choice replay (`ri = 0`) of the **3 multi-level witnesses** (`f19_d1: q=263491`, `g19_d1: q=263561`, `f20_d1: q=525947`).
- **Companion Conditional-Child Calibration Campaign (`c6_adaptive_h5.b85`, `9,870` roots, `4,887` first attempts, `1,856` conditional child checks across 6 settings):**

| Setting | Roots | 1st att / fail | Cond. att / fail | `O/E fail` | `O/E work` | Sh prev / run |
|---|---:|---:|---:|---:|---:|---:|
| `f19_d1` | `1,127` | `459 / 352` | `121 / 97` | `0.871` | `0.988` | `0 / 1` |
| `g19_d1` | `1,127` | `507 / 414` | `210 / 171` | `0.805` | `0.994` | `3 / 3` |
| `f20_d1` | `2,163` | `1,010 / 792` | `329 / 264` | `0.803` | `1.012` | `0 / 0` |
| `g20_d1` | `2,163` | `1,226 / 1,027` | `528 / 446` | `0.780` | `0.991` | `0 / 0` |
| `f19_d2` | `1,127` | `459 / 361` | `125 / 100` | `0.846` | `0.980` | `0 / 0` |
| `g20_d2` | `2,163` | `1,226 / 1,051` | `543 / 466` | `0.772` | `1.006` | `0 / 0` |
| **Pooled** | **`9,870`** | **`4,887 / 3,997`** | **`1,856 / 1,544`** | **`0.7935`** | **`0.9985`** | **`3 / 4`** |

- `6/6` settings have `O/E fail <= 0.871 < 0.90` (pooled `0.7935`), `10/10` child bins have **`max C_bday = 2.6198 << 256`**. Scoped `B_FB = 2^{14}` diagnostic evidence, not a proof.

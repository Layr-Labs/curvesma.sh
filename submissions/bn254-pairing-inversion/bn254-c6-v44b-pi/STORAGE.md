# STORAGE.md: Strict `< 16 EiB`/Node 64-Bit Single-Ball `0LP` `MODE-I` Architecture (`stor_bound.py`), Stage-3 FVS & Stage-4 32D-`K_f` Boot Working Sets, and Exact `< 2^{98}` B Machine Storage Certificate

Where all sieve, LA, split, and 32D-`K_f` boot data lives, how it is addressed strictly within 64-bit general-purpose pointers (`< 2^{64}` B `= 16 EiB` per node throughout), and how every pointwise storage bound is certified in exact rationals by `stor_bound.py` and checked live by `audit_bn_cm.py` (`ram_ok = True`, `storage_cap_abort = 0.0`). All `<= 2^{...}` ceilings round UP.

## 1. Strict `< 16 EiB`/Node Single-Ball `0LP` `MODE-I` Storage Certificate (`stor_bound.py`, `R_A_MAX = 572.0`)

Every active node holds strictly less than `2^{64}` B (`16 EiB`), so every load, store, and `STADDB`/`LDADD` executes natively on standard 64-bit ARMv8.1-A general-purpose registers (`-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`). Because production runs on the exact unmodified `0LP (61.50, 61.50)` single `E_{12,trunc}` ball (`R_a = 568.83 <= 572.0 = R_{A,\max}`, `sup_ng = 477.18 <= 478.0`), every count in `stor_bound.py` holds pointwise identically to `bn254-c6-v43b-pi`:

| Node / Working Set | Count | Holds | RAM (`log2` B) |
|---|---:|---|---:|
| Emitter `E_0` (`stor_bound.py`) | 1 | `a_tab`, 4-plane `B_TAB`, 1 TiB Q12 `isqrt` LUT, `rot_tab`, `2^{34}` x `256 MiB` epoch buffers (`1`-hop `emit8/9`), plus `< 1 GiB` code/stack/page tables (strictly enclosed inside the `5%` / `0.827 EiB` headroom below `16 EiB`) | `<= 63.9234` (`15.1725` EiB `< 16 EiB`) |
| Score Nodes `C[s]` (`stor_bound.py`) | `2^{34}` | `C[0]` dummy sink + `C[1..2^{34}-1]` `uint8_t` (`scale = 0.35`, 1 B/cell) score rows (`max_rows <= 224,263`; unused power-of-2 stride padding beyond each row's active `|I_a|` cells is never touched or zeroed, so init/scan charges only the `2^{94.8899}` active cells) + `32 MiB` Stage-B `g`-SRT (`L_g = 13`, `21,814,168` B) | `<= 62.9035` (`blocks <= 31,382,183 < 2^{25}`) |
| List Nodes `L[k]` (`stor_bound.py`) | `2^{36}` | Band 1 + dual `pos`/`neg` Lower-B2 prefix lists (`2^{46}` 16 KiB link pages, `B2_PREFIX_MULT = 2.8367 <= 2.84`) + build scratch | `<= 60.1617` (`model 60.1611`) |
| `f`-Cofac Node `F_0` (`stor_bound.py`) | 1 | Dedicated 1-owner `f`-side smoothness batch (`L_f = 50`, `2^{50}` 7-limb leaves: `P_f`, 51 `PT`/`RT` levels, `8x` NTT scratch) | `<= 62.7916` (`6.9239` EiB `<< 16 EiB`) |
| **Stage-3 `0LP` FVS Schur-BW LA (`bw=(1024,48)`)** | **1 owner** | **Unified SpMV (`2^{59.8926}` B, `2^{56}*nslots`) + PM-Basis (`2^{62.2970}` B) at `nu_fvs = 0.45, nu_d = 2.85, mdepth = 1.0`** | **`62.5466`** (`5.95` EiB `< 16 EiB`; exact-rat UB at actual `bw=(1024,48)`: **`62.6205`**) |
| **Stage-4 32D-`K_f` Boot Root (`LB = 58.0`)** | **1 owner** | **DF post-order all-primes root `|P| = 2^{56.8517}` B + `8x` Karatsuba/NTT3 peak (`<= 2^{59.8517}` B)** | **`59.8517`** (`0.90` EiB `<< 16 EiB`) |
| **Stage-2 Deterministic `ntt3_F` Split** | **1 owner** | **Bernstein prime-tree split on `(B_f, B_g) = (61.50, 61.50)` (`Q = 2^{46}`)** | **`<= 60.0`** (`<= 2^{59.5}` B actual `<< 16 EiB`) |
| **Pointwise Total Machine Storage** | — | **`stor_bound.py` exact-rational `MODE-I` machine envelope (`storage_cap_abort = 0.0`)** | **`<= 97.5798` (`< 100`)** |

## 2. Single-Worker Working Sets (`< 2^{64}` B)

1. **Stage-B `g`-side `P_min = 24` SRT (`L_g = 13`):** `21,814,168` B `<= 20.81 MiB <= 32 MiB` per Score Node. Dedicated `F_0` (`L_f = 50`): `<= 2^{62.7916}` B (`6.9239` EiB `<< 16 EiB`).
2. **Stage-3 100%-Matched `0LP` Zero-Loss Hybrid FVS Schur-BW (`bw = (1024, 48)`):** SpMV `2^{59.8926}` B (`D_rows*(512*nslots + 45*8) + 2^{56}*nslots + 65*2^{53}`) + PM-Basis `2^{62.2970}` B `=>` `2^{62.5466}` B (`5.95` EiB `< 16 EiB`); exact-rational `U_ub = 2^{54.62}` certificate at the ACTUAL `bw = (1024, 48)` and `w_eff_ub = 45` gives **`2^{62.6205}` B `>= 2^{62.5466}` B** (`6.26` EiB `< 16 EiB`), asserted live by `fvs_ladder_const_cap_and_16eib_ok`.
3. **Stage-2 Deterministic `ntt3_F` Bernstein Split:** `<= 2^{59.5}` B (reserved `2^{60}` B `<< 16 EiB`).
4. **Stage-4 Native 32D-`K_f` Boot Root (`x_EA = 59.3229`, `LB = 58.0 <= 61.50`):** Retained DF post-order root `|P| = 2^{56.8517}` B (`0.112` EiB), peak `8*|P| = 2^{59.8517}` B (`0.90` EiB `<< 16 EiB`; `v32_verify.py::boot_df_tree_ok` verifies DF peak `2^{59.4878}` B `= 6.22*|P| < 8*|P| < 2^{64}` at actual `x_EA = 59.3229` and `< 5*|P|` at the `62.9375` upper envelope; Stage-2 Kronecker-NTT3 / BLS MPE ECM `<= 601.97 MiB`).

## 3. Full Sensitivity Table (`part="sensitivity"` & `part="dag_survivors"`)

All scenarios below run live inside `audit_bn_cm.py` under strict single-ball `< 16 EiB`/node `0LP` `MODE-I`:

| Scenario | Description | `exp_run` (`log2 u`) | Claim (`log2 u`) | Delta vs Primary |
|---|---|---:|---:|---:|
| **Primary (`0LP`, `MODE-I`, `0.45/2.85`)** | Unmodified `0LP (61.50, 61.50)` sieve, 100%-matched `0LP` const-cap FVS, live-constructed 32D boot (`2^{43.5}` LLL, exact `B_main = 61.50`) | **`83.969589`** | **`82.611136`** | `0.0000` b |
| **`fvs_u19_030`** | `0LP` FVS at `nu_fvs = 0.30` (still `12.4x` above `U = 2^{22}` max `0.0242`) | `83.926928` | `82.568475` | `-0.0427` b |
| **`gate_0395_slack`** | Gate coin raised to `0.395` (`+0.005` success slack above `0.39`) | `83.969589` | `82.629515` | `+0.0184` b |
| **`zero_fvs_baseline_h3`** | Disable FVS completely (`32D-Kf` boot alone atop baseline `W_row = 100.5` H3 LA) | `84.004880` | `82.646429` | `+0.0353` b |
| **`stress_plus15pct_lat_and_cof`** | Simultaneous `+15%` stress penalty on BOTH `lat_b2` and `cof` | `84.078337` | `82.719884` | `+0.1087` b |
| **`dag_survivors` (`2x`)** | Double all DAG survivor rates (`H5-W` stress) | `83.971517` | `82.613064` | `+0.0019` b |
| **`boot_charge = 4` / `boot_surv x4`** | Halve boot yield / quadruple boot survivors (`H-EB0/H-EB1` stress) | — | `82.622408` / `82.615841` | `+0.0113` / `+0.0047` b |
| **`plain_dickman_rho`** | Unpadded plain Dickman `rho` (reference only) | `83.798690` | `82.440238` | `-0.1709` b |

# STORAGE.md: Strict `< 16 EiB`/Node 64-Bit Single-Ball `MODE-I` Architecture (`stor_bound.py`), Stage-3 FVS & Stage-4 32D-`K_f` Boot Working Sets, and Exact `< 2^{98}` B Machine Storage Certificate

Where all sieve, LA, split, and 32D-`K_f` boot data lives, how it is addressed strictly within 64-bit general-purpose pointers (`< 2^{64}` B `= 16 EiB` per node throughout), and how every pointwise storage bound is certified in exact rationals by `stor_bound.py` and checked live by `audit_bn_cm.py` (`ram_ok = True`, `storage_cap_abort = 0.0`). All `<= 2^{...}` ceilings round UP.

## 1. Strict `< 16 EiB`/Node Single-Ball `MODE-I` Storage Certificate (`stor_bound.py`, `R_A_MAX = 572.0`)

Every active node holds strictly less than `2^{64}` B (`16 EiB`), so every load, store, and `STADDB`/`LDADD` executes natively on standard 64-bit ARMv8.1-A general-purpose registers (`-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`). Because production runs on an unmodified single `E_{12,trunc}` ball (`T = 1`, zero CRT tiering) with **`R_a = 564.78 <= 572.0 = R_{A,\max}`** (`sup_ng = 477.05 <= 478.0`), every count in `stor_bound.py` (which is nondecreasing in `T_12`) holds pointwise without modification:

| Node / Working Set | Count | Holds | RAM (`log2` B) |
|---|---:|---|---:|
| Emitter `E_0` (`stor_bound.py`) | 1 | `a_tab`, 4-plane `B_TAB`, 1 TiB Q12 `isqrt` LUT, `rot_tab`, `2^{34}` x `256 MiB` epoch buffers (`1`-hop `emit8/9`), plus `< 1 GiB` code/stack/page tables (strictly enclosed inside the `5%` / `0.827 EiB` headroom below `16 EiB`) | `<= 63.9234` (`15.1725` EiB `< 16 EiB`) |
| Score Nodes `C[s]` (`stor_bound.py`) | `2^{34}` | `C[0]` dummy sink + `C[1..2^{34}-1]` `uint8_t` (`scale = 0.35`, 1 B/cell) score rows (`max_rows <= 224,263`) + `32 MiB` Stage-B `g`-SRT (`L_g = 13`, `21,814,168` B) | `<= 62.9035` (`blocks <= 31,382,183 < 2^{25}`) |
| List Nodes `L[k]` (`stor_bound.py`) | `2^{36}` | Band 1 + dual `pos`/`neg` Lower-B2 prefix lists (`2^{46}` 16 KiB link pages, `B2_PREFIX_MULT = 2.923 <= 2.94`) + build scratch | `<= 60.1617` (`model 60.1611`) |
| `f`-Cofac Node `F_0` (`stor_bound.py`) | 1 | Dedicated 1-owner `f`-side smoothness batch (`L_f = 50`, `2^{50}` 7-limb leaves: `P_f`, 51 `PT`/`RT` levels, `8x` NTT scratch) | `<= 62.7916` (`6.9239` EiB `<< 16 EiB`) |
| **Stage-3 FVS Schur-BW LA (`bw=(1024,40)`)** | **1 owner** | **Unified SpMV (`2^{59.7730}` B, `2^{56}*nslots`) + PM-Basis (`2^{62.3122}` B) at `nu_fvs = 0.45, nu_d = 2.70, mdepth = 1.7011`** | **`62.5412`** (`5.93` EiB `< 16 EiB`; exact-rat UB at actual `bw=(1024,40)`: **`62.8110`**) |
| **Stage-4 32D-`K_f` Boot Root (`LB = 58.0`)** | **1 owner** | **DF post-order all-primes root `|P| = 2^{56.8517}` B + `8x` Karatsuba/NTT3 peak (`<= 2^{59.8517}` B)** | **`59.8517`** (`0.90` EiB `<< 16 EiB`) |
| **Stage-2 Deterministic `ntt3_F` Split** | **1 owner** | **Bernstein prime-tree split on `(B_f, B_g) = (58.0, 61.5)` (`Q = 2^{46}`) + LP CZ root-finding** | **`<= 60.0`** (`<= 2^{59.5}` B actual `<< 16 EiB`) |
| **Pointwise Total Machine Storage** | — | **`stor_bound.py` exact-rational `MODE-I` machine envelope (`storage_cap_abort = 0.0`)** | **`<= 97.5798` (`< 100`)** |

## 2. Single-Worker Working Sets (`< 2^{64}` B)

1. **Stage-B `g`-side `P_min = 24` SRT (`L_g = 13`):** `21,814,168` B `<= 20.81 MiB <= 32 MiB` per Score Node. Dedicated `F_0` (`L_f = 50`): `<= 2^{62.7916}` B (`6.9239` EiB `<< 16 EiB`).
2. **Stage-3 Zero-Loss Hybrid FVS Schur-BW (`bw = (1024, 40)`):** SpMV `2^{59.7730}` B (`D_rows*(512*nslots + W_pre*mdepth*8) + 2^{56}*nslots + 65*2^{53}`) + PM-Basis `2^{62.3122}` B `=>` `2^{62.5412}` B (`5.93` EiB `< 16 EiB`); exact-rational `U_ub = 2^{54.62}` certificate at the ACTUAL `bw = (1024, 40)` and `w_eff_ub = 81` gives **`2^{62.8110}` B `>= 2^{62.5412}` B** (`7.15` EiB `< 16 EiB`), asserted live by `fvs_ladder_const_cap_and_16eib_ok`.
3. **Stage-2 Deterministic `ntt3_F` Bernstein Split + LP Ideal Labelling:** `<= 2^{59.5}` B (reserved `2^{60}` B `<< 16 EiB`).
4. **Stage-4 Native 32D-`K_f` Boot Root (`x_EA = 59.3229`, `LB = 58.0`):** Retained DF post-order root `|P| = 2^{56.8517}` B (`0.112` EiB), peak `8*|P| = 2^{59.8517}` B (`0.90` EiB `<< 16 EiB`; Stage-2 Kronecker-NTT3 / BLS MPE ECM `<= 601.97 MiB`).

## 3. Full Sensitivity Table (`part="sensitivity"` & `part="dag_survivors"`)

All scenarios below run live inside `audit_bn_cm.py` under strict single-ball `< 16 EiB`/node `MODE-I`:

| Scenario | Description | `exp_run` (`log2 u`) | Claim (`log2 u`) | Delta vs Primary |
|---|---|---:|---:|---:|
| **Primary (`MODE-I`, `T=1`, `0.45/2.70`)** | Single-ball `< 16 EiB`/node, full-domain live ECM 2LP, star-merge `mdepth=1.7011`, const-cap FVS, 256-lattice 32D boot | **`83.898822`** | **`82.540369`** | `0.0000` b |
| **`gate_0395_slack`** | Gate coin raised to `0.395` (`+0.005` success slack above `0.39`) | `83.898822` | `82.558748` | `+0.0184` b |
| **`zero_fvs_baseline_h3`** | Disable FVS completely (baseline `v43b-pi` `W_row = 100.5` H3 LA) | `83.948607` | `82.590154` | `+0.0498` b |
| **`free_1lp_only_zero_2lp`** | Drop `f`-side 2LP completely (1LP only) | `84.003666` | `82.645213` | `+0.1048` b |
| **`stress_plus15pct_lat_and_cof`** | Simultaneous `+15%` stress penalty on BOTH `lat_b2` and `cof` | `84.006458` | `82.648006` | `+0.1076` b |
| **`dag_survivors` (`2x`)** | Double all DAG survivor rates (`H5-W` stress) | `83.900766` | `82.542314` | `+0.0019` b |
| **`boot_charge = 4` / `boot_surv x4`** | Halve boot yield / quadruple boot survivors (`H-EB0/H-EB1` stress) | — | `82.543838` / `82.543606` | `+0.0035` / `+0.0032` b |
| **`plain_dickman_rho`** | Unpadded plain Dickman `rho` (reference only) | `83.718746` | `82.360294` | `-0.1801` b |

# STORAGE.md: Direct Live `uint8_t` Score Nodes (`N_sc = 2^{34}`), 1-Hop 16 KiB Link Pages, `256 MiB` Epoch Buffers, `0LP` FVS Schur-BW LA, DF Post-Order 12D Boot Root, & `< 2^{98}` B Storage Bound

Every active node is a flat 64-bit scalar ARMv8.1-A machine with `< 2^{64}` B (`16 EiB`) of RAM (`-march=armv8.1-a -mgeneral-regs-only`, `simd = 0`). Every cross-node byte uses explicit 1-hop 16 KiB link pages (`x16_bfi_copy_page + x16_bfi_copy`: `1.04984 <= 1.05` Ir/word) or `256 MiB` epoch buffers (`KERNELS.md` §1). Figures are certified in exact rationals by `stor_bound.py` (`R_A_MAX = 572.0 > 568.83`, byte-identical to `bn254-c6-v43b-pi`; its `2^{64}` LA envelope encloses Stage-3 LA, while `stor_bound.LA_LOG2` is an unused constant retained for byte-identity) and checked live by `audit_bn_cm.py` (`ram_ok = True`, `storage_cap_abort = 0.0`). All `<= 2^{...}` ceilings round UP.

## 1. Node Types and Pointwise Storage (`stor_bound.py`, `R_A_MAX = 572.0`)

| Node / Working Set | Count | Holds | RAM (`log2` B) |
|---|---:|---|---:|
| Emitter `E_0` | 1 | `a_tab` (`2^{59.908}` B), 4-plane `B_TAB` (`2^{63.355}` B), 1 TiB Direct Q12 `isqrt` LUT (`2^{40}` B), `rot_tab_ptr`/`mod3_tab`/`step1_tab`/`WHEEL_TABS` (`2^{31}` B), `2^{34}` x `256 MiB` epoch buffers (`4.0` EiB), `ptrs[]`, `E_a+` build arrays | `<= 63.9234` (`15.1725` EiB, `0.8275` EiB headroom) |
| Score Nodes `C[s]` | `2^{34}` | `C[0]` pure dummy sink (`0` real rows, `2^{49}` B sink); `C[1..2^{34}-1]` `uint8_t` (1 B/cell) score rows (`k_s in {41..48}`, `max_rows <= 224,263`) + incoming link pages + `init_a[]`/`a(ai)` + `32 MiB` Stage-B `g`-SRT (`L_g = 13`, `21,814,168` B) | `<= 62.9035` (`blocks <= 31,382,183 < 2^{25}`) |
| List Nodes `L[k]` | `2^{36}` | Band 1 lists (`[2^{24}, 2^{40.0})`) + dual `pos`/`neg` Lower-B2 prefix lists (`B2_PREFIX_MULT = 2.8367 <= 2.84`, pool `<= 2^{57.580}` B `<< 2^{60}` B) + build scratch (`<= 2^{56}` B) + outgoing link pages | `<= 60.1617` (`model 60.1611`) |
| `f`-Cofac Node `F_0` | 1 | Dedicated 1-owner `f`-side smoothness batch (`L_f = 50`, `2^{50}` 7-limb leaves: `P_f`, 51 `PT`/`RT` levels, `8x` NTT scratch, incoming link pages) | `<= 62.7916` (`6.9239` EiB `<< 16 EiB`) |
| **Stage-3 `0LP` FVS Schur-BW LA (`bw=(768,32)`)** | **1 owner** | **Unified SpMV (`2^{59.6419}` B, `D_rows*(512*nslots + 45*8) + 2^{56}*nslots + 65*2^{53}`) + PM-Basis (`2^{62.1872}` B) at `nu_fvs = 0.50, nu_d = 2.00, mdepth = 1.0`** | **`62.4153`** (`5.33` EiB `< 16 EiB`; exact-rat UB at `bw=(768,32)`: **`62.4901`** `<= 62.4902`) |
| Stage-2 Deterministic `ntt3_F` Split | 1 owner | Bernstein prime-tree split on `(61.50, 61.50)` (`Q = 2^{46}`) | `<= 60.0` (`<= 2^{59.5}` B actual `<< 16 EiB`) |
| Stage-4 12D Boot Root (`x_EA = 62.9375`) | 1 owner | DF post-order all-primes root `|P| <= 2^{60.4664}` B (`1.382` EiB), Karatsuba/NTT3 peak `4.2714*|P| <= 2^{62.5611}` B `<< 8*|P| <= 2^{63.4664}` B `< 16 EiB` (ECM `<= 601.97 MiB`) | `<= 63.4664` (`peak 62.5611` `<< 16 EiB`) |
| **Pointwise Total Machine Storage** | — | **Exact-rational `stor_bound.py` envelope (`storage_cap_abort = 0.0`)** | **`<= 97.5798`** (`< 100`) |

## 2. Data Movement & Single-Worker Working Sets (`< 2^{64}` B)

All Stage-1/Stage-2 sieve, link-page, `256 MiB` epoch-buffer, Stage-B `g`-SRT (`21,814,168` B `<= 32 MiB`), `F_0` (`<= 2^{62.7916}` B), and Stage-4 12D DF post-order boot root working sets and movement charges are unchanged from `bn254-c6-v43b-pi`.
- **Stage-3 `0LP` FVS Schur-BW LA (`bw = (768, 32)`):** Uses the same `bw = (768, 32)` blocking as `bn254-c6-v43b-pi` (`KERNELS.md` §3's 16-vector BW kernel slope is unchanged and executes in 16-vector passes), with strictly smaller 1-owner RAM (`2^{62.4153}` B = `5.33` EiB vs `2^{63.2867}` B = `9.95` EiB); exact-rational certificate in `audit_bn_cm.py::verify_fvs_ladder_and_caps` at `U_ub = 2^{54.62}` and `bw = (768, 32)` proves **`<= 2^{62.4901}` B (`<= 2^{62.4902}` B) `< 16 EiB`** (`la_1own_rat_log2 = 62.4901 >= 62.4153`), enclosed inside `stor_bound.py`'s `2^{64}` LA envelope. O(nnz) bucket-queue greedy unit-pivot FVS selection (`<= 32 * D * 45 <= 2^{54.5}` Ir, charged `2^{60}` Ir) and the final triangular lift `x_1 = -A_{11}^{-1} A_{12} v` (`<= 45 * 8.70 * D <= 2^{52.5}` Ir, charged `2^{58}` Ir) are priced directly inside `misc` in `model_bn254_v6.py`.

## 3. Full Sensitivity Table (`part="sensitivity"` & `part="dag_survivors"`)

| Scenario | Description | `exp_run` (`log2 u`) | Claim (`log2 u`) | vs Primary |
|---|---|---:|---:|---:|
| **Primary (`0LP`, `0.50/2.00` FVS)** | Unmodified `v43b-pi` sieve & 12D boot + exploratory `0LP` FVS | **`83.964263`** | **`82.605810`** | `0.0000` b |
| **`zero_fvs_baseline_h3`** | Disable FVS completely (reproduces `bn254-c6-v43b-pi` `25323c9f` verbatim) | `84.030218` | `82.671767` | `+0.0660` b |
| **`fvs_u19_040`** | `0LP` FVS near `2^{19}` max (`nu_fvs = 0.40, nu_d = 1.35`; `bw = (768, 40)`) | `83.921180` | `82.562727` | `-0.0431` b |
| **`fvs_stress_060_250`** | Cap-violation FVS stress (`nu_fvs = 0.60, nu_d = 2.50`; `bw = (768, 32)`, RAM `2^{62.64}` B) | `84.009841` | `82.651389` | `+0.0456` b |
| **`live_wmid_stress_058_245_925`** | Small-`U` `w_mid`-distorted live ratios (`0.58 / 2.45 / k5=9.25`; `bw = (768, 32)`) | `84.006795` | `82.648342` | `+0.0425` b |
| **`gate_0395_slack`** | Gate coin raised to `0.395` (`+0.005` success slack above `0.39`) | `83.964263` | `82.624189` | `+0.0184` b |
| **`dag_survivors` (`2x`)** | Double all DAG survivor rates (`H5-W` stress) | `83.974180` | `82.615728` | `+0.0099` b |
| **`boot_charge = 4` / `boot_surv x4`** | Halve boot yield / quadruple boot survivors (`H-EB0/H-EB1` stress) | — | `82.630111` / `82.664854` | `+0.0243` / `+0.0590` b |
| **`conservative_srt_030_cap`** | Conservative `srt_g_cap = 0.30` | `83.976238` | `82.617785` | `+0.0120` b |
| **`conservative_lat_b2_plus10`** | `+10` Ir/leaf stress on `C_LAT_LEAF` | `84.479101` | `83.120648` | `+0.5148` b |
| **`de_bruijn_tighter_cap`** | `cap_scale = 0.80` | `83.930425` | `82.571972` | `-0.0338` b |
| **`flat_shared_arena_ref`** | Shared-arena reference | `83.914485` | `82.556032` | `-0.0498` b |
| **`plain_dickman_rho`** | Unpadded plain Dickman `rho` (reference only) | `83.792691` | `82.434239` | `-0.1716` b |

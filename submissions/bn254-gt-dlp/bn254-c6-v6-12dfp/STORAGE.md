# `STORAGE.md` — Memory & Single-Owner RAM Accounting (`BN254 v6`)

Every working-set buffer in `BN254 v6` fits comfortably inside a single 64-bit address space (`< 2^64 B` per owner), and total global read-only/bucket storage stays strictly below `2^100 B` (matching accepted `bn15r` `c5-gt`).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At the `v6` operating point (`Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 96.9622`, `log2_Na = 43.8650`, `log2_Vb = 54.2629`, `L_max = 1010.65`, `P_cut = 44.0`, `P_mid = 53.0`, `lx2 = 59.0`, `W_SB = 128`):

1. **Active Score Arena (`uint16_t`, row-major `(t, a_slot)`):**
   - `(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table)` = **`2^54.8504 B`** (`~32.8 PiB`, **`570x` below `2^64 B`**).
   - For reference, even `V_b * 2 B` = **`2^55.2629 B`** (`<< 2^64 B`).
2. **Band-2 Single-Prime `a`-Table + 1-Byte `sb_u7` Directory (`p <= 2^53.0`, `W_SB = 128`):**
   - `a`-entries: `16 B * N_a` = `2^47.8650 B` (`~259 TiB`).
   - Super-bucket directory (`8 B` offset + `1 B` `sb_u7`): `9 B * (2^53 / 128)` = `2^49.1699 B` (`~640 TiB`).
   - Combined score arena + 1 prime Band-2 table: **`2^54.8894 B`** (`<< 2^64 B`).
   - Over `10,000` primes fit simultaneously inside a single `< 2^64 B` owner alongside the score arena.
3. **1-Owner `Z[omega]` Block Wiedemann LA (at 4-parameter hard cap `N_sw_cap = 2^44.3538`, `bw = (768, 35)`, `nslots = 7`):**
   - Two `5`-lane `F_r` vector buffers (`2 * 160 B * 7 * N_sw_cap`) + matrix rows (`W_ROW * 8 B * N_sw_cap`) + dense-band TABM tables (`2^55.04 B * 7`) + `65 * 2^53 B` scratch = **`2^59.68 B`** (**`20x` below `2^64 B`**, zero multi-owner LA split needed).

## 2. Global Storage & Strict `< 2^64 B` Per-Owner Sensitivity

- **Global read-only `a`-tables + Band-3 buckets:**
  - Band-1 `[2^14, 2^44)` `a`-tables + `256`-shell directories: `16 * N_a * n_lo + 8 * 256 * sp_lo = 2^87.9 B`.
  - Band-2 `[2^44, 2^53)` `a`-tables + `sb_u7` directories: `16 * N_a * n_hi + 9 * (sp_hi / 128) = 2^94.5 B`.
  - Band-3 `[2^53, 2^59)` 8-byte `log_p`-grouped buckets: `8 B * V * ln(59/53) = 2^96.8 B`.
  - **Total global RAM:** **`2^97.40 B < 2^100 B`** (`audit_bn_cm.py` asserts `ram_log2 < 100.0`).
- **Strict `< 2^64 B` Per-Owner Buffered Sensitivity (`part="sensitivity"`):**
  If Band-1 and Band-2 `a`-tables are sharded across independent `< 2^64 B` prime-owners communicating with the score-arena owner via explicit 8-byte `log_p`-grouped cross-owner buffers (`C_HIT_LO = 6.25`, `C_HIT_HI = 14.5`, plus `v3`'s conservative `C_BKT_B3 = 25.0`), the gated claim rises by only `+0.137` bits to **`2^86.6856`** — still **`-0.903` full bits below `#1` (`87.5884`)**.

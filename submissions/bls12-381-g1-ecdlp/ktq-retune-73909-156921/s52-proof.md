# s52: TMR full-child kernel, zero-spill batch round, pinned-stride leaf and aligned beta chunk

## Overview and Invariance Contract

This package inherits the five pure scalar AArch64 arithmetic listings (`fbody`, `bbody`, `kernel_text`, `leaf_text`, `beta_keys`) from submission b574037 (`tmr-full-s52`, `114.555831293763`, which itself starts from `e983f8e3-18b7-45d0-a618-3cdd1ddf3131`). The s52 revision modified only those listings; **this revision changes only the design point (`K_T`, `Q`) and the `N`-dependent quantities that follow**, preserving every algorithmic object, curve point, orbit key, hash table, pool, owner initialization (`DC ZVA`, s50), round schedule, stopping bound, gate construction, `F = 998/1000`, `ELL = 193/200`, `C = 2000`, memory layout, and conditional hypotheses `H1'`, `H2–H5` identically. No prior verdict or score is used as a premise; every listing is executed word-for-word by `a64.py` against exact big-integer curve arithmetic inside this package.

| Component | `e983f8e3` (Base) | This Package (`s52`) | Saving | Audit |
|---|---:|---:|---:|---|
| Forward batch element (`fbody`) | `292` | **`288`** | `-4` instr/pair | `batch_audit.run`, `range_audit` |
| Backward batch element (`bbody`) | `875` | **`863`** | `-12` instr/pair | `batch_audit.run`, `range_audit` |
| Batch round per lane (`BATCH`, 32-unroll) | `9337/8` (`1167.125`) | **`9209/8` (`1151.125`)** | `-16` instr/pair | `batch_audit.run` |
| Batch round fixed (`ROUND_FIXED`) | `172597` | **`172594`** | `-3` instr/round | `batch_audit.run` |
| Leaf kernel (`KX`, both paths) | `214` (213/214) | **`211` (211/211)** | `-3` instr/leaf | `kernel_check.run_leaf` |
| Full-child kernel (`KERNEL`, branchless) | `894` (12 paths) | **`778` (1 path)** | `-116` instr/full | `kernel_check.run` |
| Giant beta key kernel (`BETA`) | `139` | **`138`** | `-1` instr/giant | `batch_audit.beta_run` |
| Baby point price (`BABY_PT`) | `874.0835` | **`858.8944`** | `-15.189` instr/pt | `evidence.ledger` |
| Giant point price (`GIANT_PT`) | `1083.2085` | **`1067.0194`** | `-16.189` instr/pt | `evidence.ledger` |
| **Expected work `log2 W(2^20)`** | **`114.579559759048`** | **`114.555827210730`** | **`-0.023732548318` bits** | `log_bound.certify` |
| **Every-run ceiling `log2((W_run+64)/2000)`** | `117.550859154293` | **`117.530378341483`** | `-0.020480812810` bits | `log_bound.certify` |

## 1. Forward Element `fbody` (`292 -> 288`, `-4` instructions/pair)

1. **Zero carry word `c7` at `j = 4` (`-2` instr):** In `fbody`, the scalar `c_{i-1}` lies in `[0, p + Delta)` with `Delta = 2^376`, so its top limb satisfies `a_5 <= floor((p + Delta - 1) / 2^320) < 0.106 * 2^64`. Every `D_j = X_j + V_j` satisfies `0 <= D_j <= 2p`. Because `fbody` iterates `A[5]` down to `A[0]` against `D_0..D_5`, the 5-row partial sum after `j = 4` satisfies
   `S_4 = sum_{j=0}^4 a_{5-j} D_j <= (a_{5,max} + 4*(2^64 - 1)) * 2p < 0.834 * 2^448 < 2^448`.
   Hence word 7 (`c7`) is identically zero through `j = 4` and is allocated only at `j = 5` (`first7 = True` at `j = 5`). `range_audit.py` asserts this exact integer inequality (`rp7_j4_below_2pow448`).
2. **Live argument cursor in `x30` (`-2` instr):** In `e983f8e3`, `x30` was spilled and reloaded on every element solely because `x30` was reused as a temporary in the `e = v2 - x1` loop and as `c7`. Using dead `C[6]` in the `e` loop (written before `row0` initializes `C`) and dead `A[1]` as `c7` at `j = 5` (after `A[1]` is consumed at `j = 4`, unused in `mont_word` and `corr_l`) keeps `x30` live across `.Lfwd`, replacing `ldr x30,[sp,#32]; ldp x1,x2,[x30,#32]!; str x30,[sp,#32]` with a single `ldp x1,x2,[x30,#32]!`.

## 2. Backward Element `bbody` (`875 -> 863`, `-12` instructions/pair)

1. **1-temp interleaved `6x6` schoolbook (`mul6x6_1t`, `-3` instr):** Interleaving `mul T, A[j], b(k); adds/adcs W[k+j], W[k+j], T` uses a single temporary `T = x28` instead of a 6-register `TT[0..5]` array at identical instruction count (`146` instructions), freeing 5 registers so `PM` (`x5..x10` = `p`) remains live across `mul_tmr`. This deletes `ldrow(PM, sp, PO)` (`-3` instr/element) and `strow(PM, sp, PO)` in the prologue (`-3` instr in `ROUND_FIXED`: `172597 -> 172594`).
2. **Direct-target `tmr_l` into `U` (`-6` instr) + single `mov #384` (`-1` instr):** After `mul6x6_1t` finishes, `t` (`U`) is spilled once to `[sp, #TSP]`, freeing `U = x19..x24`. Chunk `i = 0` of `tmr_l` loads `T_0` directly into `U` and adds `first_into = W[6..11]` into `U`, then chunks `i = 1..6` reuse `W[6..11]` as `TT`. Keeping stride `#384` in dedicated register `x29` (`load_stride = True`, one `mov x29, #384`) deletes the second `mov #384` (`-1` instr), and writing `inv_i` straight into `U` deletes `strow(inv_i)` + `ldrow(inv_i)` (`-6` instr).
3. **Reversed `RP7(t, D)` (`rev = True`, `-2` instr):** Since `t < p + Delta` and `D_j <= 2p`, iterating `t[5]..t[0]` against `D_0..D_5` satisfies the exact same `S_4 < 2^448` bound as `fbody`, enabling `c7` only at `j = 5` (`-2` instr). `range_audit.py` verifies `rows.rowprod(..., rev=True)` across 40 boundary and random cases against `rows.corr(rows.rp_raw(a, Dr))`.

## 3. Pinned-Stride Leaf Kernel (`214 -> 211`, `-3` instructions/leaf)

In `measure/side_lanes.S`, the 94-instruction `SQR` block never writes `x1`, and `leaf_text` never reads `P`. Replacing the dead `ldr x1, [x30, #8]` in the prologue with `mov x1, #384` (keeping `FIXED = 32` unchanged) and advancing the cursor with `ldr x0, [x30, #32]!` keeps `x1 = 384` pinned across all lanes, eliminating both in-lane `mov #384` instructions (`-2` instr). Flipping `b.lo .Lx_done; b .Lx_y` to `b.hs .Lx_y` (falling through directly to `.Lx_done`) removes the redundant unconditional jump (`-1` instr), making both leaf paths strictly `211` instructions (`kernel_check.run_leaf`: exact canonical `x'` and `[out, out+48)` footprint across all 365 cases for both `-` and `+`).

## 4. TMR Full-Child Kernel (`894 -> 778`, `-116` instructions/full child)

Previously, only the leaf kernel used the already-allocated and already-charged 7-chunk `TMR` table (s47), while the full-child kernel (`kernel_text`) retained word-serial Montgomery `BY1` and `MULNR` (`894` instructions worst-case over 12 paths). `s52` applies the identical 7-chunk `TMR` table to both `x'` and `y'` in `kernel_text` with `x1 = 384` pinned in the prologue (replacing the three unused `p` frame stores in the prologue by loading `TMRB` at `.Lsk_consts + 56` and storing `(n0, TMRB)` to `[sp, #0]`, keeping `FIXED = 32` exact):

1. **`x'` via `SQR` + `tmr_leaf_block(branch=False)` + branchless `subs`/`csel`:** Produces exact canonical `x'` in `[0, p)` (reusing the s49 `0 < S_0 + e < 9p + h < R` proof).
2. **`y'` via `w = x' + v2`, `mul6x6_1t(w, mu)` + `tmr_l`, `S - n2`, `tmr_corr_l`, branchless `subs`/`csel`:** Because `0 <= mu < p + Delta` and `0 <= w < 2p`, the high half satisfies `H = floor(w * mu / R) < h_w = (p + Delta)(2p)/R < 0.212p`. With `T_0 >= p` (the s47 `+p` bias) and `T_i < p` (`i = 1..6`), the raw TMR sum satisfies `p <= S < 8p + h_w < R`. For any `1 <= n2 <= p - 1` (`y2 != 0` on `y^2 = x^3 + 4`), `0 < S - n2 < 8p + h_w < R` holds with zero borrow and zero carry word. `range_audit.py` asserts `8p + h_w < R` in exact rationals and checks 200 random/boundary `(w, mu, n2)` cases. Then `tmr_corr_l` maps `S - n2` into `[0, p + Delta)` and one branchless `subs`/`csel` yields exact canonical `y'`.
3. **Rows:** Loads `n0` into `x18`, runs the unchanged `rowbuild(y')` and `rowbuild(x')`, and restores `TMRB` into `x18`.
4. **Result:** Single branchless path of **`778`** instructions (`worst == best == 778`, zero signed comparisons, `FIXED = 32`), verified by `kernel_check.run` (509 cases each for `-` and `+`, 672-byte stride canary audit, 640-byte negative control) and end-to-end inside `batch_audit.run` (both random curve points and `table_level`).

## 5. Aligned Chunk 3 in `beta_keys` (`139 -> 138`, `-1` instruction/giant)

In `batch_audit.py`, the eight 48-bit chunks of canonical `X < p < 2^381` start at shifts `0, 48, 96, 144, 192, 240, 288, 336`. Chunk `i = 3` (`shift = 144, width = 48`) sits at bits `16..63` of limb `x13` (`divmod(144, 64) = (2, 16)`, `16 + 48 = 64`), never crossing a 64-bit word boundary. Extracting it with a single `lsr x3, x13, #16` instead of `extr + and` is bit-identical and reduces `BETA_LEN` from `139` to `138` (`batch_audit.beta_run` and `beta_composed` pass across all curve and boundary points).

## 6. Regenerated Exact Ledger Figures

At the retuned design point `K_T = 73909/10000`, `Q = 156921/125000`, `N_b = 278081993926198672399008286084263217`, `N_g = 221514324027853722891620852279381995`, `lo_run = 0.99888220577219534119...`, `K = 7202280856715389586`:
- `BABY_PT = 11257700252907/13107200000` (`858.8943674...`), `GIANT_PT = 13985636252907/13107200000` (`1067.0193674...`)
- `W_run = 479953627230568756514720629955568406370` (`every-run <= 2^117.530378341483` units)
- `E_run = 156393908636956402415832412106483204887`
- `E_num = 1126392854282846724778895998043138370365995771849475410206`
- **`W(2^20) <= 2^114.555827210730`** (`exact score = 114555827210730/10^12`; exact log2 `114.5558272107296685...`, certified by `log_bound.py` with both `±1e-12` grid shifts rejected).
- Peak memory (`< 2^125.0002` bytes) and challenge-independent preprocessing (`< 2^105.3247` units) are unchanged.

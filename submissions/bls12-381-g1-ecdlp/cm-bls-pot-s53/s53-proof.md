# s53: 1-temp SQR with live out/buf, branch-canonical full-child and beta kernels, and post-indexed batch-round updates

## Overview and Invariance Contract

This revision starts from accepted submission `b5740375-cdd6-478b-8210-f37d84414291` (`s52`, `114.555831293763`) and modifies **only five pure scalar AArch64 arithmetic listings** (`fbody`, `bbody`, `leaf_text`, `kernel_text`, `beta_keys`), preserving every algorithmic object, curve point, orbit key, hash table, pool, owner initialization (`DC ZVA`, s50), round schedule, stopping bound, design point (`K_T = 14717/2000`, `Q = 5/4`), `N_b`, `N_g`, `K`, `lo_run`, gate coin `K = 7202427322808921017`, `F = 998/1000`, `ELL = 193/200`, `C = 2000`, memory layout, preprocessing cost, and conditional hypotheses `H1'`, `H2–H5` identically. No prior verdict or score is used as a premise; every modified listing is executed word-for-word by `a64.py` against exact big-integer curve arithmetic inside this package, with every conditional branch exercised both ways.

| Component | `b5740375` (`s52`) | This Package (`s53`) | Saving | Audit |
|---|---:|---:|---:|---|
| Forward batch element (`fbody`) | `288` | **`287`** | `-1` instr/pair | `batch_audit.run` |
| Backward batch element (`bbody`) | `863` | **`861`** | `-2` instr/pair | `batch_audit.run` |
| Batch round per lane (`BATCH`, 32-unroll) | `9209/8` (`1151.125`) | **`9185/8` (`1148.125`)** | `-3` instr/pair | `batch_audit.run` |
| Batch round fixed (`ROUND_FIXED`) | `172594` | **`172594`** | `0` | `batch_audit.run` |
| Leaf kernel (`KX`, both paths) | `211` (211/211) | **`209` (209/209)** | `-2` instr/leaf | `kernel_check.run_leaf` |
| Full-child kernel (`KERNEL`, 4 paths) | `778` (1 path) | **`765` (764/765)** | `-13` instr/full | `kernel_check.run` |
| Giant beta key kernel (`BETA`, 4 paths) | `138` (1 path) | **`130` (130/130)** | `-8` instr/giant | `batch_audit.beta_run` |
| Baby point price (`BABY_PT`) | `858.8944` | **`854.9851`** | `-3.909` instr/pt | `evidence.ledger` |
| Giant point price (`GIANT_PT`) | `1067.0194` | **`1055.1101`** | `-11.909` instr/pt | `evidence.ledger` |
| **Expected work `log2 W(2^20)`** | `114.555831293763` | **`114.544507818417`** | **`-0.011323` bits** | `log_bound.certify` |
| **Every-run ceiling `log2((W_run+64)/2000)`** | `117.527179082161` | **`117.515879139431`** | `-0.011300` bits | `log_bound.certify` |

## 1. 1-Temp Interleaved `SQR` and Live `out`/`buf` in Leaf (`211 -> 209`, `-2` instr/leaf)

In `s52`, the extracted 91-instruction `SQR` block held all 10 diagonal half-products in separate temporary registers (`x0, x2, x3, x4, x29`) simultaneously, clobbering `x2` and `x29` and forcing both `leaf_text` and `kernel_text` to reload `buf = [x30, #24]` for `e` and reload `out = [x30]` before storing `x'`.

`s53` replaces `sqr_block()` with `sqr_1t()` (`kernel_check.py`), an exact 91-instruction 1-temp interleaved 384-bit squaring routine:
- All 15 off-diagonal products `A[j]*A[i]` (`i < j`) are accumulated into `W = x19..x28` using only `x0` as temporary (`A = x11..x16` remains unmodified).
- `W` is doubled (`adds`/`adcs` × 10).
- The 6 diagonal squares `A[i]^2` are computed and added on the fly using only `x0` (`mul`/`umulh` do not modify NZCV flags), writing word 0 to `x17`, words 1..10 to `x19..x28`, and word 11 to `x16`.
- Total length is identically **91 instructions**, and `x1` (`384`), `x2` (`out`), `x3`, `x4`, `x5..x10` (`p`), `x18` (`TMRB`), `x29` (`buf`), and `x30` (`args`) are never written.

In `tmr_leaf_block`, using `x0` as the `madd` address temporary and loading chunk 0 directly into `x11..x15, x0` likewise leaves `x2` (`out`) and `x29` (`buf`) untouched across all 7 TMR chunks and the `m(idx)*p` correction:
1. **Live `buf` in `x29` (`-1` instr/leaf):** `.Llane` loads `buf` once via `ldr x29, [x30, #24]`, loads `mu` from `[x29, #omu]`, and later loads `e` directly from `[x29, #OE]` (deleting `ldr x0, [x30, #24]`).
2. **Live `out` in `x2` (`-1` instr/leaf):** The prologue already loads `args[0][0]` into `x2` (`ldr x2, [x0]`), and each lane tail advances the cursor via `ldr x2, [x30, #32]!; cbnz x2, .Llane`, maintaining the loop invariant `x2 = args[i][0]` (`out`) at `.Llane` and deleting `ldr x2, [x30]` from both `.Lx_done` and `.Lx_y`.
3. **Result:** Both leaf paths execute identically **`209`** instructions (`worst == best == 209`, `FIXED = 32` unchanged), verified by `kernel_check.run_leaf` across 365 cases each for `-` and `+`.

## 2. Full-Child Kernel `kernel_text` (`778 -> 765`, `-13` instr/full child)

Applying `sqr_1t()` and live `x2 = out`, `x29 = buf` (loading `TMRB` via `x3` in the prologue so `x2 = args[0][0]` survives into `.Llane`, keeping `FIXED = 32` unchanged) yields four independent savings in `kernel_text`:
1. **Live `buf` in `x29` (`-1` instr):** Deletes `ldr x0, [x30, #24]` before adding `e`.
2. **Live `out` in `x2` (`-1` instr):** Deletes `ldr x2, [x30]` before storing canonical `x'`.
3. **Direct `[x29, #omu]` `mu` read (`-2` net instr):** Deletes the 3-instruction stack spill `stp x11..x16, [sp, #96..#128]`. Because `x29 = buf` remains live across `SQR`, `TMR`, and `w = x' + v2` (which loads `T` into `x0` and reads `v2` from `[x0, #96]` with `ldp x0, x2, [x0, #96]` last), `mul6x6_1t` reads `mu` directly from `[x29, #omu]`, followed by one `ldr x29, [x30, #16]` before reading `n2`.
4. **Branch canonicalization for `x'` and `y'` (`-4` and `-5` instr):**
   - For `x'`: replaces `6 * csel` with `b.hs .Lx_ge` around the 13-instruction store-`x'` + load-`v2` + `w = x' + v2` block (`b .Lx_join` on fallthrough, falling through from `.Lx_ge`), cutting worst-case cost from 6 to 2 (`-4` instr).
   - For `y'`: replaces `6 * csel` with `b.hs .Ly_ge` around two register-renamed `rowbuild(y')` + `rowbuild(x')` tails (using `Y6 = x19..x24` vs `X6 = x11..x16`, both exiting with `ldr x2, [x30, #32]!; cbnz x2, .Llane`, with `b .Lexit` once per call on the fallthrough tail, already covered by `CHARGED = FIXED + 1 + 5 + 2 = 40`), cutting 6 `csel`s to 1 `b.hs` (`-5` instr).
5. **Result:** 4 paths of lengths `764, 764, 765, 765` (`worst = 765`, `FIXED = 32`), verified by `kernel_check.run` (`509` cases each for `-` and `+`, both conditional branches exercised both ways, 672-byte stride canary audit, 640-byte negative control).

## 3. Four-Path Branch-Canonical `beta_keys` (`138 -> 130`, `-8` instr/giant)

In `batch_audit.py`, `BETA_BODY` previously executed `6 * csel` after `RED` (`S - q*p` vs `S - (q+1)*p`) and `4 * csel` after `X + Y - p` for `beta^2 * X = p - ((X + Y) mod p)`. Replacing both `csel` sequences with unsigned `b.hs` branches (`b.hs .Lb_ge` selecting `Y = AC` vs `Y = TM`, and `b.hs .Lb_ge2` / `b.hs .Lb_ge3` selecting `W` vs `Y` for `p[0..3] - ...`) gives four register-renamed tails of **identical per-point length `130` instructions** (`-5 - 3 = -8` instr/giant):
- Every non-final slot takes `cbnz x4, .Lbeta` at exactly `130` instructions on all four paths.
- Deleting the unused `ldr x29, [x2, #64]` from the prologue reduces prologue+epilogue from `26` to `25`, so even with the final-slot `b .Lbx` on paths 0–2, total call cost is `130 * n + 25` or `130 * n + 26` (never exceeding the prior fixed overhead `26`).
- `batch_audit.beta_run` includes explicit EGCD witnesses `pow(be, 2, P)` and `P - s1` exercising all 4 paths and asserting all 3 `b.hs` conditional branches are taken both ways, in addition to 663 standalone `RED` cases against exact `S % p`.

## 4. Batch Round `fbody` (`288 -> 287`, `-1`) and `bbody` (`863 -> 861`, `-2`)

Because batch slot stride `SB = 496` bytes is a multiple of 8 in `[-512, +504]` (valid AArch64 `ldp`/`stp` post-index immediate range):
1. **`fbody` post-indexed `x3` advance (`-1` instr/pair):** Biases `x3` by `+RW.OC` (`+288`) at prologue entry (`add x3, x1, #288` replacing `mov x3, x1`). All slot offsets shift by `-288` (all within `[-288, +192]`), and the final store of `c_i` uses `stp C[5],C[6],[x3,#32]; stp C[3],C[4],[x3,#16]; stp C[1],C[2],[x3],#496`, deleting `add x3, x3, #496`.
2. **`bbody` pinned `x30 = #384` (`-1` instr/pair):** Replaces the unused `ldr x30, [sp]` before `.Lbwd` with `mov x30, #384` (keeping `ROUND_FIXED = 172594` exact!), passes `stride_reg='x30', load_stride=False` to `tmr_l` (deleting `mov x29, #384`), and checks the `.Lbwd` end pointer using `x29` (`ldr x29, [sp]; cmp x2, x29; b.ne .Lbwd`).
3. **`bbody` post-indexed `x2` decrement (`-1` instr/pair):** Biases `x2` by `+240` (`sub x2, x2, #544` from `x3 = BUF + n*SB + 288`, end sentinel `BUF - 256`). All offsets shift by `-240` (`OC - SB - 240 = -448` through `OPT + 8 - 240 = +248`, all in `[-512, +504]`), and the final `j = 5` `D_5` load in `rowprod(..., rev=True, post=-SB)` executes `ldp W[4],W[5],[x2,#32]; ldp W[2],W[3],[x2,#16]; ldp W[0],W[1],[x2],#-496`, deleting `sub x2, x2, #496`.
4. **Result:** `BATCH = 9185/8 = 1148.125` (`-3` instr/pair), `ROUND_FIXED = 172594` (unchanged), verified end-to-end by `batch_audit.run({})` and `batch_audit.run({'table': True})` with exact register/memory preservation and `a64.py` visit counts.

## 5. Regenerated Exact Ledger Figures

At the unchanged design point `K_T = 14717/2000`, `Q = 5/4`, `N_b = 276877925431160214981668449342404407`, `N_g = 221502340344928171985334759473923525`, `lo_run = 0.99886189284600805486...`, `K = 7202427322808921017`:
- `BABY_PT = 11206461193451/13107200000` (`854.9851374...`), `GIANT_PT = 13829539593451/13107200000` (`1055.1101374...`)
- `W_run = 475154212042293596195741882931841881961` (`every-run <= 2^117.515879139431` units)
- `E_run = 155168486974433895693022642791753391971`
- `E_num = 1117589750223582856112034727111242783809830325425792257931`
- **`W(2^20) <= 2^114.544507818417`** (`exact score = 114544507818417/10^12`; exact log2 `114.5445078184162557...`, certified by `log_bound.py` with both `±1e-12` grid shifts rejected).
- Peak memory (`< 2^125.0002` bytes) and challenge-independent preprocessing (`< 2^105.3247` units) are unchanged.

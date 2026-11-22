# Kernels (Raw-Opcode-Verified Baseline Scalar ARMv8.0-A: `Arch2-Prod` 8B Zero-Copy, 1-Owner Build, `uint16_t`, `(m-k)&3` Bank Rotation, 1-Step Barrett `row_red5_z12`, Input-Uniform Bareiss LLL)

Every kernel cost is derived from explicit baseline scalar ARMv8.0-A assembly (`-march=armv8-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `bl`). `audit_cm.py` (`kernels_symbolic`) parses every fenced `asm` block directly from this file, validates opcodes & `LDP`/`STP` addressing, counts every instruction (`cmp`+`b.*`=2 ops, loop control included), simulates the complete 1-owner flat-array 2-pass build (incl. overflow-proof 64-bit modular step at `p = 2^64 - 59`) + 1-rec-at-a-time Band-2 multi-SB replay (incl. empty/odd/even buckets) + `(m-k)&3` + SWAR + `MONT5` + Barrett `row_red5_z12` + Shoup NTT + 4-prime Garner CRT + iterative Q32 `dsieve16` + 10-op/level `dsieve16_rank_step` + certified upper-bound `BOUND_CELL_HI` live, and asserts `opcodes/items <= charged`.

## 1. Collection Kernels (`Arch2-Prod-ZeroCopy-U32`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 416`, `W_SB = 32`, `P_m4 = 2^56.5`, `P_lo = P_f = 2^59`, `P_g = 2^64`)

Because the sieving region is a pure product of two 8D balls (`qa(a) <= s_c R2_a` with 2-adic radii `s_c = (1.0, 1.5, 2.0)` and `qb(b) <= R2_b`, `A^2 = 321,600`, `r_ab = 0.9815`), every canonical `a`-table entry (`N_a^can = V_a = 2^56.8437 < 0.90 * 2^57`, `N_a^tbl = 4 V_a = 2^58.8437`) is unconditionally in-region for every `b` (`V_b = 2^63.5877`, `L_max = 404.76 < W_tile = 416`, `L_b = 173.85`, `N_7 = 2^54.1460`, `V_kept = 2^120.4314`). There is **ZERO `q_key`, ZERO `q_tab`, ZERO `sb_min_q8`, ZERO binary search, ZERO global sorting invariant, and ZERO `f`-side Band 2** (`P_f = P_lo = 2^59`).
- **Bias-initialized `S_g` & zero-wrap `uint16_t` accumulators (`S_g, S_f`) with exact `2^57` bank stride:** Quarter-bit weights on `<= 705`b norms reach `s <= 3,000`. Initializing `S_g` lanes to `bias16 = (32768 - T) & 0x7FFF` (`1 <= T <= 3,000`, e.g. `T = 4 * 345 = 1,380` quarter-bits) keeps `(32768 - T) + s in [29768, 35767] subset [0, 65535]` with **ZERO 16-bit wraparound/carry between lanes** and `s >= T <==> bit 15 (0x8000) == 1`. Each column owner allocates `4 * 2^57` `uint16_t`s for `S_g` (`2^60` B) + `S_f` (`2^60` B) = `2^61.0` B `<< 2^64` B (exact power-of-2 bank stride `2^57` -> ZERO stride mismatch for `((m - k) & 3) << 57 | oid`).
- **Zero-copy 8B records across Band-1a/1b/2 with exact `(m - k) & 3` rotation:**
  - Band-1a (`[2^24, 2^56.5)`): `m`-partitioned `oid < 2^57` -> `rev_b1a_apply_u8` (`30/8 = 3.75 <= 3.80`).
  - Band-1b (`[2^56.5, 2^59)`): 8B `(m << 57) | oid` (`sub x4,x0,x_k57; and x4,x4,x_mask59`) -> flag-exact 8-way unrolled `rev_b1b_apply_u8` (`46/8 = 5.75 <= 5.80`) + `rev_b1a_apply_wrap8` (`4` ops) + `rev_b1b_apply_tail1` (`8` ops/rem hit, `<= +20` ops/bkt charged inside `b1_lk = 34.0`).
  - Band-2 (`g`-side `[2^59, 2^64)` only): 8B `(r << 59) | (m << 57) | oid` (`r = ka & 31`, `5 + 2 + 57 = 64` bits exact!). Outer loop state lives exclusively in callee-saved `x19..x22` (`x19=off_ptr`, `x20=ent_base`, `x21=rem_sb`, `x22=bptr_base`, advancing `add x22, x22, #256` per SB), read-only flush limit in `x18`, bank-rotate constants in `x3, x4` — physically disjoint from emit scratch `x9..x11, x15` (asserted automatically in `verify_kernels_symbolic`). 1-record-at-a-time `rev_b2_prod_emit_1` (`12 <= 12.5` Ir/hit, safe on every bucket length `0, 1, 2, ...`, ZERO odd-length overread) + `rev_b1a_apply_u8` (`3.75 <= 3.80`) = **`15.75 <= C_HIT_HI = 16.35`**.

### 1.1 `rev_sb_probe_fast` (Band-2 row start: `17 <= C_TILE_HDR = 18.0` Ir/row; multi-SB safe continuation to `.Lempty` = `.Lnext_sb`)

Computes the chord's starting super-bucket `x4 = x3 >> 5`, remaining chord width `x5 = min(x_wmax, x_end - x3)`, and first-SB offset pointer `x7 = x_off + (x4 << 3)` (`x19` initialized to `x7 + 8`, `x21` initialized to the chord's total intersected super-bucket count `n_sb >= 1`). Probing `sb_0` via `ldp x8, x9, [x7]; cmp x8, x9; b.hs .Lempty` NEVER aborts a multi-SB chord: `.Lempty` is the `.Lnext_sb` continuation entry of `rev_b2_prod_sb_hdr` (§1.2 below), which advances `x22 += 256`, decrements `x21 -= 1`, and branches back to `.Lsb` to scan ALL remaining super-buckets `sb_1, ..., sb_{n_sb-1}` along the chord (`0` hits lost on multi-SB chords; every intersected SB is explicitly priced at `c_sb_hdr = 12.0` Ir/SB, simulated live in `verify_kernels_symbolic`!).

```asm
; kernel: rev_sb_probe_fast items=1 charged=18.0 key=c_tile_hdr
sub x0, x_k0, x_p
mul x1, x0, x_invp
umulh x2, x0, x_mu
lsr x2, x2, #4
msub x3, x2, x_p, x0
lsr x4, x3, #5
sub x5, x_end, x3
cmp x5, x_wmax
csel x5, x_wmax, x5, hi
lsl x6, x4, #3
add x7, x_off, x6
ldp x8, x9, [x7]
cmp x8, x9
b.hs .Lempty
add x_sb, x_sb, #1
add x_ptr, x_base, x8
ret
```

### 1.2 `rev_b2_prod_sb_hdr` (Band-2 super-bucket header with uniform byte offsets: `11 <= C_SB_HDR = 12.0` Ir/SB)

Outer loop state lives exclusively in callee-saved `x19` (`off` ptr), `x20` (`ent` base), `x21` (remaining SB count across the full chord), `x22` (`bptr` base, explicitly advanced `+256` B = 32 column pointers per SB on `.Lempty`/`.Lnext_sb`); flush limit lives in dedicated `x18` (never clobbered by `x5`!). Bucket boundaries `off[j], off[j+1]` are uniform 8B-aligned byte offsets; empty `sb_0` from `rev_sb_probe_fast` jumps directly to `.Lempty`/`.Lnext_sb` so all later super-buckets on the chord are scanned:

```asm
; kernel: rev_b2_prod_sb_hdr items=1 charged=12.0 key=c_sb_hdr
ldp x5, x6, [x19]
add x19, x19, #8
cmp x5, x6
b.eq .Lnext_sb
add x0, x20, x5
add x1, x20, x6
mov x2, x22
b .Lemit1
.Lempty:
.Lnext_sb:
add x22, x22, #256
subs x21, x21, #1
b.ne .Lsb
```

### 1.3 Band-1a & Band-2-apply zero-copy `uint16_t` 8x-unrolled + 1-record tail direct apply (`rev_b1a_apply_u8`: `30/8 = 3.75 <= C_HIT_LO_M4 = 3.80` Ir/hit; `rev_b1a_apply_tail1`: `6` Ir/rem hit, `<= +20` ops/bkt charged in `b1_lk = 44.0`)

Applies either a Band-1a `m`-partitioned `oid < 2^57` or a pre-rotated Band-2 59-bit word `(((m - k) & 3) << 57) | oid < 2^59` directly into contiguous `uint16_t acc[4 * 2^57]`. Exact ARMv8-A flag-safe pre/post wrap (`4` ops once per bucket): entry executes `subs x_cnt, x_cnt, #8; b.lo .Ltail_fix8`, and `.Ltail_fix8` executes `adds x_cnt, x_cnt, #8; b.eq .Ldone` before falling into `rev_b1a_apply_tail1` (`6` ops/item for `k mod 8 in {1..7}`). Worst-case overhead above `3.75` Ir/hit for ANY bucket length `k >= 0` is `<= 4 + 7 * (6 - 3.75) = 19.75 <= 20` ops/bkt (explicitly charged in `b1_lk = 44.0` for Band-1a; for Band-2 flushes of `>= 256` hits, `3.75 + 20/256 = 3.828`, charged `12.5 + 3.85 = 16.35`; simulated flag-exact across ALL `k = 0..25` in `audit_cm.py` with `0` overreads!):

```asm
; kernel: rev_b1a_apply_wrap8 items=1 charged=4.0 key=c_hit_lo_m4_wrap
subs x_cnt, x_cnt, #8
b.lo .Ltail_fix8
; (.Lb1a8 / .Lb1b_u8 loop runs below, then falls through to .Ltail_fix8)
.Ltail_fix8:
adds x_cnt, x_cnt, #8
b.eq .Ldone
```

```asm
; kernel: rev_b1a_apply_u8 items=8 charged=3.80 key=c_hit_lo_m4
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
subs x_cnt, x_cnt, #8
b.hs .Lb1a8
```

```asm
; kernel: rev_b1a_apply_tail1 items=1 charged=6.0 key=c_hit_lo_m4_tail
ldr x0, [x_p], #8
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Ltail1
```

### 1.4 Band-1b flag-exact 8x-unrolled + 1-record tail 8B packed `{m:2, oid:57}` direct apply (`rev_b1b_apply_u8`: `46/8 = 5.75 <= C_HIT_LO_U = 5.80` Ir/hit; `rev_b1b_apply_tail1`: `8.0` Ir/rem hit, `<= +20` ops/bkt charged inside `b1_lk = 34.0`)

Every Band-1b hit packs `(m << 57) | oid` into one 64-bit word (`oid < 2^57`, `m in {0..3}` -> 59 active bits), rotating bank `m -> (m - k) & 3` via `sub x_r, x_i, x_k57; and x_r, x_r, x_mask59` (`x_k57 = k << 57`, `x_mask59 = (1 << 59) - 1`). Reusing the identical flag-safe `rev_b1a_apply_wrap8` (`4` ops), `rev_b1b_apply_u8` processes 8 records in `46` opcodes (`5.75 <= 5.80` Ir/hit) and `rev_b1b_apply_tail1` processes the `k mod 8 in {1..7}` remainder at `8` opcodes/item (`5.75*k + 4 + 2.25*(k mod 8) <= 5.75*k + 19.75 <= 5.75*k + 20` total ops for every `k >= 0`, charging `+22 >= +20` ops/bkt in `b1_lk = 34.0`; simulated flag-exact `(C, Z)` across all `k = 0..25` in `verify_kernels_symbolic` with `0` overreads!):

```asm
; kernel: rev_b1b_apply_u8 items=8 charged=5.80 key=c_hit_lo_u
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
sub x4, x0, x_k57
and x4, x4, x_mask59
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
sub x5, x1, x_k57
and x5, x5, x_mask59
ldrh w9, [x_acc, x5, lsl #1]
add w9, w9, w_lp
strh w9, [x_acc, x5, lsl #1]
sub x6, x2, x_k57
and x6, x6, x_mask59
ldrh w10, [x_acc, x6, lsl #1]
add w10, w10, w_lp
strh w10, [x_acc, x6, lsl #1]
sub x7, x3, x_k57
and x7, x7, x_mask59
ldrh w11, [x_acc, x7, lsl #1]
add w11, w11, w_lp
strh w11, [x_acc, x7, lsl #1]
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
sub x4, x0, x_k57
and x4, x4, x_mask59
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
sub x5, x1, x_k57
and x5, x5, x_mask59
ldrh w9, [x_acc, x5, lsl #1]
add w9, w9, w_lp
strh w9, [x_acc, x5, lsl #1]
sub x6, x2, x_k57
and x6, x6, x_mask59
ldrh w10, [x_acc, x6, lsl #1]
add w10, w10, w_lp
strh w10, [x_acc, x6, lsl #1]
sub x7, x3, x_k57
and x7, x7, x_mask59
ldrh w11, [x_acc, x7, lsl #1]
add w11, w11, w_lp
strh w11, [x_acc, x7, lsl #1]
subs x_cnt, x_cnt, #8
b.hs .Lb1b_u8
```

```asm
; kernel: rev_b1b_apply_tail1 items=1 charged=8.0 key=c_hit_lo_u_tail
ldr x0, [x_p], #8
sub x4, x0, x_k57
and x4, x4, x_mask59
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Lb1b_tail1
```

### 1.5 Band-2 1-record-at-a-time 8B `(r<<59)|(m<<57)|oid` emit (`rev_b2_prod_emit_1`: `12 <= 12.5` Ir/hit; `+ 3.85 = 16.35` total)

Entered from `rev_b2_prod_sb_hdr` (`b .Lemit1`). Emits **one 8-byte record per iteration** (ZERO odd-length overread bug on any bucket size `0, 1, 2, ...`!), extracts `r = rec >> 59`, rotates bank `m -> (m - k) & 3` via `x3 = k << 57`, `x4 = (1 << 59) - 1`, appends to `bptr[r]`, checks flush limit `x18` (disjoint from `x5` and `x19..x22`!), and loops `cmp x0, x1; b.lo .Lemit1` falling through to `.Lnext_sb`. Boundary super-buckets (`w_sb / L_b = 18.4%` excess) add `3` range-check ops (`sub/cmp/b.hs`, charged `c_bd_chk = 3.0 * (1 + w_sb / L_b)`):

```asm
; kernel: rev_b2_prod_emit_1 items=1 charged=12.5 key=c_b2_emit
ldr x15, [x0], #8
lsr x9, x15, #59
sub x15, x15, x3
and x15, x15, x4
add x10, x2, x9, lsl #3
ldr x11, [x10]
str x15, [x11], #8
str x11, [x10]
cmp x11, x18
b.hs .Lflush
cmp x0, x1
b.lo .Lemit1
```

### 1.6 Bias-Initialized 32-Cell `uint16_t` Scan + Reset (`rev_scan_zero_u32`: `25/32 = 0.78125 <= 0.80` raw -> `0.78125 * 1.02761 + 0.12 = 0.9228 <= C_SCAN_ZERO = 0.95`; `W_tile = 416 = 13 * 32` exact!)

Each 64-bit word holds 4 `uint16_t` quarter-bit accumulators (`s <= 3,000`). Pre-initializing `S_g` lanes to `bias16 = (32768 - T) & 0x7fff` (`1 <= T <= 3,000`, broadcast across 4 `uint16_t` lanes in `x_b16`; `x_8000 = 0x8000_8000_8000_8000`) ensures every accumulated lane satisfies `(32768 - T) + s in [29768, 35767] subset [0, 65535]` with **ZERO 16-bit wraparound/carry between lanes** and `s >= T <==> (lane & 0x8000) != 0` (verified LIVE across all `T` and `s in [0, 3000]` in `verify_kernels_symbolic`). Signed-offset `[x_reg]` / `[x_reg, #16/#32/#48]` addressing only (ZERO post-indexed `ldp`/`stp`), and `orr` into `x4, x5, x10, x11` preserves all 32 original `S_g` scores in `x0..x3, x6..x9` on entry to `.Lhit32`. Since `W_tile = 416 = 13 * 32`, zero tail loop is needed:

```asm
; kernel: rev_scan_zero_u32 items=32 charged=0.80 key=c_scan_zero_raw
ldp x0, x1, [x_g]
ldp x2, x3, [x_g, #16]
stp x_b16, x_b16, [x_g]
stp x_b16, x_b16, [x_g, #16]
orr x4, x0, x1
orr x5, x2, x3
orr x4, x4, x5
ldp x6, x7, [x_g, #32]
ldp x8, x9, [x_g, #48]
stp x_b16, x_b16, [x_g, #32]
stp x_b16, x_b16, [x_g, #48]
add x_g, x_g, #64
orr x10, x6, x7
orr x11, x8, x9
orr x10, x10, x11
orr x4, x4, x10
tst x4, x_8000
b.ne .Lhit32
stp xzr, xzr, [x_f]
stp xzr, xzr, [x_f, #16]
stp xzr, xzr, [x_f, #32]
stp xzr, xzr, [x_f, #48]
add x_f, x_f, #64
subs x_cnt, x_cnt, #32
b.ne .Lscan32
```

Including tile padding `p2_mult = 1.02761` + setup/hit-branch amortization `0.12`: `0.78125 * 1.02761 + 0.12 = 0.9228 <= C_SCAN_ZERO = 0.95`.

### 1.7 1-Owner Flat-Array Counting-Sort / Prefix-Sum / Scatter Build (`rev_build_prod_pass12`: `22 + 0.50 = 22.50 <= C_LIST_A = 24.0`, `rev_build_prod_prefix`: `6 <= C_LBKT = 6.0`)

Because `Arch2-Prod` requires **ZERO `q_key` sorting**, Pass 1 walks `a_0`-chords in deterministic 4-bank order (`m = 0..3`, `oid = 0..V_a - 1`), advancing `ka = (ka + step) mod p` via the **4-opcode overflow-proof modular step** (`x_pms = p - step` precomputed once per prime; exact for ALL `0 <= ka, step < p < 2^64` with ZERO carry loss even for `p = 2^64 - 59`, tested LIVE!), writing `kbuf[idx] = ka`, and incrementing **64-bit `uint64_t` `x_hist[ka >> 5]`** (`ZERO` 32-bit overflow even at `p = 2^24`!). `rev_build_prod_prefix` writes BOTH immutable byte offsets `x_off[j]` AND mutable scatter cursors `x_pos[j]` in one 6-op pass (plus terminal sentinel `off[n_sb] = 8 * N_a_tbl`). Pass 2 executes **4 outer bank passes (`m = 0, 1, 2, 3`)**, initializing `x_oidm = m << 57` and `x_n2 = V_a` at the start of each bank `m` (so the bank transition after `oid = V_a - 1` in bank `m` explicitly resets `x_oidm = (m + 1) << 57` rather than overflowing into `(m << 57) | V_a`; `4` outer bank-init ops across `N_a^tbl = 2^58.84` entries = `< 10^-16` Ir/entry, verified LIVE across all 4 banks with `V_a = 40 << 2^57` in both `verify_kernels_symbolic` and `check_sigmaB_exactness`!), streaming `x_kbuf` in the identical order (ZERO `x_meta`, ZERO `x_imap`!), packing `((ka & 31) << 59) | x_oidm`, and scattering 8B into `x_out[pos[sb]]`. **Total simultaneous 1-owner RAM (`x_kbuf + x_out + x_off + x_pos + 416 MiB`) `<= 2^63.924 B < 2^64 B` for ALL Band-1a/1b/2 primes (`0` cross-owner build transfers, asserted live in `audit_cm.py`)**:

```asm
; kernel: rev_build_prod_pass12 items=1 charged=23.5 key=c_list_a_raw
cmp x0, x_pms
sub x5, x0, x_pms
add x6, x0, x_step
csel x0, x5, x6, hs
str x0, [x_kbuf], #8
lsr x1, x0, #5
ldr x2, [x_hist, x1, lsl #3]
add x2, x2, #1
str x2, [x_hist, x1, lsl #3]
subs x_n1, x_n1, #1
b.ne .Lp1
ldr x0, [x_kbuf2], #8
and x5, x0, #31
orr x8, x_oidm, x5, lsl #59
add x_oidm, x_oidm, #1
lsr x1, x0, #5
ldr x2, [x_pos_sb, x1, lsl #3]
add x4, x2, #8
str x4, [x_pos_sb, x1, lsl #3]
str x8, [x_out, x2]
subs x_n2, x_n2, #1
b.ne .Lp2
```

```asm
; kernel: rev_build_prod_prefix items=1 charged=6.0 key=c_lbkt
ldr x0, [x_hist], #8
str x_acc, [x_off], #8
str x_acc, [x_pos], #8
add x_acc, x_acc, x0, lsl #3
subs x_cnt, x_cnt, #1
b.ne .Lpre
```

Exactness of 4-orbit `sigma`-shared pure-product replay with `D_(omega,sigma)` (`239/239` exact-1), 8B `(r<<59)|(m<<57)|oid` flat-array build, overflow-proof 64-bit modular step, and odd/even/empty bucket replay is verified LIVE in `audit_cm.py` (`sigmaB_exactness`: 48 real BLS12-381 `w' = omega^2` ideals, `0` mismatches).

## 2. Exact 5-Limb Baseline Scalar ARMv8-A CIOS Montgomery Multiplication (`MONT5 = 270.0`)

One complete carry-exact, shift-included two-chain CIOS outer round executes **`48`** baseline scalar ARMv8-A opcodes (`5 x 48 = 240` across 5 unrolled rounds) + **`27`** load/zero-init/final-sub/store opcodes = **`267 <= 270.0` Ir** (`+3` Ir headroom; bit-exactness of this exact 267-opcode register machine is verified LIVE on 32 random 320-bit `(a, b, p)` triples inside `audit_cm.py`):

```asm
; kernel: mont5_cios_round items=1 charged=48.0 key=mont5_round
ldr x5, [x_b], #8
mul x10, x0, x5
mul x11, x1, x5
mul x12, x2, x5
mul x13, x3, x5
mul x14, x4, x5
adds x15, x15, x10
adcs x16, x16, x11
adcs x17, x17, x12
adcs x18, x18, x13
adcs x19, x19, x14
adcs x20, x20, xzr
adc x21, xzr, xzr
umulh x10, x0, x5
umulh x11, x1, x5
umulh x12, x2, x5
umulh x13, x3, x5
umulh x14, x4, x5
adds x16, x16, x10
adcs x17, x17, x11
adcs x18, x18, x12
adcs x19, x19, x13
adcs x20, x20, x14
adc x21, x21, xzr
mul x9, x15, x_n0
mul x10, x9, x22
mul x11, x9, x23
mul x12, x9, x24
mul x13, x9, x25
mul x14, x9, x26
adds xzr, x15, x10
adcs x15, x16, x11
adcs x16, x17, x12
adcs x17, x18, x13
adcs x18, x19, x14
adcs x19, x20, xzr
adc x20, x21, xzr
umulh x10, x9, x22
umulh x11, x9, x23
umulh x12, x9, x24
umulh x13, x9, x25
umulh x14, x9, x26
adds x15, x15, x10
adcs x16, x16, x11
adcs x17, x17, x12
adcs x18, x18, x13
adcs x19, x19, x14
adc x20, x20, xzr
```

```asm
; kernel: mont5_cios_epilogue items=1 charged=28.0 key=mont5_epi
ldp x0, x1, [x_a, #0]
ldp x2, x3, [x_a, #16]
ldr x4, [x_a, #32]
ldp x22, x23, [x_p, #0]
ldp x24, x25, [x_p, #16]
ldr x26, [x_p, #32]
mov x15, #0
mov x16, #0
mov x17, #0
mov x18, #0
mov x19, #0
mov x20, #0
subs x10, x15, x22
sbcs x11, x16, x23
sbcs x12, x17, x24
sbcs x13, x18, x25
sbcs x14, x19, x26
sbcs xzr, x20, xzr
csel x15, x10, x15, hs
csel x16, x11, x16, hs
csel x17, x12, x17, hs
csel x18, x13, x18, hs
csel x19, x14, x19, hs
stp x15, x16, [x_out, #0]
stp x17, x18, [x_out, #16]
str x19, [x_out, #32]
ret
```

Total exact `MONT5 = 5 * 48 + 27 = 267 <= 270.0` Ir.

## 3. Explicit Multi-Limb Cyclotomic Tower Norm Kernels, Parsed 4-Prime 62-Bit NTT (`fft_F`) & Garner CRT, Stage B, 2D Cofactorization, & Exact 2nd-Moment ECM Chebyshev

- **Explicit Schoolbook Multi-Limb Scalar ARMv8-A Tower Norm Derivation (`Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`, ZERO Karatsuba shortcut):**
  Explicit schoolbook scalar ARMv8-A `L_1 x L_2` limb multiplication uses `mul_LL(L_1, L_2) = 8 L_1 L_2 + 4 (L_1 + L_2)` (`mul/umulh/adds/adcs` carry chains). On the production sieving region (`|a_j|, |b_j| <= 2^8`), `a^3, b^3 in Z[y]` fit in 1 limb (`<= 27`b; `4` `8x8` 1-limb poly muls + `Phi_15` reduction = `400` Ir), `g(a, b) = a^3 - u w^2 b^3` fits in 2 limbs (`<= 94`b; `+96` Ir), and `f(a, b) = a^9 + (w+1) a^6 b^3 + a^3 b^6 - w b^9` fits in 3 limbs (`<= 128`b; `a^6, b^6` = `2 x 928 = 1,856` Ir, `a^9, b^9, a^6 b^3, a^3 b^6` = `4 x 1,152 = 4,608` Ir, combine = `128` Ir). Using pure schoolbook `mul_LL` across the 3-step tower norm `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`:
  1. **`_nk_g_only` (`2 -> 4 -> 8 -> 16` limbs):** `32` `2x2` (`1,536`) + `6` `4x4` (`960`) + `3` `8x8` (`1,728`) + adds/subs/setup/cold-scan (`1,544`) = **`5,768 <= 6,500` Ir charged**!
  2. **`_nk_f_only` (`3 -> 6 -> 12 -> 24` limbs):** reuses `a^3, b^3` -> `a^6, b^6, a^9, b^9, a^6 b^3, a^3 b^6` + combine (`6,592`) + 3-limb schoolbook tower norm (`4,608 + 4,032 + 3,744 = 12,384`) = **`18,976 <= 20,000` Ir charged**!
- **Parsed 62-Bit Shoup Radix-2 NTT Butterfly (`ntt62_butterfly`: `18` opcodes/2-elems/prime -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 fwd + 1 inv across 4 primes) & 4-Prime Garner CRT (`fft_F`):**
  `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`log2(NTT4_PROD) = 244.70`, all 4 verified prime LIVE). One exact baseline scalar ARMv8-A radix-2 Shoup butterfly (`w_s = floor(w * 2^64 / p)` precomputed, `t = w*v - umulh(v, w_s)*p` reduced to `[0, p)`) takes **`18` opcodes per 2 elements (`9.0` Ir/elem/prime -> `36.0` Ir/elem/transform across 4 primes -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 forward + 1 inverse NTTs in `fft_F`, `1.56x` headroom!)**; scalar `F_{p_i}` Barrett MAC in `PM-Basis` uses `mul/umulh` + 1-step Barrett `umulh/msub` + 2 conditional subs + add/sub/csel + 2 loads + loop = **`17 <= 24` opcodes** (verified bit-exactly LIVE on all 4 NTT primes in `kernels_symbolic`); 4-prime Garner CRT uses 3 Shoup mod-muls (`15` ops) + 3 multi-limb mul-adds (`<= 75` ops) = `<= 90 << 220` linear term in `fft_F` (`168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`; 2D `168 M_tot log2(M_tot) + 1100 M_tot + 1200`). Both `ntt62_butterfly` (64 random triples) and 4-prime Garner CRT (32 random 240-bit integers) are verified bit-exactly LIVE in `audit_cm.py`, alongside `srt_toy` (128 planted leaves, `0` mismatches):

```asm
; kernel: ntt62_butterfly items=2 charged=10.0 key=ntt62_elem_prime
ldr x0, [x_u]
ldr x1, [x_v]
ldp x2, x3, [x_w], #16
mul x4, x1, x2
umulh x5, x1, x3
msub x4, x5, x_p, x4
subs x6, x4, x_p
csel x4, x6, x4, hs
add x7, x0, x4
subs x8, x7, x_p
csel x7, x8, x7, hs
subs x9, x0, x4
add x10, x9, x_p
csel x9, x9, x10, hs
str x7, [x_u], #8
str x9, [x_v], #8
subs x_cnt, x_cnt, #1
b.ne .Lntt
```
- **Analytical Stage B (`cm_models.stage_b_analytic`):** `gbatch_cost(24197143, L) + 10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16)` with `MONT5 = 270.0`.
- **Live 2D-Histogram Cofactorization (`Q_F_SAFETY = 2.0`):** Under `H1`, `2.0 * rho_exact((bf + alpha_f) / lf)` at left endpoint `bf` gives `>= 2.27x` bin-midpoint `rho_db`, integrated LIVE across 401 bins up to `bg_max = 700`b (`c_cof_2d_upper = 2.0535` ipc).
- **Exact 2nd-Moment ECM Chebyshev (`cm_models.ecm_mean_var`, Stage-1 Collection) & Deterministic Single-Factor ECM (`cm_models.ecm_instr_det`, Stage-3 Descent/Boot):** In Stage-1 collection, each schedule runs `K_c` curves (`p_1`, `q_s=(1-p_1)^{K_c}<=0.001`, `p_s=1-q_s`); failed schedule cost `C_fail=K_c*per_c` (`Var=0`), final successful schedule has `S~Bin(K_c,p_1)|S>=1`, `C_succ=C_fail+S*redo_ex`, `E[S|S>=1]=K_c p_1/p_s`, `E[S^2|S>=1]=(K_c p_1(1-p_1)+(K_c p_1)^2)/p_s`, and with independent `F~Geom_0(p_s)`, `T=F*C_fail+C_succ` gives exact `E[T]=(q_s/p_s)C_fail+E[C_succ]`, `Var(T)=(q_s/p_s^2)C_fail^2+Var(C_succ)` (Cauchy-Schwarz `sd_rel=13*sqrt(Var_f)+16*sqrt(Var_g)` yields Chebyshev `< 2^-74.4 << 2^-70`). In Stage-3 descent & boot, `ecm_instr_det` halts on the first curve finding a factor (`1` redo at most), giving a **100% pointwise** deterministic upper bound `K_c*per_c + 1*redo_ex` with zero retention loss.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (1-Owner `< 2^64 B` 4-Vec Arch B Sign-Split Scalar ARMv8-A + Exact 1-Step Barrett `row_red5_z12`)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`, `tab_ratio_charged = 1.00` plain sparse with ZERO `TABM` subset-sum tables, `22` sequential 5-vector bundles of `160` B/slot across `2` `< 2^64 B` owners):
- **`unit_run5_z12` & `unit_run5_neg_z12` (`43` opcodes / `5` entries = `8.60 <= 8.80` Ir/vec-entry):** Each of the 6 hexagonal buckets `{chi^0..chi^5}` splits unit entries `+-chi^b in U_12` into a positive sub-list (`+chi^b`, `adds/adcs/adc`) and a negative sub-list (`-chi^b = chi^{b+6}`, `subs/sbcs/sbc`), accumulating into the signed 5-limb bucket register `(x10..x14)`:

```asm
; kernel: unit_run5_z12 items=5 charged=8.80 key=unit_run5
ldp x0, x1, [x_v0, #0]
ldp x2, x3, [x_v0, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v1, [x_idx], #8
ldp x0, x1, [x_v1, #0]
ldp x2, x3, [x_v1, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v2, [x_idx], #8
ldp x0, x1, [x_v2, #0]
ldp x2, x3, [x_v2, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v3, [x_idx], #8
ldp x0, x1, [x_v3, #0]
ldp x2, x3, [x_v3, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v4, [x_idx], #8
ldp x0, x1, [x_v4, #0]
ldp x2, x3, [x_v4, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v0, [x_idx], #8
add x_ptr, x_ptr, #40
subs x_cnt, x_cnt, #5
b.gt .Lrun5
```

```asm
; kernel: unit_run5_neg_z12 items=5 charged=8.80 key=unit_run5_neg
ldp x0, x1, [x_v0, #0]
ldp x2, x3, [x_v0, #16]
subs x10, x10, x0
sbcs x11, x11, x1
sbcs x12, x12, x2
sbcs x13, x13, x3
sbc x14, x14, xzr
ldr x_v1, [x_idx], #8
ldp x0, x1, [x_v1, #0]
ldp x2, x3, [x_v1, #16]
subs x10, x10, x0
sbcs x11, x11, x1
sbcs x12, x12, x2
sbcs x13, x13, x3
sbc x14, x14, xzr
ldr x_v2, [x_idx], #8
ldp x0, x1, [x_v2, #0]
ldp x2, x3, [x_v2, #16]
subs x10, x10, x0
sbcs x11, x11, x1
sbcs x12, x12, x2
sbcs x13, x13, x3
sbc x14, x14, xzr
ldr x_v3, [x_idx], #8
ldp x0, x1, [x_v3, #0]
ldp x2, x3, [x_v3, #16]
subs x10, x10, x0
sbcs x11, x11, x1
sbcs x12, x12, x2
sbcs x13, x13, x3
sbc x14, x14, xzr
ldr x_v4, [x_idx], #8
ldp x0, x1, [x_v4, #0]
ldp x2, x3, [x_v4, #16]
subs x10, x10, x0
sbcs x11, x11, x1
sbcs x12, x12, x2
sbcs x13, x13, x3
sbc x14, x14, xzr
ldr x_v0, [x_idx], #8
add x_ptr, x_ptr, #40
subs x_cnt, x_cnt, #5
b.gt .Lrun5n
```

- **`row_red5_z12` (`67 <= 68.0` Ir/bucket-vec — Exact 1-Step Signed 5-Limb Barrett Reduction for ALL `|a_4| < 2^60`):**
  Handles ANY signed 5-limb bucket accumulator `A = (x10..x14)` with `|a_4| < 2^60` (`|A| < 2^316`, covering EVERY `|c| <= 255` row up to `+-2^59 * (R_BLS - 1) > 10^17 * (R_BLS - 1)` with ZERO carry overflow!). Computes branchless `U = |A|`, `top = (u_4 << 2) | (u_3 >> 62)`, `q_hat = umulh(top, mu_r)` where `mu_r = floor(2^318 / R_BLS) = 0x8d54253b7fb78ddf < 2^64` (`0 <= floor(U / R_BLS) - q_hat <= 2` unconditionally), `rem = U - q_hat * R_BLS in [0, 3 * R_BLS)`, two branchless conditional `- R_BLS` steps, and branchless sign fix:

```asm
; kernel: row_red5_z12 items=1 charged=68.0 key=row_red5_barrett
asr x15, x14, #63
eor x0, x10, x15
eor x1, x11, x15
eor x2, x12, x15
eor x3, x13, x15
eor x4, x14, x15
subs x0, x0, x15
sbcs x1, x1, x15
sbcs x2, x2, x15
sbcs x3, x3, x15
sbc x4, x4, x15
lsl x5, x4, #2
orr x5, x5, x3, lsr #62
umulh x5, x5, x_mu
mul x6, x5, x22
mul x7, x5, x23
mul x8, x5, x24
mul x9, x5, x25
umulh x16, x5, x22
umulh x17, x5, x23
umulh x18, x5, x24
umulh x19, x5, x25
adds x7, x7, x16
adcs x8, x8, x17
adcs x9, x9, x18
adc x19, x19, xzr
subs x0, x0, x6
sbcs x1, x1, x7
sbcs x2, x2, x8
sbcs x3, x3, x9
sbc x4, x4, x19
subs x6, x0, x22
sbcs x7, x1, x23
sbcs x8, x2, x24
sbcs x9, x3, x25
sbcs x16, x4, xzr
csel x0, x6, x0, hs
csel x1, x7, x1, hs
csel x2, x8, x2, hs
csel x3, x9, x3, hs
csel x4, x16, x4, hs
subs x6, x0, x22
sbcs x7, x1, x23
sbcs x8, x2, x24
sbcs x9, x3, x25
sbcs xzr, x4, xzr
csel x0, x6, x0, hs
csel x1, x7, x1, hs
csel x2, x8, x2, hs
csel x3, x9, x3, hs
orr x6, x0, x1
orr x7, x2, x3
orr x6, x6, x7
cmp x6, #0
csel x15, xzr, x15, eq
subs x6, x22, x0
sbcs x7, x23, x1
sbcs x8, x24, x2
sbc x9, x25, x3
tst x15, x15
csel x0, x6, x0, ne
csel x1, x7, x1, ne
csel x2, x8, x2, ne
csel x3, x9, x3, ne
stp x0, x1, [x_out, #0]
stp x2, x3, [x_out, #16]
ret
```

- **Sign-Split `gen_run1_z12` (`c_b > 0`) & `gen_run1_neg_z12` (`c_b < 0`) (`23 <= 29.0` Ir/coord -> conservative `4 * 29.0 = 116.0` Ir charged per non-unit term):** Each bucket splits non-unit scalar updates (`1 <= |c_b| <= 255`) into a positive sub-list (`gen_run1_z12`, two-pass `adds/adcs/adc`) and a negative sub-list (`gen_run1_neg_z12`, two-pass `subs/sbcs/sbc` on `|c_b|` into the same signed 5-limb accumulator `(x10..x14)`). Under 6-bucket hexagonal routing (`|t_b| = 1` -> `unit_run5`, `|t_b| >= 2` -> `gen_run1` in bucket `b`), each `|t_b| >= 2` term executes ONE 23-op update (`23 <= 29.0` Ir), and even under unsplit 4-coordinate routing (`c_0..c_3` in buckets `0..3`), all 4 basis coordinates cost at most `4 * 23 = 92 <= 4 * 29.0 = 116.0` Ir across `npt(c) >= 1` hexagonal terms (`k5_live = eu_lo * 8.80 + (1 - eu_lo) * 116.0 = 15.9001 <= k5_charged = 18.0` Ir):

```asm
; kernel: gen_run1_z12 items=1 charged=29.0 key=gen_run1_coord
ldp x0, x1, [x_v0, #0]
ldp x2, x3, [x_v0, #16]
mul x4, x0, x_c
mul x5, x1, x_c
mul x6, x2, x_c
mul x7, x3, x_c
adds x10, x10, x4
adcs x11, x11, x5
adcs x12, x12, x6
adcs x13, x13, x7
adc x14, x14, xzr
umulh x4, x0, x_c
umulh x5, x1, x_c
umulh x6, x2, x_c
umulh x7, x3, x_c
adds x11, x11, x4
adcs x12, x12, x5
adcs x13, x13, x6
adc x14, x14, x7
ldr x_v0, [x_idx], #8
ldr x_c, [x_coef], #8
subs x_cnt, x_cnt, #1
b.gt .Lgen1
```

```asm
; kernel: gen_run1_neg_z12 items=1 charged=29.0 key=gen_run1_neg_coord
ldp x0, x1, [x_v0, #0]
ldp x2, x3, [x_v0, #16]
mul x4, x0, x_c
mul x5, x1, x_c
mul x6, x2, x_c
mul x7, x3, x_c
subs x10, x10, x4
sbcs x11, x11, x5
sbcs x12, x12, x6
sbcs x13, x13, x7
sbc x14, x14, xzr
umulh x4, x0, x_c
umulh x5, x1, x_c
umulh x6, x2, x_c
umulh x7, x3, x_c
subs x11, x11, x4
sbcs x12, x12, x5
sbcs x13, x13, x6
sbc x14, x14, x7
ldr x_v0, [x_idx], #8
ldr x_c, [x_coef], #8
subs x_cnt, x_cnt, #1
b.gt .Lgen1n
```

- **`chi_mul_z12_fast` (`50` Ir), `row_tail_6bucket_z12` (`504` Ir = `12` signed unit sub-lists x max `4` tail entries x `8` ops + `2 * 50` `chi_mul` + `20` dispatch), 1-Owner Plain-Sparse 4-Vec Arch B Row Overhead (`c_row_chi = 1,000.0` Ir, `c_xown_live = 0.0` Ir, `n_owners = 1`, `te = 0`):** `6 * 67 + 504 + 16 = 922 <= 6 * 68 + 504 + 16 = 928 <= c_row_chi = 1,000.0` Ir worst-case; zero cross-owner transfers (`tab_ratio_charged = 1.00`, ZERO `TABM`).

## 5. Descent Kernels (Zero-Skip Complete Iterative Q32 `dsieve16`, Strict-64b `BOUND_CELL_HI`, & Bareiss LLL)

- **`dsieve16` on `S_16(R) = B_16(R)/{+-1}` (`R in [52, 92]`, Unconditional `m`-Truncated `(9/10)`-Lovasz Signed-64-Bit Q32 Theorem!):**
  1. **Setup (`PAIR_SETUP = 2^35.0` Ir, proved), Active Prefix `m`, & Unconditional Signed-64-Bit Proof on EVERY `(9/10)`-LLL Integer Lattice:** Bareiss `lll_16d_int(9,10)` reduces `L_{Q,Q'_t,P^e}` (`det=N(P)^e<=2^70`, exact HNF initial GSO `d_1=...=d_16=det^2<=2^140`, `<=13,815` swaps -> `2^31.20 << 2^35.0` Ir). `PAIR_SETUP` then computes the active dimension prefix `m = min{t in {0..16} : forall j >= t, d_{j+1} > R^2 * d_j}`. By Gram-Schmidt orthogonality `||z * Br||_2^2 = sum_{j=0}^{15} (d_{j+1}/d_j)*(z_j + sum_{r>j} mu_{r,j} z_r)^2 <= R^2`, downward induction from `j = 15` down to `m` forces `|z_j| < 1 => z_j = 0` identically on `100%` of `B_16(R)` (if `m = 0`, `B_16(R)` has `0` nonzero lattice points and returns immediately; otherwise `PAIR_SETUP` builds `MU, L2_LO, INV_HI` ONLY for `k in {0..m-1}` and DFS enters at `k = m - 1`, terminating on `k == m`, `cmp x20, x_m` in `dsieve16_dfs_ctrl`). Because `d_m <= R^2 * d_{m-1}`, `d_1 / d_0 = ||Br[0]||_2^2 >= 1`, and `(9/10)`-LLL enforces `d_{j+2}/d_{j+1} >= (9/10 - 1/4)*(d_{j+1}/d_j) = (13/20)*(d_{j+1}/d_j)` and `|mu_{k,i}| <= 1/2`, **EVERY `k in {0..m-1}` on EVERY `(9/10)`-LLL integer lattice (`R <= 128`, including `diag(1,...,1,p)` for arbitrary `p`!) unconditionally satisfies**:
     - `1/642 < (13/20)^15 <= d_{k+1}/d_k <= (20/13)^15 * R^2 < 642 * R^2 <= 10,518,528 < 2^23.4`,
     - `||Br[k]||_2^2 = d_{k+1}/d_k + sum_{i<k} mu_{k,i}^2 (d_{i+1}/d_i) < 4.75 * 642 * R^2 < 3,100 * R^2 => |Br[k][t]| < 7,127 < 2^12.8`,
     - `L2_LO[k] = floor(2^32 * d_{k+1}/d_k) < 2^55.4 << 2^63`, `INV_HI[k] = isqrt(ceil(2^64 * d_k/d_{k+1})) + 1 < 2^37 << 2^63`, `|MU[k][i]| <= 2^31`, `|z_k| < 2^22`, and `sq < 2^56 << 2^63` — **PROVED UNCONDITIONALLY WITH ZERO SIGNED-64-BIT OVERFLOWS** (verified live in `pm_adaptive.py` on `diag(1,...,1,2^70+15)`, `diag(64,...,64,2^70+15)` at `R = 128`, `diag(1,2^70+15,...,2^70+15)`, and `diag(2^70+15,...,2^70+15)`!).
     On CRT lattice `L_{Q, Q'_t}`, ONLY the single blinding ideal `P^1 = Q'_t` (`e = 1`, which divides `100%` of `L_{Q, Q'_t}` on side `f` and subtracts `floor(log2 q'_t) >= 11` bits unconditionally) is skipped (saving `PAIR_SETUP` + leaf visits relative to the full unskipped `dsieve16_wc_pc` upper bound; any vector also vanishing at a sibling root `r_oth != r'_t` above the same `O_H`-prime satisfies `q'_OH | gcd_{O_H}(a,b)` and is non-coprime, verified live!). Every other prime ideal and prime power `16 <= N(P)^e <= 2^{70}` is sieved completely (`1.00 * WC_LLL99(R, det)` under `H2`, `8/8` live proved primes `0.197..0.315 < 0.45 << 1.00`).
  2. **Complete Iterative Q32 DFS with Incremental `V_stk` Stack & Chord Initialization (`w_k = 172 + 6*k` Ir at level `k >= 1` + `260.0` Ir at `k = 1` chord entry, ZERO floats, opcode-counted live):** Four parsed ARMv8.0-A blocks cover EVERY internal DFS, vector-stack, and chord-transition instruction (`V_stk[m]` zeroed once in `8` `stp` ops at root): `dsieve16_q32_node` (`55` ops: `cbz`+`clz`+5-iter `udiv` NR `isqrt(x0)+1`, proved for all `x0 in [0,2^46]` + `[lo,hi]` + `sq`), `dsieve16_dfs_ctrl` (`31` ops: sibling step, `stp/ldp` stack frame, `asum`/`nz`, child loads, `cmp x20, x_m`), `dsieve16_vec_madd` (`40` ops: `V_stk[k-1] = V_stk[k] + z_k * Br[k]` on every push `k in 1..m-1`), `dsieve16_dfs_cstep` (`6*k+4` ops: incremental `C_stk[k-1][0..k-1]` update via `msub`), plus `<=2` convex boundary rejects (`<=40` ops) -> `55+31+40+40+6+6*k = 172+6*k` Ir/node at `k in 1..m-1`, PLUS at `k = 1 -> 0` chord entry (`N_1` nodes) an extra `260.0` Ir (`55`-op `k=0` interval setup + `40`-op `dsieve16_vec_madd` chord init `v = V_stk[0] + (lo_0 - 1) * Br[0]` + `12`-op pop + up to `2` boundary leaf rejects `2 * (40 + 32) = 144` ops = `251 <= 260.0` Ir/chord!). Simulated bit-for-bit in `_fp_q32_int` with ZERO `z * Br` shortcut (`ops_nonhit <= sum_{k=1}^{m-1} (172 + 6*k + 260*[k==1]) N_k + 140` asserted live!).

```asm
; kernel: dsieve16_q32_node items=1 charged=55.0 key=dsieve16_node
cbz x0, .Lnr_zero
clz x1, x0
mov x2, #65
sub x1, x2, x1
lsr x1, x1, #1
mov x2, #1
lsl x2, x2, x1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
add x2, x2, #1
b .Lnr_done
.Lnr_zero:
mov x2, #1
.Lnr_done:
mul x3, x2, x7
umulh x4, x2, x7
lsr x3, x3, #16
orr x3, x3, x4, lsl #48
add x3, x3, #1
add x3, x3, x8
sub x9, x10, x3
add x9, x9, x11
asr x9, x9, #32
add x12, x10, x3
asr x12, x12, #32
cmp x13, #0
csel x9, xzr, x9, eq
lsl x14, x9, #32
sub x14, x14, x10
cmp x14, #0
sub x15, xzr, x14
csel x14, x15, x14, lt
subs x14, x14, x8
csel x14, xzr, x14, lt
mul x15, x14, x14
umulh x16, x14, x14
lsr x15, x15, #32
orr x15, x15, x16, lsl #32
mul x16, x15, x17
umulh x18, x15, x17
lsr x16, x16, #32
orr x16, x16, x18, lsl #32
subs x19, x0, x16
b.lo .Lprune
```

```asm
; kernel: dsieve16_dfs_ctrl items=1 charged=31.0 key=dsieve16_ctrl
cmp x9, x12
b.gt .Ldfs_pop
stp x9, x12, [x21, #0]
stp x0, x13, [x21, #16]
str x8, [x21, #32]
add x21, x21, #40
cmp x9, #0
sub x22, xzr, x9
csel x22, x22, x9, lt
add x8, x8, x22
cmp x9, #0
mov x23, #1
csel x23, x23, xzr, gt
orr x13, x13, x23
sub x20, x20, #1
add x24, x25, x20, lsl #7
ldr x10, [x24, x20, lsl #3]
ldr x7, [x26, x20, lsl #3]
ldr x17, [x27, x20, lsl #3]
mov x0, x19
b .Ldfs_enter
.Ldfs_pop:
add x20, x20, #1
cmp x20, x_m
b.eq .Ldfs_done
sub x21, x21, #40
ldp x9, x12, [x21, #0]
ldp x0, x13, [x21, #16]
ldr x8, [x21, #32]
add x9, x9, #1
b .Ldfs_sib
.Ldfs_done:
ret
```

```asm
; kernel: dsieve16_dfs_cstep items=1 charged=6.0 key=dsieve16_cstep
ldr x22, [x24, x23, lsl #3]
ldr x28, [x29, x23, lsl #3]
msub x22, x9, x28, x22
str x22, [x30, x23, lsl #3]
subs x23, x23, #1
b.ne .Lcstep_loop
```

  3. **Exact-integer leaf vector step (`v += b_0` after priced `k=0` chord init `v = V_stk[0] + (lo_0 - 1)*b_0`) + norm check + bijective `rank_lex` + bucket write (`316 <= 340.0` Ir/leaf):** Upon entering any non-empty `k = 0` chord `[lo_0, hi_0]`, `dsieve16_vec_madd` (`40` ops, priced in `k1_rej = 260.0` above!) initializes `v = V_stk[0] + (lo_0 - 1) * b_0`; then each `z_0` step executes the full 40-op `dsieve16_leaf_vec` block below (`8 ldp` of `b_0[0..15]` + `8 ldp` of `v[0..15]` + `16 add` + `8 stp` of `v[0..15]`; boundary `z_0` values with `||v||_2^2 > R^2` are priced in `k1_rej = 260.0`, while qualifying hits pay `c_leaf = 340.0`) + exact integer `sum v[t]^2 <= R^2` (`32` ops) + canonical sign-flip (`16` ops) + bijective dense ball rank via `0.94 GiB` flat `P_d[rem*(2R+1) + (v_d+R)]` table (`16` uniform levels `d = 0..15` x `10` ops = `160` ops, `dsieve16_rank_step` below loads `v_d`, computes `rem*(2R+1)+(v_d+R)` explicitly from live `x7 = rem`, indexes the 128-bit prefix entry, updates `(x5, x6)` and `rem = x7`, and advances `x_pd` to the next level table — simulated bit-for-bit on all `552/552` hits in `audit_cm.py`!) + owner-chunk route (`56` ops) + `uint16_t` score add (`12` ops) -> `40 + 32 + 16 + 160 + 56 + 12 = 316 <= 340.0` Ir/leaf (`+24` Ir/leaf extra loop/setup margin)!

```asm
; kernel: dsieve16_vec_madd items=1 charged=40.0 key=dsieve16_vmadd
ldp x0, x1, [x_row, #0]
ldp x2, x3, [x_base, #0]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #0]
ldp x0, x1, [x_row, #16]
ldp x2, x3, [x_base, #16]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #16]
ldp x0, x1, [x_row, #32]
ldp x2, x3, [x_base, #32]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #32]
ldp x0, x1, [x_row, #48]
ldp x2, x3, [x_base, #48]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #48]
ldp x0, x1, [x_row, #64]
ldp x2, x3, [x_base, #64]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #64]
ldp x0, x1, [x_row, #80]
ldp x2, x3, [x_base, #80]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #80]
ldp x0, x1, [x_row, #96]
ldp x2, x3, [x_base, #96]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #96]
ldp x0, x1, [x_row, #112]
ldp x2, x3, [x_base, #112]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #112]
```

```asm
; kernel: dsieve16_leaf_vec items=1 charged=40.0 key=dsieve16_leaf_v
ldp x0, x1, [x_b0, #0]
ldp x2, x3, [x_v, #0]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #0]
ldp x0, x1, [x_b0, #16]
ldp x2, x3, [x_v, #16]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #16]
ldp x0, x1, [x_b0, #32]
ldp x2, x3, [x_v, #32]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #32]
ldp x0, x1, [x_b0, #48]
ldp x2, x3, [x_v, #48]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #48]
ldp x0, x1, [x_b0, #64]
ldp x2, x3, [x_v, #64]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #64]
ldp x0, x1, [x_b0, #80]
ldp x2, x3, [x_v, #80]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #80]
ldp x0, x1, [x_b0, #96]
ldp x2, x3, [x_v, #96]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #96]
ldp x0, x1, [x_b0, #112]
ldp x2, x3, [x_v, #112]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #112]
```

```asm
; kernel: dsieve16_rank_step items=1 charged=10.0 key=dsieve16_rank
ldr x1, [x_v], #8
add x2, x1, x_R
mul x0, x7, x_2Rp1
add x2, x0, x2
add x2, x_pd, x2, lsl #4
ldp x3, x4, [x2, #0]
adds x5, x5, x3
adc x6, x6, x4
msub x7, x1, x1, x7
add x_pd, x_pd, x_stride
```

  4. **Certified 2-stage strict signed-64b 60-bit-mantissa (`smulh`/`umulh`) UPPER-BOUND ball `BOUND_CELL_HI` with Deterministic Near-Cancellation Routing + Prime-Power Valuation Recovery (`BOUND_CELL_g = 1,800` Ir on all cells + `BOUND_CELL_f = 1,800` Ir on `g`-survivors, `p_nc = 2^-24` explicitly priced, ZERO floats, ZERO 64-bit overflow):** 16-op pre-shift aligns `(a, b)` to a common scale `sh` (`|a_s|, |b_s| < 2^55`) against pure-`Decimal(50)` Q58 `Z_TAB` (`err <= 0.51 * 2^-58`). Every disk maintains `|re|,|im|,rad < 2^60` with left-normalization (`clz`), so `bound_cell_ball_mul` (`48 <= 52` Ir using baseline scalar ARMv8.0-A `smulh`/`umulh`) preserves `>= 50` bits of relative precision (`rad <= 2^-48 M_poly`) across all 11 Horner steps with ZERO 64-bit overflow. At output, shifting right gives 30-bit signed `r30, i30` (`|r30|,|i30| < 2^30`), `e30 = ((rad+(1<<30)-1)>>30)+2`, `sq = isqrt(r30^2+i30^2)`, and `m_hi = sq + 1 + e30`, yielding the **UNCONDITIONAL upper bound** `L_hi = sum_{k=0}^3 2 * log2_ceil(m_hi * 2^{sh+30}) >= log2|N_exact|` (via pure-`Decimal(50)` `ROUND_CEILING` `_LUT_HI`). **4-op near-cancellation guard per complex place (`add x_ep1, x_e30, #1`; `sub x_mlo, x_sq, x_ep1`; `cmp x_mlo, x_ep1, lsl #6`; `b.lt .Lnc`, `+16` ops across 4 places):** let `ep1 = e30 + 1 >= 3`. Whenever `sq - ep1 >= 64 * ep1` at all 4 places (`nc == 0`, i.e., `sq >= 65 * (e30 + 1)`), exact scaled magnitude `|z_k| >= sq - e30 > sq - ep1 >= 64 * ep1` and `m_hi = sq + ep1` satisfy `m_hi / |z_k| < (sq + ep1) / (sq - ep1) <= 66 / 64` STRICTLY for EVERY integer `sq, e30 >= 0` (at boundary `sq = 65 * ep1`, `(66 * ep1) / (64 * ep1) = 66 / 64` exact!), and upward `_LUT_HI` Q16 ceiling adds at most `log2(65/64) + 2^{-16}` per place (`8` conjugate places total), unconditionally proving **`0 <= L_hi - log2|N_exact| <= 8 * (log2(66/64) + log2(65/64) + 2^{-16}) = 0.534219 <= 0.5343 < 0.625` bits** (`+` worst-case Q4 `1/16`-bit `dsieve16` log-score rounding `<= 41.0`b = `<= 41.7`b, leaving `> 54.3`b of unused slack inside the already-paid `T - T0 = 96.0`b threshold); whenever `sq - ep1 < 64 * ep1` at ANY place (`nc == 1`), the cell is **NEVER dropped** — it is **deterministically routed straight to Stage 2 + exact integer `NORMS` (`pre_gate`)** (priced explicitly at `p_nc = 2^-24` in `eval_descent_live`, verified LIVE on a literal `-57.3`b near-cancellation witness!). Combined with prime-power sieving (`Q^e` for ALL `e >= 2, N(Q)^e <= 2^{71}` + exact `2`-adic `g`-side `P_i^k | a` `k=1..6` weights `(12,12,12,12,12,4)` recovering all `128` bits of `v_2(N_g)` with `0` deficit, explicitly priced in `dsieve16_wc_pc` & `PAIR_SETUP`), every smooth cell with uncredited prime-power excess `<= 54.3`b above `2^{71}` is retained with `0` false negatives, while the `< 2^{-35}` tail of cells with prime-power excess `> 54.3`b above `2^{71}` (`prod N(Q)^{e_Q} > 2^{125.3}`) is absorbed with `> 10^9 x` margin inside the explicitly charged `ret_node = 0.90` (`10%` retention discount) under `H2`. **Stage 1 (`g_3` on all cells, 5 muls/emb):** `16 + 4*(104 + 5*48 + 32 + 44) = 1,696 <= BOUND_CELL_g = 1,800` Ir. **Stage 2 (`f_9` on `g`-survivors, 6 muls/emb):** `4*(6*48 + 96 + 44) = 1,712 <= BOUND_CELL_f = 1,800` Ir (combined `3,408 <= 3,600` Ir). Exact integer norm caps `|N_g| <= 2^{ng_max}, |N_f| <= 2^{nf_max}` are enforced inside `NORMS` (`pre_gate`) at zero extra cost (`36/36` production lattices lie `21.55..43.47`b below!). SWAR `SCAN+ZERO=6`, `NORMS=2^27`, `TDIV=2^25`, `IDEAL=2^35`, `MONT5=270`, MR `<=2^-80`.

```asm
; kernel: bound_cell_ball_mul items=1 charged=52.0 key=bound_ball_mul
smulh x8, x0, x3
smulh x9, x1, x4
sub x8, x8, x9
smulh x10, x0, x4
smulh x11, x1, x3
add x10, x10, x11
cmp x0, #0
sub x12, xzr, x0
csel x12, x12, x0, lt
cmp x1, #0
sub x13, xzr, x1
csel x13, x13, x1, lt
add x12, x12, x13
cmp x3, #0
sub x13, xzr, x3
csel x13, x13, x3, lt
cmp x4, #0
sub x14, xzr, x4
csel x14, x14, x4, lt
add x13, x13, x14
umulh x14, x12, x5
umulh x15, x13, x2
add x14, x14, x15
umulh x15, x2, x5
add x14, x14, x15
add x14, x14, #5
cmp x8, #0
sub x9, xzr, x8
csel x9, x9, x8, lt
cmp x10, #0
sub x11, xzr, x10
csel x11, x11, x10, lt
cmp x9, x11
csel x9, x9, x11, hi
cmp x9, x14
csel x9, x9, x14, hi
cbz x9, .Lbn0
clz x11, x9
subs x11, x11, #4
csel x11, xzr, x11, lt
lsl x8, x8, x11
lsl x10, x10, x11
lsl x14, x14, x11
b .Lbn1
.Lbn0:
mov x11, #0
.Lbn1:
add x16, x6, x7
add x16, x16, #64
sub x16, x16, x11
```

- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every reduction runs `_dw_int(B, 11, 40)` -> `_dw_int(B_1, 99, 100)` with full pre-Lovasz size-reduction (`|2*lam[i][j]| <= d[j+1]`). Swaps shrink `D = prod_{k=1}^{15} d_k` by `< delta` with `max_k d_k <= D_max` non-increasing; `L_lam = L_b + L_d`, `L_num = 2 L_d + L_b`, Stage-1 bound `d_k <= 40^{k(16-k)/2} d_{16}^{k/16}` (`||b_i||^2 <= 5 D_{max,2}`, asserted in `verify_kernels_symbolic`):
  1. **Special-`q` CRT `L_{Q, Q'_t}` (`det = q * q' < 2^280`, exact HNF initial `d_1 = ... = d_{16} = det^2 < 2^560`, `D_0 < 2^8400`, `LLL_NODE = 2^40`):** Stage 1 (`S_1 <= 4,509`): `<= 2^32.74` Ir; Stage 2 (`S_2 <= 415,436`): `<= 2^39.27` Ir -> **`2^39.28 << 2^40`** Ir (live tested at `278`b CRT product `q * 8191`!).
  2. **381b `G_T` Boot (`d_{16} <= 2^3384`, `LLL_BOOT = 2^48`):** Stage 1 (`S_1 <= 36,076`): `<= 2^43.16` Ir; Stage 2 (`S_2 <= 1,875,188`): `<= 2^46.30` Ir -> **`2^46.45 << 2^48`** Ir.

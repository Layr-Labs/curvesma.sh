# Compiled ARMv8.1-A Scalar Kernels (Part 1: Stage 1a Hybrid Split-Slab Sieve & Zero-Carry SWAR Merge-Scan, Kernels 1–7)

All 16 kernels across `KERNELS.md` (Kernels 1–7: **`583` actual scalar instructions**, `73` disclosed `staddb`) and `KERNELS_STAGE1B_LA.md` (Kernels 8–16: **`942` actual scalar instructions**, `0` `staddb`) are compiled verbatim from standalone C translation units with Apple `clang` targeting strict scalar ARMv8.1-A general-purpose registers:

```text
clang -O3 -march=armv8.1-a -mgeneral-regs-only -fno-vectorize -fno-slp-vectorize -fomit-frame-pointer -ffixed-x18 -S
```

**Global ISA & Register Hygiene Audit (All 16 Functions Combined = `1525` Executable Scalar Instructions):**
- SIMD/FP/SVE/crypto/CRC registers or mnemonics (`v0-v31`, `q/d/s/z/p`, `ld1/st1`, `pmull`, `aes`, `sha`, `crc`): **`0`** (`simd = 0`).
- Hardware divisions (`udiv`/`sdiv`) or out-of-line calls (`bl`/`blr`): **`0`**.
- **In-loop stack spills ANYWHERE across all 16 functions (both fast and slow paths):** **`0`** (every `[sp` access is strictly function-prologue/epilogue callee-save `stp`/`ldp` or prologue ABI stack-argument load prior to loop entry).
- ARMv8.1-A LSE atomic byte adds (`staddb`): **`73`** total — ALL inside `_emit_ptr64_x64` for lock-free banked slab log accumulation, explicitly disclosed and charged at `1.0 Ir` each (`100/64 = 1.5625 Ir/hit`).

## Part 1 Dynamic Instruction Count Summary (Kernels 1–7)

| # | Symbol | Role | Static Insts | In-Loop Spills | Compiled Dynamic `Ir` | Charged in Ledger |
|---:|---|---|---:|---:|---|---|
| 1 | `_tx_pack_u64_window` | Tier-1 (`p <= 2^{55}`, `W_B=1`) `M=248` Fixed-VA `ptr64` TX pack (`x64`) | 168 | **0** | `7` col + `10` non-empty + `69/64 = 1.078125 Ir/hit` + tail | Exact literal Poisson (`10 + 69*q64 + tail`) |
| 2 | `_tx_add_ptr64_chord_x16` | Tier-2 (`W_B=16`) & Tier-3 (`W_B=32`) `M=2` Half-Slab (`STRIDE=2^{56}`) multi-bucket `base+off` TX (`x16`+binary tail) | 98 | **0** | `6` chord + `8` empty bkt; non-empty `11 + 37*q16 + T16[r]` (`37/16 = 2.3125 Ir/hit`) | Exact literal Poisson (124-col Half-Slab `n_b = (L_124 - 1)/W + 1`) |
| 3 | `_chord_probe_step` | Per-prime branchless modular chord step + exact `r0`-aware `n_b` (`((r0+L-1)>>sh)+1`) + slice empty check | 19 | **0** | **`19 Ir/prime`** (0 branches, 0 spills) | **`20.0 Ir/prime`** (`+1 Ir` pad across all 3 tiers) |
| 4 | `_emit_ptr64_x64` | Banked `ptr64` slab log scatter (`73` `staddb`, `x64` unrolled, both `M=248` & `M=2`) | 126 | **0** | **`1.5625 Ir/hit`** (`100/64`, zero `add`) | **`1.5625 Ir/hit`** (256 KiB batched ring) |
| 5 | `_tx_copy_prebias_x256` | `M=248` `S_lo_g` 256B TX stream + fused prebias reset (`prebias_f`) | 54 | **0** | **`52/256 = 0.203125 Ir/cell`** (`LBB1_1`) | **`0.203125 Ir/cell`** |
| 6 | `_col_prebias_x256` | 256-byte column pre-bias / zeroing (`stp` x16) | 21 | **0** | **`19/256 = 0.07421875 Ir/cell/pass`** (`LBB2_1`) | `S_lo_f` (`1.0x`) + `S_hi_f` (`1.25x`) + `S_hi_g` interior overhang (`0.25x`) |
| 7 | `_merge_scan_swar_128_exact` | `M=2` Half-Slab fused zero-carry 64-bit `ADD` + 128B SWAR scan + `S_hi_g` zero | 97 | **0** | Fast: **`64/128 = 0.500000 Ir/cell`**; Slow: `171` + `29` Pass-2 cross-owner `f`-probe = `200 Ir/surv` | **`0.905451 Ir/cell`** (live `0.901451` + `0.004` pad) |

### Exact Verification of Zero-Carry 64-Bit `ADD`, Sequential `g`-then-`f` Reuse & `0.25x` Interior Overhang Zeroing
- **Deterministic Zero Inter-Byte Carry:** Because `0 <= prebias_i + score_lo_i + score_hi_i <= 255` unconditionally for every cell (`max_byte <= 255` via per-column $N_{g,\max}(b)$ enforcement), bytewise sum `lo[w] + hi[w]` in standard 64-bit integer addition (`ADD Xd, Xn, Xm`) never generates a carry across byte boundaries!
- **Sequential `g`-then-`f` Slab Reuse & `0.25x` Interior Overhang Proof:** Each Half-Slab (`186 * 2^{56} = 11.6250 EiB/owner << 16 EiB`) processes `g` first, records `g`-survivors into `gsB` (`2^{52.9235} B`), and zeroes `S_hi_g`, then reuses the exact same `11.6250 EiB` Half-Slab for `S_hi_f` (and `M=248` reuses `0.0584 EiB/owner` for `S_lo_f`, zero simultaneous residency!). Permanently unread outer guard columns `[-31..-1]` and `[124..154]` are write-only trash columns never read and never zeroed; only interior right-overhang columns `[L_half .. min(123, L_half+30)]` (dirtied when `L_half < 124` and reused as active columns on subsequent chords) are zeroed — shell-weighted exact ratio is `0.2445x <= 0.2500x` charged!
- **Fast path (`LBB2_2 -> LBB2_3 -> LBB2_4`):** `16 ldp` + `16 add` + `15 orr` + `tst` + `b.ne` + `mov` + `add` + `8 stp xzr` + `add` + `mov` + `sub` + `cbnz` = **`63` actual executable instructions per 128 bytes**, conservatively charged at **`64/128 = 0.5000000 Ir/cell`** (`+1 Ir/128B` pad).
- **Slow path (`LBB2_5 -> LBB2_8 -> LBB2_10`):** Zeroing of `hi[0..15]` occurs strictly *after* the branch in `LBB2_4`, so the cold slow path re-reads exact `lo[w] + hi[w]` with **zero in-loop stack spills** (`171 Ir` + `29 Ir` Pass-2 cross-owner `S_lo_f + S_hi_f` probe = `200 Ir/survivor`). At live `log2(s_g) = -13.933912`, `slab_live = 52/256 + 2.5*(19/256) + 64/128 + 200.0 * 2^{-13.933912} = 0.901451 Ir/cell`, charged **`0.905451 Ir/cell`** (`+0.004` pad).

## Verbatim Compiled ARMv8.1-A Assembly (Kernels 1–7, Zero Ellipses)

```asm
_tx_pack_u64_window:                    ; @tx_pack_u64_window
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB0_19
; %bb.1:
	mov	x8, #0                          ; =0x0
	b	LBB0_4
LBB0_2:                                 ;   in Loop: Header=BB0_4 Depth=1
	str	x9, [x0, x8, lsl #3]
LBB0_3:                                 ;   in Loop: Header=BB0_4 Depth=1
	add	x8, x8, #1
	cmp	x8, x3
	b.eq	LBB0_19
LBB0_4:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_6 Depth 2
	ldr	w10, [x2, x8, lsl #2]
	cbz	w10, LBB0_3
; %bb.5:                                ;   in Loop: Header=BB0_4 Depth=1
	ldr	x9, [x0, x8, lsl #3]
	cmp	w10, #64
	b.lo	LBB0_7
LBB0_6:                                 ;   Parent Loop BB0_4 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x11, x12, [x1]
	ldp	x13, x14, [x1, #16]
	ldp	x15, x16, [x1, #32]
	ldp	x17, x4, [x1, #48]
	stp	x11, x12, [x9]
	stp	x13, x14, [x9, #16]
	stp	x15, x16, [x9, #32]
	stp	x17, x4, [x9, #48]
	ldp	x11, x12, [x1, #64]
	ldp	x13, x14, [x1, #80]
	ldp	x15, x16, [x1, #96]
	ldp	x17, x4, [x1, #112]
	stp	x11, x12, [x9, #64]
	stp	x13, x14, [x9, #80]
	stp	x15, x16, [x9, #96]
	stp	x17, x4, [x9, #112]
	ldp	x11, x12, [x1, #128]
	ldp	x13, x14, [x1, #144]
	ldp	x15, x16, [x1, #160]
	ldp	x17, x4, [x1, #176]
	stp	x11, x12, [x9, #128]
	stp	x13, x14, [x9, #144]
	stp	x15, x16, [x9, #160]
	stp	x17, x4, [x9, #176]
	ldp	x11, x12, [x1, #192]
	ldp	x13, x14, [x1, #208]
	ldp	x15, x16, [x1, #224]
	ldp	x17, x4, [x1, #240]
	stp	x11, x12, [x9, #192]
	stp	x13, x14, [x9, #208]
	stp	x15, x16, [x9, #224]
	stp	x17, x4, [x9, #240]
	ldp	x11, x12, [x1, #256]
	ldp	x13, x14, [x1, #272]
	ldp	x15, x16, [x1, #288]
	ldp	x17, x4, [x1, #304]
	stp	x11, x12, [x9, #256]
	stp	x13, x14, [x9, #272]
	stp	x15, x16, [x9, #288]
	stp	x17, x4, [x9, #304]
	ldp	x11, x12, [x1, #320]
	ldp	x13, x14, [x1, #336]
	ldp	x15, x16, [x1, #352]
	ldp	x17, x4, [x1, #368]
	stp	x11, x12, [x9, #320]
	stp	x13, x14, [x9, #336]
	stp	x15, x16, [x9, #352]
	stp	x17, x4, [x9, #368]
	ldp	x11, x12, [x1, #384]
	ldp	x13, x14, [x1, #400]
	ldp	x15, x16, [x1, #416]
	ldp	x17, x4, [x1, #432]
	stp	x11, x12, [x9, #384]
	stp	x13, x14, [x9, #400]
	stp	x15, x16, [x9, #416]
	stp	x17, x4, [x9, #432]
	ldp	x11, x12, [x1, #448]
	ldp	x13, x14, [x1, #464]
	ldp	x15, x16, [x1, #480]
	ldp	x17, x4, [x1, #496]
	stp	x11, x12, [x9, #448]
	stp	x13, x14, [x9, #464]
	stp	x15, x16, [x9, #480]
	add	x1, x1, #512
	stp	x17, x4, [x9, #496]
	add	x9, x9, #512
	sub	x10, x10, #64
	cmp	x10, #63
	b.hi	LBB0_6
LBB0_7:                                 ;   in Loop: Header=BB0_4 Depth=1
	cmp	x10, #32
	b.hs	LBB0_13
; %bb.8:                                ;   in Loop: Header=BB0_4 Depth=1
	tbnz	w10, #4, LBB0_14
LBB0_9:                                 ;   in Loop: Header=BB0_4 Depth=1
	tbnz	w10, #3, LBB0_15
LBB0_10:                                ;   in Loop: Header=BB0_4 Depth=1
	tbnz	w10, #2, LBB0_16
LBB0_11:                                ;   in Loop: Header=BB0_4 Depth=1
	tbnz	w10, #1, LBB0_17
LBB0_12:                                ;   in Loop: Header=BB0_4 Depth=1
	tbz	w10, #0, LBB0_2
	b	LBB0_18
LBB0_13:                                ;   in Loop: Header=BB0_4 Depth=1
	ldp	x11, x12, [x1]
	ldp	x13, x14, [x1, #16]
	ldp	x15, x16, [x1, #32]
	ldp	x17, x4, [x1, #48]
	stp	x11, x12, [x9]
	stp	x13, x14, [x9, #16]
	stp	x15, x16, [x9, #32]
	stp	x17, x4, [x9, #48]
	ldp	x11, x12, [x1, #64]
	ldp	x13, x14, [x1, #80]
	ldp	x15, x16, [x1, #96]
	ldp	x17, x4, [x1, #112]
	stp	x11, x12, [x9, #64]
	stp	x13, x14, [x9, #80]
	stp	x15, x16, [x9, #96]
	stp	x17, x4, [x9, #112]
	ldp	x11, x12, [x1, #128]
	ldp	x13, x14, [x1, #144]
	ldp	x15, x16, [x1, #160]
	ldp	x17, x4, [x1, #176]
	stp	x11, x12, [x9, #128]
	stp	x13, x14, [x9, #144]
	stp	x15, x16, [x9, #160]
	stp	x17, x4, [x9, #176]
	ldp	x11, x12, [x1, #192]
	ldp	x13, x14, [x1, #208]
	ldp	x15, x16, [x1, #224]
	ldp	x17, x4, [x1, #240]
	stp	x11, x12, [x9, #192]
	stp	x13, x14, [x9, #208]
	stp	x15, x16, [x9, #224]
	add	x1, x1, #256
	stp	x17, x4, [x9, #240]
	add	x9, x9, #256
	tbz	w10, #4, LBB0_9
LBB0_14:                                ;   in Loop: Header=BB0_4 Depth=1
	ldp	x11, x12, [x1]
	ldp	x13, x14, [x1, #16]
	ldp	x15, x16, [x1, #32]
	ldp	x17, x4, [x1, #48]
	stp	x11, x12, [x9]
	stp	x13, x14, [x9, #16]
	stp	x15, x16, [x9, #32]
	stp	x17, x4, [x9, #48]
	ldp	x11, x12, [x1, #64]
	ldp	x13, x14, [x1, #80]
	ldp	x15, x16, [x1, #96]
	ldp	x17, x4, [x1, #112]
	stp	x11, x12, [x9, #64]
	stp	x13, x14, [x9, #80]
	stp	x15, x16, [x9, #96]
	add	x1, x1, #128
	stp	x17, x4, [x9, #112]
	add	x9, x9, #128
	tbz	w10, #3, LBB0_10
LBB0_15:                                ;   in Loop: Header=BB0_4 Depth=1
	ldp	x11, x12, [x1]
	ldp	x13, x14, [x1, #16]
	ldp	x15, x16, [x1, #32]
	ldp	x17, x4, [x1, #48]
	stp	x11, x12, [x9]
	stp	x13, x14, [x9, #16]
	stp	x15, x16, [x9, #32]
	add	x1, x1, #64
	stp	x17, x4, [x9, #48]
	add	x9, x9, #64
	tbz	w10, #2, LBB0_11
LBB0_16:                                ;   in Loop: Header=BB0_4 Depth=1
	ldp	x11, x12, [x1]
	ldp	x13, x14, [x1, #16]
	stp	x11, x12, [x9]
	stp	x13, x14, [x9, #16]
	add	x1, x1, #32
	add	x9, x9, #32
	tbz	w10, #1, LBB0_12
LBB0_17:                                ;   in Loop: Header=BB0_4 Depth=1
	ldp	x11, x12, [x1], #16
	stp	x11, x12, [x9], #16
	tbz	w10, #0, LBB0_2
LBB0_18:                                ;   in Loop: Header=BB0_4 Depth=1
	ldr	x10, [x1], #8
	str	x10, [x9], #8
	b	LBB0_2
LBB0_19:
	ret
	.cfi_endproc

_tx_add_ptr64_chord_x16:                ; @tx_add_ptr64_chord_x16
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB0_14
; %bb.1:
	mov	x8, #0                          ; =0x0
	ldr	x11, [x2]
	b	LBB0_3
LBB0_2:                                 ;   in Loop: Header=BB0_3 Depth=1
	add	x4, x4, x5
	mov	x11, x9
	cmp	x8, x3
	b.eq	LBB0_14
LBB0_3:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_5 Depth 2
	add	x8, x8, #1
	ldr	x9, [x2, x8, lsl #3]
	subs	x10, x9, x11
	b.eq	LBB0_2
; %bb.4:                                ;   in Loop: Header=BB0_3 Depth=1
	add	x11, x1, x11, lsl #3
	cmp	x10, #16
	b.lo	LBB0_6
LBB0_5:                                 ;   Parent Loop BB0_3 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x12, x13, [x11]
	ldp	x14, x15, [x11, #16]
	ldp	x16, x17, [x11, #32]
	ldp	x6, x7, [x11, #48]
	add	x12, x12, x4
	add	x13, x13, x4
	stp	x12, x13, [x0]
	add	x12, x14, x4
	add	x13, x15, x4
	stp	x12, x13, [x0, #16]
	add	x12, x16, x4
	add	x13, x17, x4
	stp	x12, x13, [x0, #32]
	add	x12, x6, x4
	add	x13, x7, x4
	stp	x12, x13, [x0, #48]
	ldp	x12, x13, [x11, #64]
	ldp	x14, x15, [x11, #80]
	ldp	x16, x17, [x11, #96]
	ldp	x6, x7, [x11, #112]
	add	x12, x12, x4
	add	x13, x13, x4
	stp	x12, x13, [x0, #64]
	add	x12, x14, x4
	add	x13, x15, x4
	stp	x12, x13, [x0, #80]
	add	x12, x16, x4
	add	x13, x17, x4
	add	x14, x6, x4
	stp	x12, x13, [x0, #96]
	add	x12, x7, x4
	add	x11, x11, #128
	stp	x14, x12, [x0, #112]
	add	x0, x0, #128
	sub	x10, x10, #16
	cmp	x10, #15
	b.hi	LBB0_5
LBB0_6:                                 ;   in Loop: Header=BB0_3 Depth=1
	cmp	x10, #8
	b.hs	LBB0_10
; %bb.7:                                ;   in Loop: Header=BB0_3 Depth=1
	tbnz	w10, #2, LBB0_11
LBB0_8:                                 ;   in Loop: Header=BB0_3 Depth=1
	tbnz	w10, #1, LBB0_12
LBB0_9:                                 ;   in Loop: Header=BB0_3 Depth=1
	tbz	w10, #0, LBB0_2
	b	LBB0_13
LBB0_10:                                ;   in Loop: Header=BB0_3 Depth=1
	ldp	x12, x13, [x11]
	ldp	x14, x15, [x11, #16]
	ldp	x16, x17, [x11, #32]
	ldp	x6, x7, [x11, #48]
	add	x12, x12, x4
	add	x13, x13, x4
	stp	x12, x13, [x0]
	add	x12, x14, x4
	add	x13, x15, x4
	stp	x12, x13, [x0, #16]
	add	x12, x16, x4
	add	x13, x17, x4
	stp	x12, x13, [x0, #32]
	add	x12, x6, x4
	add	x13, x7, x4
	add	x11, x11, #64
	stp	x12, x13, [x0, #48]
	add	x0, x0, #64
	tbz	w10, #2, LBB0_8
LBB0_11:                                ;   in Loop: Header=BB0_3 Depth=1
	ldp	x12, x13, [x11]
	ldp	x14, x15, [x11, #16]
	add	x12, x12, x4
	add	x13, x13, x4
	stp	x12, x13, [x0]
	add	x12, x14, x4
	add	x13, x15, x4
	stp	x12, x13, [x0, #16]
	add	x11, x11, #32
	add	x0, x0, #32
	tbz	w10, #1, LBB0_9
LBB0_12:                                ;   in Loop: Header=BB0_3 Depth=1
	ldp	x12, x13, [x11], #16
	add	x12, x12, x4
	add	x13, x13, x4
	stp	x12, x13, [x0], #16
	tbz	w10, #0, LBB0_2
LBB0_13:                                ;   in Loop: Header=BB0_3 Depth=1
	ldr	x10, [x11]
	add	x10, x10, x4
	str	x10, [x0], #8
	b	LBB0_2
LBB0_14:
	ret
	.cfi_endproc

_chord_probe_step:                      ; @chord_probe_step
	.cfi_startproc
; %bb.0:
                                        ; kill: def $w3 killed $w3 def $x3
	ldr	x8, [x0]
	add	x8, x8, x1
	cmp	x8, x2
	csel	x9, xzr, x2, lo
	sub	x8, x8, x9
	str	x8, [x0]
	lsr	x9, x8, x3
	lsl	x10, x9, x3
	add	x8, x8, w4, uxtw
	sub	x8, x8, x10
	lsr	x8, x8, x3
	add	x8, x8, #1
	str	x9, [x6]
	str	x8, [x7]
	add	x9, x5, x9, lsl #3
	ldr	x8, [x9, x8, lsl #3]
	ldr	x9, [x9]
	sub	x0, x8, x9
	ret
	.cfi_endproc

_emit_ptr64_x64:                        ; @emit_ptr64_x64
	.cfi_startproc
; %bb.0:
	cmp	x1, #64
	b.lo	LBB3_2
LBB3_1:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #64]
	ldp	x10, x11, [x0, #80]
	ldp	x12, x13, [x0, #96]
	ldp	x14, x15, [x0, #112]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #128]
	ldp	x10, x11, [x0, #144]
	ldp	x12, x13, [x0, #160]
	ldp	x14, x15, [x0, #176]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #192]
	ldp	x10, x11, [x0, #208]
	ldp	x12, x13, [x0, #224]
	ldp	x14, x15, [x0, #240]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #256]
	ldp	x10, x11, [x0, #272]
	ldp	x12, x13, [x0, #288]
	ldp	x14, x15, [x0, #304]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #320]
	ldp	x10, x11, [x0, #336]
	ldp	x12, x13, [x0, #352]
	ldp	x14, x15, [x0, #368]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #384]
	ldp	x10, x11, [x0, #400]
	ldp	x12, x13, [x0, #416]
	ldp	x14, x15, [x0, #432]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	ldp	x8, x9, [x0, #448]
	ldp	x10, x11, [x0, #464]
	ldp	x12, x13, [x0, #480]
	ldp	x14, x15, [x0, #496]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	add	x0, x0, #512
	sub	x1, x1, #64
	cmp	x1, #63
	b.hi	LBB3_1
LBB3_2:
	cmp	x1, #8
	b.lo	LBB3_4
LBB3_3:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x9]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x10]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x11]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x12]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x13]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x14]
	; InlineAsm End
	; InlineAsm Start
	staddb	w2, [x15]
	; InlineAsm End
	add	x0, x0, #64
	sub	x1, x1, #8
	cmp	x1, #7
	b.hi	LBB3_3
LBB3_4:
	cbz	x1, LBB3_6
LBB3_5:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x8, [x0], #8
	; InlineAsm Start
	staddb	w2, [x8]
	; InlineAsm End
	subs	x1, x1, #1
	b.ne	LBB3_5
LBB3_6:
	ret
	.cfi_endproc

_tx_copy_prebias_x256:                  ; @tx_copy_prebias_x256
	.cfi_startproc
; %bb.0:
	cbz	x2, LBB1_2
LBB1_1:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x1]
	ldp	x10, x11, [x1, #16]
	ldp	x12, x13, [x1, #32]
	ldp	x14, x15, [x1, #48]
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x3, x3, [x1]
	stp	x3, x3, [x1, #16]
	stp	x3, x3, [x1, #32]
	stp	x3, x3, [x1, #48]
	ldp	x8, x9, [x1, #64]
	ldp	x10, x11, [x1, #80]
	ldp	x12, x13, [x1, #96]
	ldp	x14, x15, [x1, #112]
	stp	x8, x9, [x0, #64]
	stp	x10, x11, [x0, #80]
	stp	x12, x13, [x0, #96]
	stp	x14, x15, [x0, #112]
	stp	x3, x3, [x1, #64]
	stp	x3, x3, [x1, #80]
	stp	x3, x3, [x1, #96]
	stp	x3, x3, [x1, #112]
	ldp	x8, x9, [x1, #128]
	ldp	x10, x11, [x1, #144]
	ldp	x12, x13, [x1, #160]
	ldp	x14, x15, [x1, #176]
	stp	x8, x9, [x0, #128]
	stp	x10, x11, [x0, #144]
	stp	x12, x13, [x0, #160]
	stp	x14, x15, [x0, #176]
	stp	x3, x3, [x1, #128]
	stp	x3, x3, [x1, #144]
	stp	x3, x3, [x1, #160]
	stp	x3, x3, [x1, #176]
	ldp	x8, x9, [x1, #192]
	ldp	x10, x11, [x1, #208]
	ldp	x12, x13, [x1, #224]
	ldp	x14, x15, [x1, #240]
	stp	x8, x9, [x0, #192]
	stp	x10, x11, [x0, #208]
	stp	x12, x13, [x0, #224]
	stp	x14, x15, [x0, #240]
	stp	x3, x3, [x1, #192]
	stp	x3, x3, [x1, #208]
	stp	x3, x3, [x1, #224]
	add	x0, x0, #256
	stp	x3, x3, [x1, #240]
	add	x1, x1, #256
	subs	x2, x2, #1
	b.ne	LBB1_1
LBB1_2:
	ret
	.cfi_endproc

_col_prebias_x256:                      ; @col_prebias_x256
	.cfi_startproc
; %bb.0:
	cbz	x1, LBB2_2
LBB2_1:                                 ; =>This Inner Loop Header: Depth=1
	stp	x2, x2, [x0]
	stp	x2, x2, [x0, #16]
	stp	x2, x2, [x0, #32]
	stp	x2, x2, [x0, #48]
	stp	x2, x2, [x0, #64]
	stp	x2, x2, [x0, #80]
	stp	x2, x2, [x0, #96]
	stp	x2, x2, [x0, #112]
	stp	x2, x2, [x0, #128]
	stp	x2, x2, [x0, #144]
	stp	x2, x2, [x0, #160]
	stp	x2, x2, [x0, #176]
	stp	x2, x2, [x0, #192]
	stp	x2, x2, [x0, #208]
	stp	x2, x2, [x0, #224]
	stp	x2, x2, [x0, #240]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB2_1
LBB2_2:
	ret
	.cfi_endproc

_merge_scan_swar_128_exact:             ; @merge_scan_swar_128_exact
	.cfi_startproc
; %bb.0:
	cbz	x2, LBB2_12
; %bb.1:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	mov	x9, #0                          ; =0x0
LBB2_2:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB2_8 Depth 2
                                        ;       Child Loop BB2_10 Depth 3
	ldp	x8, x10, [x0]
	ldp	x11, x12, [x1]
	add	x8, x11, x8
	ldp	x11, x13, [x0, #16]
	ldp	x14, x15, [x1, #16]
	add	x11, x14, x11
	ldp	x14, x16, [x0, #32]
	ldp	x17, x5, [x1, #32]
	add	x14, x17, x14
	ldp	x17, x6, [x0, #48]
	ldp	x7, x19, [x1, #48]
	add	x6, x19, x6
	ldp	x19, x20, [x0, #80]
	ldp	x21, x22, [x1, #80]
	add	x20, x22, x20
	add	x10, x12, x10
	orr	x8, x10, x8
	add	x10, x15, x13
	orr	x10, x11, x10
	orr	x8, x8, x10
	ldp	x10, x11, [x0, #64]
	add	x12, x5, x16
	orr	x12, x14, x12
	add	x13, x7, x17
	orr	x12, x12, x13
	orr	x8, x8, x12
	ldp	x12, x13, [x1, #64]
	add	x10, x12, x10
	orr	x10, x6, x10
	add	x11, x13, x11
	orr	x10, x10, x11
	ldp	x11, x12, [x0, #96]
	add	x13, x21, x19
	orr	x10, x10, x13
	orr	x8, x8, x10
	ldp	x10, x13, [x1, #96]
	add	x10, x10, x11
	orr	x10, x20, x10
	add	x11, x13, x12
	orr	x10, x10, x11
	ldp	x11, x12, [x0, #112]
	ldp	x13, x14, [x1, #112]
	add	x11, x13, x11
	orr	x10, x10, x11
	add	x11, x14, x12
	orr	x10, x10, x11
	orr	x8, x8, x10
	tst	x8, #0x8080808080808080
	b.ne	LBB2_5
; %bb.3:                                ;   in Loop: Header=BB2_2 Depth=1
	mov	x8, x9
LBB2_4:                                 ;   in Loop: Header=BB2_2 Depth=1
	add	x0, x0, #128
	stp	xzr, xzr, [x1, #112]
	stp	xzr, xzr, [x1, #96]
	stp	xzr, xzr, [x1, #80]
	stp	xzr, xzr, [x1, #64]
	stp	xzr, xzr, [x1, #48]
	stp	xzr, xzr, [x1, #32]
	stp	xzr, xzr, [x1, #16]
	stp	xzr, xzr, [x1], #128
	add	x3, x3, #128
	mov	x9, x8
	sub	x2, x2, #1
	cbnz	x2, LBB2_2
	b	LBB2_11
LBB2_5:                                 ;   in Loop: Header=BB2_2 Depth=1
	mov	x10, #0                         ; =0x0
	b	LBB2_8
LBB2_6:                                 ;   in Loop: Header=BB2_8 Depth=2
	mov	x8, x9
LBB2_7:                                 ;   in Loop: Header=BB2_8 Depth=2
	add	x10, x10, #1
	mov	x9, x8
	cmp	x10, #16
	b.eq	LBB2_4
LBB2_8:                                 ;   Parent Loop BB2_2 Depth=1
                                        ; =>  This Loop Header: Depth=2
                                        ;       Child Loop BB2_10 Depth 3
	ldr	x8, [x0, x10, lsl #3]
	ldr	x11, [x1, x10, lsl #3]
	add	x8, x11, x8
	ands	x11, x8, #0x8080808080808080
	b.eq	LBB2_6
; %bb.9:                                ;   in Loop: Header=BB2_8 Depth=2
	add	x12, x3, x10, lsl #3
LBB2_10:                                ;   Parent Loop BB2_2 Depth=1
                                        ;     Parent Loop BB2_8 Depth=2
                                        ; =>    This Inner Loop Header: Depth=3
	rbit	x8, x11
	clz	x8, x8
	add	x13, x12, x8, lsr #3
	add	x8, x9, #1
	str	x13, [x4, x9, lsl #3]
	sub	x13, x11, #1
	mov	x9, x8
	ands	x11, x13, x11
	b.ne	LBB2_10
	b	LBB2_7
LBB2_11:
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	mov	x0, x8
	ret
LBB2_12:
	mov	x0, #0                          ; =0x0
	ret
	.cfi_endproc
```

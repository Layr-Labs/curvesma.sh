# Compiled ARMv8.1-A Scalar Kernels (Part 1: Stage 1a Sieve & Banked Slab Scan, Kernels 1–6)

All 16 kernels across `KERNELS.md` (Kernels 1–6: **`615` actual scalar instructions**, `73` disclosed `staddb`) and `KERNELS_STAGE1B_LA.md` (Kernels 7–16: **987 actual scalar instructions**, `0` `staddb`) are compiled verbatim from a single standalone C translation unit with Apple `clang` targeting strict scalar ARMv8.1-A general-purpose registers:

```text
clang -O3 -march=armv8.1-a -mgeneral-regs-only -fno-vectorize -fno-slp-vectorize -fomit-frame-pointer -ffixed-x18 -S
```

**Global ISA & Register Hygiene Audit (All 16 Functions Combined = `1,602` Executable Scalar Instructions):**
- SIMD/FP/SVE/crypto/CRC registers or mnemonics (`v0-v31`, `q/d/s/z/p`, `ld1/st1`, `pmull`, `aes`, `sha`, `crc`): **`0`** (`simd = 0`).
- Hardware divisions (`udiv`/`sdiv`) or out-of-line calls (`bl`/`blr`): **`0`**.
- **In-loop stack spills ANYWHERE across all 16 functions (both fast and slow paths):** **`0`** (every `[sp` access in the entire 16-function binary is strictly function-prologue/epilogue callee-save `stp`/`ldp` or prologue ABI stack-argument load prior to loop entry).
- ARMv8.1-A LSE atomic byte adds (`staddb`): **`73`** total — ALL inside `_emit_ptr64_x64` (`LBB3_1`/`LBB3_3`/`LBB3_5`) for lock-free banked slab log accumulation, explicitly disclosed and charged at `1.0 Ir` each (`100/64 = 1.5625 Ir/hit` in `LBB3_1`).

## Part 1 Dynamic Instruction Count Summary (Kernels 1–6)

| # | Symbol | Role | Static Insts | In-Loop Spills | Compiled Dynamic `Ir` | Charged in Ledger |
|---:|---|---|---:|---:|---|---|
| 1 | `_tx_pack_u64_window` | Tier-1 (`p <= 2^{54}`) contiguous `ptr64` slice pack (`x64`) | 168 | **0** | `7` col + `10` non-empty + `69/64 = 1.078125 Ir/hit` + tail (`LBB0_6`) | Exact literal Poisson (`10 + 69*q64 + tail`) |
| 2 | `_w4_packed_tx_window` | Tier-2 (`2^{54} < p <= 2^{57}`) `W_B=4` 4-slice header pack | 193 | **0** | `5` empty bkt + `12` non-empty B4; `k=1: 6.5, k=2: 9.5, k>=3: 13.5+13*q8+tail` | Exact literal Poisson (`5/12`, `6.5/9.5/13.5+`) |
| 3 | `_demux_w16_window` | Tier-3 (`p > 2^{57}`) `W_B=16` zero-AND (`SLAB_BASE=2^{56}`) demux | 48 | **0** | `5` empty bkt, `4` non-empty, `26/4 = 6.50 Ir/hit` (`LBB2_5`), `8.00 Ir` tail (`LBB2_7`) | Exact literal Poisson over ALL `n_b*16` cols |
| 4 | `_emit_ptr64_x64` | Banked `ptr64` slab log scatter (`73` `staddb`, `x64` unrolled) | 126 | **0** | **`1.5625 Ir/hit`** (`100/64` in `LBB3_1`, zero `add`) | **`1.5625 Ir/hit`** (256 KiB batched ring) |
| 5 | `_col_prebias_x128` | 128-byte column pre-bias (`stp` x8) | 16 | **0** | `0.0859375 Ir/cell/pass` (`11/128`) | **`0.171875 Ir/cell`** (2 passes) |
| 6 | `_col_scan_swar_128_nospill` | 128-byte SWAR MSB scan + `rbit`/`clz` slow path (**0 spills**) | 64 | **0** | Fast: `30/128 = 0.234375 Ir/cell`; Slow: `140<=142 Ir` + `8 Ir` `f`-probe = `150 Ir/surv` | **`0.417621 Ir/cell`** (`+0.002` pad) |

### Exact Verification of Zero-Spill `_col_scan_swar_128_nospill` (`64` Scalar Instructions)
- Prologue/epilogue saves `x19–x22` once outside the loop (`stp/ldp` 4 insns total). **Zero `[sp` accesses inside `LBB5_2`, `LBB5_5`, `LBB5_8`, or `LBB5_10`.**
- **Fast path (`LBB5_2 -> LBB5_3 -> LBB5_4`):** `8 ldp` + `15 orr` + `tst` + `b.ne` + `mov` + `add` + `mov` + `subs` + `b.ne` = **`30` instructions per 128 bytes = `0.234375 Ir/cell`**.
- **Slow path (`LBB5_5 -> LBB5_8 -> LBB5_10`):** `3` setup (`mov/sub/b`) + `15 * 8` empty-word checks (`ldr/ands/b.eq/mov/add/mov/cmp/b.eq`) + `15` hit-word (`ldr/ands/b.eq/add` + `rbit/clz/add/add/str/sub/mov/ands/b.ne` + `b`) + `2` tail = **`140 Ir <= 142 Ir`** per 128B block with 1 survivor, plus `8 Ir` `f`-side probe = **`150.0 Ir/survivor`**. At `log2(s_g) = -13.966476`, `slab_live = 0.171875 + 0.234375 + 150.0 * 2^{-13.966476} = 0.415621 Ir/cell`, charged **`0.417621 Ir/cell`** (`+0.002` pad).

## Verbatim Compiled ARMv8.1-A Assembly (Kernels 1–6, Zero Ellipses)

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

_w4_packed_tx_window:                   ; @w4_packed_tx_window
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB1_59
; %bb.1:
	mov	x8, #0                          ; =0x0
	b	LBB1_6
LBB1_2:                                 ;   in Loop: Header=BB1_6 Depth=1
	ldr	x11, [x10], #8
	str	x11, [x12], #8
LBB1_3:                                 ;   in Loop: Header=BB1_6 Depth=1
	str	x12, [x9, #24]
LBB1_4:                                 ;   in Loop: Header=BB1_6 Depth=1
	mov	x1, x10
LBB1_5:                                 ;   in Loop: Header=BB1_6 Depth=1
	add	x8, x8, #1
	cmp	x8, x3
	b.eq	LBB1_59
LBB1_6:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB1_31 Depth 2
                                        ;     Child Loop BB1_36 Depth 2
                                        ;     Child Loop BB1_41 Depth 2
                                        ;     Child Loop BB1_46 Depth 2
	ldr	w9, [x2, x8, lsl #2]
	cbz	w9, LBB1_5
; %bb.7:                                ;   in Loop: Header=BB1_6 Depth=1
	mov	x10, x1
	ldr	x11, [x10], #8
	add	x9, x0, x8, lsl #5
	ands	x13, x11, #0xffff
	b.eq	LBB1_12
; %bb.8:                                ;   in Loop: Header=BB1_6 Depth=1
	ldr	x12, [x9]
	and	w14, w11, #0xffff
	cmp	w14, #1
	b.ne	LBB1_24
; %bb.9:                                ;   in Loop: Header=BB1_6 Depth=1
	add	x10, x1, #16
	ldr	x13, [x1, #8]
LBB1_10:                                ;   in Loop: Header=BB1_6 Depth=1
	str	x13, [x12], #8
LBB1_11:                                ;   in Loop: Header=BB1_6 Depth=1
	str	x12, [x9]
LBB1_12:                                ;   in Loop: Header=BB1_6 Depth=1
	ubfx	x13, x11, #16, #16
	cbz	x13, LBB1_16
; %bb.13:                               ;   in Loop: Header=BB1_6 Depth=1
	ldr	x12, [x9, #8]
	cmp	w13, #1
	b.ne	LBB1_26
LBB1_14:                                ;   in Loop: Header=BB1_6 Depth=1
	ldr	x13, [x10], #8
	str	x13, [x12], #8
LBB1_15:                                ;   in Loop: Header=BB1_6 Depth=1
	str	x12, [x9, #8]
LBB1_16:                                ;   in Loop: Header=BB1_6 Depth=1
	ubfx	x13, x11, #32, #16
	cbz	x13, LBB1_20
; %bb.17:                               ;   in Loop: Header=BB1_6 Depth=1
	ldr	x12, [x9, #16]
	cmp	w13, #1
	b.ne	LBB1_28
LBB1_18:                                ;   in Loop: Header=BB1_6 Depth=1
	ldr	x13, [x10], #8
	str	x13, [x12], #8
LBB1_19:                                ;   in Loop: Header=BB1_6 Depth=1
	str	x12, [x9, #16]
LBB1_20:                                ;   in Loop: Header=BB1_6 Depth=1
	lsr	x13, x11, #48
	cbz	x13, LBB1_4
; %bb.21:                               ;   in Loop: Header=BB1_6 Depth=1
	ldr	x12, [x9, #24]
	cmp	w13, #1
	b.eq	LBB1_2
; %bb.22:                               ;   in Loop: Header=BB1_6 Depth=1
	cmp	w13, #2
	b.ne	LBB1_45
; %bb.23:                               ;   in Loop: Header=BB1_6 Depth=1
	ldp	x11, x13, [x10], #16
	stp	x11, x13, [x12], #16
	b	LBB1_3
LBB1_24:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	w14, #2
	b.ne	LBB1_30
; %bb.25:                               ;   in Loop: Header=BB1_6 Depth=1
	ldp	x10, x13, [x1, #8]
	stp	x10, x13, [x12], #16
	add	x10, x1, #24
	b	LBB1_11
LBB1_26:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	w13, #2
	b.ne	LBB1_35
; %bb.27:                               ;   in Loop: Header=BB1_6 Depth=1
	ldp	x13, x14, [x10], #16
	stp	x13, x14, [x12], #16
	b	LBB1_15
LBB1_28:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	w13, #2
	b.ne	LBB1_40
; %bb.29:                               ;   in Loop: Header=BB1_6 Depth=1
	ldp	x13, x14, [x10], #16
	stp	x13, x14, [x12], #16
	b	LBB1_19
LBB1_30:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #8
	b.lo	LBB1_32
LBB1_31:                                ;   Parent Loop BB1_6 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x14, x15, [x10]
	ldp	x16, x17, [x10, #16]
	ldp	x1, x4, [x10, #32]
	ldp	x5, x6, [x10, #48]
	stp	x14, x15, [x12]
	stp	x16, x17, [x12, #16]
	stp	x1, x4, [x12, #32]
	stp	x5, x6, [x12, #48]
	add	x12, x12, #64
	add	x10, x10, #64
	sub	x13, x13, #8
	cmp	x13, #7
	b.hi	LBB1_31
LBB1_32:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #4
	b.hs	LBB1_50
; %bb.33:                               ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #1, LBB1_51
LBB1_34:                                ;   in Loop: Header=BB1_6 Depth=1
	tbz	w13, #0, LBB1_11
	b	LBB1_52
LBB1_35:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #8
	b.lo	LBB1_37
LBB1_36:                                ;   Parent Loop BB1_6 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x14, x15, [x10]
	ldp	x16, x17, [x10, #16]
	ldp	x1, x4, [x10, #32]
	ldp	x5, x6, [x10, #48]
	stp	x14, x15, [x12]
	stp	x16, x17, [x12, #16]
	stp	x1, x4, [x12, #32]
	stp	x5, x6, [x12, #48]
	add	x12, x12, #64
	add	x10, x10, #64
	sub	x13, x13, #8
	cmp	x13, #7
	b.hi	LBB1_36
LBB1_37:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #4
	b.hs	LBB1_53
; %bb.38:                               ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #1, LBB1_54
LBB1_39:                                ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #0, LBB1_14
	b	LBB1_15
LBB1_40:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #8
	b.lo	LBB1_42
LBB1_41:                                ;   Parent Loop BB1_6 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x14, x15, [x10]
	ldp	x16, x17, [x10, #16]
	ldp	x1, x4, [x10, #32]
	ldp	x5, x6, [x10, #48]
	stp	x14, x15, [x12]
	stp	x16, x17, [x12, #16]
	stp	x1, x4, [x12, #32]
	stp	x5, x6, [x12, #48]
	add	x12, x12, #64
	add	x10, x10, #64
	sub	x13, x13, #8
	cmp	x13, #7
	b.hi	LBB1_41
LBB1_42:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #4
	b.hs	LBB1_55
; %bb.43:                               ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #1, LBB1_56
LBB1_44:                                ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #0, LBB1_18
	b	LBB1_19
LBB1_45:                                ;   in Loop: Header=BB1_6 Depth=1
	lsr	x11, x11, #51
	cbz	x11, LBB1_47
LBB1_46:                                ;   Parent Loop BB1_6 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x11, x14, [x10]
	ldp	x15, x16, [x10, #16]
	ldp	x17, x1, [x10, #32]
	ldp	x4, x5, [x10, #48]
	stp	x11, x14, [x12]
	stp	x15, x16, [x12, #16]
	stp	x17, x1, [x12, #32]
	stp	x4, x5, [x12, #48]
	add	x12, x12, #64
	add	x10, x10, #64
	sub	x13, x13, #8
	cmp	x13, #7
	b.hi	LBB1_46
LBB1_47:                                ;   in Loop: Header=BB1_6 Depth=1
	cmp	x13, #4
	b.hs	LBB1_57
; %bb.48:                               ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #1, LBB1_58
LBB1_49:                                ;   in Loop: Header=BB1_6 Depth=1
	tbnz	w13, #0, LBB1_2
	b	LBB1_3
LBB1_50:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x14, x15, [x10]
	ldp	x16, x17, [x10, #16]
	stp	x14, x15, [x12]
	stp	x16, x17, [x12, #16]
	add	x12, x12, #32
	add	x10, x10, #32
	tbz	w13, #1, LBB1_34
LBB1_51:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x14, x15, [x10], #16
	stp	x14, x15, [x12], #16
	tbz	w13, #0, LBB1_11
LBB1_52:                                ;   in Loop: Header=BB1_6 Depth=1
	ldr	x13, [x10], #8
	b	LBB1_10
LBB1_53:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x14, x15, [x10]
	ldp	x16, x17, [x10, #16]
	stp	x14, x15, [x12]
	stp	x16, x17, [x12, #16]
	add	x12, x12, #32
	add	x10, x10, #32
	tbz	w13, #1, LBB1_39
LBB1_54:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x14, x15, [x10], #16
	stp	x14, x15, [x12], #16
	tbnz	w13, #0, LBB1_14
	b	LBB1_15
LBB1_55:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x14, x15, [x10]
	ldp	x16, x17, [x10, #16]
	stp	x14, x15, [x12]
	stp	x16, x17, [x12, #16]
	add	x12, x12, #32
	add	x10, x10, #32
	tbz	w13, #1, LBB1_44
LBB1_56:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x14, x15, [x10], #16
	stp	x14, x15, [x12], #16
	tbnz	w13, #0, LBB1_18
	b	LBB1_19
LBB1_57:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x11, x14, [x10]
	ldp	x15, x16, [x10, #16]
	stp	x11, x14, [x12]
	stp	x15, x16, [x12, #16]
	add	x12, x12, #32
	add	x10, x10, #32
	tbz	w13, #1, LBB1_49
LBB1_58:                                ;   in Loop: Header=BB1_6 Depth=1
	ldp	x11, x14, [x10], #16
	stp	x11, x14, [x12], #16
	tbnz	w13, #0, LBB1_2
	b	LBB1_3
LBB1_59:
	ret
	.cfi_endproc

_demux_w16_window:                      ; @demux_w16_window
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB2_8
; %bb.1:
	mov	x8, #0                          ; =0x0
	b	LBB2_3
LBB2_2:                                 ;   in Loop: Header=BB2_3 Depth=1
	add	x8, x8, #1
	cmp	x8, x3
	b.eq	LBB2_8
LBB2_3:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB2_5 Depth 2
                                        ;     Child Loop BB2_7 Depth 2
	ldr	w9, [x2, x8, lsl #2]
	cbz	w9, LBB2_2
; %bb.4:                                ;   in Loop: Header=BB2_3 Depth=1
	add	x10, x0, x8, lsl #7
	cmp	w9, #4
	b.lo	LBB2_7
LBB2_5:                                 ;   Parent Loop BB2_3 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldp	x11, x12, [x1]
	ldp	x13, x14, [x1, #16]
	lsr	x15, x11, #57
	and	x11, x11, #0xfffffffffffffff
	ldr	x16, [x10, x15]
	str	x11, [x16], #8
	str	x16, [x10, x15]
	lsr	x11, x12, #57
	and	x12, x12, #0xfffffffffffffff
	ldr	x15, [x10, x11]
	str	x12, [x15], #8
	str	x15, [x10, x11]
	lsr	x11, x13, #57
	and	x12, x13, #0xfffffffffffffff
	ldr	x13, [x10, x11]
	str	x12, [x13], #8
	str	x13, [x10, x11]
	lsr	x11, x14, #57
	and	x12, x14, #0xfffffffffffffff
	ldr	x13, [x10, x11]
	str	x12, [x13], #8
	str	x13, [x10, x11]
	add	x1, x1, #32
	sub	x9, x9, #4
	cmp	x9, #3
	b.hi	LBB2_5
; %bb.6:                                ;   in Loop: Header=BB2_3 Depth=1
	cbz	x9, LBB2_2
LBB2_7:                                 ;   Parent Loop BB2_3 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	ldr	x11, [x1], #8
	lsr	x12, x11, #57
	ldr	x13, [x10, x12]
	and	x11, x11, #0xfffffffffffffff
	str	x11, [x13], #8
	str	x13, [x10, x12]
	subs	x9, x9, #1
	b.ne	LBB2_7
	b	LBB2_2
LBB2_8:
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

_col_prebias_x128:                      ; @col_prebias_x128
	.cfi_startproc
; %bb.0:
	cbz	x1, LBB4_3
; %bb.1:
	mov	w8, w2
	mov	x9, #72340172838076673          ; =0x101010101010101
	mul	x8, x8, x9
LBB4_2:                                 ; =>This Inner Loop Header: Depth=1
	stp	x8, x8, [x0]
	stp	x8, x8, [x0, #16]
	stp	x8, x8, [x0, #32]
	stp	x8, x8, [x0, #48]
	stp	x8, x8, [x0, #64]
	stp	x8, x8, [x0, #80]
	stp	x8, x8, [x0, #96]
	stp	x8, x8, [x0, #112]
	add	x0, x0, #128
	subs	x1, x1, #1
	b.ne	LBB4_2
LBB4_3:
	ret
	.cfi_endproc

_col_scan_swar_128_nospill:             ; @col_scan_swar_128_nospill
	.cfi_startproc
; %bb.0:
	cbz	x1, LBB5_12
; %bb.1:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	mov	x10, #0                         ; =0x0
	mov	x9, x0
LBB5_2:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB5_8 Depth 2
                                        ;       Child Loop BB5_10 Depth 3
	ldp	x8, x11, [x9]
	ldp	x12, x13, [x9, #16]
	ldp	x14, x15, [x9, #32]
	ldp	x16, x17, [x9, #48]
	ldp	x3, x4, [x9, #64]
	ldp	x5, x6, [x9, #80]
	ldp	x7, x19, [x9, #96]
	ldp	x20, x21, [x9, #112]
	orr	x8, x11, x8
	orr	x11, x12, x13
	orr	x8, x8, x11
	orr	x11, x14, x15
	orr	x11, x11, x16
	orr	x8, x8, x11
	orr	x11, x17, x3
	orr	x11, x11, x4
	orr	x11, x11, x5
	orr	x8, x8, x11
	orr	x11, x6, x7
	orr	x11, x11, x19
	orr	x11, x11, x20
	orr	x11, x11, x21
	orr	x8, x8, x11
	tst	x8, #0x8080808080808080
	b.ne	LBB5_5
; %bb.3:                                ;   in Loop: Header=BB5_2 Depth=1
	mov	x8, x10
LBB5_4:                                 ;   in Loop: Header=BB5_2 Depth=1
	add	x9, x9, #128
	mov	x10, x8
	subs	x1, x1, #1
	b.ne	LBB5_2
	b	LBB5_11
LBB5_5:                                 ;   in Loop: Header=BB5_2 Depth=1
	mov	x11, #0                         ; =0x0
	sub	x12, x9, x0
	b	LBB5_8
LBB5_6:                                 ;   in Loop: Header=BB5_8 Depth=2
	mov	x8, x10
LBB5_7:                                 ;   in Loop: Header=BB5_8 Depth=2
	add	x11, x11, #1
	mov	x10, x8
	cmp	x11, #16
	b.eq	LBB5_4
LBB5_8:                                 ;   Parent Loop BB5_2 Depth=1
                                        ; =>  This Loop Header: Depth=2
                                        ;       Child Loop BB5_10 Depth 3
	ldr	x8, [x9, x11, lsl #3]
	ands	x13, x8, #0x8080808080808080
	b.eq	LBB5_6
; %bb.9:                                ;   in Loop: Header=BB5_8 Depth=2
	add	x14, x12, x11, lsl #3
LBB5_10:                                ;   Parent Loop BB5_2 Depth=1
                                        ;     Parent Loop BB5_8 Depth=2
                                        ; =>    This Inner Loop Header: Depth=3
	rbit	x8, x13
	clz	x8, x8
	add	x15, x14, x8, lsr #3
	add	x8, x10, #1
	str	x15, [x2, x10, lsl #3]
	sub	x15, x13, #1
	mov	x10, x8
	ands	x13, x15, x13
	b.ne	LBB5_10
	b	LBB5_7
LBB5_11:
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	mov	x0, x8
	ret
LBB5_12:
	mov	x0, #0                          ; =0x0
	ret
	.cfi_endproc
```

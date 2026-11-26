# Stage 1a Guard4 x16 Sieve, Replay-Build, and Sequential Single-Slab Active-Only Scan Kernels (`KERNELS_SIEVE.md`)

All assembly listings in this file were compiled with Apple Clang (`arm64`) under strict scalar general-purpose-register constraints:

```text
-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize
```

## 1. ISA Compliance, `STADDB` Disclosure & Guard-Masked Correctness

- **SIMD / FP / SVE / Platform-Register Audit:** `simd = 0`. Zero NEON/FP/SVE vector registers (`v0..v31`, `q0..q31`, `d0..d31`, `s0..s31`, `h0..h31`, `z0..z31`, `p0..p15`), zero floating-point instructions, zero `dc zva`, and zero `x18` accesses appear anywhere in this file.
- **Explicit ARMv8.1-A Scalar LSE Disclosure:** The Guard4 x16 hit-emission kernel `_emit_u64_x16_direct` (invoked by `_sieve_window_guard4`) uses the single scalar ARMv8.1-A Large System Extensions (LSE) byte atomic-add instruction **`STADDB Ws, [Xn]`** (a 64-bit general-purpose-register-only memory add with zero SIMD state, counted as `1` instruction each). The replay builder and sequential single-slab active-only scan kernels use pure ARMv8.0 scalar instructions only (`0` `staddb`).
- **Sequential Single-Slab Active-Only Masked Execution (Closing Guard-Leak & Dual-Slab Overhead):**
  - Each `a_id` row occupies 256 bytes: 4-byte left guard `[0..3]`, $L \le 248$ active chord columns `[4..4+L-1]`, and $\ge 4$-byte right guard `[4+L..255]`.
  - ONE single `14.9616 EiB` flat slab is reused sequentially:
    1. **Pass 1 ($g$-side):** `_slab_prebias_active_sl` writes only the $K_{\text{even}} = 2\lceil\lceil(L+4)/16\rceil/2\rceil$ active 16-byte pairs per row (`K_even_avg = 8.40` at `L_win_eff = 112.49`), followed by `_sieve_window_guard4` ($g$-band), followed by `_slab_scan_g_active_sl` which **explicitly masks BOTH the 4-byte left guard `[0..3]` AND the right guard `[4+L..255]`** before testing sign bits, writing the `2^{-14.41}` $g$-survivor offsets into a compact `<= 1.4 PiB` list.
    2. **Pass 2 ($f$-side, SAME slab):** `_slab_prebias_active_sl` overwrites the active region with the $f$-bias, `_sieve_window_guard4` emits $f$-band hits, and `_slab_filter_f_at_gsurv` probes ONLY at the exact $g$-survivor coordinates (`5 Ir` per non-passing $g$-survivor $\implies 5 \cdot 2^{-14.4064} = 0.00023$ Ir/active cell).
- **Live Differential Verification:** `experiment.py` verifies Guard4 x16 addressing (`W_B = 4`, `L <= 248`, forced prime wrap-around and 256-boundary splits) AND active-only masked scanning with **`0xFF` poison written into EVERY left `[0..3]` and right `[4+L..255]` guard byte** across randomized chord lengths `L in [1, 248]`: **0 false positives, 0 missed survivors, 0 outer-guard corruptions**.

## 2. Exact Basic-Block Instruction Counts

| Function | Basic Block(s) | Exact Compiled Count | Charged in Ledger (`proof.md`) | Notes |
|---|---|---:|---:|---|
| `_emit_u64_x16_direct` | `LBB0_2` (x16 unrolled hits) | `43 Ir / 16 = 2.6875 Ir/hit` | **`2.6875 Ir/hit`** | `8 ldp + 17 add + 16 staddb + cmp + b.lo` |
| `_emit_u64_x16_direct` | `LBB0_4` (x4 hits) | `13 Ir / 4 = 3.2500 Ir/hit` | **`3.2500 Ir/hit`** | `2 ldp + 5 add + 4 staddb + cmp + b.lo` |
| `_emit_u64_x16_direct` | `LBB0_6` (x1 tail) | `5 Ir / 1 = 5.0000 Ir/hit` | **`5.0000 Ir/hit`** | `ldr + add + staddb + cmp + b.lo` |
| `_emit_u64_x16_direct` | Prologue/epilogue/loop-dispatch (`%bb.0..LBB0_7`) | `16–18 Ir/call` (`17` nominal) | Included in probe | Included in `<= 56` / `<= 84` Ir probe bounds |
| `_sieve_window_guard4` | Non-split probe + call (`LBB1_4..LBB1_2..LBB1_3` + callee) | `28 Ir` (empty) / `55 Ir` (hit) | **`56.0 Ir/prime`** | Includes full `emit_u64_x16_direct` call overhead |
| `_sieve_window_guard4` | Split-256 probe + calls (`LBB1_4..LBB1_11..LBB1_3` + callees) | `32 Ir` (empty) / `81–85 Ir` (Poisson-avg `< 84.0`) | **`84.0 Ir/prime`** | Includes both slice calls' prologue/epilogue |
| `_replay_pass1_count_guard4` | `LBB0_1` (counting pass) | `11 Ir/entry` | `11.0` (+Pass 2 + pad = **`27.5 Ir/entry`**) | `Pass 1 (11) + Pass 2 (15) = 26 Ir` |
| `_replay_pass2_scatter_guard4` | `LBB1_1` (scatter pass) | `15 Ir/entry` | `15.0` (+`1.5` startup = **`27.5 Ir/entry`**) | `+1.5 Ir/entry` per-$a_0$-line startup pad |
| `_replay_prefix_sum_guard4` | `LBB2_2` (prefix sum) | `5 Ir/bucket` | **`6.5 Ir/bucket`** | `+1.5 Ir/bucket` zeroing cushion |
| `_slab_prebias_active_sl` | `LBB0_4/8/11/14/16/18/20/22` (`K_even` `stp` + 3) | `K_even + 3` Ir/row (`0.2027` for 2 passes) | Part of **`0.5500`** | `2 * (8.40 + 3) / 112.49 = 0.2027 Ir/active` |
| `_slab_scan_g_active_sl` | Masked SWAR fast path + `_extract_row_survivors` | `37.0 / 112.49 = 0.3289 Ir/active` | Part of **`0.5500`** | Masks left `[0..3]` & right `[4+L..255]` |
| `_slab_filter_f_at_gsurv` | `LBB3_3` + `LBB3_2` (non-surv `5 Ir`; surv `8 Ir`) | `5 * 2^-14.4064 = 0.00023 Ir/active` | **`0.5500 Ir/active`** total | Exact sum `0.5318 <= 0.5500 Ir/active cell` |

---

## 3. Full Assembly Listing: Guard4 x16 Chord-Window Sieve (`kernels_guard4_x16_v81.s`)

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_emit_u64_x16_direct            ; -- Begin function emit_u64_x16_direct
	.p2align	2
_emit_u64_x16_direct:                   ; @emit_u64_x16_direct
	.cfi_startproc
; %bb.0:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	and	x9, x2, #0xfffffffffffffff0
	mov	x8, x1
	cmp	x9, #1
	b.lt	LBB0_3
; %bb.1:
	add	x9, x1, x9, lsl #3
	mov	x8, x1
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x10, x11, [x8]
	ldp	x12, x13, [x8, #16]
	ldp	x14, x15, [x8, #32]
	ldp	x16, x17, [x8, #48]
	ldp	x4, x5, [x8, #64]
	ldp	x6, x7, [x8, #80]
	ldp	x19, x20, [x8, #96]
	ldp	x21, x22, [x8, #112]
	add	x8, x8, #128
	add	x10, x0, x10
	add	x11, x0, x11
	add	x12, x0, x12
	add	x13, x0, x13
	add	x14, x0, x14
	add	x15, x0, x15
	add	x16, x0, x16
	add	x17, x0, x17
	; InlineAsm Start
	staddb	w3, [x10]
	staddb	w3, [x11]
	staddb	w3, [x12]
	staddb	w3, [x13]
	staddb	w3, [x14]
	staddb	w3, [x15]
	staddb	w3, [x16]
	staddb	w3, [x17]

	; InlineAsm End
	add	x10, x0, x4
	add	x11, x0, x5
	add	x12, x0, x6
	add	x13, x0, x7
	add	x14, x0, x19
	add	x15, x0, x20
	add	x16, x0, x21
	add	x17, x0, x22
	; InlineAsm Start
	staddb	w3, [x10]
	staddb	w3, [x11]
	staddb	w3, [x12]
	staddb	w3, [x13]
	staddb	w3, [x14]
	staddb	w3, [x15]
	staddb	w3, [x16]
	staddb	w3, [x17]

	; InlineAsm End
	cmp	x8, x9
	b.lo	LBB0_2
LBB0_3:
	and	x9, x2, #0x1ffffffffffffffc
	add	x9, x1, x9, lsl #3
	cmp	x8, x9
	b.hs	LBB0_5
LBB0_4:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x10, x11, [x8]
	ldp	x12, x13, [x8, #16]
	add	x8, x8, #32
	add	x10, x0, x10
	add	x11, x0, x11
	add	x12, x0, x12
	add	x13, x0, x13
	; InlineAsm Start
	staddb	w3, [x10]
	staddb	w3, [x11]
	staddb	w3, [x12]
	staddb	w3, [x13]

	; InlineAsm End
	cmp	x8, x9
	b.lo	LBB0_4
LBB0_5:
	add	x9, x1, x2, lsl #3
	cmp	x8, x9
	b.hs	LBB0_7
LBB0_6:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x10, [x8], #8
	add	x10, x0, x10
	; InlineAsm Start
	staddb	w3, [x10]
	; InlineAsm End
	cmp	x8, x9
	b.lo	LBB0_6
LBB0_7:
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_sieve_window_guard4            ; -- Begin function sieve_window_guard4
	.p2align	2
_sieve_window_guard4:                   ; @sieve_window_guard4
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #144
	stp	x28, x27, [sp, #48]             ; 16-byte Folded Spill
	stp	x26, x25, [sp, #64]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #80]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #96]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #112]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #128]            ; 16-byte Folded Spill
	add	x29, sp, #128
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	cbz	x2, LBB1_22
; %bb.1:
	mov	x19, x3
	mov	x20, x2
	mov	x22, x1
	mov	x21, #0                         ; =0x0
	add	x8, x0, #4
	stp	x8, x0, [sp, #32]               ; 16-byte Folded Spill
	b	LBB1_4
LBB1_2:                                 ;   in Loop: Header=BB1_4 Depth=1
	and	w3, w26, #0xff
	mov	x2, x24
	bl	_emit_u64_x16_direct
	add	x21, x24, x21
LBB1_3:                                 ;   in Loop: Header=BB1_4 Depth=1
	add	x22, x22, #48
	subs	x20, x20, #1
	b.eq	LBB1_23
LBB1_4:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x22]
	ldp	x24, x28, [x22, #16]
	add	x8, x9, x8
	cmp	x8, x24
	csel	x9, xzr, x24, lo
	sub	x25, x8, x9
	add	x10, x25, x19
	str	x25, [x22]
	ldr	x27, [x22, #32]
	ldrb	w26, [x22, #40]
	lsr	x8, x25, #8
	lsl	x9, x25, #1
	cmp	x10, x24
	b.hi	LBB1_12
; %bb.5:                                ;   in Loop: Header=BB1_4 Depth=1
	sub	x10, x10, #1
	lsr	x11, x10, #2
	and	x9, x9, #0xfffffffffffffff8
	ldr	x9, [x28, x9]
	and	x25, x25, #0xff
	cmp	x8, x10, lsr #8
	add	x23, x28, x11, lsl #3
	b.ne	LBB1_8
; %bb.6:                                ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x24, x8, x9
	b.ls	LBB1_3
; %bb.7:                                ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x25
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	b	LBB1_2
LBB1_8:                                 ;   in Loop: Header=BB1_4 Depth=1
	add	x8, x28, x8, lsl #9
	ldr	x28, [x8, #512]
	subs	x24, x28, x9
	b.ls	LBB1_10
; %bb.9:                                ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x25
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	and	w3, w26, #0xff
	mov	x2, x24
	bl	_emit_u64_x16_direct
	add	x21, x24, x21
LBB1_10:                                ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x24, x8, x28
	b.ls	LBB1_3
; %bb.11:                               ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x25
	add	x1, x27, x28, lsl #3
	add	x0, x8, #260
	b	LBB1_2
LBB1_12:                                ;   in Loop: Header=BB1_4 Depth=1
	sub	x10, x24, #1
	lsr	x11, x10, #2
	and	x9, x9, #0xfffffffffffffff8
	ldr	x9, [x28, x9]
	and	x12, x25, #0xff
	cmp	x8, x10, lsr #8
	add	x23, x28, x11, lsl #3
	b.ne	LBB1_15
; %bb.13:                               ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x23, x8, x9
	b.ls	LBB1_20
; %bb.14:                               ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x12
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	b	LBB1_19
LBB1_15:                                ;   in Loop: Header=BB1_4 Depth=1
	add	x8, x28, x8, lsl #9
	ldr	x10, [x8, #512]
	subs	x11, x10, x9
	str	x12, [sp, #24]                  ; 8-byte Folded Spill
	b.ls	LBB1_17
; %bb.16:                               ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x12
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	and	w3, w26, #0xff
	mov	x2, x11
	stp	x11, x10, [sp, #8]              ; 16-byte Folded Spill
	bl	_emit_u64_x16_direct
	ldp	x8, x10, [sp, #8]               ; 16-byte Folded Reload
	add	x21, x8, x21
LBB1_17:                                ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x23, x8, x10
	b.ls	LBB1_20
; %bb.18:                               ;   in Loop: Header=BB1_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	ldr	x9, [sp, #24]                   ; 8-byte Folded Reload
	sub	x8, x8, x9
	add	x1, x27, x10, lsl #3
	add	x0, x8, #260
LBB1_19:                                ;   in Loop: Header=BB1_4 Depth=1
	and	w3, w26, #0xff
	mov	x2, x23
	bl	_emit_u64_x16_direct
	add	x21, x23, x21
LBB1_20:                                ;   in Loop: Header=BB1_4 Depth=1
	sub	x8, x24, x25
	mvn	x9, x8
	add	x10, x9, x19
	ldr	x9, [x28]
	lsr	x10, x10, #2
	add	x10, x28, x10, lsl #3
	ldr	x10, [x10, #8]
	subs	x24, x10, x9
	b.ls	LBB1_3
; %bb.21:                               ;   in Loop: Header=BB1_4 Depth=1
	add	x1, x27, x9, lsl #3
	ldr	x9, [sp, #32]                   ; 8-byte Folded Reload
	add	x0, x9, x8
	b	LBB1_2
LBB1_22:
	mov	x21, #0                         ; =0x0
LBB1_23:
	mov	x0, x21
	ldp	x29, x30, [sp, #128]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #112]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #96]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #80]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #64]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #48]             ; 16-byte Folded Reload
	add	sp, sp, #144
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
```

---

## 4. Full Assembly Listing: 2-Pass Guard4 Replay Builder (`replay_build_guard4_v81.s`)

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_replay_pass1_count_guard4      ; -- Begin function replay_pass1_count_guard4
	.p2align	2
_replay_pass1_count_guard4:             ; @replay_pass1_count_guard4
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB0_2
LBB0_1:                                 ; =>This Inner Loop Header: Depth=1
	lsl	x8, x0, #1
	and	x8, x8, #0xfffffffffffffff8
	ldr	x9, [x4, x8]
	add	x9, x9, #1
	str	x9, [x4, x8]
	add	x8, x0, x1
	cmp	x8, x2
	csel	x9, xzr, x2, lo
	sub	x0, x8, x9
	subs	x3, x3, #1
	b.ne	LBB0_1
LBB0_2:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_replay_pass2_scatter_guard4    ; -- Begin function replay_pass2_scatter_guard4
	.p2align	2
_replay_pass2_scatter_guard4:           ; @replay_pass2_scatter_guard4
	.cfi_startproc
; %bb.0:
	cbz	x4, LBB1_2
LBB1_1:                                 ; =>This Inner Loop Header: Depth=1
	lsl	x8, x0, #1
	and	x8, x8, #0xfffffffffffffff8
	ldr	x9, [x5, x8]
	add	x10, x9, #1
	str	x10, [x5, x8]
	and	x8, x0, #0xff
	orr	x8, x8, x3
	str	x8, [x6, x9, lsl #3]
	add	x3, x3, #256
	add	x8, x0, x1
	cmp	x8, x2
	csel	x9, xzr, x2, lo
	sub	x0, x8, x9
	subs	x4, x4, #1
	b.ne	LBB1_1
LBB1_2:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_replay_prefix_sum_guard4       ; -- Begin function replay_prefix_sum_guard4
	.p2align	2
_replay_prefix_sum_guard4:              ; @replay_prefix_sum_guard4
	.cfi_startproc
; %bb.0:
	str	xzr, [x1]
	cbz	x2, LBB2_3
; %bb.1:
	mov	x8, #0                          ; =0x0
	add	x9, x1, #8
LBB2_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x10, [x0], #8
	add	x8, x10, x8
	str	x8, [x9], #8
	subs	x2, x2, #1
	b.ne	LBB2_2
LBB2_3:
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
```

---

## 5. Full Assembly Listing: Sequential Single-Slab Active-Only Pre-Bias, Masked Scan & Direct Probe (`kernels_slab_compact_v81.s`)

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_slab_prebias_active_sl         ; -- Begin function slab_prebias_active_sl
	.p2align	2
_slab_prebias_active_sl:                ; @slab_prebias_active_sl
	.cfi_startproc
; %bb.0:
	add	x8, x2, #19
	lsr	x8, x8, #4
	add	x8, x8, #1
	and	x8, x8, #0x1ffffffffffffffe
	cmp	x8, #9
	b.gt	LBB0_5
; %bb.1:
	cmp	x8, #5
	b.gt	LBB0_9
; %bb.2:
	cmp	x8, #2
	b.eq	LBB0_15
; %bb.3:
	cmp	x8, #4
	ccmp	x1, #0, #4, eq
	b.eq	LBB0_23
LBB0_4:                                 ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_4
	b	LBB0_23
LBB0_5:
	cmp	x8, #13
	b.gt	LBB0_12
; %bb.6:
	cmp	x8, #10
	b.eq	LBB0_17
; %bb.7:
	cmp	x8, #12
	ccmp	x1, #0, #4, eq
	b.eq	LBB0_23
LBB0_8:                                 ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	stp	x3, x3, [x0, #64]
	stp	x3, x3, [x0, #80]
	stp	x3, x3, [x0, #96]
	stp	x3, x3, [x0, #112]
	stp	x3, x3, [x0, #128]
	stp	x3, x3, [x0, #144]
	stp	x3, x3, [x0, #160]
	stp	x3, x3, [x0, #176]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_8
	b	LBB0_23
LBB0_9:
	cmp	x8, #6
	b.eq	LBB0_19
; %bb.10:
	cmp	x8, #8
	ccmp	x1, #0, #4, eq
	b.eq	LBB0_23
LBB0_11:                                ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	stp	x3, x3, [x0, #64]
	stp	x3, x3, [x0, #80]
	stp	x3, x3, [x0, #96]
	stp	x3, x3, [x0, #112]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_11
	b	LBB0_23
LBB0_12:
	cmp	x8, #14
	b.eq	LBB0_21
; %bb.13:
	cmp	x8, #16
	ccmp	x1, #0, #4, eq
	b.eq	LBB0_23
LBB0_14:                                ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	stp	x3, x3, [x0, #64]
	stp	x3, x3, [x0, #80]
	stp	x3, x3, [x0, #96]
	stp	x3, x3, [x0, #112]
	stp	x3, x3, [x0, #128]
	stp	x3, x3, [x0, #144]
	stp	x3, x3, [x0, #160]
	stp	x3, x3, [x0, #176]
	stp	x3, x3, [x0, #192]
	stp	x3, x3, [x0, #208]
	stp	x3, x3, [x0, #224]
	stp	x3, x3, [x0, #240]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_14
	b	LBB0_23
LBB0_15:
	cbz	x1, LBB0_23
LBB0_16:                                ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_16
	b	LBB0_23
LBB0_17:
	cbz	x1, LBB0_23
LBB0_18:                                ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	stp	x3, x3, [x0, #64]
	stp	x3, x3, [x0, #80]
	stp	x3, x3, [x0, #96]
	stp	x3, x3, [x0, #112]
	stp	x3, x3, [x0, #128]
	stp	x3, x3, [x0, #144]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_18
	b	LBB0_23
LBB0_19:
	cbz	x1, LBB0_23
LBB0_20:                                ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	stp	x3, x3, [x0, #64]
	stp	x3, x3, [x0, #80]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_20
	b	LBB0_23
LBB0_21:
	cbz	x1, LBB0_23
LBB0_22:                                ; =>This Inner Loop Header: Depth=1
	stp	x3, x3, [x0]
	stp	x3, x3, [x0, #16]
	stp	x3, x3, [x0, #32]
	stp	x3, x3, [x0, #48]
	stp	x3, x3, [x0, #64]
	stp	x3, x3, [x0, #80]
	stp	x3, x3, [x0, #96]
	stp	x3, x3, [x0, #112]
	stp	x3, x3, [x0, #128]
	stp	x3, x3, [x0, #144]
	stp	x3, x3, [x0, #160]
	stp	x3, x3, [x0, #176]
	stp	x3, x3, [x0, #192]
	stp	x3, x3, [x0, #208]
	add	x0, x0, #256
	subs	x1, x1, #1
	b.ne	LBB0_22
LBB0_23:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_slab_scan_g_active_sl          ; -- Begin function slab_scan_g_active_sl
	.p2align	2
_slab_scan_g_active_sl:                 ; @slab_scan_g_active_sl
	.cfi_startproc
; %bb.0:
	stp	x28, x27, [sp, #-96]!           ; 16-byte Folded Spill
	stp	x26, x25, [sp, #16]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #32]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #48]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #64]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #80]             ; 16-byte Folded Spill
	add	x29, sp, #80
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	add	x9, x2, #4
	add	x8, x2, #19
	lsr	x8, x8, #4
	add	x8, x8, #1
	and	x8, x8, #0x1ffffffffffffffe
	lsl	x10, x8, #4
	subs	x11, x10, #32
	cset	w12, ne
	cmp	x11, x9
	cset	w11, lo
	tst	w12, w11
	mov	w11, #128                       ; =0x80
	csel	x13, x11, xzr, ne
	sub	x11, x10, #31
	cmp	x11, x9
	cset	w11, lo
	tst	w12, w11
	mov	w11, #32768                     ; =0x8000
	csel	x14, x11, xzr, ne
	sub	x11, x10, #30
	cmp	x11, x9
	cset	w11, lo
	tst	w12, w11
	mov	w11, #8388608                   ; =0x800000
	csel	x15, x11, xzr, ne
	sub	x11, x10, #29
	cmp	x11, x9
	cset	w11, lo
	tst	w12, w11
	mov	w11, #-2147483648               ; =0x80000000
	csel	x12, x11, xzr, ne
	sub	x11, x10, #28
	cmp	x11, x9
	cset	w11, lo
	lsl	x11, x11, #39
	sub	x16, x10, #27
	cmp	x16, x9
	cset	w16, lo
	sub	x17, x10, #26
	cmp	x17, x9
	cset	w17, lo
	sub	x4, x10, #25
	cmp	x4, x9
	cset	w4, lo
	orr	x11, x11, x16, lsl #47
	orr	x11, x11, x17, lsl #55
	orr	x11, x11, x4, lsl #63
	orr	x13, x13, x15
	orr	x13, x13, x14
	orr	x12, x13, x12
	sub	x13, x10, #24
	cmp	x13, x9
	cset	w13, lo
	ubfiz	x14, x13, #7, #32
	sub	x13, x10, #23
	cmp	x13, x9
	cset	w15, lo
	sub	x13, x10, #22
	cmp	x13, x9
	cset	w16, lo
	sub	x13, x10, #21
	cmp	x13, x9
	cset	w17, lo
	sub	x13, x10, #20
	cmp	x13, x9
	cset	w4, lo
	sub	x13, x10, #19
	cmp	x13, x9
	cset	w5, lo
	sub	x13, x10, #18
	cmp	x13, x9
	cset	w6, lo
	sub	x13, x10, #17
	cmp	x13, x9
	cset	w13, lo
	orr	x14, x14, x16, lsl #23
	orr	x14, x14, x15, lsl #15
	orr	x14, x14, x17, lsl #31
	orr	x14, x14, x4, lsl #39
	orr	x14, x14, x5, lsl #47
	orr	x14, x14, x6, lsl #55
	sub	x15, x10, #16
	cmp	x15, x9
	cset	w15, lo
	ubfiz	x16, x15, #7, #32
	sub	x15, x10, #15
	cmp	x15, x9
	cset	w17, lo
	sub	x15, x10, #14
	cmp	x15, x9
	cset	w4, lo
	sub	x15, x10, #13
	cmp	x15, x9
	cset	w5, lo
	sub	x15, x10, #12
	cmp	x15, x9
	cset	w6, lo
	sub	x15, x10, #11
	cmp	x15, x9
	cset	w7, lo
	sub	x15, x10, #10
	cmp	x15, x9
	cset	w19, lo
	sub	x15, x10, #9
	cmp	x15, x9
	cset	w15, lo
	orr	x16, x16, x4, lsl #23
	orr	x16, x16, x17, lsl #15
	orr	x16, x16, x5, lsl #31
	orr	x16, x16, x6, lsl #39
	orr	x16, x16, x7, lsl #47
	orr	x16, x16, x19, lsl #55
	cbz	x10, LBB1_10
; %bb.1:
	sub	x17, x10, #8
	cmp	x17, x9
	cset	w17, lo
	ubfiz	x17, x17, #7, #32
	sub	x4, x10, #7
	cmp	x4, x9
	cset	w4, lo
	sub	x5, x10, #6
	cmp	x5, x9
	cset	w5, lo
	sub	x6, x10, #5
	cmp	x6, x9
	cset	w6, lo
	sub	x7, x10, #4
	cmp	x7, x9
	cset	w7, lo
	sub	x19, x10, #3
	cmp	x19, x9
	cset	w19, lo
	sub	x20, x10, #2
	cmp	x20, x9
	cset	w20, lo
	sub	x10, x10, #1
	cmp	x10, x9
	cset	w9, lo
	orr	x10, x17, x5, lsl #23
	orr	x10, x10, x4, lsl #15
	orr	x10, x10, x6, lsl #31
	orr	x10, x10, x7, lsl #39
	orr	x10, x10, x19, lsl #47
	orr	x10, x10, x20, lsl #55
	orr	x23, x10, x9, lsl #63
	mov	x4, #0                          ; =0x0
	orr	x24, x11, x12
	orr	x25, x14, x13, lsl #63
	orr	x26, x16, x15, lsl #63
	cmp	x8, #9
	b.gt	LBB1_11
LBB1_2:
	cmp	x8, #5
	b.gt	LBB1_19
; %bb.3:
	cmp	x8, #2
	b.eq	LBB1_33
; %bb.4:
	cmp	x8, #4
	b.ne	LBB1_54
; %bb.5:
	cbz	x1, LBB1_53
; %bb.6:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x27, #141287244169216           ; =0x808000000000
	movk	x27, #32896, lsl #48
LBB1_7:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	and	x8, x8, x27
	ldp	x10, x11, [x0, #16]
	orr	x9, x10, x9
	orr	x9, x9, x11
	and	x9, x9, #0x8080808080808080
	ldp	x10, x11, [x0, #32]
	and	x10, x10, x24
	and	x11, x11, x25
	ldp	x12, x13, [x0, #48]
	and	x12, x12, x26
	and	x13, x13, x23
	orr	x8, x10, x8
	orr	x8, x8, x11
	orr	x8, x8, x9
	orr	x9, x12, x13
	orr	x8, x8, x9
	cbnz	x8, LBB1_9
; %bb.8:                                ;   in Loop: Header=BB1_7 Depth=1
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_7
	b	LBB1_54
LBB1_9:                                 ;   in Loop: Header=BB1_7 Depth=1
	mov	x20, x0
	mov	x28, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_7
	b	LBB1_54
LBB1_10:
	mov	x23, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	orr	x24, x11, x12
	orr	x25, x14, x13, lsl #63
	orr	x26, x16, x15, lsl #63
	cmp	x8, #9
	b.le	LBB1_2
LBB1_11:
	cmp	x8, #13
	b.gt	LBB1_26
; %bb.12:
	cmp	x8, #10
	b.eq	LBB1_38
; %bb.13:
	cmp	x8, #12
	b.ne	LBB1_54
; %bb.14:
	cbz	x1, LBB1_53
; %bb.15:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x27, #141287244169216           ; =0x808000000000
	movk	x27, #32896, lsl #48
LBB1_16:                                ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	and	x8, x8, x27
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #96]
	ldp	x21, x22, [x0, #112]
	orr	x9, x10, x9
	orr	x10, x11, x12
	orr	x9, x9, x10
	orr	x10, x13, x14
	ldp	x11, x12, [x0, #128]
	orr	x10, x10, x15
	orr	x9, x9, x10
	orr	x10, x16, x17
	orr	x10, x10, x5
	ldp	x13, x14, [x0, #144]
	orr	x10, x10, x6
	orr	x9, x9, x10
	orr	x10, x7, x20
	orr	x10, x10, x21
	orr	x10, x10, x22
	orr	x10, x10, x11
	orr	x9, x9, x10
	orr	x10, x12, x13
	orr	x10, x10, x14
	orr	x9, x9, x10
	and	x9, x9, #0x8080808080808080
	ldp	x10, x11, [x0, #160]
	and	x10, x10, x24
	and	x11, x11, x25
	ldp	x12, x13, [x0, #176]
	and	x12, x12, x26
	and	x13, x13, x23
	orr	x8, x10, x8
	orr	x8, x8, x11
	orr	x8, x8, x12
	orr	x8, x8, x13
	orr	x8, x8, x9
	cbnz	x8, LBB1_18
; %bb.17:                               ;   in Loop: Header=BB1_16 Depth=1
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_16
	b	LBB1_54
LBB1_18:                                ;   in Loop: Header=BB1_16 Depth=1
	mov	x20, x0
	mov	x28, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_16
	b	LBB1_54
LBB1_19:
	cmp	x8, #6
	b.eq	LBB1_43
; %bb.20:
	cmp	x8, #8
	b.ne	LBB1_54
; %bb.21:
	cbz	x1, LBB1_53
; %bb.22:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x27, #141287244169216           ; =0x808000000000
	movk	x27, #32896, lsl #48
LBB1_23:                                ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	and	x8, x8, x27
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	orr	x9, x10, x9
	orr	x10, x11, x12
	orr	x9, x9, x10
	orr	x10, x13, x14
	orr	x10, x10, x15
	orr	x9, x9, x10
	orr	x10, x16, x17
	orr	x10, x10, x5
	orr	x10, x10, x6
	orr	x9, x9, x10
	and	x9, x9, #0x8080808080808080
	ldp	x10, x11, [x0, #96]
	and	x10, x10, x24
	and	x11, x11, x25
	ldp	x12, x13, [x0, #112]
	and	x12, x12, x26
	and	x13, x13, x23
	orr	x8, x10, x8
	orr	x8, x8, x11
	orr	x8, x8, x12
	orr	x8, x8, x13
	orr	x8, x8, x9
	cbnz	x8, LBB1_25
; %bb.24:                               ;   in Loop: Header=BB1_23 Depth=1
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_23
	b	LBB1_54
LBB1_25:                                ;   in Loop: Header=BB1_23 Depth=1
	mov	x20, x0
	mov	x28, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_23
	b	LBB1_54
LBB1_26:
	cmp	x8, #14
	b.eq	LBB1_48
; %bb.27:
	cmp	x8, #16
	b.ne	LBB1_54
; %bb.28:
	cbz	x1, LBB1_53
; %bb.29:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x22, #141287244169216           ; =0x808000000000
	movk	x22, #32896, lsl #48
LBB1_30:                                ; =>This Inner Loop Header: Depth=1
	add	x8, x0, x19
	ldp	x9, x10, [x8, #8]
	orr	x9, x10, x9
	ldp	x11, x10, [x8, #24]
	orr	x10, x11, x10
	orr	x9, x9, x10
	ldp	x12, x10, [x8, #40]
	ldr	x11, [x8, #56]
	orr	x10, x12, x10
	orr	x10, x10, x11
	ldp	x11, x12, [x8, #64]
	orr	x9, x9, x10
	orr	x10, x11, x12
	ldp	x11, x12, [x8, #80]
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #96]
	orr	x9, x9, x10
	orr	x10, x11, x12
	ldp	x11, x12, [x8, #112]
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #128]
	orr	x10, x10, x11
	orr	x9, x9, x10
	ldp	x10, x11, [x8, #144]
	orr	x10, x12, x10
	orr	x10, x10, x11
	ldp	x11, x12, [x8, #160]
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #176]
	orr	x10, x10, x11
	orr	x9, x9, x10
	ldp	x10, x11, [x8, #192]
	orr	x10, x12, x10
	orr	x10, x10, x11
	ldp	x11, x12, [x8, #208]
	orr	x10, x10, x11
	ldr	x11, [x8]
	and	x11, x11, x22
	orr	x10, x10, x12
	orr	x9, x9, x10
	ldp	x10, x12, [x8, #224]
	and	x10, x10, x24
	and	x12, x12, x25
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #240]
	and	x11, x11, x26
	orr	x10, x10, x11
	and	x11, x12, x23
	orr	x10, x10, x11
	and	x9, x9, #0x8080808080808080
	orr	x9, x10, x9
	cbnz	x9, LBB1_32
; %bb.31:                               ;   in Loop: Header=BB1_30 Depth=1
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_30
	b	LBB1_54
LBB1_32:                                ;   in Loop: Header=BB1_30 Depth=1
	mov	x27, x0
	mov	x0, x8
	mov	x28, x1
	mov	x1, x19
	mov	x20, x2
	mov	x21, x3
	bl	_extract_row_survivors
	mov	x2, x20
	mov	x3, x21
	mov	x4, x0
	mov	x0, x27
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_30
	b	LBB1_54
LBB1_33:
	cbz	x1, LBB1_53
; %bb.34:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
LBB1_35:                                ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	and	x8, x8, x24
	and	x9, x9, x25
	orr	x8, x9, x8
	ldp	x9, x10, [x0, #16]
	and	x9, x9, x26
	and	x10, x10, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	cbnz	x8, LBB1_37
; %bb.36:                               ;   in Loop: Header=BB1_35 Depth=1
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_35
	b	LBB1_54
LBB1_37:                                ;   in Loop: Header=BB1_35 Depth=1
	mov	x20, x0
	mov	x27, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x27, #1
	b.ne	LBB1_35
	b	LBB1_54
LBB1_38:
	cbz	x1, LBB1_53
; %bb.39:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x27, #141287244169216           ; =0x808000000000
	movk	x27, #32896, lsl #48
LBB1_40:                                ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	and	x8, x8, x27
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #96]
	ldp	x21, x22, [x0, #112]
	orr	x9, x10, x9
	orr	x10, x11, x12
	orr	x9, x9, x10
	orr	x10, x13, x14
	orr	x10, x10, x15
	orr	x9, x9, x10
	orr	x10, x16, x17
	orr	x10, x10, x5
	orr	x10, x10, x6
	orr	x9, x9, x10
	orr	x10, x7, x20
	orr	x10, x10, x21
	orr	x10, x10, x22
	orr	x9, x9, x10
	and	x9, x9, #0x8080808080808080
	ldp	x10, x11, [x0, #128]
	and	x10, x10, x24
	and	x11, x11, x25
	ldp	x12, x13, [x0, #144]
	and	x12, x12, x26
	and	x13, x13, x23
	orr	x8, x10, x8
	orr	x8, x8, x11
	orr	x8, x8, x12
	orr	x8, x8, x13
	orr	x8, x8, x9
	cbnz	x8, LBB1_42
; %bb.41:                               ;   in Loop: Header=BB1_40 Depth=1
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_40
	b	LBB1_54
LBB1_42:                                ;   in Loop: Header=BB1_40 Depth=1
	mov	x20, x0
	mov	x28, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_40
	b	LBB1_54
LBB1_43:
	cbz	x1, LBB1_53
; %bb.44:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x27, #141287244169216           ; =0x808000000000
	movk	x27, #32896, lsl #48
LBB1_45:                                ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	and	x8, x8, x27
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	orr	x9, x10, x9
	orr	x10, x11, x12
	orr	x9, x9, x10
	orr	x10, x13, x14
	orr	x10, x10, x15
	orr	x9, x9, x10
	and	x9, x9, #0x8080808080808080
	ldp	x10, x11, [x0, #64]
	and	x10, x10, x24
	and	x11, x11, x25
	ldp	x12, x13, [x0, #80]
	and	x12, x12, x26
	and	x13, x13, x23
	orr	x8, x10, x8
	orr	x8, x8, x11
	orr	x8, x8, x12
	orr	x8, x8, x9
	orr	x8, x8, x13
	cbnz	x8, LBB1_47
; %bb.46:                               ;   in Loop: Header=BB1_45 Depth=1
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_45
	b	LBB1_54
LBB1_47:                                ;   in Loop: Header=BB1_45 Depth=1
	mov	x20, x0
	mov	x28, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_45
	b	LBB1_54
LBB1_48:
	cbz	x1, LBB1_53
; %bb.49:
	mov	x19, #0                         ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x22, #141287244169216           ; =0x808000000000
	movk	x22, #32896, lsl #48
LBB1_50:                                ; =>This Inner Loop Header: Depth=1
	add	x8, x0, x19
	ldp	x9, x10, [x8, #8]
	orr	x9, x10, x9
	ldp	x11, x10, [x8, #24]
	orr	x10, x11, x10
	orr	x9, x9, x10
	ldp	x12, x10, [x8, #40]
	ldr	x11, [x8, #56]
	orr	x10, x12, x10
	orr	x10, x10, x11
	ldp	x11, x12, [x8, #64]
	orr	x9, x9, x10
	orr	x10, x11, x12
	ldp	x11, x12, [x8, #80]
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #96]
	orr	x9, x9, x10
	orr	x10, x11, x12
	ldp	x11, x12, [x8, #112]
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #128]
	orr	x10, x10, x11
	orr	x9, x9, x10
	ldp	x10, x11, [x8, #144]
	orr	x10, x12, x10
	orr	x10, x10, x11
	ldp	x11, x12, [x8, #160]
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #176]
	orr	x10, x10, x11
	ldr	x11, [x8]
	and	x11, x11, x22
	orr	x9, x9, x10
	orr	x9, x9, x12
	ldp	x10, x12, [x8, #192]
	and	x10, x10, x24
	and	x12, x12, x25
	orr	x10, x10, x11
	orr	x10, x10, x12
	ldp	x11, x12, [x8, #208]
	and	x11, x11, x26
	orr	x10, x10, x11
	and	x11, x12, x23
	orr	x10, x10, x11
	and	x9, x9, #0x8080808080808080
	orr	x9, x10, x9
	cbnz	x9, LBB1_52
; %bb.51:                               ;   in Loop: Header=BB1_50 Depth=1
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_50
	b	LBB1_54
LBB1_52:                                ;   in Loop: Header=BB1_50 Depth=1
	mov	x27, x0
	mov	x0, x8
	mov	x28, x1
	mov	x1, x19
	mov	x20, x2
	mov	x21, x3
	bl	_extract_row_survivors
	mov	x2, x20
	mov	x3, x21
	mov	x4, x0
	mov	x0, x27
	add	x19, x19, #256
	subs	x1, x28, #1
	b.ne	LBB1_50
	b	LBB1_54
LBB1_53:
	mov	x4, #0                          ; =0x0
LBB1_54:
	mov	x0, x4
	ldp	x29, x30, [sp, #80]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #32]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #16]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp], #96             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function extract_row_survivors
_extract_row_survivors:                 ; @extract_row_survivors
	.cfi_startproc
; %bb.0:
	add	x9, x2, #4
	cmp	x9, #4
	b.ls	LBB2_13
; %bb.1:
	ldrsb	w8, [x0, #4]
	tbnz	w8, #31, LBB2_3
; %bb.2:
	cmp	x9, #5
	b.ne	LBB2_4
	b	LBB2_12
LBB2_3:
	add	x8, x1, #4
	str	x8, [x3, x4, lsl #3]
	add	x4, x4, #1
	cmp	x9, #5
	b.eq	LBB2_12
LBB2_4:
	ldrsb	w8, [x0, #5]
	tbz	w8, #31, LBB2_6
; %bb.5:
	add	x8, x1, #5
	str	x8, [x3, x4, lsl #3]
	add	x4, x4, #1
LBB2_6:
	cmp	x9, #6
	b.eq	LBB2_12
; %bb.7:
	ldrsb	w8, [x0, #6]
	tbz	w8, #31, LBB2_9
; %bb.8:
	add	x8, x1, #6
	str	x8, [x3, x4, lsl #3]
	add	x4, x4, #1
LBB2_9:
	cmp	x9, #7
	b.eq	LBB2_12
; %bb.10:
	ldrsb	w8, [x0, #7]
	tbz	w8, #31, LBB2_12
; %bb.11:
	add	x8, x1, #7
	str	x8, [x3, x4, lsl #3]
	add	x4, x4, #1
LBB2_12:
	cmp	x9, #9
	b.hs	LBB2_14
LBB2_13:
	mov	x0, x4
	ret
LBB2_14:
	cmp	x9, #16
	b.lo	LBB2_21
; %bb.15:
	lsr	x10, x9, #3
	mov	w11, #1                         ; =0x1
	b	LBB2_18
LBB2_16:                                ;   in Loop: Header=BB2_18 Depth=1
	mov	x8, x4
LBB2_17:                                ;   in Loop: Header=BB2_18 Depth=1
	add	x11, x11, #1
	mov	x4, x8
	cmp	x11, x10
	b.eq	LBB2_22
LBB2_18:                                ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB2_20 Depth 2
	ldr	x8, [x0, x11, lsl #3]
	ands	x12, x8, #0x8080808080808080
	b.eq	LBB2_16
; %bb.19:                               ;   in Loop: Header=BB2_18 Depth=1
	add	x13, x1, x11, lsl #3
LBB2_20:                                ;   Parent Loop BB2_18 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	rbit	x8, x12
	clz	x8, x8
	add	x14, x13, x8, lsr #3
	add	x8, x4, #1
	str	x14, [x3, x4, lsl #3]
	sub	x14, x12, #1
	mov	x4, x8
	ands	x12, x14, x12
	b.ne	LBB2_20
	b	LBB2_17
LBB2_21:
	mov	x8, x4
LBB2_22:
	and	x10, x9, #0xfffffffffffffff8
	cmp	x10, x9
	b.ne	LBB2_25
LBB2_23:
	mov	x0, x8
	ret
LBB2_24:                                ;   in Loop: Header=BB2_25 Depth=1
	add	x10, x10, #1
	cmp	x10, x9
	b.hs	LBB2_23
LBB2_25:                                ; =>This Inner Loop Header: Depth=1
	ldrsb	w11, [x0, x10]
	tbz	w11, #31, LBB2_24
; %bb.26:                               ;   in Loop: Header=BB2_25 Depth=1
	add	x11, x1, x10
	str	x11, [x3, x8, lsl #3]
	add	x8, x8, #1
	b	LBB2_24
	.cfi_endproc
                                        ; -- End function
	.globl	_slab_filter_f_at_gsurv         ; -- Begin function slab_filter_f_at_gsurv
	.p2align	2
_slab_filter_f_at_gsurv:                ; @slab_filter_f_at_gsurv
	.cfi_startproc
; %bb.0:
	mov	x8, #0                          ; =0x0
	cbnz	x2, LBB3_3
LBB3_1:
	mov	x0, x8
	ret
LBB3_2:                                 ;   in Loop: Header=BB3_3 Depth=1
	subs	x2, x2, #1
	b.eq	LBB3_1
LBB3_3:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x9, [x1], #8
	ldrsb	w10, [x0, x9]
	tbz	w10, #31, LBB3_2
; %bb.4:                                ;   in Loop: Header=BB3_3 Depth=1
	str	x9, [x3, x8, lsl #3]
	add	x8, x8, #1
	b	LBB3_2
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
```

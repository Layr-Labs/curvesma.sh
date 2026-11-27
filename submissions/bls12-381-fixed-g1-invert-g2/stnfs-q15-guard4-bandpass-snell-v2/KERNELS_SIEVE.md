# Stage 1a Zero-Spill Guard4 x64 Sieve, Replay-Build, and 16B-Granular Adaptive Slab Kernels (`KERNELS_SIEVE.md`)

All assembly listings in this file were compiled with Apple Clang (`arm64`) under strict scalar general-purpose-register constraints:

```text
-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize
```

## 1. ISA Compliance, Zero-Spill `STADDB` Disclosure & Guard-Masked Correctness

- **SIMD / FP / SVE / Platform-Register Audit:** `simd = 0`. Zero NEON/FP/SVE vector registers (`v0..v31`, `q0..q31`, `d0..d31`, `s0..s31`, `h0..h31`, `z0..z31`, `p0..p15`), zero floating-point instructions, zero `dc zva`, and zero `x18` accesses appear anywhere in this file.
- **Explicit ARMv8.1-A Scalar LSE Disclosure:** The zero-spill Guard4 x64 hit-emission kernel `_emit_u64_x64_direct` (invoked by `_sieve_window_guard4_x64`) uses the scalar ARMv8.1-A Large System Extensions (LSE) byte atomic-add instruction **`STADDB Ws, [Xn]`** (a 64-bit general-purpose-register-only memory add with zero SIMD state, counted at `1` instruction each). Every other kernel in this file (`_replay_*`, `_slab_prebias_adaptive_sl`, `_slab_scan_g_active_16`, `_slab_filter_f_at_gsurv`) uses pure ARMv8.0 scalar instructions only (`0` `staddb`).
- **Zero-Spill `x64 STADDB` (`_emit_u64_x64_direct`):** By emitting in 8-hit caller-saved bursts (`x10..x17`), `_emit_u64_x64_direct` spills **ZERO callee-saved registers (`sp = 0`)**, cutting callee dispatch overhead from `17` Ir (`x16`) down to **`14` Ir** and prime probe bounds from `56 / 84` Ir down to **`51 / 77` Ir** (non-split / split-256; note `_sieve_window_guard4_x64`'s outer prime-probe loop stack reloads are fully counted inside the `51 / 77` Ir probe bounds), while achieving **`163/64 = 2.546875 Ir/hit`** in the `x64` loop (`LBB1_2`), **`23/8 = 2.875000 Ir/hit`** in the `x8` loop (`LBB1_4`), and **`5.000000 Ir/hit`** in the `x1` tail (`LBB1_6`). Exact Poisson expectation over Chebotarev root counts and split-256 slices yields `hits_g = 1.9218`, `hits_f = 1.5691` (`ipc_hits = 3.4909`, `ipc_probe = 0.3209`).
- **16B-Granular 2-Row Adaptive Prebias & Active-Only Masked Scan (`kernels_slab_adaptive_v81.s`):**
  - Each `a_id` row occupies 256 bytes: 4-byte left guard `[0..3]`, $L \le 248$ active columns `[4..4+L-1]`, and $\ge 4$-byte right guard `[4+L..255]`.
  - ONE single `14.9616 EiB` flat slab is reused sequentially ($S_g \to S_f$):
    1. **Pass 1 ($g$-side):** `_slab_prebias_adaptive_sl` broadcasts a per-column 256-byte adaptive bias template `bias_g[c] = 128 - T_g(shell, N_g(b_c), nf0(b_c))` across 2 rows per iteration with **`0` stack spills across all 16 switch cases `k16 = 1..16`** (`k16 + 2.0` Ir/row for `k16 = 1..12` covering `L <= 188` with `0` in-loop loads; `2*k16 - 7.0` Ir/row for `k16 = 13..16`). After `_sieve_window_guard4_x64`, `_slab_scan_g_active_16` scans strictly active columns at 16-byte granularity (`9.0` Ir/row at `k16 = 1`; `3*k16 + 8.0` Ir/row at `k16 = 2..16`), **explicitly masking BOTH the 4-byte left guard `[0..3]` AND the right guard `[4+L..255]`** before sign-bit detection (`248/248` lengths `L in [1..248]` verified bit-exact under `0xFF` guard poison!).
    2. **Pass 2 ($f$-side, SAME slab):** `_slab_prebias_adaptive_sl` writes the per-column adaptive $f$-bias template `bias_f[c]`, `_sieve_window_guard4_x64` emits $f$-band hits, and `_slab_filter_f_at_gsurv` probes ONLY at the `2^{-14.5751}` $g$-survivor offsets (`5` Ir per non-passing $g$-survivor).
  - Exact `200,000`-chord average at `CAND2_FINAL_ETA85_7000_NL205` (`L_win_eff = 111.534`): `0.1924` (2-pass adaptive prebias) `+ 0.2775` (masked $S_g$ scan) `+ 0.00774` (`189 * 2^{-14.5751}` slow-path + $S_f$ probe) `= 0.4777 Ir/active cell`, conservatively charged **`0.4950 Ir/active cell`**.

## 2. Exact Basic-Block Instruction Counts

| Function | Basic Block(s) | Exact Compiled Count | Charged in Ledger (`proof.md`) | Notes |
|---|---|---:|---:|---|
| `_emit_u64_x64_direct` | `LBB1_2` (x64 unrolled hits) | `163 Ir / 64 = 2.546875 Ir/hit` | **`2.546875 Ir/hit`** | `32 ldp + 65 add + 64 staddb + cmp + b.lo`, `0` spills |
| `_emit_u64_x64_direct` | `LBB1_4` (x8 hits) | `23 Ir / 8 = 2.875000 Ir/hit` | **`2.875000 Ir/hit`** | `4 ldp + 9 add + 8 staddb + cmp + b.lo`, `0` spills |
| `_emit_u64_x64_direct` | `LBB1_6` (x1 tail) | `5 Ir / 1 = 5.000000 Ir/hit` | **`5.000000 Ir/hit`** | `ldr + add + staddb + cmp + b.lo` |
| `_sieve_window_guard4_x64` | Non-split probe + `x64` call | `<= 51 Ir/prime` | **`51.0 Ir/prime`** | Includes `14` Ir zero-spill callee dispatch |
| `_sieve_window_guard4_x64` | Split-256 probe + 2 calls | `<= 77 Ir/prime` | **`77.0 Ir/prime`** | `p_split = L_win_eff / 256 = 0.4357` |
| `_replay_pass1_count_guard4` | `LBB0_1` (counting pass) | `11 Ir/entry` | Part of **`27.5 Ir/entry`** | `11` (Pass 1) `+ 15` (Pass 2) `+ 1.5` startup |
| `_replay_pass2_scatter_guard4` | `LBB1_1` (scatter pass) | `15 Ir/entry` | Part of **`27.5 Ir/entry`** | `26 Ir` compiled `+ 1.5` pad `= 27.5 Ir/entry` |
| `_replay_prefix_sum_guard4` | `LBB2_2` (prefix sum) | `5 Ir/bucket` | **`6.5 Ir/bucket`** | `+1.5 Ir/bucket` zeroing cushion |
| `_slab_prebias_adaptive_sl` | 16-case 2-row switch (`k16=1..16`) | `k16+2` (`k<=12`), `2k16-7` (`k>12`) | `0.1924` (in **`0.4950`**) | `0` stack spills across all 16 cases |
| `_slab_scan_g_active_16` | 16-case masked SWAR (`k16=1..16`) | `9` (`k=1`), `3k16+8` (`k>=2`) | `0.2775` (in **`0.4950`**) | Masks left `[0..3]` & right `[4+L..255]`, `0` spills |
| `_slab_filter_f_at_gsurv` | `LBB2_2` (direct $S_f$ probe) | `5 Ir` non-surv / `8 Ir` surv | `0.0077` (in **`0.4950`**) | `189 * 2^{-14.5751} = 0.00774` slow+probe |

---

## 3. Full Assembly Listing: Zero-Spill Guard4 x64 Chord-Window Sieve (`kernels_guard4_x32_x64_v81.s`)

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_emit_u64_x64_direct            ; -- Begin function emit_u64_x64_direct
	.p2align	2
_emit_u64_x64_direct:                   ; @emit_u64_x64_direct
	.cfi_startproc
; %bb.0:
	and	x9, x2, #0xffffffffffffffc0
	mov	x8, x1
	cmp	x9, #1
	b.lt	LBB1_3
; %bb.1:
	add	x9, x1, x9, lsl #3
	mov	x8, x1
LBB1_2:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x10, x11, [x8]
	ldp	x12, x13, [x8, #16]
	ldp	x14, x15, [x8, #32]
	ldp	x16, x17, [x8, #48]
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
	ldp	x10, x11, [x8, #64]
	ldp	x12, x13, [x8, #80]
	ldp	x14, x15, [x8, #96]
	ldp	x16, x17, [x8, #112]
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
	ldp	x10, x11, [x8, #128]
	ldp	x12, x13, [x8, #144]
	ldp	x14, x15, [x8, #160]
	ldp	x16, x17, [x8, #176]
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
	ldp	x10, x11, [x8, #192]
	ldp	x12, x13, [x8, #208]
	ldp	x14, x15, [x8, #224]
	ldp	x16, x17, [x8, #240]
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
	ldp	x10, x11, [x8, #256]
	ldp	x12, x13, [x8, #272]
	ldp	x14, x15, [x8, #288]
	ldp	x16, x17, [x8, #304]
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
	ldp	x10, x11, [x8, #320]
	ldp	x12, x13, [x8, #336]
	ldp	x14, x15, [x8, #352]
	ldp	x16, x17, [x8, #368]
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
	ldp	x10, x11, [x8, #384]
	ldp	x12, x13, [x8, #400]
	ldp	x14, x15, [x8, #416]
	ldp	x16, x17, [x8, #432]
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
	ldp	x10, x11, [x8, #448]
	ldp	x12, x13, [x8, #464]
	ldp	x14, x15, [x8, #480]
	ldp	x16, x17, [x8, #496]
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
	add	x8, x8, #512
	cmp	x8, x9
	b.lo	LBB1_2
LBB1_3:
	and	x9, x2, #0x1ffffffffffffff8
	add	x9, x1, x9, lsl #3
	cmp	x8, x9
	b.hs	LBB1_5
LBB1_4:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x10, x11, [x8]
	ldp	x12, x13, [x8, #16]
	ldp	x14, x15, [x8, #32]
	ldp	x16, x17, [x8, #48]
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
	add	x8, x8, #64
	cmp	x8, x9
	b.lo	LBB1_4
LBB1_5:
	add	x9, x1, x2, lsl #3
	cmp	x8, x9
	b.hs	LBB1_7
LBB1_6:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x10, [x8], #8
	add	x10, x0, x10
	; InlineAsm Start
	staddb	w3, [x10]
	; InlineAsm End
	cmp	x8, x9
	b.lo	LBB1_6
LBB1_7:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_sieve_window_guard4_x64        ; -- Begin function sieve_window_guard4_x64
	.p2align	2
_sieve_window_guard4_x64:               ; @sieve_window_guard4_x64
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
	cbz	x2, LBB2_22
; %bb.1:
	mov	x19, x3
	mov	x20, x2
	mov	x22, x1
	mov	x21, #0                         ; =0x0
	add	x8, x0, #4
	stp	x8, x0, [sp, #32]               ; 16-byte Folded Spill
	b	LBB2_4
LBB2_2:                                 ;   in Loop: Header=BB2_4 Depth=1
	and	w3, w26, #0xff
	mov	x2, x24
	bl	_emit_u64_x64_direct
	add	x21, x24, x21
LBB2_3:                                 ;   in Loop: Header=BB2_4 Depth=1
	add	x22, x22, #48
	subs	x20, x20, #1
	b.eq	LBB2_23
LBB2_4:                                 ; =>This Inner Loop Header: Depth=1
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
	b.hi	LBB2_12
; %bb.5:                                ;   in Loop: Header=BB2_4 Depth=1
	sub	x10, x10, #1
	lsr	x11, x10, #2
	and	x9, x9, #0xfffffffffffffff8
	ldr	x9, [x28, x9]
	and	x25, x25, #0xff
	cmp	x8, x10, lsr #8
	add	x23, x28, x11, lsl #3
	b.ne	LBB2_8
; %bb.6:                                ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x24, x8, x9
	b.ls	LBB2_3
; %bb.7:                                ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x25
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	b	LBB2_2
LBB2_8:                                 ;   in Loop: Header=BB2_4 Depth=1
	add	x8, x28, x8, lsl #9
	ldr	x28, [x8, #512]
	subs	x24, x28, x9
	b.ls	LBB2_10
; %bb.9:                                ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x25
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	and	w3, w26, #0xff
	mov	x2, x24
	bl	_emit_u64_x64_direct
	add	x21, x24, x21
LBB2_10:                                ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x24, x8, x28
	b.ls	LBB2_3
; %bb.11:                               ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x25
	add	x1, x27, x28, lsl #3
	add	x0, x8, #260
	b	LBB2_2
LBB2_12:                                ;   in Loop: Header=BB2_4 Depth=1
	sub	x10, x24, #1
	lsr	x11, x10, #2
	and	x9, x9, #0xfffffffffffffff8
	ldr	x9, [x28, x9]
	and	x12, x25, #0xff
	cmp	x8, x10, lsr #8
	add	x23, x28, x11, lsl #3
	b.ne	LBB2_15
; %bb.13:                               ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x23, x8, x9
	b.ls	LBB2_20
; %bb.14:                               ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x12
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	b	LBB2_19
LBB2_15:                                ;   in Loop: Header=BB2_4 Depth=1
	add	x8, x28, x8, lsl #9
	ldr	x10, [x8, #512]
	subs	x11, x10, x9
	str	x12, [sp, #24]                  ; 8-byte Folded Spill
	b.ls	LBB2_17
; %bb.16:                               ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	sub	x8, x8, x12
	add	x1, x27, x9, lsl #3
	add	x0, x8, #4
	and	w3, w26, #0xff
	mov	x2, x11
	stp	x11, x10, [sp, #8]              ; 16-byte Folded Spill
	bl	_emit_u64_x64_direct
	ldp	x8, x10, [sp, #8]               ; 16-byte Folded Reload
	add	x21, x8, x21
LBB2_17:                                ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [x23, #8]
	subs	x23, x8, x10
	b.ls	LBB2_20
; %bb.18:                               ;   in Loop: Header=BB2_4 Depth=1
	ldr	x8, [sp, #40]                   ; 8-byte Folded Reload
	ldr	x9, [sp, #24]                   ; 8-byte Folded Reload
	sub	x8, x8, x9
	add	x1, x27, x10, lsl #3
	add	x0, x8, #260
LBB2_19:                                ;   in Loop: Header=BB2_4 Depth=1
	and	w3, w26, #0xff
	mov	x2, x23
	bl	_emit_u64_x64_direct
	add	x21, x23, x21
LBB2_20:                                ;   in Loop: Header=BB2_4 Depth=1
	sub	x8, x24, x25
	mvn	x9, x8
	add	x10, x9, x19
	ldr	x9, [x28]
	lsr	x10, x10, #2
	add	x10, x28, x10, lsl #3
	ldr	x10, [x10, #8]
	subs	x24, x10, x9
	b.ls	LBB2_3
; %bb.21:                               ;   in Loop: Header=BB2_4 Depth=1
	add	x1, x27, x9, lsl #3
	ldr	x9, [sp, #32]                   ; 8-byte Folded Reload
	add	x0, x9, x8
	b	LBB2_2
LBB2_22:
	mov	x21, #0                         ; =0x0
LBB2_23:
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

## 5. Full Assembly Listing: 16B-Granular Adaptive Prebias, Masked Active Scan & $S_f$ Probe (`kernels_slab_adaptive_v81.s`)

*(Note: only `.cfi_*` unwind metadata and trailing `;` comments were stripped to stay under the 65,536-byte file cap; 100% of executable instructions, labels, and jump tables across all `k16 = 1..16` cases are preserved verbatim.)*

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_slab_prebias_adaptive_sl
	.p2align	2
_slab_prebias_adaptive_sl:
	add	x8, x2, #19
	lsr	x8, x8, #4
	sub	x8, x8, #1
	cmp	x8, #15
	b.hi	LBB0_99
	stp	x28, x27, [sp, #-80]!
	stp	x26, x25, [sp, #16]
	stp	x24, x23, [sp, #32]
	stp	x22, x21, [sp, #48]
	stp	x20, x19, [sp, #64]
Lloh0:
	adrp	x9, lJTI0_0@PAGE
Lloh1:
	add	x9, x9, lJTI0_0@PAGEOFF
	adr	x10, LBB0_2
	ldrh	w11, [x9, x8, lsl #1]
	add	x10, x10, x11, lsl #2
	br	x10
LBB0_2:
	ldp	x8, x9, [x3]
	cmp	x1, #2
	b.lo	LBB0_77
	and	x10, x1, #0xfffffffffffffffe
	mov	w11, #1
LBB0_4:
	stp	x8, x9, [x0]
	stp	x8, x9, [x0, #256]
	add	x0, x0, #512
	add	x11, x11, #2
	cmp	x11, x1
	b.lo	LBB0_4
	b	LBB0_78
LBB0_5:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x23, [x3, #144]
	cmp	x1, #2
	b.lo	LBB0_50
	mov	x24, #0
LBB0_7:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x23, [x0, #400]
	ldr	x25, [x3, #160]
	ldr	x26, [x3, #168]
	stp	x25, x26, [x0, #160]
	stp	x25, x26, [x0, #416]
	ldr	x25, [x3, #176]
	ldr	x26, [x3, #184]
	stp	x25, x26, [x0, #176]
	stp	x25, x26, [x0, #432]
	ldr	x25, [x3, #192]
	ldr	x26, [x3, #200]
	stp	x25, x26, [x0, #192]
	stp	x25, x26, [x0, #448]
	ldr	x25, [x3, #208]
	ldr	x26, [x3, #216]
	stp	x25, x26, [x0, #208]
	stp	x25, x26, [x0, #464]
	ldr	x26, [x3, #224]
	ldr	x27, [x3, #232]
	add	x25, x24, #2
	stp	x26, x27, [x0, #224]
	stp	x26, x27, [x0, #480]
	add	x0, x0, #512
	add	x26, x24, #3
	mov	x24, x25
	cmp	x26, x1
	b.lo	LBB0_7
	b	LBB0_51
LBB0_8:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x23, [x3, #144]
	ldp	x24, x25, [x3, #160]
	ldp	x26, x3, [x3, #176]
	cmp	x1, #2
	b.lo	LBB0_53
	and	x27, x1, #0xfffffffffffffffe
	mov	w28, #1
LBB0_10:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x24, x25, [x0, #160]
	stp	x26, x3, [x0, #176]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x23, [x0, #400]
	stp	x24, x25, [x0, #416]
	stp	x26, x3, [x0, #432]
	add	x0, x0, #512
	add	x28, x28, #2
	cmp	x28, x1
	b.lo	LBB0_10
	b	LBB0_54
LBB0_11:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	cmp	x1, #2
	b.lo	LBB0_56
	and	x2, x1, #0xfffffffffffffffe
	mov	w3, #1
LBB0_13:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	add	x0, x0, #512
	add	x3, x3, #2
	cmp	x3, x1
	b.lo	LBB0_13
	b	LBB0_57
LBB0_14:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x23, [x3, #144]
	cmp	x1, #2
	b.lo	LBB0_59
	mov	x24, #0
LBB0_16:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x23, [x0, #400]
	ldr	x25, [x3, #160]
	ldr	x26, [x3, #168]
	stp	x25, x26, [x0, #160]
	stp	x25, x26, [x0, #416]
	ldr	x25, [x3, #176]
	ldr	x26, [x3, #184]
	stp	x25, x26, [x0, #176]
	stp	x25, x26, [x0, #432]
	ldr	x26, [x3, #192]
	ldr	x27, [x3, #200]
	add	x25, x24, #2
	stp	x26, x27, [x0, #192]
	stp	x26, x27, [x0, #448]
	add	x0, x0, #512
	add	x26, x24, #3
	mov	x24, x25
	cmp	x26, x1
	b.lo	LBB0_16
	b	LBB0_60
LBB0_17:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x3, [x3, #144]
	cmp	x1, #2
	b.lo	LBB0_62
	and	x23, x1, #0xfffffffffffffffe
	mov	w24, #1
LBB0_19:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x3, [x0, #144]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x3, [x0, #400]
	add	x0, x0, #512
	add	x24, x24, #2
	cmp	x24, x1
	b.lo	LBB0_19
	b	LBB0_63
LBB0_20:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	cmp	x1, #2
	b.lo	LBB0_65
	and	x14, x1, #0xfffffffffffffffe
	mov	w15, #1
LBB0_22:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	add	x0, x0, #512
	add	x15, x15, #2
	cmp	x15, x1
	b.lo	LBB0_22
	b	LBB0_66
LBB0_23:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x23, [x3, #144]
	cmp	x1, #2
	b.lo	LBB0_68
	mov	x24, #0
LBB0_25:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x23, [x0, #400]
	ldr	x25, [x3, #160]
	ldr	x26, [x3, #168]
	stp	x25, x26, [x0, #160]
	stp	x25, x26, [x0, #416]
	ldr	x25, [x3, #176]
	ldr	x26, [x3, #184]
	stp	x25, x26, [x0, #176]
	stp	x25, x26, [x0, #432]
	ldr	x25, [x3, #192]
	ldr	x26, [x3, #200]
	stp	x25, x26, [x0, #192]
	stp	x25, x26, [x0, #448]
	ldr	x25, [x3, #208]
	ldr	x26, [x3, #216]
	stp	x25, x26, [x0, #208]
	stp	x25, x26, [x0, #464]
	ldr	x25, [x3, #224]
	ldr	x26, [x3, #232]
	stp	x25, x26, [x0, #224]
	stp	x25, x26, [x0, #480]
	ldr	x26, [x3, #240]
	ldr	x27, [x3, #248]
	add	x25, x24, #2
	stp	x26, x27, [x0, #240]
	stp	x26, x27, [x0, #496]
	add	x0, x0, #512
	add	x26, x24, #3
	mov	x24, x25
	cmp	x26, x1
	b.lo	LBB0_25
	b	LBB0_69
LBB0_26:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	cmp	x1, #2
	b.lo	LBB0_71
	and	x16, x1, #0xfffffffffffffffe
	mov	w17, #1
LBB0_28:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	add	x0, x0, #512
	add	x17, x17, #2
	cmp	x17, x1
	b.lo	LBB0_28
	b	LBB0_72
LBB0_29:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x3, [x3, #112]
	cmp	x1, #2
	b.lo	LBB0_74
	and	x19, x1, #0xfffffffffffffffe
	mov	w20, #1
LBB0_31:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x3, [x0, #112]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x3, [x0, #368]
	add	x0, x0, #512
	add	x20, x20, #2
	cmp	x20, x1
	b.lo	LBB0_31
	b	LBB0_75
LBB0_32:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	cmp	x1, #2
	b.lo	LBB0_80
	and	x12, x1, #0xfffffffffffffffe
	mov	w13, #1
LBB0_34:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	add	x0, x0, #512
	add	x13, x13, #2
	cmp	x13, x1
	b.lo	LBB0_34
	b	LBB0_81
LBB0_35:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x3, [x3, #80]
	cmp	x1, #2
	b.lo	LBB0_83
	and	x4, x1, #0xfffffffffffffffe
	mov	w5, #1
LBB0_37:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x3, [x0, #80]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x3, [x0, #336]
	add	x0, x0, #512
	add	x5, x5, #2
	cmp	x5, x1
	b.lo	LBB0_37
	b	LBB0_84
LBB0_38:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x23, [x3, #144]
	cmp	x1, #2
	b.lo	LBB0_86
	mov	x24, #0
LBB0_40:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x23, [x0, #400]
	ldr	x25, [x3, #160]
	ldr	x26, [x3, #168]
	stp	x25, x26, [x0, #160]
	stp	x25, x26, [x0, #416]
	ldr	x25, [x3, #176]
	ldr	x26, [x3, #184]
	stp	x25, x26, [x0, #176]
	stp	x25, x26, [x0, #432]
	ldr	x25, [x3, #192]
	ldr	x26, [x3, #200]
	stp	x25, x26, [x0, #192]
	stp	x25, x26, [x0, #448]
	ldr	x26, [x3, #208]
	ldr	x27, [x3, #216]
	add	x25, x24, #2
	stp	x26, x27, [x0, #208]
	stp	x26, x27, [x0, #464]
	add	x0, x0, #512
	add	x26, x24, #3
	mov	x24, x25
	cmp	x26, x1
	b.lo	LBB0_40
	b	LBB0_87
LBB0_41:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x3, [x3, #96]
	cmp	x1, #2
	b.lo	LBB0_89
	and	x6, x1, #0xfffffffffffffffe
	mov	w7, #1
LBB0_43:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x3, [x0, #96]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x3, [x0, #352]
	add	x0, x0, #512
	add	x7, x7, #2
	cmp	x7, x1
	b.lo	LBB0_43
	b	LBB0_90
LBB0_44:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x21, [x3, #128]
	ldp	x22, x23, [x3, #144]
	ldp	x24, x3, [x3, #160]
	cmp	x1, #2
	b.lo	LBB0_92
	and	x25, x1, #0xfffffffffffffffe
	mov	w26, #1
LBB0_46:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x24, x3, [x0, #160]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x21, [x0, #384]
	stp	x22, x23, [x0, #400]
	stp	x24, x3, [x0, #416]
	add	x0, x0, #512
	add	x26, x26, #2
	cmp	x26, x1
	b.lo	LBB0_46
	b	LBB0_93
LBB0_47:
	ldp	x8, x9, [x3]
	ldp	x10, x11, [x3, #16]
	ldp	x12, x13, [x3, #32]
	ldp	x14, x15, [x3, #48]
	ldp	x16, x17, [x3, #64]
	ldp	x2, x4, [x3, #80]
	ldp	x5, x6, [x3, #96]
	ldp	x7, x19, [x3, #112]
	ldp	x20, x3, [x3, #128]
	cmp	x1, #2
	b.lo	LBB0_95
	and	x21, x1, #0xfffffffffffffffe
	mov	w22, #1
LBB0_49:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x3, [x0, #128]
	stp	x8, x9, [x0, #256]
	stp	x10, x11, [x0, #272]
	stp	x12, x13, [x0, #288]
	stp	x14, x15, [x0, #304]
	stp	x16, x17, [x0, #320]
	stp	x2, x4, [x0, #336]
	stp	x5, x6, [x0, #352]
	stp	x7, x19, [x0, #368]
	stp	x20, x3, [x0, #384]
	add	x0, x0, #512
	add	x22, x22, #2
	cmp	x22, x1
	b.lo	LBB0_49
	b	LBB0_96
LBB0_50:
	mov	x25, #0
LBB0_51:
	cmp	x25, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	ldp	x8, x9, [x3, #160]
	stp	x8, x9, [x0, #160]
	ldp	x8, x9, [x3, #176]
	stp	x8, x9, [x0, #176]
	ldp	x8, x9, [x3, #192]
	stp	x8, x9, [x0, #192]
	ldp	x8, x9, [x3, #208]
	stp	x8, x9, [x0, #208]
	ldp	x8, x9, [x3, #224]
	stp	x8, x9, [x0, #224]
	b	LBB0_98
LBB0_53:
	mov	x27, #0
LBB0_54:
	cmp	x27, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x24, x25, [x0, #160]
	stp	x26, x3, [x0, #176]
	b	LBB0_98
LBB0_56:
	mov	x2, #0
LBB0_57:
	cmp	x2, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	b	LBB0_98
LBB0_59:
	mov	x25, #0
LBB0_60:
	cmp	x25, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	ldp	x8, x9, [x3, #160]
	stp	x8, x9, [x0, #160]
	ldp	x8, x9, [x3, #176]
	stp	x8, x9, [x0, #176]
	ldp	x8, x9, [x3, #192]
	stp	x8, x9, [x0, #192]
	b	LBB0_98
LBB0_62:
	mov	x23, #0
LBB0_63:
	cmp	x23, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x3, [x0, #144]
	b	LBB0_98
LBB0_65:
	mov	x14, #0
LBB0_66:
	cmp	x14, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	b	LBB0_98
LBB0_68:
	mov	x25, #0
LBB0_69:
	cmp	x25, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	ldp	x8, x9, [x3, #160]
	stp	x8, x9, [x0, #160]
	ldp	x8, x9, [x3, #176]
	stp	x8, x9, [x0, #176]
	ldp	x8, x9, [x3, #192]
	stp	x8, x9, [x0, #192]
	ldp	x8, x9, [x3, #208]
	stp	x8, x9, [x0, #208]
	ldp	x8, x9, [x3, #224]
	stp	x8, x9, [x0, #224]
	ldp	x8, x9, [x3, #240]
	stp	x8, x9, [x0, #240]
	b	LBB0_98
LBB0_71:
	mov	x16, #0
LBB0_72:
	cmp	x16, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	b	LBB0_98
LBB0_74:
	mov	x19, #0
LBB0_75:
	cmp	x19, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x3, [x0, #112]
	b	LBB0_98
LBB0_77:
	mov	x10, #0
LBB0_78:
	cmp	x10, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	b	LBB0_98
LBB0_80:
	mov	x12, #0
LBB0_81:
	cmp	x12, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	b	LBB0_98
LBB0_83:
	mov	x4, #0
LBB0_84:
	cmp	x4, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x3, [x0, #80]
	b	LBB0_98
LBB0_86:
	mov	x25, #0
LBB0_87:
	cmp	x25, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	ldp	x8, x9, [x3, #160]
	stp	x8, x9, [x0, #160]
	ldp	x8, x9, [x3, #176]
	stp	x8, x9, [x0, #176]
	ldp	x8, x9, [x3, #192]
	stp	x8, x9, [x0, #192]
	ldp	x8, x9, [x3, #208]
	stp	x8, x9, [x0, #208]
	b	LBB0_98
LBB0_89:
	mov	x6, #0
LBB0_90:
	cmp	x6, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x3, [x0, #96]
	b	LBB0_98
LBB0_92:
	mov	x25, #0
LBB0_93:
	cmp	x25, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x21, [x0, #128]
	stp	x22, x23, [x0, #144]
	stp	x24, x3, [x0, #160]
	b	LBB0_98
LBB0_95:
	mov	x21, #0
LBB0_96:
	cmp	x21, x1
	b.hs	LBB0_98
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x2, x4, [x0, #80]
	stp	x5, x6, [x0, #96]
	stp	x7, x19, [x0, #112]
	stp	x20, x3, [x0, #128]
LBB0_98:
	ldp	x20, x19, [sp, #64]
	ldp	x22, x21, [sp, #48]
	ldp	x24, x23, [sp, #32]
	ldp	x26, x25, [sp, #16]
	ldp	x28, x27, [sp], #80
LBB0_99:
	ret
	.loh AdrpAdd	Lloh0, Lloh1
	.section	__TEXT,__const
	.p2align	1, 0x0
lJTI0_0:
	.short	(LBB0_2-LBB0_2)>>2
	.short	(LBB0_32-LBB0_2)>>2
	.short	(LBB0_20-LBB0_2)>>2
	.short	(LBB0_26-LBB0_2)>>2
	.short	(LBB0_11-LBB0_2)>>2
	.short	(LBB0_35-LBB0_2)>>2
	.short	(LBB0_41-LBB0_2)>>2
	.short	(LBB0_29-LBB0_2)>>2
	.short	(LBB0_47-LBB0_2)>>2
	.short	(LBB0_17-LBB0_2)>>2
	.short	(LBB0_44-LBB0_2)>>2
	.short	(LBB0_8-LBB0_2)>>2
	.short	(LBB0_14-LBB0_2)>>2
	.short	(LBB0_38-LBB0_2)>>2
	.short	(LBB0_5-LBB0_2)>>2
	.short	(LBB0_23-LBB0_2)>>2
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_slab_scan_g_active_16
	.p2align	2
_slab_scan_g_active_16:
	stp	x28, x27, [sp, #-96]!
	stp	x26, x25, [sp, #16]
	stp	x24, x23, [sp, #32]
	stp	x22, x21, [sp, #48]
	stp	x20, x19, [sp, #64]
	stp	x29, x30, [sp, #80]
	add	x29, sp, #80
	add	x13, x2, #4
	add	x9, x2, #19
	lsr	x8, x9, #4
	and	x14, x9, #0xfffffffffffffff0
	subs	x9, x14, #16
	cset	w12, ne
	cmp	x9, x13
	cset	w9, lo
	tst	w12, w9
	mov	w9, #128
	csel	x9, x9, xzr, ne
	sub	x10, x14, #15
	cmp	x10, x13
	cset	w10, lo
	tst	w12, w10
	mov	w10, #32768
	csel	x10, x10, xzr, ne
	sub	x11, x14, #14
	cmp	x11, x13
	cset	w11, lo
	tst	w12, w11
	mov	w11, #8388608
	csel	x11, x11, xzr, ne
	sub	x15, x14, #13
	cmp	x15, x13
	cset	w15, lo
	tst	w12, w15
	mov	w12, #-2147483648
	csel	x12, x12, xzr, ne
	sub	x15, x14, #12
	cmp	x15, x13
	cset	w15, lo
	sub	x16, x14, #11
	cmp	x16, x13
	cset	w16, lo
	sub	x17, x14, #10
	cmp	x17, x13
	cset	w17, lo
	sub	x4, x14, #9
	cmp	x4, x13
	cset	w4, lo
	cbz	x14, LBB1_2
	sub	x5, x14, #8
	cmp	x5, x13
	cset	w5, lo
	ubfiz	x5, x5, #7, #32
	sub	x6, x14, #7
	cmp	x6, x13
	cset	w6, lo
	sub	x7, x14, #6
	cmp	x7, x13
	cset	w7, lo
	sub	x19, x14, #5
	cmp	x19, x13
	cset	w19, lo
	sub	x20, x14, #4
	cmp	x20, x13
	cset	w20, lo
	sub	x21, x14, #3
	cmp	x21, x13
	cset	w21, lo
	sub	x22, x14, #2
	cmp	x22, x13
	cset	w22, lo
	sub	x14, x14, #1
	cmp	x14, x13
	cset	w13, lo
	orr	x14, x5, x7, lsl #23
	orr	x14, x14, x6, lsl #15
	orr	x14, x14, x19, lsl #31
	orr	x14, x14, x20, lsl #39
	orr	x14, x14, x21, lsl #47
	orr	x14, x14, x22, lsl #55
	orr	x23, x14, x13, lsl #63
	sub	x8, x8, #1
	cmp	x8, #15
	b.ls	LBB1_3
	b	LBB1_84
LBB1_2:
	mov	x23, #0
	sub	x8, x8, #1
	cmp	x8, #15
	b.hi	LBB1_84
LBB1_3:
	lsl	x13, x15, #39
	orr	x13, x13, x16, lsl #47
	orr	x9, x9, x11
	orr	x11, x13, x17, lsl #55
	orr	x9, x9, x10
	orr	x10, x11, x4, lsl #63
	orr	x9, x9, x12
	orr	x24, x10, x9
Lloh2:
	adrp	x9, lJTI1_0@PAGE
Lloh3:
	add	x9, x9, lJTI1_0@PAGEOFF
	adr	x10, LBB1_4
	ldrh	w11, [x9, x8, lsl #1]
	add	x10, x10, x11, lsl #2
	br	x10
LBB1_4:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
LBB1_6:
	ldp	x8, x9, [x0]
	and	x8, x8, x24
	and	x9, x9, x23
	orr	x8, x9, x8
	cbnz	x8, LBB1_8
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_6
	b	LBB1_85
LBB1_8:
	mov	x20, x0
	mov	x25, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x25, #1
	b.ne	LBB1_6
	b	LBB1_85
LBB1_9:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_11:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #32]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #48]
	orr	x8, x8, x9
	orr	x9, x10, x11
	ldp	x10, x11, [x0, #64]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #80]
	orr	x8, x8, x9
	orr	x9, x10, x11
	ldp	x10, x11, [x0, #96]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #112]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #128]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #144]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #160]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #176]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #192]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #208]
	orr	x9, x9, x10
	orr	x9, x9, x11
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #224]
	and	x9, x9, x24
	and	x10, x10, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_13
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_11
	b	LBB1_85
LBB1_13:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_11
	b	LBB1_85
LBB1_14:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_16:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #176]
	and	x7, x7, x24
	and	x20, x20, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #96]
	orr	x11, x11, x12
	orr	x11, x11, x13
	orr	x8, x8, x11
	orr	x11, x14, x15
	ldp	x12, x13, [x0, #112]
	orr	x11, x11, x16
	orr	x11, x11, x17
	orr	x8, x8, x11
	orr	x11, x5, x6
	ldp	x14, x15, [x0, #128]
	orr	x9, x11, x9
	orr	x9, x9, x10
	orr	x9, x9, x12
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #144]
	orr	x11, x13, x14
	orr	x11, x11, x15
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #160]
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x7
	orr	x9, x9, x20
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_18
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_16
	b	LBB1_85
LBB1_18:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_16
	b	LBB1_85
LBB1_19:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_21:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	and	x16, x16, x24
	and	x17, x17, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x9, x9, x16
	orr	x9, x9, x17
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_23
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_21
	b	LBB1_85
LBB1_23:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_21
	b	LBB1_85
LBB1_24:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_26:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #192]
	and	x5, x5, x24
	and	x6, x6, x23
	orr	x9, x9, x10
	ldp	x10, x7, [x0, #80]
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	ldp	x9, x11, [x0, #96]
	orr	x12, x14, x15
	orr	x12, x12, x16
	orr	x12, x12, x17
	orr	x8, x8, x12
	ldp	x12, x13, [x0, #112]
	orr	x10, x10, x7
	orr	x9, x10, x9
	orr	x9, x9, x11
	orr	x9, x9, x12
	ldp	x10, x11, [x0, #128]
	orr	x8, x8, x9
	orr	x9, x13, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #144]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #160]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #176]
	orr	x9, x11, x9
	orr	x9, x9, x10
	orr	x9, x9, x5
	orr	x9, x9, x6
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_28
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_26
	b	LBB1_85
LBB1_28:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_26
	b	LBB1_85
LBB1_29:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_31:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #144]
	and	x7, x7, x24
	and	x20, x20, x23
	ldp	x21, x22, [x0, #96]
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	ldp	x10, x11, [x0, #112]
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x9, x9, x16
	orr	x9, x9, x17
	ldp	x12, x13, [x0, #128]
	orr	x8, x8, x9
	orr	x9, x5, x6
	orr	x9, x9, x21
	orr	x9, x9, x22
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x9, x9, x7
	orr	x9, x9, x20
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_33
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_31
	b	LBB1_85
LBB1_33:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_31
	b	LBB1_85
LBB1_34:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_36:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	and	x12, x12, x24
	and	x13, x13, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x8, x8, x9
	orr	x8, x8, x13
	tst	x8, #0x8080808080808080
	b.ne	LBB1_38
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_36
	b	LBB1_85
LBB1_38:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_36
	b	LBB1_85
LBB1_39:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_41:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #32]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #48]
	orr	x8, x8, x9
	orr	x9, x10, x11
	ldp	x10, x11, [x0, #64]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #80]
	orr	x8, x8, x9
	orr	x9, x10, x11
	ldp	x10, x11, [x0, #96]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #112]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #128]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #144]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #160]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #176]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #192]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #208]
	orr	x9, x9, x10
	orr	x9, x9, x11
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #224]
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #240]
	and	x10, x10, x24
	orr	x9, x9, x10
	and	x10, x11, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_43
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_41
	b	LBB1_85
LBB1_43:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_41
	b	LBB1_85
LBB1_44:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_46:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	and	x14, x14, x24
	and	x15, x15, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_48
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_46
	b	LBB1_85
LBB1_48:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_46
	b	LBB1_85
LBB1_49:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_51:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #96]
	ldp	x21, x22, [x0, #112]
	and	x21, x21, x24
	and	x22, x22, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x9, x9, x16
	orr	x9, x9, x17
	orr	x8, x8, x9
	orr	x9, x5, x6
	orr	x9, x9, x7
	orr	x9, x9, x20
	orr	x9, x9, x21
	orr	x8, x8, x9
	orr	x8, x8, x22
	tst	x8, #0x8080808080808080
	b.ne	LBB1_53
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_51
	b	LBB1_85
LBB1_53:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_51
	b	LBB1_85
LBB1_54:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_56:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #16]
	and	x9, x9, x24
	and	x10, x10, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_58
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_56
	b	LBB1_85
LBB1_58:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_56
	b	LBB1_85
LBB1_59:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_61:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	and	x5, x5, x24
	and	x6, x6, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x9, x9, x16
	orr	x9, x9, x17
	orr	x8, x8, x9
	orr	x9, x5, x6
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_63
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_61
	b	LBB1_85
LBB1_63:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_61
	b	LBB1_85
LBB1_64:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_66:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #32]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #48]
	orr	x8, x8, x9
	orr	x9, x10, x11
	ldp	x10, x11, [x0, #64]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #80]
	orr	x8, x8, x9
	orr	x9, x10, x11
	ldp	x10, x11, [x0, #96]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #112]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #128]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #144]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #160]
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #176]
	orr	x9, x11, x9
	orr	x9, x9, x10
	ldp	x10, x11, [x0, #192]
	orr	x9, x9, x10
	orr	x9, x9, x11
	ldp	x10, x11, [x0, #208]
	and	x10, x10, x24
	orr	x9, x9, x10
	and	x10, x11, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_68
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_66
	b	LBB1_85
LBB1_68:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_66
	b	LBB1_85
LBB1_69:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_71:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #96]
	and	x7, x7, x24
	and	x20, x20, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x9, x9, x16
	orr	x9, x9, x17
	orr	x8, x8, x9
	orr	x9, x5, x6
	orr	x9, x9, x7
	orr	x9, x9, x20
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_73
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_71
	b	LBB1_85
LBB1_73:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_71
	b	LBB1_85
LBB1_74:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_76:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #96]
	ldp	x21, x22, [x0, #160]
	and	x21, x21, x24
	and	x22, x22, x23
	orr	x9, x9, x10
	orr	x8, x8, x9
	ldp	x9, x10, [x0, #112]
	orr	x11, x11, x12
	orr	x11, x11, x13
	orr	x8, x8, x11
	orr	x11, x14, x15
	ldp	x12, x13, [x0, #128]
	orr	x11, x11, x16
	orr	x11, x11, x17
	orr	x8, x8, x11
	orr	x11, x5, x6
	ldp	x14, x15, [x0, #144]
	orr	x11, x11, x7
	orr	x11, x11, x20
	orr	x9, x11, x9
	orr	x8, x8, x9
	orr	x9, x10, x12
	orr	x9, x9, x13
	orr	x9, x9, x14
	orr	x9, x9, x15
	orr	x9, x9, x21
	orr	x8, x8, x9
	orr	x8, x8, x22
	tst	x8, #0x8080808080808080
	b.ne	LBB1_78
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_76
	b	LBB1_85
LBB1_78:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_76
	b	LBB1_85
LBB1_79:
	cbz	x1, LBB1_84
	mov	x19, #0
	mov	x4, #0
	mov	x25, #141287244169216
	movk	x25, #32896, lsl #48
LBB1_81:
	ldp	x8, x9, [x0]
	and	x8, x8, x25
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x5, x6, [x0, #80]
	ldp	x7, x20, [x0, #96]
	ldp	x21, x22, [x0, #128]
	and	x21, x21, x24
	and	x22, x22, x23
	ldp	x26, x27, [x0, #112]
	orr	x9, x9, x10
	orr	x8, x8, x9
	orr	x9, x11, x12
	orr	x9, x9, x13
	orr	x8, x8, x9
	orr	x9, x14, x15
	orr	x9, x9, x16
	orr	x9, x9, x17
	orr	x8, x8, x9
	orr	x9, x5, x6
	orr	x9, x9, x7
	orr	x9, x9, x20
	orr	x9, x9, x26
	orr	x8, x8, x9
	orr	x9, x27, x21
	orr	x9, x9, x22
	orr	x8, x8, x9
	tst	x8, #0x8080808080808080
	b.ne	LBB1_83
	add	x0, x0, #256
	add	x19, x19, #256
	subs	x1, x1, #1
	b.ne	LBB1_81
	b	LBB1_85
LBB1_83:
	mov	x20, x0
	mov	x26, x1
	mov	x1, x19
	mov	x21, x2
	mov	x22, x3
	bl	_extract_row_survivors
	mov	x2, x21
	mov	x3, x22
	mov	x4, x0
	add	x0, x20, #256
	add	x19, x19, #256
	subs	x1, x26, #1
	b.ne	LBB1_81
	b	LBB1_85
LBB1_84:
	mov	x4, #0
LBB1_85:
	mov	x0, x4
	ldp	x29, x30, [sp, #80]
	ldp	x20, x19, [sp, #64]
	ldp	x22, x21, [sp, #48]
	ldp	x24, x23, [sp, #32]
	ldp	x26, x25, [sp, #16]
	ldp	x28, x27, [sp], #96
	ret
	.loh AdrpAdd	Lloh2, Lloh3
	.section	__TEXT,__const
	.p2align	1, 0x0
lJTI1_0:
	.short	(LBB1_4-LBB1_4)>>2
	.short	(LBB1_54-LBB1_4)>>2
	.short	(LBB1_34-LBB1_4)>>2
	.short	(LBB1_44-LBB1_4)>>2
	.short	(LBB1_19-LBB1_4)>>2
	.short	(LBB1_59-LBB1_4)>>2
	.short	(LBB1_69-LBB1_4)>>2
	.short	(LBB1_49-LBB1_4)>>2
	.short	(LBB1_79-LBB1_4)>>2
	.short	(LBB1_29-LBB1_4)>>2
	.short	(LBB1_74-LBB1_4)>>2
	.short	(LBB1_14-LBB1_4)>>2
	.short	(LBB1_24-LBB1_4)>>2
	.short	(LBB1_64-LBB1_4)>>2
	.short	(LBB1_9-LBB1_4)>>2
	.short	(LBB1_39-LBB1_4)>>2
	.section	__TEXT,__text,regular,pure_instructions
	.p2align	2
_extract_row_survivors:
	add	x9, x2, #4
	cmp	x9, #4
	b.ls	LBB2_13
	ldrsb	w8, [x0, #4]
	tbnz	w8, #31, LBB2_3
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
	add	x8, x1, #5
	str	x8, [x3, x4, lsl #3]
	add	x4, x4, #1
LBB2_6:
	cmp	x9, #6
	b.eq	LBB2_12
	ldrsb	w8, [x0, #6]
	tbz	w8, #31, LBB2_9
	add	x8, x1, #6
	str	x8, [x3, x4, lsl #3]
	add	x4, x4, #1
LBB2_9:
	cmp	x9, #7
	b.eq	LBB2_12
	ldrsb	w8, [x0, #7]
	tbz	w8, #31, LBB2_12
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
	lsr	x10, x9, #3
	mov	w11, #1
	b	LBB2_18
LBB2_16:
	mov	x8, x4
LBB2_17:
	add	x11, x11, #1
	mov	x4, x8
	cmp	x11, x10
	b.eq	LBB2_22
LBB2_18:
	ldr	x8, [x0, x11, lsl #3]
	ands	x12, x8, #0x8080808080808080
	b.eq	LBB2_16
	add	x13, x1, x11, lsl #3
LBB2_20:
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
LBB2_24:
	add	x10, x10, #1
	cmp	x10, x9
	b.hs	LBB2_23
LBB2_25:
	ldrsb	w11, [x0, x10]
	tbz	w11, #31, LBB2_24
	add	x11, x1, x10
	str	x11, [x3, x8, lsl #3]
	add	x8, x8, #1
	b	LBB2_24
	.globl	_slab_filter_f_at_gsurv
	.p2align	2
_slab_filter_f_at_gsurv:
	mov	x8, #0
	cbnz	x2, LBB3_3
LBB3_1:
	mov	x0, x8
	ret
LBB3_2:
	subs	x2, x2, #1
	b.eq	LBB3_1
LBB3_3:
	ldr	x9, [x1], #8
	ldrsb	w10, [x0, x9]
	tbz	w10, #31, LBB3_2
	str	x9, [x3, x8, lsl #3]
	add	x8, x8, #1
	b	LBB3_2
.subsections_via_symbols
```

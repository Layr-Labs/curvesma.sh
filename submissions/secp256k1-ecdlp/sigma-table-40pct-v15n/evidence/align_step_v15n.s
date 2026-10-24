	.build_version macos, 27, 0	sdk_version 27, 0
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_class_slot                     ; -- Begin function class_slot
	.p2align	2
_class_slot:                            ; @class_slot
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #128
	stp	x22, x21, [sp, #80]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #96]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #112]            ; 16-byte Folded Spill
	add	x29, sp, #112
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	mov	x19, x2
	mov	x20, x1
	mov	x21, x0
Lloh0:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh1:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh2:
	ldr	x8, [x8]
	stur	x8, [x29, #-40]
Lloh3:
	adrp	x22, _BETA@GOTPAGE
Lloh4:
	ldr	x22, [x22, _BETA@GOTPAGEOFF]
	add	x0, sp, #40
	mov	x1, x21
	mov	x2, x22
	bl	_fmul
	add	x0, sp, #8
	add	x1, sp, #40
	mov	x2, x22
	bl	_fmul
	ldr	x9, [x21, #24]
	ldr	x13, [sp, #64]
	cmp	x13, x9
	b.hs	LBB0_3
; %bb.1:
	ldr	x14, [sp, #56]
LBB0_2:
	ldp	x8, x10, [sp, #40]
	mov	w12, #1                         ; =0x1
	mov	x9, x13
	mov	x11, x14
	b	LBB0_5
LBB0_3:
	ldp	x8, x10, [x21]
	ldr	x11, [x21, #16]
	b.ls	LBB0_11
LBB0_4:
	mov	x12, #0                         ; =0x0
LBB0_5:
	ldr	x13, [sp, #32]
	cmp	x13, x9
	b.hs	LBB0_8
; %bb.6:
	ldr	x14, [sp, #24]
LBB0_7:
	ldp	x8, x10, [sp, #8]
	mov	w12, #2                         ; =0x2
	mov	x9, x13
	mov	x11, x14
	b	LBB0_9
LBB0_8:
	b.ls	LBB0_16
LBB0_9:
	stp	x8, x10, [x20]
	stp	x11, x9, [x20, #16]
	str	x12, [x19]
	ldur	x8, [x29, #-40]
Lloh5:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh6:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh7:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB0_21
; %bb.10:
	ldp	x29, x30, [sp, #112]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #96]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #128
	ret
LBB0_11:
	ldr	x14, [sp, #56]
	cmp	x14, x11
	b.lo	LBB0_2
; %bb.12:
	b.hi	LBB0_4
; %bb.13:
	ldr	x15, [sp, #48]
	cmp	x15, x10
	b.lo	LBB0_2
; %bb.14:
	mov	x12, #0                         ; =0x0
	cmp	x15, x10
	b.hi	LBB0_5
; %bb.15:
	ldr	x15, [sp, #40]
	cmp	x15, x8
	b.lo	LBB0_2
	b	LBB0_5
LBB0_16:
	ldr	x14, [sp, #24]
	cmp	x14, x11
	b.lo	LBB0_7
; %bb.17:
	b.hi	LBB0_9
; %bb.18:
	ldr	x15, [sp, #16]
	cmp	x15, x10
	b.lo	LBB0_7
; %bb.19:
	b.hi	LBB0_9
; %bb.20:
	ldr	x15, [sp, #8]
	cmp	x15, x8
	b.lo	LBB0_7
	b	LBB0_9
LBB0_21:
	bl	___stack_chk_fail
	.loh AdrpLdrGot	Lloh3, Lloh4
	.loh AdrpLdrGotLdr	Lloh0, Lloh1, Lloh2
	.loh AdrpLdrGotLdr	Lloh5, Lloh6, Lloh7
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function fmul
_fmul:                                  ; @fmul
	.cfi_startproc
; %bb.0:
	ldp	x10, x9, [x2]
	ldp	x11, x14, [x1]
	mul	x8, x10, x11
	umulh	x12, x10, x11
	umulh	x13, x9, x11
	mul	x15, x9, x11
	adds	x15, x12, x15
	cinc	x16, x13, hs
	ldp	x13, x12, [x2, #16]
	umulh	x17, x13, x11
	mul	x2, x13, x11
	adds	x16, x16, x2
	cinc	x17, x17, hs
	umulh	x2, x12, x11
	mul	x11, x12, x11
	adds	x17, x17, x11
	cinc	x2, x2, hs
	umulh	x3, x10, x14
	mul	x11, x10, x14
	adds	x11, x15, x11
	cinc	x15, x3, hs
	umulh	x3, x9, x14
	mul	x4, x9, x14
	adds	x15, x16, x15
	cset	w16, hs
	adds	x15, x15, x4
	adc	x16, x16, x3
	umulh	x3, x13, x14
	mul	x4, x13, x14
	adds	x16, x17, x16
	cset	w17, hs
	adds	x16, x16, x4
	adc	x17, x17, x3
	umulh	x3, x12, x14
	mul	x14, x12, x14
	adds	x17, x2, x17
	cset	w2, hs
	adds	x14, x17, x14
	adc	x17, x2, x3
	ldp	x2, x1, [x1, #16]
	umulh	x3, x10, x2
	mul	x4, x10, x2
	adds	x15, x15, x4
	cinc	x3, x3, hs
	umulh	x4, x9, x2
	mul	x5, x9, x2
	adds	x16, x16, x3
	cset	w3, hs
	adds	x16, x16, x5
	adc	x3, x3, x4
	umulh	x4, x13, x2
	mul	x5, x13, x2
	adds	x14, x14, x3
	cset	w3, hs
	adds	x14, x14, x5
	adc	x3, x3, x4
	umulh	x4, x12, x2
	mul	x2, x12, x2
	adds	x17, x17, x3
	cset	w3, hs
	adds	x17, x17, x2
	adc	x2, x3, x4
	umulh	x3, x10, x1
	mul	x10, x10, x1
	adds	x10, x16, x10
	cinc	x16, x3, hs
	umulh	x3, x9, x1
	mul	x9, x9, x1
	adds	x14, x14, x16
	cset	w16, hs
	adds	x14, x14, x9
	adc	x9, x16, x3
	umulh	x16, x13, x1
	mul	x13, x13, x1
	adds	x9, x17, x9
	cset	w17, hs
	adds	x13, x9, x13
	adc	x9, x17, x16
	umulh	x16, x12, x1
	mul	x12, x12, x1
	adds	x9, x2, x9
	cset	w17, hs
	adds	x12, x9, x12
	adc	x16, x17, x16
	mov	x9, #977                        ; =0x3d1
	movk	x9, #1, lsl #32
	umulh	x17, x16, x9
	mul	x16, x16, x9
	adds	x10, x16, x10
	cinc	x16, x17, hs
	umulh	x17, x12, x9
	mul	x12, x12, x9
	adds	x12, x12, x15
	adcs	x15, x10, x17
	cset	w17, hs
	umulh	x10, x13, x9
	mul	x13, x13, x9
	adds	x11, x13, x11
	adcs	x10, x12, x10
	cset	w13, hs
	umulh	x1, x14, x9
	mul	x12, x14, x9
	adds	x12, x12, x8
	adcs	x11, x11, x1
	adcs	x10, x10, xzr
	adcs	x8, x13, x15
	adc	x13, x17, x16
	cbz	x13, LBB1_2
; %bb.1:
	mul	x14, x13, x9
	adds	x12, x14, x12
	cset	w14, hs
	adcs	xzr, x11, xzr
	adcs	xzr, x10, xzr
	adcs	xzr, x8, xzr
	umulh	x13, x13, x9
	adcs	x12, x12, x13
	adcs	x11, x14, x11
	adcs	x10, x10, xzr
	cinc	x8, x8, hs
LBB1_2:
	stp	x12, x11, [x0]
	stp	x10, x8, [x0, #16]
	adds	x9, x12, x9
	mov	x12, #-1                        ; =0xffffffffffffffff
	adc	x13, xzr, x12
	asr	x14, x13, #63
	adds	x11, x11, x13
	cinc	x13, x14, hs
	adds	x11, x11, #1
	adc	x13, x13, x12
	asr	x14, x13, #63
	adds	x10, x10, x13
	cinc	x13, x14, hs
	adds	x10, x10, #1
	adc	x13, x13, x12
	asr	x14, x13, #63
	adds	x8, x8, x13
	cinc	x13, x14, hs
	adds	x8, x8, #1
	adc	x12, x13, x12
	cbnz	x12, LBB1_4
; %bb.3:
	stp	x9, x11, [x0]
	stp	x10, x8, [x0, #16]
LBB1_4:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_even_y                         ; -- Begin function even_y
	.p2align	2
_even_y:                                ; @even_y
	.cfi_startproc
; %bb.0:
	ldr	x8, [x0]
	tbz	w8, #0, LBB2_2
; %bb.1:
	ldp	x9, x10, [x0, #8]
	ldr	x11, [x0, #24]
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	neg	x8, x8
	adds	x8, x8, x12
	mvn	x9, x9
	adcs	x9, x9, xzr
	cset	w12, hs
	subs	x9, x9, #1
	mvn	x10, x10
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	cinc	x12, x12, hs
	sub	x11, x12, x11
	stp	x8, x9, [x0]
	sub	x8, x11, #2
	stp	x10, x8, [x0, #16]
LBB2_2:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_align_sweep                    ; -- Begin function align_sweep
	.p2align	2
_align_sweep:                           ; @align_sweep
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #96
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
Lloh8:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh9:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh10:
	ldr	x8, [x8]
	str	x8, [sp, #40]
Lloh11:
	adrp	x8, _NREC@GOTPAGE
Lloh12:
	ldr	x8, [x8, _NREC@GOTPAGEOFF]
Lloh13:
	ldr	x20, [x8]
	cbz	x20, LBB3_5
; %bb.1:
Lloh14:
	adrp	x19, _A@GOTPAGE
Lloh15:
	ldr	x19, [x19, _A@GOTPAGEOFF]
	mov	x21, #-977                      ; =0xfffffffffffffc2f
	movk	x21, #65534, lsl #32
	b	LBB3_3
LBB3_2:                                 ;   in Loop: Header=BB3_3 Depth=1
	ldur	q0, [sp, #8]
	ldur	q1, [sp, #24]
	stp	q0, q1, [x19]
	ldr	x9, [x19, #72]
	add	x9, x9, #1
	stp	x8, x9, [x19, #64]
	add	x19, x19, #96
	subs	x20, x20, #1
	b.eq	LBB3_5
LBB3_3:                                 ; =>This Inner Loop Header: Depth=1
	add	x1, sp, #8
	mov	x2, sp
	mov	x0, x19
	bl	_class_slot
	ldr	x8, [sp, #8]
	and	x8, x8, #0x7ff
	add	x9, x8, x8, lsl #1
	ldr	x8, [sp]
	add	x9, x9, x8
	str	x9, [x19, #88]
	ldr	x9, [x19, #32]
	tbz	w9, #0, LBB3_2
; %bb.4:                                ;   in Loop: Header=BB3_3 Depth=1
	ldp	x10, x11, [x19, #40]
	ldr	x12, [x19, #56]
	neg	x9, x9
	adds	x9, x9, x21
	mvn	x10, x10
	adcs	x10, x10, xzr
	cset	w13, hs
	subs	x10, x10, #1
	mvn	x11, x11
	adcs	x11, x11, x13
	cset	w13, hs
	subs	x11, x11, #1
	cinc	x13, x13, hs
	sub	x12, x13, x12
	stp	x9, x10, [x19, #32]
	sub	x9, x12, #2
	stp	x11, x9, [x19, #48]
	b	LBB3_2
LBB3_5:
	ldr	x8, [sp, #40]
Lloh16:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh17:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh18:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB3_7
; %bb.6:
	ldp	x29, x30, [sp, #80]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	add	sp, sp, #96
	ret
LBB3_7:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh11, Lloh12, Lloh13
	.loh AdrpLdrGotLdr	Lloh8, Lloh9, Lloh10
	.loh AdrpLdrGot	Lloh14, Lloh15
	.loh AdrpLdrGotLdr	Lloh16, Lloh17, Lloh18
	.cfi_endproc
                                        ; -- End function
	.globl	_report_proc                    ; -- Begin function report_proc
	.p2align	2
_report_proc:                           ; @report_proc
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #112
	stp	x24, x23, [sp, #48]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #64]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #80]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #96]             ; 16-byte Folded Spill
	add	x29, sp, #96
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	mov	x19, x0
Lloh19:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh20:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh21:
	ldr	x8, [x8]
	str	x8, [sp, #40]
Lloh22:
	adrp	x20, _REPORT_REC@GOTPAGE
Lloh23:
	ldr	x20, [x20, _REPORT_REC@GOTPAGEOFF]
	ldp	x21, x22, [x20]
	ldp	x23, x24, [x20, #16]
Lloh24:
	adrp	x2, _BETA@GOTPAGE
Lloh25:
	ldr	x2, [x2, _BETA@GOTPAGEOFF]
	add	x0, sp, #8
	mov	x1, x20
	bl	_fmul
	add	x0, sp, #8
	add	x1, sp, #8
	mov	x2, x20
	bl	_fmul
	ldrh	w8, [sp, #8]
	str	x8, [x19]
Lloh26:
	adrp	x9, _NODE_HEAD@GOTPAGE
Lloh27:
	ldr	x9, [x9, _NODE_HEAD@GOTPAGEOFF]
	add	x8, x9, x8, lsl #3
	ldr	x9, [x8, #24]
	cmp	x24, x9
	b.ne	LBB4_5
; %bb.1:
	ldr	x9, [x8, #16]
	cmp	x23, x9
	b.ne	LBB4_5
; %bb.2:
	ldr	x9, [x8, #8]
	cmp	x22, x9
	b.ne	LBB4_5
; %bb.3:
	ldr	x9, [x8]
	cmp	x21, x9
	b.ne	LBB4_5
; %bb.4:
	mov	w9, #1                          ; =0x1
	b	LBB4_6
LBB4_5:
	mov	x9, #0                          ; =0x0
	stp	x21, x22, [x8]
	stp	x23, x24, [x8, #16]
LBB4_6:
Lloh28:
	adrp	x8, _REPORT_STORE@GOTPAGE
Lloh29:
	ldr	x8, [x8, _REPORT_STORE@GOTPAGEOFF]
Lloh30:
	str	x9, [x8]
	ldr	x8, [sp, #40]
Lloh31:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh32:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh33:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB4_8
; %bb.7:
	ldp	x29, x30, [sp, #96]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #80]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #64]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #48]             ; 16-byte Folded Reload
	add	sp, sp, #112
	ret
LBB4_8:
	bl	___stack_chk_fail
	.loh AdrpLdrGot	Lloh26, Lloh27
	.loh AdrpLdrGot	Lloh24, Lloh25
	.loh AdrpLdrGot	Lloh22, Lloh23
	.loh AdrpLdrGotLdr	Lloh19, Lloh20, Lloh21
	.loh AdrpLdrGotLdr	Lloh31, Lloh32, Lloh33
	.loh AdrpLdrGotStr	Lloh28, Lloh29, Lloh30
	.cfi_endproc
                                        ; -- End function
	.globl	_api_fmul                       ; -- Begin function api_fmul
	.p2align	2
_api_fmul:                              ; @api_fmul
	.cfi_startproc
; %bb.0:
	b	_fmul
	.cfi_endproc
                                        ; -- End function
	.globl	_api_fsub                       ; -- Begin function api_fsub
	.p2align	2
_api_fsub:                              ; @api_fsub
	.cfi_startproc
; %bb.0:
	ldp	x8, x9, [x1]
	ldp	x10, x11, [x2]
	subs	x8, x8, x10
	ngc	x10, xzr
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	asr	x11, x10, #63
	ldp	x12, x13, [x1, #16]
	ldp	x14, x15, [x2, #16]
	subs	x12, x12, x14
	ngc	x14, xzr
	adds	x10, x10, x12
	adc	x11, x11, x14
	asr	x12, x11, #63
	subs	x13, x13, x15
	ngc	x14, xzr
	adds	x11, x11, x13
	adc	x12, x12, x14
	mov	x13, #-977                      ; =0xfffffffffffffc2f
	movk	x13, #65534, lsl #32
	adds	x13, x8, x13
	adcs	x14, x9, xzr
	cset	w15, hs
	subs	x14, x14, #1
	adcs	x15, x10, x15
	cset	w16, hs
	subs	x15, x15, #1
	adc	x16, x11, x16
	sub	x16, x16, #1
	cmp	x12, #0
	csel	x8, x8, x13, eq
	csel	x9, x9, x14, eq
	csel	x10, x10, x15, eq
	stp	x8, x9, [x0]
	csel	x8, x11, x16, eq
	stp	x10, x8, [x0, #16]
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_api_class                      ; -- Begin function api_class
	.p2align	2
_api_class:                             ; @api_class
	.cfi_startproc
; %bb.0:
	b	_class_slot
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols

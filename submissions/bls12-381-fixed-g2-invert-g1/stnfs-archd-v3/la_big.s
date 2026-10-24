	.build_version macos, 26, 0	sdk_version 27, 0
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_main                           ; -- Begin function main
	.p2align	2
_main:                                  ; @main
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #272
	stp	x28, x27, [sp, #176]            ; 16-byte Folded Spill
	stp	x26, x25, [sp, #192]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #208]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #224]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #240]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #256]            ; 16-byte Folded Spill
	add	x29, sp, #256
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
	cmp	w0, #1
	b.gt	LBB0_2
; %bb.1:
Lloh0:
	adrp	x8, ___stderrp@GOTPAGE
Lloh1:
	ldr	x8, [x8, ___stderrp@GOTPAGEOFF]
Lloh2:
	ldr	x3, [x8]
Lloh3:
	adrp	x0, l_.str@PAGE
Lloh4:
	add	x0, x0, l_.str@PAGEOFF
	mov	w19, #1                         ; =0x1
	mov	w1, #47                         ; =0x2f
	mov	w2, #1                          ; =0x1
	bl	_fwrite
	b	LBB0_60
LBB0_2:
	mov	x20, x1
	ldr	x0, [x1, #8]
Lloh5:
	adrp	x1, l_.str.1@PAGE
Lloh6:
	add	x1, x1, l_.str.1@PAGEOFF
	bl	_strcmp
	mov	x23, x0
	cbz	w0, LBB0_4
; %bb.3:
	ldr	x0, [x20, #16]
	bl	_atoi
	mov	x26, x0
	ldr	x0, [x20, #24]
	bl	_atoi
	mov	x21, x0
	ldr	x0, [x20, #32]
	bl	_atoi
	mov	x22, x0
	b	LBB0_5
LBB0_4:
	mov	w22, #1                         ; =0x1
	mov	w26, #64                        ; =0x40
	mov	w21, #48                        ; =0x30
LBB0_5:
	sxtw	x19, w26
	sbfiz	x25, x26, #5, #32
	sxtw	x20, w21
	mul	x0, x25, x20
	bl	_malloc
	mov	x24, x0
	sbfiz	x0, x21, #5, #32
	bl	_malloc
	mov	x21, x0
	add	x8, x19, w26, sxtw #2
	lsl	x0, x8, #4
	str	x0, [sp, #104]                  ; 8-byte Folded Spill
	bl	_malloc
	mov	x19, x0
	mov	x0, x25
	bl	_malloc
	cmp	x24, #0
	stur	x21, [x29, #-120]               ; 8-byte Folded Spill
	ccmp	x21, #0, #4, ne
	stp	x0, x19, [x29, #-104]           ; 16-byte Folded Spill
	ccmp	x19, #0, #4, ne
	ccmp	x0, #0, #4, ne
	b.ne	LBB0_7
; %bb.6:
	mov	w19, #1                         ; =0x1
	b	LBB0_60
LBB0_7:
	mov	x8, #32072                      ; =0x7d48
	movk	x8, #10653, lsl #16
	movk	x8, #42835, lsl #32
	movk	x8, #29677, lsl #48
	mov	x9, #23550                      ; =0x5bfe
	movk	x9, #65534, lsl #16
	movk	x9, #41986, lsl #32
	movk	x9, #21437, lsl #48
	mov	x10, #55301                     ; =0xd805
	movk	x10, #2465, lsl #16
	movk	x10, #55304, lsl #32
	movk	x10, #13113, lsl #48
	mul	w12, w26, w20
	adrp	x11, _rngs@PAGE
	cmp	w12, #1
	b.lt	LBB0_17
; %bb.8:
	mov	x13, #0                         ; =0x0
	ldr	x14, [x11, _rngs@PAGEOFF]
	b	LBB0_10
LBB0_9:                                 ;   in Loop: Header=BB0_10 Depth=1
	stp	x16, x17, [x15]
	stp	x0, x1, [x15, #16]
	add	x13, x13, #1
	cmp	x13, x12
	b.eq	LBB0_16
LBB0_10:                                ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_12 Depth 2
	add	x15, x24, x13, lsl #5
	b	LBB0_12
LBB0_11:                                ;   in Loop: Header=BB0_12 Depth=2
	cmp	x2, x3
	b.ls	LBB0_9
LBB0_12:                                ;   Parent Loop BB0_10 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	eor	x14, x14, x14, lsl #13
	eor	x14, x14, x14, lsr #7
	eor	x16, x14, x14, lsl #17
	eor	x14, x16, x16, lsl #13
	eor	x14, x14, x14, lsr #7
	eor	x17, x14, x14, lsl #17
	eor	x14, x17, x17, lsl #13
	eor	x14, x14, x14, lsr #7
	eor	x0, x14, x14, lsl #17
	eor	x14, x0, x0, lsl #13
	eor	x14, x14, x14, lsr #7
	eor	x14, x14, x14, lsl #17
	and	x1, x14, #0x7fffffffffffffff
	mov	x2, x1
	mov	x3, x8
	cmp	x1, x8
	b.ne	LBB0_11
; %bb.13:                               ;   in Loop: Header=BB0_12 Depth=2
	mov	x2, x0
	mov	x3, x10
	cmp	x0, x10
	b.ne	LBB0_11
; %bb.14:                               ;   in Loop: Header=BB0_12 Depth=2
	mov	x2, x17
	mov	x3, x9
	cmp	x17, x9
	b.ne	LBB0_11
; %bb.15:                               ;   in Loop: Header=BB0_12 Depth=2
	mov	x3, #-4294967295                ; =0xffffffff00000001
	mov	x2, x16
	cmp	x16, x3
	b.ne	LBB0_11
	b	LBB0_9
LBB0_16:
	str	x14, [x11, _rngs@PAGEOFF]
LBB0_17:
	cmp	w20, #1
	b.lt	LBB0_27
; %bb.18:
	mov	x12, #0                         ; =0x0
	ldr	x13, [x11, _rngs@PAGEOFF]
	b	LBB0_20
LBB0_19:                                ;   in Loop: Header=BB0_20 Depth=1
	stp	x15, x16, [x14]
	stp	x17, x0, [x14, #16]
	add	x12, x12, #1
	cmp	x12, x20
	b.eq	LBB0_26
LBB0_20:                                ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_22 Depth 2
	ldur	x14, [x29, #-120]               ; 8-byte Folded Reload
	add	x14, x14, x12, lsl #5
	b	LBB0_22
LBB0_21:                                ;   in Loop: Header=BB0_22 Depth=2
	cmp	x1, x2
	b.ls	LBB0_19
LBB0_22:                                ;   Parent Loop BB0_20 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	eor	x13, x13, x13, lsl #13
	eor	x13, x13, x13, lsr #7
	eor	x15, x13, x13, lsl #17
	eor	x13, x15, x15, lsl #13
	eor	x13, x13, x13, lsr #7
	eor	x16, x13, x13, lsl #17
	eor	x13, x16, x16, lsl #13
	eor	x13, x13, x13, lsr #7
	eor	x17, x13, x13, lsl #17
	eor	x13, x17, x17, lsl #13
	eor	x13, x13, x13, lsr #7
	eor	x13, x13, x13, lsl #17
	and	x0, x13, #0x7fffffffffffffff
	mov	x1, x0
	mov	x2, x8
	cmp	x0, x8
	b.ne	LBB0_21
; %bb.23:                               ;   in Loop: Header=BB0_22 Depth=2
	mov	x1, x17
	mov	x2, x10
	cmp	x17, x10
	b.ne	LBB0_21
; %bb.24:                               ;   in Loop: Header=BB0_22 Depth=2
	mov	x1, x16
	mov	x2, x9
	cmp	x16, x9
	b.ne	LBB0_21
; %bb.25:                               ;   in Loop: Header=BB0_22 Depth=2
	mov	x2, #-4294967295                ; =0xffffffff00000001
	mov	x1, x15
	cmp	x15, x2
	b.ne	LBB0_21
	b	LBB0_19
LBB0_26:
	str	x13, [x11, _rngs@PAGEOFF]
LBB0_27:
	str	x24, [sp, #120]                 ; 8-byte Folded Spill
	str	w23, [sp, #92]                  ; 4-byte Folded Spill
	cmp	w22, #1
	stur	x20, [x29, #-112]               ; 8-byte Folded Spill
	str	x26, [sp, #112]                 ; 8-byte Folded Spill
	str	x22, [sp, #96]                  ; 8-byte Folded Spill
	b.lt	LBB0_36
; %bb.28:
	mov	w9, #0                          ; =0x0
	mov	w23, w20
	lsl	x24, x20, #5
	mov	w22, #80                        ; =0x50
	mov	w20, w26
	b	LBB0_30
LBB0_29:                                ;   in Loop: Header=BB0_30 Depth=1
	ldur	w9, [x29, #-124]                ; 4-byte Folded Reload
	add	w9, w9, #1
	ldr	x8, [sp, #96]                   ; 8-byte Folded Reload
	cmp	w9, w8
	b.eq	LBB0_42
LBB0_30:                                ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_33 Depth 2
                                        ;       Child Loop BB0_35 Depth 3
	stur	w9, [x29, #-124]                ; 4-byte Folded Spill
	ldur	x0, [x29, #-96]                 ; 8-byte Folded Reload
	ldr	x1, [sp, #104]                  ; 8-byte Folded Reload
	bl	_bzero
	ldr	x8, [sp, #112]                  ; 8-byte Folded Reload
	cmp	w8, #1
	b.lt	LBB0_29
; %bb.31:                               ;   in Loop: Header=BB0_30 Depth=1
	mov	x21, #0                         ; =0x0
	ldr	x19, [sp, #120]                 ; 8-byte Folded Reload
	b	LBB0_33
LBB0_32:                                ;   in Loop: Header=BB0_33 Depth=2
	ldp	x8, x9, [x29, #-104]            ; 16-byte Folded Reload
	madd	x0, x21, x22, x9
	add	x1, x8, x21, lsl #5
	bl	_row_finish
	add	x21, x21, #1
	add	x19, x19, x24
	cmp	x21, x20
	b.eq	LBB0_29
LBB0_33:                                ;   Parent Loop BB0_30 Depth=1
                                        ; =>  This Loop Header: Depth=2
                                        ;       Child Loop BB0_35 Depth 3
	ldur	x8, [x29, #-112]                ; 8-byte Folded Reload
	cmp	w8, #1
	b.lt	LBB0_32
; %bb.34:                               ;   in Loop: Header=BB0_33 Depth=2
	ldur	x8, [x29, #-96]                 ; 8-byte Folded Reload
	madd	x27, x21, x22, x8
	mov	x26, x23
	mov	x28, x19
	ldur	x25, [x29, #-120]               ; 8-byte Folded Reload
LBB0_35:                                ;   Parent Loop BB0_30 Depth=1
                                        ;     Parent Loop BB0_33 Depth=2
                                        ; =>    This Inner Loop Header: Depth=3
	mov	x0, x27
	mov	x1, x28
	mov	x2, x25
	bl	_big_mac
	add	x25, x25, #32
	add	x28, x28, #32
	subs	x26, x26, #1
	b.ne	LBB0_35
	b	LBB0_32
LBB0_36:
	cbnz	w22, LBB0_42
; %bb.37:
	ldr	x8, [sp, #112]                  ; 8-byte Folded Reload
	cmp	w8, #1
	b.lt	LBB0_42
; %bb.38:
	mov	x19, #0                         ; =0x0
	ldur	x8, [x29, #-112]                ; 8-byte Folded Reload
	mov	w20, w8
	lsl	x22, x8, #5
	mov	w24, #80                        ; =0x50
	movi.2d	v0, #0000000000000000
	ldr	x8, [sp, #112]                  ; 8-byte Folded Reload
	mov	w25, w8
	b	LBB0_40
LBB0_39:                                ;   in Loop: Header=BB0_40 Depth=1
	ldur	x8, [x29, #-104]                ; 8-byte Folded Reload
	add	x1, x8, x19, lsl #5
	mov	x0, x26
	bl	_row_finish
	movi.2d	v0, #0000000000000000
	stp	q0, q0, [x26, #48]
	stp	q0, q0, [x26, #16]
	add	x19, x19, #1
	ldr	x8, [sp, #120]                  ; 8-byte Folded Reload
	add	x8, x8, x22
	str	x8, [sp, #120]                  ; 8-byte Folded Spill
	str	q0, [x26]
	cmp	x19, x25
	b.eq	LBB0_42
LBB0_40:                                ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_41 Depth 2
	ldur	x8, [x29, #-96]                 ; 8-byte Folded Reload
	madd	x26, x19, x24, x8
	stp	q0, q0, [x26, #48]
	stp	q0, q0, [x26, #16]
	str	q0, [x26]
	mov	x21, x20
	ldr	x27, [sp, #120]                 ; 8-byte Folded Reload
	ldp	x28, x8, [x29, #-120]           ; 16-byte Folded Reload
	cmp	w8, #1
	b.lt	LBB0_39
LBB0_41:                                ;   Parent Loop BB0_40 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	mov	x0, x26
	mov	x1, x27
	mov	x2, x28
	bl	_big_mac
	add	x28, x28, #32
	add	x27, x27, #32
	subs	x21, x21, #1
	b.ne	LBB0_41
	b	LBB0_39
LBB0_42:
	ldr	w8, [sp, #92]                   ; 4-byte Folded Reload
	cbz	w8, LBB0_46
; %bb.43:
	ldr	x17, [sp, #112]                 ; 8-byte Folded Reload
	cmp	w17, #1
	ldur	x16, [x29, #-112]               ; 8-byte Folded Reload
	b.lt	LBB0_52
; %bb.44:
	mov	w8, w17
	cmp	w17, #9
	ldr	x0, [sp, #96]                   ; 8-byte Folded Reload
	b.hs	LBB0_53
; %bb.45:
	mov	x9, #0                          ; =0x0
	mov	x10, #0                         ; =0x0
	b	LBB0_56
LBB0_46:
	ldr	x8, [sp, #112]                  ; 8-byte Folded Reload
	cmp	w8, #0
	ldur	x23, [x29, #-112]               ; 8-byte Folded Reload
	b.le	LBB0_61
; %bb.47:
	mov	x20, #0                         ; =0x0
	mov	w22, w8
	ldur	x8, [x29, #-104]                ; 8-byte Folded Reload
	add	x21, x8, #16
Lloh7:
	adrp	x19, l_.str.2@PAGE
Lloh8:
	add	x19, x19, l_.str.2@PAGEOFF
LBB0_48:                                ; =>This Inner Loop Header: Depth=1
	ldp	x9, x8, [x21]
	ldp	x11, x10, [x21, #-16]
	stp	x10, x11, [sp, #32]
	stp	x8, x9, [sp, #16]
	stp	x20, x23, [sp]
	mov	x0, x19
	bl	_printf
	add	x20, x20, #1
	add	x21, x21, #32
	cmp	x22, x20
	b.ne	LBB0_48
; %bb.49:
Lloh9:
	adrp	x0, l_.str.3@PAGE
Lloh10:
	add	x0, x0, l_.str.3@PAGEOFF
	bl	_getenv
	cbz	x0, LBB0_59
; %bb.50:
	mov	x20, #0                         ; =0x0
	ldur	x8, [x29, #-96]                 ; 8-byte Folded Reload
	add	x21, x8, #40
Lloh11:
	adrp	x19, l_.str.4@PAGE
Lloh12:
	add	x19, x19, l_.str.4@PAGEOFF
LBB0_51:                                ; =>This Inner Loop Header: Depth=1
	ldp	x9, x8, [x21, #24]
	ldp	x11, x10, [x21, #8]
	ldp	x13, x12, [x21, #-8]
	ldp	x15, x14, [x21, #-24]
	ldp	x17, x16, [x21, #-40]
	stp	x16, x17, [sp, #72]
	stp	x14, x15, [sp, #56]
	stp	x12, x13, [sp, #40]
	stp	x10, x11, [sp, #24]
	stp	x8, x9, [sp, #8]
	str	x20, [sp]
	mov	x0, x19
	bl	_printf
	add	x20, x20, #1
	add	x21, x21, #80
	cmp	x22, x20
	b.ne	LBB0_51
	b	LBB0_59
LBB0_52:
	mov	x10, #0                         ; =0x0
	ldr	x0, [sp, #96]                   ; 8-byte Folded Reload
	b	LBB0_58
LBB0_53:
	ands	x9, x8, #0x7
	mov	w10, #8                         ; =0x8
	csel	x9, x10, x9, eq
	sub	x9, x8, x9
	ldur	x10, [x29, #-104]               ; 8-byte Folded Reload
	add	x10, x10, #128
	movi.2d	v0, #0000000000000000
	mov	x11, x9
	movi.2d	v1, #0000000000000000
	movi.2d	v2, #0000000000000000
	movi.2d	v3, #0000000000000000
LBB0_54:                                ; =>This Inner Loop Header: Depth=1
	sub	x12, x10, #96
	sub	x13, x10, #32
	add	x14, x10, #32
	add	x15, x10, #96
	ldur	d4, [x10, #-128]
	ld1.d	{ v4 }[1], [x12]
	ldur	d5, [x10, #-64]
	ld1.d	{ v5 }[1], [x13]
	ldr	d6, [x10]
	ld1.d	{ v6 }[1], [x14]
	ldr	d7, [x10, #64]
	ld1.d	{ v7 }[1], [x15]
	eor.16b	v0, v4, v0
	eor.16b	v1, v5, v1
	eor.16b	v2, v6, v2
	eor.16b	v3, v7, v3
	add	x10, x10, #256
	subs	x11, x11, #8
	b.ne	LBB0_54
; %bb.55:
	eor.16b	v0, v1, v0
	eor3.16b	v0, v2, v0, v3
	ext.16b	v1, v0, v0, #8
	eor.8b	v0, v0, v1
	fmov	x10, d0
LBB0_56:
	ldur	x11, [x29, #-104]               ; 8-byte Folded Reload
	add	x11, x11, x9, lsl #5
	sub	x8, x8, x9
LBB0_57:                                ; =>This Inner Loop Header: Depth=1
	ldr	x9, [x11], #32
	eor	x10, x9, x10
	subs	x8, x8, #1
	b.ne	LBB0_57
LBB0_58:
	stp	x0, x10, [sp, #16]
	mov	x8, x17
	stp	x8, x16, [sp]
Lloh13:
	adrp	x0, l_.str.5@PAGE
Lloh14:
	add	x0, x0, l_.str.5@PAGEOFF
	bl	_printf
LBB0_59:
	mov	w19, #0                         ; =0x0
LBB0_60:
	mov	x0, x19
	ldp	x29, x30, [sp, #256]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #240]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #224]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #208]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #192]            ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #176]            ; 16-byte Folded Reload
	add	sp, sp, #272
	ret
LBB0_61:
Lloh15:
	adrp	x0, l_.str.3@PAGE
Lloh16:
	add	x0, x0, l_.str.3@PAGEOFF
	bl	_getenv
	b	LBB0_59
	.loh AdrpAdd	Lloh3, Lloh4
	.loh AdrpLdrGotLdr	Lloh0, Lloh1, Lloh2
	.loh AdrpAdd	Lloh5, Lloh6
	.loh AdrpAdd	Lloh7, Lloh8
	.loh AdrpAdd	Lloh9, Lloh10
	.loh AdrpAdd	Lloh11, Lloh12
	.loh AdrpAdd	Lloh13, Lloh14
	.loh AdrpAdd	Lloh15, Lloh16
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function big_mac
_big_mac:                               ; @big_mac
	.cfi_startproc
; %bb.0:
	stp	x20, x19, [sp, #-16]!           ; 16-byte Folded Spill
	.cfi_def_cfa_offset 16
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	ldp	x9, x8, [x2, #16]
	ldp	x14, x10, [x2]
	ldp	x11, x13, [x1]
	mul	x15, x14, x11
	umulh	x12, x14, x11
	umulh	x16, x10, x11
	mul	x17, x10, x11
	adds	x12, x17, x12
	cinc	x16, x16, hs
	umulh	x17, x9, x11
	mul	x2, x9, x11
	adds	x2, x16, x2
	cinc	x16, x17, hs
	umulh	x17, x8, x11
	mul	x11, x8, x11
	adds	x11, x16, x11
	cinc	x16, x17, hs
	umulh	x17, x14, x13
	mul	x3, x14, x13
	adds	x3, x12, x3
	cinc	x12, x17, hs
	umulh	x17, x10, x13
	mul	x4, x10, x13
	adds	x2, x2, x4
	cset	w4, hs
	adds	x2, x2, x12
	adc	x12, x4, x17
	umulh	x17, x9, x13
	mul	x4, x9, x13
	adds	x11, x11, x4
	cset	w4, hs
	adds	x5, x11, x12
	adc	x4, x4, x17
	ldp	x12, x11, [x1, #16]
	umulh	x17, x14, x12
	mul	x1, x14, x12
	adds	x6, x2, x1
	cinc	x17, x17, hs
	umulh	x1, x10, x12
	mul	x2, x10, x12
	adds	x2, x5, x2
	cset	w5, hs
	adds	x17, x2, x17
	adc	x5, x5, x1
	umulh	x1, x14, x11
	mul	x14, x14, x11
	adds	x14, x17, x14
	cinc	x7, x1, hs
	ldp	x19, x2, [x0, #16]
	adds	x1, x6, x19
	ldp	x17, x20, [x0]
	adds	x15, x17, x15
	adcs	x17, x3, x20
	adcs	x1, x1, xzr
	cset	w3, hs
	adds	x20, x14, x2
	cmn	x6, x19
	umulh	x6, x8, x13
	mul	x19, x8, x13
	adcs	x13, x20, x3
	cset	w3, hs
	adds	x16, x16, x19
	cset	w19, hs
	adds	x16, x16, x4
	umulh	x4, x9, x12
	mul	x20, x9, x12
	adc	x6, x19, x6
	adds	x16, x16, x20
	cset	w19, hs
	adds	x16, x16, x5
	umulh	x5, x10, x11
	mul	x10, x10, x11
	adc	x4, x19, x4
	adds	x10, x16, x10
	cset	w16, hs
	adds	x10, x10, x7
	adc	x16, x16, x5
	ldp	x5, x7, [x0, #32]
	adds	x19, x10, x5
	cmn	x14, x2
	umulh	x14, x8, x12
	mul	x12, x8, x12
	adcs	x2, x19, x3
	cset	w3, hs
	adds	x12, x6, x12
	cset	w6, hs
	adds	x12, x12, x4
	umulh	x4, x9, x11
	mul	x9, x9, x11
	adc	x14, x6, x14
	adds	x9, x12, x9
	cset	w12, hs
	adds	x9, x9, x16
	stp	x15, x17, [x0]
	stp	x1, x13, [x0, #16]
	adc	x12, x12, x4
	adds	x13, x9, x7
	cmn	x10, x5
	umulh	x10, x8, x11
	mul	x8, x8, x11
	adcs	x11, x13, x3
	cset	w13, hs
	adds	x8, x14, x8
	cset	w14, hs
	adds	x8, x8, x12
	adc	x10, x14, x10
	stp	x2, x11, [x0, #32]
	ldp	x11, x12, [x0, #48]
	adds	x14, x8, x11
	cmn	x9, x7
	adcs	x9, x14, x13
	cset	w13, hs
	adds	x14, x10, x12
	cmn	x8, x11
	adcs	x8, x14, x13
	cset	w11, hs
	cmn	x10, x12
	stp	x9, x8, [x0, #48]
	cinc	x8, x11, hs
	cbz	x8, LBB1_2
; %bb.1:
	ldr	x9, [x0, #64]
	adds	x8, x9, x8
	str	x8, [x0, #64]
	b.eq	LBB1_3
LBB1_2:
	ldp	x20, x19, [sp], #16             ; 16-byte Folded Reload
	ret
LBB1_3:
	ldr	x8, [x0, #72]
	add	x8, x8, #1
	str	x8, [x0, #72]
	ldp	x20, x19, [sp], #16             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function row_finish
_row_finish:                            ; @row_finish
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #336
	stp	x26, x25, [sp, #256]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #272]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #288]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #304]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #320]            ; 16-byte Folded Spill
	add	x29, sp, #320
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
	mov	x14, #0                         ; =0x0
Lloh17:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh18:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh19:
	ldr	x8, [x8]
	stur	x8, [x29, #-72]
	ldp	x13, x12, [x0, #32]
	extr	x8, x12, x13, #2
	ldp	x11, x10, [x0, #48]
	extr	x9, x11, x12, #2
	extr	x15, x10, x11, #2
	stp	x8, x9, [x29, #-112]
	ldp	x9, x8, [x0, #64]
	extr	x16, x9, x10, #2
	stp	x15, x16, [x29, #-96]
	extr	x15, x8, x9, #2
	stur	x15, [x29, #-80]
	movi.2d	v0, #0000000000000000
	stp	q0, q0, [sp, #176]
	stp	q0, q0, [sp, #144]
	str	q0, [sp, #128]
	sub	x15, x29, #112
	add	x16, sp, #128
	mov	x17, #14641                     ; =0x3931
	movk	x17, #49366, lsl #16
	movk	x17, #40992, lsl #32
	movk	x17, #10039, lsl #48
	ldr	x6, [sp, #160]
	mov	x2, #55060                      ; =0xd714
	movk	x2, #58554, lsl #16
	movk	x2, #60235, lsl #32
	movk	x2, #20547, lsl #48
	mov	x3, #60726                      ; =0xed36
	movk	x3, #32264, lsl #16
	movk	x3, #52080, lsl #32
	movk	x3, #35677, lsl #48
	mov	x4, #30659                      ; =0x77c3
	movk	x4, #60899, lsl #16
	movk	x4, #20191, lsl #32
	movk	x4, #21769, lsl #48
	mov	w5, #35                         ; =0x23
LBB2_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x7, [x15, x14]
	add	x19, x16, x14
	mul	x20, x7, x17
	ldp	x21, x22, [x19]
	adds	x20, x21, x20
	umulh	x21, x7, x17
	adcs	x21, x22, x21
	cset	w22, hs
	mul	x23, x7, x2
	adds	x21, x21, x23
	stp	x20, x21, [x19]
	umulh	x20, x7, x2
	adc	x20, x20, x22
	ldp	x21, x22, [x19, #16]
	adds	x20, x21, x20
	cset	w21, hs
	mul	x23, x7, x3
	adds	x20, x20, x23
	umulh	x23, x7, x3
	adc	x21, x23, x21
	adds	x21, x22, x21
	cset	w22, hs
	mul	x23, x7, x4
	adds	x21, x21, x23
	stp	x20, x21, [x19, #16]
	umulh	x20, x7, x4
	adc	x20, x20, x22
	adds	x6, x6, x20
	cset	w20, hs
	mul	x21, x7, x5
	adds	x21, x6, x21
	umulh	x6, x7, x5
	adc	x6, x6, x20
	ldr	x7, [x19, #40]
	add	x6, x7, x6
	stp	x21, x6, [x19, #32]
	add	x14, x14, #8
	cmp	x14, #40
	b.ne	LBB2_1
; %bb.2:
	mov	x14, #0                         ; =0x0
	ldp	x15, x16, [sp, #160]
	extr	x15, x16, x15, #2
	ldp	x17, x2, [sp, #176]
	extr	x16, x17, x16, #2
	stp	x15, x16, [sp, #88]
	extr	x15, x2, x17, #2
	ldp	x16, x17, [sp, #192]
	extr	x2, x16, x2, #2
	stp	x15, x2, [sp, #104]
	extr	x15, x17, x16, #2
	str	x15, [sp, #120]
	movi.2d	v0, #0000000000000000
	stp	q0, q0, [sp, #48]
	stp	q0, q0, [sp, #16]
	add	x15, sp, #88
	mov	x16, sp
	mov	x17, #-4294967295               ; =0xffffffff00000001
	str	q0, [sp]
	mov	x2, #23550                      ; =0x5bfe
	movk	x2, #65534, lsl #16
	movk	x2, #41986, lsl #32
	movk	x2, #21437, lsl #48
	mov	x3, #55301                      ; =0xd805
	movk	x3, #2465, lsl #16
	movk	x3, #55304, lsl #32
	movk	x3, #13113, lsl #48
	mov	x4, #32072                      ; =0x7d48
	movk	x4, #10653, lsl #16
	movk	x4, #42835, lsl #32
	movk	x4, #29677, lsl #48
	ldr	x5, [sp, #24]
LBB2_3:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x6, [x15, x14]
	add	x7, x16, x14
	sub	x19, x6, x6, lsl #32
	ldp	x20, x21, [x7]
	adds	x19, x20, x19
	umulh	x20, x6, x17
	cinc	x20, x20, hs
	adds	x20, x21, x20
	cset	w21, hs
	mul	x22, x6, x2
	adds	x20, x20, x22
	stp	x19, x20, [x7]
	umulh	x19, x6, x2
	adc	x19, x19, x21
	ldr	x20, [x7, #16]
	adds	x19, x20, x19
	cset	w20, hs
	mul	x21, x6, x3
	adds	x19, x19, x21
	umulh	x21, x6, x3
	adc	x20, x21, x20
	adds	x5, x5, x20
	cset	w20, hs
	mul	x21, x6, x4
	adds	x5, x5, x21
	stp	x19, x5, [x7, #16]
	umulh	x5, x6, x4
	adc	x5, x5, x20
	ldr	x6, [x7, #32]
	add	x5, x6, x5
	str	x5, [x7, #32]
	add	x14, x14, #8
	cmp	x14, #40
	b.ne	LBB2_3
; %bb.4:
	mov	x14, #10234                     ; =0x27fa
	movk	x14, #63070, lsl #16
	movk	x14, #10231, lsl #32
	movk	x14, #52422, lsl #48
	mov	x15, #41985                     ; =0xa401
	movk	x15, #1, lsl #16
	movk	x15, #23549, lsl #32
	movk	x15, #44098, lsl #48
	mov	x16, #33463                     ; =0x82b7
	movk	x16, #54882, lsl #16
	movk	x16, #22700, lsl #32
	movk	x16, #35858, lsl #48
	mov	x17, #32072                     ; =0x7d48
	movk	x17, #10653, lsl #16
	movk	x17, #42835, lsl #32
	movk	x17, #29677, lsl #48
	mov	x2, #23550                      ; =0x5bfe
	movk	x2, #65534, lsl #16
	movk	x2, #41986, lsl #32
	movk	x2, #21437, lsl #48
	mov	x3, #55301                      ; =0xd805
	movk	x3, #2465, lsl #16
	movk	x3, #55304, lsl #32
	movk	x3, #13113, lsl #48
	ldp	x4, x6, [sp]
	ldp	x5, x7, [x0]
	subs	x4, x5, x4
	cset	w19, lo
	cinc	x5, x6, lo
	sub	x5, x7, x5
	cmp	x7, x6
	csel	w6, wzr, w19, ne
	csinc	w7, w6, wzr, hs
	ldp	x6, x20, [sp, #16]
	ldp	x19, x0, [x0, #16]
	cmp	x19, x6
	sub	x6, x19, x6
	sub	x6, x6, x7
	csel	w7, wzr, w7, ne
	csinc	w7, w7, wzr, hs
	sub	x19, x0, x20
	sub	x19, x19, x7
	cmp	x0, x20
	csel	w0, wzr, w7, ne
	csinc	w7, w0, wzr, hs
	ldp	x0, x20, [sp, #32]
	cmp	x13, x0
	sub	x13, x13, x0
	sub	x0, x13, x7
	csel	w13, wzr, w7, ne
	csinc	w7, w13, wzr, hs
	sub	x13, x12, x20
	sub	x13, x13, x7
	cmp	x12, x20
	csel	w12, wzr, w7, ne
	csinc	w12, w12, wzr, hs
	ldp	x7, x20, [sp, #48]
	cmp	x11, x7
	sub	x11, x11, x7
	sub	x11, x11, x12
	csel	w12, wzr, w12, ne
	csinc	w7, w12, wzr, hs
	sub	x12, x10, x20
	sub	x12, x12, x7
	cmp	x10, x20
	csel	w10, wzr, w7, ne
	csinc	w7, w10, wzr, hs
	ldr	x20, [sp, #64]
	sub	x10, x9, x20
	sub	x10, x10, x7
	cmp	x9, x20
	csel	w9, wzr, w7, ne
	csinc	w9, w9, wzr, hs
	sub	x8, x8, x9
	mov	w9, #-1                         ; =0xffffffff
	mov	x7, #-4294967295                ; =0xffffffff00000001
	b	LBB2_6
LBB2_5:                                 ;   in Loop: Header=BB2_6 Depth=1
	sub	x20, x0, #1
	sub	x21, x13, #1
	sub	x22, x11, #1
	sub	x23, x12, #1
	sub	x24, x10, #1
	sub	x25, x8, #1
	cmp	x10, #0
	csel	x25, x25, x8, eq
	cmp	x12, #0
	csel	x24, x10, x24, ne
	csel	x25, x8, x25, ne
	cmp	x11, #0
	csel	x24, x10, x24, ne
	csel	x23, x12, x23, ne
	csel	x25, x8, x25, ne
	cmp	x13, #0
	csel	x24, x10, x24, ne
	csel	x23, x12, x23, ne
	csel	x22, x11, x22, ne
	csel	x25, x8, x25, ne
	cmp	x0, #0
	csel	x10, x10, x24, ne
	csel	x12, x12, x23, ne
	csel	x11, x11, x22, ne
	csel	x13, x13, x21, ne
	csel	x8, x8, x25, ne
	mov	x0, x20
LBB2_6:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB2_7 Depth 2
	orr	x20, x8, x10
	orr	x21, x12, x11
	orr	x20, x20, x21
	orr	x21, x13, x0
	orr	x20, x20, x21
	cmp	x20, #0
	cset	w20, eq
	mov	x21, x19
	mov	x22, x6
	mov	x23, x5
	mov	x24, x4
LBB2_7:                                 ;   Parent Loop BB2_6 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	cbz	w20, LBB2_13
; %bb.8:                                ;   in Loop: Header=BB2_7 Depth=2
	mov	x4, x21
	mov	x5, x17
	cmp	x21, x17
	b.ne	LBB2_12
; %bb.9:                                ;   in Loop: Header=BB2_7 Depth=2
	mov	x4, x22
	mov	x5, x3
	cmp	x22, x3
	b.ne	LBB2_12
; %bb.10:                               ;   in Loop: Header=BB2_7 Depth=2
	mov	x4, x23
	mov	x5, x2
	cmp	x23, x2
	b.ne	LBB2_12
; %bb.11:                               ;   in Loop: Header=BB2_7 Depth=2
	mov	x5, #-4294967295                ; =0xffffffff00000001
	mov	x4, x24
	cmp	x24, x5
	b.eq	LBB2_13
LBB2_12:                                ;   in Loop: Header=BB2_7 Depth=2
	cmp	x4, x5
	b.ls	LBB2_15
LBB2_13:                                ;   in Loop: Header=BB2_7 Depth=2
	add	x4, x24, x9
	cmp	x24, x7
	cset	w6, lo
	cinc	x5, x15, hs
	add	x5, x5, x23
	cmp	x23, x2
	csel	w6, wzr, w6, ne
	csinc	w19, w6, wzr, hs
	cmp	w19, #0
	cinc	x6, x14, eq
	add	x6, x6, x22
	cmp	x22, x3
	csel	w19, wzr, w19, ne
	csinc	w22, w19, wzr, hs
	cmp	w22, #0
	cinc	x19, x16, eq
	add	x19, x19, x21
	cmp	x21, x17
	csel	w25, wzr, w22, ne
	b.lo	LBB2_5
; %bb.14:                               ;   in Loop: Header=BB2_7 Depth=2
	mov	x21, x19
	mov	x22, x6
	mov	x23, x5
	mov	x24, x4
	cbz	w25, LBB2_7
	b	LBB2_5
LBB2_15:
	stp	x24, x23, [x1]
	stp	x22, x21, [x1, #16]
	ldur	x8, [x29, #-72]
Lloh20:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh21:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh22:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB2_17
; %bb.16:
	ldp	x29, x30, [sp, #320]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #304]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #288]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #272]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #256]            ; 16-byte Folded Reload
	add	sp, sp, #336
	ret
LBB2_17:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh17, Lloh18, Lloh19
	.loh AdrpLdrGotLdr	Lloh20, Lloh21, Lloh22
	.cfi_endproc
                                        ; -- End function
	.section	__TEXT,__cstring,cstring_literals
l_.str:                                 ; @.str
	.asciz	"usage: la_big check | la_big bench rows W reps\n"

l_.str.1:                               ; @.str.1
	.asciz	"check"

l_.str.2:                               ; @.str.2
	.asciz	"%d %d %016llx%016llx%016llx%016llx\n"

l_.str.3:                               ; @.str.3
	.asciz	"LA_BIG_RAW"

l_.str.4:                               ; @.str.4
	.asciz	"acc %d %016llx%016llx%016llx%016llx%016llx%016llx%016llx%016llx%016llx%016llx\n"

l_.str.5:                               ; @.str.5
	.asciz	"bench rows=%d W=%d reps=%d sink=%llu\n"

	.section	__DATA,__data
	.p2align	3, 0x0                          ; @rngs
_rngs:
	.quad	88172645463325252               ; 0x139408dcbbf7a44

.subsections_via_symbols

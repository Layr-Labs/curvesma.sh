	.build_version macos, 26, 0	sdk_version 27, 0
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_mac_run                        ; -- Begin function mac_run
	.p2align	2
_mac_run:                               ; @mac_run
	ldp	x8, x12, [x0]
	ldp	x11, x10, [x0, #16]
	ldr	x9, [x0, #32]
	cmp	x4, #1
	b.lt	LBB0_4
; %bb.1:
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x13, [x2], #8
	add	x13, x1, x13, lsl #5
	ldrb	w14, [x3], #1
	ldp	x15, x16, [x13]
	mul	x17, x15, x14
	umulh	x15, x15, x14
	umulh	x5, x16, x14
	mul	x16, x16, x14
	ldp	x6, x13, [x13, #16]
	umulh	x7, x6, x14
	mul	x6, x6, x14
	umulh	x19, x13, x14
	mul	x13, x13, x14
	adds	x12, x15, x12
	cset	w14, hs
	adds	x12, x12, x16
	cinc	x14, x14, hs
	adds	x11, x5, x11
	cset	w15, hs
	adds	x11, x11, x6
	cinc	x15, x15, hs
	adds	x10, x7, x10
	cset	w16, hs
	adds	x10, x10, x13
	cinc	x13, x16, hs
	adds	x8, x17, x8
	adcs	x12, x12, xzr
	adcs	x11, x11, x14
	adcs	x10, x10, x15
	add	x9, x9, x19
	adc	x9, x9, x13
	subs	x4, x4, #1
	b.ne	LBB0_2
; %bb.3:
LBB0_4:
	stp	x8, x12, [x0]
	stp	x11, x10, [x0, #16]
	str	x9, [x0, #32]
	ret
	.globl	_row_dot                        ; -- Begin function row_dot
	.p2align	2
_row_dot:                               ; @row_dot
	add	x29, sp, #48
	mov	x19, x6
	mov	x20, x5
	mov	x21, x4
	mov	x22, x3
	mov	x23, x2
	mov	x24, x1
	mov	x1, x2
	mov	x2, x3
	mov	x3, x4
	mov	x4, x5
	bl	_mac_run
	add	x2, x22, x20, lsl #3
	sub	x4, x19, x20
	add	x3, x21, x20
	mov	x0, x24
	mov	x1, x23
	b	_mac_run
	.globl	_group_run                      ; -- Begin function group_run
	.p2align	2
_group_run:                             ; @group_run
; %bb.0:
	stp	x20, x19, [sp, #-16]!           ; 16-byte Folded Spill
	ldp	x8, x12, [x0]
	ldp	x11, x10, [x0, #16]
	ldr	x9, [x0, #32]
	cmp	x3, #1
	b.lt	LBB9_2
LBB9_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x13, [x2], #8
	add	x13, x1, x13, lsl #5
	ldp	x14, x15, [x13]
	adds	x8, x14, x8
	adcs	x12, x15, x12
	ldp	x14, x13, [x13, #16]
	adcs	x11, x14, x11
	adcs	x10, x13, x10
	cinc	x9, x9, hs
	subs	x3, x3, #1
	b.ne	LBB9_1
LBB9_2:
	cmp	x6, #1
	b.lt	LBB9_4
LBB9_3:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x13, [x4], #8
	add	x13, x1, x13, lsl #5
	ldrb	w14, [x5], #1
	ldp	x15, x16, [x13]
	mul	x17, x15, x14
	umulh	x15, x15, x14
	umulh	x2, x16, x14
	mul	x16, x16, x14
	ldp	x3, x13, [x13, #16]
	umulh	x7, x3, x14
	mul	x3, x3, x14
	umulh	x19, x13, x14
	mul	x13, x13, x14
	adds	x12, x15, x12
	cset	w14, hs
	adds	x12, x12, x16
	cinc	x14, x14, hs
	adds	x11, x2, x11
	cset	w15, hs
	adds	x11, x11, x3
	cinc	x15, x15, hs
	adds	x10, x7, x10
	cset	w16, hs
	adds	x10, x10, x13
	cinc	x13, x16, hs
	adds	x8, x17, x8
	adcs	x12, x12, xzr
	adcs	x11, x11, x14
	adcs	x10, x10, x15
	add	x9, x9, x19
	adc	x9, x9, x13
	subs	x6, x6, #1
	b.ne	LBB9_3
LBB9_4:
	stp	x8, x12, [x0]
	stp	x11, x10, [x0, #16]
	str	x9, [x0, #32]
	ldp	x20, x19, [sp], #16             ; 16-byte Folded Reload
	ret
                                        ; -- End function
.subsections_via_symbols

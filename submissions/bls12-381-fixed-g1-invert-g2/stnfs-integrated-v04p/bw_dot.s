	.build_version macos, 26, 0	sdk_version 27, 0
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_bw_dot                         ; -- Begin function bw_dot
	.p2align	2
_bw_dot:                                ; @bw_dot
	ldp	x8, x14, [x0]
	ldp	x13, x12, [x0, #16]
	ldp	x11, x9, [x0, #32]
	ldr	x10, [x0, #48]
	cmp	x3, #1
	b.lt	LBB0_4
; %bb.1:
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x15, [x1], #8
	ldp	x16, x17, [x2]
	mul	x4, x16, x15
	umulh	x16, x16, x15
	umulh	x5, x17, x15
	mul	x17, x17, x15
	ldp	x6, x7, [x2, #16]
	umulh	x19, x6, x15
	mul	x6, x6, x15
	umulh	x20, x7, x15
	mul	x15, x7, x15
	adds	x14, x16, x14
	cset	w16, hs
	adds	x14, x14, x17
	cinc	x16, x16, hs
	adds	x13, x5, x13
	cset	w17, hs
	adds	x13, x13, x6
	cinc	x17, x17, hs
	adds	x12, x19, x12
	cset	w5, hs
	adds	x12, x12, x15
	cinc	x15, x5, hs
	adds	x11, x20, x11
	cset	w5, hs
	adds	x8, x4, x8
	adcs	x14, x14, xzr
	adcs	x13, x13, x16
	adcs	x12, x12, x17
	adcs	x11, x11, x15
	adcs	x9, x9, x5
	cinc	x10, x10, hs
	add	x2, x2, #32
	subs	x3, x3, #1
	b.ne	LBB0_2
; %bb.3:
LBB0_4:
	stp	x8, x14, [x0]
	stp	x13, x12, [x0, #16]
	stp	x11, x9, [x0, #32]
	str	x10, [x0, #48]
	ret
.subsections_via_symbols
	.globl	_bw_full
	.p2align	2
_bw_full:
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldr	x16, [x0, #64]
	cmp	x3, #1
	b.lt	LBB1_4
; %bb.1:
LBB1_2:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x5, x4, [x2]
	ldp	x7, x20, [x1]
	mul	x17, x5, x7
	umulh	x6, x5, x7
	umulh	x19, x4, x7
	mul	x21, x4, x7
	adds	x21, x6, x21
	cinc	x22, x19, hs
	ldp	x19, x6, [x2, #16]
	umulh	x23, x19, x7
	mul	x24, x19, x7
	adds	x22, x22, x24
	cinc	x23, x23, hs
	umulh	x24, x6, x7
	mul	x7, x6, x7
	adds	x23, x23, x7
	umulh	x25, x5, x20
	mul	x7, x5, x20
	cinc	x24, x24, hs
	adds	x7, x21, x7
	cinc	x21, x25, hs
	umulh	x25, x4, x20
	mul	x26, x4, x20
	adds	x21, x21, x22
	cset	w22, hs
	adds	x21, x21, x26
	adc	x22, x22, x25
	umulh	x25, x19, x20
	mul	x26, x19, x20
	adds	x22, x22, x23
	cset	w23, hs
	adds	x22, x22, x26
	adc	x23, x23, x25
	umulh	x25, x6, x20
	mul	x20, x6, x20
	adds	x23, x23, x24
	cset	w24, hs
	adds	x20, x23, x20
	ldp	x23, x26, [x1, #16]
	umulh	x27, x5, x23
	mul	x28, x5, x23
	adc	x24, x24, x25
	adds	x21, x21, x28
	cinc	x25, x27, hs
	mul	x27, x4, x23
	adds	x22, x25, x22
	cset	w25, hs
	adds	x22, x22, x27
	umulh	x27, x4, x23
	adc	x25, x25, x27
	mul	x27, x19, x23
	adds	x20, x25, x20
	cset	w25, hs
	adds	x20, x20, x27
	umulh	x27, x19, x23
	adc	x25, x25, x27
	umulh	x27, x6, x23
	mul	x23, x6, x23
	adds	x24, x25, x24
	cset	w25, hs
	adds	x23, x24, x23
	adc	x24, x25, x27
	umulh	x25, x5, x26
	mul	x5, x5, x26
	adds	x5, x22, x5
	cinc	x22, x25, hs
	umulh	x25, x4, x26
	mul	x4, x4, x26
	adds	x20, x22, x20
	cset	w22, hs
	adds	x4, x20, x4
	adc	x20, x22, x25
	umulh	x22, x19, x26
	mul	x19, x19, x26
	adds	x20, x20, x23
	cset	w23, hs
	adds	x19, x20, x19
	adc	x20, x23, x22
	umulh	x22, x6, x26
	mul	x6, x6, x26
	adds	x20, x20, x24
	cset	w23, hs
	adds	x6, x20, x6
	adc	x20, x23, x22
	adds	x8, x8, x17
	adcs	x9, x7, x9
	adcs	x10, x21, x10
	adcs	x11, x5, x11
	adcs	x12, x4, x12
	adcs	x13, x19, x13
	adcs	x14, x6, x14
	adcs	x15, x20, x15
	cinc	x16, x16, hs
	add	x1, x1, #32
	add	x2, x2, #32
	subs	x3, x3, #1
	b.ne	LBB1_2
; %bb.3:
LBB1_4:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	str	x16, [x0, #64]
	ret
                                        ; -- End function
.subsections_via_symbols

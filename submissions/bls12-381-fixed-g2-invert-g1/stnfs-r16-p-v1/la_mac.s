	.section	__TEXT,__text,regular,pure_instructions
	.globl	_mac_run                        ; -- Begin function mac_run
	.p2align	2
_mac_run:
	ldp	x8, x12, [x0]
	ldp	x11, x10, [x0, #16]
	ldr	x9, [x0, #32]
	cmp	x4, #1
	b.lt	LBB0_4
LBB0_2:
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
;bb
LBB0_4:
	stp	x8, x12, [x0]
	stp	x11, x10, [x0, #16]
	str	x9, [x0, #32]
	ret
.subsections_via_symbols

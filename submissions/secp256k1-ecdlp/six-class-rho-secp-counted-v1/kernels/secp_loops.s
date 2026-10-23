pass1:
	sub	sp, sp, #128
	ldr	x8, [x2]
	ldp	x9, x12, [x1]
	ldp	x13, x14, [x1, #16]
	mov	x15, x2
	ldp	x10, x11, [x8]
	mov	x2, #977
	movk	x2, #1, lsl #32
	str	x30, [sp, #32]
	stp	x28, x27, [sp, #48]
	subs	x9, x10, x9
	stp	x26, x25, [sp, #64]
	sbcs	x11, x11, x12
	ldp	x10, x12, [x8, #16]
	stp	x24, x23, [sp, #80]
	stp	x22, x21, [sp, #96]
	sbcs	x10, x10, x13
	stp	x20, x19, [sp, #112]
	sbcs	x12, x12, x14
	csel	x13, x2, xzr, lo
	subs	x9, x9, x13
	cset	w13, lo
	subs	x11, x11, x13
	cset	w13, lo
	stp	x9, x11, [x3]
	subs	x10, x10, x13
	cset	w9, lo
	cmp	x0, #2
	sub	x9, x12, x9
	stp	x10, x9, [x3, #16]
	b.lo	.LBB0_3
	mov	x11, #977
	sub	x10, x0, #1
	add	x8, x15, #8
	add	x22, x1, #88
	add	x12, x3, #32
	movk	x11, #1, lsl #32
.LBB0_2:
	ldr	x13, [x8], #8
	str	x8, [sp, #24]
	ldp	x15, x16, [x22, #-24]
	str	x10, [sp, #40]
	ldp	x14, x17, [x13]
	ldp	x18, x13, [x13, #16]
	ldp	x25, x26, [x12, #-16]
	subs	x14, x14, x15
	ldp	x15, x0, [x22, #-8]
	sbcs	x16, x17, x16
	add	x22, x22, #64
	sbcs	x15, x18, x15
	sbcs	x0, x13, x0
	ldp	x13, x1, [x12, #-32]
	csel	x17, x2, xzr, lo
	subs	x14, x14, x17
	cset	w17, lo
	umulh	x18, x14, x13
	subs	x17, x16, x17
	cset	w16, lo
	mul	x9, x1, x14
	subs	x16, x15, x16
	mul	x3, x17, x13
	cset	w15, lo
	sub	x15, x0, x15
	umulh	x8, x16, x13
	umulh	x4, x17, x13
	adds	x18, x3, x18
	mul	x5, x16, x13
	stp	x8, x9, [sp, #8]
	mul	x0, x15, x13
	umulh	x7, x15, x13
	adcs	x3, x5, x4
	ldr	x5, [sp, #8]
	umulh	x6, x1, x14
	adcs	x0, x0, x5
	mul	x20, x1, x17
	cinc	x5, x7, hs
	umulh	x21, x1, x17
	mul	x19, x1, x16
	adds	x6, x20, x6
	umulh	x23, x1, x16
	mul	x24, x1, x15
	adcs	x19, x19, x21
	umulh	x1, x1, x15
	umulh	x28, x25, x14
	adcs	x20, x24, x23
	ldr	x23, [sp, #16]
	mul	x30, x25, x17
	cinc	x1, x1, hs
	adds	x18, x18, x23
	umulh	x8, x25, x17
	adcs	x3, x3, x6
	adcs	x0, x0, x19
	mul	x9, x25, x16
	adcs	x5, x5, x20
	cinc	x1, x1, hs
	adds	x19, x30, x28
	umulh	x2, x25, x16
	mul	x10, x25, x15
	adcs	x8, x9, x8
	mul	x27, x25, x14
	umulh	x25, x25, x15
	adcs	x10, x10, x2
	umulh	x7, x26, x14
	mul	x21, x26, x17
	umulh	x17, x26, x17
	mul	x6, x26, x16
	umulh	x9, x26, x16
	cinc	x16, x25, hs
	adds	x2, x3, x27
	adcs	x0, x0, x19
	mul	x3, x26, x15
	adcs	x8, x5, x8
	adcs	x10, x1, x10
	mul	x4, x26, x14
	cinc	x16, x16, hs
	adds	x1, x21, x7
	adcs	x17, x6, x17
	umulh	x15, x26, x15
	adcs	x9, x3, x9
	mul	x13, x14, x13
	cinc	x14, x15, hs
	adds	x15, x0, x4
	adcs	x8, x8, x1
	adcs	x10, x10, x17
	mul	x0, x8, x11
	adcs	x9, x16, x9
	umulh	x8, x8, x11
	cinc	x14, x14, hs
	mul	x17, x10, x11
	umulh	x10, x10, x11
	mul	x16, x9, x11
	adds	x8, x17, x8
	umulh	x9, x9, x11
	mul	x1, x14, x11
	adcs	x10, x16, x10
	umulh	x14, x14, x11
	adcs	x9, x1, x9
	cinc	x14, x14, hs
	adds	x13, x13, x0
	adcs	x8, x18, x8
	adcs	x10, x2, x10
	mov	x2, #977
	cset	w16, hs
	adcs	x17, x15, x9
	movk	x2, #1, lsl #32
	cinc	x14, x14, hs
	add	x9, x15, x9
	mul	x18, x14, x2
	umulh	x14, x14, x2
	adds	x13, x13, x18
	adcs	x8, x8, x14
	cset	w14, hs
	adcs	xzr, x10, xzr
	adcs	xzr, x17, xzr
	csel	x17, x2, xzr, hs
	adds	x13, x13, x17
	adcs	x8, x8, xzr
	stp	x13, x8, [x12]
	adcs	x8, x10, x14
	ldr	x10, [sp, #40]
	adc	x9, x9, x16
	stp	x8, x9, [x12, #16]
	ldr	x8, [sp, #24]
	subs	x10, x10, #1
	add	x12, x12, #32
	b.ne	.LBB0_2
.LBB0_3:
	ldp	x20, x19, [sp, #112]
	ldr	x30, [sp, #32]
	ldp	x22, x21, [sp, #96]
	ldp	x24, x23, [sp, #80]
	ldp	x26, x25, [sp, #64]
	ldp	x28, x27, [sp, #48]
	add	sp, sp, #128
	ret
pass2:
	sub	sp, sp, #304
	ldp	x9, x15, [x4]
	subs	x11, x0, #1
	ldp	x14, x10, [x4, #16]
	str	x30, [sp, #208]
	stp	x28, x27, [sp, #224]
	stp	x26, x25, [sp, #240]
	stp	x24, x23, [sp, #256]
	stp	x22, x21, [sp, #272]
	stp	x20, x19, [sp, #288]
	str	x1, [sp, #8]
	str	x10, [sp, #192]
	b.eq	.LBB1_3
	add	x12, x1, x0, lsl #6
	add	x13, x3, x0, lsl #5
	str	x2, [sp, #184]
	sub	x8, x12, #32
	sub	x12, x13, #64
	str	x8, [sp, #216]
.LBB1_2:
	ldp	x24, x17, [x12]
	mov	x22, x9
	ldp	x18, x16, [x12, #16]
	str	x11, [sp, #200]
	ldr	x10, [sp, #192]
	umulh	x0, x24, x9
	mul	x1, x17, x9
	umulh	x3, x17, x9
	mul	x4, x18, x9
	adds	x0, x1, x0
	umulh	x5, x18, x9
	mul	x6, x16, x9
	adcs	x1, x4, x3
	umulh	x7, x16, x9
	umulh	x20, x15, x24
	adcs	x4, x6, x5
	mul	x21, x15, x17
	cinc	x5, x7, hs
	mul	x8, x15, x24
	umulh	x19, x15, x17
	adds	x6, x21, x20
	mul	x23, x15, x18
	stp	x8, x12, [sp, #168]
	mov	x12, x14
	umulh	x25, x15, x18
	mul	x26, x15, x16
	adcs	x20, x23, x19
	umulh	x27, x15, x16
	umulh	x30, x14, x24
	adcs	x21, x26, x25
	ldr	x25, [sp, #168]
	mul	x11, x14, x17
	cinc	x23, x27, hs
	adds	x0, x0, x25
	ldr	x27, [sp, #216]
	umulh	x2, x14, x17
	str	x0, [sp, #168]
	adcs	x0, x1, x6
	adcs	x1, x4, x20
	str	x27, [sp, #216]
	mul	x13, x14, x18
	adcs	x4, x5, x21
	cinc	x6, x23, hs
	adds	x11, x11, x30
	umulh	x9, x14, x18
	mul	x8, x14, x16
	adcs	x13, x13, x2
	mul	x28, x14, x24
	umulh	x14, x14, x16
	adcs	x9, x8, x9
	umulh	x7, x10, x24
	mul	x19, x10, x17
	cinc	x8, x14, hs
	adds	x14, x0, x28
	umulh	x17, x10, x17
	str	x14, [sp, #160]
	adcs	x11, x1, x11
	adcs	x13, x4, x13
	mul	x5, x10, x18
	adcs	x9, x6, x9
	mov	x6, #977
	cinc	x8, x8, hs
	adds	x0, x19, x7
	movk	x6, #1, lsl #32
	umulh	x18, x10, x18
	mul	x14, x10, x16
	adcs	x17, x5, x17
	mul	x3, x10, x24
	umulh	x16, x10, x16
	adcs	x14, x14, x18
	cinc	x18, x16, hs
	adds	x11, x11, x3
	ldr	x16, [sp, #184]
	str	x11, [sp, #136]
	ldr	x11, [sp, #200]
	ldr	x21, [x16, x11, lsl #3]
	adcs	x11, x13, x0
	adcs	x9, x9, x17
	stp	x11, x9, [sp, #48]
	adcs	x8, x8, x14
	ldp	x9, x14, [x21]
	str	x8, [sp, #64]
	cinc	x11, x18, hs
	str	x11, [sp, #80]
	ldp	x8, x11, [x27, #-32]
	stp	x14, x9, [sp, #144]
	subs	x9, x9, x8
	ldp	x13, x17, [x21, #16]
	stp	x11, x8, [sp, #104]
	sbcs	x16, x14, x11
	ldp	x8, x14, [x27, #-16]
	stp	x14, x8, [sp, #88]
	sbcs	x11, x13, x8
	mul	x8, x24, x22
	stp	x17, x13, [sp, #120]
	sbcs	x13, x17, x14
	str	x8, [sp, #72]
	mov	x8, #977
	movk	x8, #1, lsl #32
	csel	x14, x8, xzr, lo
	subs	x24, x9, x14
	mul	x8, x15, x24
	cset	w9, lo
	subs	x9, x16, x9
	cset	w1, lo
	mul	x14, x9, x22
	subs	x11, x11, x1
	cset	w30, lo
	umulh	x28, x15, x24
	str	x8, [sp, #40]
	sub	x30, x13, x30
	mul	x8, x12, x24
	mul	x19, x15, x9
	umulh	x4, x15, x9
	str	x8, [sp, #32]
	umulh	x8, x12, x24
	mul	x3, x15, x11
	umulh	x5, x15, x11
	str	x8, [sp, #24]
	umulh	x8, x12, x9
	mul	x0, x15, x30
	umulh	x1, x15, x30
	str	x8, [sp, #16]
	umulh	x8, x24, x22
	umulh	x15, x9, x22
	mul	x13, x11, x22
	adds	x8, x14, x8
	mul	x23, x12, x9
	mul	x20, x12, x11
	adcs	x13, x13, x15
	umulh	x26, x12, x11
	mul	x7, x12, x30
	umulh	x25, x12, x30
	mul	x17, x10, x9
	umulh	x18, x10, x9
	umulh	x12, x11, x22
	mul	x9, x30, x22
	mul	x16, x10, x11
	umulh	x2, x10, x11
	adcs	x9, x9, x12
	umulh	x11, x30, x22
	mul	x14, x10, x30
	umulh	x12, x10, x30
	cinc	x11, x11, hs
	adds	x15, x19, x28
	adcs	x3, x3, x4
	mul	x4, x10, x24
	adcs	x0, x0, x5
	ldr	x5, [sp, #40]
	cinc	x1, x1, hs
	umulh	x28, x10, x24
	adds	x8, x8, x5
	adcs	x13, x13, x15
	adcs	x15, x9, x3
	ldr	x9, [sp, #48]
	adcs	x10, x11, x0
	ldr	x11, [sp, #24]
	mul	x3, x9, x6
	cinc	x0, x1, hs
	adds	x1, x23, x11
	ldr	x11, [sp, #56]
	umulh	x23, x9, x6
	ldr	x9, [sp, #16]
	adcs	x19, x20, x9
	ldr	x9, [sp, #32]
	mul	x20, x11, x6
	adcs	x5, x7, x26
	cinc	x7, x25, hs
	adds	x9, x13, x9
	adcs	x13, x15, x1
	umulh	x15, x11, x6
	adcs	x1, x10, x19
	ldr	x10, [sp, #64]
	adcs	x0, x0, x5
	cinc	x7, x7, hs
	adds	x17, x17, x28
	mul	x5, x10, x6
	adcs	x16, x16, x18
	adcs	x11, x14, x2
	ldr	x2, [sp, #80]
	umulh	x18, x10, x6
	cinc	x12, x12, hs
	adds	x10, x13, x4
	mul	x13, x2, x6
	adcs	x14, x1, x17
	adcs	x16, x0, x16
	umulh	x17, x2, x6
	adcs	x11, x7, x11
	cinc	x12, x12, hs
	adds	x0, x20, x23
	ldr	x20, [sp, #136]
	adcs	x15, x5, x15
	ldp	x1, x2, [sp, #160]
	adcs	x13, x13, x18
	ldr	x18, [sp, #72]
	mul	x5, x14, x6
	cinc	x17, x17, hs
	adds	x18, x18, x3
	umulh	x14, x14, x6
	adcs	x0, x2, x0
	adcs	x15, x1, x15
	mul	x19, x16, x6
	cset	w2, hs
	adcs	x3, x20, x13
	add	x13, x20, x13
	cinc	x17, x17, hs
	umulh	x16, x16, x6
	mul	x4, x17, x6
	umulh	x17, x17, x6
	mul	x1, x24, x22
	adds	x18, x18, x4
	mul	x4, x11, x6
	adcs	x17, x0, x17
	cset	w0, hs
	adcs	xzr, x15, xzr
	umulh	x11, x11, x6
	adcs	xzr, x3, xzr
	csel	x3, x6, xzr, hs
	mul	x7, x12, x6
	adds	x22, x18, x3
	adcs	x25, x17, xzr
	umulh	x12, x12, x6
	adcs	x24, x15, x0
	adc	x23, x13, x2
	adds	x13, x19, x14
	adcs	x14, x4, x16
	adcs	x11, x7, x11
	cinc	x12, x12, hs
	adds	x15, x1, x5
	adcs	x8, x8, x13
	adcs	x9, x9, x14
	cset	w13, hs
	adcs	x14, x10, x11
	add	x10, x10, x11
	cinc	x12, x12, hs
	mul	x16, x12, x6
	umulh	x12, x12, x6
	adds	x15, x15, x16
	adcs	x8, x8, x12
	cset	w12, hs
	adcs	xzr, x9, xzr
	adcs	xzr, x14, xzr
	csel	x14, x6, xzr, hs
	adds	x14, x15, x14
	adcs	x8, x8, xzr
	stp	x8, x14, [sp, #160]
	adcs	x8, x9, x12
	adc	x9, x10, x13
	str	x8, [sp, #136]
	ldp	x11, x12, [x27]
	str	x9, [sp, #192]
	ldp	x8, x9, [x21, #32]
	stp	x12, x11, [sp, #72]
	subs	x8, x8, x11
	ldp	x10, x11, [x21, #48]
	sbcs	x9, x9, x12
	ldp	x13, x12, [x27, #16]
	sbcs	x10, x10, x13
	stp	x12, x13, [sp, #56]
	sbcs	x20, x11, x12
	csel	x12, x6, xzr, lo
	subs	x26, x8, x12
	umulh	x12, x26, x22
	cset	w8, lo
	subs	x8, x9, x8
	mul	x0, x26, x25
	cset	w5, lo
	subs	x10, x10, x5
	umulh	x1, x26, x25
	cset	w21, lo
	sub	x20, x20, x21
	mul	x2, x26, x24
	adds	x12, x0, x12
	umulh	x3, x26, x24
	mul	x4, x26, x23
	adcs	x0, x2, x1
	mul	x9, x8, x22
	umulh	x19, x8, x22
	adcs	x2, x4, x3
	mul	x14, x8, x25
	str	x9, [sp, #48]
	umulh	x11, x26, x23
	umulh	x15, x8, x25
	mul	x16, x8, x24
	cinc	x3, x11, hs
	adds	x14, x14, x19
	umulh	x17, x8, x24
	mul	x18, x8, x23
	adcs	x15, x16, x15
	umulh	x13, x8, x23
	umulh	x27, x10, x22
	adcs	x16, x18, x17
	mul	x28, x10, x25
	cinc	x13, x13, hs
	adds	x12, x12, x9
	umulh	x30, x10, x25
	adcs	x14, x0, x14
	adcs	x15, x2, x15
	mul	x8, x10, x24
	adcs	x16, x3, x16
	cinc	x13, x13, hs
	adds	x2, x28, x27
	umulh	x5, x10, x24
	mul	x9, x10, x23
	adcs	x8, x8, x30
	mul	x7, x10, x22
	umulh	x10, x10, x23
	adcs	x9, x9, x5
	mul	x1, x20, x25
	umulh	x0, x20, x22
	cinc	x10, x10, hs
	adds	x11, x14, x7
	umulh	x4, x20, x25
	adcs	x15, x15, x2
	adcs	x8, x16, x8
	mul	x17, x20, x24
	adcs	x9, x13, x9
	cinc	x10, x10, hs
	adds	x13, x1, x0
	umulh	x18, x20, x24
	mul	x3, x20, x23
	adcs	x17, x17, x4
	umulh	x14, x20, x23
	mul	x16, x20, x22
	adcs	x18, x3, x18
	mul	x0, x26, x22
	cinc	x14, x14, hs
	adds	x15, x15, x16
	adcs	x8, x8, x13
	adcs	x9, x9, x17
	mul	x13, x8, x6
	adcs	x10, x10, x18
	umulh	x8, x8, x6
	cinc	x14, x14, hs
	mul	x16, x9, x6
	umulh	x9, x9, x6
	mul	x17, x10, x6
	adds	x8, x16, x8
	umulh	x10, x10, x6
	mul	x18, x14, x6
	adcs	x9, x17, x9
	umulh	x14, x14, x6
	adcs	x10, x18, x10
	cinc	x14, x14, hs
	adds	x13, x0, x13
	adcs	x8, x12, x8
	adcs	x9, x11, x9
	cset	w11, hs
	adcs	x12, x15, x10
	add	x10, x15, x10
	cinc	x14, x14, hs
	mul	x16, x14, x6
	umulh	x14, x14, x6
	adds	x13, x13, x16
	adcs	x8, x8, x14
	cset	w14, hs
	adcs	xzr, x9, xzr
	adcs	xzr, x12, xzr
	csel	x12, x6, xzr, hs
	adds	x22, x13, x12
	adcs	x25, x8, xzr
	umulh	x2, x22, x22
	adcs	x24, x9, x14
	umulh	x12, x25, x22
	adc	x23, x10, x11
	mul	x9, x24, x22
	umulh	x13, x24, x22
	mul	x10, x23, x22
	adds	x9, x12, x9
	umulh	x11, x23, x22
	umulh	x15, x24, x25
	adcs	x10, x13, x10
	mul	x16, x23, x25
	cinc	x11, x11, hs
	mul	x14, x24, x25
	umulh	x17, x23, x25
	adds	x13, x15, x16
	mul	x8, x25, x22
	mul	x18, x23, x24
	cinc	x15, x17, hs
	adds	x10, x10, x14
	umulh	x0, x23, x24
	adcs	x11, x11, x13
	extr	x17, x9, x8, #63
	extr	x9, x10, x9, #63
	extr	x10, x11, x10, #63
	mul	x3, x25, x25
	adcs	x13, x15, x18
	umulh	x12, x25, x25
	extr	x11, x13, x11, #63
	cinc	x16, x0, hs
	adds	x8, x2, x8, lsl #1
	mul	x14, x24, x24
	extr	x13, x16, x13, #63
	lsr	x16, x16, #63
	adcs	x17, x17, x3
	umulh	x15, x24, x24
	adcs	x9, x9, x12
	mul	x18, x23, x23
	adcs	x10, x10, x14
	umulh	x12, x23, x23
	adcs	x11, x11, x15
	mul	x14, x10, x6
	adcs	x13, x13, x18
	umulh	x10, x10, x6
	adc	x12, x16, x12
	mul	x15, x11, x6
	umulh	x11, x11, x6
	mul	x18, x13, x6
	adds	x10, x15, x10
	umulh	x13, x13, x6
	mul	x16, x12, x6
	adcs	x11, x18, x11
	ldp	x0, x18, [sp, #88]
	mul	x1, x22, x22
	umulh	x12, x12, x6
	adcs	x13, x16, x13
	cinc	x12, x12, hs
	adds	x14, x1, x14
	adcs	x8, x8, x10
	adcs	x10, x17, x11
	cset	w11, hs
	adcs	x15, x9, x13
	add	x9, x9, x13
	cinc	x12, x12, hs
	mul	x16, x12, x6
	umulh	x12, x12, x6
	adds	x14, x14, x16
	ldp	x17, x16, [sp, #104]
	adcs	x8, x8, x12
	cset	w12, hs
	adcs	xzr, x10, xzr
	adcs	xzr, x15, xzr
	csel	x15, x6, xzr, hs
	adds	x14, x14, x15
	adcs	x8, x8, xzr
	adcs	x10, x10, x12
	adc	x9, x9, x11
	subs	x11, x14, x16
	sbcs	x8, x8, x17
	ldp	x13, x14, [sp, #144]
	sbcs	x10, x10, x18
	sbcs	x9, x9, x0
	csel	x12, x6, xzr, lo
	subs	x11, x11, x12
	cset	w12, lo
	subs	x8, x8, x12
	cset	w12, lo
	subs	x10, x10, x12
	cset	w12, lo
	subs	x11, x11, x14
	sbcs	x8, x8, x13
	sub	x9, x9, x12
	ldp	x12, x13, [sp, #120]
	sbcs	x10, x10, x13
	sbcs	x9, x9, x12
	csel	x12, x6, xzr, lo
	subs	x11, x11, x12
	cset	w12, lo
	subs	x8, x8, x12
	cset	w12, lo
	subs	x10, x10, x12
	cset	w12, lo
	adds	x13, x11, x6
	adcs	x14, x8, xzr
	sub	x9, x9, x12
	adcs	x12, x10, xzr
	adcs	x15, x9, xzr
	csel	x13, x13, x11, hs
	csel	x11, x12, x10, hs
	csel	x10, x14, x8, hs
	csel	x12, x15, x9, hs
	subs	x8, x16, x13
	stp	x13, x10, [sp, #120]
	sbcs	x9, x17, x10
	stp	x11, x12, [sp, #144]
	sbcs	x10, x18, x11
	sbcs	x11, x0, x12
	csel	x12, x6, xzr, lo
	subs	x16, x8, x12
	cset	w8, lo
	umulh	x18, x16, x22
	subs	x8, x9, x8
	cset	w13, lo
	mul	x15, x8, x22
	subs	x17, x10, x13
	cset	w7, lo
	umulh	x14, x8, x22
	sub	x11, x11, x7
	mul	x13, x17, x22
	adds	x15, x15, x18
	umulh	x4, x17, x22
	mul	x21, x11, x22
	adcs	x13, x13, x14
	mul	x5, x25, x16
	umulh	x3, x25, x16
	adcs	x18, x21, x4
	mul	x0, x25, x8
	umulh	x1, x25, x8
	mul	x2, x25, x17
	umulh	x19, x25, x17
	mul	x7, x25, x11
	umulh	x20, x25, x11
	umulh	x25, x11, x22
	umulh	x27, x24, x16
	mul	x28, x24, x8
	cinc	x4, x25, hs
	adds	x0, x0, x3
	adcs	x1, x2, x1
	umulh	x30, x24, x8
	adcs	x2, x7, x19
	mul	x9, x24, x17
	cinc	x7, x20, hs
	adds	x15, x15, x5
	adcs	x13, x13, x0
	umulh	x12, x24, x17
	adcs	x18, x18, x1
	adcs	x0, x4, x2
	mul	x10, x24, x11
	cinc	x2, x7, hs
	adds	x4, x28, x27
	adcs	x9, x9, x30
	mul	x26, x24, x16
	umulh	x24, x24, x11
	adcs	x10, x10, x12
	mul	x14, x23, x8
	umulh	x1, x23, x16
	cinc	x12, x24, hs
	adds	x13, x13, x26
	umulh	x8, x23, x8
	adcs	x18, x18, x4
	adcs	x9, x0, x9
	mul	x3, x23, x17
	adcs	x10, x2, x10
	cinc	x12, x12, hs
	adds	x14, x14, x1
	umulh	x17, x23, x17
	mul	x5, x23, x11
	adcs	x8, x3, x8
	umulh	x11, x23, x11
	mul	x0, x23, x16
	adcs	x17, x5, x17
	mul	x16, x16, x22
	cinc	x11, x11, hs
	adds	x18, x18, x0
	adcs	x9, x9, x14
	adcs	x8, x10, x8
	mul	x14, x9, x6
	adcs	x12, x12, x17
	umulh	x9, x9, x6
	cinc	x11, x11, hs
	mul	x10, x8, x6
	umulh	x8, x8, x6
	mul	x17, x12, x6
	adds	x9, x10, x9
	umulh	x12, x12, x6
	mul	x0, x11, x6
	adcs	x8, x17, x8
	umulh	x11, x11, x6
	adcs	x10, x0, x12
	cinc	x11, x11, hs
	adds	x12, x16, x14
	ldr	x16, [sp, #216]
	adcs	x9, x15, x9
	adcs	x8, x13, x8
	cset	w13, hs
	adcs	x14, x18, x10
	add	x10, x18, x10
	cinc	x11, x11, hs
	mul	x15, x11, x6
	umulh	x11, x11, x6
	adds	x12, x12, x15
	ldr	x15, [sp, #120]
	adcs	x9, x9, x11
	cset	w11, hs
	adcs	xzr, x8, xzr
	adcs	xzr, x14, xzr
	csel	x14, x6, xzr, hs
	adds	x12, x12, x14
	adcs	x9, x9, xzr
	adcs	x8, x8, x11
	ldr	x11, [sp, #80]
	adc	x10, x10, x13
	subs	x11, x12, x11
	ldp	x12, x13, [sp, #64]
	sbcs	x9, x9, x13
	sbcs	x8, x8, x12
	ldr	x12, [sp, #56]
	sbcs	x10, x10, x12
	csel	x12, x6, xzr, lo
	subs	x11, x11, x12
	cset	w12, lo
	subs	x9, x9, x12
	cset	w12, lo
	subs	x8, x8, x12
	cset	w12, lo
	adds	x13, x11, x6
	sub	x10, x10, x12
	ldr	x12, [sp, #128]
	adcs	x14, x9, xzr
	stp	x15, x12, [x16, #-32]
	ldr	x15, [sp, #144]
	adcs	x12, x8, xzr
	stur	x15, [x16, #-16]
	adcs	x15, x10, xzr
	csel	x11, x13, x11, hs
	ldr	x13, [sp, #152]
	csel	x9, x14, x9, hs
	csel	x8, x12, x8, hs
	ldp	x12, x2, [sp, #176]
	stp	x9, x8, [x16, #8]
	csel	x9, x15, x10, hs
	ldr	x14, [sp, #136]
	stp	x13, x11, [x16, #-8]
	ldr	x11, [sp, #200]
	str	x9, [x16, #24]
	ldp	x15, x9, [sp, #160]
	sub	x11, x11, #1
	sub	x16, x16, #64
	sub	x12, x12, #32
	str	x16, [sp, #216]
	cbnz	x11, .LBB1_2
.LBB1_3:
	ldr	x12, [x2]
	ldr	x27, [sp, #8]
	ldp	x13, x18, [x12, #32]
	ldp	x8, x10, [x27, #32]
	ldp	x0, x3, [x12, #48]
	str	x8, [sp, #216]
	subs	x2, x13, x8
	ldp	x8, x11, [x27, #48]
	sbcs	x4, x18, x10
	str	x10, [sp, #184]
	stp	x12, x8, [sp, #168]
	sbcs	x0, x0, x8
	mov	x8, #977
	sbcs	x10, x3, x11
	movk	x8, #1, lsl #32
	str	x11, [sp, #200]
	csel	x18, x8, xzr, lo
	subs	x18, x2, x18
	umulh	x8, x18, x14
	cset	w21, lo
	subs	x4, x4, x21
	umulh	x2, x18, x9
	cset	w30, lo
	subs	x0, x0, x30
	mul	x5, x18, x15
	cset	w3, lo
	stp	x8, x10, [sp, #152]
	ldr	x10, [sp, #192]
	umulh	x6, x18, x15
	ldr	x20, [sp, #160]
	mul	x8, x18, x10
	sub	x3, x20, x3
	adds	x2, x5, x2
	mul	x7, x18, x14
	umulh	x22, x18, x10
	str	x8, [sp, #144]
	mov	x8, #977
	umulh	x23, x4, x9
	movk	x8, #1, lsl #32
	adcs	x5, x7, x6
	ldp	x20, x7, [sp, #144]
	mul	x24, x4, x15
	umulh	x25, x4, x15
	adcs	x7, x20, x7
	cinc	x20, x22, hs
	mul	x26, x4, x14
	adds	x22, x24, x23
	umulh	x19, x4, x14
	mul	x28, x4, x10
	adcs	x24, x26, x25
	mul	x21, x4, x9
	umulh	x4, x4, x10
	adcs	x25, x28, x19
	umulh	x1, x0, x9
	mul	x13, x0, x15
	cinc	x4, x4, hs
	adds	x2, x2, x21
	umulh	x17, x0, x15
	adcs	x5, x5, x22
	adcs	x7, x7, x24
	mul	x11, x0, x14
	adcs	x20, x20, x25
	cinc	x4, x4, hs
	adds	x13, x13, x1
	umulh	x16, x0, x14
	mul	x12, x0, x10
	adcs	x11, x11, x17
	mul	x30, x0, x9
	umulh	x0, x0, x10
	adcs	x12, x12, x16
	umulh	x23, x3, x9
	mul	x26, x3, x15
	cinc	x16, x0, hs
	adds	x17, x5, x30
	umulh	x15, x3, x15
	adcs	x13, x7, x13
	adcs	x11, x20, x11
	mul	x21, x3, x14
	adcs	x12, x4, x12
	cinc	x16, x16, hs
	adds	x1, x26, x23
	umulh	x14, x3, x14
	mul	x0, x3, x10
	adcs	x15, x21, x15
	mul	x6, x3, x9
	umulh	x10, x3, x10
	adcs	x14, x0, x14
	mul	x9, x18, x9
	cinc	x10, x10, hs
	adds	x13, x13, x6
	adcs	x11, x11, x1
	adcs	x12, x12, x15
	mul	x18, x11, x8
	adcs	x14, x16, x14
	umulh	x11, x11, x8
	cinc	x10, x10, hs
	mul	x15, x12, x8
	umulh	x12, x12, x8
	mul	x16, x14, x8
	adds	x11, x15, x11
	umulh	x14, x14, x8
	mul	x0, x10, x8
	adcs	x12, x16, x12
	umulh	x10, x10, x8
	adcs	x16, x0, x14
	cinc	x10, x10, hs
	adds	x9, x9, x18
	adcs	x11, x2, x11
	adcs	x12, x17, x12
	cset	w17, hs
	adcs	x14, x13, x16
	cinc	x10, x10, hs
	mul	x15, x10, x8
	umulh	x10, x10, x8
	adds	x9, x9, x15
	adcs	x10, x11, x10
	cset	w11, hs
	adcs	xzr, x12, xzr
	adcs	xzr, x14, xzr
	csel	x14, x8, xzr, hs
	adds	x9, x9, x14
	adcs	x15, x10, xzr
	add	x10, x13, x16
	umulh	x7, x9, x9
	adcs	x14, x12, x11
	umulh	x0, x15, x9
	adc	x10, x10, x17
	mul	x11, x14, x9
	umulh	x12, x14, x9
	mul	x13, x10, x9
	adds	x11, x0, x11
	umulh	x16, x10, x9
	umulh	x1, x14, x15
	adcs	x12, x12, x13
	mul	x2, x10, x15
	cinc	x16, x16, hs
	mul	x17, x14, x15
	umulh	x3, x10, x15
	adds	x0, x1, x2
	mul	x18, x15, x9
	mul	x4, x10, x14
	cinc	x1, x3, hs
	adds	x12, x12, x17
	umulh	x5, x10, x14
	adcs	x16, x16, x0
	extr	x3, x11, x18, #63
	extr	x11, x12, x11, #63
	extr	x12, x16, x12, #63
	mul	x20, x15, x15
	adcs	x0, x1, x4
	umulh	x13, x15, x15
	extr	x16, x0, x16, #63
	cinc	x2, x5, hs
	adds	x18, x7, x18, lsl #1
	mul	x17, x14, x14
	extr	x0, x2, x0, #63
	lsr	x2, x2, #63
	adcs	x3, x3, x20
	umulh	x1, x14, x14
	adcs	x11, x11, x13
	mul	x4, x10, x10
	adcs	x12, x12, x17
	umulh	x5, x10, x10
	adcs	x16, x16, x1
	mul	x13, x12, x8
	adcs	x0, x0, x4
	umulh	x12, x12, x8
	adc	x2, x2, x5
	mul	x17, x16, x8
	umulh	x16, x16, x8
	mul	x1, x0, x8
	adds	x12, x17, x12
	umulh	x0, x0, x8
	mul	x4, x2, x8
	adcs	x16, x1, x16
	mul	x6, x9, x9
	umulh	x2, x2, x8
	adcs	x17, x4, x0
	ldp	x4, x5, [x27, #16]
	cinc	x0, x2, hs
	adds	x13, x6, x13
	adcs	x12, x18, x12
	adcs	x16, x3, x16
	cset	w18, hs
	adcs	x1, x11, x17
	add	x11, x11, x17
	cinc	x0, x0, hs
	mul	x2, x0, x8
	umulh	x0, x0, x8
	adds	x13, x13, x2
	ldr	x2, [sp, #168]
	adcs	x12, x12, x0
	cset	w0, hs
	adcs	xzr, x16, xzr
	adcs	xzr, x1, xzr
	csel	x1, x8, xzr, hs
	adds	x13, x13, x1
	ldp	x1, x3, [x27]
	adcs	x12, x12, xzr
	adcs	x16, x16, x0
	adc	x11, x11, x18
	subs	x13, x13, x1
	sbcs	x12, x12, x3
	sbcs	x16, x16, x4
	sbcs	x11, x11, x5
	csel	x17, x8, xzr, lo
	subs	x13, x13, x17
	cset	w17, lo
	subs	x12, x12, x17
	cset	w17, lo
	subs	x16, x16, x17
	ldp	x18, x17, [x2]
	cset	w0, lo
	sub	x11, x11, x0
	subs	x13, x13, x18
	sbcs	x12, x12, x17
	ldp	x18, x17, [x2, #16]
	sbcs	x16, x16, x18
	sbcs	x11, x11, x17
	csel	x17, x8, xzr, lo
	subs	x13, x13, x17
	cset	w17, lo
	subs	x12, x12, x17
	cset	w17, lo
	subs	x16, x16, x17
	cset	w17, lo
	adds	x18, x13, x8
	adcs	x0, x12, xzr
	sub	x11, x11, x17
	adcs	x17, x16, xzr
	adcs	x6, x11, xzr
	csel	x13, x18, x13, hs
	csel	x16, x17, x16, hs
	csel	x12, x0, x12, hs
	csel	x17, x6, x11, hs
	subs	x11, x1, x13
	stp	x13, x12, [sp, #152]
	sbcs	x12, x3, x12
	str	x16, [sp, #168]
	sbcs	x13, x4, x16
	str	x17, [sp, #192]
	sbcs	x16, x5, x17
	csel	x1, x8, xzr, lo
	subs	x3, x11, x1
	cset	w11, lo
	umulh	x1, x3, x9
	subs	x11, x12, x11
	cset	w4, lo
	mul	x18, x11, x9
	subs	x5, x13, x4
	cset	w4, lo
	umulh	x6, x11, x9
	sub	x4, x16, x4
	mul	x17, x5, x9
	adds	x18, x18, x1
	umulh	x7, x5, x9
	mul	x2, x4, x9
	adcs	x17, x17, x6
	umulh	x19, x4, x9
	umulh	x21, x15, x3
	adcs	x2, x2, x7
	mul	x22, x15, x11
	cinc	x6, x19, hs
	umulh	x23, x15, x11
	mul	x24, x15, x5
	adds	x7, x22, x21
	umulh	x25, x15, x5
	mul	x26, x15, x4
	adcs	x19, x24, x23
	ldp	x24, x23, [sp, #256]
	mul	x20, x15, x3
	umulh	x0, x15, x4
	adcs	x21, x26, x25
	ldp	x26, x25, [sp, #240]
	mul	x12, x14, x3
	umulh	x28, x14, x3
	cinc	x0, x0, hs
	adds	x18, x18, x20
	mul	x30, x14, x11
	adcs	x17, x17, x7
	str	x12, [sp, #144]
	adcs	x2, x2, x19
	umulh	x16, x14, x11
	adcs	x6, x6, x21
	cinc	x0, x0, hs
	mul	x12, x14, x5
	adds	x20, x30, x28
	ldr	x30, [sp, #208]
	umulh	x15, x14, x5
	mul	x13, x14, x4
	adcs	x12, x12, x16
	umulh	x14, x14, x4
	mul	x1, x10, x11
	adcs	x13, x13, x15
	ldr	x15, [sp, #144]
	umulh	x7, x10, x3
	cinc	x14, x14, hs
	adds	x15, x17, x15
	umulh	x11, x10, x11
	adcs	x2, x2, x20
	adcs	x12, x6, x12
	mul	x19, x10, x5
	adcs	x13, x0, x13
	cinc	x14, x14, hs
	adds	x0, x1, x7
	umulh	x16, x10, x5
	mul	x17, x10, x4
	adcs	x11, x19, x11
	ldp	x20, x19, [sp, #288]
	mul	x22, x10, x3
	umulh	x10, x10, x4
	adcs	x16, x17, x16
	mul	x9, x3, x9
	cinc	x10, x10, hs
	adds	x17, x2, x22
	adcs	x12, x12, x0
	ldp	x22, x21, [sp, #272]
	adcs	x11, x13, x11
	mul	x0, x12, x8
	adcs	x14, x14, x16
	umulh	x12, x12, x8
	cinc	x10, x10, hs
	mul	x13, x11, x8
	umulh	x11, x11, x8
	mul	x16, x14, x8
	adds	x12, x13, x12
	umulh	x14, x14, x8
	mul	x1, x10, x8
	adcs	x11, x16, x11
	umulh	x10, x10, x8
	adcs	x13, x1, x14
	cinc	x10, x10, hs
	adds	x9, x9, x0
	adcs	x12, x18, x12
	adcs	x11, x15, x11
	cset	w14, hs
	adcs	x15, x17, x13
	add	x13, x17, x13
	cinc	x10, x10, hs
	mul	x16, x10, x8
	umulh	x10, x10, x8
	adds	x9, x9, x16
	adcs	x10, x12, x10
	cset	w12, hs
	adcs	xzr, x11, xzr
	adcs	xzr, x15, xzr
	csel	x15, x8, xzr, hs
	adds	x9, x9, x15
	ldr	x15, [sp, #152]
	adcs	x10, x10, xzr
	adcs	x11, x11, x12
	str	x15, [x27]
	adc	x12, x13, x14
	ldr	x13, [sp, #216]
	ldp	x15, x16, [sp, #160]
	subs	x9, x9, x13
	ldp	x13, x14, [sp, #176]
	stp	x15, x16, [x27, #8]
	sbcs	x10, x10, x14
	sbcs	x11, x11, x13
	ldr	x13, [sp, #200]
	sbcs	x12, x12, x13
	csel	x13, x8, xzr, lo
	subs	x9, x9, x13
	cset	w13, lo
	subs	x10, x10, x13
	cset	w13, lo
	subs	x11, x11, x13
	cset	w13, lo
	adds	x8, x9, x8
	adcs	x14, x10, xzr
	sub	x12, x12, x13
	adcs	x13, x11, xzr
	adcs	x15, x12, xzr
	csel	x8, x8, x9, hs
	ldr	x9, [sp, #192]
	stp	x9, x8, [x27, #24]
	csel	x9, x14, x10, hs
	csel	x8, x13, x11, hs
	stp	x9, x8, [x27, #40]
	csel	x9, x15, x12, hs
	str	x9, [x27, #56]
	ldp	x28, x27, [sp, #224]
	add	sp, sp, #304
	ret

pass1:
	sub	sp, sp, #112
	ldr	x8, [x2]
	stp	x20, x19, [sp, #96]
	mov	x19, #977
	ldp	x9, x10, [x1]
	movk	x19, #1, lsl #32
	ldp	x11, x12, [x8]
	str	x30, [sp, #16]
	ldp	x13, x14, [x8, #16]
	stp	x28, x27, [sp, #32]
	stp	x26, x25, [sp, #48]
	ldp	x15, x16, [x1, #16]
	subs	x17, x11, x9
	sbcs	x18, x12, x10
	sbcs	x4, x13, x15
	sbcs	x5, x14, x16
	csel	x6, x19, xzr, lo
	subs	x17, x17, x6
	sbcs	x18, x18, xzr
	sbcs	x4, x4, xzr
	sbc	x5, x5, xzr
	stp	x24, x23, [sp, #64]
	cmp	x0, #2
	stp	x22, x21, [sp, #80]
	stp	x17, x18, [x3]
	stp	x4, x5, [x3, #16]
	b.lo	.LBB0_3
	sub	x10, x0, #1
	add	x8, x2, #8
	add	x11, x1, #88
	add	x12, x3, #32
.LBB0_2:
	ldr	x13, [x8], #8
	str	x8, [sp, #8]
	ldp	x3, x4, [x11, #-24]
	str	x10, [sp, #24]
	ldp	x18, x0, [x13]
	ldp	x1, x2, [x13, #16]
	ldp	x5, x6, [x11, #-8]
	subs	x14, x18, x3
	sbcs	x17, x0, x4
	sbcs	x16, x1, x5
	sbcs	x15, x2, x6
	csel	x7, x19, xzr, lo
	subs	x14, x14, x7
	sbcs	x17, x17, xzr
	sbcs	x16, x16, xzr
	sbc	x15, x15, xzr
	ldp	x13, x18, [x12, #-32]
	add	x11, x11, #64
	ldp	x25, x26, [x12, #-16]
	umulh	x7, x13, x14
	mul	x1, x13, x17
	umulh	x2, x13, x17
	mul	x3, x13, x16
	adds	x1, x1, x7
	umulh	x4, x13, x16
	mul	x5, x13, x15
	adcs	x2, x3, x2
	umulh	x6, x13, x15
	mul	x8, x18, x14
	adcs	x4, x5, x4
	umulh	x19, x18, x14
	cinc	x5, x6, hs
	mul	x20, x18, x17
	str	x8, [sp]
	umulh	x21, x18, x17
	mul	x22, x18, x16
	adds	x6, x20, x19
	umulh	x23, x18, x16
	mul	x24, x18, x15
	adcs	x19, x22, x21
	umulh	x0, x18, x15
	umulh	x28, x25, x14
	adcs	x20, x24, x23
	mul	x30, x25, x17
	cinc	x0, x0, hs
	adds	x1, x1, x8
	umulh	x18, x25, x17
	adcs	x2, x2, x6
	adcs	x4, x4, x19
	mul	x9, x25, x16
	adcs	x5, x5, x20
	cinc	x0, x0, hs
	adds	x19, x30, x28
	umulh	x10, x25, x16
	mul	x8, x25, x15
	adcs	x9, x9, x18
	mul	x27, x25, x14
	umulh	x25, x25, x15
	adcs	x8, x8, x10
	umulh	x7, x26, x14
	mul	x21, x26, x17
	cinc	x10, x25, hs
	adds	x18, x2, x27
	umulh	x17, x26, x17
	adcs	x4, x4, x19
	mov	x19, #977
	adcs	x9, x5, x9
	movk	x19, #1, lsl #32
	mul	x6, x26, x16
	adcs	x8, x0, x8
	cinc	x10, x10, hs
	adds	x0, x21, x7
	umulh	x16, x26, x16
	mul	x2, x26, x15
	adcs	x17, x6, x17
	mul	x3, x26, x14
	umulh	x15, x26, x15
	adcs	x16, x2, x16
	mul	x13, x13, x14
	cinc	x14, x15, hs
	adds	x15, x4, x3
	adcs	x9, x9, x0
	adcs	x8, x8, x17
	mul	x0, x9, x19
	adcs	x10, x10, x16
	umulh	x9, x9, x19
	cinc	x14, x14, hs
	mul	x17, x8, x19
	umulh	x8, x8, x19
	mul	x16, x10, x19
	adds	x9, x17, x9
	umulh	x10, x10, x19
	mul	x2, x14, x19
	adcs	x8, x16, x8
	umulh	x14, x14, x19
	adcs	x10, x2, x10
	cinc	x14, x14, hs
	adds	x13, x13, x0
	adcs	x9, x1, x9
	adcs	x8, x18, x8
	cset	w16, hs
	adcs	x17, x15, x10
	add	x10, x15, x10
	cinc	x14, x14, hs
	mul	x18, x14, x19
	umulh	x14, x14, x19
	adds	x13, x13, x18
	adcs	x9, x9, x14
	cset	w14, hs
	adcs	xzr, x8, xzr
	adcs	xzr, x17, xzr
	csel	x17, x19, xzr, hs
	adds	x13, x13, x17
	adcs	x9, x9, xzr
	adcs	x8, x8, x14
	stp	x13, x9, [x12]
	adc	x9, x10, x16
	ldr	x10, [sp, #24]
	stp	x8, x9, [x12, #16]
	ldr	x8, [sp, #8]
	add	x12, x12, #32
	subs	x10, x10, #1
	b.ne	.LBB0_2
.LBB0_3:
	ldp	x20, x19, [sp, #96]
	ldr	x30, [sp, #16]
	ldp	x22, x21, [sp, #80]
	ldp	x24, x23, [sp, #64]
	ldp	x26, x25, [sp, #48]
	ldp	x28, x27, [sp, #32]
	add	sp, sp, #112
	ret
pass2:
	sub	sp, sp, #304
	ldp	x9, x15, [x4]
	subs	x10, x0, #1
	ldp	x14, x13, [x4, #16]
	str	x30, [sp, #208]
	stp	x28, x27, [sp, #224]
	stp	x26, x25, [sp, #240]
	stp	x24, x23, [sp, #256]
	stp	x22, x21, [sp, #272]
	stp	x20, x19, [sp, #288]
	str	x1, [sp, #8]
	str	x13, [sp, #192]
	b.eq	.LBB1_3
	add	x11, x1, x0, lsl #6
	add	x12, x3, x0, lsl #5
	str	x2, [sp, #184]
	sub	x8, x11, #32
	sub	x11, x12, #64
	str	x8, [sp, #216]
.LBB1_2:
	ldp	x22, x17, [x11]
	mov	x7, x9
	ldp	x18, x16, [x11, #16]
	str	x11, [sp, #176]
	str	x10, [sp, #200]
	ldr	x10, [sp, #192]
	umulh	x0, x22, x9
	mul	x1, x17, x9
	umulh	x3, x17, x9
	mul	x4, x18, x9
	adds	x0, x1, x0
	umulh	x5, x18, x9
	mul	x6, x16, x9
	adcs	x1, x4, x3
	umulh	x11, x16, x9
	mul	x8, x15, x22
	adcs	x4, x6, x5
	umulh	x20, x15, x22
	cinc	x5, x11, hs
	mul	x21, x15, x17
	str	x8, [sp, #168]
	umulh	x23, x15, x17
	mul	x24, x15, x18
	adds	x6, x21, x20
	mov	x20, x8
	umulh	x25, x15, x18
	mul	x26, x15, x16
	adcs	x21, x24, x23
	umulh	x27, x15, x16
	umulh	x30, x14, x22
	adcs	x23, x26, x25
	mul	x8, x14, x17
	cinc	x25, x27, hs
	adds	x0, x0, x20
	umulh	x19, x14, x17
	str	x0, [sp, #104]
	adcs	x0, x1, x6
	adcs	x1, x4, x21
	mul	x12, x14, x18
	adcs	x4, x5, x23
	cinc	x6, x25, hs
	adds	x8, x8, x30
	umulh	x9, x14, x18
	mul	x13, x14, x16
	adcs	x12, x12, x19
	mul	x28, x14, x22
	umulh	x2, x14, x16
	adcs	x9, x13, x9
	umulh	x11, x10, x22
	mul	x24, x10, x17
	cinc	x13, x2, hs
	adds	x0, x0, x28
	umulh	x17, x10, x17
	adcs	x8, x1, x8
	str	x0, [sp, #96]
	adcs	x12, x4, x12
	ldr	x1, [sp, #184]
	mul	x5, x10, x18
	adcs	x9, x6, x9
	cinc	x13, x13, hs
	adds	x0, x24, x11
	ldr	x11, [sp, #200]
	umulh	x18, x10, x18
	ldr	x21, [x1, x11, lsl #3]
	mul	x2, x10, x16
	adcs	x17, x5, x17
	mul	x3, x10, x22
	umulh	x16, x10, x16
	adcs	x18, x2, x18
	cinc	x1, x16, hs
	adds	x8, x8, x3
	str	x8, [sp, #88]
	adcs	x8, x12, x0
	str	x8, [sp, #48]
	adcs	x8, x9, x17
	mul	x9, x22, x7
	str	x8, [sp, #64]
	ldr	x8, [sp, #216]
	ldp	x0, x2, [x21]
	ldp	x5, x6, [x8, #-32]
	ldp	x11, x17, [x8, #-16]
	adcs	x8, x13, x18
	str	x8, [sp, #56]
	cinc	x8, x1, hs
	mov	x1, #977
	stp	x8, x9, [sp, #72]
	mov	x8, #977
	movk	x1, #1, lsl #32
	movk	x8, #1, lsl #32
	ldp	x3, x4, [x21, #16]
	stp	x11, x17, [sp, #112]
	subs	x23, x0, x5
	sbcs	x9, x2, x6
	sbcs	x16, x3, x11
	sbcs	x28, x4, x17
	csel	x12, x8, xzr, lo
	subs	x23, x23, x12
	sbcs	x9, x9, xzr
	sbcs	x16, x16, xzr
	sbc	x28, x28, xzr
	mul	x11, x15, x23
	stp	x5, x6, [sp, #128]
	mul	x8, x14, x23
	stp	x3, x0, [sp, #144]
	stp	x4, x2, [sp, #160]
	mul	x13, x9, x7
	umulh	x6, x15, x23
	stp	x11, x8, [sp, #32]
	umulh	x11, x14, x23
	umulh	x8, x14, x9
	mul	x3, x15, x9
	umulh	x4, x15, x9
	stp	x8, x11, [sp, #16]
	umulh	x8, x23, x7
	mul	x2, x15, x16
	umulh	x30, x15, x16
	adds	x8, x13, x8
	mul	x0, x15, x28
	umulh	x5, x15, x28
	mul	x22, x14, x9
	mul	x20, x14, x16
	umulh	x27, x14, x16
	mul	x19, x14, x28
	umulh	x26, x14, x28
	umulh	x14, x9, x7
	mul	x15, x16, x7
	mul	x18, x10, x9
	umulh	x24, x10, x9
	adcs	x13, x15, x14
	ldr	x14, [sp, #48]
	umulh	x12, x16, x7
	mul	x9, x28, x7
	umulh	x11, x28, x7
	mul	x17, x10, x16
	adcs	x9, x9, x12
	umulh	x25, x10, x16
	cinc	x12, x11, hs
	adds	x3, x3, x6
	ldr	x11, [sp, #32]
	adcs	x2, x2, x4
	mul	x16, x10, x28
	adcs	x0, x0, x30
	cinc	x6, x5, hs
	adds	x11, x8, x11
	umulh	x15, x10, x28
	adcs	x8, x13, x3
	adcs	x9, x9, x2
	mul	x4, x10, x23
	umulh	x5, x10, x23
	adcs	x10, x12, x0
	ldr	x12, [sp, #24]
	cinc	x13, x6, hs
	adds	x0, x22, x12
	ldr	x12, [sp, #16]
	mul	x2, x14, x1
	umulh	x22, x14, x1
	adcs	x3, x20, x12
	ldr	x14, [sp, #40]
	adcs	x6, x19, x27
	ldr	x12, [sp, #64]
	cinc	x19, x26, hs
	adds	x8, x8, x14
	adcs	x9, x9, x0
	mul	x20, x12, x1
	adcs	x10, x10, x3
	adcs	x3, x13, x6
	umulh	x0, x12, x1
	ldr	x12, [sp, #56]
	cinc	x19, x19, hs
	adds	x18, x18, x5
	adcs	x14, x17, x24
	mul	x6, x12, x1
	umulh	x17, x12, x1
	adcs	x12, x16, x25
	ldr	x16, [sp, #72]
	cinc	x15, x15, hs
	adds	x13, x9, x4
	mul	x9, x16, x1
	adcs	x10, x10, x18
	adcs	x14, x3, x14
	umulh	x16, x16, x1
	adcs	x12, x19, x12
	cinc	x15, x15, hs
	adds	x18, x20, x22
	adcs	x0, x6, x0
	mul	x20, x23, x7
	adcs	x9, x9, x17
	ldp	x17, x7, [sp, #80]
	mul	x5, x10, x1
	cinc	x16, x16, hs
	adds	x17, x17, x2
	ldp	x2, x3, [sp, #96]
	umulh	x10, x10, x1
	mul	x6, x14, x1
	adcs	x18, x3, x18
	adcs	x0, x2, x0
	umulh	x14, x14, x1
	cset	w2, hs
	adcs	x3, x7, x9
	add	x9, x7, x9
	cinc	x16, x16, hs
	mul	x19, x15, x1
	mul	x4, x16, x1
	umulh	x16, x16, x1
	umulh	x15, x15, x1
	adds	x17, x17, x4
	mul	x4, x12, x1
	adcs	x16, x18, x16
	cset	w18, hs
	adcs	xzr, x0, xzr
	umulh	x12, x12, x1
	adcs	xzr, x3, xzr
	csel	x3, x1, xzr, hs
	adds	x22, x17, x3
	adcs	x25, x16, xzr
	adcs	x24, x0, x18
	adc	x23, x9, x2
	adds	x9, x6, x10
	adcs	x10, x4, x14
	adcs	x12, x19, x12
	cinc	x14, x15, hs
	adds	x15, x20, x5
	adcs	x9, x11, x9
	adcs	x10, x8, x10
	cset	w11, hs
	adcs	x8, x13, x12
	add	x12, x13, x12
	cinc	x14, x14, hs
	mul	x16, x14, x1
	umulh	x14, x14, x1
	adds	x15, x15, x16
	ldp	x16, x17, [x21, #48]
	adcs	x9, x9, x14
	cset	w14, hs
	adcs	xzr, x10, xzr
	adcs	xzr, x8, xzr
	csel	x8, x1, xzr, hs
	adds	x8, x15, x8
	str	x8, [sp, #104]
	adcs	x8, x9, xzr
	str	x8, [sp, #96]
	ldr	x8, [sp, #216]
	ldp	x9, x15, [x21, #32]
	ldp	x18, x0, [x8]
	ldp	x2, x3, [x8, #16]
	adcs	x8, x10, x14
	str	x8, [sp, #88]
	adc	x8, x12, x11
	subs	x26, x9, x18
	sbcs	x10, x15, x0
	sbcs	x11, x16, x2
	sbcs	x27, x17, x3
	csel	x12, x1, xzr, lo
	subs	x26, x26, x12
	sbcs	x10, x10, xzr
	sbcs	x11, x11, xzr
	sbc	x27, x27, xzr
	umulh	x13, x26, x22
	str	x8, [sp, #192]
	stp	x3, x0, [sp, #72]
	mul	x0, x26, x25
	stp	x2, x18, [sp, #56]
	umulh	x21, x26, x25
	mul	x2, x26, x24
	adds	x13, x0, x13
	umulh	x3, x26, x24
	mul	x4, x26, x23
	adcs	x0, x2, x21
	umulh	x7, x10, x22
	mul	x14, x10, x25
	adcs	x2, x4, x3
	umulh	x5, x26, x23
	umulh	x19, x10, x25
	mul	x16, x10, x24
	cinc	x3, x5, hs
	adds	x14, x14, x7
	umulh	x17, x10, x24
	mul	x18, x10, x23
	adcs	x16, x16, x19
	mul	x20, x10, x22
	umulh	x15, x10, x23
	adcs	x17, x18, x17
	umulh	x28, x11, x22
	mul	x30, x11, x25
	cinc	x15, x15, hs
	adds	x13, x13, x20
	umulh	x12, x11, x25
	adcs	x14, x0, x14
	adcs	x16, x2, x16
	mul	x8, x11, x24
	adcs	x17, x3, x17
	cinc	x15, x15, hs
	adds	x3, x30, x28
	umulh	x10, x11, x24
	mul	x9, x11, x23
	adcs	x8, x8, x12
	mul	x6, x11, x22
	umulh	x11, x11, x23
	adcs	x9, x9, x10
	mul	x21, x27, x25
	umulh	x2, x27, x22
	cinc	x10, x11, hs
	adds	x11, x14, x6
	umulh	x4, x27, x25
	adcs	x16, x16, x3
	adcs	x8, x17, x8
	ldp	x6, x5, [sp, #144]
	mul	x18, x27, x24
	adcs	x9, x15, x9
	cinc	x10, x10, hs
	adds	x15, x21, x2
	umulh	x0, x27, x24
	mul	x12, x27, x23
	adcs	x18, x18, x4
	umulh	x14, x27, x23
	mul	x17, x27, x22
	adcs	x12, x12, x0
	mul	x0, x26, x22
	cinc	x14, x14, hs
	adds	x16, x16, x17
	adcs	x8, x8, x15
	adcs	x9, x9, x18
	mul	x15, x8, x1
	adcs	x10, x10, x12
	umulh	x8, x8, x1
	cinc	x14, x14, hs
	mul	x17, x9, x1
	umulh	x9, x9, x1
	mul	x12, x10, x1
	adds	x8, x17, x8
	umulh	x10, x10, x1
	mul	x18, x14, x1
	adcs	x9, x12, x9
	umulh	x14, x14, x1
	adcs	x10, x18, x10
	cinc	x12, x14, hs
	adds	x14, x0, x15
	adcs	x8, x13, x8
	adcs	x9, x11, x9
	cset	w11, hs
	adcs	x13, x16, x10
	add	x10, x16, x10
	cinc	x12, x12, hs
	mul	x15, x12, x1
	umulh	x12, x12, x1
	adds	x14, x14, x15
	adcs	x8, x8, x12
	cset	w12, hs
	adcs	xzr, x9, xzr
	adcs	xzr, x13, xzr
	csel	x13, x1, xzr, hs
	adds	x22, x14, x13
	adcs	x25, x8, xzr
	umulh	x2, x22, x22
	adcs	x24, x9, x12
	umulh	x13, x25, x22
	adc	x23, x10, x11
	mul	x9, x24, x22
	umulh	x12, x24, x22
	mul	x10, x23, x22
	adds	x9, x13, x9
	umulh	x11, x23, x22
	umulh	x15, x24, x25
	adcs	x10, x12, x10
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
	ldp	x2, x0, [sp, #112]
	umulh	x15, x24, x24
	adcs	x9, x9, x12
	mul	x18, x23, x23
	adcs	x10, x10, x14
	umulh	x12, x23, x23
	adcs	x11, x11, x15
	mul	x14, x10, x1
	adcs	x13, x13, x18
	umulh	x10, x10, x1
	adc	x12, x16, x12
	mul	x15, x11, x1
	umulh	x11, x11, x1
	mul	x18, x13, x1
	adds	x10, x15, x10
	umulh	x13, x13, x1
	mul	x16, x12, x1
	adcs	x11, x18, x11
	mul	x4, x22, x22
	umulh	x12, x12, x1
	adcs	x13, x16, x13
	cinc	x12, x12, hs
	adds	x14, x4, x14
	adcs	x8, x8, x10
	ldp	x4, x3, [sp, #160]
	adcs	x10, x17, x11
	ldp	x17, x18, [sp, #128]
	cset	w11, hs
	adcs	x15, x9, x13
	add	x9, x9, x13
	cinc	x12, x12, hs
	mul	x16, x12, x1
	umulh	x12, x12, x1
	adds	x14, x14, x16
	adcs	x8, x8, x12
	cset	w12, hs
	adcs	xzr, x10, xzr
	adcs	xzr, x15, xzr
	csel	x15, x1, xzr, hs
	adds	x14, x14, x15
	adcs	x8, x8, xzr
	adcs	x10, x10, x12
	adc	x9, x9, x11
	subs	x11, x14, x17
	sbcs	x12, x8, x18
	sbcs	x13, x10, x2
	sbcs	x15, x9, x0
	csel	x16, x1, xzr, lo
	subs	x11, x11, x16
	sbcs	x12, x12, xzr
	sbcs	x13, x13, xzr
	sbc	x15, x15, xzr
	subs	x8, x11, x5
	sbcs	x9, x12, x3
	sbcs	x10, x13, x6
	sbcs	x14, x15, x4
	csel	x16, x1, xzr, lo
	subs	x8, x8, x16
	sbcs	x9, x9, xzr
	sbcs	x10, x10, xzr
	sbc	x14, x14, xzr
	adds	x11, x8, x1
	adcs	x12, x9, xzr
	adcs	x13, x10, xzr
	adcs	x15, x14, xzr
	csel	x16, x11, x8, hs
	csel	x12, x12, x9, hs
	csel	x11, x13, x10, hs
	stp	x16, x12, [sp, #144]
	csel	x10, x15, x14, hs
	subs	x26, x17, x16
	sbcs	x8, x18, x12
	sbcs	x28, x2, x11
	sbcs	x27, x0, x10
	csel	x9, x1, xzr, lo
	subs	x26, x26, x9
	sbcs	x8, x8, xzr
	sbcs	x28, x28, xzr
	sbc	x27, x27, xzr
	umulh	x16, x22, x26
	stp	x11, x10, [sp, #160]
	mul	x15, x22, x8
	umulh	x14, x22, x8
	mul	x13, x22, x28
	adds	x15, x15, x16
	umulh	x5, x22, x28
	mul	x7, x22, x27
	adcs	x13, x13, x14
	umulh	x6, x25, x26
	mul	x17, x25, x8
	adcs	x16, x7, x5
	umulh	x19, x22, x27
	umulh	x4, x25, x8
	mul	x18, x25, x28
	cinc	x5, x19, hs
	adds	x17, x17, x6
	umulh	x3, x25, x28
	mul	x0, x25, x27
	adcs	x18, x18, x4
	mul	x21, x25, x26
	umulh	x2, x25, x27
	adcs	x0, x0, x3
	umulh	x25, x24, x26
	mul	x30, x24, x8
	cinc	x2, x2, hs
	adds	x15, x15, x21
	umulh	x12, x24, x8
	adcs	x13, x13, x17
	adcs	x16, x16, x18
	mul	x9, x24, x28
	adcs	x18, x5, x0
	cinc	x2, x2, hs
	adds	x3, x30, x25
	umulh	x11, x24, x28
	mul	x10, x24, x27
	adcs	x9, x9, x12
	mul	x20, x24, x26
	umulh	x24, x24, x27
	adcs	x10, x10, x11
	mul	x14, x23, x8
	umulh	x0, x23, x26
	cinc	x11, x24, hs
	adds	x13, x13, x20
	umulh	x8, x23, x8
	adcs	x16, x16, x3
	adcs	x9, x18, x9
	mul	x6, x23, x28
	adcs	x10, x2, x10
	cinc	x11, x11, hs
	adds	x14, x14, x0
	umulh	x17, x23, x28
	mul	x12, x23, x27
	adcs	x8, x6, x8
	umulh	x4, x23, x27
	mul	x18, x23, x26
	adcs	x12, x12, x17
	mul	x17, x22, x26
	cinc	x0, x4, hs
	adds	x16, x16, x18
	adcs	x9, x9, x14
	adcs	x8, x10, x8
	mul	x14, x9, x1
	adcs	x11, x11, x12
	umulh	x9, x9, x1
	cinc	x18, x0, hs
	mul	x10, x8, x1
	umulh	x8, x8, x1
	mul	x12, x11, x1
	adds	x9, x10, x9
	umulh	x11, x11, x1
	mul	x0, x18, x1
	adcs	x8, x12, x8
	umulh	x18, x18, x1
	adcs	x10, x0, x11
	ldp	x2, x0, [sp, #56]
	cinc	x11, x18, hs
	adds	x12, x17, x14
	adcs	x9, x15, x9
	ldp	x18, x17, [sp, #72]
	adcs	x8, x13, x8
	cset	w13, hs
	adcs	x14, x16, x10
	add	x10, x16, x10
	cinc	x11, x11, hs
	mul	x15, x11, x1
	umulh	x11, x11, x1
	adds	x12, x12, x15
	adcs	x9, x9, x11
	cset	w11, hs
	adcs	xzr, x8, xzr
	adcs	xzr, x14, xzr
	csel	x14, x1, xzr, hs
	adds	x12, x12, x14
	adcs	x9, x9, xzr
	adcs	x8, x8, x11
	adc	x10, x10, x13
	subs	x11, x12, x0
	sbcs	x13, x9, x17
	sbcs	x14, x8, x2
	sbcs	x15, x10, x18
	csel	x16, x1, xzr, lo
	subs	x11, x11, x16
	sbcs	x13, x13, xzr
	sbcs	x14, x14, xzr
	sbc	x15, x15, xzr
	ldp	x9, x10, [sp, #144]
	ldr	x16, [sp, #216]
	adds	x8, x11, x1
	ldp	x12, x17, [sp, #160]
	stp	x9, x10, [x16, #-32]
	adcs	x9, x13, xzr
	adcs	x10, x14, xzr
	stp	x12, x17, [x16, #-16]
	adcs	x12, x15, xzr
	csel	x8, x8, x11, hs
	ldp	x11, x2, [sp, #176]
	csel	x9, x9, x13, hs
	stp	x8, x9, [x16]
	csel	x8, x10, x14, hs
	ldr	x10, [sp, #200]
	csel	x9, x12, x15, hs
	ldp	x14, x15, [sp, #88]
	stp	x8, x9, [x16, #16]
	ldr	x9, [sp, #104]
	sub	x10, x10, #1
	sub	x16, x16, #64
	sub	x11, x11, #32
	str	x16, [sp, #216]
	cbnz	x10, .LBB1_2
.LBB1_3:
	ldr	x26, [sp, #8]
	ldr	x30, [x2]
	mov	x8, #977
	movk	x8, #1, lsl #32
	ldp	x12, x13, [x26, #48]
	ldp	x2, x3, [x30, #32]
	ldp	x4, x5, [x30, #48]
	ldp	x10, x11, [x26, #32]
	str	x13, [sp, #200]
	subs	x18, x2, x10
	sbcs	x6, x3, x11
	sbcs	x7, x4, x12
	sbcs	x0, x5, x13
	csel	x19, x8, xzr, lo
	subs	x18, x18, x19
	sbcs	x6, x6, xzr
	sbcs	x7, x7, xzr
	sbc	x0, x0, xzr
	umulh	x2, x18, x9
	ldr	x13, [sp, #192]
	mov	x8, #977
	mul	x3, x18, x15
	str	x11, [sp, #216]
	movk	x8, #1, lsl #32
	stp	x12, x10, [sp, #176]
	umulh	x4, x18, x15
	mul	x5, x18, x14
	adds	x2, x3, x2
	umulh	x22, x18, x14
	mul	x19, x18, x13
	adcs	x3, x5, x4
	umulh	x21, x18, x13
	mul	x11, x6, x9
	adcs	x5, x19, x22
	umulh	x23, x6, x9
	cinc	x19, x21, hs
	mul	x24, x6, x15
	mul	x10, x7, x9
	umulh	x25, x6, x15
	adds	x21, x24, x23
	mul	x20, x6, x14
	stp	x10, x11, [sp, #160]
	umulh	x27, x6, x14
	mul	x28, x6, x13
	adcs	x23, x20, x25
	ldr	x20, [sp, #168]
	umulh	x6, x6, x13
	umulh	x1, x7, x9
	adcs	x24, x28, x27
	mul	x12, x7, x15
	cinc	x6, x6, hs
	adds	x2, x2, x20
	umulh	x17, x7, x15
	adcs	x3, x3, x21
	adcs	x5, x5, x23
	mul	x10, x7, x14
	adcs	x19, x19, x24
	cinc	x6, x6, hs
	adds	x12, x12, x1
	umulh	x16, x7, x14
	mul	x11, x7, x13
	adcs	x10, x10, x17
	ldr	x17, [sp, #160]
	umulh	x7, x7, x13
	umulh	x22, x0, x9
	adcs	x11, x11, x16
	mul	x25, x0, x15
	cinc	x16, x7, hs
	adds	x17, x3, x17
	umulh	x15, x0, x15
	adcs	x12, x5, x12
	adcs	x10, x19, x10
	mul	x21, x0, x14
	adcs	x11, x6, x11
	cinc	x16, x16, hs
	umulh	x14, x0, x14
	mul	x1, x0, x13
	mul	x4, x0, x9
	umulh	x13, x0, x13
	adds	x0, x25, x22
	adcs	x15, x21, x15
	adcs	x14, x1, x14
	mul	x9, x18, x9
	cinc	x13, x13, hs
	adds	x12, x12, x4
	adcs	x10, x10, x0
	adcs	x11, x11, x15
	mul	x18, x10, x8
	adcs	x14, x16, x14
	umulh	x10, x10, x8
	cinc	x13, x13, hs
	mul	x15, x11, x8
	umulh	x11, x11, x8
	mul	x16, x14, x8
	adds	x10, x15, x10
	umulh	x14, x14, x8
	mul	x0, x13, x8
	adcs	x11, x16, x11
	umulh	x13, x13, x8
	adcs	x16, x0, x14
	cinc	x13, x13, hs
	adds	x9, x9, x18
	adcs	x10, x2, x10
	adcs	x11, x17, x11
	cset	w17, hs
	adcs	x14, x12, x16
	add	x12, x12, x16
	cinc	x13, x13, hs
	mul	x15, x13, x8
	umulh	x13, x13, x8
	adds	x9, x9, x15
	adcs	x10, x10, x13
	cset	w13, hs
	adcs	xzr, x11, xzr
	adcs	xzr, x14, xzr
	csel	x14, x8, xzr, hs
	adds	x9, x9, x14
	adcs	x15, x10, xzr
	umulh	x7, x9, x9
	adcs	x14, x11, x13
	umulh	x18, x15, x9
	adc	x13, x12, x17
	mul	x11, x14, x9
	umulh	x16, x14, x9
	mul	x12, x13, x9
	adds	x11, x18, x11
	umulh	x17, x13, x9
	umulh	x1, x14, x15
	adcs	x12, x16, x12
	mul	x2, x13, x15
	cinc	x17, x17, hs
	mul	x0, x14, x15
	umulh	x3, x13, x15
	adds	x18, x1, x2
	mul	x10, x15, x9
	mul	x4, x13, x14
	cinc	x1, x3, hs
	adds	x12, x12, x0
	umulh	x5, x13, x14
	adcs	x17, x17, x18
	extr	x3, x11, x10, #63
	extr	x11, x12, x11, #63
	extr	x12, x17, x12, #63
	mul	x19, x15, x15
	adcs	x18, x1, x4
	umulh	x16, x15, x15
	extr	x17, x18, x17, #63
	cinc	x2, x5, hs
	adds	x10, x7, x10, lsl #1
	mul	x0, x14, x14
	extr	x18, x2, x18, #63
	lsr	x2, x2, #63
	adcs	x3, x3, x19
	ldp	x7, x19, [x26, #16]
	umulh	x1, x14, x14
	adcs	x11, x11, x16
	mul	x4, x13, x13
	adcs	x12, x12, x0
	umulh	x5, x13, x13
	adcs	x17, x17, x1
	mul	x16, x12, x8
	adcs	x18, x18, x4
	umulh	x12, x12, x8
	adc	x2, x2, x5
	mul	x0, x17, x8
	umulh	x17, x17, x8
	mul	x1, x18, x8
	adds	x12, x0, x12
	umulh	x18, x18, x8
	mul	x4, x2, x8
	adcs	x17, x1, x17
	mul	x6, x9, x9
	umulh	x2, x2, x8
	adcs	x18, x4, x18
	cinc	x0, x2, hs
	adds	x16, x6, x16
	adcs	x10, x10, x12
	adcs	x12, x3, x17
	cset	w17, hs
	adcs	x1, x11, x18
	add	x11, x11, x18
	cinc	x0, x0, hs
	mul	x2, x0, x8
	umulh	x0, x0, x8
	adds	x16, x16, x2
	adcs	x10, x10, x0
	cset	w0, hs
	adcs	xzr, x12, xzr
	adcs	xzr, x1, xzr
	csel	x1, x8, xzr, hs
	adds	x16, x16, x1
	ldp	x1, x6, [x26]
	adcs	x10, x10, xzr
	adcs	x12, x12, x0
	ldp	x0, x2, [x30, #16]
	adc	x11, x11, x17
	ldp	x17, x18, [x30]
	subs	x4, x16, x1
	sbcs	x5, x10, x6
	sbcs	x20, x12, x7
	sbcs	x21, x11, x19
	csel	x3, x8, xzr, lo
	subs	x4, x4, x3
	sbcs	x5, x5, xzr
	sbcs	x20, x20, xzr
	sbc	x21, x21, xzr
	subs	x10, x4, x17
	sbcs	x11, x5, x18
	sbcs	x12, x20, x0
	sbcs	x16, x21, x2
	csel	x3, x8, xzr, lo
	subs	x10, x10, x3
	sbcs	x11, x11, xzr
	sbcs	x12, x12, xzr
	sbc	x16, x16, xzr
	adds	x17, x10, x8
	adcs	x18, x11, xzr
	adcs	x2, x12, xzr
	adcs	x3, x16, xzr
	csel	x17, x17, x10, hs
	csel	x18, x18, x11, hs
	csel	x12, x2, x12, hs
	stp	x17, x18, [sp, #152]
	csel	x16, x3, x16, hs
	subs	x3, x1, x17
	sbcs	x10, x6, x18
	sbcs	x5, x7, x12
	sbcs	x4, x19, x16
	csel	x11, x8, xzr, lo
	subs	x3, x3, x11
	sbcs	x10, x10, xzr
	sbcs	x5, x5, xzr
	sbc	x4, x4, xzr
	umulh	x0, x9, x3
	str	x16, [sp, #192]
	str	x12, [sp, #168]
	mul	x18, x9, x10
	umulh	x16, x9, x10
	mul	x1, x9, x5
	adds	x18, x18, x0
	umulh	x6, x9, x5
	mul	x7, x9, x4
	adcs	x16, x1, x16
	umulh	x19, x9, x4
	umulh	x21, x15, x3
	adcs	x1, x7, x6
	mul	x22, x15, x10
	cinc	x6, x19, hs
	umulh	x23, x15, x10
	mul	x24, x15, x5
	adds	x7, x22, x21
	mul	x11, x15, x3
	umulh	x25, x15, x5
	adcs	x19, x24, x23
	ldp	x24, x23, [sp, #256]
	mul	x2, x15, x4
	str	x11, [sp, #144]
	umulh	x20, x15, x4
	umulh	x28, x14, x3
	adcs	x21, x2, x25
	mul	x30, x14, x10
	cinc	x2, x20, hs
	umulh	x15, x14, x10
	adds	x18, x18, x11
	mul	x11, x14, x5
	adcs	x16, x16, x7
	adcs	x1, x1, x19
	umulh	x17, x14, x5
	adcs	x6, x6, x21
	cinc	x2, x2, hs
	adds	x20, x30, x28
	ldr	x30, [sp, #208]
	mul	x12, x14, x4
	adcs	x11, x11, x15
	mul	x27, x14, x3
	umulh	x14, x14, x4
	adcs	x12, x12, x17
	mul	x0, x13, x10
	umulh	x7, x13, x3
	cinc	x14, x14, hs
	adds	x16, x16, x27
	umulh	x10, x13, x10
	adcs	x1, x1, x20
	adcs	x11, x6, x11
	ldp	x28, x27, [sp, #224]
	mul	x19, x13, x5
	adcs	x12, x2, x12
	cinc	x14, x14, hs
	adds	x0, x0, x7
	umulh	x15, x13, x5
	mul	x17, x13, x4
	adcs	x10, x19, x10
	ldp	x20, x19, [sp, #288]
	mul	x22, x13, x3
	umulh	x13, x13, x4
	adcs	x15, x17, x15
	mul	x9, x9, x3
	cinc	x13, x13, hs
	adds	x17, x1, x22
	adcs	x11, x11, x0
	ldp	x22, x21, [sp, #272]
	adcs	x10, x12, x10
	mul	x0, x11, x8
	adcs	x14, x14, x15
	umulh	x11, x11, x8
	cinc	x13, x13, hs
	mul	x12, x10, x8
	umulh	x10, x10, x8
	mul	x15, x14, x8
	adds	x11, x12, x11
	umulh	x14, x14, x8
	mul	x1, x13, x8
	adcs	x10, x15, x10
	umulh	x13, x13, x8
	adcs	x12, x1, x14
	ldp	x2, x1, [sp, #176]
	cinc	x13, x13, hs
	adds	x9, x9, x0
	ldr	x0, [sp, #200]
	adcs	x11, x18, x11
	ldr	x18, [sp, #216]
	adcs	x10, x16, x10
	cset	w14, hs
	adcs	x15, x17, x12
	add	x12, x17, x12
	cinc	x13, x13, hs
	mul	x16, x13, x8
	umulh	x13, x13, x8
	adds	x9, x9, x16
	adcs	x11, x11, x13
	cset	w13, hs
	adcs	xzr, x10, xzr
	adcs	xzr, x15, xzr
	csel	x15, x8, xzr, hs
	adds	x9, x9, x15
	adcs	x11, x11, xzr
	adcs	x10, x10, x13
	adc	x12, x12, x14
	subs	x13, x9, x1
	sbcs	x14, x11, x18
	sbcs	x15, x10, x2
	sbcs	x16, x12, x0
	csel	x17, x8, xzr, lo
	subs	x13, x13, x17
	sbcs	x14, x14, xzr
	sbcs	x15, x15, xzr
	sbc	x16, x16, xzr
	ldp	x10, x11, [sp, #152]
	adds	x8, x13, x8
	adcs	x9, x14, xzr
	ldr	x12, [sp, #168]
	stp	x10, x11, [x26]
	adcs	x10, x15, xzr
	adcs	x11, x16, xzr
	str	x12, [x26, #16]
	ldr	x12, [sp, #192]
	csel	x8, x8, x13, hs
	csel	x9, x9, x14, hs
	stp	x12, x8, [x26, #24]
	csel	x8, x10, x15, hs
	stp	x9, x8, [x26, #40]
	csel	x9, x11, x16, hs
	str	x9, [x26, #56]
	ldp	x26, x25, [sp, #240]
	add	sp, sp, #304
	ret

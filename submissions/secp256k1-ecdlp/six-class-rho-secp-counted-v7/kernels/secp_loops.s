	.text
	.p2align 2
	.globl pass1
	.type pass1, %function
pass1:
	stp x19, x20, [sp, #-80]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	mov x4, #977
	movk x4, #1, lsl #32
	ldr x5, [x2], #8
	ldp x6, x7, [x5, #0]
	ldp x8, x9, [x5, #16]
	ldp x12, x13, [x1, #16]
	ldp x10, x11, [x1], #64
	subs x6, x6, x10
	sbcs x7, x7, x11
	sbcs x8, x8, x12
	sbcs x9, x9, x13
	csel x14, x4, xzr, lo
	subs x6, x6, x14
	sbcs x7, x7, xzr
	sbcs x8, x8, xzr
	sbc x9, x9, xzr
	stp x8, x9, [x3, #16]
	stp x6, x7, [x3], #32
	sub x0, x0, #1
1:
	ldr x5, [x2], #8
	ldp x6, x7, [x5, #0]
	ldp x8, x9, [x5, #16]
	ldp x12, x13, [x1, #16]
	ldp x10, x11, [x1], #64
	ldp x14, x15, [x3, #-32]
	ldp x16, x17, [x3, #-16]
	subs x6, x6, x10
	sbcs x7, x7, x11
	sbcs x8, x8, x12
	sbcs x9, x9, x13
	csel x19, x4, xzr, lo
	subs x6, x6, x19
	sbcs x7, x7, xzr
	sbcs x8, x8, xzr
	sbc x9, x9, xzr
	mul x13, x14, x6
	umulh x12, x14, x6
	mul x5, x14, x7
	umulh x11, x14, x7
	adds x12, x12, x5
	mul x5, x14, x8
	umulh x10, x14, x8
	adcs x11, x11, x5
	mul x5, x14, x9
	umulh x19, x14, x9
	adcs x10, x10, x5
	adc x19, x19, xzr
	mul x14, x15, x6
	umulh x20, x15, x6
	mul x5, x15, x7
	umulh x21, x15, x7
	adds x20, x20, x5
	mul x5, x15, x8
	umulh x22, x15, x8
	adcs x21, x21, x5
	mul x5, x15, x9
	umulh x23, x15, x9
	adcs x22, x22, x5
	adc x23, x23, xzr
	adds x12, x12, x14
	adcs x11, x11, x20
	adcs x10, x10, x21
	adcs x19, x19, x22
	adc x23, x23, xzr
	mul x15, x16, x6
	umulh x22, x16, x6
	mul x5, x16, x7
	umulh x21, x16, x7
	adds x22, x22, x5
	mul x5, x16, x8
	umulh x20, x16, x8
	adcs x21, x21, x5
	mul x5, x16, x9
	umulh x14, x16, x9
	adcs x20, x20, x5
	adc x14, x14, xzr
	adds x11, x11, x15
	adcs x10, x10, x22
	adcs x19, x19, x21
	adcs x23, x23, x20
	adc x14, x14, xzr
	mul x16, x17, x6
	umulh x20, x17, x6
	mul x5, x17, x7
	umulh x21, x17, x7
	adds x20, x20, x5
	mul x5, x17, x8
	umulh x22, x17, x8
	adcs x21, x21, x5
	mul x5, x17, x9
	umulh x15, x17, x9
	adcs x22, x22, x5
	adc x15, x15, xzr
	adds x10, x10, x16
	adcs x19, x19, x20
	adcs x23, x23, x21
	adcs x14, x14, x22
	adc x15, x15, xzr
	mul x9, x19, x4
	umulh x5, x19, x4
	mul x8, x23, x4
	umulh x17, x23, x4
	mul x7, x14, x4
	umulh x22, x14, x4
	mul x6, x15, x4
	umulh x21, x15, x4
	adds x8, x8, x5
	adcs x7, x7, x17
	adcs x6, x6, x22
	adc x21, x21, xzr
	adds x13, x13, x9
	adcs x12, x12, x8
	adcs x11, x11, x7
	adcs x10, x10, x6
	adc x21, x21, xzr
	mul x6, x21, x4
	umulh x7, x21, x4
	adds x13, x13, x6
	adcs x12, x12, x7
	adcs x11, x11, xzr
	adcs x10, x10, xzr
	csel x6, x4, xzr, cs
	adds x13, x13, x6
	adcs x12, x12, xzr
	adcs x11, x11, xzr
	adc x10, x10, xzr
	stp x11, x10, [x3, #16]
	stp x13, x12, [x3], #32
	subs x0, x0, #1
	b.ne 1b
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #80
	ret
	.text
	.p2align 2
	.globl pass2
	.type pass2, %function
pass2:
	stp x19, x20, [sp, #-112]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	ldp x5, x6, [x4]
	ldp x7, x8, [x4, #16]
	stp x5, x6, [sp, #80]
	stp x7, x8, [sp, #96]
	mov x4, #977
	movk x4, #1, lsl #32
	sub x5, x0, #1
	add x1, x1, x5, lsl #6
	add x1, x1, #32
	add x2, x2, x5, lsl #3
	sub x6, x5, #1
	add x3, x3, x6, lsl #5
	mov x0, x5
1:
	ldr x5, [x2], #-8
	ldp x6, x7, [sp, #80]
	ldp x8, x9, [sp, #96]
	ldp x12, x13, [x3, #16]
	ldp x10, x11, [x3], #-32
	mul x15, x6, x10
	umulh x16, x6, x10
	mul x14, x6, x11
	umulh x17, x6, x11
	adds x16, x16, x14
	mul x14, x6, x12
	umulh x19, x6, x12
	adcs x17, x17, x14
	mul x14, x6, x13
	umulh x20, x6, x13
	adcs x19, x19, x14
	adc x20, x20, xzr
	mul x21, x7, x10
	umulh x22, x7, x10
	mul x14, x7, x11
	umulh x23, x7, x11
	adds x22, x22, x14
	mul x14, x7, x12
	umulh x24, x7, x12
	adcs x23, x23, x14
	mul x14, x7, x13
	umulh x25, x7, x13
	adcs x24, x24, x14
	adc x25, x25, xzr
	adds x16, x16, x21
	adcs x17, x17, x22
	adcs x19, x19, x23
	adcs x20, x20, x24
	adc x25, x25, xzr
	mul x24, x8, x10
	umulh x23, x8, x10
	mul x14, x8, x11
	umulh x22, x8, x11
	adds x23, x23, x14
	mul x14, x8, x12
	umulh x21, x8, x12
	adcs x22, x22, x14
	mul x14, x8, x13
	umulh x26, x8, x13
	adcs x21, x21, x14
	adc x26, x26, xzr
	adds x17, x17, x24
	adcs x19, x19, x23
	adcs x20, x20, x22
	adcs x25, x25, x21
	adc x26, x26, xzr
	mul x21, x9, x10
	umulh x22, x9, x10
	mul x14, x9, x11
	umulh x23, x9, x11
	adds x22, x22, x14
	mul x14, x9, x12
	umulh x24, x9, x12
	adcs x23, x23, x14
	mul x14, x9, x13
	umulh x27, x9, x13
	adcs x24, x24, x14
	adc x27, x27, xzr
	adds x19, x19, x21
	adcs x20, x20, x22
	adcs x25, x25, x23
	adcs x26, x26, x24
	adc x27, x27, xzr
	mul x13, x20, x4
	umulh x14, x20, x4
	mul x12, x25, x4
	umulh x24, x25, x4
	mul x11, x26, x4
	umulh x23, x26, x4
	mul x10, x27, x4
	umulh x22, x27, x4
	adds x12, x12, x14
	adcs x11, x11, x24
	adcs x10, x10, x23
	adc x22, x22, xzr
	adds x15, x15, x13
	adcs x16, x16, x12
	adcs x17, x17, x11
	adcs x19, x19, x10
	adc x22, x22, xzr
	mul x10, x22, x4
	umulh x11, x22, x4
	adds x15, x15, x10
	adcs x16, x16, x11
	adcs x17, x17, xzr
	adcs x19, x19, xzr
	csel x10, x4, xzr, cs
	adds x15, x15, x10
	adcs x16, x16, xzr
	adcs x17, x17, xzr
	adc x19, x19, xzr
	ldp x11, x10, [x5, #0]
	ldp x22, x12, [x5, #16]
	ldp x13, x23, [x1, #-32]
	ldp x24, x14, [x1, #-16]
	subs x11, x11, x13
	sbcs x10, x10, x23
	sbcs x22, x22, x24
	sbcs x12, x12, x14
	csel x27, x4, xzr, lo
	subs x11, x11, x27
	sbcs x10, x10, xzr
	sbcs x22, x22, xzr
	sbc x12, x12, xzr
	mul x24, x6, x11
	umulh x23, x6, x11
	mul x14, x6, x10
	umulh x13, x6, x10
	adds x23, x23, x14
	mul x14, x6, x22
	umulh x27, x6, x22
	adcs x13, x13, x14
	mul x14, x6, x12
	umulh x26, x6, x12
	adcs x27, x27, x14
	adc x26, x26, xzr
	mul x6, x7, x11
	umulh x25, x7, x11
	mul x14, x7, x10
	umulh x20, x7, x10
	adds x25, x25, x14
	mul x14, x7, x22
	umulh x21, x7, x22
	adcs x20, x20, x14
	mul x14, x7, x12
	umulh x28, x7, x12
	adcs x21, x21, x14
	adc x28, x28, xzr
	adds x23, x23, x6
	adcs x13, x13, x25
	adcs x27, x27, x20
	adcs x26, x26, x21
	adc x28, x28, xzr
	mul x7, x8, x11
	umulh x21, x8, x11
	mul x14, x8, x10
	umulh x20, x8, x10
	adds x21, x21, x14
	mul x14, x8, x22
	umulh x25, x8, x22
	adcs x20, x20, x14
	mul x14, x8, x12
	umulh x6, x8, x12
	adcs x25, x25, x14
	adc x6, x6, xzr
	adds x13, x13, x7
	adcs x27, x27, x21
	adcs x26, x26, x20
	adcs x28, x28, x25
	adc x6, x6, xzr
	mul x8, x9, x11
	umulh x25, x9, x11
	mul x14, x9, x10
	umulh x20, x9, x10
	adds x25, x25, x14
	mul x14, x9, x22
	umulh x21, x9, x22
	adcs x20, x20, x14
	mul x14, x9, x12
	umulh x7, x9, x12
	adcs x21, x21, x14
	adc x7, x7, xzr
	adds x27, x27, x8
	adcs x26, x26, x25
	adcs x28, x28, x20
	adcs x6, x6, x21
	adc x7, x7, xzr
	mul x12, x26, x4
	umulh x14, x26, x4
	mul x22, x28, x4
	umulh x9, x28, x4
	mul x10, x6, x4
	umulh x21, x6, x4
	mul x11, x7, x4
	umulh x20, x7, x4
	adds x22, x22, x14
	adcs x10, x10, x9
	adcs x11, x11, x21
	adc x20, x20, xzr
	adds x24, x24, x12
	adcs x23, x23, x22
	adcs x13, x13, x10
	adcs x27, x27, x11
	adc x20, x20, xzr
	mul x11, x20, x4
	umulh x10, x20, x4
	adds x24, x24, x11
	adcs x23, x23, x10
	adcs x13, x13, xzr
	adcs x27, x27, xzr
	csel x11, x4, xzr, cs
	adds x24, x24, x11
	adcs x23, x23, xzr
	adcs x13, x13, xzr
	adc x27, x27, xzr
	stp x24, x23, [sp, #80]
	stp x13, x27, [sp, #96]
	ldp x27, x13, [x5, #32]
	ldp x23, x24, [x5, #48]
	ldp x10, x11, [x1, #0]
	ldp x20, x22, [x1, #16]
	subs x27, x27, x10
	sbcs x13, x13, x11
	sbcs x23, x23, x20
	sbcs x24, x24, x22
	csel x12, x4, xzr, lo
	subs x27, x27, x12
	sbcs x13, x13, xzr
	sbcs x23, x23, xzr
	sbc x24, x24, xzr
	mul x20, x27, x15
	umulh x11, x27, x15
	mul x22, x27, x16
	umulh x10, x27, x16
	adds x11, x11, x22
	mul x22, x27, x17
	umulh x12, x27, x17
	adcs x10, x10, x22
	mul x22, x27, x19
	umulh x21, x27, x19
	adcs x12, x12, x22
	adc x21, x21, xzr
	mul x27, x13, x15
	umulh x9, x13, x15
	mul x22, x13, x16
	umulh x14, x13, x16
	adds x9, x9, x22
	mul x22, x13, x17
	umulh x7, x13, x17
	adcs x14, x14, x22
	mul x22, x13, x19
	umulh x6, x13, x19
	adcs x7, x7, x22
	adc x6, x6, xzr
	adds x11, x11, x27
	adcs x10, x10, x9
	adcs x12, x12, x14
	adcs x21, x21, x7
	adc x6, x6, xzr
	mul x13, x23, x15
	umulh x7, x23, x15
	mul x22, x23, x16
	umulh x14, x23, x16
	adds x7, x7, x22
	mul x22, x23, x17
	umulh x9, x23, x17
	adcs x14, x14, x22
	mul x22, x23, x19
	umulh x27, x23, x19
	adcs x9, x9, x22
	adc x27, x27, xzr
	adds x10, x10, x13
	adcs x12, x12, x7
	adcs x21, x21, x14
	adcs x6, x6, x9
	adc x27, x27, xzr
	mul x23, x24, x15
	umulh x9, x24, x15
	mul x22, x24, x16
	umulh x14, x24, x16
	adds x9, x9, x22
	mul x22, x24, x17
	umulh x7, x24, x17
	adcs x14, x14, x22
	mul x22, x24, x19
	umulh x13, x24, x19
	adcs x7, x7, x22
	adc x13, x13, xzr
	adds x12, x12, x23
	adcs x21, x21, x9
	adcs x6, x6, x14
	adcs x27, x27, x7
	adc x13, x13, xzr
	mul x19, x21, x4
	umulh x22, x21, x4
	mul x17, x6, x4
	umulh x24, x6, x4
	mul x16, x27, x4
	umulh x7, x27, x4
	mul x15, x13, x4
	umulh x14, x13, x4
	adds x17, x17, x22
	adcs x16, x16, x24
	adcs x15, x15, x7
	adc x14, x14, xzr
	adds x20, x20, x19
	adcs x11, x11, x17
	adcs x10, x10, x16
	adcs x12, x12, x15
	adc x14, x14, xzr
	mul x15, x14, x4
	umulh x16, x14, x4
	adds x20, x20, x15
	adcs x11, x11, x16
	adcs x10, x10, xzr
	adcs x12, x12, xzr
	csel x15, x4, xzr, cs
	adds x20, x20, x15
	adcs x11, x11, xzr
	adcs x10, x10, xzr
	adc x12, x12, xzr
	mul x16, x20, x11
	umulh x15, x20, x11
	mul x14, x20, x10
	umulh x17, x20, x10
	mul x19, x20, x12
	umulh x7, x20, x12
	mul x24, x11, x10
	umulh x22, x11, x10
	mul x13, x11, x12
	umulh x27, x11, x12
	mul x6, x10, x12
	umulh x21, x10, x12
	adds x15, x15, x14
	adcs x17, x17, x19
	adc x7, x7, xzr
	adds x22, x22, x13
	adc x27, x27, xzr
	adds x17, x17, x24
	adcs x7, x7, x22
	adcs x27, x27, x6
	adc x21, x21, xzr
	lsr x6, x21, #63
	extr x21, x21, x27, #63
	extr x27, x27, x7, #63
	extr x7, x7, x17, #63
	extr x17, x17, x15, #63
	extr x15, x15, x16, #63
	lsl x16, x16, #1
	mul x22, x20, x20
	umulh x24, x20, x20
	mul x9, x11, x11
	umulh x13, x11, x11
	mul x23, x10, x10
	umulh x19, x10, x10
	mul x28, x12, x12
	umulh x14, x12, x12
	adds x16, x16, x24
	adcs x15, x15, x9
	adcs x17, x17, x13
	adcs x7, x7, x23
	adcs x27, x27, x19
	adcs x21, x21, x28
	adc x6, x6, x14
	mul x28, x7, x4
	umulh x19, x7, x4
	mul x23, x27, x4
	umulh x13, x27, x4
	mul x9, x21, x4
	umulh x24, x21, x4
	mul x14, x6, x4
	umulh x26, x6, x4
	adds x23, x23, x19
	adcs x9, x9, x13
	adcs x14, x14, x24
	adc x26, x26, xzr
	adds x22, x22, x28
	adcs x16, x16, x23
	adcs x15, x15, x9
	adcs x17, x17, x14
	adc x26, x26, xzr
	mul x14, x26, x4
	umulh x9, x26, x4
	adds x22, x22, x14
	adcs x16, x16, x9
	adcs x15, x15, xzr
	adcs x17, x17, xzr
	csel x14, x4, xzr, cs
	adds x22, x22, x14
	adcs x16, x16, xzr
	adcs x15, x15, xzr
	adc x17, x17, xzr
	ldp x9, x14, [x1, #-32]
	ldp x26, x23, [x1, #-16]
	subs x22, x22, x9
	sbcs x16, x16, x14
	sbcs x15, x15, x26
	sbcs x17, x17, x23
	csel x28, x4, xzr, lo
	subs x22, x22, x28
	sbcs x16, x16, xzr
	sbcs x15, x15, xzr
	sbc x17, x17, xzr
	ldp x23, x26, [x5, #0]
	ldp x14, x9, [x5, #16]
	subs x22, x22, x23
	sbcs x16, x16, x26
	sbcs x15, x15, x14
	sbcs x17, x17, x9
	csel x28, x4, xzr, lo
	subs x22, x22, x28
	sbcs x16, x16, xzr
	sbcs x15, x15, xzr
	sbc x17, x17, xzr
	adds x9, x22, x4
	adcs x14, x16, xzr
	adcs x26, x15, xzr
	adcs x23, x17, xzr
	csel x22, x9, x22, cs
	csel x16, x14, x16, cs
	csel x15, x26, x15, cs
	csel x17, x23, x17, cs
	ldp x23, x26, [x1, #-32]
	ldp x14, x9, [x1, #-16]
	subs x23, x23, x22
	sbcs x26, x26, x16
	sbcs x14, x14, x15
	sbcs x9, x9, x17
	csel x28, x4, xzr, lo
	subs x23, x23, x28
	sbcs x26, x26, xzr
	sbcs x14, x14, xzr
	sbc x9, x9, xzr
	stp x22, x16, [x1, #-32]
	stp x15, x17, [x1, #-16]
	mul x15, x20, x23
	umulh x16, x20, x23
	mul x17, x20, x26
	umulh x22, x20, x26
	adds x16, x16, x17
	mul x17, x20, x14
	umulh x28, x20, x14
	adcs x22, x22, x17
	mul x17, x20, x9
	umulh x24, x20, x9
	adcs x28, x28, x17
	adc x24, x24, xzr
	mul x20, x11, x23
	umulh x13, x11, x23
	mul x17, x11, x26
	umulh x19, x11, x26
	adds x13, x13, x17
	mul x17, x11, x14
	umulh x6, x11, x14
	adcs x19, x19, x17
	mul x17, x11, x9
	umulh x21, x11, x9
	adcs x6, x6, x17
	adc x21, x21, xzr
	adds x16, x16, x20
	adcs x22, x22, x13
	adcs x28, x28, x19
	adcs x24, x24, x6
	adc x21, x21, xzr
	mul x11, x10, x23
	umulh x6, x10, x23
	mul x17, x10, x26
	umulh x19, x10, x26
	adds x6, x6, x17
	mul x17, x10, x14
	umulh x13, x10, x14
	adcs x19, x19, x17
	mul x17, x10, x9
	umulh x20, x10, x9
	adcs x13, x13, x17
	adc x20, x20, xzr
	adds x22, x22, x11
	adcs x28, x28, x6
	adcs x24, x24, x19
	adcs x21, x21, x13
	adc x20, x20, xzr
	mul x10, x12, x23
	umulh x13, x12, x23
	mul x17, x12, x26
	umulh x19, x12, x26
	adds x13, x13, x17
	mul x17, x12, x14
	umulh x6, x12, x14
	adcs x19, x19, x17
	mul x17, x12, x9
	umulh x11, x12, x9
	adcs x6, x6, x17
	adc x11, x11, xzr
	adds x28, x28, x10
	adcs x24, x24, x13
	adcs x21, x21, x19
	adcs x20, x20, x6
	adc x11, x11, xzr
	mul x9, x24, x4
	umulh x17, x24, x4
	mul x14, x21, x4
	umulh x12, x21, x4
	mul x26, x20, x4
	umulh x6, x20, x4
	mul x23, x11, x4
	umulh x19, x11, x4
	adds x14, x14, x17
	adcs x26, x26, x12
	adcs x23, x23, x6
	adc x19, x19, xzr
	adds x15, x15, x9
	adcs x16, x16, x14
	adcs x22, x22, x26
	adcs x28, x28, x23
	adc x19, x19, xzr
	mul x23, x19, x4
	umulh x26, x19, x4
	adds x15, x15, x23
	adcs x16, x16, x26
	adcs x22, x22, xzr
	adcs x28, x28, xzr
	csel x23, x4, xzr, cs
	adds x15, x15, x23
	adcs x16, x16, xzr
	adcs x22, x22, xzr
	adc x28, x28, xzr
	ldp x26, x23, [x1, #0]
	ldp x19, x14, [x1, #16]
	subs x15, x15, x26
	sbcs x16, x16, x23
	sbcs x22, x22, x19
	sbcs x28, x28, x14
	csel x9, x4, xzr, lo
	subs x15, x15, x9
	sbcs x16, x16, xzr
	sbcs x22, x22, xzr
	sbc x28, x28, xzr
	adds x14, x15, x4
	adcs x19, x16, xzr
	adcs x23, x22, xzr
	adcs x26, x28, xzr
	csel x15, x14, x15, cs
	csel x16, x19, x16, cs
	csel x22, x23, x22, cs
	csel x28, x26, x28, cs
	stp x22, x28, [x1, #16]
	stp x15, x16, [x1], #-64
	subs x0, x0, #1
	b.ne 1b
	ldr x5, [x2]
	ldp x6, x7, [sp, #80]
	ldp x8, x9, [sp, #96]
	ldp x10, x11, [x5, #32]
	ldp x12, x13, [x5, #48]
	ldp x14, x15, [x1, #0]
	ldp x16, x17, [x1, #16]
	subs x10, x10, x14
	sbcs x11, x11, x15
	sbcs x12, x12, x16
	sbcs x13, x13, x17
	csel x19, x4, xzr, lo
	subs x10, x10, x19
	sbcs x11, x11, xzr
	sbcs x12, x12, xzr
	sbc x13, x13, xzr
	mul x16, x10, x6
	umulh x15, x10, x6
	mul x17, x10, x7
	umulh x14, x10, x7
	adds x15, x15, x17
	mul x17, x10, x8
	umulh x19, x10, x8
	adcs x14, x14, x17
	mul x17, x10, x9
	umulh x20, x10, x9
	adcs x19, x19, x17
	adc x20, x20, xzr
	mul x10, x11, x6
	umulh x21, x11, x6
	mul x17, x11, x7
	umulh x22, x11, x7
	adds x21, x21, x17
	mul x17, x11, x8
	umulh x23, x11, x8
	adcs x22, x22, x17
	mul x17, x11, x9
	umulh x24, x11, x9
	adcs x23, x23, x17
	adc x24, x24, xzr
	adds x15, x15, x10
	adcs x14, x14, x21
	adcs x19, x19, x22
	adcs x20, x20, x23
	adc x24, x24, xzr
	mul x11, x12, x6
	umulh x23, x12, x6
	mul x17, x12, x7
	umulh x22, x12, x7
	adds x23, x23, x17
	mul x17, x12, x8
	umulh x21, x12, x8
	adcs x22, x22, x17
	mul x17, x12, x9
	umulh x10, x12, x9
	adcs x21, x21, x17
	adc x10, x10, xzr
	adds x14, x14, x11
	adcs x19, x19, x23
	adcs x20, x20, x22
	adcs x24, x24, x21
	adc x10, x10, xzr
	mul x12, x13, x6
	umulh x21, x13, x6
	mul x17, x13, x7
	umulh x22, x13, x7
	adds x21, x21, x17
	mul x17, x13, x8
	umulh x23, x13, x8
	adcs x22, x22, x17
	mul x17, x13, x9
	umulh x11, x13, x9
	adcs x23, x23, x17
	adc x11, x11, xzr
	adds x19, x19, x12
	adcs x20, x20, x21
	adcs x24, x24, x22
	adcs x10, x10, x23
	adc x11, x11, xzr
	mul x9, x20, x4
	umulh x17, x20, x4
	mul x8, x24, x4
	umulh x13, x24, x4
	mul x7, x10, x4
	umulh x23, x10, x4
	mul x6, x11, x4
	umulh x22, x11, x4
	adds x8, x8, x17
	adcs x7, x7, x13
	adcs x6, x6, x23
	adc x22, x22, xzr
	adds x16, x16, x9
	adcs x15, x15, x8
	adcs x14, x14, x7
	adcs x19, x19, x6
	adc x22, x22, xzr
	mul x6, x22, x4
	umulh x7, x22, x4
	adds x16, x16, x6
	adcs x15, x15, x7
	adcs x14, x14, xzr
	adcs x19, x19, xzr
	csel x6, x4, xzr, cs
	adds x16, x16, x6
	adcs x15, x15, xzr
	adcs x14, x14, xzr
	adc x19, x19, xzr
	mul x7, x16, x15
	umulh x6, x16, x15
	mul x22, x16, x14
	umulh x8, x16, x14
	mul x9, x16, x19
	umulh x23, x16, x19
	mul x13, x15, x14
	umulh x17, x15, x14
	mul x11, x15, x19
	umulh x10, x15, x19
	mul x24, x14, x19
	umulh x20, x14, x19
	adds x6, x6, x22
	adcs x8, x8, x9
	adc x23, x23, xzr
	adds x17, x17, x11
	adc x10, x10, xzr
	adds x8, x8, x13
	adcs x23, x23, x17
	adcs x10, x10, x24
	adc x20, x20, xzr
	lsr x24, x20, #63
	extr x20, x20, x10, #63
	extr x10, x10, x23, #63
	extr x23, x23, x8, #63
	extr x8, x8, x6, #63
	extr x6, x6, x7, #63
	lsl x7, x7, #1
	mul x17, x16, x16
	umulh x13, x16, x16
	mul x21, x15, x15
	umulh x11, x15, x15
	mul x12, x14, x14
	umulh x9, x14, x14
	mul x25, x19, x19
	umulh x22, x19, x19
	adds x7, x7, x13
	adcs x6, x6, x21
	adcs x8, x8, x11
	adcs x23, x23, x12
	adcs x10, x10, x9
	adcs x20, x20, x25
	adc x24, x24, x22
	mul x25, x23, x4
	umulh x9, x23, x4
	mul x12, x10, x4
	umulh x11, x10, x4
	mul x21, x20, x4
	umulh x13, x20, x4
	mul x22, x24, x4
	umulh x26, x24, x4
	adds x12, x12, x9
	adcs x21, x21, x11
	adcs x22, x22, x13
	adc x26, x26, xzr
	adds x17, x17, x25
	adcs x7, x7, x12
	adcs x6, x6, x21
	adcs x8, x8, x22
	adc x26, x26, xzr
	mul x22, x26, x4
	umulh x21, x26, x4
	adds x17, x17, x22
	adcs x7, x7, x21
	adcs x6, x6, xzr
	adcs x8, x8, xzr
	csel x22, x4, xzr, cs
	adds x17, x17, x22
	adcs x7, x7, xzr
	adcs x6, x6, xzr
	adc x8, x8, xzr
	ldp x21, x22, [x1, #-32]
	ldp x26, x12, [x1, #-16]
	subs x17, x17, x21
	sbcs x7, x7, x22
	sbcs x6, x6, x26
	sbcs x8, x8, x12
	csel x25, x4, xzr, lo
	subs x17, x17, x25
	sbcs x7, x7, xzr
	sbcs x6, x6, xzr
	sbc x8, x8, xzr
	ldp x12, x26, [x5, #0]
	ldp x22, x21, [x5, #16]
	subs x17, x17, x12
	sbcs x7, x7, x26
	sbcs x6, x6, x22
	sbcs x8, x8, x21
	csel x25, x4, xzr, lo
	subs x17, x17, x25
	sbcs x7, x7, xzr
	sbcs x6, x6, xzr
	sbc x8, x8, xzr
	adds x21, x17, x4
	adcs x22, x7, xzr
	adcs x26, x6, xzr
	adcs x12, x8, xzr
	csel x17, x21, x17, cs
	csel x7, x22, x7, cs
	csel x6, x26, x6, cs
	csel x8, x12, x8, cs
	ldp x12, x26, [x1, #-32]
	ldp x22, x21, [x1, #-16]
	subs x12, x12, x17
	sbcs x26, x26, x7
	sbcs x22, x22, x6
	sbcs x21, x21, x8
	csel x25, x4, xzr, lo
	subs x12, x12, x25
	sbcs x26, x26, xzr
	sbcs x22, x22, xzr
	sbc x21, x21, xzr
	stp x17, x7, [x1, #-32]
	stp x6, x8, [x1, #-16]
	mul x6, x16, x12
	umulh x7, x16, x12
	mul x8, x16, x26
	umulh x17, x16, x26
	adds x7, x7, x8
	mul x8, x16, x22
	umulh x25, x16, x22
	adcs x17, x17, x8
	mul x8, x16, x21
	umulh x13, x16, x21
	adcs x25, x25, x8
	adc x13, x13, xzr
	mul x16, x15, x12
	umulh x11, x15, x12
	mul x8, x15, x26
	umulh x9, x15, x26
	adds x11, x11, x8
	mul x8, x15, x22
	umulh x24, x15, x22
	adcs x9, x9, x8
	mul x8, x15, x21
	umulh x20, x15, x21
	adcs x24, x24, x8
	adc x20, x20, xzr
	adds x7, x7, x16
	adcs x17, x17, x11
	adcs x25, x25, x9
	adcs x13, x13, x24
	adc x20, x20, xzr
	mul x15, x14, x12
	umulh x24, x14, x12
	mul x8, x14, x26
	umulh x9, x14, x26
	adds x24, x24, x8
	mul x8, x14, x22
	umulh x11, x14, x22
	adcs x9, x9, x8
	mul x8, x14, x21
	umulh x16, x14, x21
	adcs x11, x11, x8
	adc x16, x16, xzr
	adds x17, x17, x15
	adcs x25, x25, x24
	adcs x13, x13, x9
	adcs x20, x20, x11
	adc x16, x16, xzr
	mul x14, x19, x12
	umulh x11, x19, x12
	mul x8, x19, x26
	umulh x9, x19, x26
	adds x11, x11, x8
	mul x8, x19, x22
	umulh x24, x19, x22
	adcs x9, x9, x8
	mul x8, x19, x21
	umulh x15, x19, x21
	adcs x24, x24, x8
	adc x15, x15, xzr
	adds x25, x25, x14
	adcs x13, x13, x11
	adcs x20, x20, x9
	adcs x16, x16, x24
	adc x15, x15, xzr
	mul x21, x13, x4
	umulh x8, x13, x4
	mul x22, x20, x4
	umulh x19, x20, x4
	mul x26, x16, x4
	umulh x24, x16, x4
	mul x12, x15, x4
	umulh x9, x15, x4
	adds x22, x22, x8
	adcs x26, x26, x19
	adcs x12, x12, x24
	adc x9, x9, xzr
	adds x6, x6, x21
	adcs x7, x7, x22
	adcs x17, x17, x26
	adcs x25, x25, x12
	adc x9, x9, xzr
	mul x12, x9, x4
	umulh x26, x9, x4
	adds x6, x6, x12
	adcs x7, x7, x26
	adcs x17, x17, xzr
	adcs x25, x25, xzr
	csel x12, x4, xzr, cs
	adds x6, x6, x12
	adcs x7, x7, xzr
	adcs x17, x17, xzr
	adc x25, x25, xzr
	ldp x26, x12, [x1, #0]
	ldp x9, x22, [x1, #16]
	subs x6, x6, x26
	sbcs x7, x7, x12
	sbcs x17, x17, x9
	sbcs x25, x25, x22
	csel x21, x4, xzr, lo
	subs x6, x6, x21
	sbcs x7, x7, xzr
	sbcs x17, x17, xzr
	sbc x25, x25, xzr
	adds x22, x6, x4
	adcs x9, x7, xzr
	adcs x12, x17, xzr
	adcs x26, x25, xzr
	csel x6, x22, x6, cs
	csel x7, x9, x7, cs
	csel x17, x12, x17, cs
	csel x25, x26, x25, cs
	stp x6, x7, [x1, #0]
	stp x17, x25, [x1, #16]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #112
	ret

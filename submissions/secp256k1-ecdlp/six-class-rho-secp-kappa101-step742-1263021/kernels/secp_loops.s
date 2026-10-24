	.text
	.p2align 2
	.globl pass1
	.type pass1, %function
pass1:
	stp x19, x20, [sp, #-96]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	mov x4, #977
	movk x4, #1, lsl #32
	ldr x9, [x2], #8
	ldp x5, x6, [x9, #0]
	ldp x7, x8, [x9, #16]
	ldp x12, x13, [x1, #16]
	ldp x10, x11, [x1], #64
	subs x5, x5, x10
	sbcs x6, x6, x11
	sbcs x7, x7, x12
	sbcs x8, x8, x13
	csel x14, x4, xzr, lo
	subs x5, x5, x14
	sbcs x6, x6, xzr
	sbcs x7, x7, xzr
	sbc x8, x8, xzr
	add x3, x3, #64
	sub x0, x0, #1
1:
	ldr x9, [x2], #8
	ldp x10, x11, [x9, #0]
	ldp x12, x13, [x9, #16]
	ldp x16, x17, [x1, #16]
	ldp x14, x15, [x1], #64
	subs x10, x10, x14
	sbcs x11, x11, x15
	sbcs x12, x12, x16
	sbcs x13, x13, x17
	csel x19, x4, xzr, lo
	subs x10, x10, x19
	sbcs x11, x11, xzr
	sbcs x12, x12, xzr
	sbc x13, x13, xzr
	stp x5, x6, [x3, #32]
	stp x7, x8, [x3, #48]
	stp x12, x13, [x3, #16]
	stp x10, x11, [x3], #64
	mul x17, x10, x5
	umulh x16, x10, x5
	mul x9, x10, x6
	umulh x15, x10, x6
	adds x16, x16, x9
	mul x9, x10, x7
	umulh x14, x10, x7
	adcs x15, x15, x9
	mul x9, x10, x8
	umulh x19, x10, x8
	adcs x14, x14, x9
	adc x19, x19, xzr
	mul x10, x11, x5
	umulh x20, x11, x5
	mul x9, x11, x6
	umulh x21, x11, x6
	adds x20, x20, x9
	mul x9, x11, x7
	umulh x22, x11, x7
	adcs x21, x21, x9
	mul x9, x11, x8
	umulh x23, x11, x8
	adcs x22, x22, x9
	adc x23, x23, xzr
	adds x16, x16, x10
	adcs x15, x15, x20
	adcs x14, x14, x21
	adcs x19, x19, x22
	adc x23, x23, xzr
	mul x11, x12, x5
	umulh x22, x12, x5
	mul x9, x12, x6
	umulh x21, x12, x6
	adds x22, x22, x9
	mul x9, x12, x7
	umulh x20, x12, x7
	adcs x21, x21, x9
	mul x9, x12, x8
	umulh x10, x12, x8
	adcs x20, x20, x9
	adc x10, x10, xzr
	adds x15, x15, x11
	adcs x14, x14, x22
	adcs x19, x19, x21
	adcs x23, x23, x20
	adc x10, x10, xzr
	mul x12, x13, x5
	umulh x20, x13, x5
	mul x9, x13, x6
	umulh x21, x13, x6
	adds x20, x20, x9
	mul x9, x13, x7
	umulh x22, x13, x7
	adcs x21, x21, x9
	mul x9, x13, x8
	umulh x11, x13, x8
	adcs x22, x22, x9
	adc x11, x11, xzr
	adds x14, x14, x12
	adcs x19, x19, x20
	adcs x23, x23, x21
	adcs x10, x10, x22
	adc x11, x11, xzr
	mul x9, x19, x4
	umulh x20, x19, x4
	mul x13, x23, x4
	umulh x12, x23, x4
	mul x22, x10, x4
	umulh x24, x10, x4
	mul x21, x11, x4
	umulh x25, x11, x4
	adds x13, x13, x20
	adcs x22, x22, x12
	adcs x21, x21, x24
	adc x25, x25, xzr
	adds x17, x17, x9
	adcs x16, x16, x13
	adcs x15, x15, x22
	adcs x14, x14, x21
	adc x25, x25, xzr
	mul x21, x25, x4
	umulh x22, x25, x4
	adds x17, x17, x21
	adcs x16, x16, x22
	adcs x15, x15, xzr
	adcs x14, x14, xzr
	csel x21, x4, xzr, cs
	adds x5, x17, x21
	adcs x6, x16, xzr
	adcs x7, x15, xzr
	adc x8, x14, xzr
	subs x0, x0, #1
	b.ne 1b
	stp x5, x6, [x3, #32]
	stp x7, x8, [x3, #48]
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #96
	ret
	.text
	.p2align 2
	.globl pass2
	.type pass2, %function
pass2:
	stp x19, x20, [sp, #-128]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	ldp x5, x6, [x4]
	ldp x7, x8, [x4, #16]
	stp x5, x6, [sp, #96]
	stp x7, x8, [sp, #112]
	mov x4, #977
	movk x4, #1, lsl #32
	sub x5, x0, #1
	add x1, x1, x5, lsl #6
	add x1, x1, #32
	add x2, x2, x5, lsl #3
	add x3, x3, x5, lsl #6
	mov x0, x5
1:
	ldr x5, [x2], #-8
	ldp x6, x7, [sp, #96]
	ldp x8, x9, [sp, #112]
	ldp x10, x11, [x3, #32]
	ldp x12, x13, [x3, #48]
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
	ldp x22, x12, [x3, #16]
	ldp x11, x10, [x3], #-64
	mul x23, x6, x11
	umulh x24, x6, x11
	mul x13, x6, x10
	umulh x14, x6, x10
	adds x24, x24, x13
	mul x13, x6, x22
	umulh x27, x6, x22
	adcs x14, x14, x13
	mul x13, x6, x12
	umulh x26, x6, x12
	adcs x27, x27, x13
	adc x26, x26, xzr
	mul x6, x7, x11
	umulh x25, x7, x11
	mul x13, x7, x10
	umulh x20, x7, x10
	adds x25, x25, x13
	mul x13, x7, x22
	umulh x21, x7, x22
	adcs x20, x20, x13
	mul x13, x7, x12
	umulh x28, x7, x12
	adcs x21, x21, x13
	adc x28, x28, xzr
	adds x24, x24, x6
	adcs x14, x14, x25
	adcs x27, x27, x20
	adcs x26, x26, x21
	adc x28, x28, xzr
	mul x7, x8, x11
	umulh x21, x8, x11
	mul x13, x8, x10
	umulh x20, x8, x10
	adds x21, x21, x13
	mul x13, x8, x22
	umulh x25, x8, x22
	adcs x20, x20, x13
	mul x13, x8, x12
	umulh x6, x8, x12
	adcs x25, x25, x13
	adc x6, x6, xzr
	adds x14, x14, x7
	adcs x27, x27, x21
	adcs x26, x26, x20
	adcs x28, x28, x25
	adc x6, x6, xzr
	mul x8, x9, x11
	umulh x25, x9, x11
	mul x13, x9, x10
	umulh x20, x9, x10
	adds x25, x25, x13
	mul x13, x9, x22
	umulh x21, x9, x22
	adcs x20, x20, x13
	mul x13, x9, x12
	umulh x7, x9, x12
	adcs x21, x21, x13
	adc x7, x7, xzr
	adds x27, x27, x8
	adcs x26, x26, x25
	adcs x28, x28, x20
	adcs x6, x6, x21
	adc x7, x7, xzr
	mul x12, x26, x4
	umulh x13, x26, x4
	mul x22, x28, x4
	umulh x9, x28, x4
	mul x10, x6, x4
	umulh x21, x6, x4
	mul x11, x7, x4
	umulh x20, x7, x4
	adds x22, x22, x13
	adcs x10, x10, x9
	adcs x11, x11, x21
	adc x20, x20, xzr
	adds x23, x23, x12
	adcs x24, x24, x22
	adcs x14, x14, x10
	adcs x27, x27, x11
	adc x20, x20, xzr
	mul x11, x20, x4
	umulh x10, x20, x4
	adds x23, x23, x11
	adcs x24, x24, x10
	adcs x14, x14, xzr
	adcs x27, x27, xzr
	csel x11, x4, xzr, cs
	adds x23, x23, x11
	adcs x24, x24, xzr
	adcs x14, x14, xzr
	adc x27, x27, xzr
	stp x23, x24, [sp, #96]
	stp x14, x27, [sp, #112]
	ldp x27, x14, [x5, #32]
	ldp x24, x23, [x5, #48]
	ldp x10, x11, [x1, #0]
	ldp x20, x22, [x1, #16]
	subs x27, x27, x10
	sbcs x14, x14, x11
	sbcs x24, x24, x20
	sbcs x23, x23, x22
	csel x12, x4, xzr, lo
	subs x27, x27, x12
	sbcs x14, x14, xzr
	sbcs x24, x24, xzr
	sbc x23, x23, xzr
	mul x21, x27, x15
	umulh x9, x27, x15
	mul x12, x27, x16
	umulh x13, x27, x16
	adds x9, x9, x12
	mul x12, x27, x17
	umulh x7, x27, x17
	adcs x13, x13, x12
	mul x12, x27, x19
	umulh x6, x27, x19
	adcs x7, x7, x12
	adc x6, x6, xzr
	mul x27, x14, x15
	umulh x28, x14, x15
	mul x12, x14, x16
	umulh x26, x14, x16
	adds x28, x28, x12
	mul x12, x14, x17
	umulh x25, x14, x17
	adcs x26, x26, x12
	mul x12, x14, x19
	umulh x8, x14, x19
	adcs x25, x25, x12
	adc x8, x8, xzr
	adds x9, x9, x27
	adcs x13, x13, x28
	adcs x7, x7, x26
	adcs x6, x6, x25
	adc x8, x8, xzr
	mul x14, x24, x15
	umulh x25, x24, x15
	mul x12, x24, x16
	umulh x26, x24, x16
	adds x25, x25, x12
	mul x12, x24, x17
	umulh x28, x24, x17
	adcs x26, x26, x12
	mul x12, x24, x19
	umulh x27, x24, x19
	adcs x28, x28, x12
	adc x27, x27, xzr
	adds x13, x13, x14
	adcs x7, x7, x25
	adcs x6, x6, x26
	adcs x8, x8, x28
	adc x27, x27, xzr
	mul x24, x23, x15
	umulh x28, x23, x15
	mul x12, x23, x16
	umulh x26, x23, x16
	adds x28, x28, x12
	mul x12, x23, x17
	umulh x25, x23, x17
	adcs x26, x26, x12
	mul x12, x23, x19
	umulh x14, x23, x19
	adcs x25, x25, x12
	adc x14, x14, xzr
	adds x7, x7, x24
	adcs x6, x6, x28
	adcs x8, x8, x26
	adcs x27, x27, x25
	adc x14, x14, xzr
	mul x19, x6, x4
	umulh x12, x6, x4
	mul x17, x8, x4
	umulh x23, x8, x4
	mul x16, x27, x4
	umulh x25, x27, x4
	mul x15, x14, x4
	umulh x26, x14, x4
	adds x17, x17, x12
	adcs x16, x16, x23
	adcs x15, x15, x25
	adc x26, x26, xzr
	adds x21, x21, x19
	adcs x9, x9, x17
	adcs x13, x13, x16
	adcs x7, x7, x15
	adc x26, x26, xzr
	mul x15, x26, x4
	umulh x16, x26, x4
	adds x21, x21, x15
	adcs x9, x9, x16
	adcs x13, x13, xzr
	adcs x7, x7, xzr
	csel x15, x4, xzr, cs
	adds x21, x21, x15
	adcs x9, x9, xzr
	adcs x13, x13, xzr
	adc x7, x7, xzr
	mul x16, x21, x9
	umulh x15, x21, x9
	mul x26, x21, x13
	umulh x17, x21, x13
	mul x19, x21, x7
	umulh x25, x21, x7
	mul x23, x9, x13
	umulh x12, x9, x13
	mul x14, x9, x7
	umulh x27, x9, x7
	mul x8, x13, x7
	umulh x6, x13, x7
	adds x15, x15, x26
	adcs x17, x17, x19
	adc x25, x25, xzr
	adds x12, x12, x14
	adc x27, x27, xzr
	adds x17, x17, x23
	adcs x25, x25, x12
	adcs x27, x27, x8
	adc x6, x6, xzr
	lsr x8, x6, #63
	extr x6, x6, x27, #63
	extr x27, x27, x25, #63
	extr x25, x25, x17, #63
	extr x17, x17, x15, #63
	extr x15, x15, x16, #63
	lsl x16, x16, #1
	mul x12, x21, x21
	umulh x23, x21, x21
	mul x28, x9, x9
	umulh x14, x9, x9
	mul x24, x13, x13
	umulh x19, x13, x13
	mul x29, x7, x7
	umulh x26, x7, x7
	adds x16, x16, x23
	adcs x15, x15, x28
	adcs x17, x17, x14
	adcs x25, x25, x24
	adcs x27, x27, x19
	adcs x6, x6, x29
	adc x8, x8, x26
	mul x29, x25, x4
	umulh x19, x25, x4
	mul x24, x27, x4
	umulh x14, x27, x4
	mul x28, x6, x4
	umulh x23, x6, x4
	mul x26, x8, x4
	umulh x30, x8, x4
	adds x24, x24, x19
	adcs x28, x28, x14
	adcs x26, x26, x23
	adc x30, x30, xzr
	adds x12, x12, x29
	adcs x16, x16, x24
	adcs x15, x15, x28
	adcs x17, x17, x26
	adc x30, x30, xzr
	mul x26, x30, x4
	umulh x28, x30, x4
	adds x12, x12, x26
	adcs x16, x16, x28
	adcs x15, x15, xzr
	adcs x17, x17, xzr
	csel x26, x4, xzr, cs
	adds x12, x12, x26
	adcs x16, x16, xzr
	adcs x15, x15, xzr
	adc x17, x17, xzr
	ldp x28, x26, [x1, #-32]
	ldp x30, x24, [x1, #-16]
	subs x12, x12, x28
	sbcs x16, x16, x26
	sbcs x15, x15, x30
	sbcs x17, x17, x24
	csel x29, x4, xzr, lo
	subs x12, x12, x29
	sbcs x16, x16, xzr
	sbcs x15, x15, xzr
	sbc x17, x17, xzr
	ldp x29, x23, [x5, #0]
	ldp x14, x19, [x5, #16]
	subs x12, x12, x29
	sbcs x16, x16, x23
	sbcs x15, x15, x14
	sbcs x17, x17, x19
	csel x8, x4, xzr, lo
	subs x12, x12, x8
	sbcs x16, x16, xzr
	sbcs x15, x15, xzr
	sbc x17, x17, xzr
	adds x19, x12, x4
	adcs x14, x16, xzr
	adcs x23, x15, xzr
	adcs x29, x17, xzr
	csel x12, x19, x12, cs
	csel x16, x14, x16, cs
	csel x15, x23, x15, cs
	csel x17, x29, x17, cs
	subs x28, x28, x12
	sbcs x26, x26, x16
	sbcs x30, x30, x15
	sbcs x24, x24, x17
	csel x29, x4, xzr, lo
	subs x28, x28, x29
	sbcs x26, x26, xzr
	sbcs x30, x30, xzr
	sbc x24, x24, xzr
	stp x12, x16, [x1, #-32]
	stp x15, x17, [x1, #-16]
	mul x15, x21, x28
	umulh x16, x21, x28
	mul x17, x21, x26
	umulh x12, x21, x26
	adds x16, x16, x17
	mul x17, x21, x30
	umulh x29, x21, x30
	adcs x12, x12, x17
	mul x17, x21, x24
	umulh x23, x21, x24
	adcs x29, x29, x17
	adc x23, x23, xzr
	mul x21, x9, x28
	umulh x14, x9, x28
	mul x17, x9, x26
	umulh x19, x9, x26
	adds x14, x14, x17
	mul x17, x9, x30
	umulh x8, x9, x30
	adcs x19, x19, x17
	mul x17, x9, x24
	umulh x6, x9, x24
	adcs x8, x8, x17
	adc x6, x6, xzr
	adds x16, x16, x21
	adcs x12, x12, x14
	adcs x29, x29, x19
	adcs x23, x23, x8
	adc x6, x6, xzr
	mul x9, x13, x28
	umulh x8, x13, x28
	mul x17, x13, x26
	umulh x19, x13, x26
	adds x8, x8, x17
	mul x17, x13, x30
	umulh x14, x13, x30
	adcs x19, x19, x17
	mul x17, x13, x24
	umulh x21, x13, x24
	adcs x14, x14, x17
	adc x21, x21, xzr
	adds x12, x12, x9
	adcs x29, x29, x8
	adcs x23, x23, x19
	adcs x6, x6, x14
	adc x21, x21, xzr
	mul x13, x7, x28
	umulh x14, x7, x28
	mul x17, x7, x26
	umulh x19, x7, x26
	adds x14, x14, x17
	mul x17, x7, x30
	umulh x8, x7, x30
	adcs x19, x19, x17
	mul x17, x7, x24
	umulh x9, x7, x24
	adcs x8, x8, x17
	adc x9, x9, xzr
	adds x29, x29, x13
	adcs x23, x23, x14
	adcs x6, x6, x19
	adcs x21, x21, x8
	adc x9, x9, xzr
	mul x24, x23, x4
	umulh x17, x23, x4
	mul x30, x6, x4
	umulh x7, x6, x4
	mul x26, x21, x4
	umulh x8, x21, x4
	mul x28, x9, x4
	umulh x19, x9, x4
	adds x30, x30, x17
	adcs x26, x26, x7
	adcs x28, x28, x8
	adc x19, x19, xzr
	adds x15, x15, x24
	adcs x16, x16, x30
	adcs x12, x12, x26
	adcs x29, x29, x28
	adc x19, x19, xzr
	mul x28, x19, x4
	umulh x26, x19, x4
	adds x15, x15, x28
	adcs x16, x16, x26
	adcs x12, x12, xzr
	adcs x29, x29, xzr
	csel x28, x4, xzr, cs
	adds x15, x15, x28
	adcs x16, x16, xzr
	adcs x12, x12, xzr
	adc x29, x29, xzr
	subs x15, x15, x10
	sbcs x16, x16, x11
	sbcs x12, x12, x20
	sbcs x29, x29, x22
	csel x26, x4, xzr, lo
	subs x15, x15, x26
	sbcs x16, x16, xzr
	sbcs x12, x12, xzr
	sbc x29, x29, xzr
	adds x22, x15, x4
	adcs x20, x16, xzr
	adcs x11, x12, xzr
	adcs x10, x29, xzr
	csel x15, x22, x15, cs
	csel x16, x20, x16, cs
	csel x12, x11, x12, cs
	csel x29, x10, x29, cs
	stp x12, x29, [x1, #16]
	stp x15, x16, [x1], #-64
	subs x0, x0, #1
	b.ne 1b
	ldr x5, [x2]
	ldp x6, x7, [sp, #96]
	ldp x8, x9, [sp, #112]
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
	mul x20, x10, x6
	umulh x21, x10, x6
	mul x19, x10, x7
	umulh x22, x10, x7
	adds x21, x21, x19
	mul x19, x10, x8
	umulh x23, x10, x8
	adcs x22, x22, x19
	mul x19, x10, x9
	umulh x24, x10, x9
	adcs x23, x23, x19
	adc x24, x24, xzr
	mul x10, x11, x6
	umulh x25, x11, x6
	mul x19, x11, x7
	umulh x26, x11, x7
	adds x25, x25, x19
	mul x19, x11, x8
	umulh x27, x11, x8
	adcs x26, x26, x19
	mul x19, x11, x9
	umulh x28, x11, x9
	adcs x27, x27, x19
	adc x28, x28, xzr
	adds x21, x21, x10
	adcs x22, x22, x25
	adcs x23, x23, x26
	adcs x24, x24, x27
	adc x28, x28, xzr
	mul x11, x12, x6
	umulh x27, x12, x6
	mul x19, x12, x7
	umulh x26, x12, x7
	adds x27, x27, x19
	mul x19, x12, x8
	umulh x25, x12, x8
	adcs x26, x26, x19
	mul x19, x12, x9
	umulh x10, x12, x9
	adcs x25, x25, x19
	adc x10, x10, xzr
	adds x22, x22, x11
	adcs x23, x23, x27
	adcs x24, x24, x26
	adcs x28, x28, x25
	adc x10, x10, xzr
	mul x12, x13, x6
	umulh x25, x13, x6
	mul x19, x13, x7
	umulh x26, x13, x7
	adds x25, x25, x19
	mul x19, x13, x8
	umulh x27, x13, x8
	adcs x26, x26, x19
	mul x19, x13, x9
	umulh x11, x13, x9
	adcs x27, x27, x19
	adc x11, x11, xzr
	adds x23, x23, x12
	adcs x24, x24, x25
	adcs x28, x28, x26
	adcs x10, x10, x27
	adc x11, x11, xzr
	mul x9, x24, x4
	umulh x19, x24, x4
	mul x8, x28, x4
	umulh x13, x28, x4
	mul x7, x10, x4
	umulh x27, x10, x4
	mul x6, x11, x4
	umulh x26, x11, x4
	adds x8, x8, x19
	adcs x7, x7, x13
	adcs x6, x6, x27
	adc x26, x26, xzr
	adds x20, x20, x9
	adcs x21, x21, x8
	adcs x22, x22, x7
	adcs x23, x23, x6
	adc x26, x26, xzr
	mul x6, x26, x4
	umulh x7, x26, x4
	adds x20, x20, x6
	adcs x21, x21, x7
	adcs x22, x22, xzr
	adcs x23, x23, xzr
	csel x6, x4, xzr, cs
	adds x20, x20, x6
	adcs x21, x21, xzr
	adcs x22, x22, xzr
	adc x23, x23, xzr
	mul x7, x20, x21
	umulh x6, x20, x21
	mul x26, x20, x22
	umulh x8, x20, x22
	mul x9, x20, x23
	umulh x27, x20, x23
	mul x13, x21, x22
	umulh x19, x21, x22
	mul x11, x21, x23
	umulh x10, x21, x23
	mul x28, x22, x23
	umulh x24, x22, x23
	adds x6, x6, x26
	adcs x8, x8, x9
	adc x27, x27, xzr
	adds x19, x19, x11
	adc x10, x10, xzr
	adds x8, x8, x13
	adcs x27, x27, x19
	adcs x10, x10, x28
	adc x24, x24, xzr
	lsr x28, x24, #63
	extr x24, x24, x10, #63
	extr x10, x10, x27, #63
	extr x27, x27, x8, #63
	extr x8, x8, x6, #63
	extr x6, x6, x7, #63
	lsl x7, x7, #1
	mul x19, x20, x20
	umulh x13, x20, x20
	mul x25, x21, x21
	umulh x11, x21, x21
	mul x12, x22, x22
	umulh x9, x22, x22
	mul x29, x23, x23
	umulh x26, x23, x23
	adds x7, x7, x13
	adcs x6, x6, x25
	adcs x8, x8, x11
	adcs x27, x27, x12
	adcs x10, x10, x9
	adcs x24, x24, x29
	adc x28, x28, x26
	mul x29, x27, x4
	umulh x9, x27, x4
	mul x12, x10, x4
	umulh x11, x10, x4
	mul x25, x24, x4
	umulh x13, x24, x4
	mul x26, x28, x4
	umulh x30, x28, x4
	adds x12, x12, x9
	adcs x25, x25, x11
	adcs x26, x26, x13
	adc x30, x30, xzr
	adds x19, x19, x29
	adcs x7, x7, x12
	adcs x6, x6, x25
	adcs x8, x8, x26
	adc x30, x30, xzr
	mul x26, x30, x4
	umulh x25, x30, x4
	adds x19, x19, x26
	adcs x7, x7, x25
	adcs x6, x6, xzr
	adcs x8, x8, xzr
	csel x26, x4, xzr, cs
	adds x19, x19, x26
	adcs x7, x7, xzr
	adcs x6, x6, xzr
	adc x8, x8, xzr
	ldp x25, x26, [x1, #-32]
	ldp x30, x12, [x1, #-16]
	subs x19, x19, x25
	sbcs x7, x7, x26
	sbcs x6, x6, x30
	sbcs x8, x8, x12
	csel x29, x4, xzr, lo
	subs x19, x19, x29
	sbcs x7, x7, xzr
	sbcs x6, x6, xzr
	sbc x8, x8, xzr
	ldp x29, x13, [x5, #0]
	ldp x11, x9, [x5, #16]
	subs x19, x19, x29
	sbcs x7, x7, x13
	sbcs x6, x6, x11
	sbcs x8, x8, x9
	csel x28, x4, xzr, lo
	subs x19, x19, x28
	sbcs x7, x7, xzr
	sbcs x6, x6, xzr
	sbc x8, x8, xzr
	adds x9, x19, x4
	adcs x11, x7, xzr
	adcs x13, x6, xzr
	adcs x29, x8, xzr
	csel x19, x9, x19, cs
	csel x7, x11, x7, cs
	csel x6, x13, x6, cs
	csel x8, x29, x8, cs
	subs x25, x25, x19
	sbcs x26, x26, x7
	sbcs x30, x30, x6
	sbcs x12, x12, x8
	csel x29, x4, xzr, lo
	subs x25, x25, x29
	sbcs x26, x26, xzr
	sbcs x30, x30, xzr
	sbc x12, x12, xzr
	stp x19, x7, [x1, #-32]
	stp x6, x8, [x1, #-16]
	mul x6, x20, x25
	umulh x7, x20, x25
	mul x8, x20, x26
	umulh x19, x20, x26
	adds x7, x7, x8
	mul x8, x20, x30
	umulh x29, x20, x30
	adcs x19, x19, x8
	mul x8, x20, x12
	umulh x13, x20, x12
	adcs x29, x29, x8
	adc x13, x13, xzr
	mul x20, x21, x25
	umulh x11, x21, x25
	mul x8, x21, x26
	umulh x9, x21, x26
	adds x11, x11, x8
	mul x8, x21, x30
	umulh x28, x21, x30
	adcs x9, x9, x8
	mul x8, x21, x12
	umulh x24, x21, x12
	adcs x28, x28, x8
	adc x24, x24, xzr
	adds x7, x7, x20
	adcs x19, x19, x11
	adcs x29, x29, x9
	adcs x13, x13, x28
	adc x24, x24, xzr
	mul x21, x22, x25
	umulh x28, x22, x25
	mul x8, x22, x26
	umulh x9, x22, x26
	adds x28, x28, x8
	mul x8, x22, x30
	umulh x11, x22, x30
	adcs x9, x9, x8
	mul x8, x22, x12
	umulh x20, x22, x12
	adcs x11, x11, x8
	adc x20, x20, xzr
	adds x19, x19, x21
	adcs x29, x29, x28
	adcs x13, x13, x9
	adcs x24, x24, x11
	adc x20, x20, xzr
	mul x22, x23, x25
	umulh x11, x23, x25
	mul x8, x23, x26
	umulh x9, x23, x26
	adds x11, x11, x8
	mul x8, x23, x30
	umulh x28, x23, x30
	adcs x9, x9, x8
	mul x8, x23, x12
	umulh x21, x23, x12
	adcs x28, x28, x8
	adc x21, x21, xzr
	adds x29, x29, x22
	adcs x13, x13, x11
	adcs x24, x24, x9
	adcs x20, x20, x28
	adc x21, x21, xzr
	mul x12, x13, x4
	umulh x8, x13, x4
	mul x30, x24, x4
	umulh x23, x24, x4
	mul x26, x20, x4
	umulh x28, x20, x4
	mul x25, x21, x4
	umulh x9, x21, x4
	adds x30, x30, x8
	adcs x26, x26, x23
	adcs x25, x25, x28
	adc x9, x9, xzr
	adds x6, x6, x12
	adcs x7, x7, x30
	adcs x19, x19, x26
	adcs x29, x29, x25
	adc x9, x9, xzr
	mul x25, x9, x4
	umulh x26, x9, x4
	adds x6, x6, x25
	adcs x7, x7, x26
	adcs x19, x19, xzr
	adcs x29, x29, xzr
	csel x25, x4, xzr, cs
	adds x6, x6, x25
	adcs x7, x7, xzr
	adcs x19, x19, xzr
	adc x29, x29, xzr
	subs x6, x6, x14
	sbcs x7, x7, x15
	sbcs x19, x19, x16
	sbcs x29, x29, x17
	csel x26, x4, xzr, lo
	subs x6, x6, x26
	sbcs x7, x7, xzr
	sbcs x19, x19, xzr
	sbc x29, x29, xzr
	adds x17, x6, x4
	adcs x16, x7, xzr
	adcs x15, x19, xzr
	adcs x14, x29, xzr
	csel x6, x17, x6, cs
	csel x7, x16, x7, cs
	csel x19, x15, x19, cs
	csel x29, x14, x29, cs
	stp x6, x7, [x1, #0]
	stp x19, x29, [x1, #16]
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #128
	ret

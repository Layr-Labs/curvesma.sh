	.text
	.p2align 2
	.globl t_mul
	.type t_mul, %function
t_mul:
	stp x19, x20, [sp, #-96]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	mov x4, #977
	movk x4, #1, lsl #32
	ldp x5, x6, [x1, #0]
	ldp x7, x8, [x1, #16]
	ldp x9, x10, [x2, #0]
	ldp x11, x12, [x2, #16]
	mul x13, x5, x9
	umulh x14, x5, x9
	mul x15, x5, x10
	adds x14, x14, x15
	umulh x15, x5, x10
	mul x16, x5, x11
	adcs x15, x15, x16
	umulh x16, x5, x11
	mul x17, x5, x12
	adcs x16, x16, x17
	umulh x17, x5, x12
	mul x5, x6, x12
	adcs x17, x17, x5
	umulh x5, x6, x12
	adc x5, x5, xzr
	mul x19, x6, x9
	adds x14, x14, x19
	umulh x19, x6, x9
	adcs x15, x15, x19
	umulh x19, x7, x9
	adcs x16, x16, x19
	umulh x19, x8, x9
	adcs x17, x17, x19
	umulh x19, x7, x11
	adc x19, x19, xzr
	mul x20, x7, x9
	adds x15, x15, x20
	mul x20, x8, x9
	adcs x16, x16, x20
	mul x20, x7, x11
	adcs x17, x17, x20
	adcs x5, x5, x19
	umulh x19, x7, x12
	adc x19, x19, xzr
	mul x20, x6, x10
	adds x15, x15, x20
	umulh x20, x6, x10
	adcs x16, x16, x20
	umulh x20, x6, x11
	adcs x17, x17, x20
	umulh x20, x8, x10
	adcs x5, x5, x20
	umulh x20, x8, x11
	adc x20, x20, xzr
	mul x9, x6, x11
	adds x16, x16, x9
	umulh x9, x7, x10
	adcs x17, x17, x9
	mul x9, x7, x12
	adcs x5, x5, x9
	adcs x19, x19, x20
	umulh x20, x8, x12
	adc x20, x20, xzr
	mul x9, x7, x10
	adds x16, x16, x9
	mul x9, x8, x10
	adcs x17, x17, x9
	mul x9, x8, x11
	adcs x5, x5, x9
	mul x9, x8, x12
	adcs x19, x19, x9
	adc x20, x20, xzr
	mul x9, x17, x4
	umulh x12, x17, x4
	mul x17, x5, x4
	adds x17, x17, x12
	umulh x12, x5, x4
	mul x5, x19, x4
	adcs x5, x5, x12
	umulh x12, x19, x4
	mul x19, x20, x4
	adcs x19, x19, x12
	umulh x12, x20, x4
	adc x12, x12, xzr
	adds x13, x13, x9
	adcs x14, x14, x17
	adcs x15, x15, x5
	adcs x16, x16, x19
	adc x12, x12, xzr
	mul x19, x12, x4
	umulh x5, x12, x4
	adds x13, x13, x19
	adcs x14, x14, x5
	adcs x15, x15, xzr
	adc x16, x16, xzr
	stp x13, x14, [x0, #0]
	stp x15, x16, [x0, #16]
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #96
	ret
	.text
	.p2align 2
	.globl t_sqr
	.type t_sqr, %function
t_sqr:
	stp x19, x20, [sp, #-96]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	mov x4, #977
	movk x4, #1, lsl #32
	ldp x5, x6, [x1, #0]
	ldp x7, x8, [x1, #16]
	mul x9, x5, x8
	mul x10, x6, x7
	adds x9, x9, x10
	umulh x10, x5, x8
	umulh x11, x6, x7
	adcs x10, x10, x11
	umulh x11, x6, x8
	adc x11, x11, xzr
	umulh x12, x5, x6
	mul x13, x5, x7
	adds x12, x12, x13
	umulh x13, x5, x7
	adcs x9, x9, x13
	mul x13, x6, x8
	adcs x10, x10, x13
	mul x13, x7, x8
	adcs x11, x11, x13
	umulh x13, x7, x8
	adc x13, x13, xzr
	mul x14, x5, x6
	lsr x15, x13, #63
	extr x13, x13, x11, #63
	extr x11, x11, x10, #63
	extr x10, x10, x9, #63
	extr x9, x9, x12, #63
	extr x12, x12, x14, #63
	mul x16, x5, x5
	umulh x17, x5, x5
	adds x14, x17, x14, lsl #1
	mul x17, x6, x6
	adcs x12, x12, x17
	umulh x17, x6, x6
	adcs x9, x9, x17
	mul x17, x7, x7
	adcs x10, x10, x17
	umulh x17, x7, x7
	adcs x11, x11, x17
	mul x17, x8, x8
	adcs x13, x13, x17
	umulh x17, x8, x8
	adc x15, x15, x17
	mul x17, x10, x4
	umulh x8, x10, x4
	mul x10, x11, x4
	adds x10, x10, x8
	umulh x8, x11, x4
	mul x11, x13, x4
	adcs x11, x11, x8
	umulh x8, x13, x4
	mul x13, x15, x4
	adcs x13, x13, x8
	umulh x8, x15, x4
	adc x8, x8, xzr
	adds x16, x16, x17
	adcs x14, x14, x10
	adcs x12, x12, x11
	adcs x9, x9, x13
	adc x8, x8, xzr
	mul x13, x8, x4
	umulh x11, x8, x4
	adds x16, x16, x13
	adcs x14, x14, x11
	adcs x12, x12, xzr
	adc x9, x9, xzr
	stp x16, x14, [x0, #0]
	stp x12, x9, [x0, #16]
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #96
	ret
	.text
	.p2align 2
	.globl t_sub
	.type t_sub, %function
t_sub:
	stp x19, x20, [sp, #-96]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	mov x4, #977
	movk x4, #1, lsl #32
	ldp x5, x6, [x1, #0]
	ldp x7, x8, [x1, #16]
	ldp x9, x10, [x2, #0]
	ldp x11, x12, [x2, #16]
	subs x5, x5, x9
	sbcs x6, x6, x10
	sbcs x7, x7, x11
	sbcs x8, x8, x12
	csel x12, x4, xzr, lo
	subs x5, x5, x12
	sbcs x6, x6, xzr
	sbcs x7, x7, xzr
	sbc x8, x8, xzr
	stp x5, x6, [x0, #0]
	stp x7, x8, [x0, #16]
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #96
	ret
	.text
	.p2align 2
	.globl t_reduce
	.type t_reduce, %function
t_reduce:
	stp x19, x20, [sp, #-96]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	mov x4, #977
	movk x4, #1, lsl #32
	ldp x5, x6, [x1, #0]
	ldp x7, x8, [x1, #16]
	ldp x9, x10, [x1, #32]
	ldp x11, x12, [x1, #48]
	mul x13, x9, x4
	umulh x14, x9, x4
	mul x9, x10, x4
	adds x9, x9, x14
	umulh x14, x10, x4
	mul x10, x11, x4
	adcs x10, x10, x14
	umulh x14, x11, x4
	mul x11, x12, x4
	adcs x11, x11, x14
	umulh x14, x12, x4
	adc x14, x14, xzr
	adds x5, x5, x13
	adcs x6, x6, x9
	adcs x7, x7, x10
	adcs x8, x8, x11
	adc x14, x14, xzr
	mul x11, x14, x4
	umulh x10, x14, x4
	adds x5, x5, x11
	adcs x6, x6, x10
	adcs x7, x7, xzr
	adc x8, x8, xzr
	stp x5, x6, [x0, #0]
	stp x7, x8, [x0, #16]
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #96
	ret

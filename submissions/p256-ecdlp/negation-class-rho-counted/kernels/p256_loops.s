// clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S p256_walk.c
// Apple clang 21.0.0 (LLVM). Comments and assembler directives stripped; functions: pass1, pass2, pass1_dbl, pass2_dbl
pass1:
	cbz x5, .LBB4_4
	sub sp, sp, #112
	mov x8, x0
	mov x0, xzr
	add x9, x2, #16
	stp x1, x30, [sp, #8]
	stp x28, x27, [sp, #32]
	stp x26, x25, [sp, #48]
	stp x24, x23, [sp, #64]
	stp x22, x21, [sp, #80]
	stp x20, x19, [sp, #96]
	str x4, [sp]
.LBB4_2:
	ldp x13, x12, [x8, #32]
	mov w10, #-1
	ldp x15, x16, [x8, #48]
	mov x11, #-4294967295
	str x5, [sp, #24]
	subs x14, x10, x12
	mvn x18, x13
	mov w10, #-1
	ngc x17, xzr
	subs x2, x17, x15
	sbc x17, x17, x16
	cmp x18, x13
	ngc x13, xzr
	subs x12, x13, x12
	sbc x13, x13, xzr
	cmn x12, x14
	cinc x12, x13, hs
	subs x13, x2, x15
	asr x14, x12, #63
	ngc x15, xzr
	cmn x13, x12
	add x12, x17, x11
	adc x13, x15, x14
	subs x12, x12, x16
	asr x14, x13, #63
	ngc x15, xzr
	cmn x12, x13
	adc x13, x15, x14
	ldp x12, x14, [x8]
	and x13, x13, #0x1
	bfi x13, x12, #1, #11
	add x15, x1, x13, lsl #6
	str x13, [x8, #160]
	ldp x13, x16, [x15]
	ldp x17, x15, [x15, #16]
	subs x12, x12, x13
	ldp x13, x18, [x8, #16]
	sbcs x14, x14, x16
	add x8, x8, #168
	sbcs x13, x13, x17
	sbcs x16, x18, x15
	cset w15, lo
	csel x17, x10, xzr, lo
	csel x18, x11, xzr, lo
	sbfx x15, x15, #0, #1
	adds x15, x12, x15
	adcs x14, x14, x17
	adcs x13, x13, xzr
	stp x15, x14, [x9, #-16]
	adc x12, x18, x16
	stp x13, x12, [x9], #32
	ldp x16, x6, [x4]
	ldp x17, x2, [x4, #16]
	umulh x5, x15, x16
	mul x21, x14, x16
	umulh x22, x14, x16
	mul x26, x15, x6
	adds x5, x5, x21
	umulh x24, x15, x6
	cinc x22, x22, hs
	mul x7, x13, x16
	adds x21, x5, x26
	umulh x23, x13, x16
	adcs x22, x22, x24
	mul x19, x14, x6
	umulh x20, x14, x6
	cinc x23, x23, hs
	adds x7, x22, x7
	mul x27, x15, x17
	adcs x22, x23, xzr
	mrs x23, NZCV
	umulh x28, x15, x17
	mul x30, x12, x16
	mul x1, x13, x6
	umulh x11, x13, x6
	mul x24, x12, x6
	umulh x26, x12, x6
	adds x6, x7, x19
	adcs x7, x22, x20
	umulh x10, x12, x16
	cset w20, hs
	msr NZCV, x23
	cinc x20, x20, hs
	adds x6, x6, x27
	adcs x7, x7, x28
	mul x4, x14, x17
	cinc x20, x20, hs
	adds x7, x7, x30
	mul x16, x15, x16
	adcs x10, x20, x10
	mrs x20, NZCV
	umulh x25, x14, x17
	adds x1, x7, x1
	mul x18, x15, x2
	adcs x10, x10, x11
	umulh x5, x15, x2
	orr x15, x14, x15
	orr x15, x15, x13
	orr x11, x15, x12
	cset w15, hs
	msr NZCV, x20
	cinc x15, x15, hs
	cmp x11, #0
	mul x19, x13, x17
	cinc x0, x0, eq
	adds x4, x1, x4
	adcs x10, x10, x25
	umulh x22, x13, x17
	cset w1, hs
	adds x18, x4, x18
	mul x23, x14, x2
	adcs x10, x10, x5
	adc x15, x15, x1
	adds x10, x10, x24
	umulh x14, x14, x2
	adcs x15, x15, x26
	cset w1, hs
	adds x10, x10, x19
	mul x11, x12, x17
	adcs x15, x15, x22
	cset w7, hs
	adds x10, x10, x23
	umulh x17, x12, x17
	adcs x14, x15, x14
	mul x4, x13, x2
	adc x15, x7, x1
	mov x7, #-4294967295
	adds x11, x14, x11
	umulh x5, x12, x2
	adcs x14, x15, x17
	lsl x17, x16, #32
	umulh x13, x13, x2
	mul x12, x12, x2
	cinc x15, x5, hs
	adds x11, x11, x4
	ldr x5, [sp, #24]
	umulh x2, x16, x7
	ldr x4, [sp]
	adcs x13, x14, x13
	lsr x14, x16, #32
	sub x16, x16, x17
	cset w1, hs
	adds x12, x13, x12
	adc x13, x15, x1
	adds x15, x21, x17
	adcs x14, x6, x14
	lsl x17, x15, #32
	mov x6, #-4294967295
	adcs x16, x18, x16
	umulh x18, x15, x7
	adcs x10, x10, x2
	lsr x2, x15, #32
	sub x15, x15, x17
	adcs x11, x11, xzr
	adcs x12, x12, xzr
	cset w1, hs
	adds x14, x14, x17
	adcs x16, x16, x2
	umulh x17, x14, x6
	adcs x10, x10, x15
	lsl x15, x14, #32
	adcs x11, x11, x18
	lsr x18, x14, #32
	adcs x12, x12, xzr
	sub x14, x14, x15
	adcs x13, x13, x1
	ldr x1, [sp, #8]
	cset w2, hs
	adds x16, x16, x15
	adcs x10, x10, x18
	umulh x15, x16, x6
	adcs x11, x11, x14
	lsl x14, x16, #32
	adcs x12, x12, x17
	lsr x17, x16, #32
	adcs x13, x13, xzr
	cset w18, hs
	adds x10, x10, x14
	sub x14, x16, x14
	adcs x11, x11, x17
	mov w16, #-1
	eor w18, w2, w18
	adcs x12, x12, x14
	adcs x13, x13, x15
	cset w14, hs
	adds x15, x10, #1
	sbcs x16, x11, x16
	eor w14, w18, w14
	cset w17, lo
	eor w14, w14, #0x1
	subs x17, x12, x17
	sbcs x18, x13, x6
	cset w2, lo
	tst w2, w14
	csel x10, x10, x15, ne
	csel x11, x11, x16, ne
	csel x12, x12, x17, ne
	csel x13, x13, x18, ne
	subs x5, x5, #1
	stp x10, x11, [x4]
	stp x12, x13, [x4, #16]
	stp x12, x13, [x3, #16]
	stp x10, x11, [x3], #32
	b.ne .LBB4_2
	ldp x20, x19, [sp, #96]
	ldr x30, [sp, #16]
	ldp x22, x21, [sp, #80]
	ldp x24, x23, [sp, #64]
	ldp x26, x25, [sp, #48]
	ldp x28, x27, [sp, #32]
	add sp, sp, #112
	ret
.LBB4_4:
	mov x0, xzr
	ret
.Lfunc_end4:
pass2:
	sub sp, sp, #288
	stp x29, x30, [sp, #192]
	add x29, sp, #192
	stp x28, x27, [sp, #208]
	stp x26, x25, [sp, #224]
	stp x24, x23, [sp, #240]
	stp x22, x21, [sp, #256]
	stp x20, x19, [sp, #272]
	str x7, [sp, #24]
	stur x5, [x29, #-24]
	stp x1, x4, [sp, #8]
	cbz x6, .LBB5_8
	mov w9, #168
	mov w10, #1
	stur xzr, [x29, #-16]
	madd x8, x6, x9, x0
	lsl x9, x6, #5
	sub x11, x10, x6
	add x10, x9, x3
	add x9, x9, x2
	sub x10, x10, #64
	sub x12, x9, #16
	sub x27, x8, #168
.LBB5_2:
	ldr x8, [sp, #16]
	cmp x11, #0
	stp x10, x11, [x29, #-40]
	ldp x16, x5, [x12, #-16]
	csel x9, x8, x10, eq
	ldur x10, [x29, #-24]
	ldp x17, x3, [x9, #16]
	ldp x11, x8, [x10]
	ldp x25, x2, [x9]
	ldp x13, x14, [x10, #16]
	mul x9, x17, x11
	ldp x4, x7, [x12]
	umulh x15, x25, x11
	stp x4, x9, [x29, #-72]
	mul x6, x2, x11
	umulh x9, x16, x11
	umulh x20, x2, x11
	adds x19, x15, x6
	mul x22, x25, x8
	stp x9, x12, [x29, #-56]
	umulh x9, x5, x11
	cinc x26, x20, hs
	umulh x12, x3, x11
	umulh x21, x25, x8
	str x9, [sp, #88]
	mul x9, x16, x8
	ldr x15, [sp, #88]
	str x12, [sp, #48]
	adds x12, x19, x22
	umulh x28, x17, x11
	stur x12, [x29, #-80]
	ldur x12, [x29, #-56]
	adcs x26, x26, x21
	mul x24, x5, x11
	str x9, [sp, #80]
	umulh x9, x16, x8
	cinc x28, x28, hs
	mul x30, x4, x11
	adds x24, x12, x24
	cinc x20, x15, hs
	ldr x15, [sp, #80]
	umulh x23, x4, x11
	str x9, [sp, #72]
	mul x9, x5, x8
	adds x15, x24, x15
	stur x15, [x29, #-56]
	ldr x15, [sp, #72]
	umulh x10, x2, x8
	adcs x20, x20, x15
	mov x15, x4
	cinc x23, x23, hs
	adds x20, x20, x30
	mul x1, x16, x13
	str x9, [sp, #64]
	adcs x23, x23, xzr
	umulh x9, x5, x8
	umulh x0, x16, x13
	mul x19, x7, x11
	str x9, [sp, #56]
	mul x9, x25, x13
	umulh x12, x7, x11
	mul x4, x4, x8
	str x9, [sp, #40]
	umulh x9, x25, x13
	umulh x15, x15, x8
	mul x18, x25, x11
	str x9, [sp, #32]
	mul x9, x25, x14
	mul x24, x17, x8
	stur x9, [x29, #-88]
	umulh x9, x25, x14
	mul x25, x3, x11
	mul x11, x16, x11
	str x9, [sp, #96]
	mrs x30, NZCV
	ldr x6, [sp, #64]
	mul x9, x2, x8
	adds x20, x20, x6
	ldr x6, [sp, #56]
	adcs x23, x23, x6
	ldur x6, [x29, #-64]
	cset w21, hs
	adds x26, x26, x6
	umulh x6, x17, x8
	adcs x28, x28, xzr
	mrs x22, NZCV
	msr NZCV, x30
	cinc x21, x21, hs
	umulh x30, x2, x13
	adds x9, x26, x9
	adcs x10, x28, x10
	mul x26, x2, x13
	cset w28, hs
	msr NZCV, x22
	cinc x22, x28, hs
	ldr x28, [sp, #40]
	adds x9, x9, x28
	ldr x28, [sp, #32]
	adcs x28, x10, x28
	mul x10, x3, x8
	cinc x22, x22, hs
	adds x1, x20, x1
	adcs x0, x23, x0
	stur x1, [x29, #-64]
	umulh x1, x3, x8
	cinc x21, x21, hs
	adds x0, x0, x19
	adcs x12, x21, x12
	mul x21, x17, x13
	mrs x19, NZCV
	adds x0, x0, x4
	umulh x4, x17, x13
	adcs x20, x12, x15
	ldr x12, [sp, #48]
	cset w15, hs
	adds x23, x28, x25
	mul x28, x2, x14
	adcs x22, x22, x12
	mrs x25, NZCV
	umulh x2, x2, x14
	msr NZCV, x19
	cinc x15, x15, hs
	adds x19, x23, x24
	mul x23, x3, x13
	adcs x6, x22, x6
	cset w22, hs
	msr NZCV, x25
	ldur x12, [x29, #-88]
	cinc x22, x22, hs
	adds x19, x19, x26
	umulh x24, x3, x13
	adcs x6, x6, x30
	mov x30, #-4294967295
	cset w25, hs
	adds x19, x19, x12
	ldr x12, [sp, #96]
	mul x26, x17, x14
	adcs x6, x6, x12
	adc x22, x22, x25
	adds x10, x6, x10
	umulh x17, x17, x14
	adcs x1, x22, x1
	mov x22, #-4294967295
	mov x25, x7
	cset w6, hs
	adds x10, x10, x21
	adcs x1, x1, x4
	umulh x4, x3, x14
	cset w21, hs
	adds x10, x10, x28
	ldp x12, x28, [x29, #-80]
	adcs x1, x1, x2
	mul x2, x3, x14
	adc x3, x21, x6
	adds x1, x1, x23
	mov x6, #-4294967295
	adcs x3, x3, x24
	umulh x6, x18, x6
	cinc x4, x4, hs
	adds x1, x1, x26
	adcs x17, x3, x17
	lsl x3, x18, #32
	mul x24, x16, x14
	cset w21, hs
	adds x17, x17, x2
	lsr x2, x18, #32
	adc x4, x4, x21
	adds x21, x12, x3
	sub x18, x18, x3
	adcs x9, x9, x2
	lsl x2, x21, #32
	umulh x3, x21, x22
	adcs x18, x19, x18
	adcs x10, x10, x6
	lsr x6, x21, #32
	mul x12, x5, x13
	adcs x1, x1, xzr
	adcs x17, x17, xzr
	umulh x26, x28, x13
	cset w19, hs
	adds x22, x9, x2
	sub x9, x21, x2
	adcs x18, x18, x6
	mov x6, #-4294967295
	lsl x2, x22, #32
	adcs x10, x10, x9
	lsr x23, x22, #32
	umulh x21, x5, x13
	adcs x1, x1, x3
	adcs x17, x17, xzr
	adcs x3, x4, x19
	umulh x4, x22, x6
	cset w9, hs
	adds x18, x18, x2
	sub x2, x22, x2
	adcs x10, x10, x23
	umulh x6, x18, x6
	adcs x1, x1, x2
	umulh x22, x16, x14
	adcs x16, x17, x4
	lsl x17, x18, #32
	adcs x2, x3, xzr
	lsr x3, x18, #32
	mul x19, x7, x8
	cset w23, hs
	adds x10, x10, x17
	sub x17, x18, x17
	adcs x18, x1, x3
	mov w3, #-1
	umulh x7, x7, x8
	adcs x16, x16, x17
	eor w9, w9, w23
	adcs x1, x2, x6
	mul x4, x28, x13
	cset w17, hs
	adds x2, x10, #1
	sbcs x8, x18, x3
	eor w9, w9, w17
	mul x6, x5, x14
	cset w3, lo
	eor w9, w9, #0x1
	subs x23, x16, x3
	umulh x5, x5, x14
	sbcs x30, x1, x30
	cset w17, lo
	tst w17, w9
	mul x9, x25, x13
	csel x3, x10, x2, ne
	csel x2, x18, x8, ne
	csel x17, x16, x23, ne
	csel x16, x1, x30, ne
	adds x10, x0, x12
	umulh x13, x25, x13
	adcs x12, x20, x21
	cset w18, hs
	adds x10, x10, x24
	mul x0, x28, x14
	adcs x12, x12, x22
	adc x15, x15, x18
	adds x12, x12, x19
	umulh x1, x25, x14
	adcs x8, x15, x7
	cset w18, hs
	adds x12, x12, x4
	umulh x15, x28, x14
	adcs x8, x8, x26
	cset w4, hs
	adds x12, x12, x6
	mul x14, x25, x14
	adcs x8, x8, x5
	mov x5, #-4294967295
	adc x18, x4, x18
	adds x8, x8, x9
	adcs x9, x18, x13
	cinc x13, x1, hs
	adds x8, x8, x0
	mov x1, #-4294967295
	adcs x9, x9, x15
	lsl x15, x11, #32
	umulh x18, x11, x1
	cset w0, hs
	adds x9, x9, x14
	lsr x14, x11, #32
	adc x13, x13, x0
	ldur x0, [x29, #-56]
	sub x11, x11, x15
	adds x0, x0, x15
	ldur x15, [x29, #-64]
	adcs x14, x15, x14
	umulh x15, x0, x1
	adcs x10, x10, x11
	adcs x11, x12, x18
	lsl x12, x0, #32
	lsr x18, x0, #32
	adcs x8, x8, xzr
	adcs x9, x9, xzr
	cset w1, hs
	adds x14, x14, x12
	sub x12, x0, x12
	adcs x10, x10, x18
	lsr x18, x14, #32
	adcs x11, x11, x12
	lsl x12, x14, #32
	adcs x8, x8, x15
	umulh x15, x14, x5
	adcs x9, x9, xzr
	adcs x13, x13, x1
	mov w1, #-1
	cset w0, hs
	adds x10, x10, x12
	sub x12, x14, x12
	adcs x11, x11, x18
	umulh x14, x10, x5
	adcs x8, x8, x12
	lsl x12, x10, #32
	adcs x9, x9, x15
	lsr x15, x10, #32
	adcs x13, x13, xzr
	sub x10, x10, x12
	cset w18, hs
	adds x11, x11, x12
	adcs x8, x8, x15
	eor w18, w0, w18
	adcs x9, x9, x10
	adcs x10, x13, x14
	cset w12, hs
	adds x13, x11, #1
	sbcs x14, x8, x1
	eor w12, w18, w12
	cset w15, lo
	eor w12, w12, #0x1
	subs x15, x9, x15
	sbcs x0, x10, x5
	cset w18, lo
	tst w18, w12
	ldr x12, [x27, #160]
	csel x11, x11, x13, ne
	ldur x13, [x29, #-24]
	csel x8, x8, x14, ne
	csel x9, x9, x15, ne
	stp x11, x8, [x13]
	ldr x8, [sp, #8]
	add x12, x8, x12, lsl #6
	csel x8, x10, x0, ne
	stp x9, x8, [x13, #16]
	ldp x10, x11, [x12, #32]
	ldp x8, x9, [x27, #32]
	ldp x14, x15, [x12, #48]
	stp x9, x8, [x29, #-64]
	subs x13, x8, x10
	sbcs x18, x9, x11
	ldp x8, x9, [x27, #48]
	sbcs x14, x8, x14
	sbcs x15, x9, x15
	stp x9, x8, [x29, #-80]
	mov x9, #-4294967295
	cset w0, lo
	csel x4, x1, xzr, lo
	csel x1, x5, xzr, lo
	sbfx x0, x0, #0, #1
	mov w8, #-1
	adds x5, x13, x0
	umulh x13, x3, x5
	adcs x18, x18, x4
	mrs x4, NZCV
	mul x0, x2, x5
	umulh x6, x2, x5
	mul x7, x3, x18
	adds x13, x13, x0
	umulh x19, x3, x18
	cinc x0, x6, hs
	mul x20, x17, x5
	adds x13, x13, x7
	umulh x21, x17, x5
	adcs x0, x0, x19
	mul x22, x2, x18
	umulh x6, x2, x18
	cinc x7, x21, hs
	adds x0, x0, x20
	adcs x7, x7, xzr
	mul x21, x16, x5
	mrs x19, NZCV
	adds x20, x0, x22
	umulh x22, x16, x5
	adcs x6, x7, x6
	cset w0, hs
	msr NZCV, x19
	mul x23, x17, x18
	cinc x7, x0, hs
	msr NZCV, x4
	adcs x0, x14, xzr
	umulh x24, x17, x18
	adc x1, x1, x15
	mul x4, x3, x0
	umulh x19, x3, x0
	mul x14, x3, x5
	adds x15, x20, x4
	mul x5, x2, x0
	adcs x6, x6, x19
	cinc x7, x7, hs
	adds x6, x6, x21
	umulh x25, x2, x0
	adcs x7, x7, x22
	mrs x20, NZCV
	mul x26, x3, x1
	adds x6, x6, x23
	umulh x3, x3, x1
	adcs x7, x7, x24
	cset w21, hs
	msr NZCV, x20
	mul x28, x16, x18
	cinc x20, x21, hs
	adds x5, x6, x5
	adcs x6, x7, x25
	umulh x18, x16, x18
	cset w7, hs
	adds x5, x5, x26
	adcs x3, x6, x3
	mul x30, x17, x0
	adc x7, x20, x7
	adds x3, x3, x28
	umulh x4, x17, x0
	adcs x18, x7, x18
	mul x19, x2, x1
	cset w7, hs
	adds x3, x3, x30
	umulh x2, x2, x1
	adcs x18, x18, x4
	mul x22, x16, x0
	cset w20, hs
	adds x3, x3, x19
	umulh x0, x16, x0
	adcs x18, x18, x2
	mul x6, x17, x1
	umulh x4, x16, x1
	umulh x17, x17, x1
	mul x16, x16, x1
	adc x1, x20, x7
	adds x18, x18, x22
	adcs x0, x1, x0
	cinc x1, x4, hs
	adds x18, x18, x6
	umulh x2, x14, x9
	adcs x17, x0, x17
	lsl x0, x14, #32
	cset w4, hs
	adds x16, x17, x16
	lsr x17, x14, #32
	sub x14, x14, x0
	adc x1, x1, x4
	adds x13, x13, x0
	adcs x15, x15, x17
	lsl x0, x13, #32
	adcs x14, x5, x14
	adcs x17, x3, x2
	umulh x2, x13, x9
	lsr x3, x13, #32
	adcs x18, x18, xzr
	sub x13, x13, x0
	adcs x16, x16, xzr
	cset w4, hs
	adds x15, x15, x0
	adcs x14, x14, x3
	umulh x0, x15, x9
	adcs x13, x17, x13
	adcs x17, x18, x2
	lsl x18, x15, #32
	lsr x2, x15, #32
	adcs x16, x16, xzr
	adcs x1, x1, x4
	sub x15, x15, x18
	cset w3, hs
	adds x14, x14, x18
	adcs x13, x13, x2
	umulh x18, x14, x9
	adcs x15, x17, x15
	lsl x17, x14, #32
	adcs x16, x16, x0
	adcs x0, x1, xzr
	lsr x1, x14, #32
	sub x14, x14, x17
	cset w2, hs
	adds x13, x13, x17
	adcs x15, x15, x1
	eor w2, w3, w2
	adcs x14, x16, x14
	adcs x16, x0, x18
	cset w17, hs
	adds x18, x13, #1
	sbcs x0, x15, x8
	eor w17, w2, w17
	cset w1, lo
	eor w17, w17, #0x1
	subs x1, x14, x1
	sbcs x3, x16, x9
	cset w2, lo
	tst w2, w17
	csel x18, x13, x18, ne
	csel x15, x15, x0, ne
	csel x14, x14, x1, ne
	umulh x17, x18, x18
	csel x13, x16, x3, ne
	mul x0, x18, x15
	umulh x2, x18, x15
	mul x1, x18, x14
	adds x16, x17, x0
	umulh x4, x18, x14
	cinc x17, x2, hs
	adds x16, x16, x0
	mul x5, x15, x15
	adcs x17, x17, x2
	umulh x6, x15, x15
	cinc x2, x4, hs
	adds x17, x17, x1
	adcs x2, x2, xzr
	mul x3, x18, x13
	mrs x7, NZCV
	adds x17, x17, x5
	umulh x0, x18, x13
	adcs x5, x2, x6
	cset w2, hs
	msr NZCV, x7
	mul x19, x15, x14
	cinc x7, x2, hs
	adds x17, x17, x1
	umulh x6, x15, x14
	adcs x1, x5, x4
	cinc x4, x7, hs
	adds x1, x1, x3
	adcs x4, x4, x0
	mul x5, x15, x13
	mrs x7, NZCV
	adds x1, x1, x19
	umulh x20, x15, x13
	adcs x4, x4, x6
	cset w21, hs
	msr NZCV, x7
	mul x22, x14, x14
	cinc x7, x21, hs
	adds x1, x1, x19
	adcs x4, x4, x6
	umulh x6, x14, x14
	cset w19, hs
	adds x1, x1, x3
	adcs x0, x4, x0
	mul x3, x14, x13
	adc x4, x7, x19
	adds x0, x0, x5
	adcs x4, x4, x20
	umulh x7, x14, x13
	cset w19, hs
	adds x0, x0, x22
	adcs x4, x4, x6
	umulh x6, x13, x13
	cset w21, hs
	adds x0, x0, x5
	mul x2, x18, x18
	adcs x4, x4, x20
	adc x19, x21, x19
	adds x4, x4, x3
	mul x5, x13, x13
	adcs x19, x19, x7
	cinc x6, x6, hs
	adds x3, x4, x3
	adcs x4, x19, x7
	lsl x7, x2, #32
	umulh x19, x2, x9
	cset w20, hs
	adds x4, x4, x5
	lsr x5, x2, #32
	sub x2, x2, x7
	adc x6, x6, x20
	adds x16, x16, x7
	adcs x17, x17, x5
	umulh x5, x16, x9
	lsr x7, x16, #32
	adcs x1, x1, x2
	adcs x0, x0, x19
	adcs x2, x3, xzr
	lsl x3, x16, #32
	adcs x4, x4, xzr
	cset w19, hs
	adds x17, x17, x3
	sub x16, x16, x3
	adcs x1, x1, x7
	lsl x3, x17, #32
	adcs x16, x0, x16
	adcs x0, x2, x5
	adcs x2, x4, xzr
	umulh x4, x17, x9
	adcs x5, x6, x19
	lsr x6, x17, #32
	sub x17, x17, x3
	cset w7, hs
	adds x1, x1, x3
	adcs x16, x16, x6
	umulh x3, x1, x9
	adcs x17, x0, x17
	ldp x19, x20, [x27]
	adcs x0, x2, x4
	lsl x2, x1, #32
	adcs x4, x5, xzr
	lsr x5, x1, #32
	cset w6, hs
	adds x16, x16, x2
	sub x1, x1, x2
	adcs x17, x17, x5
	eor w6, w7, w6
	adcs x0, x0, x1
	adcs x1, x4, x3
	cset w2, hs
	adds x3, x16, #1
	sbcs x4, x17, x8
	eor w2, w6, w2
	mov w8, #-1
	cset w5, lo
	eor w2, w2, #0x1
	subs x5, x0, x5
	sbcs x6, x1, x9
	cset w7, lo
	tst w7, w2
	csel x16, x16, x3, ne
	csel x17, x17, x4, ne
	csel x0, x0, x5, ne
	ldp x3, x4, [x27, #16]
	csel x1, x1, x6, ne
	subs x16, x16, x19
	sbcs x17, x17, x20
	sbcs x0, x0, x3
	sbcs x1, x1, x4
	cset w2, lo
	csel x5, x8, xzr, lo
	csel x6, x9, xzr, lo
	sbfx x2, x2, #0, #1
	adds x16, x16, x2
	adcs x17, x17, x5
	ldp x2, x5, [x12]
	adcs x0, x0, xzr
	adc x1, x6, x1
	subs x16, x16, x2
	ldp x2, x12, [x12, #16]
	sbcs x17, x17, x5
	sbcs x0, x0, x2
	sbcs x1, x1, x12
	cset w12, lo
	csel x2, x8, xzr, lo
	csel x5, x9, xzr, lo
	sbfx x12, x12, #0, #1
	adds x12, x16, x12
	adcs x2, x17, x2
	adcs x17, x0, xzr
	adc x16, x5, x1
	subs x5, x19, x12
	sbcs x21, x20, x2
	sbcs x1, x3, x17
	sbcs x0, x4, x16
	cset w3, lo
	csel x22, x8, xzr, lo
	csel x6, x9, xzr, lo
	sbfx x3, x3, #0, #1
	adds x5, x5, x3
	adcs x21, x21, x22
	umulh x20, x5, x18
	adcs x1, x1, xzr
	mul x22, x21, x18
	adc x6, x6, x0
	mul x23, x5, x15
	umulh x25, x21, x18
	adds x20, x20, x22
	umulh x24, x5, x15
	mul x30, x1, x18
	cinc x25, x25, hs
	adds x20, x20, x23
	umulh x8, x1, x18
	adcs x23, x25, x24
	mul x26, x21, x15
	umulh x28, x21, x15
	cinc x8, x8, hs
	adds x23, x23, x30
	mul x19, x5, x14
	adcs x8, x8, xzr
	mrs x25, NZCV
	umulh x7, x5, x14
	mul x9, x1, x15
	umulh x10, x1, x15
	mul x0, x6, x18
	mul x24, x6, x15
	umulh x30, x6, x15
	adds x15, x23, x26
	adcs x8, x8, x28
	mul x3, x5, x18
	cset w26, hs
	msr NZCV, x25
	cinc x25, x26, hs
	adds x15, x15, x19
	umulh x18, x6, x18
	adcs x8, x8, x7
	cinc x7, x25, hs
	adds x8, x8, x0
	mul x11, x21, x14
	umulh x22, x21, x14
	adcs x18, x7, x18
	mrs x7, NZCV
	mul x4, x5, x13
	adds x8, x8, x9
	umulh x5, x5, x13
	adcs x10, x18, x10
	cset w18, hs
	msr NZCV, x7
	stp x12, x2, [x27]
	cinc x18, x18, hs
	adds x8, x8, x11
	mul x23, x1, x14
	adcs x10, x10, x22
	umulh x19, x1, x14
	stp x17, x16, [x27, #16]
	mul x0, x21, x13
	umulh x9, x21, x13
	mul x21, x6, x14
	umulh x11, x6, x14
	cset w14, hs
	adds x8, x8, x4
	adcs x10, x10, x5
	adc x14, x18, x14
	adds x10, x10, x24
	mul x4, x1, x13
	adcs x14, x14, x30
	umulh x18, x1, x13
	cset w1, hs
	adds x10, x10, x23
	adcs x14, x14, x19
	umulh x5, x6, x13
	cset w7, hs
	adds x10, x10, x0
	adcs x9, x14, x9
	mul x13, x6, x13
	adc x14, x7, x1
	adds x9, x9, x21
	adcs x11, x14, x11
	cinc x14, x5, hs
	adds x9, x9, x4
	mov x4, #-4294967295
	adcs x11, x11, x18
	lsl x18, x3, #32
	umulh x0, x3, x4
	cset w1, hs
	adds x11, x11, x13
	lsr x13, x3, #32
	adc x14, x14, x1
	adds x1, x20, x18
	sub x18, x3, x18
	adcs x13, x15, x13
	lsl x15, x1, #32
	adcs x8, x8, x18
	umulh x18, x1, x4
	adcs x10, x10, x0
	lsr x0, x1, #32
	adcs x9, x9, xzr
	adcs x11, x11, xzr
	cset w3, hs
	adds x13, x13, x15
	sub x15, x1, x15
	adcs x8, x8, x0
	lsr x0, x13, #32
	adcs x10, x10, x15
	lsl x15, x13, #32
	adcs x9, x9, x18
	umulh x18, x13, x4
	adcs x11, x11, xzr
	sub x13, x13, x15
	adcs x14, x14, x3
	mov w3, #-1
	cset w1, hs
	adds x8, x8, x15
	adcs x10, x10, x0
	umulh x15, x8, x4
	adcs x9, x9, x13
	lsl x13, x8, #32
	adcs x11, x11, x18
	lsr x18, x8, #32
	adcs x14, x14, xzr
	sub x8, x8, x13
	cset w0, hs
	adds x10, x10, x13
	adcs x9, x9, x18
	eor w0, w1, w0
	adcs x8, x11, x8
	adcs x11, x14, x15
	cset w13, hs
	adds x14, x10, #1
	sbcs x15, x9, x3
	eor w13, w0, w13
	cset w18, lo
	eor w13, w13, #0x1
	subs x18, x8, x18
	sbcs x1, x11, x4
	cset w0, lo
	tst w0, w13
	csel x10, x10, x14, ne
	ldp x13, x14, [x29, #-64]
	csel x9, x9, x15, ne
	csel x8, x8, x18, ne
	csel x11, x11, x1, ne
	subs x10, x10, x14
	sbcs x9, x9, x13
	ldp x13, x14, [x29, #-80]
	sbcs x8, x8, x14
	sbcs x11, x11, x13
	cset w13, lo
	csel x14, x3, xzr, lo
	csel x15, x4, xzr, lo
	sbfx x13, x13, #0, #1
	adds x10, x10, x13
	adcs x9, x9, x14
	stp x10, x9, [x27, #32]
	ldp x13, x10, [x27, #136]
	adcs x14, x8, xzr
	ldr x9, [sp, #24]
	adc x11, x15, x11
	adds x8, x13, #1
	ldr x13, [x27, #152]
	stp x14, x11, [x27, #48]
	cinc x10, x10, hs
	cmp x10, x9
	ldp x9, x14, [x27, #88]
	ccmp x2, #0, #4, lo
	and x13, x13, #0x1
	stp x8, x10, [x27, #136]
	csinc x13, x13, xzr, ne
	cmp x12, x9
	ldp x9, x11, [x27, #104]
	ccmp x17, x14, #0, eq
	ccmp x16, x9, #0, eq
	csinc x10, x13, xzr, ne
	cmp x12, x11
	ldp x9, x11, [x27, #120]
	ccmp x17, x9, #0, eq
	ccmp x16, x11, #0, eq
	csinc x9, x10, xzr, ne
	cbnz x9, .LBB5_7
	tst x8, #0x3f
	b.ne .LBB5_5
	stp x12, x17, [x27, #112]
	str x16, [x27, #128]
.LBB5_5:
	ldp x9, x10, [x27, #72]
	mov x8, xzr
	stp x9, x10, [x27, #96]
	ldr x9, [x27, #64]
	stp x12, x17, [x27, #64]
	stp x16, x9, [x27, #80]
.LBB5_6:
	ldp x10, x11, [x29, #-40]
	sub x27, x27, #168
	ldur x12, [x29, #-48]
	ldur x9, [x29, #-16]
	add x11, x11, #1
	sub x12, x12, #32
	add x9, x8, x9
	cmp x11, #1
	sub x10, x10, #32
	stur x9, [x29, #-16]
	b.ne .LBB5_2
	b .LBB5_9
.LBB5_7:
	mov x0, x27
	bl slow_path
	mov w8, #1
	b .LBB5_6
.LBB5_8:
	stur xzr, [x29, #-16]
.LBB5_9:
	ldur x0, [x29, #-16]
	ldp x20, x19, [sp, #272]
	ldp x22, x21, [sp, #256]
	ldp x24, x23, [sp, #240]
	ldp x26, x25, [sp, #224]
	ldp x28, x27, [sp, #208]
	ldp x29, x30, [sp, #192]
	add sp, sp, #288
	ret
.Lfunc_end5:
pass1_dbl:
	cbz x4, .LBB6_4
	str x30, [sp, #-96]!
	mov x8, x0
	mov x0, xzr
	add x9, x1, #16
	stp x28, x27, [sp, #16]
	stp x26, x25, [sp, #32]
	stp x24, x23, [sp, #48]
	stp x22, x21, [sp, #64]
	stp x20, x19, [sp, #80]
.LBB6_2:
	ldr x12, [x8], #8
	str x8, [sp, #8]
	mov w8, #-1
	ldp x13, x14, [x12, #32]
	ldp x15, x12, [x12, #48]
	adds x13, x13, x13
	adcs x14, x14, x14
	adcs x16, x15, x15
	adcs x12, x12, x12
	cset w15, hs
	adds x17, x13, #1
	sbcs x18, x14, x8
	mov x8, #-4294967295
	eor w15, w15, #0x1
	cset w1, lo
	subs x1, x16, x1
	sbcs x5, x12, x8
	cset w6, lo
	tst w6, w15
	csel x15, x13, x17, ne
	csel x14, x14, x18, ne
	csel x13, x16, x1, ne
	csel x12, x12, x5, ne
	stp x15, x14, [x9, #-16]
	stp x13, x12, [x9], #32
	ldp x16, x5, [x3]
	ldp x17, x1, [x3, #16]
	umulh x19, x15, x16
	mul x20, x14, x16
	umulh x18, x14, x16
	mul x21, x15, x5
	adds x19, x19, x20
	umulh x22, x15, x5
	cinc x18, x18, hs
	mul x23, x13, x16
	adds x19, x19, x21
	umulh x24, x13, x16
	adcs x18, x18, x22
	mul x6, x14, x5
	umulh x7, x14, x5
	cinc x22, x24, hs
	adds x18, x18, x23
	mul x25, x15, x17
	adcs x22, x22, xzr
	mrs x23, NZCV
	umulh x26, x15, x17
	adds x18, x18, x6
	adcs x6, x22, x7
	mul x27, x12, x16
	umulh x28, x12, x16
	mul x30, x13, x5
	umulh x10, x13, x5
	mul x24, x12, x5
	umulh x7, x12, x5
	cset w5, hs
	msr NZCV, x23
	cinc x22, x5, hs
	adds x5, x18, x25
	adcs x6, x6, x26
	mul x11, x14, x17
	cinc x22, x22, hs
	adds x6, x6, x27
	adcs x22, x22, x28
	umulh x8, x14, x17
	mrs x25, NZCV
	adds x6, x6, x30
	mul x16, x15, x16
	adcs x10, x22, x10
	orr x22, x13, x12
	mul x20, x15, x1
	umulh x21, x15, x1
	orr x15, x15, x14
	orr x15, x15, x22
	cset w22, hs
	msr NZCV, x25
	cinc x22, x22, hs
	cmp x15, #0
	mul x18, x13, x17
	cinc x0, x0, eq
	adds x11, x6, x11
	adcs x8, x10, x8
	umulh x23, x13, x17
	mul x15, x12, x17
	umulh x10, x12, x17
	cset w17, hs
	adds x11, x11, x20
	adcs x8, x8, x21
	mul x26, x14, x1
	adc x17, x22, x17
	adds x8, x8, x24
	adcs x17, x17, x7
	umulh x14, x14, x1
	cset w7, hs
	adds x8, x8, x18
	adcs x17, x17, x23
	mul x6, x13, x1
	cset w20, hs
	adds x8, x8, x26
	umulh x18, x12, x1
	adcs x14, x17, x14
	umulh x13, x13, x1
	adc x17, x20, x7
	adds x14, x14, x15
	adcs x10, x17, x10
	mul x12, x12, x1
	mov x1, #-4294967295
	cinc x15, x18, hs
	adds x14, x14, x6
	umulh x17, x16, x1
	adcs x10, x10, x13
	lsl x13, x16, #32
	cset w18, hs
	adds x10, x10, x12
	lsr x12, x16, #32
	adc x15, x15, x18
	adds x18, x19, x13
	sub x13, x16, x13
	adcs x12, x5, x12
	umulh x16, x18, x1
	mov x5, #-4294967295
	adcs x11, x11, x13
	adcs x8, x8, x17
	lsr x17, x18, #32
	adcs x13, x14, xzr
	lsl x14, x18, #32
	adcs x10, x10, xzr
	cset w1, hs
	adds x12, x12, x14
	sub x14, x18, x14
	adcs x11, x11, x17
	lsr x17, x12, #32
	adcs x8, x8, x14
	lsl x14, x12, #32
	adcs x13, x13, x16
	umulh x16, x12, x5
	adcs x10, x10, xzr
	sub x12, x12, x14
	adcs x15, x15, x1
	cset w18, hs
	adds x11, x11, x14
	adcs x8, x8, x17
	umulh x14, x11, x5
	adcs x12, x13, x12
	lsl x13, x11, #32
	adcs x10, x10, x16
	lsr x16, x11, #32
	adcs x15, x15, xzr
	sub x11, x11, x13
	cset w17, hs
	adds x8, x8, x13
	adcs x12, x12, x16
	eor w17, w18, w17
	adcs x10, x10, x11
	adcs x11, x15, x14
	mov w15, #-1
	cset w13, hs
	adds x14, x8, #1
	sbcs x15, x12, x15
	eor w13, w17, w13
	cset w16, lo
	eor w13, w13, #0x1
	subs x16, x10, x16
	sbcs x17, x11, x5
	cset w18, lo
	tst w18, w13
	csel x8, x8, x14, ne
	csel x12, x12, x15, ne
	csel x10, x10, x16, ne
	csel x11, x11, x17, ne
	stp x8, x12, [x3]
	subs x4, x4, #1
	stp x10, x11, [x3, #16]
	stp x10, x11, [x2, #16]
	stp x8, x12, [x2], #32
	ldr x8, [sp, #8]
	b.ne .LBB6_2
	ldp x20, x19, [sp, #80]
	ldp x22, x21, [sp, #64]
	ldp x24, x23, [sp, #48]
	ldp x26, x25, [sp, #32]
	ldp x28, x27, [sp, #16]
	ldr x30, [sp], #96
	ret
.LBB6_4:
	mov x0, xzr
	ret
.Lfunc_end6:
pass2_dbl:
	sub sp, sp, #336
	stp x29, x30, [sp, #240]
	add x29, sp, #240
	stp x28, x27, [sp, #256]
	stp x26, x25, [sp, #272]
	stp x24, x23, [sp, #288]
	stp x22, x21, [sp, #304]
	stp x20, x19, [sp, #320]
	stp x3, x7, [sp, #24]
	stp x5, x4, [sp, #8]
	stur x6, [x29, #-16]
	cbz x6, .LBB7_10
	ldur x9, [x29, #-16]
	sub x8, x0, #8
	stur xzr, [x29, #-24]
	str x8, [sp]
	lsl x9, x9, #5
	add x10, x9, x2
	add x8, x9, x1
	sub x13, x10, #64
	sub x11, x8, #16
.LBB7_2:
	ldr x12, [sp, #8]
	ldp x9, x8, [x11, #-16]
	stur x11, [x29, #-40]
	ldp x10, x11, [x11]
	ldp x16, x15, [x12]
	ldp x28, x25, [x12, #16]
	umulh x18, x9, x16
	mul x0, x8, x16
	umulh x1, x8, x16
	mul x3, x9, x15
	adds x18, x18, x0
	umulh x2, x9, x15
	cinc x0, x1, hs
	mul x4, x10, x16
	adds x18, x18, x3
	umulh x5, x10, x16
	mul x6, x8, x15
	umulh x7, x8, x15
	mul x19, x9, x28
	umulh x20, x9, x28
	mul x17, x9, x16
	mul x27, x9, x25
	umulh x1, x9, x25
	adcs x9, x0, x2
	cinc x0, x5, hs
	adds x9, x9, x4
	adcs x0, x0, xzr
	mul x21, x11, x16
	mrs x3, NZCV
	adds x9, x9, x6
	umulh x22, x11, x16
	adcs x5, x0, x7
	cset w0, hs
	msr NZCV, x3
	mul x23, x10, x15
	cinc x3, x0, hs
	adds x0, x9, x19
	umulh x24, x10, x15
	adcs x9, x5, x20
	cinc x3, x3, hs
	adds x9, x9, x21
	adcs x3, x3, x22
	mul x14, x8, x28
	mrs x7, NZCV
	adds x20, x9, x23
	umulh x26, x8, x28
	adcs x3, x3, x24
	mul x19, x8, x25
	cset w9, hs
	umulh x21, x8, x25
	ldur x8, [x29, #-16]
	msr NZCV, x7
	cinc x7, x9, hs
	stur x13, [x29, #-32]
	subs x9, x8, #1
	mul x2, x11, x15
	stur x9, [x29, #-16]
	ldr x9, [sp, #24]
	umulh x4, x11, x15
	csel x9, x9, x13, eq
	adds x20, x20, x14
	mov x13, #-4294967295
	mul x6, x10, x28
	adcs x3, x3, x26
	mov w14, #-1
	cset w24, hs
	adds x20, x20, x27
	umulh x5, x10, x28
	adcs x1, x3, x1
	adc x7, x7, x24
	adds x1, x1, x2
	mul x22, x11, x28
	adcs x2, x7, x4
	cset w4, hs
	adds x1, x1, x6
	umulh x23, x11, x28
	ldp x26, x24, [x9, #16]
	adcs x2, x2, x5
	mul x3, x10, x25
	cset w6, hs
	adds x1, x1, x19
	adcs x2, x2, x21
	mov x19, #-4294967295
	umulh x5, x11, x25
	adc x4, x6, x4
	adds x2, x2, x22
	adcs x4, x4, x23
	umulh x10, x10, x25
	mul x11, x11, x25
	cinc x5, x5, hs
	adds x2, x2, x3
	lsl x3, x17, #32
	adcs x10, x4, x10
	umulh x4, x17, x13
	cset w6, hs
	adds x10, x10, x11
	lsr x11, x17, #32
	sub x17, x17, x3
	adc x5, x5, x6
	adds x18, x18, x3
	adcs x11, x0, x11
	umulh x3, x18, x19
	adcs x17, x20, x17
	adcs x0, x1, x4
	lsr x4, x18, #32
	adcs x1, x2, xzr
	lsl x2, x18, #32
	adcs x10, x10, xzr
	cset w6, hs
	adds x11, x11, x2
	sub x18, x18, x2
	adcs x17, x17, x4
	umulh x2, x11, x19
	lsr x4, x11, #32
	adcs x18, x0, x18
	adcs x0, x1, x3
	lsl x1, x11, #32
	adcs x10, x10, xzr
	adcs x3, x5, x6
	sub x11, x11, x1
	cset w5, hs
	adds x17, x17, x1
	adcs x18, x18, x4
	umulh x4, x17, x19
	lsr x1, x17, #32
	adcs x11, x0, x11
	lsl x0, x17, #32
	adcs x10, x10, x2
	adcs x3, x3, xzr
	sub x17, x17, x0
	cset w7, hs
	adds x2, x18, x0
	adcs x1, x11, x1
	eor w7, w5, w7
	adcs x13, x10, x17
	adcs x10, x3, x4
	cset w11, hs
	adds x6, x2, #1
	sbcs x3, x1, x14
	eor w7, w7, w11
	cset w17, lo
	subs x4, x13, x17
	ldp x17, x0, [x9]
	sbcs x5, x10, x19
	eor w19, w7, #0x1
	cset w11, lo
	umulh x14, x17, x16
	tst w11, w19
	csel x11, x2, x6, ne
	csel x9, x1, x3, ne
	csel x10, x10, x5, ne
	stp x11, x9, [x12]
	umulh x11, x0, x16
	csel x9, x13, x4, ne
	stp x9, x10, [x12, #16]
	ldr x9, [sp, #16]
	mov x6, x28
	mul x23, x17, x16
	stur x14, [x29, #-96]
	mul x14, x0, x16
	mul x21, x17, x15
	umulh x20, x17, x15
	stp x11, x14, [sp, #112]
	ldr x11, [sp]
	ldr x1, [x11, x8, lsl #3]
	mul x8, x26, x16
	umulh x11, x26, x16
	ldp x18, x13, [x1]
	stur x1, [x29, #-80]
	stp x18, x13, [x29, #-56]
	stp x11, x8, [sp, #96]
	mul x8, x24, x16
	ldp x10, x11, [x9]
	subs x19, x18, x10
	stur x8, [x29, #-104]
	mov x8, x13
	sbcs x7, x13, x11
	ldp x14, x13, [x1, #16]
	ldp x12, x1, [x9, #16]
	stp x13, x14, [x29, #-72]
	sbcs x9, x14, x12
	str x9, [sp, #72]
	umulh x9, x24, x16
	mov x16, #-4294967295
	stur x9, [x29, #-112]
	sbcs x9, x13, x1
	str x9, [sp, #88]
	mov w9, #-1
	cset w30, lo
	csel x27, x9, xzr, lo
	csel x22, x16, xzr, lo
	adds x5, x18, x10
	adcs x2, x8, x11
	mul x10, x17, x28
	adcs x12, x14, x12
	adcs x8, x13, x1
	cset w1, hs
	adds x14, x5, #1
	sbcs x13, x2, x9
	umulh x9, x17, x28
	eor w1, w1, #0x1
	cset w4, lo
	subs x11, x12, x4
	umulh x4, x26, x15
	stp x9, x10, [sp, #40]
	sbcs x10, x8, x16
	mul x16, x17, x25
	cset w9, lo
	tst w9, w1
	sbfx x9, x30, #0, #1
	csel x2, x2, x13, ne
	umulh x13, x17, x25
	csel x30, x5, x14, ne
	csel x1, x12, x11, ne
	csel x3, x8, x10, ne
	adds x19, x19, x9
	umulh x8, x30, x19
	adcs x27, x7, x27
	str x16, [sp, #80]
	mul x9, x2, x19
	stp x25, x13, [sp, #56]
	mrs x16, NZCV
	umulh x10, x2, x19
	mul x12, x27, x30
	adds x8, x8, x9
	umulh x13, x27, x30
	cinc x9, x10, hs
	mul x14, x1, x19
	adds x8, x8, x12
	umulh x7, x1, x19
	stur x8, [x29, #-88]
	adcs x9, x9, x13
	mul x11, x27, x2
	umulh x8, x27, x2
	cinc x10, x7, hs
	adds x9, x9, x14
	adcs x10, x10, xzr
	mul x13, x0, x15
	mrs x14, NZCV
	adds x18, x9, x11
	umulh x7, x0, x15
	adcs x17, x10, x8
	cset w8, hs
	msr NZCV, x14
	ldr x9, [sp, #120]
	cinc x12, x8, hs
	ldur x8, [x29, #-96]
	umulh x11, x3, x19
	adds x10, x8, x9
	ldr x8, [sp, #112]
	mul x5, x26, x15
	cinc x8, x8, hs
	adds x9, x10, x21
	str x9, [sp, #120]
	ldp x10, x9, [sp, #96]
	adcs x8, x8, x20
	cinc x20, x10, hs
	adds x8, x8, x9
	mul x9, x24, x15
	adcs x20, x20, xzr
	umulh x15, x24, x15
	str x9, [sp, #112]
	mrs x9, NZCV
	msr NZCV, x16
	ldr x10, [sp, #72]
	ldr x14, [sp, #88]
	adcs x28, x10, xzr
	mul x10, x3, x19
	adc x25, x22, x14
	adds x8, x8, x13
	adcs x7, x20, x7
	mul x16, x28, x30
	cset w20, hs
	msr NZCV, x9
	ldr x14, [sp, #48]
	umulh x21, x28, x30
	cinc x9, x20, hs
	adds x20, x8, x14
	ldr x8, [sp, #40]
	mul x13, x27, x1
	adcs x8, x7, x8
	cinc x9, x9, hs
	adds x16, x18, x16
	umulh x22, x27, x1
	stur x16, [x29, #-96]
	adcs x16, x17, x21
	cinc x12, x12, hs
	adds x10, x16, x10
	mul x14, x0, x6
	adcs x11, x12, x11
	mrs x16, NZCV
	ldr x21, [sp, #56]
	umulh x17, x0, x6
	adds x12, x10, x13
	mul x18, x0, x21
	adcs x7, x11, x22
	mov x10, x21
	cset w13, hs
	umulh x11, x0, x21
	ldur x0, [x29, #-104]
	adds x8, x8, x0
	ldur x0, [x29, #-112]
	mul x21, x26, x6
	adcs x9, x9, x0
	mrs x0, NZCV
	msr NZCV, x16
	cinc x16, x13, hs
	umulh x13, x26, x6
	adds x8, x8, x5
	adcs x9, x9, x4
	mul x5, x24, x6
	cset w4, hs
	msr NZCV, x0
	cinc x0, x4, hs
	adds x8, x8, x14
	ldr x4, [sp, #80]
	adcs x9, x9, x17
	umulh x14, x24, x6
	mov x6, #-4294967295
	cset w17, hs
	adds x8, x8, x4
	ldr x4, [sp, #64]
	adcs x9, x9, x4
	mul x4, x26, x10
	adc x17, x0, x17
	ldr x0, [sp, #112]
	adds x9, x9, x0
	adcs x15, x17, x15
	umulh x17, x26, x10
	cset w0, hs
	adds x9, x9, x21
	adcs x13, x15, x13
	umulh x15, x24, x10
	cset w21, hs
	adds x9, x9, x18
	adcs x11, x13, x11
	mul x13, x24, x10
	mov x10, #-4294967295
	adc x18, x21, x0
	adds x11, x11, x5
	adcs x14, x18, x14
	umulh x18, x23, x10
	cinc x15, x15, hs
	adds x11, x11, x4
	adcs x14, x14, x17
	lsl x17, x23, #32
	umulh x24, x25, x2
	cset w0, hs
	adds x13, x14, x13
	lsr x14, x23, #32
	adc x15, x15, x0
	ldr x0, [sp, #120]
	mul x21, x25, x30
	adds x0, x0, x17
	sub x17, x23, x17
	umulh x23, x28, x1
	adcs x14, x20, x14
	lsr x4, x0, #32
	adcs x8, x8, x17
	lsl x17, x0, #32
	adcs x9, x9, x18
	umulh x18, x0, x10
	adcs x11, x11, xzr
	adcs x13, x13, xzr
	mul x10, x28, x2
	cset w5, hs
	adds x14, x14, x17
	sub x17, x0, x17
	adcs x8, x8, x4
	lsl x0, x14, #32
	umulh x4, x14, x6
	adcs x17, x9, x17
	lsr x20, x14, #32
	adcs x18, x11, x18
	sub x14, x14, x0
	mul x9, x30, x19
	adcs x13, x13, xzr
	adcs x15, x15, x5
	umulh x19, x28, x2
	cset w11, hs
	adds x8, x8, x0
	adcs x17, x17, x20
	umulh x5, x8, x6
	adcs x14, x18, x14
	lsl x18, x8, #32
	adcs x13, x13, x4
	lsr x4, x8, #32
	umulh x0, x25, x30
	adcs x15, x15, xzr
	sub x8, x8, x18
	cset w22, hs
	adds x17, x17, x18
	mul x20, x27, x3
	adcs x14, x14, x4
	eor w11, w11, w22
	adcs x8, x13, x8
	umulh x4, x27, x3
	adcs x13, x15, x5
	mul x15, x25, x2
	cset w18, hs
	adds x5, x17, #1
	mov w2, #-1
	eor w11, w11, w18
	sbcs x26, x14, x2
	eor w11, w11, #0x1
	mul x27, x28, x1
	cset w2, lo
	subs x22, x8, x2
	sbcs x30, x13, x6
	cset w18, lo
	tst w18, w11
	mul x11, x28, x3
	csel x2, x17, x5, ne
	csel x18, x14, x26, ne
	csel x14, x8, x22, ne
	csel x13, x13, x30, ne
	adds x8, x12, x10
	umulh x12, x28, x3
	adcs x10, x7, x19
	mov x28, #-4294967295
	mov w30, #-1
	cset w17, hs
	adds x8, x8, x21
	adcs x10, x10, x0
	mul x0, x25, x1
	adc x16, x16, x17
	adds x10, x10, x20
	adcs x16, x16, x4
	umulh x17, x25, x1
	cset w1, hs
	adds x10, x10, x27
	adcs x16, x16, x23
	umulh x4, x25, x3
	cset w5, hs
	adds x10, x10, x15
	adcs x15, x16, x24
	mul x16, x25, x3
	mov x3, #-4294967295
	adc x1, x5, x1
	adds x11, x15, x11
	adcs x12, x1, x12
	cinc x15, x4, hs
	adds x11, x11, x0
	umulh x0, x9, x3
	adcs x12, x12, x17
	lsl x17, x9, #32
	cset w1, hs
	adds x12, x12, x16
	lsr x16, x9, #32
	adc x15, x15, x1
	ldur x1, [x29, #-88]
	sub x9, x9, x17
	adds x1, x1, x17
	ldur x17, [x29, #-96]
	adcs x16, x17, x16
	umulh x17, x1, x3
	adcs x8, x8, x9
	adcs x9, x10, x0
	lsr x0, x1, #32
	adcs x10, x11, xzr
	lsl x11, x1, #32
	adcs x12, x12, xzr
	cset w3, hs
	adds x16, x16, x11
	sub x11, x1, x11
	adcs x8, x8, x0
	lsr x0, x16, #32
	adcs x9, x9, x11
	adcs x10, x10, x17
	umulh x17, x16, x28
	adcs x11, x12, xzr
	lsl x12, x16, #32
	adcs x15, x15, x3
	mov w3, #-1
	cset w1, hs
	adds x8, x8, x12
	sub x12, x16, x12
	adcs x9, x9, x0
	umulh x16, x8, x28
	adcs x10, x10, x12
	lsl x12, x8, #32
	adcs x11, x11, x17
	lsr x17, x8, #32
	adcs x15, x15, xzr
	sub x8, x8, x12
	cset w0, hs
	adds x9, x9, x12
	adcs x10, x10, x17
	eor w0, w1, w0
	adcs x8, x11, x8
	adcs x11, x15, x16
	cset w12, hs
	adds x15, x9, #1
	sbcs x16, x10, x3
	eor w12, w0, w12
	cset w17, lo
	eor w12, w12, #0x1
	subs x17, x8, x17
	sbcs x1, x11, x28
	cset w0, lo
	tst w0, w12
	csel x9, x9, x15, ne
	csel x10, x10, x16, ne
	csel x8, x8, x17, ne
	csel x11, x11, x1, ne
	adds x12, x9, x9
	adcs x15, x10, x10
	adcs x16, x8, x8
	adcs x17, x11, x11
	cset w0, hs
	adds x1, x12, #1
	sbcs x3, x15, x3
	eor w0, w0, #0x1
	cset w4, lo
	subs x4, x16, x4
	sbcs x5, x17, x28
	cset w7, lo
	tst w7, w0
	csel x12, x12, x1, ne
	csel x15, x15, x3, ne
	csel x16, x16, x4, ne
	csel x17, x17, x5, ne
	adds x9, x12, x9
	adcs x12, x15, x10
	adcs x8, x16, x8
	adcs x10, x17, x11
	cset w11, hs
	adds x15, x9, #1
	sbcs x16, x12, x30
	eor w11, w11, #0x1
	cset w17, lo
	subs x17, x8, x17
	sbcs x0, x10, x28
	cset w1, lo
	tst w1, w11
	csel x9, x9, x15, ne
	csel x15, x12, x16, ne
	csel x11, x8, x17, ne
	umulh x1, x9, x2
	csel x12, x10, x0, ne
	mul x3, x9, x18
	umulh x4, x9, x18
	mul x16, x15, x2
	adds x10, x1, x3
	umulh x5, x15, x2
	mul x7, x9, x14
	umulh x19, x9, x14
	mul x20, x15, x18
	mul x8, x11, x2
	umulh x17, x11, x2
	mul x22, x9, x13
	umulh x23, x9, x13
	mul x9, x9, x2
	mul x0, x12, x2
	umulh x1, x12, x2
	cinc x2, x4, hs
	adds x10, x10, x16
	adcs x16, x2, x5
	umulh x21, x15, x18
	cinc x3, x19, hs
	adds x16, x16, x7
	adcs x3, x3, xzr
	mrs x4, NZCV
	mul x24, x15, x14
	umulh x25, x15, x14
	mul x2, x15, x13
	umulh x5, x15, x13
	adds x15, x16, x20
	adcs x16, x3, x21
	cset w7, hs
	msr NZCV, x4
	mul x26, x11, x18
	cinc x4, x7, hs
	adds x15, x15, x8
	adcs x16, x16, x17
	umulh x27, x11, x18
	cinc x17, x4, hs
	adds x16, x16, x22
	adcs x17, x17, x23
	mul x3, x11, x14
	mrs x7, NZCV
	adds x16, x16, x24
	umulh x8, x11, x14
	adcs x17, x17, x25
	cset w19, hs
	msr NZCV, x7
	mul x4, x12, x18
	cinc x7, x19, hs
	adds x16, x16, x26
	adcs x17, x17, x27
	umulh x18, x12, x18
	cset w19, hs
	adds x16, x16, x0
	adcs x17, x17, x1
	mul x20, x11, x13
	adc x1, x7, x19
	adds x17, x17, x2
	mov w19, #-1
	adcs x1, x1, x5
	umulh x11, x11, x13
	cset w2, hs
	adds x17, x17, x3
	mul x0, x12, x14
	adcs x8, x1, x8
	cset w3, hs
	adds x17, x17, x4
	umulh x1, x12, x13
	adcs x8, x8, x18
	umulh x14, x12, x14
	mul x12, x12, x13
	adc x13, x3, x2
	adds x8, x8, x20
	adcs x11, x13, x11
	cinc x13, x1, hs
	adds x8, x8, x0
	umulh x18, x9, x28
	adcs x11, x11, x14
	lsl x14, x9, #32
	cset w0, hs
	adds x11, x11, x12
	lsr x12, x9, #32
	sub x9, x9, x14
	adc x13, x13, x0
	adds x10, x10, x14
	adcs x12, x15, x12
	lsl x15, x10, #32
	adcs x9, x16, x9
	umulh x16, x10, x28
	adcs x14, x17, x18
	lsr x17, x10, #32
	sub x10, x10, x15
	adcs x8, x8, xzr
	adcs x11, x11, xzr
	cset w18, hs
	adds x12, x12, x15
	adcs x9, x9, x17
	umulh x15, x12, x28
	adcs x10, x14, x10
	lsl x14, x12, #32
	adcs x8, x8, x16
	lsr x16, x12, #32
	adcs x11, x11, xzr
	sub x12, x12, x14
	adcs x13, x13, x18
	cset w17, hs
	adds x9, x9, x14
	adcs x10, x10, x16
	umulh x14, x9, x28
	adcs x8, x8, x12
	lsl x12, x9, #32
	adcs x11, x11, x15
	lsr x15, x9, #32
	adcs x13, x13, xzr
	sub x9, x9, x12
	cset w16, hs
	adds x10, x10, x12
	adcs x8, x8, x15
	eor w16, w17, w16
	adcs x9, x11, x9
	adcs x11, x13, x14
	cset w12, hs
	adds x13, x10, #1
	sbcs x14, x8, x30
	eor w12, w16, w12
	cset w15, lo
	eor w12, w12, #0x1
	subs x17, x9, x15
	sbcs x16, x11, x28
	cset w15, lo
	tst w15, w12
	csel x18, x10, x13, ne
	csel x15, x8, x14, ne
	csel x14, x9, x17, ne
	umulh x10, x18, x18
	csel x13, x11, x16, ne
	mul x8, x18, x15
	umulh x0, x18, x15
	mul x17, x18, x14
	adds x9, x10, x8
	umulh x1, x18, x14
	cinc x10, x0, hs
	adds x9, x9, x8
	mul x2, x15, x15
	adcs x8, x10, x0
	umulh x3, x15, x15
	cinc x10, x1, hs
	adds x8, x8, x17
	adcs x10, x10, xzr
	mul x16, x18, x13
	mrs x11, NZCV
	adds x8, x8, x2
	umulh x12, x18, x13
	adcs x3, x10, x3
	cset w4, hs
	msr NZCV, x11
	mul x0, x15, x14
	cinc x4, x4, hs
	adds x11, x8, x17
	umulh x2, x15, x14
	adcs x8, x3, x1
	cinc x1, x4, hs
	adds x8, x8, x16
	adcs x1, x1, x12
	mul x17, x15, x13
	mrs x3, NZCV
	adds x8, x8, x0
	umulh x4, x15, x13
	adcs x1, x1, x2
	cset w7, hs
	msr NZCV, x3
	mul x5, x14, x14
	cinc x3, x7, hs
	adds x8, x8, x0
	mov x7, #-4294967295
	adcs x1, x1, x2
	umulh x0, x14, x14
	cset w2, hs
	adds x8, x8, x16
	adcs x12, x1, x12
	mul x16, x14, x13
	adc x1, x3, x2
	adds x12, x12, x17
	adcs x1, x1, x4
	umulh x2, x14, x13
	cset w3, hs
	adds x12, x12, x5
	umulh x5, x13, x13
	adcs x0, x1, x0
	cset w1, hs
	adds x12, x12, x17
	mul x10, x18, x18
	adcs x0, x0, x4
	adc x1, x1, x3
	adds x0, x0, x16
	mul x17, x13, x13
	adcs x1, x1, x2
	cinc x3, x5, hs
	adds x16, x0, x16
	adcs x1, x1, x2
	lsl x0, x10, #32
	umulh x4, x10, x28
	lsr x2, x10, #32
	cset w5, hs
	adds x17, x1, x17
	sub x10, x10, x0
	adc x1, x3, x5
	adds x9, x9, x0
	adcs x11, x11, x2
	umulh x0, x9, x28
	lsr x2, x9, #32
	adcs x8, x8, x10
	adcs x10, x12, x4
	lsl x12, x9, #32
	adcs x16, x16, xzr
	ldp x5, x4, [x29, #-72]
	adcs x17, x17, xzr
	sub x9, x9, x12
	cset w3, hs
	adds x11, x11, x12
	adcs x8, x8, x2
	lsl x12, x11, #32
	adcs x9, x10, x9
	adcs x10, x16, x0
	umulh x16, x11, x28
	lsr x0, x11, #32
	adcs x17, x17, xzr
	sub x11, x11, x12
	adcs x1, x1, x3
	cset w2, hs
	adds x8, x8, x12
	adcs x9, x9, x0
	umulh x12, x8, x28
	adcs x10, x10, x11
	lsl x11, x8, #32
	adcs x16, x17, x16
	lsr x17, x8, #32
	adcs x0, x1, xzr
	sub x8, x8, x11
	cset w1, hs
	adds x9, x9, x11
	adcs x10, x10, x17
	eor w1, w2, w1
	adcs x8, x16, x8
	adcs x11, x0, x12
	cset w12, hs
	adds x16, x9, #1
	sbcs x17, x10, x30
	eor w12, w1, w12
	cset w0, lo
	eor w12, w12, #0x1
	subs x0, x8, x0
	sbcs x1, x11, x28
	cset w2, lo
	tst w2, w12
	ldp x2, x3, [x29, #-56]
	csel x9, x9, x16, ne
	csel x10, x10, x17, ne
	csel x8, x8, x0, ne
	csel x11, x11, x1, ne
	mov w1, #-1
	mov x0, #-4294967295
	subs x9, x9, x2
	sbcs x10, x10, x3
	sbcs x8, x8, x4
	sbcs x11, x11, x5
	cset w12, lo
	csel x16, x1, xzr, lo
	csel x17, x0, xzr, lo
	sbfx x12, x12, #0, #1
	adds x9, x9, x12
	adcs x10, x10, x16
	adcs x8, x8, xzr
	adc x11, x17, x11
	subs x9, x9, x2
	sbcs x10, x10, x3
	sbcs x8, x8, x4
	sbcs x11, x11, x5
	cset w12, lo
	csel x17, x1, xzr, lo
	csel x0, x0, xzr, lo
	sbfx x12, x12, #0, #1
	adds x16, x9, x12
	adcs x6, x10, x17
	adcs x1, x8, xzr
	adc x17, x0, x11
	subs x8, x2, x16
	sbcs x20, x3, x6
	sbcs x11, x4, x1
	sbcs x10, x5, x17
	cset w9, lo
	csel x0, x19, xzr, lo
	csel x12, x7, xzr, lo
	sbfx x9, x9, #0, #1
	adds x7, x8, x9
	adcs x20, x20, x0
	umulh x8, x7, x18
	adcs x28, x11, xzr
	mul x24, x20, x18
	adc x12, x12, x10
	mul x4, x7, x15
	umulh x25, x20, x18
	adds x8, x8, x24
	umulh x3, x7, x15
	mul x11, x28, x18
	cinc x25, x25, hs
	adds x10, x8, x4
	umulh x30, x28, x18
	adcs x3, x25, x3
	mul x26, x20, x15
	umulh x27, x20, x15
	cinc x4, x30, hs
	adds x11, x3, x11
	mul x19, x7, x14
	adcs x4, x4, xzr
	mrs x25, NZCV
	umulh x5, x7, x14
	adds x11, x11, x26
	adcs x4, x4, x27
	mul x3, x12, x18
	cset w26, hs
	msr NZCV, x25
	cinc x25, x26, hs
	adds x11, x11, x19
	mul x9, x7, x18
	adcs x4, x4, x5
	umulh x18, x12, x18
	cinc x5, x25, hs
	adds x3, x4, x3
	mul x2, x28, x15
	umulh x24, x28, x15
	adcs x18, x5, x18
	mul x0, x20, x14
	mrs x5, NZCV
	adds x2, x3, x2
	umulh x21, x20, x14
	adcs x18, x18, x24
	mul x23, x7, x13
	cset w24, hs
	msr NZCV, x5
	cinc x5, x24, hs
	adds x0, x2, x0
	umulh x7, x7, x13
	adcs x18, x18, x21
	mul x3, x12, x15
	cset w21, hs
	adds x0, x0, x23
	umulh x15, x12, x15
	adcs x18, x18, x7
	mul x8, x28, x14
	adc x5, x5, x21
	adds x18, x18, x3
	umulh x19, x28, x14
	adcs x15, x5, x15
	mul x22, x20, x13
	cset w3, hs
	adds x8, x18, x8
	umulh x20, x20, x13
	adcs x15, x15, x19
	mul x7, x12, x14
	cset w5, hs
	adds x8, x8, x22
	umulh x14, x12, x14
	mul x4, x28, x13
	umulh x18, x12, x13
	umulh x2, x28, x13
	mul x12, x12, x13
	adcs x13, x15, x20
	adc x15, x5, x3
	adds x13, x13, x7
	adcs x14, x15, x14
	cinc x15, x18, hs
	adds x13, x13, x4
	mov x4, #-4294967295
	lsl x18, x9, #32
	umulh x3, x9, x4
	adcs x14, x14, x2
	cset w2, hs
	adds x12, x14, x12
	lsr x14, x9, #32
	adc x15, x15, x2
	adds x10, x10, x18
	sub x9, x9, x18
	adcs x11, x11, x14
	lsl x14, x10, #32
	umulh x18, x10, x4
	adcs x9, x0, x9
	lsr x0, x10, #32
	adcs x8, x8, x3
	sub x10, x10, x14
	mov x3, #-4294967295
	adcs x13, x13, xzr
	adcs x12, x12, xzr
	cset w2, hs
	adds x11, x11, x14
	adcs x9, x9, x0
	umulh x14, x11, x3
	adcs x8, x8, x10
	adcs x10, x13, x18
	lsl x13, x11, #32
	lsr x18, x11, #32
	adcs x12, x12, xzr
	adcs x15, x15, x2
	sub x11, x11, x13
	mov w2, #-1
	cset w0, hs
	adds x9, x9, x13
	adcs x8, x8, x18
	umulh x13, x9, x3
	adcs x10, x10, x11
	lsl x11, x9, #32
	adcs x12, x12, x14
	lsr x14, x9, #32
	adcs x15, x15, xzr
	sub x9, x9, x11
	cset w18, hs
	adds x8, x8, x11
	adcs x10, x10, x14
	eor w18, w0, w18
	adcs x9, x12, x9
	adcs x11, x15, x13
	cset w12, hs
	adds x13, x8, #1
	sbcs x14, x10, x2
	eor w12, w18, w12
	cset w15, lo
	eor w12, w12, #0x1
	subs x15, x9, x15
	sbcs x18, x11, x3
	cset w0, lo
	tst w0, w12
	ldur x0, [x29, #-80]
	csel x8, x8, x13, ne
	csel x10, x10, x14, ne
	csel x9, x9, x15, ne
	ldp x12, x13, [x0, #32]
	csel x11, x11, x18, ne
	stp xzr, xzr, [x0, #64]
	stp xzr, xzr, [x0, #96]
	subs x8, x8, x12
	ldp x12, x14, [x0, #48]
	sbcs x10, x10, x13
	stp xzr, xzr, [x0, #80]
	stp xzr, xzr, [x0, #112]
	sbcs x9, x9, x12
	stp x16, x6, [x0]
	sbcs x11, x11, x14
	stp x1, x17, [x0, #16]
	cset w12, lo
	csel x13, x2, xzr, lo
	csel x14, x3, xzr, lo
	sbfx x12, x12, #0, #1
	adds x8, x8, x12
	adcs x10, x10, x13
	stp x8, x10, [x0, #32]
	adcs x8, x9, xzr
	adc x9, x14, x11
	stp x8, x9, [x0, #48]
	ldp x8, x10, [x0, #136]
	adds x9, x8, #1
	cinc x8, x10, hs
	ldr x10, [sp, #32]
	stp xzr, x9, [x0, #128]
	stp x8, xzr, [x0, #144]
	cmp x8, x10
	b.hs .LBB7_9
	cbz x6, .LBB7_9
	orr x8, x1, x16
	orr x8, x8, x17
	cbz x8, .LBB7_9
	tst x9, #0x3f
	b.ne .LBB7_7
	stp x16, x1, [x0, #112]
	str x17, [x0, #128]
.LBB7_7:
	mov x9, xzr
	stp xzr, xzr, [x0, #96]
	stp x16, x1, [x0, #64]
	stp x17, xzr, [x0, #80]
.LBB7_8:
	ldur x8, [x29, #-24]
	ldp x11, x13, [x29, #-40]
	add x8, x9, x8
	stur x8, [x29, #-24]
	ldur x8, [x29, #-16]
	sub x11, x11, #32
	sub x13, x13, #32
	cbnz x8, .LBB7_2
	b .LBB7_11
.LBB7_9:
	bl slow_path
	mov w9, #1
	b .LBB7_8
.LBB7_10:
	stur xzr, [x29, #-24]
.LBB7_11:
	ldur x0, [x29, #-24]
	ldp x20, x19, [sp, #320]
	ldp x22, x21, [sp, #304]
	ldp x24, x23, [sp, #288]
	ldp x26, x25, [sp, #272]
	ldp x28, x27, [sp, #256]
	ldp x29, x30, [sp, #240]
	add sp, sp, #336
	ret
.Lfunc_end7:

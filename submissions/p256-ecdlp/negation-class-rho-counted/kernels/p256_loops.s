// clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S p256_walk.c
// Apple clang 21.0.0 (LLVM). Comments and assembler directives stripped; functions: pass1, pass2, pass1_dbl, pass2_dbl
pass1:
	cbz x5, .LBB4_4
	sub sp, sp, #128
	mov x8, x0
	mov x0, xzr
	add x9, x2, #16
	str x30, [sp, #32]
	stp x28, x27, [sp, #48]
	stp x26, x25, [sp, #64]
	stp x24, x23, [sp, #80]
	stp x22, x21, [sp, #96]
	stp x20, x19, [sp, #112]
	str x1, [sp, #8]
.LBB4_2:
	ldp x14, x13, [x8, #32]
	mov w10, #-1
	ldp x16, x17, [x8, #48]
	str x5, [sp, #40]
	stp x9, x0, [sp, #16]
	subs x15, x10, x13
	mvn x2, x14
	mov x10, #-4294967295
	ngc x18, xzr
	subs x6, x18, x16
	sbc x18, x18, x17
	cmp x2, x14
	ngc x14, xzr
	subs x13, x14, x13
	sbc x14, x14, xzr
	cmn x13, x15
	cinc x13, x14, hs
	subs x14, x6, x16
	asr x15, x13, #63
	ngc x16, xzr
	cmn x14, x13
	add x13, x18, x10
	ldp x2, x6, [x8, #16]
	adc x14, x16, x15
	subs x13, x13, x17
	asr x15, x14, #63
	ngc x16, xzr
	cmn x13, x14
	adc x14, x16, x15
	ldp x13, x15, [x8]
	and x14, x14, #0x1
	bfi x14, x13, #1, #11
	add x16, x1, x14, lsl #6
	str x14, [x8, #152]
	add x8, x8, #160
	ldp x14, x17, [x16]
	subs x13, x13, x14
	ldp x14, x16, [x16, #16]
	ngc x18, xzr
	subs x15, x15, x17
	ngc x17, xzr
	adds x15, x15, x18
	adc x17, x17, x18
	subs x14, x2, x14
	asr x18, x17, #63
	ngc x2, xzr
	adds x17, x17, x14
	adc x14, x18, x2
	subs x16, x6, x16
	asr x18, x14, #63
	ngc x2, xzr
	adds x6, x14, x16
	adc x16, x18, x2
	adds x13, x16, x13
	mov w14, w16
	adcs x14, x14, x15
	and x15, x16, #0xffffffff00000001
	adcs x16, x17, xzr
	stp x13, x14, [x9, #-16]
	adc x15, x15, x6
	stp x16, x15, [x9]
	ldp x18, x2, [x4]
	umulh x6, x13, x18
	mul x19, x14, x18
	umulh x7, x14, x18
	mul x21, x16, x18
	adds x6, x19, x6
	umulh x20, x16, x18
	cinc x7, x7, hs
	mul x23, x15, x18
	adds x7, x21, x7
	umulh x22, x15, x18
	cinc x20, x20, hs
	mul x25, x13, x2
	adds x20, x23, x20
	umulh x24, x13, x2
	cinc x22, x22, hs
	mul x27, x14, x2
	adds x6, x6, x25
	umulh x26, x14, x2
	cinc x23, x24, hs
	mul x30, x16, x2
	adds x7, x7, x27
	mul x17, x13, x18
	ldp x1, x18, [x4, #16]
	cinc x24, x26, hs
	adds x7, x7, x23
	umulh x28, x16, x2
	cinc x24, x24, hs
	adds x20, x20, x30
	mul x10, x15, x2
	umulh x12, x15, x2
	cinc x25, x28, hs
	adds x20, x20, x24
	mul x2, x13, x1
	cinc x25, x25, hs
	adds x10, x22, x10
	umulh x11, x13, x1
	cinc x12, x12, hs
	adds x10, x10, x25
	mul x5, x14, x1
	cinc x12, x12, hs
	adds x2, x7, x2
	umulh x9, x14, x1
	cinc x11, x11, hs
	mul x19, x16, x1
	adds x11, x11, x5
	umulh x0, x16, x1
	cinc x9, x9, hs
	adds x11, x11, x20
	umulh x21, x15, x1
	cinc x9, x9, hs
	adds x10, x10, x19
	mul x1, x15, x1
	cinc x0, x0, hs
	adds x9, x10, x9
	mul x24, x13, x18
	cinc x0, x0, hs
	umulh x23, x13, x18
	adds x12, x12, x1
	mul x7, x14, x18
	cinc x1, x21, hs
	adds x12, x12, x0
	cinc x1, x1, hs
	adds x11, x11, x24
	umulh x22, x14, x18
	cinc x19, x23, hs
	mul x10, x16, x18
	adds x9, x9, x7
	umulh x5, x16, x18
	cinc x7, x22, hs
	adds x9, x9, x19
	mov x19, #-4294967295
	mul x0, x15, x18
	cinc x7, x7, hs
	adds x10, x12, x10
	lsl x12, x17, #32
	umulh x18, x15, x18
	cinc x5, x5, hs
	adds x10, x10, x7
	umulh x7, x17, x19
	cinc x5, x5, hs
	adds x0, x1, x0
	lsr x1, x17, #32
	cinc x18, x18, hs
	adds x0, x0, x5
	cinc x18, x18, hs
	adds x5, x12, x6
	sub x12, x17, x12
	adcs x17, x2, x1
	lsr x2, x5, #32
	adcs x12, x12, xzr
	cinc x1, x7, hs
	adds x11, x12, x11
	cinc x12, x1, hs
	umulh x1, x5, x19
	adds x9, x9, x12
	lsl x12, x5, #32
	adcs x10, x10, xzr
	adcs x0, x0, xzr
	adcs x18, x18, xzr
	cset w6, hs
	adds x17, x12, x17
	sub x12, x5, x12
	adcs x11, x11, x2
	lsr x2, x17, #32
	adcs x12, x12, xzr
	cinc x1, x1, hs
	adds x9, x12, x9
	cinc x12, x1, hs
	umulh x1, x17, x19
	adds x10, x12, x10
	lsl x12, x17, #32
	adcs x0, x0, xzr
	adcs x18, x18, xzr
	cset w5, hs
	adds x11, x12, x11
	sub x12, x17, x12
	adcs x9, x9, x2
	lsr x2, x11, #32
	adcs x12, x12, xzr
	cinc x17, x1, hs
	adds x10, x12, x10
	lsl x12, x11, #32
	umulh x1, x11, x19
	cinc x17, x17, hs
	adds x7, x17, x0
	sub x11, x11, x12
	adds x9, x12, x9
	orr x12, x14, x13
	adcs x10, x10, x2
	orr x12, x12, x16
	mov w2, #-1
	adcs x11, x11, xzr
	orr x12, x12, x15
	cinc x13, x1, hs
	adds x11, x11, x7
	ldr x1, [sp, #8]
	cinc x13, x13, hs
	cmn x17, x0
	ldr x0, [sp, #24]
	adcs x13, x13, x18
	mov x17, #-1
	adc w14, w5, w6
	cmp x12, #0
	ldr x5, [sp, #40]
	cinc x0, x0, eq
	adds x12, x9, #1
	adc x15, xzr, x17
	asr x16, x15, #63
	adds x15, x10, x15
	cinc x16, x16, hs
	adds x15, x15, x19
	adc x16, x16, x17
	asr x17, x16, #63
	adds x16, x16, x11
	cinc x17, x17, hs
	asr x18, x17, #63
	adds x17, x13, x17
	cinc w18, w18, hs
	adds x17, x17, x2
	adc w18, w18, w2
	bic w14, w18, w14
	tst x14, #0x1
	csel x9, x12, x9, eq
	csel x10, x15, x10, eq
	csel x11, x16, x11, eq
	csel x12, x17, x13, eq
	stp x9, x10, [x4]
	subs x5, x5, #1
	stp x11, x12, [x4, #16]
	stp x11, x12, [x3, #16]
	stp x9, x10, [x3], #32
	ldr x9, [sp, #16]
	add x9, x9, #32
	b.ne .LBB4_2
	ldp x20, x19, [sp, #112]
	ldr x30, [sp, #32]
	ldp x22, x21, [sp, #96]
	ldp x24, x23, [sp, #80]
	ldp x26, x25, [sp, #64]
	ldp x28, x27, [sp, #48]
	add sp, sp, #128
	ret
.LBB4_4:
	mov x0, xzr
	ret
.Lfunc_end4:
pass2:
	sub sp, sp, #208
	stp x29, x30, [sp, #112]
	add x29, sp, #112
	stp x28, x27, [sp, #128]
	stp x26, x25, [sp, #144]
	stp x24, x23, [sp, #160]
	stp x22, x21, [sp, #176]
	stp x20, x19, [sp, #192]
	str x7, [sp, #24]
	stur x5, [x29, #-32]
	stp x1, x4, [sp, #8]
	cbz x6, .LBB5_8
	mov w9, #160
	mov w10, #1
	stur xzr, [x29, #-24]
	madd x8, x6, x9, x0
	lsl x9, x6, #5
	sub x11, x10, x6
	add x10, x9, x3
	add x9, x9, x2
	sub x10, x10, #64
	sub x15, x9, #16
	sub x25, x8, #160
.LBB5_2:
	ldur x8, [x29, #-32]
	ldp x13, x12, [x15, #-16]
	stp x10, x11, [x29, #-16]
	ldp x11, x10, [x15]
	stur x15, [x29, #-40]
	ldp x16, x14, [x8]
	ldp x15, x8, [x8, #16]
	umulh x1, x13, x16
	mul x2, x12, x16
	umulh x18, x12, x16
	mul x4, x11, x16
	adds x1, x2, x1
	umulh x3, x11, x16
	cinc x18, x18, hs
	mul x7, x10, x16
	adds x4, x18, x4
	umulh x6, x10, x16
	cinc x18, x3, hs
	mul x21, x13, x14
	adds x7, x18, x7
	umulh x20, x13, x14
	cinc x6, x6, hs
	mul x24, x12, x14
	adds x18, x1, x21
	umulh x22, x12, x14
	cinc x20, x20, hs
	mul x26, x11, x14
	adds x4, x4, x24
	umulh x0, x11, x14
	cinc x21, x22, hs
	adds x4, x4, x20
	mul x28, x10, x14
	cinc x21, x21, hs
	adds x7, x7, x26
	umulh x27, x10, x14
	cinc x0, x0, hs
	adds x7, x7, x21
	mul x23, x13, x15
	cinc x0, x0, hs
	adds x6, x6, x28
	umulh x30, x13, x15
	cinc x21, x27, hs
	adds x0, x6, x0
	mul x19, x12, x15
	cinc x21, x21, hs
	adds x4, x4, x23
	umulh x17, x12, x15
	cinc x22, x30, hs
	mul x2, x11, x15
	adds x7, x7, x19
	umulh x5, x11, x15
	cinc x17, x17, hs
	adds x7, x7, x22
	mul x1, x10, x15
	cinc x17, x17, hs
	adds x0, x0, x2
	umulh x3, x10, x15
	cinc x5, x5, hs
	adds x17, x0, x17
	mul x9, x13, x16
	cinc x0, x5, hs
	adds x1, x21, x1
	umulh x20, x13, x8
	cinc x3, x3, hs
	adds x0, x1, x0
	mul x13, x13, x8
	cinc x1, x3, hs
	umulh x6, x12, x8
	mul x12, x12, x8
	adds x13, x7, x13
	umulh x2, x11, x8
	cinc x5, x20, hs
	mul x11, x11, x8
	adds x12, x17, x12
	mul x3, x10, x8
	cinc x17, x6, hs
	adds x12, x12, x5
	cinc x17, x17, hs
	mov x5, #-4294967295
	umulh x10, x10, x8
	adds x11, x0, x11
	cinc x0, x2, hs
	adds x11, x11, x17
	mov x2, #-4294967295
	cinc x17, x0, hs
	adds x0, x1, x3
	lsl x1, x9, #32
	umulh x2, x9, x2
	cinc x10, x10, hs
	adds x17, x0, x17
	lsr x0, x9, #32
	cinc x10, x10, hs
	adds x18, x1, x18
	sub x9, x9, x1
	adcs x0, x4, x0
	adcs x9, x9, xzr
	cinc x1, x2, hs
	adds x9, x9, x13
	lsr x2, x18, #32
	cinc x13, x1, hs
	umulh x1, x18, x5
	adds x12, x12, x13
	adcs x11, x11, xzr
	adcs x13, x17, xzr
	lsl x17, x18, #32
	adcs x10, x10, xzr
	cset w3, hs
	adds x0, x17, x0
	sub x17, x18, x17
	adcs x9, x9, x2
	adcs x17, x17, xzr
	cinc x18, x1, hs
	adds x12, x17, x12
	cinc x17, x18, hs
	umulh x18, x0, x5
	adds x11, x17, x11
	lsl x17, x0, #32
	adcs x13, x13, xzr
	adcs x1, x10, xzr
	lsr x10, x0, #32
	cset w2, hs
	adds x9, x17, x9
	sub x17, x0, x17
	adcs x10, x12, x10
	umulh x0, x9, x5
	lsr x5, x9, #32
	adcs x12, x17, xzr
	cinc x17, x18, hs
	adds x11, x12, x11
	lsl x18, x9, #32
	cinc x17, x17, hs
	adds x4, x17, x13
	sub x9, x9, x18
	adds x12, x18, x10
	adcs x10, x11, x5
	adcs x9, x9, xzr
	cinc x18, x0, hs
	adds x9, x9, x4
	stp x10, x9, [sp, #48]
	cinc x9, x18, hs
	cmn x17, x13
	adcs x13, x9, x1
	ldr x17, [sp, #16]
	adc w9, w2, w3
	stur w9, [x29, #-48]
	ldp x9, x10, [x29, #-16]
	cmp x10, #0
	csel x18, x17, x9, eq
	ldp x6, x3, [x18]
	ldp x2, x18, [x18, #16]
	umulh x0, x6, x16
	mul x4, x3, x16
	umulh x1, x3, x16
	mul x7, x2, x16
	adds x0, x4, x0
	umulh x5, x2, x16
	cinc x1, x1, hs
	mul x17, x6, x16
	adds x1, x1, x7
	umulh x19, x18, x16
	cinc x5, x5, hs
	mul x16, x18, x16
	mul x21, x6, x14
	umulh x20, x6, x14
	adds x5, x5, x16
	mul x23, x3, x14
	cinc x19, x19, hs
	adds x16, x0, x21
	umulh x22, x3, x14
	cinc x0, x20, hs
	mul x26, x2, x14
	adds x1, x1, x23
	umulh x24, x2, x14
	cinc x20, x22, hs
	adds x0, x1, x0
	mul x28, x18, x14
	cinc x20, x20, hs
	adds x5, x5, x26
	umulh x27, x18, x14
	cinc x21, x24, hs
	adds x5, x5, x20
	mul x9, x6, x15
	cinc x20, x21, hs
	adds x19, x19, x28
	umulh x30, x6, x15
	cinc x21, x27, hs
	adds x19, x19, x20
	mul x14, x3, x15
	cinc x21, x21, hs
	adds x9, x0, x9
	umulh x10, x3, x15
	cinc x0, x30, hs
	mul x4, x2, x15
	adds x14, x5, x14
	umulh x11, x2, x15
	cinc x10, x10, hs
	adds x14, x14, x0
	umulh x7, x18, x15
	cinc x10, x10, hs
	adds x0, x19, x4
	mul x15, x18, x15
	cinc x11, x11, hs
	adds x10, x0, x10
	umulh x1, x6, x8
	cinc x11, x11, hs
	mul x6, x6, x8
	adds x15, x21, x15
	umulh x20, x3, x8
	mul x3, x3, x8
	mul x0, x2, x8
	umulh x4, x2, x8
	cinc x2, x7, hs
	adds x11, x15, x11
	cinc x15, x2, hs
	adds x14, x14, x6
	mul x2, x18, x8
	cinc x1, x1, hs
	adds x10, x10, x3
	cinc x3, x20, hs
	adds x10, x10, x1
	umulh x8, x18, x8
	cinc x18, x3, hs
	adds x11, x11, x0
	cinc x0, x4, hs
	adds x11, x11, x18
	mov x4, #-4294967295
	cinc x18, x0, hs
	lsl x0, x17, #32
	mov w3, #-1
	adds x15, x15, x2
	mov x2, #-4294967295
	umulh x1, x17, x2
	cinc x8, x8, hs
	adds x15, x15, x18
	lsr x18, x17, #32
	cinc x8, x8, hs
	adds x16, x0, x16
	sub x17, x17, x0
	adcs x9, x9, x18
	lsr x0, x16, #32
	adcs x17, x17, xzr
	cinc x18, x1, hs
	adds x14, x17, x14
	cinc x17, x18, hs
	umulh x18, x16, x2
	adds x10, x10, x17
	lsl x17, x16, #32
	adcs x11, x11, xzr
	adcs x15, x15, xzr
	sub x16, x16, x17
	adcs x8, x8, xzr
	cset w1, hs
	adds x9, x17, x9
	adcs x14, x14, x0
	adcs x16, x16, xzr
	cinc x17, x18, hs
	adds x10, x16, x10
	lsr x18, x9, #32
	cinc x16, x17, hs
	umulh x17, x9, x4
	adds x11, x16, x11
	lsl x16, x9, #32
	adcs x15, x15, xzr
	adcs x8, x8, xzr
	sub x9, x9, x16
	cset w0, hs
	adds x14, x16, x14
	adcs x10, x10, x18
	lsr x18, x14, #32
	adcs x9, x9, xzr
	cinc x16, x17, hs
	adds x9, x9, x11
	lsl x11, x14, #32
	cinc x16, x16, hs
	umulh x17, x14, x4
	adds x2, x16, x15
	adds x10, x11, x10
	sub x11, x14, x11
	adcs x9, x9, x18
	adcs x11, x11, xzr
	cinc x14, x17, hs
	adds x11, x11, x2
	mov x2, #-1
	cinc x14, x14, hs
	cmn x16, x15
	adcs x8, x14, x8
	adc w14, w0, w1
	adds x15, x10, #1
	adc x16, xzr, x2
	asr x17, x16, #63
	adds x16, x9, x16
	cinc x17, x17, hs
	adds x16, x16, x4
	adc x17, x17, x2
	asr x18, x17, #63
	adds x0, x17, x11
	cinc x17, x18, hs
	asr x18, x17, #63
	adds x17, x8, x17
	cinc w18, w18, hs
	adds x1, x17, x3
	adc w17, w18, w3
	bic w14, w17, w14
	tst x14, #0x1
	csel x17, x15, x10, eq
	csel x16, x16, x9, eq
	csel x15, x0, x11, eq
	csel x14, x1, x8, eq
	adds x8, x12, #1
	ldp x0, x1, [sp, #48]
	adc x9, xzr, x2
	asr x10, x9, #63
	adds x9, x0, x9
	cinc x10, x10, hs
	adds x9, x9, x4
	adc x10, x10, x2
	ldur w2, [x29, #-48]
	asr x11, x10, #63
	adds x10, x10, x1
	cinc x11, x11, hs
	asr x18, x11, #63
	adds x11, x13, x11
	cinc w18, w18, hs
	adds x11, x11, x3
	adc w18, w18, w3
	bic w18, w18, w2
	tst x18, #0x1
	ldur x18, [x29, #-32]
	csel x8, x8, x12, eq
	csel x9, x9, x0, eq
	ldr x12, [x25, #152]
	stp x8, x9, [x18]
	ldr x8, [sp, #8]
	csel x10, x10, x1, eq
	add x7, x8, x12, lsl #6
	csel x8, x11, x13, eq
	stp x10, x8, [x18, #16]
	ldp x8, x9, [x7, #32]
	ldp x10, x11, [x25, #32]
	subs x13, x10, x8
	stur x10, [x29, #-48]
	ngc x18, xzr
	subs x9, x11, x9
	str x11, [sp, #56]
	ngc x0, xzr
	adds x9, x9, x18
	ldp x10, x11, [x25, #48]
	adc x18, x0, x18
	ldp x8, x0, [x7, #48]
	asr x1, x18, #63
	stp x10, x11, [sp, #40]
	subs x8, x10, x8
	ngc x2, xzr
	adds x18, x18, x8
	adc x8, x1, x2
	subs x0, x11, x0
	asr x1, x8, #63
	ngc x2, xzr
	adds x8, x8, x0
	adc x0, x1, x2
	adds x1, x0, x13
	mov w21, w0
	and x0, x0, #0xffffffff00000001
	umulh x2, x17, x1
	adcs x9, x21, x9
	adcs x30, x18, xzr
	mul x5, x16, x1
	adc x8, x0, x8
	umulh x3, x16, x1
	mul x19, x15, x1
	adds x0, x2, x5
	umulh x6, x15, x1
	cinc x2, x3, hs
	mul x13, x17, x1
	adds x2, x2, x19
	umulh x20, x14, x1
	cinc x5, x6, hs
	mul x1, x14, x1
	mul x22, x17, x9
	umulh x21, x17, x9
	adds x1, x5, x1
	mul x24, x16, x9
	cinc x5, x20, hs
	adds x0, x0, x22
	umulh x23, x16, x9
	cinc x20, x21, hs
	mul x27, x15, x9
	adds x2, x2, x24
	umulh x26, x15, x9
	cinc x21, x23, hs
	adds x2, x2, x20
	mul x10, x14, x9
	cinc x21, x21, hs
	adds x1, x1, x27
	umulh x28, x14, x9
	cinc x22, x26, hs
	adds x1, x1, x21
	mul x18, x17, x30
	cinc x21, x22, hs
	mov x22, #-4294967295
	adds x10, x5, x10
	ldp x26, x24, [x25, #16]
	umulh x12, x17, x30
	cinc x5, x28, hs
	adds x10, x10, x21
	mul x9, x16, x30
	cinc x5, x5, hs
	adds x18, x2, x18
	umulh x11, x16, x30
	cinc x12, x12, hs
	mul x3, x15, x30
	adds x9, x1, x9
	umulh x4, x15, x30
	cinc x11, x11, hs
	adds x9, x9, x12
	mul x19, x14, x30
	cinc x11, x11, hs
	adds x10, x10, x3
	umulh x6, x14, x30
	cinc x1, x4, hs
	adds x10, x10, x11
	umulh x20, x17, x8
	cinc x11, x1, hs
	adds x1, x5, x19
	mul x17, x17, x8
	cinc x2, x6, hs
	adds x11, x1, x11
	umulh x21, x16, x8
	cinc x1, x2, hs
	mul x16, x16, x8
	adds x9, x9, x17
	umulh x12, x15, x8
	cinc x2, x20, hs
	mul x15, x15, x8
	adds x10, x10, x16
	mul x17, x14, x8
	cinc x16, x21, hs
	adds x10, x10, x2
	umulh x8, x14, x8
	cinc x14, x16, hs
	adds x11, x11, x15
	lsl x15, x13, #32
	cinc x12, x12, hs
	adds x11, x11, x14
	cinc x12, x12, hs
	adds x14, x1, x17
	mov x1, #-4294967295
	umulh x16, x13, x1
	cinc x8, x8, hs
	adds x12, x14, x12
	lsr x14, x13, #32
	cinc x8, x8, hs
	adds x17, x15, x0
	sub x13, x13, x15
	adcs x14, x18, x14
	adcs x13, x13, xzr
	cinc x15, x16, hs
	adds x9, x13, x9
	lsr x16, x17, #32
	cinc x13, x15, hs
	umulh x15, x17, x1
	adds x10, x10, x13
	lsl x13, x17, #32
	adcs x11, x11, xzr
	adcs x12, x12, xzr
	adcs x8, x8, xzr
	cset w18, hs
	adds x14, x13, x14
	sub x13, x17, x13
	adcs x9, x9, x16
	lsr x16, x14, #32
	adcs x13, x13, xzr
	cinc x15, x15, hs
	adds x10, x13, x10
	cinc x13, x15, hs
	umulh x15, x14, x22
	adds x11, x13, x11
	lsl x13, x14, #32
	adcs x12, x12, xzr
	adcs x8, x8, xzr
	cset w17, hs
	adds x9, x13, x9
	sub x13, x14, x13
	adcs x10, x10, x16
	lsr x16, x9, #32
	adcs x13, x13, xzr
	cinc x14, x15, hs
	adds x11, x13, x11
	lsl x13, x9, #32
	umulh x15, x9, x22
	cinc x14, x14, hs
	adds x0, x14, x12
	sub x9, x9, x13
	adds x10, x13, x10
	adcs x11, x11, x16
	mov x16, #-1
	adcs x9, x9, xzr
	cinc x13, x15, hs
	adds x9, x9, x0
	mov w0, #-1
	cinc x13, x13, hs
	cmn x14, x12
	adcs x8, x13, x8
	adc w12, w17, w18
	adds x13, x10, #1
	adc x14, xzr, x16
	asr x15, x14, #63
	adds x14, x11, x14
	cinc x15, x15, hs
	adds x14, x14, x22
	adc x15, x15, x16
	asr x16, x15, #63
	adds x17, x15, x9
	cinc x15, x16, hs
	asr x16, x15, #63
	adds x15, x8, x15
	cinc w16, w16, hs
	adds x18, x15, x0
	adc w15, w16, w0
	bic w12, w15, w12
	tst x12, #0x1
	csel x16, x13, x10, eq
	csel x15, x14, x11, eq
	csel x14, x17, x9, eq
	umulh x10, x16, x16
	csel x13, x18, x8, eq
	mul x12, x16, x15
	umulh x11, x16, x15
	mul x6, x16, x14
	adds x10, x10, x12
	umulh x3, x16, x14
	cinc x8, x11, hs
	mul x2, x16, x13
	adds x19, x8, x6
	umulh x17, x16, x13
	cinc x9, x3, hs
	mul x1, x15, x15
	adds x20, x9, x2
	umulh x0, x15, x15
	cinc x21, x17, hs
	adds x18, x10, x12
	mul x5, x15, x14
	cinc x10, x11, hs
	adds x11, x19, x1
	umulh x4, x15, x14
	cinc x0, x0, hs
	adds x10, x11, x10
	mul x9, x15, x13
	cinc x11, x0, hs
	adds x0, x20, x5
	umulh x8, x15, x13
	cinc x19, x4, hs
	adds x11, x0, x11
	mul x1, x14, x14
	cinc x0, x19, hs
	adds x19, x21, x9
	umulh x12, x14, x14
	cinc x21, x8, hs
	adds x0, x19, x0
	mul x20, x14, x13
	cinc x19, x21, hs
	adds x10, x10, x6
	cinc x3, x3, hs
	adds x11, x11, x5
	umulh x6, x14, x13
	cinc x4, x4, hs
	adds x11, x11, x3
	cinc x3, x4, hs
	adds x0, x0, x1
	cinc x12, x12, hs
	adds x0, x0, x3
	mul x1, x16, x16
	cinc x12, x12, hs
	adds x3, x19, x20
	mul x4, x13, x13
	cinc x5, x6, hs
	adds x12, x3, x12
	cinc x3, x5, hs
	adds x11, x11, x2
	umulh x2, x13, x13
	cinc x17, x17, hs
	adds x9, x0, x9
	cinc x8, x8, hs
	adds x9, x9, x17
	lsl x17, x1, #32
	cinc x8, x8, hs
	adds x12, x12, x20
	cinc x0, x6, hs
	adds x8, x12, x8
	umulh x12, x1, x22
	cinc x0, x0, hs
	adds x3, x3, x4
	lsr x4, x1, #32
	cinc x2, x2, hs
	adds x0, x3, x0
	cinc x2, x2, hs
	adds x18, x17, x18
	sub x17, x1, x17
	adcs x10, x10, x4
	lsr x1, x18, #32
	adcs x17, x17, xzr
	cinc x12, x12, hs
	adds x11, x17, x11
	umulh x17, x18, x22
	cinc x12, x12, hs
	adds x9, x9, x12
	lsl x12, x18, #32
	adcs x8, x8, xzr
	adcs x0, x0, xzr
	adcs x2, x2, xzr
	cset w3, hs
	adds x10, x12, x10
	sub x12, x18, x12
	adcs x11, x11, x1
	adcs x12, x12, xzr
	cinc x17, x17, hs
	adds x9, x12, x9
	cinc x12, x17, hs
	umulh x17, x10, x22
	adds x8, x12, x8
	lsl x12, x10, #32
	adcs x18, x0, xzr
	lsr x0, x10, #32
	adcs x1, x2, xzr
	sub x10, x10, x12
	cset w2, hs
	adds x11, x12, x11
	adcs x9, x9, x0
	lsr x0, x11, #32
	adcs x10, x10, xzr
	cinc x12, x17, hs
	adds x8, x10, x8
	lsl x10, x11, #32
	umulh x17, x11, x22
	cinc x12, x12, hs
	adds x4, x12, x18
	adds x9, x10, x9
	sub x10, x11, x10
	adcs x8, x8, x0
	adcs x10, x10, xzr
	cinc x11, x17, hs
	adds x10, x10, x4
	cinc x11, x11, hs
	cmn x12, x18
	adcs x11, x11, x1
	mov x1, #-1
	adc w12, w2, w3
	adds x17, x9, #1
	mov w3, #-1
	adc x18, xzr, x1
	asr x0, x18, #63
	adds x18, x8, x18
	cinc x0, x0, hs
	adds x18, x18, x22
	adc x0, x0, x1
	asr x1, x0, #63
	adds x0, x0, x10
	cinc x1, x1, hs
	asr x2, x1, #63
	adds x1, x11, x1
	cinc w2, w2, hs
	adds x1, x1, x3
	adc w2, w2, w3
	ldp x3, x6, [x25]
	bic w12, w2, w12
	tst x12, #0x1
	csel x9, x17, x9, eq
	csel x8, x18, x8, eq
	csel x10, x0, x10, eq
	csel x11, x1, x11, eq
	subs x9, x9, x3
	ngc x12, xzr
	subs x8, x8, x6
	ngc x17, xzr
	adds x8, x12, x8
	adc x12, x12, x17
	subs x10, x10, x26
	asr x17, x12, #63
	ngc x18, xzr
	adds x10, x12, x10
	adc x12, x17, x18
	subs x11, x11, x24
	asr x17, x12, #63
	ngc x18, xzr
	adds x11, x12, x11
	adc x12, x17, x18
	adds x9, x12, x9
	mov w17, w12
	and x12, x12, #0xffffffff00000001
	adcs x8, x17, x8
	ldp x17, x18, [x7]
	adcs x10, x10, xzr
	adc x11, x12, x11
	subs x9, x9, x17
	ngc x12, xzr
	subs x8, x8, x18
	ngc x18, xzr
	adds x8, x8, x12
	adc x12, x18, x12
	ldp x17, x18, [x7, #16]
	subs x10, x10, x17
	asr x17, x12, #63
	ngc x0, xzr
	adds x10, x12, x10
	adc x12, x17, x0
	subs x11, x11, x18
	asr x17, x12, #63
	ngc x18, xzr
	adds x11, x11, x12
	adc x17, x18, x17
	adds x2, x17, x9
	mov w9, w17
	adcs x12, x9, x8
	and x8, x17, #0xffffffff00000001
	adcs x18, x10, xzr
	str x12, [sp, #32]
	adc x17, x8, x11
	subs x8, x3, x2
	ngc x9, xzr
	subs x10, x6, x12
	stp x18, x17, [x25, #16]
	ngc x11, xzr
	adds x27, x10, x9
	adc x9, x11, x9
	subs x10, x26, x18
	asr x11, x9, #63
	ngc x0, xzr
	adds x6, x9, x10
	adc x9, x11, x0
	subs x10, x24, x17
	asr x11, x9, #63
	ngc x0, xzr
	adds x26, x10, x9
	adc x28, x0, x11
	adds x9, x28, x8
	mov w10, w28
	and x23, x28, #0xffffffff00000001
	adcs x10, x10, x27
	umulh x30, x9, x16
	mul x0, x10, x16
	umulh x11, x10, x16
	umulh x4, x10, x15
	mul x5, x10, x15
	umulh x7, x10, x14
	mul x19, x10, x14
	umulh x27, x10, x13
	mul x22, x10, x13
	adcs x10, x6, xzr
	adc x23, x23, x26
	adds x0, x0, x30
	mul x28, x10, x16
	cinc x11, x11, hs
	umulh x6, x10, x16
	mul x3, x9, x16
	adds x11, x28, x11
	umulh x26, x16, x23
	cinc x6, x6, hs
	mul x16, x16, x23
	mul x8, x9, x15
	umulh x1, x9, x15
	adds x16, x16, x6
	mul x30, x10, x15
	cinc x6, x26, hs
	adds x8, x0, x8
	umulh x12, x10, x15
	cinc x0, x1, hs
	adds x11, x11, x5
	umulh x28, x15, x23
	cinc x1, x4, hs
	adds x11, x11, x0
	cinc x1, x1, hs
	adds x16, x16, x30
	mul x15, x15, x23
	cinc x12, x12, hs
	adds x16, x16, x1
	mul x20, x9, x14
	cinc x12, x12, hs
	umulh x21, x9, x14
	adds x15, x6, x15
	mul x1, x10, x14
	cinc x4, x28, hs
	adds x12, x15, x12
	cinc x4, x4, hs
	adds x11, x11, x20
	umulh x0, x10, x14
	cinc x5, x21, hs
	umulh x15, x14, x23
	adds x5, x5, x19
	cinc x6, x7, hs
	adds x16, x5, x16
	mul x14, x14, x23
	cinc x6, x6, hs
	adds x12, x12, x1
	cinc x0, x0, hs
	adds x12, x12, x6
	umulh x24, x9, x13
	cinc x0, x0, hs
	mul x9, x9, x13
	adds x14, x4, x14
	mov x4, #-4294967295
	umulh x5, x10, x13
	cinc x15, x15, hs
	adds x14, x14, x0
	cinc x15, x15, hs
	mul x10, x10, x13
	adds x9, x16, x9
	mul x0, x13, x23
	cinc x16, x24, hs
	adds x12, x12, x22
	cinc x1, x27, hs
	adds x12, x12, x16
	umulh x13, x13, x23
	cinc x16, x1, hs
	adds x10, x14, x10
	lsl x14, x3, #32
	cinc x1, x5, hs
	adds x10, x10, x16
	umulh x16, x3, x4
	cinc x1, x1, hs
	adds x15, x15, x0
	lsr x0, x3, #32
	cinc x13, x13, hs
	adds x15, x15, x1
	cinc x13, x13, hs
	adds x8, x14, x8
	sub x14, x3, x14
	adcs x11, x11, x0
	lsr x0, x8, #32
	adcs x14, x14, xzr
	cinc x16, x16, hs
	adds x9, x14, x9
	cinc x14, x16, hs
	umulh x16, x8, x4
	adds x12, x12, x14
	lsl x14, x8, #32
	adcs x10, x10, xzr
	adcs x15, x15, xzr
	sub x8, x8, x14
	adcs x13, x13, xzr
	cset w1, hs
	adds x11, x14, x11
	adcs x9, x9, x0
	adcs x8, x8, xzr
	cinc x14, x16, hs
	adds x8, x8, x12
	lsr x16, x11, #32
	cinc x12, x14, hs
	umulh x14, x11, x4
	adds x10, x12, x10
	lsl x12, x11, #32
	adcs x15, x15, xzr
	adcs x13, x13, xzr
	sub x11, x11, x12
	cset w0, hs
	adds x9, x12, x9
	adcs x8, x8, x16
	lsr x16, x9, #32
	adcs x11, x11, xzr
	cinc x12, x14, hs
	adds x10, x11, x10
	lsl x11, x9, #32
	umulh x14, x9, x4
	cinc x12, x12, hs
	adds x3, x12, x15
	sub x9, x9, x11
	adds x8, x11, x8
	adcs x10, x10, x16
	mov x16, #-1
	adcs x9, x9, xzr
	cinc x11, x14, hs
	adds x9, x9, x3
	cinc x11, x11, hs
	cmn x12, x15
	adcs x11, x11, x13
	adc w12, w0, w1
	adds x13, x8, #1
	mov w1, #-1
	adc x14, xzr, x16
	asr x15, x14, #63
	adds x14, x10, x14
	cinc x15, x15, hs
	adds x14, x14, x4
	adc x15, x15, x16
	asr x16, x15, #63
	adds x15, x15, x9
	cinc x16, x16, hs
	asr x0, x16, #63
	adds x16, x11, x16
	cinc w0, w0, hs
	adds x16, x16, x1
	adc w0, w0, w1
	bic w12, w0, w12
	tst x12, #0x1
	ldur x12, [x29, #-48]
	csel x8, x13, x8, eq
	ldr x13, [sp, #56]
	csel x10, x14, x10, eq
	csel x9, x15, x9, eq
	csel x11, x16, x11, eq
	subs x8, x8, x12
	ngc x12, xzr
	subs x10, x10, x13
	ngc x13, xzr
	adds x10, x12, x10
	adc x12, x12, x13
	ldr x13, [sp, #40]
	subs x9, x9, x13
	asr x13, x12, #63
	ngc x14, xzr
	adds x9, x12, x9
	adc x12, x13, x14
	ldr x13, [sp, #48]
	subs x11, x11, x13
	asr x13, x12, #63
	ngc x14, xzr
	adds x11, x12, x11
	adc x12, x13, x14
	adds x8, x12, x8
	mov w13, w12
	and x12, x12, #0xffffffff00000001
	adcs x10, x13, x10
	adcs x14, x9, xzr
	stp x8, x10, [x25, #32]
	adc x11, x12, x11
	ldp x13, x12, [x25, #136]
	stp x14, x11, [x25, #48]
	add x8, x13, #1
	ldp x9, x13, [sp, #24]
	and x10, x12, #0x1
	str x8, [x25, #136]
	cmp x8, x9
	ldp x9, x12, [x25, #88]
	ccmp x13, #0, #4, lo
	stp x2, x13, [x25]
	csinc x10, x10, xzr, ne
	cmp x9, x2
	ccmp x12, x18, #0, eq
	ldp x9, x12, [x25, #104]
	ccmp x17, x9, #0, eq
	ldp x9, x11, [x25, #120]
	csinc x10, x10, xzr, ne
	cmp x12, x2
	ccmp x9, x18, #0, eq
	ccmp x17, x11, #0, eq
	csinc x9, x10, xzr, ne
	cbnz x9, .LBB5_7
	tst x8, #0x3f
	b.ne .LBB5_5
	stp x2, x18, [x25, #112]
	str x17, [x25, #128]
.LBB5_5:
	ldp x9, x10, [x25, #72]
	mov x8, xzr
	stp x9, x10, [x25, #96]
	ldr x9, [x25, #64]
	stp x2, x18, [x25, #64]
	stp x17, x9, [x25, #80]
.LBB5_6:
	ldp x10, x11, [x29, #-16]
	sub x25, x25, #160
	ldur x15, [x29, #-40]
	ldur x9, [x29, #-24]
	add x11, x11, #1
	sub x15, x15, #32
	add x9, x8, x9
	cmp x11, #1
	sub x10, x10, #32
	stur x9, [x29, #-24]
	b.ne .LBB5_2
	b .LBB5_9
.LBB5_7:
	mov x0, x25
	bl slow_path
	mov w8, #1
	b .LBB5_6
.LBB5_8:
	stur xzr, [x29, #-24]
.LBB5_9:
	ldur x0, [x29, #-24]
	ldp x20, x19, [sp, #192]
	ldp x22, x21, [sp, #176]
	ldp x24, x23, [sp, #160]
	ldp x26, x25, [sp, #144]
	ldp x28, x27, [sp, #128]
	ldp x29, x30, [sp, #112]
	add sp, sp, #208
	ret
.Lfunc_end5:
pass1_dbl:
	cbz x4, .LBB6_4
	sub sp, sp, #112
	mov x8, x0
	mov x0, xzr
	add x9, x1, #16
	str x30, [sp, #16]
	stp x28, x27, [sp, #32]
	stp x26, x25, [sp, #48]
	stp x24, x23, [sp, #64]
	stp x22, x21, [sp, #80]
	stp x20, x19, [sp, #96]
.LBB6_2:
	ldr x14, [x8], #8
	stp x8, x0, [sp]
	mov x8, #-4294967296
	mov x10, #-1
	str x4, [sp, #24]
	ldp x15, x16, [x14, #32]
	ldp x18, x14, [x14, #48]
	extr x17, x16, x15, #63
	extr x16, x18, x16, #63
	extr x18, x14, x18, #63
	adds x1, x17, x8
	mov x8, #-1
	adc x5, xzr, x8
	mov w8, #-1
	asr x6, x5, #63
	adds x5, x5, x16
	cinc x6, x6, hs
	asr x7, x6, #63
	adds x6, x18, x6
	cinc x7, x7, hs
	adds x6, x6, x8
	adc x7, x7, x10
	cmp x14, #0
	and x14, x7, #0x1
	csel x14, xzr, x14, lt
	orr x7, x14, x15, lsl #1
	cmp x14, #0
	csel x15, x5, x16, eq
	csel x14, x6, x18, eq
	csel x17, x1, x17, eq
	eor x16, x7, #0x1
	stp x15, x14, [x9]
	stp x16, x17, [x9, #-16]
	add x9, x9, #32
	ldp x1, x5, [x3]
	umulh x6, x16, x1
	mul x19, x17, x1
	umulh x7, x17, x1
	mul x21, x15, x1
	adds x6, x6, x19
	umulh x20, x15, x1
	cinc x7, x7, hs
	mul x23, x14, x1
	adds x7, x7, x21
	umulh x22, x14, x1
	cinc x20, x20, hs
	mul x25, x16, x5
	adds x20, x20, x23
	umulh x24, x16, x5
	cinc x22, x22, hs
	mul x27, x17, x5
	adds x6, x6, x25
	umulh x26, x17, x5
	cinc x23, x24, hs
	mul x30, x15, x5
	adds x7, x7, x27
	mul x18, x16, x1
	ldp x11, x1, [x3, #16]
	cinc x24, x26, hs
	adds x7, x7, x23
	umulh x28, x15, x5
	cinc x24, x24, hs
	adds x20, x20, x30
	mul x12, x14, x5
	umulh x10, x14, x5
	cinc x25, x28, hs
	adds x20, x20, x24
	mul x8, x16, x11
	cinc x25, x25, hs
	adds x12, x22, x12
	umulh x13, x16, x11
	cinc x10, x10, hs
	adds x12, x12, x25
	mul x5, x17, x11
	cinc x10, x10, hs
	adds x8, x7, x8
	umulh x4, x17, x11
	cinc x13, x13, hs
	mul x19, x15, x11
	adds x5, x20, x5
	umulh x0, x15, x11
	cinc x4, x4, hs
	adds x13, x5, x13
	umulh x21, x14, x11
	cinc x4, x4, hs
	adds x12, x12, x19
	mul x11, x14, x11
	cinc x0, x0, hs
	adds x12, x12, x4
	mul x24, x16, x1
	cinc x0, x0, hs
	umulh x23, x16, x1
	adds x10, x10, x11
	mul x7, x17, x1
	cinc x4, x21, hs
	adds x10, x10, x0
	cinc x0, x4, hs
	adds x13, x13, x24
	umulh x22, x17, x1
	cinc x19, x23, hs
	mul x11, x15, x1
	adds x12, x12, x7
	umulh x5, x15, x1
	cinc x7, x22, hs
	adds x12, x12, x19
	mov x19, #-4294967295
	mul x4, x14, x1
	cinc x7, x7, hs
	adds x10, x10, x11
	umulh x1, x14, x1
	cinc x11, x5, hs
	adds x10, x10, x7
	cinc x11, x11, hs
	umulh x5, x18, x19
	adds x0, x0, x4
	lsl x4, x18, #32
	cinc x1, x1, hs
	adds x11, x0, x11
	lsr x0, x18, #32
	cinc x1, x1, hs
	adds x6, x4, x6
	sub x18, x18, x4
	adcs x8, x8, x0
	lsr x4, x6, #32
	adcs x18, x18, xzr
	cinc x0, x5, hs
	adds x13, x18, x13
	cinc x18, x0, hs
	umulh x0, x6, x19
	adds x12, x12, x18
	lsl x18, x6, #32
	adcs x10, x10, xzr
	adcs x11, x11, xzr
	adcs x1, x1, xzr
	cset w5, hs
	adds x8, x18, x8
	sub x18, x6, x18
	adcs x13, x13, x4
	lsr x4, x8, #32
	adcs x18, x18, xzr
	cinc x0, x0, hs
	adds x12, x18, x12
	cinc x18, x0, hs
	umulh x0, x8, x19
	adds x10, x18, x10
	lsl x18, x8, #32
	adcs x11, x11, xzr
	adcs x1, x1, xzr
	sub x8, x8, x18
	cset w6, hs
	adds x13, x18, x13
	adcs x12, x12, x4
	lsr x4, x13, #32
	adcs x8, x8, xzr
	cinc x18, x0, hs
	adds x8, x8, x10
	lsl x10, x13, #32
	umulh x0, x13, x19
	cinc x18, x18, hs
	adds x7, x18, x11
	adds x12, x10, x12
	sub x10, x13, x10
	orr x13, x17, x15
	adcs x8, x8, x4
	orr x13, x16, x13
	mov x17, #-1
	adcs x10, x10, xzr
	ldr x4, [sp, #24]
	cinc x15, x0, hs
	adds x10, x10, x7
	ldr x0, [sp, #8]
	cinc x15, x15, hs
	cmn x18, x11
	orr x11, x13, x14
	adcs x13, x15, x1
	mov w1, #-1
	adc w14, w6, w5
	cmp x11, #0
	cinc x0, x0, eq
	adds x11, x12, #1
	adc x15, xzr, x17
	asr x16, x15, #63
	adds x15, x8, x15
	cinc x16, x16, hs
	adds x15, x15, x19
	adc x16, x16, x17
	asr x17, x16, #63
	adds x16, x16, x10
	cinc x17, x17, hs
	asr x18, x17, #63
	adds x17, x13, x17
	cinc w18, w18, hs
	adds x17, x17, x1
	adc w18, w18, w1
	bic w14, w18, w14
	tst x14, #0x1
	csel x11, x11, x12, eq
	csel x8, x15, x8, eq
	csel x10, x16, x10, eq
	csel x12, x17, x13, eq
	stp x11, x8, [x3]
	subs x4, x4, #1
	stp x10, x12, [x3, #16]
	stp x10, x12, [x2, #16]
	stp x11, x8, [x2], #32
	ldr x8, [sp]
	b.ne .LBB6_2
	ldp x20, x19, [sp, #96]
	ldr x30, [sp, #16]
	ldp x22, x21, [sp, #80]
	ldp x24, x23, [sp, #64]
	ldp x26, x25, [sp, #48]
	ldp x28, x27, [sp, #32]
	add sp, sp, #112
	ret
.LBB6_4:
	mov x0, xzr
	ret
.Lfunc_end6:
pass2_dbl:
	sub sp, sp, #320
	stp x29, x30, [sp, #224]
	add x29, sp, #224
	stp x28, x27, [sp, #240]
	stp x26, x25, [sp, #256]
	stp x24, x23, [sp, #272]
	stp x22, x21, [sp, #288]
	stp x20, x19, [sp, #304]
	stp x3, x7, [sp, #32]
	stp x5, x4, [sp, #16]
	stur x6, [x29, #-16]
	cbz x6, .LBB7_10
	ldur x9, [x29, #-16]
	sub x8, x0, #8
	stur xzr, [x29, #-24]
	str x8, [sp, #8]
	lsl x9, x9, #5
	add x8, x9, x2
	add x9, x9, x1
	sub x13, x8, #64
	sub x14, x9, #16
.LBB7_2:
	ldr x0, [sp, #16]
	ldp x12, x11, [x14, #-16]
	ldp x10, x9, [x14]
	ldp x16, x23, [x0]
	umulh x17, x12, x16
	mul x1, x11, x16
	umulh x18, x11, x16
	mul x3, x10, x16
	adds x17, x1, x17
	umulh x2, x10, x16
	cinc x18, x18, hs
	mul x5, x9, x16
	adds x18, x18, x3
	umulh x8, x10, x23
	cinc x2, x2, hs
	umulh x4, x9, x16
	adds x2, x2, x5
	mul x7, x12, x23
	stp x14, x8, [x29, #-40]
	ldp x14, x8, [x0, #16]
	umulh x6, x12, x23
	cinc x4, x4, hs
	mul x20, x11, x23
	stur x8, [x29, #-88]
	adds x17, x17, x7
	umulh x19, x11, x23
	cinc x6, x6, hs
	mul x22, x10, x23
	adds x18, x18, x20
	mul x24, x9, x23
	cinc x7, x19, hs
	adds x18, x18, x6
	ldur x19, [x29, #-32]
	umulh x21, x9, x23
	cinc x7, x7, hs
	stur x13, [x29, #-32]
	adds x2, x2, x22
	mul x26, x12, x14
	cinc x19, x19, hs
	adds x2, x2, x7
	cinc x7, x19, hs
	adds x4, x4, x24
	umulh x25, x12, x14
	cinc x19, x21, hs
	adds x4, x4, x7
	mul x28, x11, x14
	cinc x19, x19, hs
	adds x18, x18, x26
	umulh x27, x11, x14
	cinc x20, x25, hs
	mul x1, x10, x14
	adds x2, x2, x28
	umulh x30, x10, x14
	cinc x21, x27, hs
	adds x2, x2, x20
	mul x5, x9, x14
	cinc x20, x21, hs
	adds x1, x4, x1
	umulh x3, x9, x14
	cinc x21, x30, hs
	adds x1, x1, x20
	mul x15, x12, x16
	cinc x20, x21, hs
	adds x5, x19, x5
	umulh x6, x12, x8
	cinc x3, x3, hs
	adds x5, x5, x20
	mul x12, x12, x8
	cinc x3, x3, hs
	umulh x7, x11, x8
	mul x11, x11, x8
	adds x12, x2, x12
	umulh x4, x10, x8
	cinc x6, x6, hs
	mul x10, x10, x8
	adds x11, x1, x11
	mul x2, x9, x8
	cinc x1, x7, hs
	adds x11, x11, x6
	cinc x1, x1, hs
	umulh x9, x9, x8
	adds x10, x5, x10
	mov x5, #-4294967295
	cinc x4, x4, hs
	adds x10, x10, x1
	cinc x1, x4, hs
	mov x4, #-4294967295
	adds x2, x3, x2
	lsl x3, x15, #32
	umulh x4, x15, x4
	cinc x9, x9, hs
	adds x1, x2, x1
	lsr x2, x15, #32
	cinc x9, x9, hs
	adds x17, x3, x17
	sub x8, x15, x3
	adcs x18, x18, x2
	lsr x3, x17, #32
	mov x15, #-4294967295
	adcs x8, x8, xzr
	cinc x2, x4, hs
	adds x8, x8, x12
	cinc x12, x2, hs
	umulh x2, x17, x5
	adds x11, x11, x12
	adcs x10, x10, xzr
	adcs x12, x1, xzr
	lsl x1, x17, #32
	adcs x9, x9, xzr
	cset w4, hs
	adds x18, x1, x18
	sub x17, x17, x1
	adcs x8, x8, x3
	adcs x17, x17, xzr
	cinc x1, x2, hs
	adds x11, x17, x11
	lsr x2, x18, #32
	cinc x17, x1, hs
	umulh x1, x18, x15
	adds x10, x17, x10
	lsl x17, x18, #32
	adcs x12, x12, xzr
	adcs x9, x9, xzr
	cset w3, hs
	adds x8, x17, x8
	sub x17, x18, x17
	adcs x11, x11, x2
	lsr x5, x8, #32
	adcs x17, x17, xzr
	cinc x18, x1, hs
	adds x10, x17, x10
	umulh x1, x8, x15
	cinc x17, x18, hs
	lsl x18, x8, #32
	adds x2, x17, x12
	adds x11, x18, x11
	sub x8, x8, x18
	adcs x10, x10, x5
	adcs x8, x8, xzr
	cinc x18, x1, hs
	adds x5, x8, x2
	cinc x8, x18, hs
	cmn x17, x12
	ldur x12, [x29, #-16]
	adcs x9, x8, x9
	adc w17, w3, w4
	subs x8, x12, #1
	stur x8, [x29, #-16]
	ldr x8, [sp, #32]
	csel x8, x8, x13, eq
	adds x18, x11, #1
	mov x13, #-1
	adc x1, xzr, x13
	asr x2, x1, #63
	adds x1, x10, x1
	cinc x2, x2, hs
	adds x3, x1, x15
	adc x1, x2, x13
	mov w13, #-1
	asr x2, x1, #63
	adds x4, x1, x5
	cinc x1, x2, hs
	asr x2, x1, #63
	adds x1, x9, x1
	cinc w2, w2, hs
	adds x6, x1, x13
	adc w1, w2, w13
	ldr x13, [x8]
	mov w2, #-1
	bic w17, w1, w17
	str x13, [sp, #104]
	ldr x13, [x8, #8]
	tst x17, #0x1
	csel x11, x18, x11, eq
	csel x9, x6, x9, eq
	str x13, [sp, #112]
	ldp x13, x8, [x8, #16]
	stp x13, x8, [x29, #-104]
	csel x8, x3, x10, eq
	csel x10, x4, x5, eq
	stp x11, x8, [x0]
	ldr x8, [sp, #8]
	stp x10, x9, [x0, #16]
	ldr x9, [sp, #24]
	ldr x8, [x8, x12, lsl #3]
	ldp x10, x11, [x9]
	ldp x18, x17, [x8]
	stur x8, [x29, #-48]
	ldp x13, x8, [x8, #16]
	subs x12, x18, x10
	stp x18, x17, [x29, #-80]
	ngc x3, xzr
	subs x4, x17, x11
	stp x13, x8, [x29, #-64]
	ngc x5, xzr
	adds x4, x4, x3
	adc x3, x5, x3
	ldp x0, x5, [x9, #16]
	asr x6, x3, #63
	subs x9, x13, x0
	ngc x19, xzr
	adds x9, x3, x9
	adc x3, x6, x19
	subs x6, x8, x5
	asr x19, x3, #63
	ngc x20, xzr
	adds x3, x3, x6
	adc x6, x19, x20
	adds x19, x6, x12
	mov w12, w6
	adcs x12, x12, x4
	and x4, x6, #0xffffffff00000001
	adcs x9, x9, xzr
	adc x26, x4, x3
	adds x10, x10, x18
	adcs x3, x11, x17
	adcs x0, x0, x13
	adcs x4, x5, x8
	mov x8, #-1
	cset w11, hs
	adds x5, x10, #1
	adc x6, xzr, x8
	eor x11, x11, #0x1
	asr x20, x6, #63
	adds x6, x6, x3
	cinc x20, x20, hs
	adds x6, x6, x15
	adc x20, x20, x8
	asr x21, x20, #63
	adds x20, x20, x0
	cinc x21, x21, hs
	asr x22, x21, #63
	adds x21, x4, x21
	cinc x22, x22, hs
	adds x21, x21, x2
	adc x22, x22, x8
	tst x11, x22
	csel x11, x5, x10, eq
	csel x10, x6, x3, eq
	csel x28, x20, x0, eq
	csel x30, x21, x4, eq
	mul x8, x19, x11
	umulh x5, x19, x11
	umulh x3, x10, x19
	mul x6, x10, x19
	umulh x7, x28, x19
	mul x20, x28, x19
	adds x6, x5, x6
	umulh x22, x30, x19
	mul x15, x30, x19
	umulh x19, x28, x9
	mul x25, x12, x11
	umulh x18, x12, x11
	stp x19, x8, [sp, #88]
	umulh x19, x30, x9
	mul x24, x10, x12
	umulh x1, x9, x11
	str x19, [sp, #80]
	mul x13, x9, x11
	umulh x17, x10, x9
	mul x8, x10, x9
	mul x21, x28, x9
	mul x5, x30, x9
	cinc x9, x3, hs
	adds x9, x9, x20
	umulh x0, x10, x12
	cinc x3, x7, hs
	adds x3, x3, x15
	mul x27, x28, x12
	umulh x2, x28, x12
	umulh x4, x30, x12
	mul x12, x30, x12
	umulh x19, x26, x11
	mul x15, x26, x11
	cinc x11, x22, hs
	adds x25, x6, x25
	cinc x18, x18, hs
	adds x9, x9, x24
	ldr x6, [sp, #104]
	cinc x0, x0, hs
	adds x9, x9, x18
	umulh x20, x10, x26
	cinc x18, x0, hs
	adds x0, x3, x27
	mul x7, x10, x26
	cinc x10, x2, hs
	adds x18, x0, x18
	cinc x10, x10, hs
	adds x11, x11, x12
	cinc x12, x4, hs
	adds x10, x11, x10
	umulh x3, x28, x26
	cinc x12, x12, hs
	adds x22, x9, x13
	cinc x11, x1, hs
	adds x8, x18, x8
	mul x13, x28, x26
	cinc x9, x17, hs
	adds x11, x8, x11
	cinc x8, x9, hs
	adds x9, x10, x21
	ldr x10, [sp, #88]
	umulh x18, x30, x26
	cinc x10, x10, hs
	adds x8, x9, x8
	umulh x9, x6, x16
	mul x28, x6, x16
	str x9, [sp, #88]
	cinc x9, x10, hs
	adds x10, x12, x5
	ldr x12, [sp, #80]
	cinc x12, x12, hs
	adds x9, x10, x9
	mul x10, x30, x26
	cinc x12, x12, hs
	adds x11, x11, x15
	ldr x15, [sp, #112]
	cinc x0, x19, hs
	adds x8, x8, x7
	umulh x17, x15, x16
	cinc x1, x20, hs
	adds x8, x8, x0
	cinc x0, x1, hs
	adds x9, x9, x13
	ldur x20, [x29, #-96]
	cinc x1, x3, hs
	adds x0, x9, x0
	mov x9, #-4294967295
	str x17, [sp, #80]
	ldr x17, [sp, #96]
	umulh x2, x17, x9
	cinc x9, x1, hs
	adds x10, x12, x10
	cinc x18, x18, hs
	adds x1, x10, x9
	lsl x13, x17, #32
	mul x9, x15, x16
	lsr x12, x17, #32
	cinc x18, x18, hs
	adds x3, x13, x25
	sub x13, x17, x13
	mov x10, x15
	adcs x4, x22, x12
	lsl x5, x3, #32
	umulh x15, x6, x23
	umulh x22, x10, x23
	str x9, [sp, #72]
	ldur x9, [x29, #-104]
	mul x25, x10, x23
	umulh x12, x9, x16
	mul x17, x9, x16
	umulh x21, x9, x23
	str x12, [sp, #64]
	adcs x12, x13, xzr
	cinc x13, x2, hs
	adds x11, x12, x11
	cinc x12, x13, hs
	adds x13, x8, x12
	umulh x8, x20, x16
	adcs x2, x0, xzr
	adcs x1, x1, xzr
	mul x12, x20, x16
	lsr x16, x3, #32
	adcs x18, x18, xzr
	cset w30, hs
	adds x19, x5, x4
	mul x0, x6, x23
	stp x15, x8, [sp, #48]
	mov x8, #-4294967295
	adcs x16, x11, x16
	umulh x7, x3, x8
	sub x3, x3, x5
	lsr x15, x19, #32
	adcs x11, x3, xzr
	umulh x24, x19, x8
	mul x5, x20, x23
	cinc x3, x7, hs
	adds x13, x11, x13
	cinc x3, x3, hs
	mul x11, x9, x23
	adds x2, x3, x2
	lsl x3, x19, #32
	adcs x1, x1, xzr
	umulh x4, x20, x23
	mov x23, x20
	adcs x18, x18, xzr
	cset w7, hs
	adds x16, x3, x16
	sub x3, x19, x3
	adcs x13, x13, x15
	umulh x19, x16, x8
	lsl x26, x16, #32
	adcs x15, x3, xzr
	lsr x8, x16, #32
	cinc x24, x24, hs
	adds x2, x15, x2
	umulh x3, x6, x14
	cinc x24, x24, hs
	adds x27, x24, x1
	mul x15, x6, x14
	adds x13, x26, x13
	str x13, [sp, #96]
	sub x13, x16, x26
	adcs x16, x2, x8
	adcs x8, x13, xzr
	umulh x26, x10, x14
	mov x2, x10
	cinc x13, x19, hs
	adds x19, x8, x27
	cinc x8, x13, hs
	cmn x24, x1
	mul x27, x10, x14
	adcs x24, x8, x18
	ldr x8, [sp, #88]
	ldr x10, [sp, #72]
	adc w7, w7, w30
	umulh x13, x9, x14
	adds x8, x10, x8
	ldr x10, [sp, #80]
	mul x1, x9, x14
	cinc x18, x10, hs
	ldr x10, [sp, #64]
	adds x18, x18, x17
	umulh x30, x20, x14
	cinc x17, x10, hs
	ldr x10, [sp, #56]
	adds x17, x17, x12
	cinc x12, x10, hs
	mul x10, x20, x14
	ldr x14, [sp, #48]
	adds x8, x8, x0
	ldur x20, [x29, #-88]
	cinc x14, x14, hs
	adds x18, x18, x25
	cinc x0, x22, hs
	adds x14, x18, x14
	umulh x18, x6, x20
	cinc x0, x0, hs
	adds x11, x17, x11
	cinc x17, x21, hs
	adds x11, x11, x0
	mul x0, x6, x20
	cinc x17, x17, hs
	adds x12, x12, x5
	cinc x4, x4, hs
	adds x12, x12, x17
	umulh x17, x2, x20
	cinc x4, x4, hs
	adds x15, x14, x15
	cinc x14, x3, hs
	adds x11, x11, x27
	mul x3, x2, x20
	cinc x5, x26, hs
	adds x11, x11, x14
	cinc x14, x5, hs
	adds x12, x12, x1
	umulh x1, x9, x20
	cinc x13, x13, hs
	adds x12, x12, x14
	mov w5, #-1
	cinc x13, x13, hs
	adds x10, x4, x10
	mul x14, x9, x20
	cinc x4, x30, hs
	adds x10, x10, x13
	mov x30, #-4294967295
	cinc x13, x4, hs
	adds x11, x11, x0
	mul x0, x23, x20
	cinc x18, x18, hs
	adds x12, x12, x3
	mov x4, #-1
	cinc x17, x17, hs
	adds x12, x12, x18
	umulh x18, x23, x20
	cinc x17, x17, hs
	adds x10, x10, x14
	cinc x14, x1, hs
	adds x10, x10, x17
	mov x1, #-4294967295
	cinc x14, x14, hs
	adds x13, x13, x0
	lsl x17, x28, #32
	umulh x0, x28, x1
	cinc x18, x18, hs
	adds x13, x13, x14
	lsr x14, x28, #32
	cinc x18, x18, hs
	adds x8, x17, x8
	sub x17, x28, x17
	adcs x9, x15, x14
	ldr x15, [sp, #96]
	adcs x14, x17, xzr
	cinc x17, x0, hs
	adds x11, x14, x11
	lsr x0, x8, #32
	cinc x14, x17, hs
	umulh x17, x8, x1
	adds x12, x12, x14
	lsl x14, x8, #32
	adcs x10, x10, xzr
	adcs x13, x13, xzr
	sub x8, x8, x14
	adcs x18, x18, xzr
	cset w1, hs
	adds x9, x14, x9
	adcs x11, x11, x0
	adcs x8, x8, xzr
	cinc x14, x17, hs
	adds x8, x8, x12
	lsr x17, x9, #32
	cinc x12, x14, hs
	umulh x14, x9, x30
	adds x10, x12, x10
	lsl x12, x9, #32
	adcs x13, x13, xzr
	adcs x18, x18, xzr
	sub x9, x9, x12
	cset w0, hs
	adds x11, x12, x11
	adcs x8, x8, x17
	lsr x17, x11, #32
	adcs x9, x9, xzr
	cinc x12, x14, hs
	adds x9, x9, x10
	lsl x10, x11, #32
	umulh x14, x11, x30
	cinc x12, x12, hs
	adds x3, x12, x13
	adds x8, x10, x8
	sub x10, x11, x10
	adcs x9, x9, x17
	adcs x10, x10, xzr
	cinc x11, x14, hs
	adds x10, x10, x3
	mov w3, #-1
	cinc x11, x11, hs
	cmn x12, x13
	adcs x11, x11, x18
	mov x18, #-1
	adc w12, w0, w1
	adds x13, x8, #1
	adc x14, xzr, x18
	asr x17, x14, #63
	adds x14, x9, x14
	cinc x17, x17, hs
	adds x14, x14, x30
	adc x17, x17, x18
	asr x18, x17, #63
	adds x0, x17, x10
	cinc x17, x18, hs
	asr x18, x17, #63
	adds x17, x11, x17
	cinc w18, w18, hs
	adds x1, x17, x3
	adc w17, w18, w3
	bic w12, w17, w12
	tst x12, #0x1
	csel x18, x13, x8, eq
	csel x17, x14, x9, eq
	csel x14, x0, x10, eq
	csel x13, x1, x11, eq
	adds x8, x15, #1
	adc x9, xzr, x4
	asr x10, x9, #63
	adds x9, x16, x9
	cinc x10, x10, hs
	adds x9, x9, x30
	adc x10, x10, x4
	asr x11, x10, #63
	adds x10, x10, x19
	cinc x11, x11, hs
	asr x12, x11, #63
	adds x11, x24, x11
	cinc w12, w12, hs
	adds x11, x11, x5
	adc w12, w12, w5
	bic w12, w12, w7
	tst x12, #0x1
	csel x8, x8, x15, eq
	csel x9, x9, x16, eq
	mov x15, #-4294967296
	extr x12, x9, x8, #63
	csel x10, x10, x19, eq
	csel x11, x11, x24, eq
	extr x16, x10, x9, #63
	extr x2, x11, x10, #63
	adds x15, x12, x15
	adc x0, xzr, x4
	asr x1, x0, #63
	adds x0, x0, x16
	cinc x1, x1, hs
	asr x3, x1, #63
	adds x1, x2, x1
	cinc x3, x3, hs
	adds x1, x1, x5
	adc x3, x3, x4
	tst x11, #0x8000000000000000
	and x3, x3, #0x1
	csel x3, x3, xzr, eq
	orr x4, x3, x8, lsl #1
	cmp x3, #0
	mov x3, #-1
	csel x12, x15, x12, eq
	csel x16, x0, x16, eq
	csel x0, x1, x2, eq
	eor x15, x4, #0x1
	mov x1, #-1
	adds x8, x15, x8
	adcs x12, x12, x9
	adcs x10, x16, x10
	adcs x9, x0, x11
	cset w11, hs
	adds x15, x8, #1
	adc x16, xzr, x1
	eor x11, x11, #0x1
	asr x0, x16, #63
	adds x16, x16, x12
	cinc x0, x0, hs
	adds x16, x16, x30
	adc x0, x0, x1
	asr x1, x0, #63
	adds x0, x0, x10
	cinc x1, x1, hs
	asr x2, x1, #63
	adds x1, x9, x1
	cinc x2, x2, hs
	adds x1, x1, x5
	adc x2, x2, x3
	tst x11, x2
	csel x8, x15, x8, eq
	csel x12, x16, x12, eq
	csel x0, x0, x10, eq
	umulh x11, x18, x8
	csel x10, x1, x9, eq
	mul x3, x17, x8
	umulh x2, x17, x8
	mul x5, x14, x8
	adds x9, x11, x3
	umulh x4, x14, x8
	cinc x1, x2, hs
	mul x15, x18, x8
	adds x1, x1, x5
	umulh x6, x13, x8
	cinc x3, x4, hs
	mul x8, x13, x8
	mul x7, x12, x18
	umulh x16, x12, x18
	adds x8, x3, x8
	mul x20, x12, x17
	cinc x3, x6, hs
	adds x9, x9, x7
	umulh x19, x12, x17
	cinc x16, x16, hs
	mul x22, x12, x14
	adds x1, x1, x20
	umulh x21, x12, x14
	cinc x5, x19, hs
	adds x16, x1, x16
	umulh x23, x12, x13
	cinc x1, x5, hs
	adds x8, x8, x22
	mul x12, x12, x13
	cinc x5, x21, hs
	adds x8, x8, x1
	mul x25, x0, x18
	cinc x5, x5, hs
	umulh x24, x0, x18
	adds x12, x3, x12
	mul x27, x0, x17
	cinc x3, x23, hs
	adds x12, x12, x5
	cinc x3, x3, hs
	adds x16, x16, x25
	umulh x26, x0, x17
	cinc x5, x24, hs
	mul x11, x0, x14
	adds x8, x8, x27
	umulh x28, x0, x14
	cinc x7, x26, hs
	adds x8, x8, x5
	umulh x2, x0, x13
	cinc x5, x7, hs
	adds x11, x12, x11
	mul x0, x0, x13
	umulh x4, x10, x18
	mul x18, x10, x18
	umulh x1, x10, x17
	mul x17, x10, x17
	umulh x6, x10, x14
	mul x12, x10, x14
	cinc x14, x28, hs
	adds x11, x11, x5
	cinc x14, x14, hs
	adds x0, x3, x0
	cinc x2, x2, hs
	adds x14, x0, x14
	mul x3, x10, x13
	cinc x0, x2, hs
	adds x8, x8, x18
	umulh x10, x10, x13
	cinc x13, x4, hs
	adds x11, x11, x17
	cinc x17, x1, hs
	adds x11, x11, x13
	cinc x13, x17, hs
	adds x12, x14, x12
	lsl x14, x15, #32
	cinc x17, x6, hs
	adds x12, x12, x13
	umulh x13, x15, x30
	cinc x17, x17, hs
	adds x18, x0, x3
	lsr x0, x15, #32
	cinc x10, x10, hs
	adds x17, x18, x17
	cinc x10, x10, hs
	adds x9, x14, x9
	sub x14, x15, x14
	adcs x15, x16, x0
	adcs x14, x14, xzr
	cinc x13, x13, hs
	adds x8, x14, x8
	umulh x14, x9, x30
	cinc x13, x13, hs
	adds x11, x11, x13
	lsl x13, x9, #32
	adcs x12, x12, xzr
	adcs x16, x17, xzr
	lsr x17, x9, #32
	sub x9, x9, x13
	adcs x10, x10, xzr
	cset w18, hs
	adds x15, x13, x15
	adcs x8, x8, x17
	adcs x9, x9, xzr
	cinc x13, x14, hs
	adds x9, x9, x11
	cinc x11, x13, hs
	umulh x13, x15, x30
	adds x11, x11, x12
	lsl x12, x15, #32
	adcs x14, x16, xzr
	lsr x16, x15, #32
	adcs x10, x10, xzr
	cset w17, hs
	adds x8, x12, x8
	sub x12, x15, x12
	adcs x9, x9, x16
	umulh x15, x8, x30
	lsr x16, x8, #32
	adcs x12, x12, xzr
	cinc x13, x13, hs
	adds x11, x12, x11
	lsl x12, x8, #32
	cinc x13, x13, hs
	adds x0, x13, x14
	sub x8, x8, x12
	adds x9, x12, x9
	adcs x11, x11, x16
	mov x16, #-1
	adcs x8, x8, xzr
	cinc x12, x15, hs
	adds x8, x8, x0
	cinc x12, x12, hs
	cmn x13, x14
	adcs x10, x12, x10
	adc w12, w17, w18
	adds x13, x9, #1
	mov w18, #-1
	adc x14, xzr, x16
	asr x15, x14, #63
	adds x14, x11, x14
	cinc x15, x15, hs
	adds x14, x14, x30
	adc x15, x15, x16
	asr x16, x15, #63
	adds x15, x15, x8
	cinc x16, x16, hs
	asr x17, x16, #63
	adds x16, x10, x16
	cinc w17, w17, hs
	adds x0, x16, x18
	adc w16, w17, w18
	bic w12, w16, w12
	tst x12, #0x1
	csel x1, x13, x9, eq
	csel x16, x14, x11, eq
	csel x14, x15, x8, eq
	umulh x12, x1, x1
	csel x13, x0, x10, eq
	mul x3, x1, x16
	umulh x2, x1, x16
	mul x18, x1, x14
	adds x10, x12, x3
	umulh x17, x1, x14
	cinc x12, x2, hs
	mul x11, x1, x13
	adds x6, x12, x18
	umulh x9, x1, x13
	cinc x15, x17, hs
	mul x0, x16, x16
	adds x7, x15, x11
	umulh x8, x16, x16
	cinc x19, x9, hs
	adds x10, x10, x3
	mul x5, x16, x14
	cinc x2, x2, hs
	adds x0, x6, x0
	umulh x4, x16, x14
	cinc x8, x8, hs
	adds x0, x0, x2
	mul x15, x16, x13
	cinc x8, x8, hs
	adds x2, x7, x5
	umulh x12, x16, x13
	cinc x7, x4, hs
	adds x8, x2, x8
	mul x6, x14, x14
	cinc x2, x7, hs
	adds x7, x19, x15
	umulh x3, x14, x14
	cinc x20, x12, hs
	adds x2, x7, x2
	mul x19, x14, x13
	cinc x7, x20, hs
	adds x18, x0, x18
	cinc x17, x17, hs
	adds x8, x8, x5
	umulh x0, x14, x13
	cinc x4, x4, hs
	adds x8, x8, x17
	cinc x17, x4, hs
	adds x2, x2, x6
	cinc x3, x3, hs
	adds x17, x2, x17
	mul x4, x1, x1
	cinc x2, x3, hs
	adds x3, x7, x19
	mul x5, x13, x13
	cinc x6, x0, hs
	adds x2, x3, x2
	cinc x3, x6, hs
	adds x8, x8, x11
	umulh x11, x13, x13
	cinc x9, x9, hs
	adds x15, x17, x15
	lsl x17, x4, #32
	cinc x12, x12, hs
	adds x9, x15, x9
	cinc x12, x12, hs
	adds x15, x2, x19
	cinc x0, x0, hs
	adds x12, x15, x12
	umulh x15, x4, x30
	cinc x0, x0, hs
	adds x2, x3, x5
	lsr x3, x4, #32
	cinc x11, x11, hs
	adds x0, x2, x0
	cinc x11, x11, hs
	adds x10, x17, x10
	sub x17, x4, x17
	adcs x18, x18, x3
	lsr x2, x10, #32
	adcs x17, x17, xzr
	cinc x15, x15, hs
	adds x8, x17, x8
	umulh x17, x10, x30
	cinc x15, x15, hs
	adds x9, x9, x15
	lsl x15, x10, #32
	adcs x12, x12, xzr
	adcs x0, x0, xzr
	sub x10, x10, x15
	adcs x11, x11, xzr
	cset w3, hs
	adds x18, x15, x18
	adcs x8, x8, x2
	adcs x10, x10, xzr
	cinc x15, x17, hs
	adds x9, x10, x9
	cinc x10, x15, hs
	umulh x15, x18, x30
	adds x10, x10, x12
	lsl x12, x18, #32
	adcs x17, x0, xzr
	lsr x0, x18, #32
	adcs x11, x11, xzr
	cset w2, hs
	adds x8, x12, x8
	sub x12, x18, x12
	adcs x9, x9, x0
	umulh x18, x8, x30
	lsr x0, x8, #32
	adcs x12, x12, xzr
	cinc x15, x15, hs
	adds x10, x12, x10
	lsl x12, x8, #32
	cinc x15, x15, hs
	adds x4, x15, x17
	sub x8, x8, x12
	adds x9, x12, x9
	adcs x10, x10, x0
	mov x0, #-1
	adcs x8, x8, xzr
	cinc x12, x18, hs
	adds x8, x8, x4
	cinc x12, x12, hs
	cmn x15, x17
	adcs x11, x12, x11
	adc w12, w2, w3
	adds x15, x9, #1
	mov w3, #-1
	adc x17, xzr, x0
	asr x18, x17, #63
	adds x17, x10, x17
	cinc x18, x18, hs
	adds x17, x17, x30
	adc x18, x18, x0
	asr x0, x18, #63
	adds x18, x18, x8
	cinc x0, x0, hs
	asr x2, x0, #63
	adds x0, x11, x0
	cinc w2, w2, hs
	adds x0, x0, x3
	adc w2, w2, w3
	ldp x4, x3, [x29, #-64]
	bic w12, w2, w12
	tst x12, #0x1
	csel x11, x0, x11, eq
	ldp x0, x2, [x29, #-80]
	csel x9, x15, x9, eq
	csel x10, x17, x10, eq
	csel x8, x18, x8, eq
	subs x9, x9, x0
	ngc x12, xzr
	subs x10, x10, x2
	ngc x15, xzr
	adds x10, x12, x10
	adc x12, x12, x15
	subs x8, x8, x4
	asr x15, x12, #63
	ngc x17, xzr
	adds x8, x12, x8
	adc x12, x15, x17
	subs x11, x11, x3
	asr x15, x12, #63
	ngc x17, xzr
	adds x11, x12, x11
	adc x12, x15, x17
	adds x9, x12, x9
	mov w15, w12
	and x12, x12, #0xffffffff00000001
	adcs x10, x15, x10
	adcs x8, x8, xzr
	adc x11, x12, x11
	subs x9, x9, x0
	ngc x12, xzr
	subs x10, x10, x2
	ngc x15, xzr
	adds x10, x10, x12
	adc x12, x15, x12
	subs x8, x8, x4
	asr x15, x12, #63
	ngc x17, xzr
	adds x12, x12, x8
	adc x8, x15, x17
	subs x11, x11, x3
	asr x15, x8, #63
	ngc x17, xzr
	adds x11, x11, x8
	adc x17, x17, x15
	adds x15, x17, x9
	mov w8, w17
	and x9, x17, #0xffffffff00000001
	adcs x8, x8, x10
	adcs x18, x12, xzr
	adc x17, x9, x11
	subs x9, x0, x15
	ngc x10, xzr
	subs x11, x2, x8
	ngc x12, xzr
	adds x0, x11, x10
	adc x10, x12, x10
	subs x11, x4, x18
	asr x12, x10, #63
	ngc x2, xzr
	adds x10, x10, x11
	adc x11, x12, x2
	subs x12, x3, x17
	asr x2, x11, #63
	ngc x3, xzr
	adds x11, x12, x11
	adc x12, x3, x2
	adds x7, x12, x9
	mov w2, w12
	and x12, x12, #0xffffffff00000001
	adcs x2, x2, x0
	umulh x3, x7, x1
	mul x23, x2, x1
	umulh x22, x2, x1
	umulh x24, x2, x16
	mul x25, x2, x16
	umulh x26, x2, x14
	mul x27, x2, x14
	umulh x0, x2, x13
	mul x21, x2, x13
	adcs x2, x10, xzr
	adc x11, x12, x11
	adds x3, x23, x3
	mul x28, x2, x1
	cinc x22, x22, hs
	umulh x10, x2, x1
	mul x9, x7, x1
	adds x22, x28, x22
	umulh x12, x1, x11
	cinc x10, x10, hs
	mul x1, x1, x11
	mul x5, x7, x16
	umulh x6, x7, x16
	adds x1, x1, x10
	mul x23, x2, x16
	cinc x12, x12, hs
	adds x10, x3, x5
	umulh x30, x2, x16
	cinc x3, x6, hs
	adds x5, x22, x25
	umulh x28, x16, x11
	cinc x6, x24, hs
	adds x3, x5, x3
	cinc x6, x6, hs
	adds x1, x1, x23
	mul x16, x16, x11
	cinc x22, x30, hs
	adds x1, x1, x6
	mul x20, x7, x14
	cinc x22, x22, hs
	umulh x19, x7, x14
	adds x12, x12, x16
	mul x6, x2, x14
	cinc x16, x28, hs
	adds x12, x12, x22
	cinc x16, x16, hs
	adds x3, x3, x20
	umulh x5, x2, x14
	cinc x19, x19, hs
	umulh x22, x14, x11
	adds x19, x19, x27
	cinc x20, x26, hs
	adds x1, x19, x1
	mul x14, x14, x11
	cinc x20, x20, hs
	adds x12, x12, x6
	cinc x5, x5, hs
	adds x12, x12, x20
	mov x6, #-4294967295
	umulh x4, x7, x13
	cinc x5, x5, hs
	mul x7, x7, x13
	adds x14, x16, x14
	umulh x19, x2, x13
	cinc x16, x22, hs
	adds x14, x14, x5
	cinc x16, x16, hs
	mul x2, x2, x13
	adds x1, x1, x7
	mul x5, x13, x11
	umulh x11, x13, x11
	cinc x13, x4, hs
	adds x12, x12, x21
	cinc x0, x0, hs
	adds x12, x12, x13
	lsr x4, x9, #32
	cinc x13, x0, hs
	adds x14, x14, x2
	lsl x0, x9, #32
	cinc x2, x19, hs
	adds x13, x14, x13
	umulh x14, x9, x6
	cinc x2, x2, hs
	adds x16, x16, x5
	sub x9, x9, x0
	cinc x11, x11, hs
	adds x16, x16, x2
	mov x5, #-4294967295
	cinc x11, x11, hs
	adds x10, x0, x10
	adcs x0, x3, x4
	lsr x2, x10, #32
	adcs x9, x9, xzr
	cinc x14, x14, hs
	adds x9, x9, x1
	umulh x1, x10, x6
	cinc x14, x14, hs
	adds x12, x12, x14
	lsl x14, x10, #32
	adcs x13, x13, xzr
	adcs x16, x16, xzr
	sub x10, x10, x14
	adcs x11, x11, xzr
	cset w3, hs
	adds x0, x14, x0
	adcs x9, x9, x2
	adcs x10, x10, xzr
	cinc x14, x1, hs
	adds x10, x10, x12
	lsr x1, x0, #32
	cinc x12, x14, hs
	umulh x14, x0, x5
	adds x12, x12, x13
	lsl x13, x0, #32
	adcs x16, x16, xzr
	adcs x11, x11, xzr
	cset w2, hs
	adds x9, x13, x9
	sub x13, x0, x13
	adcs x10, x10, x1
	umulh x0, x9, x5
	lsr x1, x9, #32
	adcs x13, x13, xzr
	cinc x14, x14, hs
	adds x12, x13, x12
	lsl x13, x9, #32
	cinc x14, x14, hs
	adds x4, x14, x16
	sub x9, x9, x13
	adds x10, x13, x10
	adcs x12, x12, x1
	mov x1, #-1
	adcs x9, x9, xzr
	cinc x13, x0, hs
	adds x9, x9, x4
	cinc x13, x13, hs
	cmn x14, x16
	adcs x11, x13, x11
	adc w13, w2, w3
	adds x14, x10, #1
	mov w3, #-1
	adc x16, xzr, x1
	asr x0, x16, #63
	adds x16, x12, x16
	cinc x0, x0, hs
	adds x16, x16, x5
	adc x0, x0, x1
	asr x1, x0, #63
	adds x0, x0, x9
	cinc x1, x1, hs
	asr x2, x1, #63
	adds x1, x11, x1
	cinc w2, w2, hs
	adds x1, x1, x3
	adc w2, w2, w3
	bic w13, w2, w13
	ldur x2, [x29, #-48]
	tst x13, #0x1
	csel x10, x14, x10, eq
	ldp x13, x14, [x2, #32]
	csel x12, x16, x12, eq
	csel x9, x0, x9, eq
	csel x11, x1, x11, eq
	stp xzr, xzr, [x2, #64]
	subs x10, x10, x13
	stp xzr, xzr, [x2, #96]
	ngc x13, xzr
	subs x12, x12, x14
	stp xzr, xzr, [x2, #80]
	ngc x16, xzr
	adds x12, x13, x12
	stp xzr, xzr, [x2, #112]
	adc x13, x13, x16
	ldp x14, x16, [x2, #48]
	stp x15, x8, [x2]
	stp x18, x17, [x2, #16]
	subs x9, x9, x14
	asr x14, x13, #63
	str xzr, [x2, #128]
	ngc x0, xzr
	adds x9, x13, x9
	adc x13, x14, x0
	subs x11, x11, x16
	mov x0, x2
	asr x14, x13, #63
	ngc x16, xzr
	adds x11, x13, x11
	adc x13, x14, x16
	adds x10, x13, x10
	mov w14, w13
	adcs x12, x14, x12
	stp x10, x12, [x2, #32]
	and x10, x13, #0xffffffff00000001
	adcs x9, x9, xzr
	adc x10, x10, x11
	stp x9, x10, [x2, #48]
	ldr x9, [x2, #136]
	ldr x10, [sp, #40]
	add x9, x9, #1
	cmp x9, x10
	stp x9, xzr, [x2, #136]
	b.hs .LBB7_9
	cbz x8, .LBB7_9
	orr x8, x18, x15
	orr x8, x8, x17
	cbz x8, .LBB7_9
	tst x9, #0x3f
	b.ne .LBB7_7
	stp x15, x18, [x0, #112]
	str x17, [x0, #128]
.LBB7_7:
	mov x9, xzr
	stp xzr, xzr, [x0, #96]
	stp x15, x18, [x0, #64]
	stp x17, xzr, [x0, #80]
.LBB7_8:
	ldur x8, [x29, #-24]
	ldp x14, x13, [x29, #-40]
	add x8, x9, x8
	stur x8, [x29, #-24]
	ldur x8, [x29, #-16]
	sub x14, x14, #32
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
	ldp x20, x19, [sp, #304]
	ldp x22, x21, [sp, #288]
	ldp x24, x23, [sp, #272]
	ldp x26, x25, [sp, #256]
	ldp x28, x27, [sp, #240]
	ldp x29, x30, [sp, #224]
	add sp, sp, #320
	ret
.Lfunc_end7:

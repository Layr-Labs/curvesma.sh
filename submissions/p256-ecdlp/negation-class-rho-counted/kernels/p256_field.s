// clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S p256_walk.c
// Apple clang 21.0.0 (LLVM). Comments and assembler directives stripped; functions: fe_mul, fe_add, fe_sub, fe_inv
fe_mul:
	stp x28, x27, [sp, #-80]!
	ldp x13, x11, [x2]
	stp x24, x23, [sp, #32]
	ldp x12, x14, [x1]
	stp x22, x21, [sp, #48]
	stp x20, x19, [sp, #64]
	ldp x10, x9, [x2, #16]
	stp x26, x25, [sp, #16]
	umulh x15, x13, x12
	mul x17, x11, x12
	umulh x16, x11, x12
	mul x2, x10, x12
	umulh x18, x10, x12
	mul x4, x9, x12
	umulh x3, x9, x12
	umulh x5, x13, x14
	mul x6, x13, x14
	umulh x7, x11, x14
	mul x19, x11, x14
	umulh x20, x10, x14
	mul x21, x10, x14
	umulh x22, x9, x14
	mul x23, x9, x14
	adds x14, x17, x15
	cinc x15, x16, hs
	adds x15, x15, x2
	mul x8, x13, x12
	cinc x17, x18, hs
	ldp x24, x12, [x1, #16]
	adds x17, x17, x4
	cinc x2, x3, hs
	adds x14, x14, x6
	cinc x4, x5, hs
	adds x15, x15, x19
	mul x25, x13, x24
	cinc x5, x7, hs
	adds x15, x15, x4
	umulh x1, x13, x24
	cinc x5, x5, hs
	adds x17, x17, x21
	cinc x6, x20, hs
	adds x17, x17, x5
	mul x27, x11, x24
	cinc x5, x6, hs
	adds x2, x2, x23
	cinc x6, x22, hs
	adds x2, x2, x5
	umulh x26, x11, x24
	cinc x6, x6, hs
	adds x15, x15, x25
	cinc x1, x1, hs
	ldp x20, x19, [sp, #64]
	mul x16, x10, x24
	ldp x22, x21, [sp, #48]
	adds x17, x17, x27
	umulh x28, x10, x24
	cinc x7, x26, hs
	adds x17, x17, x1
	mul x3, x9, x24
	cinc x1, x7, hs
	adds x16, x2, x16
	ldp x26, x25, [sp, #16]
	umulh x18, x9, x24
	ldp x24, x23, [sp, #32]
	cinc x7, x28, hs
	adds x16, x16, x1
	umulh x4, x13, x12
	cinc x1, x7, hs
	adds x3, x6, x3
	mul x13, x13, x12
	cinc x18, x18, hs
	adds x1, x3, x1
	umulh x5, x11, x12
	cinc x18, x18, hs
	mul x11, x11, x12
	adds x13, x17, x13
	umulh x2, x10, x12
	cinc x3, x4, hs
	mul x10, x10, x12
	adds x11, x16, x11
	mul x17, x9, x12
	cinc x16, x5, hs
	adds x11, x11, x3
	umulh x12, x9, x12
	cinc x9, x16, hs
	adds x10, x1, x10
	cinc x16, x2, hs
	adds x10, x10, x9
	mov x9, #-4294967295
	cinc x16, x16, hs
	adds x17, x18, x17
	lsl x18, x8, #32
	umulh x1, x8, x9
	cinc x12, x12, hs
	adds x16, x17, x16
	lsr x17, x8, #32
	cinc x12, x12, hs
	adds x14, x18, x14
	sub x8, x8, x18
	adcs x15, x15, x17
	lsr x18, x14, #32
	adcs x8, x8, xzr
	cinc x17, x1, hs
	adds x8, x8, x13
	cinc x13, x17, hs
	umulh x17, x14, x9
	adds x11, x11, x13
	adcs x10, x10, xzr
	adcs x13, x16, xzr
	lsl x16, x14, #32
	adcs x12, x12, xzr
	cset w1, hs
	adds x15, x16, x15
	sub x14, x14, x16
	adcs x8, x8, x18
	adcs x14, x14, xzr
	cinc x16, x17, hs
	adds x11, x14, x11
	lsr x17, x15, #32
	cinc x14, x16, hs
	umulh x16, x15, x9
	adds x10, x14, x10
	lsl x14, x15, #32
	adcs x13, x13, xzr
	adcs x12, x12, xzr
	cset w18, hs
	adds x8, x14, x8
	sub x14, x15, x14
	adcs x11, x11, x17
	lsr x2, x8, #32
	adcs x14, x14, xzr
	cinc x15, x16, hs
	adds x10, x14, x10
	umulh x16, x8, x9
	cinc x14, x15, hs
	lsl x15, x8, #32
	adds x17, x14, x13
	adds x11, x15, x11
	sub x8, x8, x15
	adcs x10, x10, x2
	adcs x8, x8, xzr
	cinc x15, x16, hs
	adds x8, x8, x17
	cinc x15, x15, hs
	cmn x14, x13
	adcs x12, x15, x12
	mov x15, #-1
	adc w13, w18, w1
	adds x14, x11, #1
	mov w18, #-1
	adc x16, xzr, x15
	asr x17, x16, #63
	adds x16, x10, x16
	cinc x17, x17, hs
	adds x9, x16, x9
	adc x15, x17, x15
	asr x16, x15, #63
	adds x15, x15, x8
	cinc x16, x16, hs
	asr x17, x16, #63
	adds x16, x12, x16
	cinc w17, w17, hs
	adds x16, x16, x18
	adc w17, w17, w18
	bic w13, w17, w13
	tst x13, #0x1
	csel x11, x14, x11, eq
	csel x9, x9, x10, eq
	csel x8, x15, x8, eq
	stp x11, x9, [x0]
	csel x9, x16, x12, eq
	stp x8, x9, [x0, #16]
	ldp x28, x27, [sp], #80
	ret
.Lfunc_end0:
fe_add:
	ldp x8, x11, [x2]
	mov x17, #-4294967295
	ldp x9, x10, [x1]
	ldp x12, x14, [x2, #16]
	adds x8, x8, x9
	ldp x9, x13, [x1, #16]
	adcs x10, x11, x10
	mov w1, #-1
	adcs x9, x12, x9
	mov x12, #-1
	adcs x11, x14, x13
	cset w13, hs
	adds x14, x8, #1
	adc x15, xzr, x12
	eor x13, x13, #0x1
	asr x16, x15, #63
	adds x15, x15, x10
	cinc x16, x16, hs
	adds x15, x15, x17
	adc x16, x16, x12
	asr x17, x16, #63
	adds x16, x16, x9
	cinc x17, x17, hs
	asr x18, x17, #63
	adds x17, x11, x17
	cinc x18, x18, hs
	adds x17, x17, x1
	adc x12, x18, x12
	tst x13, x12
	csel x8, x14, x8, eq
	csel x10, x15, x10, eq
	stp x8, x10, [x0]
	csel x8, x16, x9, eq
	csel x9, x17, x11, eq
	stp x8, x9, [x0, #16]
	ret
.Lfunc_end1:
fe_sub:
	ldp x8, x11, [x2]
	ldp x9, x10, [x1]
	ldp x13, x14, [x2, #16]
	subs x8, x9, x8
	ngc x9, xzr
	subs x10, x10, x11
	ngc x11, xzr
	adds x10, x10, x9
	adc x9, x11, x9
	ldp x12, x11, [x1, #16]
	subs x12, x12, x13
	asr x13, x9, #63
	ngc x15, xzr
	adds x9, x9, x12
	adc x12, x13, x15
	subs x11, x11, x14
	asr x13, x12, #63
	ngc x14, xzr
	adds x11, x12, x11
	adc x12, x13, x14
	adds x8, x12, x8
	mov w13, w12
	adcs x10, x13, x10
	stp x8, x10, [x0]
	adcs x8, x9, xzr
	and x9, x12, #0xffffffff00000001
	adc x9, x9, x11
	stp x8, x9, [x0, #16]
	ret
.Lfunc_end2:
fe_inv:
	str x30, [sp, #-96]!
	ldp x10, x11, [x1, #16]
	mov w8, #254
	ldp x13, x14, [x1]
	mov x9, #-4294967295
	mov x12, #-1
	mov w30, #-1
	stp x28, x27, [sp, #16]
	mov x17, x11
	mov x18, x10
	stp x26, x25, [sp, #32]
	mov x3, x14
	mov x1, x13
	stp x24, x23, [sp, #48]
	stp x22, x21, [sp, #64]
	stp x20, x19, [sp, #80]
	str x0, [sp, #8]
	b .LBB3_2
.LBB3_1:
	sub w8, w8, #1
	cmn w8, #1
	b.eq .LBB3_4
.LBB3_2:
	umulh x4, x1, x1
	adrp x15, fe_inv.E
	add x15, x15, :lo12:fe_inv.E
	mul x21, x3, x1
	umulh x20, x3, x1
	mul x19, x18, x1
	adds x4, x21, x4
	umulh x7, x18, x1
	cinc x6, x20, hs
	mul x5, x17, x1
	adds x26, x6, x19
	umulh x2, x17, x1
	cinc x27, x7, hs
	mul x23, x3, x3
	adds x27, x27, x5
	umulh x22, x3, x3
	cinc x28, x2, hs
	adds x4, x4, x21
	mul x25, x18, x3
	cinc x20, x20, hs
	adds x21, x26, x23
	umulh x24, x18, x3
	cinc x22, x22, hs
	adds x20, x21, x20
	umulh x6, x17, x3
	cinc x22, x22, hs
	adds x23, x27, x25
	mul x3, x17, x3
	cinc x26, x24, hs
	adds x22, x23, x22
	mul x23, x18, x18
	cinc x26, x26, hs
	umulh x21, x18, x18
	adds x27, x28, x3
	cinc x28, x6, hs
	adds x26, x27, x26
	mul x27, x17, x18
	cinc x28, x28, hs
	adds x19, x20, x19
	cinc x7, x7, hs
	adds x20, x22, x25
	umulh x18, x17, x18
	cinc x22, x24, hs
	adds x7, x20, x7
	cinc x20, x22, hs
	adds x22, x26, x23
	mul x1, x1, x1
	cinc x21, x21, hs
	adds x20, x22, x20
	cinc x21, x21, hs
	adds x22, x28, x27
	cinc x23, x18, hs
	adds x21, x22, x21
	cinc x22, x23, hs
	adds x5, x7, x5
	mul x7, x17, x17
	cinc x2, x2, hs
	adds x3, x20, x3
	umulh x17, x17, x17
	cinc x6, x6, hs
	adds x2, x3, x2
	cinc x3, x6, hs
	adds x6, x21, x27
	cinc x18, x18, hs
	adds x3, x6, x3
	umulh x20, x1, x9
	cinc x18, x18, hs
	adds x6, x22, x7
	lsl x7, x1, #32
	cinc x17, x17, hs
	adds x18, x6, x18
	lsr x6, x1, #32
	cinc x17, x17, hs
	adds x4, x7, x4
	sub x1, x1, x7
	adcs x6, x19, x6
	lsr x19, x4, #32
	adcs x1, x1, xzr
	cinc x7, x20, hs
	adds x1, x1, x5
	cinc x5, x7, hs
	umulh x7, x4, x9
	adds x2, x2, x5
	lsl x5, x4, #32
	adcs x3, x3, xzr
	adcs x18, x18, xzr
	sub x4, x4, x5
	adcs x17, x17, xzr
	cset w20, hs
	adds x6, x5, x6
	adcs x1, x1, x19
	adcs x4, x4, xzr
	cinc x5, x7, hs
	adds x2, x4, x2
	lsr x7, x6, #32
	cinc x4, x5, hs
	umulh x5, x6, x9
	adds x3, x4, x3
	lsl x4, x6, #32
	adcs x18, x18, xzr
	adcs x17, x17, xzr
	cset w19, hs
	adds x1, x4, x1
	sub x4, x6, x4
	adcs x2, x2, x7
	umulh x6, x1, x9
	lsr x7, x1, #32
	adcs x4, x4, xzr
	cinc x5, x5, hs
	adds x3, x4, x3
	lsl x4, x1, #32
	cinc x5, x5, hs
	adds x21, x5, x18
	sub x1, x1, x4
	adds x2, x4, x2
	adcs x3, x3, x7
	adcs x1, x1, xzr
	cinc x4, x6, hs
	adds x6, x1, x21
	cinc x1, x4, hs
	cmn x5, x18
	cinc x18, x17, hs
	adcs xzr, x1, x17
	adc w17, w19, w20
	adds x4, x1, x18
	lsr w20, w8, #6
	adds x18, x2, #1
	adc x1, xzr, x12
	asr x5, x1, #63
	adds x1, x3, x1
	cinc x5, x5, hs
	adds x7, x1, x9
	adc x1, x5, x12
	asr x5, x1, #63
	adds x19, x1, x6
	cinc x1, x5, hs
	asr x5, x1, #63
	adds x1, x4, x1
	cinc x5, x5, hs
	adds x21, x1, x30
	ldr x1, [x15, w20, uxtw #3]
	adc x5, x5, x12
	bic w17, w5, w17
	lsr x5, x1, x8
	tst x17, #0x1
	csel x1, x18, x2, eq
	csel x3, x7, x3, eq
	csel x18, x19, x6, eq
	csel x17, x21, x4, eq
	tbz w5, #0, .LBB3_1
	umulh x4, x1, x13
	mul x6, x1, x14
	umulh x5, x1, x14
	mul x19, x1, x10
	adds x4, x4, x6
	umulh x7, x1, x10
	cinc x5, x5, hs
	mul x2, x1, x13
	adds x5, x5, x19
	umulh x20, x1, x11
	cinc x7, x7, hs
	mul x1, x1, x11
	mul x22, x3, x13
	umulh x21, x3, x13
	adds x7, x7, x1
	mul x24, x3, x14
	cinc x20, x20, hs
	adds x1, x4, x22
	umulh x23, x3, x14
	cinc x4, x21, hs
	mul x26, x3, x10
	adds x5, x5, x24
	umulh x25, x3, x10
	cinc x21, x23, hs
	adds x4, x5, x4
	umulh x27, x3, x11
	cinc x21, x21, hs
	adds x7, x7, x26
	mul x3, x3, x11
	cinc x22, x25, hs
	adds x7, x7, x21
	mul x30, x18, x13
	cinc x22, x22, hs
	umulh x28, x18, x13
	adds x3, x20, x3
	mul x0, x18, x14
	cinc x20, x27, hs
	adds x3, x3, x22
	cinc x20, x20, hs
	adds x4, x4, x30
	mov w30, #-1
	umulh x16, x18, x14
	cinc x23, x28, hs
	mul x6, x18, x10
	adds x0, x7, x0
	umulh x15, x18, x10
	cinc x16, x16, hs
	adds x0, x0, x23
	umulh x19, x18, x11
	cinc x16, x16, hs
	adds x3, x3, x6
	mul x18, x18, x11
	cinc x15, x15, hs
	adds x16, x3, x16
	mul x21, x17, x13
	cinc x15, x15, hs
	umulh x5, x17, x13
	adds x18, x20, x18
	mul x7, x17, x14
	cinc x19, x19, hs
	adds x15, x18, x15
	cinc x18, x19, hs
	adds x0, x0, x21
	umulh x22, x17, x14
	cinc x5, x5, hs
	mul x3, x17, x10
	adds x16, x16, x7
	umulh x6, x17, x10
	cinc x7, x22, hs
	adds x16, x16, x5
	mul x19, x17, x11
	cinc x5, x7, hs
	adds x15, x15, x3
	umulh x17, x17, x11
	cinc x3, x6, hs
	adds x15, x15, x5
	lsl x5, x2, #32
	cinc x3, x3, hs
	umulh x6, x2, x9
	adds x18, x18, x19
	cinc x17, x17, hs
	adds x18, x18, x3
	lsr x3, x2, #32
	cinc x17, x17, hs
	adds x1, x5, x1
	sub x2, x2, x5
	adcs x3, x4, x3
	lsr x5, x1, #32
	adcs x2, x2, xzr
	cinc x4, x6, hs
	adds x0, x2, x0
	cinc x2, x4, hs
	umulh x4, x1, x9
	adds x16, x16, x2
	lsl x2, x1, #32
	adcs x15, x15, xzr
	adcs x18, x18, xzr
	sub x1, x1, x2
	adcs x17, x17, xzr
	cset w6, hs
	adds x3, x2, x3
	adcs x0, x0, x5
	adcs x1, x1, xzr
	cinc x2, x4, hs
	adds x16, x1, x16
	lsr x4, x3, #32
	cinc x1, x2, hs
	umulh x2, x3, x9
	adds x15, x1, x15
	lsl x1, x3, #32
	adcs x18, x18, xzr
	adcs x17, x17, xzr
	cset w5, hs
	adds x0, x1, x0
	sub x1, x3, x1
	adcs x16, x16, x4
	umulh x3, x0, x9
	lsr x7, x0, #32
	adcs x1, x1, xzr
	cinc x2, x2, hs
	adds x15, x1, x15
	cinc x1, x2, hs
	lsl x2, x0, #32
	adds x4, x1, x18
	adds x16, x2, x16
	sub x0, x0, x2
	adcs x15, x15, x7
	adcs x0, x0, xzr
	cinc x2, x3, hs
	adds x0, x0, x4
	cinc x2, x2, hs
	cmn x1, x18
	adcs x17, x2, x17
	adc w18, w5, w6
	adds x1, x16, #1
	adc x2, xzr, x12
	asr x3, x2, #63
	adds x2, x15, x2
	cinc x3, x3, hs
	adds x2, x2, x9
	adc x3, x3, x12
	asr x4, x3, #63
	adds x5, x3, x0
	cinc x3, x4, hs
	asr x4, x3, #63
	adds x3, x17, x3
	cinc w4, w4, hs
	adds x6, x3, x30
	adc w3, w4, w30
	bic w18, w3, w18
	tst x18, #0x1
	csel x1, x1, x16, eq
	csel x3, x2, x15, eq
	csel x18, x5, x0, eq
	csel x17, x6, x17, eq
	b .LBB3_1
.LBB3_4:
	ldr x8, [sp, #8]
	ldp x20, x19, [sp, #80]
	ldp x22, x21, [sp, #64]
	ldp x24, x23, [sp, #48]
	stp x1, x3, [x8]
	stp x18, x17, [x8, #16]
	ldp x26, x25, [sp, #32]
	ldp x28, x27, [sp, #16]
	ldr x30, [sp], #96
	ret
.Lfunc_end3:

// clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S p256_walk.c
// Apple clang 21.0.0 (LLVM). Comments and assembler directives stripped; functions: fe_mul, fe_add, fe_sub, fe_inv
fe_mul:
	stp x26, x25, [sp, #-64]!
	ldp x13, x10, [x2]
	stp x22, x21, [sp, #32]
	ldp x11, x15, [x1]
	stp x20, x19, [sp, #48]
	stp x24, x23, [sp, #16]
	ldp x8, x9, [x2, #16]
	ldp x14, x12, [x1, #16]
	umulh x16, x13, x11
	mul x18, x10, x11
	umulh x17, x10, x11
	mul x3, x13, x15
	adds x16, x16, x18
	umulh x4, x13, x15
	cinc x17, x17, hs
	mul x2, x8, x11
	umulh x5, x8, x11
	mul x6, x10, x15
	mul x1, x13, x14
	umulh x19, x13, x14
	mul x20, x9, x11
	umulh x21, x9, x11
	mul x11, x13, x11
	mul x26, x13, x12
	umulh x18, x13, x12
	adds x13, x16, x3
	adcs x16, x17, x4
	umulh x7, x10, x15
	cinc x3, x5, hs
	adds x16, x16, x2
	adcs x2, x3, xzr
	mrs x3, NZCV
	mul x22, x8, x15
	umulh x23, x8, x15
	mul x17, x9, x15
	umulh x4, x9, x15
	adds x15, x16, x6
	adcs x16, x2, x7
	cset w5, hs
	msr NZCV, x3
	mul x24, x10, x14
	cinc x3, x5, hs
	adds x15, x15, x1
	adcs x16, x16, x19
	umulh x25, x10, x14
	cinc x3, x3, hs
	adds x16, x16, x20
	adcs x3, x3, x21
	mul x2, x8, x14
	mrs x6, NZCV
	adds x16, x16, x22
	umulh x1, x8, x14
	adcs x3, x3, x23
	cset w7, hs
	msr NZCV, x6
	mul x5, x10, x12
	cinc x6, x7, hs
	adds x16, x16, x24
	adcs x3, x3, x25
	umulh x10, x10, x12
	cset w7, hs
	adds x16, x16, x26
	adcs x18, x3, x18
	mul x19, x9, x14
	adc x6, x6, x7
	adds x17, x18, x17
	adcs x18, x6, x4
	umulh x14, x9, x14
	ldp x22, x21, [sp, #32]
	mul x3, x8, x12
	ldp x24, x23, [sp, #16]
	umulh x4, x8, x12
	cset w8, hs
	adds x17, x17, x2
	adcs x18, x18, x1
	umulh x1, x9, x12
	cset w2, hs
	adds x17, x17, x5
	adcs x10, x18, x10
	lsl x18, x11, #32
	mul x9, x9, x12
	adc x8, x2, x8
	adds x10, x10, x19
	adcs x12, x8, x14
	mov x8, #-4294967295
	ldp x20, x19, [sp, #48]
	cinc x14, x1, hs
	adds x10, x10, x3
	umulh x1, x11, x8
	adcs x12, x12, x4
	cset w2, hs
	adds x9, x12, x9
	lsr x12, x11, #32
	adc x14, x14, x2
	adds x13, x13, x18
	sub x11, x11, x18
	adcs x12, x15, x12
	lsr x18, x13, #32
	adcs x11, x16, x11
	lsl x16, x13, #32
	adcs x15, x17, x1
	umulh x17, x13, x8
	adcs x10, x10, xzr
	sub x13, x13, x16
	adcs x9, x9, xzr
	cset w1, hs
	adds x12, x12, x16
	adcs x11, x11, x18
	umulh x16, x12, x8
	adcs x13, x15, x13
	lsl x15, x12, #32
	adcs x10, x10, x17
	lsr x17, x12, #32
	adcs x9, x9, xzr
	sub x12, x12, x15
	adcs x14, x14, x1
	cset w18, hs
	adds x11, x11, x15
	adcs x13, x13, x17
	umulh x15, x11, x8
	adcs x10, x10, x12
	lsl x12, x11, #32
	adcs x9, x9, x16
	lsr x16, x11, #32
	adcs x14, x14, xzr
	sub x11, x11, x12
	cset w17, hs
	adds x13, x13, x12
	adcs x10, x10, x16
	eor w17, w18, w17
	adcs x9, x9, x11
	adcs x11, x14, x15
	mov w15, #-1
	cset w12, hs
	adds x14, x13, #1
	sbcs x15, x10, x15
	eor w12, w17, w12
	cset w16, lo
	eor w12, w12, #0x1
	subs x16, x9, x16
	sbcs x8, x11, x8
	cset w17, lo
	tst w17, w12
	csel x12, x13, x14, ne
	csel x10, x10, x15, ne
	csel x9, x9, x16, ne
	csel x8, x11, x8, ne
	stp x12, x10, [x0]
	stp x9, x8, [x0, #16]
	ldp x26, x25, [sp], #64
	ret
.Lfunc_end0:
fe_add:
	ldp x8, x11, [x2]
	mov x16, #-4294967295
	ldp x9, x10, [x1]
	ldp x12, x14, [x2, #16]
	adds x8, x9, x8
	ldp x9, x13, [x1, #16]
	adcs x10, x10, x11
	adcs x9, x9, x12
	mov w12, #-1
	adcs x11, x13, x14
	cset w13, hs
	adds x14, x8, #1
	sbcs x12, x10, x12
	eor w13, w13, #0x1
	cset w15, lo
	subs x15, x9, x15
	sbcs x16, x11, x16
	cset w17, lo
	tst w17, w13
	csel x8, x8, x14, ne
	csel x10, x10, x12, ne
	stp x8, x10, [x0]
	csel x8, x9, x15, ne
	csel x9, x11, x16, ne
	stp x8, x9, [x0, #16]
	ret
.Lfunc_end1:
fe_sub:
	ldp x8, x11, [x2]
	mov x14, #-4294967295
	ldp x9, x10, [x1]
	ldp x12, x13, [x2, #16]
	subs x8, x9, x8
	sbcs x10, x10, x11
	ldp x9, x11, [x1, #16]
	sbcs x9, x9, x12
	sbcs x11, x11, x13
	mov w13, #-1
	cset w12, lo
	csel x13, x13, xzr, lo
	csel x14, x14, xzr, lo
	sbfx x12, x12, #0, #1
	adds x8, x8, x12
	adcs x10, x10, x13
	stp x8, x10, [x0]
	adcs x8, x9, xzr
	adc x9, x14, x11
	stp x8, x9, [x0, #16]
	ret
.Lfunc_end2:
fe_inv:
	stp x30, x0, [sp, #-96]!
	ldp x10, x11, [x1, #16]
	mov w8, #254
	ldp x13, x14, [x1]
	mov x9, #-4294967295
	mov w12, #-1
	adrp x0, fe_inv.E
	add x0, x0, :lo12:fe_inv.E
	mov x16, x11
	mov x17, x10
	mov x18, x14
	mov x1, x13
	stp x28, x27, [sp, #16]
	stp x26, x25, [sp, #32]
	stp x24, x23, [sp, #48]
	stp x22, x21, [sp, #64]
	stp x20, x19, [sp, #80]
	b .LBB3_2
.LBB3_1:
	sub w8, w8, #1
	cmn w8, #1
	b.eq .LBB3_4
.LBB3_2:
	umulh x2, x1, x1
	mul x3, x18, x1
	umulh x4, x18, x1
	mul x5, x17, x1
	adds x2, x2, x3
	umulh x6, x17, x1
	cinc x21, x4, hs
	adds x2, x2, x3
	mul x7, x18, x18
	adcs x3, x21, x4
	umulh x19, x18, x18
	cinc x21, x6, hs
	adds x3, x3, x5
	adcs x21, x21, xzr
	mul x20, x16, x1
	mrs x22, NZCV
	adds x3, x3, x7
	umulh x4, x16, x1
	adcs x7, x21, x19
	cset w21, hs
	msr NZCV, x22
	mul x23, x17, x18
	cinc x21, x21, hs
	adds x3, x3, x5
	umulh x19, x17, x18
	adcs x5, x7, x6
	cinc x6, x21, hs
	adds x5, x5, x20
	adcs x6, x6, x4
	mul x7, x16, x18
	mrs x21, NZCV
	adds x5, x5, x23
	umulh x18, x16, x18
	adcs x6, x6, x19
	cset w22, hs
	msr NZCV, x21
	mul x24, x17, x17
	cinc x21, x22, hs
	adds x5, x5, x23
	adcs x6, x6, x19
	umulh x19, x17, x17
	cset w22, hs
	adds x5, x5, x20
	adcs x4, x6, x4
	mul x6, x16, x17
	adc x20, x21, x22
	adds x4, x4, x7
	adcs x20, x20, x18
	umulh x17, x16, x17
	cset w21, hs
	adds x4, x4, x24
	adcs x19, x20, x19
	umulh x20, x16, x16
	cset w22, hs
	adds x4, x4, x7
	mul x1, x1, x1
	adcs x18, x19, x18
	adc x7, x22, x21
	adds x18, x18, x6
	mul x16, x16, x16
	adcs x7, x7, x17
	cinc x19, x20, hs
	adds x18, x18, x6
	adcs x17, x7, x17
	lsl x6, x1, #32
	umulh x7, x1, x9
	cset w20, hs
	adds x16, x17, x16
	lsr x17, x1, #32
	sub x1, x1, x6
	adc x19, x19, x20
	adds x2, x2, x6
	adcs x17, x3, x17
	lsr x6, x2, #32
	adcs x1, x5, x1
	umulh x5, x2, x9
	adcs x3, x4, x7
	lsl x4, x2, #32
	adcs x18, x18, xzr
	adcs x16, x16, xzr
	sub x2, x2, x4
	cset w7, hs
	adds x17, x17, x4
	adcs x1, x1, x6
	umulh x4, x17, x9
	adcs x2, x3, x2
	lsl x3, x17, #32
	adcs x18, x18, x5
	lsr x5, x17, #32
	adcs x16, x16, xzr
	sub x17, x17, x3
	adcs x6, x19, x7
	cset w7, hs
	adds x1, x1, x3
	adcs x2, x2, x5
	umulh x3, x1, x9
	adcs x17, x18, x17
	lsl x18, x1, #32
	adcs x16, x16, x4
	lsr x4, x1, #32
	adcs x5, x6, xzr
	cset w6, hs
	adds x2, x2, x18
	sub x18, x1, x18
	adcs x17, x17, x4
	eor w6, w7, w6
	adcs x16, x16, x18
	adcs x3, x5, x3
	lsr w5, w8, #6
	cset w18, hs
	adds x1, x2, #1
	sbcs x4, x17, x12
	ldr x5, [x0, w5, uxtw #3]
	eor w18, w6, w18
	cset w19, lo
	eor w18, w18, #0x1
	subs x7, x16, x19
	lsr x5, x5, x8
	sbcs x6, x3, x9
	cset w19, lo
	tst w19, w18
	csel x1, x2, x1, ne
	csel x18, x17, x4, ne
	csel x17, x16, x7, ne
	csel x16, x3, x6, ne
	tbz w5, #0, .LBB3_1
	umulh x2, x1, x13
	mul x3, x1, x14
	umulh x4, x1, x14
	mul x5, x18, x13
	adds x2, x2, x3
	umulh x6, x18, x13
	cinc x4, x4, hs
	mul x7, x1, x10
	adds x2, x2, x5
	umulh x19, x1, x10
	adcs x4, x4, x6
	mul x20, x18, x14
	umulh x21, x18, x14
	cinc x6, x19, hs
	adds x4, x4, x7
	mul x22, x17, x13
	adcs x6, x6, xzr
	mrs x7, NZCV
	umulh x23, x17, x13
	mul x24, x1, x11
	mul x26, x18, x10
	umulh x27, x18, x10
	mul x5, x18, x11
	umulh x19, x18, x11
	adds x18, x4, x20
	adcs x4, x6, x21
	umulh x25, x1, x11
	cset w20, hs
	msr NZCV, x7
	cinc x7, x20, hs
	adds x18, x18, x22
	adcs x4, x4, x23
	mul x28, x17, x14
	cinc x7, x7, hs
	adds x4, x4, x24
	umulh x30, x17, x14
	adcs x7, x7, x25
	mrs x22, NZCV
	mul x15, x16, x13
	adds x4, x4, x26
	umulh x3, x16, x13
	adcs x7, x7, x27
	cset w24, hs
	msr NZCV, x22
	cinc x22, x24, hs
	adds x4, x4, x28
	mul x6, x17, x10
	adcs x7, x7, x30
	umulh x20, x17, x10
	cset w24, hs
	adds x15, x4, x15
	adcs x3, x7, x3
	mul x21, x16, x14
	adc x7, x22, x24
	adds x3, x3, x5
	adcs x5, x7, x19
	umulh x23, x16, x14
	cset w19, hs
	adds x3, x3, x6
	adcs x5, x5, x20
	mul x25, x17, x11
	cset w20, hs
	adds x3, x3, x21
	umulh x17, x17, x11
	adcs x5, x5, x23
	mul x4, x16, x10
	adc x19, x20, x19
	adds x5, x5, x25
	umulh x6, x16, x11
	adcs x17, x19, x17
	mul x1, x1, x13
	umulh x7, x16, x10
	cinc x6, x6, hs
	adds x4, x5, x4
	mul x16, x16, x11
	lsl x5, x1, #32
	adcs x17, x17, x7
	umulh x7, x1, x9
	cset w19, hs
	adds x16, x17, x16
	lsr x17, x1, #32
	sub x1, x1, x5
	adc x6, x6, x19
	adds x2, x2, x5
	adcs x17, x18, x17
	lsr x5, x2, #32
	adcs x15, x15, x1
	adcs x18, x3, x7
	lsl x3, x2, #32
	adcs x1, x4, xzr
	umulh x4, x2, x9
	adcs x16, x16, xzr
	sub x2, x2, x3
	cset w7, hs
	adds x17, x17, x3
	adcs x15, x15, x5
	umulh x3, x17, x9
	lsr x5, x17, #32
	adcs x18, x18, x2
	lsl x2, x17, #32
	adcs x1, x1, x4
	adcs x16, x16, xzr
	sub x17, x17, x2
	adcs x4, x6, x7
	cset w6, hs
	adds x15, x15, x2
	adcs x18, x18, x5
	umulh x2, x15, x9
	adcs x17, x1, x17
	lsl x1, x15, #32
	adcs x16, x16, x3
	adcs x3, x4, xzr
	lsr x4, x15, #32
	sub x15, x15, x1
	cset w5, hs
	adds x18, x18, x1
	adcs x17, x17, x4
	eor w5, w6, w5
	adcs x15, x16, x15
	adcs x16, x3, x2
	cset w1, hs
	adds x2, x18, #1
	sbcs x3, x17, x12
	eor w1, w5, w1
	cset w4, lo
	eor w1, w1, #0x1
	subs x4, x15, x4
	sbcs x6, x16, x9
	cset w5, lo
	tst w5, w1
	csel x1, x18, x2, ne
	csel x18, x17, x3, ne
	csel x17, x15, x4, ne
	csel x16, x16, x6, ne
	b .LBB3_1
.LBB3_4:
	ldr x8, [sp, #8]
	ldp x20, x19, [sp, #80]
	ldp x22, x21, [sp, #64]
	ldp x24, x23, [sp, #48]
	stp x1, x18, [x8]
	stp x17, x16, [x8, #16]
	ldp x26, x25, [sp, #32]
	ldp x28, x27, [sp, #16]
	ldr x30, [sp], #96
	ret
.Lfunc_end3:

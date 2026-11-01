; clang 18 output (flags in STORAGE.md), assembler directives and comments stripped
_mac_run:
	ldp x8, x12, [x0]
	ldp x11, x10, [x0, #16]
	ldr x9, [x0, #32]
	cmp x4, #1
	b.lt LBB0_4
; %bb.1:
LBB0_2:
	ldr x13, [x2], #8
	add x13, x1, x13, lsl #5
	ldrb w14, [x3], #1
	ldp x15, x16, [x13]
	mul x17, x15, x14
	umulh x15, x15, x14
	umulh x5, x16, x14
	mul x16, x16, x14
	ldp x6, x13, [x13, #16]
	umulh x7, x6, x14
	mul x6, x6, x14
	umulh x19, x13, x14
	mul x13, x13, x14
	adds x12, x15, x12
	cset w14, hs
	adds x12, x12, x16
	cinc x14, x14, hs
	adds x11, x5, x11
	cset w15, hs
	adds x11, x11, x6
	cinc x15, x15, hs
	adds x10, x7, x10
	cset w16, hs
	adds x10, x10, x13
	cinc x13, x16, hs
	adds x8, x17, x8
	adcs x12, x12, xzr
	adcs x11, x11, x14
	adcs x10, x10, x15
	add x9, x9, x19
	adc x9, x9, x13
	subs x4, x4, #1
	b.ne LBB0_2
; %bb.3:
LBB0_4:
	stp x8, x12, [x0]
	stp x11, x10, [x0, #16]
	str x9, [x0, #32]
	ret
_group_run:
; %bb.0:
	stp x20, x19, [sp, #-16]!
	ldp x8, x12, [x0]
	ldp x11, x10, [x0, #16]
	ldr x9, [x0, #32]
	cmp x3, #1
	b.lt LBB9_2
LBB9_1:
	ldr x13, [x2], #8
	add x13, x1, x13, lsl #5
	ldp x14, x15, [x13]
	adds x8, x14, x8
	adcs x12, x15, x12
	ldp x14, x13, [x13, #16]
	adcs x11, x14, x11
	adcs x10, x13, x10
	cinc x9, x9, hs
	subs x3, x3, #1
	b.ne LBB9_1
LBB9_2:
	cmp x6, #1
	b.lt LBB9_4
LBB9_3:
	ldr x13, [x4], #8
	add x13, x1, x13, lsl #5
	ldrb w14, [x5], #1
	ldp x15, x16, [x13]
	mul x17, x15, x14
	umulh x15, x15, x14
	umulh x2, x16, x14
	mul x16, x16, x14
	ldp x3, x13, [x13, #16]
	umulh x7, x3, x14
	mul x3, x3, x14
	umulh x19, x13, x14
	mul x13, x13, x14
	adds x12, x15, x12
	cset w14, hs
	adds x12, x12, x16
	cinc x14, x14, hs
	adds x11, x2, x11
	cset w15, hs
	adds x11, x11, x3
	cinc x15, x15, hs
	adds x10, x7, x10
	cset w16, hs
	adds x10, x10, x13
	cinc x13, x16, hs
	adds x8, x17, x8
	adcs x12, x12, xzr
	adcs x11, x11, x14
	adcs x10, x10, x15
	add x9, x9, x19
	adc x9, x9, x13
	subs x6, x6, #1
	b.ne LBB9_3
LBB9_4:
	stp x8, x12, [x0]
	stp x11, x10, [x0, #16]
	str x9, [x0, #32]
	ldp x20, x19, [sp], #16
	ret
; per-row tail (row_tail.c, same clang flags): send a (pos, neg) pair, combine an incoming pair, finish (pos - neg) mod r
_send_pair:
	mov x8, x0
	ldp x9, x10, [x1]
	stp x9, x10, [x0]
	ldp x9, x10, [x1, #16]
	stp x9, x10, [x0, #16]
	ldp x9, x10, [x1, #32]
	stp x9, x10, [x0, #32]
	ldp x9, x10, [x1, #48]
	stp x9, x10, [x0, #48]
	ldp x9, x10, [x1, #64]
	add x0, x0, #80
	stp x9, x10, [x8, #64]
	ret
_combine:
	ldr x8, [x1]
	ldp x9, x10, [x0]
	adds x8, x9, x8
	str x8, [x0]
	ldr x8, [x1, #8]
	adcs x9, x10, xzr
	cset w10, hs
	adds x8, x9, x8
	str x8, [x0, #8]
	ldr x8, [x1, #16]
	ldp x9, x11, [x0, #16]
	adcs x9, x9, x10
	cset w10, hs
	adds x8, x9, x8
	str x8, [x0, #16]
	ldr x8, [x1, #24]
	adcs x9, x11, x10
	cset w10, hs
	adds x8, x9, x8
	str x8, [x0, #24]
	ldr x8, [x1, #32]
	ldp x9, x11, [x0, #32]
	adc x9, x9, x10
	add x8, x9, x8
	str x8, [x0, #32]
	ldr x8, [x1, #40]
	adds x8, x11, x8
	str x8, [x0, #40]
	ldr x8, [x1, #48]
	ldp x9, x10, [x0, #48]
	adcs x9, x9, xzr
	cset w11, hs
	adds x8, x9, x8
	str x8, [x0, #48]
	ldr x8, [x1, #56]
	adcs x9, x10, x11
	cset w10, hs
	adds x8, x9, x8
	str x8, [x0, #56]
	ldr x8, [x1, #64]
	ldp x9, x11, [x0, #64]
	adcs x9, x9, x10
	cset w10, hs
	adds x8, x9, x8
	str x8, [x0, #64]
	ldr x8, [x1, #72]
	adc x9, x11, x10
	add x8, x9, x8
	str x8, [x0, #72]
	add x0, x1, #80
	ret
_finish:
	ldp x8, x9, [x2]
	ldp x11, x10, [x1, #32]
	subs x8, x8, x10
	ngc x10, xzr
	and x10, x10, #0x1
	ldp x12, x13, [x1, #48]
	adds x10, x10, x12
	cset w12, hs
	subs x9, x9, x10
	ngc x10, x12
	ldp x12, x14, [x1]
	adds x15, x8, x12
	adcs x8, x14, xzr
	cset w12, hs
	and x10, x10, #0x1
	adds x10, x10, x13
	cset w13, hs
	ldp x14, x16, [x2, #16]
	subs x10, x14, x10
	ngc x13, x13
	adds x14, x8, x9
	ldp x8, x9, [x1, #16]
	adcs x8, x8, x12
	cset w17, hs
	and x12, x13, #0x1
	ldp x13, x1, [x1, #64]
	adds x12, x12, x13
	cset w13, hs
	subs x16, x16, x12
	ngc x13, x13
	adds x12, x8, x10
	adcs x8, x9, x17
	cset w9, hs
	ldr x10, [x2, #32]
	and x13, x13, #0x1
	sub x10, x10, x1
	sub x13, x10, x13
	adds x10, x8, x16
	adc x8, x11, x9
	add x8, x8, x13
	ldp x9, x11, [x3]
	umulh x13, x10, x9
	umulh x16, x10, x11
	mul x17, x10, x11
	umulh x1, x8, x9
	mul x9, x8, x9
	adds x13, x17, x13
	cset w17, hs
	cmn x13, x9
	madd x9, x8, x11, x16
	add x9, x9, x1
	adc x13, x9, x17
	ldp x16, x17, [x4]
	umulh x11, x16, x13
	mul x9, x16, x13
	subs x9, x15, x9
	ngc x15, xzr
	umulh x1, x17, x13
	mul x2, x17, x13
	adds x11, x2, x11
	cinc x1, x1, hs
	and x15, x15, #0x1
	adds x11, x11, x15
	cset w15, hs
	subs x11, x14, x11
	ngc x14, x15
	ldp x15, x2, [x4, #16]
	umulh x3, x15, x13
	mul x4, x15, x13
	adds x1, x4, x1
	cinc x3, x3, hs
	and x14, x14, #0x1
	adds x14, x1, x14
	cset w1, hs
	subs x12, x12, x14
	ngc x14, x1
	mul x1, x2, x13
	adds x4, x1, x3
	and x14, x14, #0x1
	adds x14, x4, x14
	cset w4, hs
	subs x10, x10, x14
	ngc x14, x4
	cmn x1, x3
	umulh x13, x2, x13
	and x14, x14, #0x1
	adc x1, x14, x13
	subs x13, x9, x16
	ngc x14, xzr
	and x14, x14, #0x1
	adds x14, x14, x17
	cset w16, hs
	subs x14, x11, x14
	ngc x16, x16
	and x16, x16, #0x1
	adds x15, x16, x15
	cset w17, hs
	subs x16, x12, x15
	ngc x15, x17
	and x15, x15, #0x1
	adds x15, x15, x2
	cset w17, hs
	subs x15, x10, x15
	ngc x17, x17
	and x17, x17, #0x1
	cmp x8, x1
	csel x8, x17, xzr, eq
	cbz x8, LBB12_2
	neg x8, x8
	and x9, x9, x8
	and x11, x11, x8
	stp x9, x11, [x0]
	and x9, x12, x8
	and x15, x10, x8
	stp x9, x15, [x0, #16]
	ret
LBB12_2:
	stp x13, x14, [x0]
	stp x16, x15, [x0, #16]
	ret

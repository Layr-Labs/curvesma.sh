	.macro FM a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p
	mul \a,\b,\c
	umulh \d,\b,\c
	mul \e,\b,\f
	adds \d,\d,\e
	umulh \e,\b,\f
	mul \g,\b,\h
	adcs \e,\e,\g
	umulh \g,\b,\h
	mul \i,\b,\j
	adcs \g,\g,\i
	umulh \i,\b,\j
	mul \b,\k,\j
	adcs \i,\i,\b
	umulh \b,\k,\j
	adc \b,\b,xzr
	mul \l,\k,\c
	adds \d,\d,\l
	umulh \l,\k,\c
	adcs \e,\e,\l
	umulh \l,\m,\c
	adcs \g,\g,\l
	umulh \l,\n,\c
	adcs \i,\i,\l
	umulh \l,\m,\h
	adc \l,\l,xzr
	mul \o,\m,\c
	adds \e,\e,\o
	mul \o,\n,\c
	adcs \g,\g,\o
	mul \o,\m,\h
	adcs \i,\i,\o
	adcs \b,\b,\l
	umulh \l,\m,\j
	adc \l,\l,xzr
	mul \o,\k,\f
	adds \e,\e,\o
	umulh \o,\k,\f
	adcs \g,\g,\o
	umulh \o,\k,\h
	adcs \i,\i,\o
	umulh \o,\n,\f
	adcs \b,\b,\o
	umulh \o,\n,\h
	adc \o,\o,xzr
	mul \c,\k,\h
	adds \g,\g,\c
	umulh \c,\m,\f
	adcs \i,\i,\c
	mul \c,\m,\j
	adcs \b,\b,\c
	adcs \l,\l,\o
	umulh \o,\n,\j
	adc \o,\o,xzr
	mul \c,\m,\f
	adds \g,\g,\c
	mul \c,\n,\f
	adcs \i,\i,\c
	mul \c,\n,\h
	adcs \b,\b,\c
	mul \c,\n,\j
	adcs \l,\l,\c
	adc \o,\o,xzr
	mul \c,\i,\p
	umulh \j,\i,\p
	mul \i,\b,\p
	adds \i,\i,\j
	umulh \j,\b,\p
	mul \b,\l,\p
	adcs \b,\b,\j
	umulh \j,\l,\p
	mul \l,\o,\p
	adcs \l,\l,\j
	umulh \j,\o,\p
	adc \j,\j,xzr
	adds \a,\a,\c
	adcs \d,\d,\i
	adcs \e,\e,\b
	adcs \g,\g,\l
	adc \j,\j,xzr
	.endm
	.macro PR
	stp x19,x20,[sp,#-96]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	mov x4,#977
	movk x4,#1,lsl #32
	.endm
	.macro EP
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	ret
	.endm
	.macro SQB
	ldp x5,x6,[x1,#0]
	ldp x7,x8,[x1,#16]
	mul x9,x5,x8
	mul x10,x6,x7
	adds x9,x9,x10
	umulh x10,x5,x8
	umulh x11,x6,x7
	adcs x10,x10,x11
	umulh x11,x6,x8
	adc x11,x11,xzr
	umulh x12,x5,x6
	mul x13,x5,x7
	adds x12,x12,x13
	umulh x13,x5,x7
	adcs x9,x9,x13
	mul x13,x6,x8
	adcs x10,x10,x13
	mul x13,x7,x8
	adcs x11,x11,x13
	umulh x13,x7,x8
	adc x13,x13,xzr
	mul x14,x5,x6
	lsr x15,x13,#63
	extr x13,x13,x11,#63
	extr x11,x11,x10,#63
	extr x10,x10,x9,#63
	extr x9,x9,x12,#63
	extr x12,x12,x14,#63
	mul x16,x5,x5
	umulh x17,x5,x5
	adds x14,x17,x14,lsl #1
	mul x17,x6,x6
	adcs x12,x12,x17
	umulh x17,x6,x6
	adcs x9,x9,x17
	mul x17,x7,x7
	adcs x10,x10,x17
	umulh x17,x7,x7
	adcs x11,x11,x17
	mul x17,x8,x8
	adcs x13,x13,x17
	umulh x17,x8,x8
	adc x15,x15,x17
	mul x17,x10,x4
	umulh x8,x10,x4
	mul x10,x11,x4
	adds x10,x10,x8
	umulh x8,x11,x4
	mul x11,x13,x4
	adcs x11,x11,x8
	umulh x8,x13,x4
	mul x13,x15,x4
	adcs x13,x13,x8
	umulh x8,x15,x4
	adc x8,x8,xzr
	adds x16,x16,x17
	adcs x14,x14,x10
	adcs x12,x12,x11
	adcs x9,x9,x13
	adc x8,x8,xzr
	.endm
	.text
	.p2align 2
	.globl t_mul, t_sqr, t_sub, t_reduce, t_sqr_sub2, t_mul_sub
t_mul:
	PR
	ldp x5,x6,[x1,#0]
	ldp x7,x8,[x1,#16]
	ldp x9,x10,[x2,#0]
	ldp x11,x12,[x2,#16]
	FM x13,x5,x9,x14,x15,x10,x16,x11,x17,x12,x6,x19,x7,x8,x20,x4
	mul x19,x12,x4
	umulh x5,x12,x4
	adds x13,x13,x19
	adcs x14,x14,x5
	adcs x15,x15,xzr
	adc x16,x16,xzr
	stp x13,x14,[x0,#0]
	stp x15,x16,[x0,#16]
	EP
t_sqr:
	PR
	SQB
	mul x13,x8,x4
	umulh x11,x8,x4
	adds x16,x16,x13
	adcs x14,x14,x11
	adcs x12,x12,xzr
	adc x9,x9,xzr
	stp x16,x14,[x0,#0]
	stp x12,x9,[x0,#16]
	EP
t_sub:
	PR
	ldp x5,x6,[x1,#0]
	ldp x7,x8,[x1,#16]
	ldp x9,x10,[x2,#0]
	ldp x11,x12,[x2,#16]
	subs x5,x5,x9
	sbcs x6,x6,x10
	sbcs x7,x7,x11
	sbcs x8,x8,x12
	csel x12,x4,xzr,lo
	subs x5,x5,x12
	sbcs x6,x6,xzr
	sbcs x7,x7,xzr
	sbc x8,x8,xzr
	stp x5,x6,[x0,#0]
	stp x7,x8,[x0,#16]
	EP
t_reduce:
	PR
	ldp x5,x6,[x1,#0]
	ldp x7,x8,[x1,#16]
	ldp x9,x10,[x1,#32]
	ldp x11,x12,[x1,#48]
	mul x13,x9,x4
	umulh x14,x9,x4
	mul x9,x10,x4
	adds x9,x9,x14
	umulh x14,x10,x4
	mul x10,x11,x4
	adcs x10,x10,x14
	umulh x14,x11,x4
	mul x11,x12,x4
	adcs x11,x11,x14
	umulh x14,x12,x4
	adc x14,x14,xzr
	adds x5,x5,x13
	adcs x6,x6,x9
	adcs x7,x7,x10
	adcs x8,x8,x11
	adc x14,x14,xzr
	mul x11,x14,x4
	umulh x10,x14,x4
	adds x5,x5,x11
	adcs x6,x6,x10
	adcs x7,x7,xzr
	adc x8,x8,xzr
	stp x5,x6,[x0,#0]
	stp x7,x8,[x0,#16]
	EP
t_sqr_sub2:
	PR
	SQB
	ldp x13,x11,[x2,#0]
	subs x16,x16,x13
	sbcs x14,x14,x11
	ldp x13,x11,[x2,#16]
	sbcs x12,x12,x13
	sbcs x9,x9,x11
	sbc x8,x8,xzr
	ldp x11,x13,[x3,#0]
	subs x16,x16,x11
	sbcs x14,x14,x13
	ldp x11,x13,[x3,#16]
	sbcs x12,x12,x11
	sbcs x9,x9,x13
	sbc x8,x8,xzr
	asr x13,x8,#63
	mul x11,x8,x4
	smulh x10,x8,x4
	adds x16,x16,x11
	adcs x14,x14,x10
	adcs x12,x12,x13
	adc x9,x9,x13
	stp x16,x14,[x0,#0]
	stp x12,x9,[x0,#16]
	EP
t_mul_sub:
	PR
	ldp x5,x6,[x1,#0]
	ldp x7,x8,[x1,#16]
	ldp x9,x10,[x2,#0]
	ldp x11,x12,[x2,#16]
	FM x13,x5,x9,x14,x15,x10,x16,x11,x17,x12,x6,x19,x7,x8,x20,x4
	ldp x19,x5,[x3,#0]
	subs x13,x13,x19
	sbcs x14,x14,x5
	ldp x19,x5,[x3,#16]
	sbcs x15,x15,x19
	sbcs x16,x16,x5
	sbc x12,x12,xzr
	asr x5,x12,#63
	mul x19,x12,x4
	smulh x17,x12,x4
	adds x13,x13,x19
	adcs x14,x14,x17
	adcs x15,x15,x5
	adc x16,x16,x5
	stp x13,x14,[x0,#0]
	stp x15,x16,[x0,#16]
	EP

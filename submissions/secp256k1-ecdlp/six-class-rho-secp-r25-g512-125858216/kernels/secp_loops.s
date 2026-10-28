	.macro FS d,a,b,e,c,f,g,h,i,j,k,l,t,s
	subs \d,\a,\b
	sbcs \e,\c,\f
	sbcs \g,\h,\i
	sbcs \j,\k,\l
	csel \t, x4, xzr, lo
	subs \d,\d,\t
	b.lo .Lr\s
.Lb\s:
	.endm
	.macro FO l,s,h,d,z,e,y,c
	mul \l,\s,x4
	umulh \h,\s,x4
	adds \d,\z,\l
	adcs \e,\y,\h
	b.cs .Lr\c
.Lb\c:
	.endm
	.macro OD a,b,c,d,e,f,g,h,i,j,k
	mul \a,\b,\c
	mul \d,\e,\f
	adds \a,\a,\d
	umulh \d,\b,\c
	umulh \g,\e,\f
	adcs \d,\d,\g
	umulh \g,\e,\c
	adc \g,\g,xzr
	umulh \h,\b,\e
	mul \i,\b,\f
	adds \h,\h,\i
	umulh \i,\b,\f
	adcs \a,\a,\i
	mul \i,\e,\c
	adcs \d,\d,\i
	mul \i,\f,\c
	adcs \g,\g,\i
	umulh \i,\f,\c
	adc \i,\i,xzr
	mul \j,\b,\e
	lsr \k,\i,#63
	extr \i,\i,\g,#63
	extr \g,\g,\d,#63
	extr \d,\d,\a,#63
	extr \a,\a,\h,#63
	extr \h,\h,\j,#63
	.endm
	.macro SQ a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r
	mul \a,\b,\b
	umulh \c,\b,\b
	adds \d,\c,\d,lsl #1
	mul \c,\e,\e
	adcs \f,\f,\c
	umulh \c,\e,\e
	adcs \g,\g,\c
	mul \c,\h,\h
	adcs \i,\i,\c
	umulh \c,\h,\h
	adcs \j,\j,\c
	mul \c,\k,\k
	adcs \l,\l,\c
	umulh \c,\k,\k
	adc \m,\m,\c
	mul \c,\i,\n
	umulh \o,\i,\n
	mul \i,\j,\n
	adds \i,\i,\o
	umulh \o,\j,\n
	mul \j,\l,\n
	adcs \j,\j,\o
	umulh \o,\l,\n
	mul \l,\m,\n
	adcs \l,\l,\o
	umulh \o,\m,\n
	adc \o,\o,xzr
	adds \a,\a,\c
	adcs \d,\d,\i
	adcs \f,\f,\j
	adcs \g,\g,\l
	adc \o,\o,xzr
	ldp \l,\j,[\p,#-32]
	ldp \i,\c,[\p,#-16]
	subs \a,\a,\l
	sbcs \d,\d,\j
	sbcs \f,\f,\i
	sbcs \g,\g,\c
	adc \o,\o,xzr
	ldp \m,\q,[\r,#192]
	subs \a,\a,\m
	sbcs \d,\d,\q
	ldp \m,\q,[\r,#208]
	sbcs \f,\f,\m
	sbcs \g,\g,\q
	adc \o,\o,xzr
	mul \q,\o,\n
	umulh \m,\o,\n
	adds \a,\a,\q
	adcs \d,\d,\m
	.endm
	.macro FM a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u
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
	mul \k,\l,\j
	adcs \i,\i,\k
	umulh \k,\l,\j
	adc \k,\k,xzr
	mul \m,\l,\c
	adds \d,\d,\m
	umulh \m,\l,\c
	adcs \e,\e,\m
	umulh \m,\n,\c
	adcs \g,\g,\m
	umulh \m,\o,\c
	adcs \i,\i,\m
	umulh \m,\n,\h
	adc \m,\m,xzr
	mul \p,\n,\c
	adds \e,\e,\p
	mul \p,\o,\c
	adcs \g,\g,\p
	mul \p,\n,\h
	adcs \i,\i,\p
	adcs \k,\k,\m
	umulh \m,\n,\j
	adc \m,\m,xzr
	mul \p,\l,\f
	adds \e,\e,\p
	umulh \p,\l,\f
	adcs \g,\g,\p
	umulh \p,\l,\h
	adcs \i,\i,\p
	umulh \p,\o,\f
	adcs \k,\k,\p
	umulh \p,\o,\h
	adc \p,\p,xzr
	mul \q,\l,\h
	adds \g,\g,\q
	umulh \q,\n,\f
	adcs \i,\i,\q
	mul \q,\n,\j
	adcs \k,\k,\q
	adcs \m,\m,\p
	umulh \p,\o,\j
	adc \p,\p,xzr
	mul \q,\n,\f
	adds \g,\g,\q
	mul \q,\o,\f
	adcs \i,\i,\q
	mul \q,\o,\h
	adcs \k,\k,\q
	mul \q,\o,\j
	adcs \m,\m,\q
	adc \p,\p,xzr
	mul \q,\i,\r
	umulh \s,\i,\r
	mul \i,\k,\r
	adds \i,\i,\s
	umulh \s,\k,\r
	mul \k,\m,\r
	adcs \k,\k,\s
	umulh \s,\m,\r
	mul \m,\p,\r
	adcs \m,\m,\s
	umulh \s,\p,\r
	adc \s,\s,xzr
	adds \a,\a,\q
	adcs \d,\d,\i
	adcs \t,\e,\k
	adcs \u,\g,\m
	adc \s,\s,xzr
	.endm
	.macro RCHG
	str x0,[sp,#-16]!
	ldr x0,[sp,#16]
	sub x0,x0,#20
	str x0,[sp,#16]
	ldr x0,[sp],#16
	.endm
	.macro SB a,b,c,s
.Lr\s:
	RCHG
	sbcs \a,\a,xzr
	sbcs \b,\b,xzr
	sbc \c,\c,xzr
	b .Lb\s
	.endm
	.macro SC a,b,s
.Lr\s:
	RCHG
	adcs \a,\a,xzr
	adc \b,\b,xzr
	b .Lb\s
	.endm
	.macro PR
	stp x19,x20,[sp,#-96]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	.endm
	.macro EP
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	.endm
	.macro W1 k,p
	ldp x10,x11,[\p,#0]
	ldp x12,x13,[\p,#16]
	ldp x16,x17,[x1,#16]
	ldp x14,x15,[x1],#64
	FS x10,x10,x14,x11,x11,x15,x12,x12,x16,x13,x13,x17,x17,1_0_\k
	stp x5,x6,[x3,#32]
	stp x7,x8,[x3,#48]
	stp x12,x13,[x3,#16]
	stp x10,x11,[x3],#64
	FM x9,x10,x5,x17,x16,x6,x15,x7,x14,x8,x10,x11,x19,x12,x13,x20,x21,x4,x13,x7,x8
	FO x19,x13,x10,x5,x9,x6,x17,1_1_\k
	.endm
	.macro S1 k
	SB x11,x12,x13,1_0_\k
	SC x7,x8,1_1_\k
	.endm
	.macro P1 i,j
	ldp x9,x26,[x2],#16
	W1 \i,x9
	W1 \j,x26
	.endm
	.macro Q1 v
	P1 0\v,1\v
	P1 2\v,3\v
	P1 4\v,5\v
	P1 6\v,7\v
	.endm
	.macro H1 v
	Q1 0\v
	Q1 1\v
	Q1 2\v
	Q1 3\v
	.endm
	.text
	.p2align 2
	.globl pass1
	.type pass1,%function
pass1:
	PR
	mov x4,#977
	movk x4,#1,lsl #32
	ldr x9,[x2],#8
	ldp x5,x6,[x9,#0]
	ldp x7,x8,[x9,#16]
	ldp x12,x13,[x1,#16]
	ldp x10,x11,[x1],#64
	FS x5,x5,x10,x6,x6,x11,x7,x7,x12,x8,x8,x13,x13,1_0_w0
	add x3,x3,#64
	sub x0,x0,#1
	lsr x29,x0,#9
	and x0,x0,#511
	cbz x0, 3f
2:
	ldr x9,[x2],#8
	W1 256,x9
	subs x0,x0,#1
	b.ne 2b
3:
	cbz x29, 4f
1:
	H1 0
	H1 1
	H1 2
	H1 3
	H1 4
	H1 5
	H1 6
	H1 7
	H1 8
	H1 9
	H1 10
	H1 11
	H1 12
	H1 13
	H1 14
	H1 15
	subs x29,x29,#1
	b.ne 1b
4:
	stp x5,x6,[x3,#32]
	stp x7,x8,[x3,#48]
	EP
	ret
	.macro Q3 v
	S1 0\v
	S1 1\v
	S1 2\v
	S1 3\v
	S1 4\v
	S1 5\v
	S1 6\v
	S1 7\v
	.endm
	.macro H3 v
	Q3 0\v
	Q3 1\v
	Q3 2\v
	Q3 3\v
	.endm
pass1_rare:
	H3 0
	H3 1
	H3 2
	H3 3
	H3 4
	H3 5
	H3 6
	H3 7
	H3 8
	H3 9
	H3 10
	H3 11
	H3 12
	H3 13
	H3 14
	H3 15
	S1 256
	SB x6,x7,x8,1_0_w0
	.macro RF o,a,b,c,d,e,f,g,h,i,t,u
	mul \e,\o,\a
	umulh \t,\o,\a
	mul \u,\o,\b
	adds \f,\t,\u
	umulh \t,\o,\b
	mul \u,\o,\c
	adcs \g,\t,\u
	umulh \t,\o,\c
	mul \u,\o,\d
	adcs \h,\t,\u
	umulh \t,\o,\d
	adc \i,\t,xzr
	.endm
	.macro RA o,a,b,c,d,r,s,t,v,u,e,f,g,h,i,j,y,m,n
	mul \r,\o,\a
	umulh \s,\o,\a
	mul \u,\o,\b
	adds \s,\s,\u
	umulh \t,\o,\b
	mul \u,\o,\c
	adcs \t,\t,\u
	umulh \v,\o,\c
	mul \u,\o,\d
	adcs \v,\v,\u
	umulh \u,\o,\d
	adc \u,\u,xzr
	adds \e,\e,\r
	adcs \f,\f,\s
	adcs \m,\g,\t
	adcs \n,\h,\v
	adcs \i,\i,\u
	adc \j,\y,xzr
	.endm
	.macro SH l,a,d,s
	mul \l,\d,x4
	umulh \d,\d,x4
	adds \a,\a,\d
	b.cs .Lr2_\s
.Lb2_\s:
	.endm
	.macro FD e,f,i,j,o,p,l,h,s
	mul \l,\i,x4
	umulh \h,\i,x4
	madd \h,\j,x4,\h
	adds \o,\e,\l
	adcs \p,\f,\h
	b.cs .Lr2_\s
.Lb2_\s:
	.endm
	.macro ST a,b,c,d,s
.Lr2_\s:
	RCHG
	adcs \c,\c,xzr
	adcs \d,\d,xzr
	b.cc .Lb2_\s
	adds \a,\a,x4
	adcs \b,\b,xzr
	adcs \c,\c,xzr
	adc \d,\d,xzr
	b .Lb2_\s
	.endm
	.macro W2P k
	ldp x5,x6,[x3,#32]
	RF x5,x25,x26,x27,x28,x10,x11,x12,x13,x7,x8,x9
	ldp x5,x8,[x3],#-64
	RF x5,x25,x26,x27,x28,x9,x14,x15,x16,x17,x19,x20
	SH x5,x25,x28,7_\k
	RA x6,x5,x25,x26,x27,x19,x20,x21,x22,x23,x10,x11,x12,x13,x7,x24,xzr,x12,x13
	RA x8,x5,x25,x26,x27,x6,x19,x20,x21,x22,x9,x14,x15,x16,x17,x23,xzr,x15,x16
	SH x28,x5,x27,8_\k
	ldp x6,x8,[x3,#112]
	RA x6,x28,x5,x25,x26,x19,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	ldp x6,x19,[x3,#80]
	RA x6,x28,x5,x25,x26,x20,x21,x22,x29,x30,x9,x14,x15,x16,x17,x23,x23,x15,x16
	SH x27,x28,x26,9_\k
	RA x8,x27,x28,x5,x25,x6,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	RA x19,x27,x28,x5,x25,x6,x8,x20,x21,x22,x9,x14,x15,x16,x17,x23,x23,x27,x28
	FD x10,x11,x7,x24,x10,x11,x5,x6,0_\k
	FD x9,x14,x17,x23,x25,x26,x5,x6,1_\k
	.endm
	.macro W2 k,p
	ldp x7,x8,[\p,#32]
	ldp x17,x6,[\p,#48]
	ldp x15,x20,[x1,#0]
	ldp x14,x19,[x1,#16]
	FS x7,x7,x15,x8,x8,x20,x17,x17,x14,x6,x6,x19,x22,2_2_\k
	FM x22,x7,x10,x21,x9,x11,x16,x12,x23,x13,x7,x8,x24,x17,x6,x29,x10,x4,x13,x9,x16
	FO x24,x13,x7,x22,x22,x21,x21,2_3_\k
	OD x7,x22,x16,x24,x21,x9,x13,x23,x10,x29,x6
	SQ x12,x22,x11,x29,x21,x23,x7,x9,x24,x13,x16,x10,x6,x4,x17,x1,x8,\p
	b.cs .Lr2_4_\k
.Lb2_4_\k:
	FS x6,x12,x10,x8,x29,x13,x17,x23,x24,x30,x7,x11,x11,2_5_\k
	stp x12,x29,[x1,#-32]
	stp x23,x7,[x1,#-16]
	FM x7,x22,x6,x23,x29,x8,x12,x17,x11,x30,x22,x21,x24,x9,x16,x13,x6,x4,x30,x29,x12
	adds x7,x7,x15
	adcs x23,x23,x20
	adcs x29,x29,x14
	adcs x12,x12,x19
	adc x30,x30,xzr
	FO x19,x30,x14,x7,x7,x23,x23,2_6_\k
	stp x29,x12,[x1,#16]
	stp x7,x23,[x1],#-64
	.endm
	.macro S2 k
	ST x5,x25,x26,x27,7_\k
	ST x28,x5,x25,x26,8_\k
	ST x27,x28,x5,x25,9_\k
	ST x10,x11,x12,x13,0_\k
	ST x25,x26,x27,x28,1_\k
	SB x8,x17,x6,2_2_\k
	SC x9,x16,2_3_\k
	SC x23,x7,2_4_\k
	SB x8,x17,x30,2_5_\k
	SC x29,x12,2_6_\k
	.endm
	.macro P2 i,j
	W2P \i
	ldp x0,x5,[x2,#-16]!
	W2 \i,x5
	W2P \j
	W2 \j,x0
	.endm
	.macro Q2 v
	P2 0\v,1\v
	P2 2\v,3\v
	P2 4\v,5\v
	P2 6\v,7\v
	.endm
	.macro H2 v
	Q2 0\v
	Q2 1\v
	Q2 2\v
	Q2 3\v
	.endm
	.macro Q4 v
	S2 0\v
	S2 1\v
	S2 2\v
	S2 3\v
	S2 4\v
	S2 5\v
	S2 6\v
	S2 7\v
	.endm
	.macro H4 v
	Q4 0\v
	Q4 1\v
	Q4 2\v
	Q4 3\v
	.endm
	.text
	.p2align 2
	.globl pass2
	.type pass2,%function
pass2:
	PR
	ldp x25,x26,[x4]
	ldp x27,x28,[x4,#16]
	mov x4,#977
	movk x4,#1,lsl #32
	sub x5,x0,#1
	add x1,x1,x5,lsl #6
	add x1,x1,#32
	add x2,x2,x0,lsl #3
	add x3,x3,x5,lsl #6
	and x0,x5,#511
	lsr x18,x5,#9
	cbz x0, 3f
2:
	W2P 256
	ldr x5,[x2,#-8]!
	W2 256,x5
	subs x0,x0,#1
	b.ne 2b
3:
	cbz x18, 4f
1:
	H2 0
	H2 1
	H2 2
	H2 3
	H2 4
	H2 5
	H2 6
	H2 7
	H2 8
	H2 9
	H2 10
	H2 11
	H2 12
	H2 13
	H2 14
	H2 15
	subs x18,x18,#1
	b.ne 1b
4:
	ldr x5,[x2,#-8]
	ldp x6,x7,[x5,#32]
	ldp x8,x9,[x5,#48]
	ldp x10,x11,[x1,#0]
	ldp x12,x13,[x1,#16]
	FS x6,x6,x10,x7,x7,x11,x8,x8,x12,x9,x9,x13,x14,2_0_w0
	FM x14,x6,x25,x15,x16,x26,x17,x27,x19,x28,x6,x7,x20,x8,x9,x21,x25,x4,x28,x16,x17
	FO x20,x28,x6,x14,x14,x15,x15,2_1_w0
	OD x6,x14,x17,x20,x15,x16,x28,x19,x25,x21,x9
	SQ x27,x14,x26,x21,x15,x19,x6,x16,x20,x28,x17,x25,x9,x4,x8,x1,x7,x5
	b.cs .Lr2_2_w0
.Lb2_2_w0:
	FS x9,x27,x25,x7,x21,x28,x8,x19,x20,x22,x6,x26,x26,2_3_w0
	stp x27,x21,[x1,#-32]
	stp x19,x6,[x1,#-16]
	FM x6,x14,x9,x19,x21,x7,x27,x8,x26,x22,x14,x15,x20,x16,x17,x28,x9,x4,x22,x21,x27
	adds x6,x6,x10
	adcs x19,x19,x11
	adcs x21,x21,x12
	adcs x27,x27,x13
	adc x22,x22,xzr
	FO x13,x22,x12,x6,x6,x19,x19,2_4_w0
	stp x6,x19,[x1,#0]
	stp x21,x27,[x1,#16]
	EP
	ret
pass2_rare:
	S2 256
	H4 0
	H4 1
	H4 2
	H4 3
	H4 4
	H4 5
	H4 6
	H4 7
	H4 8
	H4 9
	H4 10
	H4 11
	H4 12
	H4 13
	H4 14
	H4 15
	SB x7,x8,x9,2_0_w0
	SC x16,x17,2_1_w0
	SC x19,x6,2_2_w0
	SB x7,x8,x22,2_3_w0
	SC x21,x27,2_4_w0

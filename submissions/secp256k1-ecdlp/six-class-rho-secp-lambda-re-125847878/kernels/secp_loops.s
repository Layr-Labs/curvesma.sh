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
	.macro Q m,v
	\m 0\v,1\v
	\m 2\v,3\v
	\m 4\v,5\v
	\m 6\v,7\v
	.endm
	.macro E8 m,v
	\m 0\v
	\m 1\v
	\m 2\v
	\m 3\v
	\m 4\v
	\m 5\v
	\m 6\v
	\m 7\v
	.endm
	.macro H m,w,v
	\w \m,0\v
	\w \m,1\v
	\w \m,2\v
	\w \m,3\v
	.endm
	.macro G m,w
	H \m,\w,0
	H \m,\w,1
	H \m,\w,2
	H \m,\w,3
	H \m,\w,4
	H \m,\w,5
	H \m,\w,6
	H \m,\w,7
	.endm
	.macro W1 k,p
	ldp x10,x11,[\p,#0]
	ldp x12,x13,[\p,#16]
	ldp x16,x17,[x1,#16]
	ldp x14,x15,[x1],#64
	FS x10,x10,x14,x11,x11,x15,x12,x12,x16,x13,x13,x17,x17,1_0_\k
	stp x12,x13,[x3,#16]
	stp x10,x11,[x3],#64
	ldp x10,x11,[\p,#32]
	ldp x12,x13,[\p,#48]
	ldp x14,x15,[x1,#-32]
	ldp x16,x17,[x1,#-16]
	FS x10,x10,x14,x11,x11,x15,x12,x12,x16,x13,x13,x17,x17,1_1_\k
	stp x10,x11,[x3,#-32]
	stp x12,x13,[x3,#-16]
	.endm
	.macro S1 k
	SB x11,x12,x13,1_0_\k
	SB x11,x12,x13,1_1_\k
	.endm
	.macro P1 i,j
	ldp x9,x5,[x2],#16
	W1 \i,x9
	W1 \j,x5
	.endm
	.macro W2P k,o,a,b
	ldp x5,x6,[x3,#32]
	RF x5,x25,x26,x27,x28,x10,x11,x12,x13,x7,x8,x9
	ldp x5,x8,[x3],#\o
	RF x5,x25,x26,x27,x28,x9,x14,x15,x16,x17,x19,x20
	SH x5,x25,x28,7_\k
	RA x6,x5,x25,x26,x27,x19,x20,x21,x22,x23,x10,x11,x12,x13,x7,x24,xzr,x12,x13
	RA x8,x5,x25,x26,x27,x6,x19,x20,x21,x22,x9,x14,x15,x16,x17,x23,xzr,x15,x16
	SH x28,x5,x27,8_\k
	ldp x6,x8,[x3,#\a]
	RA x6,x28,x5,x25,x26,x19,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	ldp x6,x19,[x3,#\b]
	RA x6,x28,x5,x25,x26,x20,x21,x22,x29,x30,x9,x14,x15,x16,x17,x23,x23,x15,x16
	SH x27,x28,x26,9_\k
	RA x8,x27,x28,x5,x25,x6,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	RA x19,x27,x28,x5,x25,x6,x8,x20,x21,x22,x9,x14,x15,x16,x17,x23,x23,x27,x28
	FD x10,x11,x7,x24,x10,x11,x5,x6,0_\k
	FD x9,x14,x17,x23,x25,x26,x5,x6,1_\k
	.endm
	.macro W2 k,p
	ldp x15,x20,[x1,#0]
	ldp x14,x19,[x1,#16]
	OD x7,x10,x13,x24,x11,x12,x16,x23,x22,x29,x6
	SQ x9,x10,x21,x29,x11,x23,x7,x12,x24,x16,x13,x22,x6,x4,x17,x1,x8,\p
	b.cs .Lr2_4_\k
.Lb2_4_\k:
	FS x6,x9,x22,x8,x29,x16,x17,x23,x24,x30,x7,x21,x21,2_5_\k
	stp x9,x29,[x1,#-32]
	stp x23,x7,[x1,#-16]
	FM x7,x10,x6,x23,x29,x8,x9,x17,x21,x30,x10,x11,x24,x12,x13,x16,x6,x4,x30,x29,x9
	adds x7,x7,x15
	adcs x23,x23,x20
	adcs x29,x29,x14
	adcs x9,x9,x19
	adc x30,x30,xzr
	FO x19,x30,x14,x7,x7,x23,x23,2_6_\k
	stp x29,x9,[x1,#16]
	stp x7,x23,[x1],#-64
	.endm
	.macro SR k
	ST x5,x25,x26,x27,7_\k
	ST x28,x5,x25,x26,8_\k
	ST x27,x28,x5,x25,9_\k
	ST x10,x11,x12,x13,0_\k
	ST x25,x26,x27,x28,1_\k
	.endm
	.macro S2 k
	SR \k
	SC x23,x7,2_4_\k
	SB x8,x17,x30,2_5_\k
	SC x29,x9,2_6_\k
	.endm
	.macro SRB k
	SR b\k
	.endm
	.macro PB k
	W2P b\k,64,-16,-48
	stp x10,x11,[x3,#-32]
	stp x12,x13,[x3,#-16]
	.endm
	.macro P2 i,j
	W2P \i,-64,112,80
	ldp x0,x5,[x2,#-16]!
	W2 \i,x5
	W2P \j,-64,112,80
	W2 \j,x0
	.endm
	.macro CNT
	lsr x18,x0,#8
	and x0,x0,#255
	cbz x0,3f
	.endm
	.macro LP m,w
	subs x0,x0,#1
	b.ne 2b
3:
	cbz x18,4f
1:
	G \m,\w
	subs x18,x18,#1
	b.ne 1b
4:
	.endm
	.text
	.p2align 2
	.globl pass1
pass1:
	str x19,[sp,#-16]!
	mov x4,#977
	movk x4,#1,lsl #32
	CNT
2:
	ldr x9,[x2],#8
	W1 t,x9
	LP P1,Q
	ldr x19,[sp],#16
	ret
pass1_rare:
	S1 t
	G S1,E8
	.globl prod
prod:
	PR
	mov x25,#1
	mov x26,xzr
	mov x27,xzr
	mov x28,xzr
	mov x4,#977
	movk x4,#1,lsl #32
	CNT
2:
	PB t
	LP PB,E8
	stp x25,x26,[x3,#32]
	stp x27,x28,[x3,#48]
	EP
	ret
prod_rare:
	SRB t
	G SRB,E8
	.globl pass2
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
	CNT
2:
	W2P t,-64,112,80
	ldr x5,[x2,#-8]!
	W2 t,x5
	LP P2,Q
	EP
	ret
pass2_rare:
	S2 t
	G S2,E8

// Leaf kernels on 56-bit-VA (FEAT_LVA3 + FEAT_D128) M_G=1024 Y-cluster + Phase-0 birthday X-quads (MX=8, 6-window tables).
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
	.macro RA4 o,b0,b1,b2,b3,a0,a1,a2,a3,a4,a5,y,r,s,t,v
	mul \r,\o,\b0
	adds \a0,\a0,\r
	umulh \s,\o,\b0
	mul \r,\o,\b1
	adcs \s,\s,\r
	umulh \t,\o,\b1
	mul \r,\o,\b2
	adcs \t,\t,\r
	umulh \v,\o,\b2
	mul \r,\o,\b3
	adcs \v,\v,\r
	umulh \r,\o,\b3
	adc \r,\r,xzr
	adds \a1,\a1,\s
	adcs \a2,\a2,\t
	adcs \a3,\a3,\v
	adcs \a4,\a4,\r
	adc \a5,\y,xzr
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
	.macro SHX l,a,d,p,q
	mul \l,\d,x4
	umulh \d,\d,x4
	adds \a,\a,\d
	b.lo 9f
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	stadd x4,[sp]
9:
	.endm
	.macro SHXT a,d,p,q
	umulh \d,\d,x4
	adds \a,\a,\d
	b.lo 9f
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	stadd x4,[sp]
9:
	.endm
	.macro FDX e,f,i,j,o,p,l,h,y,z
	mul \l,\i,x4
	umulh \h,\i,x4
	madd \h,\j,x4,\h
	adds \o,\e,\l
	adcs \p,\f,\h
	b.lo 9f
	adcs \y,\y,xzr
	adc \z,\z,xzr
	stadd x4,[sp]
9:
	.endm
	.macro FS d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
	subs \d0,\a0,\b0
	sbcs \d1,\a1,\b1
	sbcs \d2,\a2,\b2
	sbcs \d3,\a3,\b3
	csel \t,x4,xzr,lo
	subs \d0,\d0,\t
	b.hs 9f
	sbcs \d1,\d1,xzr
	sbcs \d2,\d2,xzr
	sbc \d3,\d3,xzr
	stadd x4,[sp]
9:
	.endm
	.macro FADDW d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
	adds \d0,\a0,\b0
	adcs \d1,\a1,\b1
	adcs \d2,\a2,\b2
	adcs \d3,\a3,\b3
	csel \t,x4,xzr,hs
	adds \d0,\d0,\t
	b.lo 9f
	adcs \d1,\d1,xzr
	adcs \d2,\d2,xzr
	adc \d3,\d3,xzr
	stadd x4,[sp]
9:
	.endm
	.macro SQ13 a0,a1,a2,a3,w0,c1,c2,c3,c4,c5,c6,c7,t
	mul \c3,\a0,\a3
	mul \t,\a1,\a2
	adds \c3,\c3,\t
	umulh \c4,\a0,\a3
	umulh \t,\a1,\a2
	adcs \c4,\c4,\t
	umulh \c5,\a1,\a3
	adc \c5,\c5,xzr
	umulh \c2,\a0,\a1
	mul \t,\a0,\a2
	adds \c2,\c2,\t
	umulh \t,\a0,\a2
	adcs \c3,\c3,\t
	mul \t,\a1,\a3
	adcs \c4,\c4,\t
	mul \t,\a2,\a3
	adcs \c5,\c5,\t
	umulh \c6,\a2,\a3
	adc \c6,\c6,xzr
	mul \c1,\a0,\a1
	lsr \c7,\c6,#63
	extr \c6,\c6,\c5,#63
	extr \c5,\c5,\c4,#63
	extr \c4,\c4,\c3,#63
	extr \c3,\c3,\c2,#63
	extr \c2,\c2,\c1,#63
	mul \w0,\a0,\a0
	umulh \t,\a0,\a0
	adds \c1,\t,\c1,lsl #1
	mul \t,\a1,\a1
	adcs \c2,\c2,\t
	umulh \t,\a1,\a1
	adcs \c3,\c3,\t
	mul \t,\a2,\a2
	adcs \c4,\c4,\t
	umulh \t,\a2,\a2
	adcs \c5,\c5,\t
	mul \t,\a3,\a3
	adcs \c6,\c6,\t
	umulh \t,\a3,\a3
	adc \c7,\c7,\t
	mul \t,\c4,x4
	adds \w0,\w0,\t
	umulh \c4,\c4,x4
	mul \t,\c5,x4
	adcs \c4,\c4,\t
	umulh \c5,\c5,x4
	mul \t,\c6,x4
	adcs \c5,\c5,\t
	umulh \c6,\c6,x4
	mul \t,\c7,x4
	adcs \c6,\c6,\t
	umulh \c7,\c7,x4
	adc \c7,\c7,xzr
	adds \c1,\c1,\c4
	adcs \c2,\c2,\c5
	adcs \c3,\c3,\c6
	adc \c7,\c7,xzr
	.endm
	.macro F5 w0,w1,w2,w3,w4,o0,o1,ta,tb
	mul \ta,\w4,x4
	umulh \tb,\w4,x4
	adds \o0,\w0,\ta
	adcs \o1,\w1,\tb
	b.lo 9f
	adcs \w2,\w2,xzr
	adc \w3,\w3,xzr
	stadd x4,[sp]
9:
	.endm
	.macro LKFOLD4 o0,o1,o2,o3,c,t
	mul \t,\c,x4
	adds \o0,\o0,\t
	b.lo 9f
	adcs \o1,\o1,xzr
	adcs \o2,\o2,xzr
	adcs \o3,\o3,xzr
	csel \t,x4,xzr,hs
	add \o0,\o0,\t
	stadd x4,[sp]
9:
	.endm
	.macro LKP_ADD5 a0,a1,a2,a3,c,b0,b1,b2,b3
	adds \a0,\a0,\b0
	adcs \a1,\a1,\b1
	adcs \a2,\a2,\b2
	adcs \a3,\a3,\b3
	adc \c,\c,xzr
	.endm
	.macro LKP6_PAIR tb_off
	ldr x2,[x3,#\tb_off]
	ldp x28,x29,[x1,#32]
	add x5,x2,x28
	ld64b x16,[x5]
	add x5,x2,x29
	ld64b x8,[x5]
	adds x16,x16,x8
	adcs x17,x17,x9
	adcs x18,x18,x10
	adcs x19,x19,x11
	cset x28,hs
	adds x20,x20,x12
	adcs x21,x21,x13
	adcs x22,x22,x14
	adcs x23,x23,x15
	cset x29,hs
	ldr x5,[x1,#48]
	add x5,x2,x5
	ld64b x8,[x5]
	LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
	ldr x5,[x1,#56]
	add x5,x2,x5
	ld64b x8,[x5]
	LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
	add x5,x2,x6
	ld64b x8,[x5]
	LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
	add x2,x2,x7
	ld64b x8,[x2]
	LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
	ld64b x8,[x1]
	LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x29,x8,x9,x10,x11
	LKFOLD4 x16,x17,x18,x19,x28,x8
	LKFOLD4 x20,x21,x22,x23,x29,x8
	.endm
	.macro MUL4 a0,a1,a2,a3,b0,b1,b2,b3,o0,o1,o2,o3,e0,e1,e4,e5,t0,t1,t2,t3,sh
	RF \a0,\b0,\b1,\b2,\b3,\e0,\e1,\o2,\o3,\e4,\t0,\t1
	SHX \sh,\b0,\b3,\b1,\b2
	RA4 \a1,\sh,\b0,\b1,\b2,\e0,\e1,\o2,\o3,\e4,\e5,xzr,\t0,\t1,\t2,\t3
	SHX \b3,\sh,\b2,\b0,\b1
	RA4 \a2,\b3,\sh,\b0,\b1,\e0,\e1,\o2,\o3,\e4,\e5,\e5,\t0,\t1,\t2,\t3
	SHX \b2,\b3,\b1,\sh,\b0
	RA4 \a3,\b2,\b3,\sh,\b0,\e0,\e1,\o2,\o3,\e4,\e5,\e5,\t0,\t1,\t2,\t3
	FDX \e0,\e1,\e4,\e5,\o0,\o1,\t0,\e4,\o2,\o3
	.endm
	.macro P1X_A
	RF x20,x24,x25,x26,x27,x12,x13,x14,x15,x11,x28,x29
	SHX x5,x24,x27,x25,x26
	RA x21,x5,x24,x25,x26,x20,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,xzr,x14,x15
	SHX x27,x5,x26,x24,x25
	RA x22,x27,x5,x24,x25,x20,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,x29,x14,x15
	SHX x26,x27,x25,x5,x24
	RA x23,x26,x27,x5,x24,x20,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,x29,x22,x23
	FDX x12,x13,x11,x29,x20,x21,x5,x28,x22,x23
	.endm
	.macro P1X_B
	RF x16,x20,x21,x22,x23,x12,x13,x14,x15,x11,x28,x29
	SHX x5,x20,x23,x21,x22
	RA x17,x5,x20,x21,x22,x16,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,xzr,x14,x15
	SHX x23,x5,x22,x20,x21
	RA x18,x23,x5,x20,x21,x16,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,x29,x14,x15
	SHX x22,x23,x21,x5,x20
	RA x19,x22,x23,x5,x20,x16,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,x29,x26,x27
	FDX x12,x13,x11,x29,x24,x25,x5,x28,x26,x27
	.endm
	.macro PAIR7C k0,k1,k2,k3,s0,s1,s2,s3
	ld64b x16,[x1]
	RF \k0,x0,x2,x14,x15,x25,x26,x27,x28,x12,x10,x11
	RF \s0,x0,x2,x14,x15,x6,x7,x8,x9,x29,x10,x11
	SHX x5,x0,x15,x2,x14
	RA4 \k1,x5,x0,x2,x14,x25,x26,x27,x28,x12,x13,xzr,x10,x11,x24,x16
	RA4 \s1,x5,x0,x2,x14,x6,x7,x8,x9,x29,x30,xzr,x10,x11,x24,x16
	SHX x15,x5,x14,x0,x2
	RA4 \k2,x15,x5,x0,x2,x25,x26,x27,x28,x12,x13,x13,x10,x11,x24,x16
	RA4 \s2,x15,x5,x0,x2,x6,x7,x8,x9,x29,x30,x30,x10,x11,x24,x16
	SHX x14,x15,x2,x5,x0
	RA4 \k3,x14,x15,x5,x0,x25,x26,x27,x28,x12,x13,x13,x10,x11,x24,x16
	RA4 \s3,x14,x15,x5,x0,x6,x7,x8,x9,x29,x30,x30,x10,x11,x24,x16
	FDX x25,x26,x12,x13,x25,x26,x10,x12,x27,x28
	FDX x6,x7,x29,x30,x6,x7,x10,x29,x8,x9
	.endm
	.macro C1BSPLIT_BDAY p0,p1,p2,p3
	ldp x16,x17,[sp,#288]
	ldp x18,x19,[sp,#304]
	ldp x21,x22,[x3,#224]
	ldp x23,x24,[x3,#240]
	FS x21,x22,x23,x24,x16,x17,x18,x19,x21,x22,x23,x24,x5
	MUL4 \p0,\p1,\p2,\p3,x21,x22,x23,x24,x0,x2,x5,x6,x0,x2,x16,x17,x9,x10,x11,x12,x13
	stp x0,x2,[x3,#-64]
	stp x5,x6,[x3,#-48]
	ldr x24,[x3,#304]
	RF x24,\p0,\p1,\p2,\p3,x16,x17,x18,x19,x20,x21,x22
	F5 x16,x17,x18,x19,x20,x16,x17,x21,x22
	FADDW x8,x9,x10,x11,x0,x2,x5,x6,x16,x17,x18,x19,x20
	.endm
	.macro LK6P_CORE tb_off
	ldr x0,[x2,#\tb_off]
	add x1,x0,x5
	ld64b x16,[x1]
	add x1,x0,x6
	ld64b x8,[x1]
	adds x16,x16,x8
	adcs x17,x17,x9
	adcs x18,x18,x10
	adcs x19,x19,x11
	cset x24,hs
	adds x20,x20,x12
	adcs x21,x21,x13
	adcs x22,x22,x14
	adcs x23,x23,x15
	cset x25,hs
	add x1,x0,x7
	ld64b x8,[x1]
	LKP_ADD5 x16,x17,x18,x19,x24,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x25,x12,x13,x14,x15
	add x1,x0,x26
	ld64b x8,[x1]
	LKP_ADD5 x16,x17,x18,x19,x24,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x25,x12,x13,x14,x15
	add x1,x0,x27
	ld64b x8,[x1]
	LKP_ADD5 x16,x17,x18,x19,x24,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x25,x12,x13,x14,x15
	add x1,x0,x28
	ld64b x8,[x1]
	LKP_ADD5 x16,x17,x18,x19,x24,x8,x9,x10,x11
	LKP_ADD5 x20,x21,x22,x23,x25,x12,x13,x14,x15
	.endm
	.macro EDELTA_G e0,e1,e2,e3,dreg,b0,b1,b2,b3,o0,o1,o2,o3,p0,p1,p2,p3,p4,t0,t1
	RF \dreg,\b0,\b1,\b2,\b3,\p0,\p1,\p2,\p3,\p4,\t0,\t1
	F5 \p0,\p1,\p2,\p3,\p4,\p0,\p1,\t0,\t1
	FADDW \o0,\o1,\o2,\o3,\e0,\e1,\e2,\e3,\p0,\p1,\p2,\p3,\p4
	.endm
	.macro EQUAD6_ONE step_ins
	LK6P_CORE 312
	add x1,sp,#256
	ld64b x8,[x1]
	LKP_ADD5 x16,x17,x18,x19,x24,x8,x9,x10,x11
	adds x8,x20,x8
	adcs x9,x21,x9
	adcs x10,x22,x10
	adcs x11,x23,x11
	adc x25,x25,xzr
	LKFOLD4 x16,x17,x18,x19,x24,x12
	LKFOLD4 x8,x9,x10,x11,x25,x12
	ldr x24,[x2,#320]
	ldr x1,[sp,#376]
	ldp x14,x15,[sp,#496]
	EDELTA_G x16,x17,x18,x19,x24,x30,x1,x14,x15,x20,x21,x22,x23,x20,x21,x22,x23,x12,x13,x25
	st64b x16,[x2]
	ldr x25,[x2,#328]
	add x2,x2,#1024
	EDELTA_G x8,x9,x10,x11,x25,x30,x1,x14,x15,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22
	st64b x8,[x2]
	\step_ins
	.endm
	.macro FGPAIR6_ONE tb_off,step_ins
	LK6P_CORE \tb_off
	LKFOLD4 x16,x17,x18,x19,x24,x8
	LKFOLD4 x20,x21,x22,x23,x25,x8
	st64b x16,[x2]
	\step_ins
	.endm
	.macro STEP1024
	add x2,x2,#1024
	.endm
	.macro REWIND8_E2F
	sub x2,x3,#3520
	.endm
	.macro REWIND8_F2G
	sub x2,x3,#3456
	.endm
	.macro NOSTEP
	.endm
	.macro EFG_QUAD_2
	stp x0,x2,[sp,#8]
	stp x14,x15,[sp,#24]
	sub x2,x3,#3584
	ldp x5,x6,[sp,#384]
	ldp x7,x26,[sp,#400]
	ldp x27,x28,[sp,#416]
	lsr x30,x5,#6
	lsr x1,x6,#6
	bfi x30,x1,#43,#21
	EQUAD6_ONE STEP1024
	EQUAD6_ONE REWIND8_E2F
	ldp x5,x6,[sp,#96]
	ldp x7,x26,[sp,#112]
	ldp x27,x28,[sp,#352]
	FGPAIR6_ONE 272,STEP1024
	FGPAIR6_ONE 272,STEP1024
	FGPAIR6_ONE 272,STEP1024
	FGPAIR6_ONE 272,REWIND8_F2G
	ldp x5,x6,[sp,#440]
	ldp x7,x26,[sp,#456]
	ldp x27,x28,[sp,#472]
	FGPAIR6_ONE 216,STEP1024
	FGPAIR6_ONE 216,STEP1024
	FGPAIR6_ONE 216,STEP1024
	FGPAIR6_ONE 216,NOSTEP
	ldp x0,x2,[sp,#8]
	ldp x14,x15,[sp,#24]
	.endm
	.macro P1W6_ACC
	adds x8,x8,x18
	adcs x9,x9,x19
	adcs x10,x10,x20
	adcs x11,x11,x21
	adc x0,x0,xzr
	adds x26,x26,x22
	adcs x27,x27,x23
	adcs x28,x28,x24
	adcs x29,x29,x25
	adc x2,x2,xzr
	.endm
	.macro P1_6WIN_BFI
	SQ13 x8,x9,x10,x11,x13,x14,x15,x16,x5,x6,x7,x12,x17
	F5 x13,x14,x15,x16,x12,x13,x14,x17,x12
	ldp x6,x7,[x3,#256]
	ldp x12,x17,[x3,#272]
	ldp x1,x30,[x3,#288]
	bfi x6,x13,#6,#43
	ld64b x22,[x6]
	adds x8,x8,x22
	adcs x9,x9,x23
	adcs x10,x10,x24
	adcs x11,x11,x25
	cset x0,hs
	extr x5,x14,x13,#43
	bfi x7,x5,#6,#43
	ld64b x18,[x7]
	adds x8,x8,x18
	adcs x9,x9,x19
	adcs x10,x10,x20
	adcs x11,x11,x21
	adc x0,x0,xzr
	adds x26,x26,x22
	adcs x27,x27,x23
	adcs x28,x28,x24
	adcs x29,x29,x25
	cset x2,hs
	lsr x5,x14,#22
	add x5,x12,x5,lsl #6
	ld64b x18,[x5]
	P1W6_ACC
	bfi x17,x15,#6,#43
	ld64b x18,[x17]
	P1W6_ACC
	extr x5,x16,x15,#43
	bfi x1,x5,#6,#43
	ld64b x18,[x1]
	P1W6_ACC
	lsr x5,x16,#22
	add x5,x30,x5,lsl #6
	ld64b x18,[x5]
	P1W6_ACC
	LKFOLD4 x8,x9,x10,x11,x0,x18
	LKFOLD4 x26,x27,x28,x29,x2,x18
	add x5,sp,#128
	ld64b x18,[x5]
	RF x26,x18,x19,x20,x21,x13,x14,x15,x16,x17,x6,x7
	RA4 x27,x22,x23,x24,x25,x13,x14,x15,x16,x17,x12,xzr,x6,x7,x30,x0
	add x5,x5,#64
	ld64b x18,[x5]
	RA4 x28,x18,x19,x20,x21,x13,x14,x15,x16,x17,x12,x12,x6,x7,x30,x0
	RA4 x29,x22,x23,x24,x25,x13,x14,x15,x16,x17,x12,x12,x6,x7,x30,x0
	FDX x13,x14,x17,x12,x13,x14,x6,x17,x15,x16
	FS x28,x29,x30,x1,x8,x9,x10,x11,x13,x14,x15,x16,x6
	FADDW x0,x2,x5,x6,x8,x9,x10,x11,x13,x14,x15,x16,x20
	.endm
	.macro BLK6_REST s1,s2,s3,o0,o1,o2,o3
	extr x10,\s2,\s1,#22
	orr x10,x8,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	adc x11,x11,xzr
	orr x10,x8,\s2,lsl #7
	ld64b x12,[x10]
	adds \o0,\o0,x12
	adcs \o1,\o1,x13
	adcs \o2,\o2,x14
	adcs \o3,\o3,x15
	adc x11,x11,xzr
	extr x10,\s3,\s2,#44
	orr x10,x9,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x12
	adcs \o1,\o1,x13
	adcs \o2,\o2,x14
	adcs \o3,\o3,x15
	adc x11,x11,xzr
	lsr x10,\s3,#23
	orr x10,x9,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	adc x11,x11,xzr
	LKFOLD4 \o0,\o1,\o2,\o3,x11,x10
	.endm
	.macro PAIRB6
	mov x7,#0xfff8000000000000
	orr x8,x7,#64
	orr x9,x7,#128
	orr x10,x7,x28,lsl #8
	ld64b x20,[x10]
	extr x10,x29,x28,#43
	orr x10,x7,x10,lsl #8
	ld64b x12,[x10]
	adds x24,x20,x16
	adcs x25,x21,x17
	adcs x26,x22,x18
	adcs x27,x23,x19
	cset x11,hs
	BLK6_REST x29,x30,x1,x24,x25,x26,x27
	orr x10,x7,x0,lsl #8
	ld64b x12,[x10]
	extr x10,x2,x0,#43
	orr x10,x7,x10,lsl #8
	ld64b x16,[x10]
	adds x20,x12,x20
	adcs x21,x13,x21
	adcs x22,x14,x22
	adcs x23,x15,x23
	cset x11,hs
	BLK6_REST x2,x5,x6,x20,x21,x22,x23
	.endm
	.macro TR0 m,f,s,r3,r2,r1,t
	umulh \f,\m,\r2
	madd \f,\m,\r3,\f
	mul \t,\m,\r2
	lsr \s,\t,#3
	umulh \t,\m,\r1
	add \s,\s,\t,lsr #3
	.endm
	.macro TRK m,f,s,r3,r2,r1,t
	umulh \t,\m,\r2
	add \f,\f,\t
	madd \f,\m,\r3,\f
	mul \t,\m,\r2
	add \s,\s,\t,lsr #3
	umulh \t,\m,\r1
	add \s,\s,\t,lsr #3
	.endm
	.macro TROW7 op,up,bp,um,bm,r3,r2,r1,t
	\op \up,x13,x12,\r3,\r2,\r1,\t
	\op \bp,x15,x14,\r3,\r2,\r1,\t
	\op \um,x17,x16,\r3,\r2,\r1,\t
	\op \bm,x19,x18,\r3,\r2,\r1,\t
	.endm
	.macro TFP4_7
	TROW7 TR0,x28,x24,x0,x20,x11,x10,x9,x7
	SHX x7,x8,x11,x9,x10
	TROW7 TRK,x29,x25,x2,x21,x10,x9,x8,x11
	SHX x11,x7,x10,x8,x9
	TROW7 TRK,x30,x26,x5,x22,x9,x8,x7,x10
	SHXT x11,x9,x7,x8
	TROW7 TRK,x1,x27,x6,x23,x8,x7,x11,x9
	.endm
	.macro TFINC f,s,cl,ch
	add \f,\f,\s,lsr #61
	adds \s,\cl,\s,lsl #3
	adc \f,\f,\ch
	.endm
	.macro TFINS f,s
	add \f,\f,\s,lsr #61
	.endm
	.macro BWDONE6_J e0,e1,f0,f1,g0,g1
	P1_6WIN_BFI
	PAIRB6
	ldp x8,x9,[x3,#\e0]
	ldp x10,x11,[x3,#\e1]
	TFP4_7
	ldp x7,x8,[sp,#320]
	ldp x9,x10,[sp,#336]
	TFINC x13,x12,x7,x8
	TFINC x15,x14,x9,x10
	TFINC x17,x16,x7,x8
	TFINC x19,x18,x9,x10
	st64b x12,[x3]
	ldp x8,x9,[x3,#\f0]
	ldp x10,x11,[x3,#\f1]
	TFP4_7
	ldp x7,x8,[x3,#192]
	ldp x9,x10,[x3,#208]
	TFINC x13,x12,x7,x8
	TFINC x15,x14,x9,x10
	TFINC x17,x16,x7,x8
	TFINC x19,x18,x9,x10
	add x11,x3,#64
	st64b x12,[x11]
	ldp x8,x9,[x3,#\g0]
	ldp x10,x11,[x3,#\g1]
	TFP4_7
	TFINS x13,x12
	TFINS x15,x14
	TFINS x17,x16
	TFINS x19,x18
	add x11,x3,#128
	st64b x12,[x11]
	.endm
	.macro BWDONE6_J0
	P1_6WIN_BFI
	PAIRB6
	ld64b x8,[x3]
	TFP4_7
	ldp x7,x8,[sp,#320]
	ldp x9,x10,[sp,#336]
	TFINC x13,x12,x7,x8
	TFINC x15,x14,x9,x10
	TFINC x17,x16,x7,x8
	TFINC x19,x18,x9,x10
	st64b x12,[x3]
	add x3,x3,#64
	ld64b x8,[x3]
	TFP4_7
	ldp x7,x8,[x3,#128]
	ldp x9,x10,[x3,#144]
	TFINC x13,x12,x7,x8
	TFINC x15,x14,x9,x10
	TFINC x17,x16,x7,x8
	TFINC x19,x18,x9,x10
	st64b x12,[x3]
	add x3,x3,#64
	ld64b x8,[x3]
	TFP4_7
	TFINS x13,x12
	TFINS x15,x14
	TFINS x17,x16
	TFINS x19,x18
	st64b x12,[x3]
	.endm
	.macro BWDPAIR p0,p1,p2,p3
	C1BSPLIT_BDAY \p0,\p1,\p2,\p3
	BWDONE6_J -480,-464,-416,-400,-352,-336
	ldp x8,x9,[x3,#-64]
	ldp x10,x11,[x3,#-48]
	sub x3,x3,#512
	BWDONE6_J0
	.endm
	.macro FLUSH7 s,r
	sub \s,x5,#2016
	str \s,[x0]
	ldr \r,[\s,#-16]
	add \r,\r,#1
	mov x5,#2048
	stp x5,\r,[\s,#-24]
	stlur x0,[\s,#-32]
	dmb ish
	.endm
	.macro LJ0_7 k,w0l,w0h,w1l,w1h
	adds x0,\w0l,\w1l
	adc x5,\w0h,\w1h
	cmp \w0h,\w1h
	b.hs kh\k
	cmn \w0h,x5
	b.hs kt\k
	and x0,x1,\w0l,lsr #20
	ldadd x4,x5,[x0]
	stp \w0h,\w1h,[x5]
	stp x2,x26,[x5,#16]
	tbz x5,#11,kn\k
	FLUSH7 \w0h,\w1h
kn\k:
	.endm
	.macro JHB7 k,w0l,w0h,w1l,w1h
kh\k:
	cmn \w1h,x5
	b.hs kt\k
	mvn \w0h,x5
	and x0,x1,\w1l,lsr #20
	ldadd x4,x5,[x0]
	stp \w1h,\w0h,[x5]
	stp x2,x26,[x5,#16]
	tbz x5,#11,km\k
	FLUSH7 \w0h,\w1h
km\k:
	.endm
	.macro JTB7 k,w0l,w0h,w1l,w1h
kt\k:
	mvn \w1h,x5
	bic x0,x1,x0,lsr #20
	ldadd x4,x5,[x0]
	stp \w1h,\w0h,[x5]
	stp x2,x26,[x5,#16]
	tbz x5,#11,ku\k
	FLUSH7 \w0h,\w1h
ku\k:
	.endm
	.macro S1_7 a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h,g1l,g1h
	adds x27,\a0l,\b0l
	adc x28,\a0h,\b0h
	adds x27,x27,\g0l,lsl #3
	adc x28,x28,\g0h
	adds x29,\a1l,\b1l
	adc x30,\a1h,\b1h
	adds x29,x29,\g1l,lsl #3
	adc x30,x30,\g1h
	.endm
	.macro S2_7 a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h,g1l,g1h
	subs x27,\b0l,\g0l,lsl #3
	sbc x28,\b0h,\g0h
	subs x29,\b1l,\g1l,lsl #3
	sbc x30,\b1h,\g1h
	subs \a0l,\a0l,\b0l
	sbc \a0h,\a0h,\b0h
	subs \a1l,\a1l,\b1l
	sbc \a1h,\a1h,\b1h
	adds \b0l,x27,x29
	adc \b0h,x28,x30
	adds \b1l,\a0l,\b0l
	adc \b1h,\a0h,\b0h
	subs \g0l,\a1l,x27
	sbc \g0h,\a1h,x28
	.endm
	.macro S3_7 a0l,a0h,a1l,a1h,b0l,b0h
	subs \a0l,\a0l,x29
	sbc \a0h,\a0h,x30
	adds \a1l,\a1l,\b0l
	adc \a1h,\a1h,\b0h
	.endm
	.macro TSUM7 k,a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h,g1l,g1h
	S1_7 \a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h,\g1l,\g1h
	LJ0_7 \k,x27,x28,x29,x30
	S2_7 \a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h,\g1l,\g1h
	LJ0_7 d\k,\b1l,\b1h,\g0l,\g0h
	S3_7 \a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	LJ0_7 e\k,\a0l,\a0h,\a1l,\a1h
	.endm
	.macro C3_7 k,j,a0l,a0h,a1l,a1h,b0l,b0h
	S3_7 \a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	LJ0_7 \k,\a0l,\a0h,\a1l,\a1h
	b \j
	JHB7 \k,\a0l,\a0h,\a1l,\a1h
	b \j
	JTB7 \k,\a0l,\a0h,\a1l,\a1h
	b \j
	.endm
	.macro C2_7 k,j,a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h,g1l,g1h
	S2_7 \a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h,\g1l,\g1h
	LJ0_7 \k,\b1l,\b1h,\g0l,\g0h
	C3_7 a\k,\j,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	JHB7 \k,\b1l,\b1h,\g0l,\g0h
	C3_7 b\k,\j,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	JTB7 \k,\b1l,\b1h,\g0l,\g0h
	C3_7 c\k,\j,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	.endm
	.macro TJX7 k,a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h,g1l,g1h
	JHB7 e\k,\a0l,\a0h,\a1l,\a1h
	b kne\k
	JTB7 e\k,\a0l,\a0h,\a1l,\a1h
	b kne\k
	JHB7 d\k,\b1l,\b1h,\g0l,\g0h
	C3_7 g\k,kne\k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	JTB7 d\k,\b1l,\b1h,\g0l,\g0h
	C3_7 i\k,kne\k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h
	JHB7 \k,x27,x28,x29,x30
	C2_7 l\k,kne\k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h,\g1l,\g1h
	JTB7 \k,x27,x28,x29,x30
	C2_7 o\k,kne\k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h,\g1l,\g1h
	.endm
	.macro KEYPAIR7 kp0,km0,kp1,km1
	ld64b x6,[x3]
	add x0,x3,#64
	ld64b x14,[x0]
	ldp x22,x23,[x3,#128]
	ldp x24,x25,[x3,#144]
	TSUM7 \kp0,x6,x7,x8,x9,x14,x15,x16,x17,x22,x23,x24,x25
	ldp x22,x23,[x3,#160]
	ldp x24,x25,[x3,#176]
	TSUM7 \km0,x10,x11,x12,x13,x18,x19,x20,x21,x22,x23,x24,x25
	add x26,x26,#4
	add x3,x3,#512
	ld64b x6,[x3]
	add x0,x3,#64
	ld64b x14,[x0]
	ldp x22,x23,[x3,#128]
	ldp x24,x25,[x3,#144]
	TSUM7 \kp1,x6,x7,x8,x9,x14,x15,x16,x17,x22,x23,x24,x25
	ldp x22,x23,[x3,#160]
	ldp x24,x25,[x3,#176]
	TSUM7 \km1,x10,x11,x12,x13,x18,x19,x20,x21,x22,x23,x24,x25
	add x26,x26,#4
	add x3,x3,#512
	.endm
	.macro FWDP6_PAIR tb_off,eo0,eo1
	LKP6_PAIR \tb_off
	add x2,x3,#\eo0
	st64b x20,[x2]
	P1X_A
	add x2,x3,#\eo1
	st64b x16,[x2]
	P1X_B
	.endm
	.macro FWD8_ONE
	ldp x6,x7,[x1,#288]
	LKP6_PAIR 416
	st64b x20,[x3]
	P1X_A
	add x2,x3,#1024
	st64b x16,[x2]
	P1X_B
	FWDP6_PAIR 424,2048,3072
	add x1,x1,#512
	add x3,x2,#1024
	.endm
	.macro BWDP8
	sub x1,x3,#576
	PAIR7C x20,x21,x22,x23,x16,x17,x18,x19
	st64b x6,[x1]
	BWDPAIR x25,x26,x27,x28
	ldp x0,x2,[x3,#-192]
	ldp x14,x15,[x3,#-176]
	sub x3,x3,#640
	.endm
	.macro BWD8_MID
	EFG_QUAD_2
	BWDP8
	BWDP8
	BWDP8
	BWDP8
	ldr x6,[sp]
	stp x6,xzr,[sp,#-512]!
	.endm
	.text
	.p2align 2
	.globl fwd
fwd:
	mov x4,#977
	movk x4,#1,lsl #32
	stp x29,x30,[sp,#-16]!
	str xzr,[sp,#-16]!
1:
	FWD8_ONE
	FWD8_ONE
	FWD8_ONE
	FWD8_ONE
	subs x0,x0,#32
	b.ne 1b
	ldr x0,[sp],#16
	ldp x29,x30,[sp],#16
	ret
	.globl bwd
bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	mov x5,sp
	stp x5,x30,[x1]
	str xzr,[sp]
	sub x3,x3,#448
1:
	BWD8_MID
	EFG_QUAD_2
	BWDP8
	BWDP8
	BWDP8
	sub x1,x3,#576
	PAIR7C x20,x21,x22,x23,x16,x17,x18,x19
	st64b x6,[x1]
	BWDPAIR x25,x26,x27,x28
	ldp x0,x2,[x3,#-192]
	ldp x14,x15,[x3,#-176]
	mov x5,sp
	tst x5,#0x7fffffff
	b.eq 2f
	ldr x6,[sp]
	stp x6,xzr,[sp,#-512]!
	sub x3,x3,#640
	b 1b
2:
	ldr x6,[sp]
	ldp x5,x30,[sp,#-32]
	mov sp,x5
	mov x5,x6
	ret
	.globl keys
keys:
	str x30,[sp,#-16]!
	mov x1,#0xffffffffff8
	mov x4,#32
1:
	KEYPAIR7 p00,m00,p10,m10
	KEYPAIR7 p20,m20,p30,m30
	tbz x3,#34,1b
	ldr x30,[sp],#16
	ret
	TJX7 p00,x6,x7,x8,x9,x14,x15,x16,x17,x22,x23,x24,x25
	TJX7 m00,x10,x11,x12,x13,x18,x19,x20,x21,x22,x23,x24,x25
	TJX7 p10,x6,x7,x8,x9,x14,x15,x16,x17,x22,x23,x24,x25
	TJX7 m10,x10,x11,x12,x13,x18,x19,x20,x21,x22,x23,x24,x25
	TJX7 p20,x6,x7,x8,x9,x14,x15,x16,x17,x22,x23,x24,x25
	TJX7 m20,x10,x11,x12,x13,x18,x19,x20,x21,x22,x23,x24,x25
	TJX7 p30,x6,x7,x8,x9,x14,x15,x16,x17,x22,x23,x24,x25
	TJX7 m30,x10,x11,x12,x13,x18,x19,x20,x21,x22,x23,x24,x25

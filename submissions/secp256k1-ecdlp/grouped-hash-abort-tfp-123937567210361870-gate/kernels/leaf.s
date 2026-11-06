// Leaf kernels (proof.md 3). fwd: e of each X pair by lookups (EFWD), ST64B of e, acc, acc *= e (P1X). bwd per pair
// (BWDP): W2S, C1BS split, two INC (C2UV, u^2, v^2, PAIRB, E/F/G by lookups, truncated products TFP4X). keys as
// 5aacf299 (3.2, 3.3). x4 = C0; SHX, SHXT, FDX, F5, FS, FADDW, BLK0, LF branch over their fixup (3.4).

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

	.macro P1X
	RF x20,x24,x25,x26,x27,x13,x14,x15,x16,x17,x28,x19
	SHX x5,x24,x27,x25,x26
	RA x21,x5,x24,x25,x26,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,xzr,x15,x16
	SHX x27,x5,x26,x24,x25
	RA x22,x27,x5,x24,x25,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,x29,x15,x16
	SHX x26,x27,x25,x5,x24
	RA x23,x26,x27,x5,x24,x6,x7,x8,x20,x21,x13,x14,x15,x16,x17,x29,x29,x26,x27
	FDX x13,x14,x17,x29,x24,x25,x5,x6,x26,x27
	.endm


	// count by C0 (BLK-type site: the rare fixup finishes limbs 2, 3 and adds C0 on their carry)
	.macro AC4F a0,a1,a2,a3,b0,b1,b2,b3,c
	adds \a0,\a0,\b0
	adcs \a1,\a1,\b1
	adcs \a2,\a2,\b2
	adcs \a3,\a3,\b3
	cset \c,hs
	.endm
	.macro AC4 a0,a1,a2,a3,b0,b1,b2,b3,c
	adds \a0,\a0,\b0
	adcs \a1,\a1,\b1
	adcs \a2,\a2,\b2
	adcs \a3,\a3,\b3
	adc \c,\c,xzr
	.endm
	.macro LF a0,a1,a2,a3,c,t
	mul \t,\c,x4
	adds \a0,\a0,\t
	adcs \a1,\a1,xzr
	b.lo 9f
	adcs \a2,\a2,xzr
	adcs \a3,\a3,xzr
	csel \t,x4,xzr,hs
	add \a0,\a0,\t
	stadd x4,[sp]
9:
	.endm
	// two windows of one 64-bit word v by BFI into the line base b (row t of window 2h at bits 8..39), one LD64B each
	.macro EFWD
	ldp x7,x8,[x3,#-280]
	ldp x20,x21,[x1,#48]
	ldp x22,x23,[x1,#64]
	ldp x8,x9,[x1,#336]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4F x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	ldp x8,x9,[x1,#352]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	ldp x8,x9,[x1,#368]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	ldp x8,x9,[x1,#384]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	LF x20,x21,x22,x23,x11,x10
	.endm
	.macro PFWD
	EFWD
	st64b x20,[x3]
	add x3,x3,#1024
	P1X
	.endm




	// beta*w by six 43-bit windows: entry t of window i at x15 + 256t + 32i holds beta*t*2^(43i) mod p (proof.md 3.1)
	.macro BLK0 s0,s1,s2,s3,o0,o1,o2,o3
	orr x10,x9,\s0,lsl #8
	ldp \o0,\o1,[x10]
	ldp \o2,\o3,[x10,#16]
	extr x10,\s1,\s0,#43
	orr x10,x9,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	cset x11,hs
	extr x10,\s2,\s1,#22
	orr x10,x5,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x12
	adcs \o1,\o1,x13
	adcs \o2,\o2,x14
	adcs \o3,\o3,x15
	adc x11,x11,xzr
	lsr x10,\s2,#1
	orr x10,x5,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	adc x11,x11,xzr
	extr x10,\s3,\s2,#44
	orr x10,x6,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x12
	adcs \o1,\o1,x13
	adcs \o2,\o2,x14
	adcs \o3,\o3,x15
	adc x11,x11,xzr
	lsr x10,\s3,#23
	orr x10,x6,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	adc x11,x11,xzr
	mul x10,x11,x4
	adds \o0,\o0,x10
	adcs \o1,\o1,xzr
	b.lo 9f
	adcs \o2,\o2,xzr
	adcs \o3,\o3,xzr
	csel x10,x4,xzr,hs
	add \o0,\o0,x10
	stadd x4,[sp]
9:
	.endm
	.macro PAIRB
	mov x9,#0xfff8000000000000
	orr x5,x9,#64
	orr x6,x9,#128
	BLK0 x21,x22,x23,x24,x29,x27,x7,x1
	BLK0 x25,x26,x0,x2,x30,x28,x8,x20
	.endm


	// limb m times row words r3..r1: f += lo(m r3)+hi(m r2), s += (lo(m r2)>>3)+(hi(m r1)>>3)
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


	// last row shift: words 1..3 only (no product reads word 0)
	.macro SHXT a,d,p,q
	umulh \d,\d,x4
	adds \a,\a,\d
	b.lo 9f
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	stadd x4,[sp]
9:
	.endm



	.macro TFINC f,s,cl,ch
	add \f,\f,\s,lsr #61
	adds \s,\cl,\s,lsl #3
	adc \f,\f,\ch
	.endm

	.macro TFINS f,s
	add \f,\f,\s,lsr #61
	.endm


	.macro W2S
	RF x20,x25,x26,x27,x28,x9,x10,x11,x12,x13,x5,x6
	RF x16,x25,x26,x27,x28,x0,x1,x2,x24,x29,x5,x6
	SHX x5,x25,x28,x26,x27
	RA4 x21,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
	RA4 x17,x5,x25,x26,x27,x0,x1,x2,x24,x29,x30,xzr,x6,x7,x8,x28
	SHX x6,x5,x27,x25,x26
	RA4 x22,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
	RA4 x18,x6,x5,x25,x26,x0,x1,x2,x24,x29,x30,x30,x7,x8,x28,x27
	SHX x7,x6,x26,x5,x25
	RA4 x23,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
	RA4 x19,x7,x6,x5,x25,x0,x1,x2,x24,x29,x30,x30,x8,x28,x27,x26
	FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12
	FDX x0,x1,x29,x30,x0,x1,x6,x29,x2,x24
	.endm
	.macro C1BS
	ldp x13,x14,[sp,#16]
	ldp x15,x16,[sp,#32]
	ldp x17,x18,[x3,#192]
	ldp x19,x20,[x3,#208]
	FS x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x21
	RF x9,x13,x14,x15,x16,x0,x1,x2,x5,x6,x7,x8
	SHX x17,x13,x16,x14,x15
	RA4 x10,x17,x13,x14,x15,x0,x1,x2,x5,x6,x18,xzr,x7,x8,x19,x20
	SHX x16,x17,x15,x13,x14
	RA4 x11,x16,x17,x13,x14,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
	SHX x15,x16,x14,x17,x13
	RA4 x12,x15,x16,x17,x13,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
	FDX x0,x1,x6,x18,x0,x1,x7,x8,x2,x5
	stp x0,x1,[x3,#-480]
	stp x2,x5,[x3,#-464]
	ldp x13,x14,[x3,#-128]
	ldp x15,x16,[x3,#-112]
	bfi x13,x9,#8,#32
	ld64b x20,[x13]
	AC4F x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x9,#32
	bfi x13,x7,#8,#32
	ld64b x20,[x13]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x14,x10,#8,#32
	ld64b x20,[x14]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x10,#32
	bfi x14,x7,#8,#32
	ld64b x20,[x14]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x15,x11,#8,#32
	ld64b x20,[x15]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x11,#32
	bfi x15,x7,#8,#32
	ld64b x20,[x15]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x16,x12,#8,#32
	ld64b x20,[x16]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x12,#32
	bfi x16,x7,#8,#32
	ld64b x20,[x16]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	LF x0,x1,x2,x5,x6,x7
	.endm
	.macro C2UV
	ldp x13,x14,[sp,#112]
	ldp x15,x16,[sp,#128]
	RF x0,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23
	ldp x13,x14,[sp,#144]
	ldp x15,x16,[sp,#160]
	RA4 x1,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,xzr,x22,x23,x25,x26
	ldp x13,x14,[sp,#176]
	ldp x15,x16,[sp,#192]
	RA4 x2,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,x24,x22,x23,x25,x26
	ldp x13,x14,[sp,#208]
	ldp x15,x16,[sp,#224]
	RA4 x5,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,x24,x22,x23,x25,x26
	FDX x17,x18,x21,x24,x17,x18,x22,x23,x19,x20
	ldp x9,x10,[x3,#352]
	ldp x11,x12,[x3,#368]
	bfi x9,x0,#8,#32
	ldp x13,x14,[x9]
	ldp x15,x16,[x9,#16]
	lsr x7,x0,#32
	bfi x9,x7,#8,#32
	ld64b x22,[x9]
	AC4F x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x10,x1,#8,#32
	ld64b x22,[x10]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x1,#32
	bfi x10,x7,#8,#32
	ld64b x22,[x10]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x11,x2,#8,#32
	ld64b x22,[x11]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x2,#32
	bfi x11,x7,#8,#32
	ld64b x22,[x11]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x12,x5,#8,#32
	ld64b x22,[x12]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x5,#32
	bfi x12,x7,#8,#32
	ld64b x22,[x12]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	LF x13,x14,x15,x16,x6,x7
	FADDW x27,x28,x29,x30,x17,x18,x19,x20,x13,x14,x15,x16,x5
	FS x17,x18,x19,x20,x17,x18,x19,x20,x13,x14,x15,x16,x5
	.endm
	.macro TROWX op,a,b,c,d,r3,r2,r1,t
	\op \a,x11,x12,\r3,\r2,\r1,\t
	\op \b,x13,x14,\r3,\r2,\r1,\t
	\op \c,x15,x16,\r3,\r2,\r1,\t
	\op \d,x17,x18,\r3,\r2,\r1,\t
	.endm
	.macro TFP4X
	TROWX TR0,x21,x25,x29,x30,x10,x9,x6,x19
	SHX x19,x5,x10,x6,x9
	TROWX TRK,x22,x26,x27,x28,x9,x6,x5,x10
	SHX x10,x19,x9,x5,x6
	TROWX TRK,x23,x0,x7,x8,x6,x5,x19,x9
	SHXT x10,x6,x19,x5
	TROWX TRK,x24,x2,x1,x20,x5,x19,x10,x6
	.endm
	.macro TOUTCX r,c,d,p,q,m,n
	ldp x9,x10,[\r,#\c]
	ldp x5,x6,[\r,#\d]
	TFINC x11,x12,x9,x10
	TFINC x13,x14,x9,x10
	TFINC x15,x16,x5,x6
	TFINC x17,x18,x5,x6
	stp x12,x11,[x3,#\p]
	stp x16,x15,[x3,#\q]
	stp x14,x13,[x3,#\m]
	stp x18,x17,[x3,#\n]
	.endm
	.macro LKE
	ldr x1,[x3,#224]
	ldp x5,x6,[sp,#80]
	ldp x9,x10,[sp,#96]
	ldp x0,x2,[sp,#272]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4F x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#288]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#304]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#320]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
	.macro LKF
	ldr x1,[x3,#224]
	ldp x0,x2,[sp,#336]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	AC4F x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#352]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#368]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#384]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
	.macro LKG
	ldr x1,[x3,#224]
	ldp x0,x2,[sp,#400]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	AC4F x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#416]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#432]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#448]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
	.macro INC
	C2UV
	SQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
	F5 x21,x22,x23,x24,x8,x21,x22,x9,x8
	SQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
	F5 x25,x26,x0,x2,x8,x25,x26,x9,x8
	PAIRB
	stp x0,x2,[sp,#464]
	str x1,[sp,#480]
	LKE
	TFP4X
	TOUTCX sp,240,256,320,336,128,144
	LKF
	TFP4X
	TOUTCX x3,160,176,288,304,96,112
	LKG
	TFP4X
	TFINS x11,x12
	TFINS x13,x14
	TFINS x15,x16
	TFINS x17,x18
	stp x12,x11,[x3,#256]
	stp x16,x15,[x3,#272]
	stp x14,x13,[x3,#64]
	stp x18,x17,[x3,#80]
	.endm
	.macro BWDP
	ld64b x16,[x3]
	W2S
	stp x0,x1,[x3,#-496]
	stp x2,x24,[x3,#-512]
	C1BS
	INC
	sub x3,x3,#512
	ldp x0,x1,[x3,#32]
	ldp x2,x5,[x3,#48]
	INC
	ldp x25,x26,[x3,#16]
	ldp x27,x28,[x3],#-512
	.endm


// After slot k (x16; x15 = 8*d1): LC back to slot 1, RC cursor (2^43 + 8*d1) += 32k (x0);
// the header sends the block iff the old cursor is below the slot end x12, else drops it.
// s, r: two leaf registers, dead once its record is stored.
	.macro FLUSH s,r
	sub \s,x16,#4080
	sub \s,\s,#4080
	str \s,[x15]
	orr \r,x15,#0x80000000000
	ldadd x0,x16,[\r]
	cmp x16,x12
	b.hs 6f
	str x0,[\s,#-16]
	stp x15,x16,[\s,#-32]
6:
	.endm

	// is (~sigma, ~S) on JH and JT: one MVN of sigma, ~S read by the emit's EON (JH) or BIC (JT) of x14
	.macro LJ0 k,w0l,w0h,w1l,w1h
	adds x14,\w0l,\w1l
	adc x15,\w0h,\w1h
	cmp \w0h,\w1h
	b.hs kh\k
	cmn \w0h,x15
	b.hs kt\k
	and x15,x1,\w0l,lsr #25
	ldadd x4,x16,[x15]
	add \w1l,x18,\w1l,lsr #28
	stp \w0h,\w1h,[x16]
	stp x17,\w1l,[x16,#16]
	tbz x16,#13,kn\k
	FLUSH \w0h,\w1h
kn\k:
	.endm

	// j=1 (JHB) and j=2 (JTB) blocks: N', emit, then on into copies of the sign's later leaves (C2, C3)
	.macro JHB k,w0l,w0h,w1l,w1h
kh\k:
	cmn \w1h,x15
	b.hs kt\k
	mvn \w0h,x15
	and x15,x1,\w1l,lsr #25
	ldadd x4,x16,[x15]
	eon \w0l,x18,x14,lsr #28
	stp \w1h,\w0h,[x16]
	stp x17,\w0l,[x16,#16]
	tbz x16,#13,km\k
	FLUSH \w0h,\w1h
km\k:
	.endm

	.macro JTB k,w0l,w0h,w1l,w1h
kt\k:
	mvn \w1h,x15
	bic x15,x1,x14,lsr #25
	ldadd x4,x16,[x15]
	add \w0l,x18,\w0l,lsr #28
	stp \w1h,\w0h,[x16]
	stp x17,\w0l,[x16,#16]
	tbz x16,#13,ku\k
	FLUSH \w0h,\w1h
ku\k:
	.endm

	// x24..x27), alpha x19,x28 and x29,x2 by LDP; x3 ends at entry+128 (sign +) or entry-256 (sign -)
	.macro TLDP
	ld64b x20,[x3]
	ldp x19,x28,[x3,#64]
	ldp x29,x2,[x3,#80]
	sub x3,x3,#128
	.endm

	.macro TLDM
	ldp x19,x28,[x3],#-64
	ldp x29,x2,[x3,#80]
	ld64b x20,[x3]
	sub x3,x3,#320
	.endm

	// W_m(t)=U_m(alpha)+U_(m+t)(beta')+U_(m+2t)(gamma), U_2=-U_0-U_1; per sign x18=x3<<30 (JH's w3 by EON with it)
	.macro S1
	lsl x18,x3,#30
	adds x5,x19,x24
	adc x6,x28,x25
	adds x5,x5,x20,lsl #3
	adc x6,x6,x21
	adds x7,x29,x26
	adc x8,x2,x27
	adds x7,x7,x22,lsl #3
	adc x8,x8,x23
	.endm

	// W(1) = (i3+i6, i5-i1), W(2) = (i3-i2, i5+i6)
	.macro S2
	subs x9,x24,x20,lsl #3
	sbc x10,x25,x21
	subs x11,x26,x22,lsl #3
	sbc x13,x27,x23
	subs x19,x19,x24
	sbc x28,x28,x25
	subs x29,x29,x26
	sbc x2,x2,x27
	adds x24,x9,x11
	adc x25,x10,x13
	adds x26,x19,x24
	adc x27,x28,x25
	subs x20,x29,x9
	sbc x21,x2,x10
	.endm

	.macro S3
	subs x19,x19,x11
	sbc x28,x28,x13
	adds x29,x29,x24
	adc x2,x2,x25
	.endm

	// of the sign's later leaves (C2, C3) and ends in one B to the sign's end (kne\k is the inline third leaf's join).
	.macro TSUM k
	S1
	LJ0 \k,x5,x6,x7,x8
	S2
	LJ0 d\k,x26,x27,x20,x21
	S3
	LJ0 e\k,x19,x28,x29,x2
	.endm

	.macro C3 k,j
	S3
	LJ0 \k,x19,x28,x29,x2
	b \j
	JHB \k,x19,x28,x29,x2
	b \j
	JTB \k,x19,x28,x29,x2
	b \j
	.endm

	.macro C2 k,j
	S2
	LJ0 \k,x26,x27,x20,x21
	C3 a\k,\j
	JHB \k,x26,x27,x20,x21
	C3 b\k,\j
	JTB \k,x26,x27,x20,x21
	C3 c\k,\j
	.endm

	.macro TJX k
	JHB e\k,x19,x28,x29,x2
	b kne\k
	JTB e\k,x19,x28,x29,x2
	b kne\k
	JHB d\k,x26,x27,x20,x21
	C3 g\k,kne\k
	JTB d\k,x26,x27,x20,x21
	C3 i\k,kne\k
	JHB \k,x5,x6,x7,x8
	C2 l\k,kne\k
	JTB \k,x5,x6,x7,x8
	C2 o\k,kne\k
	.endm


	.text
	.p2align 2
	.globl fwd
fwd:
	mov x4,#977
	movk x4,#1,lsl #32
	str xzr,[sp,#-16]!
1:
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	PFWD
	tst x3,#0x1fc00
	b.ne 2f
	add x1,x1,#512
2:
	subs x0,x0,#32
	b.ne 1b
	ldr x0,[sp],#16
	ret
	.globl bwd
bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	str x30,[sp,#8]
1:
	BWDP
	BWDP
	BWDP
	BWDP
	BWDP
	BWDP
	BWDP
	BWDP
	add x5,x3,#512
	tst x5,#0x1fe00
	b.ne 1b
	tbnz x3,#34,2f
	ldp x5,x6,[sp]
	sub sp,sp,#512
	stp x5,x6,[sp]
	b 1b
2:
	ldp x0,x30,[sp]
	ret
	.globl keys
keys:
	mov x1,#0x7ffffffff8
	mov x4,#32
	mov x0,#8192
1:
	TLDP
	TSUM p00
	TLDM
	TSUM m00
	TLDP
	TSUM p10
	TLDM
	TSUM m10
	tbz x3,#34,1b
	ret
	TJX p00
	TJX m00
	TJX p10
	TJX m10

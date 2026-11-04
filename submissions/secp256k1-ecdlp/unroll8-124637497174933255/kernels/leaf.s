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
	.macro AD4 a0,a1,a2,a3,a4,a5,k0,k1,k2,k3
	adds \a0,\a0,\k0
	adcs \a1,\a1,\k1
	adcs \a2,\a2,\k2
	adcs \a3,\a3,\k3
	adcs \a4,\a4,xzr
	adc \a5,\a5,xzr
	.endm
	.macro PAIR a0,a1,a2,a3,b0,b1,b2,b3,b0_,b1_,b2_,b3_,r0,r1,r2,r3,f0,f1,f2,f3,x0,x1,x2,x3,x4_,x5,y0,y1,y2,y3,y4,y5
	RF \a0,\r0,\r1,\r2,\r3,\x0,\x1,\x2,\x3,\x4_,\f0,\f1
	RF \b0,\r0,\r1,\r2,\r3,\y0,\y1,\y2,\y3,\y4,\f0,\f1
	SHX \f0,\r0,\r3,\r1,\r2
	RA4 \a1,\f0,\r0,\r1,\r2,\x0,\x1,\x2,\x3,\x4_,\x5,xzr,\f1,\f2,\f3,\r3
	RA4 \b1,\f0,\r0,\r1,\r2,\y0,\y1,\y2,\y3,\y4,\y5,xzr,\f1,\f2,\f3,\r3
	SHX \f1,\f0,\r2,\r0,\r1
	RA4 \a2,\f1,\f0,\r0,\r1,\x0,\x1,\x2,\x3,\x4_,\x5,\x5,\f2,\f3,\r3,\r2
	RA4 \b2,\f1,\f0,\r0,\r1,\y0,\y1,\y2,\y3,\y4,\y5,\y5,\f2,\f3,\r3,\r2
	SHX \f2,\f1,\r1,\f0,\r0
	RA4 \a3,\f2,\f1,\f0,\r0,\x0,\x1,\x2,\x3,\x4_,\x5,\x5,\f3,\r3,\r2,\r1
	RA4 \b3,\f2,\f1,\f0,\r0,\y0,\y1,\y2,\y3,\y4,\y5,\y5,\f3,\r3,\r2,\r1
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
	RF x9,x25,x26,x27,x28,x13,x14,x15,x16,x17,x23,x24
	SHX x5,x25,x28,x26,x27
	RA x10,x5,x25,x26,x27,x6,x7,x8,x23,x24,x13,x14,x15,x16,x17,x29,xzr,x15,x16
	SHX x28,x5,x27,x25,x26
	RA x11,x28,x5,x25,x26,x6,x7,x8,x23,x24,x13,x14,x15,x16,x17,x29,x29,x15,x16
	SHX x27,x28,x26,x5,x25
	RA x12,x27,x28,x5,x25,x6,x7,x8,x9,x10,x13,x14,x15,x16,x17,x29,x29,x27,x28
	FDX x13,x14,x17,x29,x25,x26,x5,x6,x27,x28
	.endm
	.macro FWDI
	ldp x7,x8,[x3,#208]
	ldp x5,x6,[x3,#192]
	FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
	stp x11,x12,[x3,#16]
	stp x9,x10,[x3],#256
	stp x27,x28,[x3,#16]
	stp x25,x26,[x3],#256
	P1X
	.endm
	.macro CP mo,mp,r0,r1,r2,r3,r4,r5,r6,r7
	ldp x5,x6,[x3,#\mo]
	ldp x9,x10,[sp,#\r0]
	ldp x11,x12,[sp,#\r1]
	RF x5,x9,x10,x11,x12,x13,x14,x15,x16,x17,x9,x19
	ldp x9,x10,[sp,#\r2]
	ldp x11,x12,[sp,#\r3]
	RA4 x6,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,xzr,x19,x9,x10,x11
	ldp x5,x6,[x3,#\mp]
	ldp x9,x10,[sp,#\r4]
	ldp x11,x12,[sp,#\r5]
	RA4 x5,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x9,x10,x11
	ldp x9,x10,[sp,#\r6]
	ldp x11,x12,[sp,#\r7]
	RA4 x6,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x9,x10,x11
	.endm
	.macro PAIRB
	ldp x5,x6,[sp,#-128]
	ldp x16,x17,[sp,#-112]
	RF x21,x5,x6,x16,x17,x29,x27,x7,x1,x13,x9,x10
	RF x25,x5,x6,x16,x17,x30,x28,x8,x20,x19,x9,x10
	ldp x5,x6,[sp,#-96]
	ldp x16,x17,[sp,#-80]
	RA4 x22,x5,x6,x16,x17,x29,x27,x7,x1,x13,x14,xzr,x9,x10,x11,x15
	RA4 x26,x5,x6,x16,x17,x30,x28,x8,x20,x19,x18,xzr,x9,x10,x11,x15
	ldp x5,x6,[sp,#-64]
	ldp x16,x17,[sp,#-48]
	RA4 x23,x5,x6,x16,x17,x29,x27,x7,x1,x13,x14,x14,x9,x10,x11,x15
	RA4 x0,x5,x6,x16,x17,x30,x28,x8,x20,x19,x18,x18,x9,x10,x11,x15
	ldp x5,x6,[sp,#-32]
	ldp x16,x17,[sp,#-16]
	RA4 x24,x5,x6,x16,x17,x29,x27,x7,x1,x13,x14,x14,x9,x10,x11,x15
	RA4 x2,x5,x6,x16,x17,x30,x28,x8,x20,x19,x18,x18,x9,x10,x11,x15
	FDX x29,x27,x13,x14,x29,x27,x5,x13,x7,x1
	FDX x30,x28,x19,x18,x30,x28,x6,x19,x8,x20
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
	.macro TROW op,a,b,c,d,r3,r2,r1,t
	\op \a,x5,x6,\r3,\r2,\r1,\t
	\op \b,x9,x10,\r3,\r2,\r1,\t
	\op \c,x11,x12,\r3,\r2,\r1,\t
	\op \d,x17,x18,\r3,\r2,\r1,\t
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
	.macro TFP4
	TROW TR0,x21,x25,x29,x30,x16,x15,x14,x19
	SHX x19,x13,x16,x14,x15
	TROW TRK,x22,x26,x27,x28,x15,x14,x13,x16
	SHX x16,x19,x15,x13,x14
	TROW TRK,x23,x0,x7,x8,x14,x13,x19,x15
	SHXT x16,x14,x19,x13
	TROW TRK,x24,x2,x1,x20,x13,x19,x16,x14
	.endm
	.macro TFIN f,s
	add \f,\f,\s,lsr #61
	lsl \s,\s,#3
	.endm
	.macro TFINC f,s,cl,ch
	add \f,\f,\s,lsr #61
	adds \s,\cl,\s,lsl #3
	adc \f,\f,\ch
	.endm
	.macro TFINS f,s
	add \f,\f,\s,lsr #61
	.endm
	.macro TOUTC r,c,d,p,q,m,n
	ldp x15,x16,[\r,#\c]
	ldp x13,x14,[\r,#\d]
	TFINC x5,x6,x15,x16
	TFINC x9,x10,x15,x16
	TFINC x11,x12,x13,x14
	TFINC x17,x18,x13,x14
	stp x6,x5,[x3,#\p]
	stp x12,x11,[x3,#\q]
	stp x10,x9,[x3,#\m]
	stp x18,x17,[x3,#\n]
	.endm
	.macro BWDI
	ldp x21,x22,[x3,#256]
	ldp x23,x24,[x3]
	RF x21,x25,x26,x27,x28,x9,x10,x11,x12,x13,x5,x6
	RF x23,x25,x26,x27,x28,x15,x16,x17,x18,x19,x5,x6
	SHX x5,x25,x28,x26,x27
	RA4 x22,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
	RA4 x24,x5,x25,x26,x27,x15,x16,x17,x18,x19,x20,xzr,x6,x7,x8,x28
	ldp x21,x22,[x3,#272]
	ldp x23,x24,[x3,#16]
	SHX x6,x5,x27,x25,x26
	RA4 x21,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
	RA4 x23,x6,x5,x25,x26,x15,x16,x17,x18,x19,x20,x20,x7,x8,x28,x27
	SHX x7,x6,x26,x5,x25
	RA4 x22,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
	RA4 x24,x7,x6,x5,x25,x15,x16,x17,x18,x19,x20,x20,x8,x28,x27,x26
	FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12
	FDX x15,x16,x19,x20,x15,x16,x6,x19,x17,x18
	stp x15,x16,[x3,#16]
	stp x17,x18,[x3]
	ldp x13,x14,[x3,#224]
	ldp x15,x16,[x3,#240]
	ldp x17,x18,[sp,#400]
	ldp x19,x20,[sp,#416]
	FS x21,x22,x23,x24,x13,x14,x15,x16,x17,x18,x19,x20,x5
	FADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
	PAIR x21,x22,x23,x24,x13,x14,x15,x16,x13,x14,x15,x16,x9,x10,x11,x12,x5,x6,x7,x8,x17,x18,x19,x20,x25,x26,x27,x28,x29,x30,x0,x2
	FDX x17,x18,x25,x26,x17,x18,x5,x25,x19,x20
	FDX x27,x28,x0,x2,x27,x28,x6,x0,x29,x30
	SQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
	F5 x21,x22,x23,x24,x8,x21,x22,x9,x8
	SQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
	F5 x25,x26,x0,x2,x8,x25,x26,x9,x8
	PAIRB
	CP 352,368,16,32,48,64,80,96,112,128
	ldp x9,x10,[x3,#448]
	ldp x11,x12,[x3,#464]
	AD4 x13,x14,x15,x16,x17,x18,x9,x10,x11,x12
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	TFP4
	TOUTC x3,160,176,320,336,64,80
	CP 384,400,144,160,176,192,208,224,240,256
	ldp x9,x10,[sp,#432]
	ldp x11,x12,[sp,#448]
	AD4 x13,x14,x15,x16,x17,x18,x9,x10,x11,x12
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	TFP4
	TOUTC sp,464,480,288,304,96,112
	CP 416,432,272,288,304,320,336,352,368,384
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	TFP4
	TFIN x5,x6
	TFIN x9,x10
	TFINS x11,x12
	TFINS x17,x18
	stp x6,x5,[x3,#256]
	stp x12,x11,[x3,#272]
	stp x10,x9,[x3,#128]
	stp x18,x17,[x3,#144]
	ldp x25,x26,[x3,#16]
	ldp x27,x28,[x3],#-512
	.endm
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
	.macro LJ0 k,w0l,w0h,w1l,w1h
	add x15,\w0h,\w1h
	cmp \w0h,\w1h
	b.hs kh\k
	cmn \w0h,x15
	b.hs kt\k
	and x15,x1,\w0l,lsr #25
	ldadd x4,x16,[x15]
	lsr \w1l,\w1l,#28
	add \w1l,\w1l,x3,lsl #30
	stp \w0h,\w1h,[x16]
	stp x17,\w1l,[x16,#16]
	tbz x16,#13,kn\k
	FLUSH \w0h,\w1h
kn\k:
	.endm
	.macro LJX k,w0l,w0h,w1l,w1h
kh\k:
	cmn \w1h,x15
	b.hs kt\k
	adds \w0l,\w0l,\w1l
	adc \w0h,\w0h,\w1h
	mvn \w0h,\w0h
	and x15,x1,\w1l,lsr #25
	ldadd x4,x16,[x15]
	eor \w0l,x14,\w0l,lsr #28
	add \w0l,\w0l,x3,lsl #30
	stp \w1h,\w0h,[x16]
	stp x17,\w0l,[x16,#16]
	tbz x16,#13,km\k
	FLUSH \w0h,\w1h
km\k:
	b kn\k
kt\k:
	adds \w1l,\w0l,\w1l
	adc \w1h,\w0h,\w1h
	mvn \w1h,\w1h
	bic x15,x1,\w1l,lsr #25
	ldadd x4,x16,[x15]
	lsr \w0l,\w0l,#28
	add \w0l,\w0l,x3,lsl #30
	stp \w1h,\w0h,[x16]
	stp x17,\w0l,[x16,#16]
	tbz x16,#13,ku\k
	FLUSH \w0h,\w1h
ku\k:
	b kn\k
	.endm
	.macro TLDP
	ldp x19,x20,[x3],#-192
	ldp x21,x22,[x3,#208]
	ldp x23,x24,[x3,#160]
	ldp x25,x26,[x3,#176]
	ldp x27,x28,[x3,#128]
	ldp x29,x2,[x3,#144]
	.endm
	.macro TLDM
	ldp x27,x28,[x3],#-320
	ldp x29,x2,[x3,#336]
	ldp x19,x20,[x3,#256]
	ldp x21,x22,[x3,#272]
	ldp x23,x24,[x3,#288]
	ldp x25,x26,[x3,#304]
	.endm
	.macro TSUM k0,k1,k2
	adds x5,x19,x23
	adc x6,x20,x24
	adds x5,x5,x27
	adc x6,x6,x28
	adds x7,x21,x25
	adc x8,x22,x26
	adds x7,x7,x29,lsl #3
	adc x8,x8,x2
	LJ0 \k0,x5,x6,x7,x8
	adds x9,x23,x25
	adc x10,x24,x26
	adds x11,x27,x29,lsl #3
	adc x13,x28,x2
	adds x25,x19,x25
	adc x26,x20,x26
	subs x25,x25,x11
	sbc x26,x26,x13
	adds x27,x21,x27
	adc x28,x22,x28
	subs x27,x27,x9
	sbc x28,x28,x10
	LJ0 \k1,x25,x26,x27,x28
	adds x19,x19,x29,lsl #3
	adc x20,x20,x2
	subs x19,x19,x9
	sbc x20,x20,x10
	adds x21,x21,x23
	adc x22,x22,x24
	subs x21,x21,x11
	sbc x22,x22,x13
	LJ0 \k2,x19,x20,x21,x22
	.endm
	.macro TJX k0,k1,k2
	LJX \k0,x5,x6,x7,x8
	LJX \k1,x25,x26,x27,x28
	LJX \k2,x19,x20,x21,x22
	.endm
	.text
	.p2align 2
	.globl fwd
fwd:
	mov x4,#977
	movk x4,#1,lsl #32
	str xzr,[sp,#-16]!
	subs x0,x0,#8
	b.lo 3f
1:
	FWDI
	FWDI
	FWDI
	FWDI
	FWDI
	FWDI
	FWDI
	FWDI
	subs x0,x0,#8
	b.hs 1b
3:
	add x0,x0,#8
	cbz x0,5f
4:
	FWDI
	FWDI
	subs x0,x0,#2
	b.ne 4b
5:
	ldr x0,[sp],#16
	ret
	.globl bwd
bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	stp xzr,x30,[sp,#-16]!
1:
	BWDI
	BWDI
	tbz x3,#34,1b
	ldp x0,x30,[sp],#16
	ret
	.globl bwd8
bwd8:
	mov x4,#977
	movk x4,#1,lsl #32
	stp xzr,x30,[sp,#-16]!
1:
	BWDI
	BWDI
	BWDI
	BWDI
	BWDI
	BWDI
	BWDI
	BWDI
	tbz x3,#34,1b
	ldp x0,x30,[sp],#16
	ret
	.globl keys
keys:
	mov x1,#0x7ffffffff8
	mov x4,#32
	mov x0,#8192
	mov x14,#0xfffffffff
1:
	TLDP
	TSUM p00,p01,p02
	TLDM
	TSUM m00,m01,m02
	TLDP
	TSUM p10,p11,p12
	TLDM
	TSUM m10,m11,m12
	tbz x3,#34,1b
	ret
	TJX p00,p01,p02
	TJX m00,m01,m02
	TJX p10,p11,p12
	TJX m10,m11,m12

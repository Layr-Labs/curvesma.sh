// Phase-0 batched +/- affine X-pair mining kernels (xmine_fwd, xmine_bwd).
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
	adcs \o1,\o1,xzr
	b.lo 9f
	adcs \o2,\o2,xzr
	adcs \o3,\o3,xzr
	csel \t,x4,xzr,hs
	add \o0,\o0,\t
	stadd x4,[sp]
9:
	.endm
	.macro MUL4_M a0,a1,a2,a3,b0,b1,b2,b3,o0,o1,o2,o3,e0,e1,e4,e5,t0,t1,t2,t3,sh
	RF \a0,\b0,\b1,\b2,\b3,\e0,\e1,\o2,\o3,\e4,\t0,\t1
	SHX \sh,\b0,\b3,\b1,\b2
	RA4 \a1,\sh,\b0,\b1,\b2,\e0,\e1,\o2,\o3,\e4,\e5,xzr,\t0,\t1,\t2,\t3
	SHX \b3,\sh,\b2,\b0,\b1
	RA4 \a2,\b3,\sh,\b0,\b1,\e0,\e1,\o2,\o3,\e4,\e5,\e5,\t0,\t1,\t2,\t3
	SHX \b2,\b3,\b1,\sh,\b0
	RA4 \a3,\b2,\b3,\sh,\b0,\e0,\e1,\o2,\o3,\e4,\e5,\e5,\t0,\t1,\t2,\t3
	FDX \e0,\e1,\e4,\e5,\o0,\o1,\t0,\e4,\o2,\o3
	.endm
	.macro LAM_TO_A
	SQ13 x8,x9,x10,x11,x13,x14,x15,x16,x5,x6,x7,x12,x17
	F5 x13,x14,x15,x16,x12,x13,x14,x17,x12
	ldp x8,x9,[sp,#64]
	ldp x10,x11,[sp,#80]
	FS x20,x21,x22,x23,x13,x14,x15,x16,x8,x9,x10,x11,x5
	SQ13 x20,x21,x22,x23,x13,x14,x15,x16,x5,x6,x7,x12,x17
	F5 x13,x14,x15,x16,x12,x13,x14,x17,x12
	MUL4_M x20,x21,x22,x23,x13,x14,x15,x16,x8,x9,x10,x11,x5,x6,x7,x12,x17,x18,x19,x24,x25
	mov x13,#7
	FADDW x8,x9,x10,x11,x8,x9,x10,x11,x13,xzr,xzr,xzr,x5
	stp x8,x9,[x28],#16
	stp x10,x11,[x28],#16
	.endm
	.text
	.p2align 2
	.globl xmine_fwd
xmine_fwd:
	mov x4,#977
	movk x4,#1,lsl #32
	sub sp,sp,#32
	stp x29,x30,[sp,#16]
	str xzr,[sp]
	ldp x16,x17,[x1]
	ldp x18,x19,[x1,#16]
1:
	ldp x20,x21,[x2],#64
	ldp x22,x23,[x2,#-48]
	FS x20,x21,x22,x23,x20,x21,x22,x23,x16,x17,x18,x19,x5
	stp x20,x21,[x3]
	stp x22,x23,[x3,#16]
	MUL4_M x20,x21,x22,x23,x24,x25,x26,x27,x8,x9,x10,x11,x6,x7,x12,x13,x14,x15,x28,x29,x5
	mov x24,x8
	mov x25,x9
	mov x26,x10
	mov x27,x11
	stp x24,x25,[x3,#32]
	stp x26,x27,[x3,#48]
	add x3,x3,#64
	subs x0,x0,#1
	b.ne 1b
	ldr x0,[sp]
	ldp x29,x30,[sp,#16]
	add sp,sp,#32
	ret
	.globl xmine_bwd
xmine_bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	sub sp,sp,#208
	stp x29,x30,[sp,#192]
	str xzr,[sp]
	ldp x8,x9,[x1]
	ldp x10,x11,[x1,#16]
	stp x8,x9,[sp,#16]
	stp x10,x11,[sp,#32]
	ldp x8,x9,[x1,#32]
	ldp x10,x11,[x1,#48]
	stp x8,x9,[sp,#96]
	stp x10,x11,[sp,#112]
1:
	ldp x8,x9,[x3,#-32]
	ldp x10,x11,[x3,#-16]
	stp x20,x21,[sp,#128]
	stp x22,x23,[sp,#144]
	MUL4_M x8,x9,x10,x11,x20,x21,x22,x23,x24,x25,x26,x27,x5,x6,x7,x12,x13,x14,x15,x16,x17
	stp x24,x25,[sp,#160]
	stp x26,x27,[sp,#176]
	ldp x20,x21,[sp,#128]
	ldp x22,x23,[sp,#144]
	ldp x8,x9,[x3]
	ldp x10,x11,[x3,#16]
	MUL4_M x8,x9,x10,x11,x20,x21,x22,x23,x24,x25,x26,x27,x5,x6,x7,x12,x13,x14,x15,x16,x17
	stp x24,x25,[sp,#128]
	stp x26,x27,[sp,#144]
	ldp x8,x9,[sp,#16]
	ldp x10,x11,[sp,#32]
	ldp x12,x13,[x2]
	ldp x14,x15,[x2,#16]
	FADDW x8,x9,x10,x11,x8,x9,x10,x11,x12,x13,x14,x15,x5
	stp x8,x9,[sp,#64]
	stp x10,x11,[sp,#80]
	ldp x12,x13,[x2,#32]
	ldp x14,x15,[x2,#48]
	ldp x16,x17,[sp,#96]
	ldp x18,x19,[sp,#112]
	FS x12,x13,x14,x15,x12,x13,x14,x15,x16,x17,x18,x19,x5
	ldp x24,x25,[sp,#160]
	ldp x26,x27,[sp,#176]
	MUL4_M x12,x13,x14,x15,x24,x25,x26,x27,x8,x9,x10,x11,x5,x6,x7,x16,x17,x18,x19,x20,x21
	LAM_TO_A
	ldp x12,x13,[x2,#32]
	ldp x14,x15,[x2,#48]
	ldp x16,x17,[sp,#96]
	ldp x18,x19,[sp,#112]
	FADDW x12,x13,x14,x15,x12,x13,x14,x15,x16,x17,x18,x19,x5
	FS x12,x13,x14,x15,xzr,xzr,xzr,xzr,x12,x13,x14,x15,x5
	ldp x24,x25,[sp,#160]
	ldp x26,x27,[sp,#176]
	MUL4_M x12,x13,x14,x15,x24,x25,x26,x27,x8,x9,x10,x11,x5,x6,x7,x16,x17,x18,x19,x20,x21
	LAM_TO_A
	ldp x20,x21,[sp,#128]
	ldp x22,x23,[sp,#144]
	sub x2,x2,#64
	sub x3,x3,#64
	subs x0,x0,#1
	b.ne 1b
	ldr x0,[sp]
	ldp x29,x30,[sp,#192]
	add sp,sp,#208
	ret

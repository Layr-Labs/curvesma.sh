// Stored-leaf batch kernels (proof.md section 3), six leaves per inverted element; field macros follow package F's
// secp_loops.s (RF, RA; SHX and FDX are its SH and FD with the carry taken inline, no exits). x4 holds C0 = 2^32 + 977
// throughout. Header H (x2): x_X +0, y_X +32, beta x_X +64, beta^2 x_X +96, x_X/2 +128, beta x_X/2 +160, 3x_X^2/4 +192,
// 3beta^2 x_X^2/4 +224, staging +256..+415, G1..G5 +416, rows of beta +456, rows of beta^2 +584. Y entry (160 bytes):
// y^2, y, x, 3x^2, x + C0. Batch entry (256 bytes; the 2^25 entries fill the lower half of a 2^34-aligned block): fwd
// writes D | acc, and bwd overwrites the entry with the leaves x(phi^t(X) + Y), x(phi^t(X) - Y), t = 0, 1, 2. fwd walks
// the entries upwards for x0 of them; bwd and keys walk them downwards and end when the entry pointer leaves the block
// (bit 33 set). Calls may overwrite x0-x29 and the flags; bwd and keys save and restore x30; sp is preserved. keys
// writes each record through the column port, a 48-byte device register block at x0 (destination +0, address +8, K2 K3
// +16, w2 w3 +32), at its destination's cursor; the cursor table follows at x0 + 64 (16 bytes per destination: cursor,
// slot end), and a record whose slot is full is dropped.
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

	.macro SHX l,a,d,p,q
	mul \l,\d,x4
	umulh \d,\d,x4
	adds \a,\a,\d
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	.endm

	.macro FDX e,f,i,j,o,p,l,h,y,z
	mul \l,\i,x4
	umulh \h,\i,x4
	madd \h,\j,x4,\h
	adds \o,\e,\l
	adcs \p,\f,\h
	adcs \y,\y,xzr
	adc \z,\z,xzr
	.endm

	.macro FS d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
	subs \d0,\a0,\b0
	sbcs \d1,\a1,\b1
	sbcs \d2,\a2,\b2
	sbcs \d3,\a3,\b3
	csel \t,x4,xzr,lo
	subs \d0,\d0,\t
	sbcs \d1,\d1,xzr
	sbcs \d2,\d2,xzr
	sbc \d3,\d3,xzr
	.endm

	.macro FADDW d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
	adds \d0,\a0,\b0
	adcs \d1,\a1,\b1
	adcs \d2,\a2,\b2
	adcs \d3,\a3,\b3
	csel \t,x4,xzr,hs
	adds \d0,\d0,\t
	adcs \d1,\d1,xzr
	adcs \d2,\d2,xzr
	adc \d3,\d3,xzr
	.endm

	.macro W2S6 
	ldp x5,x6,[x3,#32]
	RF x5,x25,x26,x27,x28,x10,x11,x12,x13,x7,x8,x9
	ldp x5,x8,[x3],#-256
	RF x5,x25,x26,x27,x28,x9,x14,x15,x16,x17,x19,x20
	SHX x5,x25,x28,x26,x27
	RA x6,x5,x25,x26,x27,x19,x20,x21,x22,x23,x10,x11,x12,x13,x7,x24,xzr,x12,x13
	RA x8,x5,x25,x26,x27,x6,x19,x20,x21,x22,x9,x14,x15,x16,x17,x23,xzr,x15,x16
	SHX x28,x5,x27,x25,x26
	ldp x6,x8,[x3,#304]
	RA x6,x28,x5,x25,x26,x19,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	ldp x6,x19,[x3,#272]
	RA x6,x28,x5,x25,x26,x20,x21,x22,x29,x30,x9,x14,x15,x16,x17,x23,x23,x15,x16
	SHX x27,x28,x26,x5,x25
	RA x8,x27,x28,x5,x25,x6,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	RA x19,x27,x28,x5,x25,x6,x8,x20,x21,x22,x9,x14,x15,x16,x17,x23,x23,x27,x28
	FDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
	FDX x9,x14,x17,x23,x25,x26,x5,x6,x27,x28
	.endm

	.macro SQ5 ia,ib,ic,id,ta,tb,tc,td,te,tf,tg,th,ti,tj,tk,tl,tm,tn,to
	mul \ta,\ia,\ib
	umulh \tb,\ia,\ib
	mul \tc,\ia,\ic
	umulh \td,\ia,\ic
	mul \te,\ia,\id
	umulh \tf,\ia,\id
	mul \tg,\ib,\ic
	umulh \th,\ib,\ic
	mul \ti,\ib,\id
	umulh \tj,\ib,\id
	mul \tk,\ic,\id
	umulh \tl,\ic,\id
	adds \te,\te,\tg
	adcs \tf,\tf,\th
	adc \tj,\tj,xzr
	adds \tb,\tb,\tc
	adcs \te,\te,\td
	adcs \tf,\tf,\ti
	adcs \tj,\tj,\tk
	adc \tl,\tl,xzr
	lsr \tc,\tl,#63
	extr \tl,\tl,\tj,#63
	extr \tj,\tj,\tf,#63
	extr \tf,\tf,\te,#63
	extr \te,\te,\tb,#63
	extr \tb,\tb,\ta,#63
	mul \td,\ia,\ia
	umulh \tg,\ia,\ia
	mul \th,\ib,\ib
	umulh \ti,\ib,\ib
	mul \tk,\ic,\ic
	umulh \tm,\ic,\ic
	mul \tn,\id,\id
	umulh \to,\id,\id
	adds \ta,\tg,\ta,lsl #1
	adcs \tb,\tb,\th
	adcs \te,\te,\ti
	adcs \tf,\tf,\tk
	adcs \tj,\tj,\tm
	adcs \tl,\tl,\tn
	adc \tc,\tc,\to
	mul \tg,\tf,x4
	umulh \th,\tf,x4
	mul \ti,\tj,x4
	umulh \tk,\tj,x4
	mul \tm,\tl,x4
	umulh \tn,\tl,x4
	mul \to,\tc,x4
	umulh \tf,\tc,x4
	adds \ti,\ti,\th
	adcs \tm,\tm,\tk
	adcs \to,\to,\tn
	adc \tf,\tf,xzr
	adds \td,\td,\tg
	adcs \ta,\ta,\ti
	adcs \tb,\tb,\tm
	adcs \te,\te,\to
	adc \tf,\tf,xzr
	.endm

	.macro F5 w0,w1,w2,w3,w4,o0,o1,o2,o3,ta,tb
	mul \ta,\w4,x4
	umulh \tb,\w4,x4
	adds \o0,\w0,\ta
	adcs \o1,\w1,\tb
	adcs \o2,\w2,xzr
	adc \o3,\w3,xzr
	.endm

	.macro AC4 w0,w1,w2,w3,w4,c0,c1,c2,c3
	adds \w0,\w0,\c0
	adcs \w1,\w1,\c1
	adcs \w2,\w2,\c2
	adcs \w3,\w3,\c3
	adc \w4,\w4,xzr
	.endm

	.macro SB6 e,f,g,h,i,j,s0,s1,s2,s3
	subs \e,\e,\s0
	sbcs \f,\f,\s1
	sbcs \g,\g,\s2
	sbcs \h,\h,\s3
	adcs \i,\i,xzr
	adc \j,\j,xzr
	.endm

	.macro LT e0,e1,l0,l1,l2,l3
	SQ5 x10,x11,x12,x13,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
	F5 x8,x5,x6,x9,x14,x10,x11,x12,x13,x15,x16
	ldp x19,x23,[x2,#288]
	RF x19,x10,x11,x12,x13,x14,x15,x16,x17,x5,x6,x7
	RF x23,x10,x11,x12,x13,x8,x9,x21,x22,x0,x6,x7
	SHX x6,x10,x13,x11,x12
	ldp x20,x24,[x2,#304]
	RA x20,x6,x10,x11,x12,x7,x13,x18,x19,x23,x14,x15,x16,x17,x5,x29,xzr,x16,x17
	RA x24,x6,x10,x11,x12,x7,x13,x18,x19,x20,x8,x9,x21,x22,x0,x23,xzr,x21,x22
	SHX x7,x6,x12,x10,x11
	ldp x12,x13,[x2,#320]
	RA x12,x7,x6,x10,x11,x18,x19,x20,x24,x30,x14,x15,x16,x17,x5,x29,x29,x16,x17
	RA x13,x7,x6,x10,x11,x18,x19,x20,x24,x30,x8,x9,x21,x22,x0,x23,x23,x21,x22
	SHX x12,x7,x11,x6,x10
	ldp x11,x13,[x2,#336]
	RA x11,x12,x7,x6,x10,x18,x19,x20,x24,x30,x14,x15,x16,x17,x5,x29,x29,x16,x17
	RA x13,x12,x7,x6,x10,x18,x19,x20,x24,x30,x8,x9,x21,x22,x0,x23,x23,x12,x13
	ldp x6,x7,[x2,#\e0]
	ldp x10,x11,[x2,#\e1]
	ldp x18,x19,[x1,#96]
	ldp x20,x21,[x1,#112]
	FADDW x6,x7,x10,x11,x6,x7,x10,x11,x18,x19,x20,x21,x22
	SB6 x14,x15,x16,x17,x5,x29,x6,x7,x10,x11
	SB6 x8,x9,x12,x13,x0,x23,x6,x7,x10,x11
	FDX x14,x15,x5,x29,x14,x15,x18,x19,x16,x17
	FDX x8,x9,x0,x23,x8,x9,x18,x19,x12,x13
	stp x14,x15,[x3,#\l0]
	stp x16,x17,[x3,#\l1]
	stp x8,x9,[x3,#\l2]
	stp x12,x13,[x3,#\l3]
	.endm

	.macro W2L6 
	RF x19,x10,x11,x12,x13,x14,x15,x16,x17,x5,x6,x7
	RF x23,x10,x11,x12,x13,x8,x9,x21,x22,x0,x6,x7
	SHX x6,x10,x13,x11,x12
	RA x20,x6,x10,x11,x12,x7,x13,x18,x19,x23,x14,x15,x16,x17,x5,x29,xzr,x16,x17
	RA x24,x6,x10,x11,x12,x7,x13,x18,x19,x20,x8,x9,x21,x22,x0,x23,xzr,x21,x22
	SHX x7,x6,x12,x10,x11
	ldp x12,x13,[x2,#256]
	RA x12,x7,x6,x10,x11,x18,x19,x20,x24,x30,x14,x15,x16,x17,x5,x29,x29,x16,x17
	RA x13,x7,x6,x10,x11,x18,x19,x20,x24,x30,x8,x9,x21,x22,x0,x23,x23,x21,x22
	SHX x12,x7,x11,x6,x10
	ldp x11,x13,[x2,#272]
	RA x11,x12,x7,x6,x10,x18,x19,x20,x24,x30,x14,x15,x16,x17,x5,x29,x29,x16,x17
	RA x13,x12,x7,x6,x10,x18,x19,x20,x24,x30,x8,x9,x21,x22,x0,x23,x23,x12,x13
	FDX x14,x15,x5,x29,x14,x15,x6,x7,x16,x17
	FDX x8,x9,x0,x23,x10,x11,x6,x7,x12,x13
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

	.macro KEY 
	umulh x24,x20,x10
	umulh x25,x21,x9
	umulh x26,x22,x8
	umulh x27,x23,x7
	add x24,x24,x25
	add x26,x26,x27
	add x24,x24,x26
	madd x24,x20,x11,x24
	madd x24,x21,x10,x24
	madd x24,x22,x9,x24
	madd x24,x23,x8,x24
	add x25,x23,x24
	subs x26,x23,x24
	csel x27,x24,x23,hs
	csel x29,x14,xzr,hs
	adds x28,x27,x25
	csel x29,x15,x29,hs
	cmn x24,#16
	ccmn x25,#17,#2,lo
	ccmn x28,x16,#2,lo
	ccmp x26,#16,#0,lo
	.endm

	.macro KMUL 
	ldp x1,x2,[x29,#0]
	mul x24,x20,x1
	umulh x5,x20,x1
	mul x6,x20,x2
	adds x25,x5,x6
	umulh x5,x20,x2
	ldp x1,x2,[x29,#16]
	mul x6,x20,x1
	adcs x26,x5,x6
	umulh x5,x20,x1
	mul x6,x20,x2
	adcs x27,x5,x6
	umulh x5,x20,x2
	adc x28,x5,xzr
	ldp x1,x2,[x29,#32]
	mul x5,x21,x1
	umulh x6,x21,x1
	mul x20,x21,x2
	adds x6,x6,x20
	umulh x13,x21,x2
	ldp x1,x2,[x29,#48]
	mul x20,x21,x1
	adcs x13,x13,x20
	umulh x19,x21,x1
	mul x20,x21,x2
	adcs x19,x19,x20
	umulh x20,x21,x2
	adc x20,x20,xzr
	adds x24,x24,x5
	adcs x25,x25,x6
	adcs x26,x26,x13
	adcs x27,x27,x19
	adcs x28,x28,x20
	adc x30,xzr,xzr
	ldp x1,x2,[x29,#64]
	mul x5,x22,x1
	umulh x6,x22,x1
	mul x20,x22,x2
	adds x6,x6,x20
	umulh x13,x22,x2
	ldp x1,x2,[x29,#80]
	mul x20,x22,x1
	adcs x13,x13,x20
	umulh x19,x22,x1
	mul x20,x22,x2
	adcs x19,x19,x20
	umulh x20,x22,x2
	adc x20,x20,xzr
	adds x24,x24,x5
	adcs x25,x25,x6
	adcs x26,x26,x13
	adcs x27,x27,x19
	adcs x28,x28,x20
	adc x30,x30,xzr
	ldp x1,x2,[x29,#96]
	mul x5,x23,x1
	umulh x6,x23,x1
	mul x20,x23,x2
	adds x6,x6,x20
	umulh x13,x23,x2
	ldp x1,x2,[x29,#112]
	mul x20,x23,x1
	adcs x13,x13,x20
	umulh x19,x23,x1
	mul x20,x23,x2
	adcs x19,x19,x20
	umulh x20,x23,x2
	adc x20,x20,xzr
	adds x24,x24,x5
	adcs x25,x25,x6
	adcs x26,x26,x13
	adcs x27,x27,x19
	adcs x28,x28,x20
	adc x30,x30,xzr
	FDX x24,x25,x28,x30,x24,x25,x5,x6,x26,x27
	.endm

	.macro LEAF6 
	ldp x22,x23,[x3,#16]
	ldp x20,x21,[x3],#32
	KEY
	b.lo 3f
	cbz x29,2f
	KMUL
	add x18,x18,#1
	and x5,x24,#0xfffffffff
	add x19,x0,x5,lsl #4
	ldp x13,x6,[x19,#64]
	cmp x13,x6
	b.hs 3f
	and x6,x25,#0xfffffffff
	add x6,x6,x3,lsl #31
	stp x5,x13,[x0]
	stp x26,x27,[x0,#16]
	stp x17,x6,[x0,#32]
	add x13,x13,#32
	str x13,[x19,#64]
	b 3f
2:
	and x5,x20,#0xfffffffff
	add x19,x0,x5,lsl #4
	ldp x13,x6,[x19,#64]
	cmp x13,x6
	b.hs 3f
	and x6,x21,#0xfffffffff
	add x6,x6,x3,lsl #31
	stp x5,x13,[x0]
	stp x22,x23,[x0,#16]
	stp x17,x6,[x0,#32]
	add x13,x13,#32
	str x13,[x19,#64]
3:
	.endm

	.text
	.p2align 2
	.globl fwd
fwd:
	mov x4,#977
	movk x4,#1,lsl #32
1:
	ldp x7,x8,[x1,#16]
	ldp x5,x6,[x1],#160
	FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
	stp x11,x12,[x3,#16]
	stp x25,x26,[x3,#32]
	stp x27,x28,[x3,#48]
	stp x9,x10,[x3],#256
	P1X
	subs x0,x0,#1
	b.ne 1b
	ret
	.globl bwd
bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	str x30,[sp,#-16]!
1:
	W2S6
	ldp x5,x6,[x1,#-160]!
	ldp x7,x8,[x1,#16]
	ldp x14,x15,[x2,#32]
	ldp x16,x17,[x2,#48]
	FS x19,x20,x21,x22,x5,x6,x7,x8,x14,x15,x16,x17,x9
	FADDW x23,x24,x0,x18,x5,x6,x7,x8,x14,x15,x16,x17,x9
	stp x21,x0,[x2,#256]
	stp x22,x18,[x2,#272]
	W2L6
	SQ5 x14,x15,x16,x17,x5,x6,x7,x8,x9,x0,x18,x29,x30,x19,x20,x21,x22,x23,x24
	F5 x8,x5,x6,x9,x0,x19,x20,x21,x22,x18,x29
	SQ5 x10,x11,x12,x13,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
	F5 x8,x5,x6,x9,x14,x10,x11,x12,x13,x15,x16
	stp x19,x10,[x2,#288]
	stp x20,x11,[x2,#304]
	stp x21,x12,[x2,#320]
	stp x22,x13,[x2,#336]
	ldp x14,x15,[x1,#32]
	ldp x16,x17,[x1,#48]
	ldp x10,x11,[x2,#128]
	ldp x12,x13,[x2,#144]
	ldp x19,x20,[x2,#160]
	ldp x21,x22,[x2,#176]
	FADDW x10,x11,x12,x13,x14,x15,x16,x17,x10,x11,x12,x13,x9
	FADDW x19,x20,x21,x22,x14,x15,x16,x17,x19,x20,x21,x22,x9
	SQ5 x10,x11,x12,x13,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
	ldp x15,x16,[x2,#192]
	ldp x17,x18,[x2,#208]
	AC4 x8,x5,x6,x9,x14,x15,x16,x17,x18
	F5 x8,x5,x6,x9,x14,x10,x11,x12,x13,x15,x16
	SQ5 x19,x20,x21,x22,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
	ldp x15,x16,[x2,#224]
	ldp x17,x18,[x2,#240]
	AC4 x8,x5,x6,x9,x14,x15,x16,x17,x18
	F5 x8,x5,x6,x9,x14,x19,x20,x21,x22,x15,x16
	ldp x5,x6,[x1,#64]
	ldp x7,x8,[x1,#80]
	FS x5,x6,x7,x8,x5,x6,x7,x8,x10,x11,x12,x13,x9
	FS x14,x15,x16,x17,x5,x6,x7,x8,x19,x20,x21,x22,x9
	stp x19,x20,[x2,#352]
	stp x21,x22,[x2,#368]
	stp x14,x15,[x2,#384]
	stp x16,x17,[x2,#400]
	LT 0,16,256,272,288,304
	ldp x10,x11,[x2,#352]
	ldp x12,x13,[x2,#368]
	LT 64,80,320,336,352,368
	ldp x10,x11,[x2,#384]
	ldp x12,x13,[x2,#400]
	LT 96,112,384,400,416,432
	tbz x3,#33,1b
	ldr x30,[sp],#16
	ret
	.globl keys
keys:
	mov x4,#977
	movk x4,#1,lsl #32
	ldp x7,x8,[x2,#416]
	ldp x9,x10,[x2,#432]
	ldr x11,[x2,#448]
	add x14,x2,#456
	add x15,x2,#584
	mov x16,#33
	str x30,[sp,#-16]!
1:
	LEAF6
	LEAF6
	LEAF6
	LEAF6
	LEAF6
	LEAF6
	sub x3,x3,#448
	tbz x3,#33,1b
	ldr x30,[sp],#16
	ret

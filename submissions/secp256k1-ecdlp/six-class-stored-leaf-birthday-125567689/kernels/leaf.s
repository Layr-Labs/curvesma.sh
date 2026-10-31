// Stored-leaf batch kernels (proof.md section 3), emitted by the package author's generator; field macros follow
// package F's secp_loops.s (RF, RA; SHX and FDX are its SH and FD with the carry taken inline, no exits).
// x4 holds C0 = 2^32 + 977 throughout. Header H (x2): x_X +0, y_X +32, R spill +64, s- +96, s+ +128,
// G1..G5 +160, rows of beta +208, rows of beta^2 +336. Y entry (96 bytes): x, y, x + 2C0. Slot (64 bytes): dx | acc.
// Calls may overwrite x0-x29 and the flags; bwd and keys save and restore x30; sp is preserved. keys writes each
// record through the column port, a 48-byte device register block at x0 (destination +0, address +8, K2 K3 +16,
// w2 w3 +32), at its destination's cursor; the cursor table follows at x0 + 64 (16 bytes per destination: cursor,
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

	.macro W2S 
	ldp x5,x6,[x3,#32]
	RF x5,x25,x26,x27,x28,x10,x11,x12,x13,x7,x8,x9
	ldp x5,x8,[x3],#-64
	RF x5,x25,x26,x27,x28,x9,x14,x15,x16,x17,x19,x20
	SHX x5,x25,x28,x26,x27
	RA x6,x5,x25,x26,x27,x19,x20,x21,x22,x23,x10,x11,x12,x13,x7,x24,xzr,x12,x13
	RA x8,x5,x25,x26,x27,x6,x19,x20,x21,x22,x9,x14,x15,x16,x17,x23,xzr,x15,x16
	SHX x28,x5,x27,x25,x26
	ldp x6,x8,[x3,#112]
	RA x6,x28,x5,x25,x26,x19,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	ldp x6,x19,[x3,#80]
	RA x6,x28,x5,x25,x26,x20,x21,x22,x29,x30,x9,x14,x15,x16,x17,x23,x23,x15,x16
	SHX x27,x28,x26,x5,x25
	RA x8,x27,x28,x5,x25,x6,x20,x21,x22,x29,x10,x11,x12,x13,x7,x24,x24,x12,x13
	RA x19,x27,x28,x5,x25,x6,x8,x20,x21,x22,x9,x14,x15,x16,x17,x23,x23,x27,x28
	FDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
	FDX x9,x14,x17,x23,x25,x26,x5,x6,x27,x28
	.endm

	.macro W2L 
	ldp x5,x6,[x2,#128]
	RF x5,x10,x11,x12,x13,x25,x26,x27,x28,x7,x8,x9
	ldp x5,x8,[x2,#96]
	RF x5,x10,x11,x12,x13,x9,x14,x15,x16,x17,x19,x20
	SHX x5,x10,x13,x11,x12
	RA x6,x5,x10,x11,x12,x19,x20,x21,x22,x23,x25,x26,x27,x28,x7,x24,xzr,x27,x28
	RA x8,x5,x10,x11,x12,x6,x19,x20,x21,x22,x9,x14,x15,x16,x17,x23,xzr,x15,x16
	SHX x13,x5,x12,x10,x11
	ldp x6,x8,[x2,#144]
	RA x6,x13,x5,x10,x11,x19,x20,x21,x22,x29,x25,x26,x27,x28,x7,x24,x24,x27,x28
	ldp x6,x19,[x2,#112]
	RA x6,x13,x5,x10,x11,x20,x21,x22,x29,x30,x9,x14,x15,x16,x17,x23,x23,x15,x16
	SHX x12,x13,x11,x5,x10
	RA x8,x12,x13,x5,x10,x6,x20,x21,x22,x29,x25,x26,x27,x28,x7,x24,x24,x27,x28
	RA x19,x12,x13,x5,x10,x6,x8,x20,x21,x22,x9,x14,x15,x16,x17,x23,x23,x12,x13
	FDX x25,x26,x7,x24,x25,x26,x5,x6,x27,x28
	FDX x9,x14,x17,x23,x10,x11,x5,x6,x12,x13
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

	.macro SQP 
	mul x5,x25,x26
	umulh x6,x25,x26
	mul x7,x25,x27
	umulh x8,x25,x27
	mul x9,x25,x28
	umulh x14,x25,x28
	mul x15,x26,x27
	umulh x16,x26,x27
	mul x17,x26,x28
	umulh x19,x26,x28
	mul x20,x27,x28
	umulh x21,x27,x28
	adds x9,x9,x15
	adcs x14,x14,x16
	adc x19,x19,xzr
	adds x6,x6,x7
	adcs x9,x9,x8
	adcs x14,x14,x17
	adcs x19,x19,x20
	adc x21,x21,xzr
	lsr x7,x21,#63
	extr x21,x21,x19,#63
	extr x19,x19,x14,#63
	extr x14,x14,x9,#63
	extr x9,x9,x6,#63
	extr x6,x6,x5,#63
	mul x8,x25,x25
	umulh x15,x25,x25
	mul x16,x26,x26
	umulh x17,x26,x26
	mul x20,x27,x27
	umulh x22,x27,x27
	mul x23,x28,x28
	umulh x24,x28,x28
	adds x5,x15,x5,lsl #1
	adcs x6,x6,x16
	adcs x9,x9,x17
	adcs x14,x14,x20
	adcs x19,x19,x22
	adcs x21,x21,x23
	adc x7,x7,x24
	mul x15,x14,x4
	umulh x16,x14,x4
	mul x17,x19,x4
	umulh x20,x19,x4
	mul x22,x21,x4
	umulh x23,x21,x4
	mul x24,x7,x4
	umulh x14,x7,x4
	adds x17,x17,x16
	adcs x22,x22,x20
	adcs x24,x24,x23
	adc x14,x14,xzr
	adds x8,x8,x15
	adcs x5,x5,x17
	adcs x6,x6,x22
	adcs x9,x9,x24
	adc x14,x14,xzr
	ldp x15,x16,[x2,#0]
	ldp x17,x19,[x2,#16]
	subs x8,x8,x15
	sbcs x5,x5,x16
	sbcs x6,x6,x17
	sbcs x9,x9,x19
	adc x14,x14,xzr
	ldp x15,x16,[x1,#64]
	ldp x17,x19,[x1,#80]
	subs x8,x8,x15
	sbcs x5,x5,x16
	sbcs x6,x6,x17
	sbcs x9,x9,x19
	adc x14,x14,xzr
	mul x15,x14,x4
	umulh x16,x14,x4
	adds x19,x8,x15
	adcs x20,x5,x16
	adcs x21,x6,xzr
	adc x22,x9,xzr
	.endm

	.macro SQM 
	mul x5,x10,x11
	umulh x6,x10,x11
	mul x7,x10,x12
	umulh x8,x10,x12
	mul x9,x10,x13
	umulh x14,x10,x13
	mul x15,x11,x12
	umulh x16,x11,x12
	mul x17,x11,x13
	umulh x23,x11,x13
	mul x24,x12,x13
	umulh x25,x12,x13
	adds x9,x9,x15
	adcs x14,x14,x16
	adc x23,x23,xzr
	adds x6,x6,x7
	adcs x9,x9,x8
	adcs x14,x14,x17
	adcs x23,x23,x24
	adc x25,x25,xzr
	lsr x7,x25,#63
	extr x25,x25,x23,#63
	extr x23,x23,x14,#63
	extr x14,x14,x9,#63
	extr x9,x9,x6,#63
	extr x6,x6,x5,#63
	mul x8,x10,x10
	umulh x15,x10,x10
	mul x16,x11,x11
	umulh x17,x11,x11
	mul x24,x12,x12
	umulh x26,x12,x12
	mul x27,x13,x13
	umulh x28,x13,x13
	adds x5,x15,x5,lsl #1
	adcs x6,x6,x16
	adcs x9,x9,x17
	adcs x14,x14,x24
	adcs x23,x23,x26
	adcs x25,x25,x27
	adc x7,x7,x28
	mul x15,x14,x4
	umulh x16,x14,x4
	mul x17,x23,x4
	umulh x24,x23,x4
	mul x26,x25,x4
	umulh x27,x25,x4
	mul x28,x7,x4
	umulh x14,x7,x4
	adds x17,x17,x16
	adcs x26,x26,x24
	adcs x28,x28,x27
	adc x14,x14,xzr
	adds x8,x8,x15
	adcs x5,x5,x17
	adcs x6,x6,x26
	adcs x9,x9,x28
	adc x14,x14,xzr
	ldp x15,x16,[x2,#0]
	ldp x17,x23,[x2,#16]
	subs x8,x8,x15
	sbcs x5,x5,x16
	sbcs x6,x6,x17
	sbcs x9,x9,x23
	adc x14,x14,xzr
	ldp x15,x16,[x1,#64]
	ldp x17,x23,[x1,#80]
	subs x8,x8,x15
	sbcs x5,x5,x16
	sbcs x6,x6,x17
	sbcs x9,x9,x23
	adc x14,x14,xzr
	mul x15,x14,x4
	umulh x16,x14,x4
	adds x19,x8,x15
	adcs x20,x5,x16
	adcs x21,x6,xzr
	adc x22,x9,xzr
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

	.macro LEAFP 
	ldp x20,x21,[x3,#0]
	ldp x22,x23,[x3,#16]
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
	add x6,x6,x3,lsl #30
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
	add x6,x6,x3,lsl #30
	stp x5,x13,[x0]
	stp x22,x23,[x0,#16]
	stp x17,x6,[x0,#32]
	add x13,x13,#32
	str x13,[x19,#64]
3:
	.endm

	.macro LEAFM 
	ldp x20,x21,[x3,#32]
	ldp x22,x23,[x3,#48]
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
	add x6,x6,x3,lsl #30
	orr x6,x6,#0x8000000000000000
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
	add x6,x6,x3,lsl #30
	orr x6,x6,#0x8000000000000000
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
	ldp x5,x6,[x1],#96
	FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
	stp x11,x12,[x3,#16]
	stp x25,x26,[x3,#32]
	stp x27,x28,[x3,#48]
	stp x9,x10,[x3],#64
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
	W2S
	stp x25,x26,[x2,#64]
	stp x27,x28,[x2,#80]
	ldp x5,x6,[x1,#32]
	ldp x7,x8,[x1,#48]
	ldp x14,x15,[x2,#32]
	ldp x16,x17,[x2,#48]
	FS x19,x20,x21,x22,x5,x6,x7,x8,x14,x15,x16,x17,x9
	FADDW x23,x24,x25,x26,x5,x6,x7,x8,x14,x15,x16,x17,x9
	stp x19,x20,[x2,#128]
	stp x21,x22,[x2,#144]
	stp x23,x24,[x2,#96]
	stp x25,x26,[x2,#112]
	W2L
	SQP
	stp x19,x20,[x3,#64]
	stp x21,x22,[x3,#80]
	SQM
	stp x19,x20,[x3,#96]
	stp x21,x22,[x3,#112]
	ldp x25,x26,[x2,#64]
	ldp x27,x28,[x2,#80]
	sub x1,x1,#96
	subs x0,x0,#1
	b.ne 1b
	ldr x30,[sp],#16
	ret
	.globl keys
keys:
	mov x4,#977
	movk x4,#1,lsl #32
	ldp x7,x8,[x2,#160]
	ldp x9,x10,[x2,#176]
	ldr x11,[x2,#192]
	add x14,x2,#208
	add x15,x2,#336
	mov x16,#33
	str x30,[sp,#-16]!
1:
	LEAFP
	LEAFM
	add x3,x3,#64
	subs x12,x12,#1
	b.ne 1b
	ldr x30,[sp],#16
	ret

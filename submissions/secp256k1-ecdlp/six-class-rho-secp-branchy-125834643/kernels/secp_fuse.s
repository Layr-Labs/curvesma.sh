// Loop A: each walk's arrival, then dx; both parities BFI y0. See proof.md 2.5 and 3.2. PASS f emits pass1f,
// PASS p pass1p (no K_cp load or filter; B.HS for B.NE; a pair's addresses in one STP; x0 counts the tail).
// TBZ x2,#26 ends a group; test stubs .Ls follow RET. A dx borrow adds 1 to x17, returned in x0.
	.macro KLf
	ldr x13,[x0],#8
	.endm
	.macro KLp
	.endm
	.macro FPf
	ccmp x29,x13,#4,hs
	.endm
	.macro FPp
	.endm
	.macro BRf l
	b.ne \l
	.endm
	.macro BRp l
	b.hs \l
	.endm
	.macro Ef
	cmp x29,#0
	orr x29,x28,#1
	.endm
	.macro Ep
	cmp x29,#0
	.endm
	.macro Of
	cmn x29,#2
	sub x29,x12,x28
	.endm
	.macro Op
	cmn x29,#2
	.endm
	.macro ST0
	.endm
	.macro ST1
	str x24,[x2],#8
	.endm
	.macro ST2
	stp x12,x24,[x2],#16
	.endm
	.macro CNf o,l
	ldr x13,[sp,#\o]
	cbz x13,\l
	.endm
	.macro CNp o,l
	ldr x0,[sp,#\o]
	cbz x0,\l
	.endm
	.macro DCf o,l
	ldr x13,[sp,#\o]
	subs x13,x13,#1
	str x13,[sp,#\o]
	b.ne \l
	.endm
	.macro DCp o,l
	subs x0,x0,#1
	b.ne \l
	.endm
	.macro Pf a,b
	BODY \a,f,x24,1
	BODY \b,f,x24,1
	.endm
	.macro Pp a,b
	BODY \a,p,x12,0
	BODY \b,p,x24,2
	.endm
	.macro B8 v,k,a,b,c,d,e,f,g,h
	P\v \k\a,\k\b
	P\v \k\c,\k\d
	P\v \k\e,\k\f
	P\v \k\g,\k\h
	.endm
	.macro D8 v,k
	B8 \v,\k,0,1,2,3,4,5,6,7
	.endm
	.macro RCHG
	ldr x13,[sp]
	sub x13,x13,#20
	str x13,[sp]
	.endm
	.macro CHAIN v,c
	ccmn x25,#16,#2,\c
	ccmn x26,#17,#2,lo
	ccmn x30,x18,#2,lo
	ccmp x27,#16,#0,lo
	FP\v
	.endm
	.macro STUB
	ldr x13,[sp]
	sub x13,x13,#25
	str x13,[sp]
	ldr x13,[sp,#120]
	stp x13,x1,[x13]
	stp x28,x29,[x13,#16]
	add x13,x13,#32
	str xzr,[x13]
	str x13,[sp,#120]
	.endm
	.macro ARR v,k,r
	KL\v
	ldp x22,x23,[x1,#16]
	ldp x28,x29,[x1,#32]
	ldp x20,x21,[x1],#64
	umulh x25,x20,x10
	umulh x26,x21,x9
	add x25,x25,x26
	umulh x26,x22,x8
	add x25,x25,x26
	umulh x26,x23,x7
	add x25,x25,x26
	madd x25,x20,x11,x25
	madd x25,x21,x10,x25
	madd x25,x22,x9,x25
	madd x25,x23,x8,x25
	add x26,x23,x25
	subs x27,x23,x25
	csel x30,x25,x23,hs
	csel \r,x15,x14,hs
	adds x30,x30,x26
	csel \r,x16,\r,hs
	tbnz x28,#0,.Lo\v\k
	E\v
	CHAIN \v,ne
	BR\v .Lb\v\k
	b .Ls\v\k
.Lo\v\k:
	O\v
	CHAIN \v,lo
	BR\v .Lb\v\k
	b .Ls\v\k
.Lb\v\k:
	bfi \r,x28,#9,#27
	.endm
	.macro D r,a,b,s
	ldp x26,x27,[\r,#\a]
	ldp x29,x30,[\r,#\b]
	subs x20,x26,x20
	sbcs x21,x27,x21
	sbcs x22,x29,x22
	sbcs x23,x30,x23
	b.hs .Lb\s
	add x17,x17,#1
	subs x20,x20,x4
	b.lo .Lr\s
.Lb\s:
	.endm
	.macro BODY k,v,r,s
	ARR \v,\k,\r
	ST\s
	D \r,0,16,1_0_\v\k
	stp x22,x23,[x3,#16]
	stp x20,x21,[x3],#64
	.endm
	.macro SS s
.Lr\s:
	RCHG
	sbcs x21,x21,xzr
	sbcs x22,x22,xzr
	sbc x23,x23,xzr
	b .Lb\s
	.endm
	.macro S1 k,v
	SS 1_0_\v\k
.Ls\v\k:
	STUB
	b .Lb\v\k
	.endm
	.macro S8 v,k,a,b,c,d,e,f,g,h
	S1 \k\a,\v
	S1 \k\b,\v
	S1 \k\c,\v
	S1 \k\d,\v
	S1 \k\e,\v
	S1 \k\f,\v
	S1 \k\g,\v
	S1 \k\h,\v
	.endm
	.macro E8 v,k
	S8 \v,\k,0,1,2,3,4,5,6,7
	.endm
	.macro PASS v
	.globl pass1\v
	.type pass1\v,%function
// x0 n, x1 walks, x2 addresses, x3 slots (dx), x4 slot base, x5 K_cp, x6 G1..G5 p0+1 2^64-1, x7 rare list
pass1\v:
	stp x19,x20,[sp,#-192]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	str x7,[sp,#120]
	str xzr,[x7]
	and x13,x0,#127
	str x13,[sp,#96]
	mov x0,x5
	mov x14,x4
	add x15,x4,#64
	add x16,x4,#128
	ldp x7,x8,[x6]
	ldp x9,x10,[x6,#16]
	ldr x11,[x6,#32]
	ldp x12,x13,[x6,#40]
	mov x18,#33
	mov x17,xzr
	mov x4,#977
	movk x4,#1,lsl #32
	CN\v 96,.Lg\v
2:
	BODY 128,\v,x24,1
	DC\v 96,2b
.Lg\v:
	tbnz x2,#26,.Ldone\v
.Lg1\v:
	D8 \v,0
	D8 \v,1
	D8 \v,2
	D8 \v,3
	D8 \v,4
	D8 \v,5
	D8 \v,6
	D8 \v,7
	D8 \v,8
	D8 \v,9
	D8 \v,10
	D8 \v,11
	D8 \v,12
	D8 \v,13
	D8 \v,14
	D8 \v,15
	tbz x2,#26,.Lg1\v
.Ldone\v:
	mov x0,x17
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#192
	ret
pass1_rare\v:
	E8 \v,0
	E8 \v,1
	E8 \v,2
	E8 \v,3
	E8 \v,4
	E8 \v,5
	E8 \v,6
	E8 \v,7
	E8 \v,8
	E8 \v,9
	E8 \v,10
	E8 \v,11
	E8 \v,12
	E8 \v,13
	E8 \v,14
	E8 \v,15
	S1 128,\v
	.endm
	.text
	.p2align 2
	PASS f
	PASS p

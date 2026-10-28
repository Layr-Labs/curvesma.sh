// Pass 1 with each walk's arrival fused in; both parities BFI y0. See proof.md 2.5 and 3.2. PASS f emits pass1f,
// PASS p pass1p (no K_cp load or filter; B.HS for B.NE; a pair's addresses in one STP; x0 counts loops).
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
	STUB
	b .Lb\v\k
.Lo\v\k:
	O\v
	CHAIN \v,lo
	BR\v .Lb\v\k
	STUB
.Lb\v\k:
	bfi \r,x28,#9,#26
	.endm
	.macro BODY k,v,r,s
	ARR \v,\k,\r
	ST\s
	ldp x26,x27,[\r]
	ldp x29,x30,[\r,#16]
	subs x20,x26,x20
	sbcs x21,x27,x21
	sbcs x22,x29,x22
	sbcs x23,x30,x23
	csel x26,x4,xzr,lo
	subs x20,x20,x26
	b.lo .Lr1_0_\v\k
.Lb1_0_\v\k:
	stp x5,x6,[x3,#32]
	stp x17,x19,[x3,#48]
	stp x22,x23,[x3,#16]
	stp x20,x21,[x3],#64
	mul x13,x20,x5
	umulh x27,x20,x5
	mul x26,x20,x6
	adds x27,x27,x26
	umulh x26,x20,x6
	mul x25,x20,x17
	adcs x26,x26,x25
	umulh x25,x20,x17
	mul x24,x20,x19
	adcs x25,x25,x24
	umulh x24,x20,x19
	mul x20,x21,x19
	adcs x24,x24,x20
	umulh x20,x21,x19
	adc x20,x20,xzr
	mul x28,x21,x5
	adds x27,x27,x28
	umulh x28,x21,x5
	adcs x26,x26,x28
	umulh x28,x22,x5
	adcs x25,x25,x28
	umulh x28,x23,x5
	adcs x24,x24,x28
	umulh x28,x22,x17
	adc x28,x28,xzr
	mul x29,x22,x5
	adds x26,x26,x29
	mul x29,x23,x5
	adcs x25,x25,x29
	mul x29,x22,x17
	adcs x24,x24,x29
	adcs x20,x20,x28
	umulh x28,x22,x19
	adc x28,x28,xzr
	mul x29,x21,x6
	adds x26,x26,x29
	umulh x29,x21,x6
	adcs x25,x25,x29
	umulh x29,x21,x17
	adcs x24,x24,x29
	umulh x29,x23,x6
	adcs x20,x20,x29
	umulh x29,x23,x17
	adc x29,x29,xzr
	mul x30,x21,x17
	adds x25,x25,x30
	umulh x30,x22,x6
	adcs x24,x24,x30
	mul x30,x22,x19
	adcs x20,x20,x30
	adcs x28,x28,x29
	umulh x29,x23,x19
	adc x29,x29,xzr
	mul x30,x22,x6
	adds x25,x25,x30
	mul x30,x23,x6
	adcs x24,x24,x30
	mul x30,x23,x17
	adcs x20,x20,x30
	mul x30,x23,x19
	adcs x28,x28,x30
	adc x29,x29,xzr
	mul x30,x24,x4
	umulh x23,x24,x4
	mul x24,x20,x4
	adds x24,x24,x23
	umulh x23,x20,x4
	mul x20,x28,x4
	adcs x20,x20,x23
	umulh x23,x28,x4
	mul x28,x29,x4
	adcs x28,x28,x23
	umulh x23,x29,x4
	adc x23,x23,xzr
	adds x13,x13,x30
	adcs x27,x27,x24
	adcs x17,x26,x20
	adcs x19,x25,x28
	adc x23,x23,xzr
	mul x28,x23,x4
	umulh x20,x23,x4
	adds x5,x13,x28
	adcs x6,x27,x20
	b.cs .Lr1_1_\v\k
.Lb1_1_\v\k:
	.endm
	.macro BODYW0 v
	ARR \v,w,x24
	str x24,[x2],#8
	ldp x26,x27,[x24]
	ldp x29,x30,[x24,#16]
	subs x5,x26,x20
	sbcs x6,x27,x21
	sbcs x17,x29,x22
	sbcs x19,x30,x23
	csel x26,x4,xzr,lo
	subs x5,x5,x26
	b.lo .Lr1_w\v
.Lb1_w\v:
	.endm
	.macro S1 k,v
.Lr1_0_\v\k:
	RCHG
	sbcs x21,x21,xzr
	sbcs x22,x22,xzr
	sbc x23,x23,xzr
	b .Lb1_0_\v\k
.Lr1_1_\v\k:
	RCHG
	adcs x17,x17,xzr
	adc x19,x19,xzr
	b .Lb1_1_\v\k
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
// x0 n, x1 walks, x2 addresses, x3 products, x4 slot base, x5 K_cp, x6 G1..G5 p0+1 2^64-1, x7 rare list
pass1\v:
	stp x19,x20,[sp,#-192]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	str x7,[sp,#120]
	str xzr,[x7]
	subs x13,x0,#1
	lsr x13,x13,#9
	str x13,[sp,#104]
	subs x13,x0,#1
	and x13,x13,#511
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
	mov x4,#977
	movk x4,#1,lsl #32
	BODYW0 \v
	add x3,x3,#64
	CN\v 96,.Lg\v
2:
	BODY 999,\v,x24,1
	DC\v 96,2b
.Lg\v:
	CN\v 104,.Ldone\v
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
	D8 \v,16
	D8 \v,17
	D8 \v,18
	D8 \v,19
	D8 \v,20
	D8 \v,21
	D8 \v,22
	D8 \v,23
	D8 \v,24
	D8 \v,25
	D8 \v,26
	D8 \v,27
	D8 \v,28
	D8 \v,29
	D8 \v,30
	D8 \v,31
	D8 \v,32
	D8 \v,33
	D8 \v,34
	D8 \v,35
	D8 \v,36
	D8 \v,37
	D8 \v,38
	D8 \v,39
	D8 \v,40
	D8 \v,41
	D8 \v,42
	D8 \v,43
	D8 \v,44
	D8 \v,45
	D8 \v,46
	D8 \v,47
	D8 \v,48
	D8 \v,49
	D8 \v,50
	D8 \v,51
	D8 \v,52
	D8 \v,53
	D8 \v,54
	D8 \v,55
	D8 \v,56
	D8 \v,57
	D8 \v,58
	D8 \v,59
	D8 \v,60
	D8 \v,61
	D8 \v,62
	D8 \v,63
	DC\v 104,.Lg1\v
.Ldone\v:
	stp x5,x6,[x3,#32]
	stp x17,x19,[x3,#48]
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#192
	ret
.Lr1_w\v:
	RCHG
	sbcs x6,x6,xzr
	sbcs x17,x17,xzr
	sbc x19,x19,xzr
	b .Lb1_w\v
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
	E8 \v,16
	E8 \v,17
	E8 \v,18
	E8 \v,19
	E8 \v,20
	E8 \v,21
	E8 \v,22
	E8 \v,23
	E8 \v,24
	E8 \v,25
	E8 \v,26
	E8 \v,27
	E8 \v,28
	E8 \v,29
	E8 \v,30
	E8 \v,31
	E8 \v,32
	E8 \v,33
	E8 \v,34
	E8 \v,35
	E8 \v,36
	E8 \v,37
	E8 \v,38
	E8 \v,39
	E8 \v,40
	E8 \v,41
	E8 \v,42
	E8 \v,43
	E8 \v,44
	E8 \v,45
	E8 \v,46
	E8 \v,47
	E8 \v,48
	E8 \v,49
	E8 \v,50
	E8 \v,51
	E8 \v,52
	E8 \v,53
	E8 \v,54
	E8 \v,55
	E8 \v,56
	E8 \v,57
	E8 \v,58
	E8 \v,59
	E8 \v,60
	E8 \v,61
	E8 \v,62
	E8 \v,63
	S1 999,\v
	.endm
	.text
	.p2align 2
	PASS f
	PASS p

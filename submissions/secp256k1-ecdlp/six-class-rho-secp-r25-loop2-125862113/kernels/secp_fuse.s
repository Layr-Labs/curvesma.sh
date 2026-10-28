// Pass 1 with each walk's arrival fused in. See proof.md section 2.5. PASS f emits pass1f (compare rounds),
// PASS p pass1p (plain rounds: no K_cp load or compare; B.HS for B.NE).
	.macro KLf
	ldr x13,[x0],#8
	.endm
	.macro KLp
	.endm
	.macro FPf f,kc
	ccmp \f,\kc,#4,hs
	.endm
	.macro FPp f,kc
	.endm
	.macro BRf l
	b.ne \l
	.endm
	.macro BRp l
	b.hs \l
	.endm
	.macro GSTARTf
	ldr x13,[sp,#104]
	cbz x13,.Ldonef
	.endm
	.macro GSTARTp
	ldr x0,[sp,#104]
	cbz x0,.Ldonep
	.endm
	.macro GLENDf
	ldr x13,[sp,#104]
	subs x13,x13,#1
	str x13,[sp,#104]
	b.ne .Lg1f
	.endm
	.macro GLENDp
	subs x0,x0,#1
	b.ne .Lg1p
	.endm
	.macro B8 v,a,b,c,d,e,f,g,h
	BODY \a,\v
	BODY \b,\v
	BODY \c,\v
	BODY \d,\v
	BODY \e,\v
	BODY \f,\v
	BODY \g,\v
	BODY \h,\v
	.endm
	.macro RCHG
	ldr x13,[sp]
	sub x13,x13,#20
	str x13,[sp]
	.endm
	.macro CHAIN kc,f,v
	ccmn x25,#16,#2,ne
	ccmn x26,#17,#2,lo
	ccmn x30,x18,#2,lo
	ccmp x27,#16,#0,lo
	FP\v \f,\kc
	.endm
	.macro BODY k,v
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
	csel x24,x15,x14,hs
	adds x30,x30,x26
	csel x24,x16,x24,hs
	tbnz x28,#0,.Lo\v\k
	cmp x29,#0
	orr x29,x28,#1
	CHAIN x13,x29,\v
	BR\v .Lb\v\k
	ldr x13,[sp]
	sub x13,x13,#25
	str x13,[sp]
	ldr x13,[sp,#120]
	stp x13,x1,[x13]
	stp x28,x29,[x13,#16]
	add x13,x13,#32
	str xzr,[x13]
	str x13,[sp,#120]
	b .Lb\v\k
.Lo\v\k:
	subs x28,x12,x28
	ccmn x29,#1,#4,cs
	CHAIN x13,x28,\v
	BR\v .Lb\v\k
	ldr x13,[sp]
	sub x13,x13,#25
	str x13,[sp]
	ldr x13,[sp,#120]
	stp x13,x1,[x13]
	stp x28,x29,[x13,#16]
	add x13,x13,#32
	str xzr,[x13]
	str x13,[sp,#120]
.Lb\v\k:
	bfi x24,x28,#9,#26
	str x24,[x2],#8
	ldp x26,x27,[x24]
	ldp x29,x30,[x24,#16]
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
	csel x24,x15,x14,hs
	adds x30,x30,x26
	csel x24,x16,x24,hs
	tbnz x28,#0,.Low\v
	cmp x29,#0
	orr x29,x28,#1
	CHAIN x13,x29,\v
	BR\v .Lbw\v
	ldr x13,[sp]
	sub x13,x13,#25
	str x13,[sp]
	ldr x13,[sp,#120]
	stp x13,x1,[x13]
	stp x28,x29,[x13,#16]
	add x13,x13,#32
	str xzr,[x13]
	str x13,[sp,#120]
	b .Lbw\v
.Low\v:
	subs x28,x12,x28
	ccmn x29,#1,#4,cs
	CHAIN x13,x28,\v
	BR\v .Lbw\v
	ldr x13,[sp]
	sub x13,x13,#25
	str x13,[sp]
	ldr x13,[sp,#120]
	stp x13,x1,[x13]
	stp x28,x29,[x13,#16]
	add x13,x13,#32
	str xzr,[x13]
	str x13,[sp,#120]
.Lbw\v:
	bfi x24,x28,#9,#26
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
	.macro S8 v,a,b,c,d,e,f,g,h
	S1 \a,\v
	S1 \b,\v
	S1 \c,\v
	S1 \d,\v
	S1 \e,\v
	S1 \f,\v
	S1 \g,\v
	S1 \h,\v
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
	lsr x13,x13,#7
	str x13,[sp,#104]
	subs x13,x0,#1
	and x13,x13,#127
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
	ldr x13,[sp,#96]
	cbz x13,.Lg\v
2:
	BODY 128,\v
	ldr x13,[sp,#96]
	subs x13,x13,#1
	str x13,[sp,#96]
	b.ne 2b
.Lg\v:
	GSTART\v
.Lg1\v:
	B8 \v,0,1,2,3,4,5,6,7
	B8 \v,8,9,10,11,12,13,14,15
	B8 \v,16,17,18,19,20,21,22,23
	B8 \v,24,25,26,27,28,29,30,31
	B8 \v,32,33,34,35,36,37,38,39
	B8 \v,40,41,42,43,44,45,46,47
	B8 \v,48,49,50,51,52,53,54,55
	B8 \v,56,57,58,59,60,61,62,63
	B8 \v,64,65,66,67,68,69,70,71
	B8 \v,72,73,74,75,76,77,78,79
	B8 \v,80,81,82,83,84,85,86,87
	B8 \v,88,89,90,91,92,93,94,95
	B8 \v,96,97,98,99,100,101,102,103
	B8 \v,104,105,106,107,108,109,110,111
	B8 \v,112,113,114,115,116,117,118,119
	B8 \v,120,121,122,123,124,125,126,127
	GLEND\v
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
	S8 \v,0,1,2,3,4,5,6,7
	S8 \v,8,9,10,11,12,13,14,15
	S8 \v,16,17,18,19,20,21,22,23
	S8 \v,24,25,26,27,28,29,30,31
	S8 \v,32,33,34,35,36,37,38,39
	S8 \v,40,41,42,43,44,45,46,47
	S8 \v,48,49,50,51,52,53,54,55
	S8 \v,56,57,58,59,60,61,62,63
	S8 \v,64,65,66,67,68,69,70,71
	S8 \v,72,73,74,75,76,77,78,79
	S8 \v,80,81,82,83,84,85,86,87
	S8 \v,88,89,90,91,92,93,94,95
	S8 \v,96,97,98,99,100,101,102,103
	S8 \v,104,105,106,107,108,109,110,111
	S8 \v,112,113,114,115,116,117,118,119
	S8 \v,120,121,122,123,124,125,126,127
	S1 128,\v
	.endm
	.text
	.p2align 2
	PASS f
	PASS p

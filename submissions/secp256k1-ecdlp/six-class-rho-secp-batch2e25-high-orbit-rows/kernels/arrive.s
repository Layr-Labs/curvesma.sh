	.macro COMMON s
	ldp x24,x25,[x1,#16]
	ldp x17,x19,[x1,#32]
	ldp x22,x23,[x1],#64
	umulh x26,x22,x10
	umulh x27,x23,x9
	umulh x28,x24,x8
	umulh x29,x25,x7
	add x26,x26,x27
	add x28,x28,x29
	add x26,x26,x28
	madd x26,x22,x11,x26
	madd x26,x23,x10,x26
	madd x26,x24,x9,x26
	madd x26,x25,x8,x26
	add x27,x25,x26
	subs x28,x25,x26
	csel x29,x26,x25,hs
	csel \s,x14,x4,hs
	adds x22,x29,x27
	csel \s,x15,\s,hs
	.endm
	.macro CHAIN kc, f
	ccmn x26,#16,#2,ne
	ccmn x27,#17,#2,lo
	ccmn x22,x16,#2,lo
	ccmp x28,#16,#0,lo
	ccmp \f,\kc,#4,hs
	.endm
	.macro STUB kc, o
	sub x5,x2,x6
	ldr x29,[sp,#96]
	add x5,x29,x5
	mov x29,#1
	tbnz x17,#0,5f
	cmp x17,x12
	adc x24,x19,xzr
	cbz x24,1f
	b 7f
5:
	cmn x19,#1
	b.eq 1f
7:
	mov x29,#2
	cmn x26,#16
	ccmn x27,#17,#2,lo
	ccmn x22,x16,#2,lo
	ccmp x28,#16,#0,lo
	b.lo 1f
	mov x29,#3
	cmp x23,\kc
	b.eq 1f
	mov x29,#4
1:
	str x29,[x5,#\o]
	.endm
	.macro WALK s, kc, o
	COMMON \s
	tbnz x17,#0,6f
	cmp x19,x19,asr #63
	orr x23,x17,#1
	CHAIN \kc,x23
	b.ne 8f
	STUB \kc,\o
	b 8f
6:
	sub x23,x12,x17
	cmp x19,x19,asr #63
	CHAIN \kc,x23
	b.ne 8f
	STUB \kc,\o
8:
	bfi \s,x17,#9,#27
	.endm
	.macro PAIR
	ldp x20,x0,[x3],#16
	WALK x21,x20,0
	WALK x30,x0,8
	stp x21,x30,[x2],#16
	.endm
	.macro PAIR4
	PAIR
	PAIR
	PAIR
	PAIR
	.endm
	.text
	.globl arrive
	.p2align 2
arrive:
	stp x19,x20,[sp,#-96]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	ldp x7,x8,[x5]
	ldp x9,x10,[x5,#16]
	ldr x11,[x5,#32]
	ldp x12,x13,[x5,#40]
	ldr x16,[x5,#56]
	add x14,x4,#64
	add x15,x4,#128
	mov x6,x2
	tbz x0,#0,pairs
	ldr x20,[x3],#8
	WALK x30,x20,0
	str x30,[x2],#8
	tbnz x2,#28,done
pairs:
	and x18,x0,#62
	cbz x18,main
2:
	PAIR
	subs x18,x18,#2
	b.ne 2b
	tbnz x2,#28,done
main:
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	tbz x2,#28,main
done:
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	ret

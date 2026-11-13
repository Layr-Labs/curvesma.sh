	.text
	.p2align 2
	.globl dlv
dlv:
	ldp x1,x2,[x0]
	ldp x3,x4,[x0,#16]
	ldp x5,x6,[x0,#32]
	stp x3,x4,[x2]
	stp x5,x6,[x2,#16]
	ret
	.globl tinit
tinit:
1:
	str x2,[x1],#8
	subs x0,x0,#1
	b.ne 1b
	ret
	.globl fin
fin:
	add x3,x0,#64
	mov x4,#0
1:
	ldr x5,[x3],#8
	cmp x5,x6
	csel x5,x5,x6,lo
	stp x4,x2,[x0]
	stp x5,xzr,[x0,#16]
	stp xzr,xzr,[x0,#32]
	add x4,x4,#1
	subs x1,x1,#1
	b.ne 1b
	ret
	.macro HREC k2,k3,w2,w3
	bfi x23,\w3,#3,#36
	ldadd x24,x4,[x23]
	bfi \w3,x17,#0,#36
	stp \k2,\k3,[x4]
	stp \w2,\w3,[x4,#16]
	tbz x4,#12,0f
	bl h2flush
0:
	.endm
	.macro H2L
	ldp x6,x7,[x1],#16
	ldp x8,x9,[x1],#16
	HREC x6,x7,x8,x9
	.endm
	.macro H2P
	ld64b x6,[x1]
	add x1,x1,#64
	HREC x6,x7,x8,x9
	HREC x10,x11,x12,x13
	.endm
	.macro H2P2
	H2P
	H2P
	.endm
	.macro H2P4
	H2P2
	H2P2
	.endm
	.macro H2P8
	H2P4
	H2P4
	.endm
	.macro H2P16
	H2P8
	H2P8
	.endm
	.globl h2b
h2b:
	mov x28,x30
	mov x24,#32
1:
	ldr x2,[x19],#32
	mov x1,x20
	add x20,x20,x21
	cmp x1,x2
	b.hs 6f
	H2L
	sub x16,x2,x1
	and x22,x2,#0xffffffffffffffc0
	tbz x16,#6,2f
	H2P
2:
	tbz x16,#7,3f
	H2P2
3:
	tbz x16,#8,4f
	H2P4
4:
	tbz x16,#9,5f
	H2P8
5:
	tbz x16,#10,7f
	H2P16
7:
	cmp x1,x22
	b.hs 8f
9:
	H2P16
	H2P16
	cmp x1,x22
	b.lo 9b
8:
	tbz x16,#5,6f
	H2L
6:
	subs x3,x3,#1
	b.ne 1b
	mov x30,x28
	ret
h2flush:
	sub x5,x4,#4064
	str x5,[x23]
	ubfx x15,x23,#3,#36
	mov x14,#128
	stp x15,x5,[x0]
	stp x14,xzr,[x0,#16]
	ret
	.globl h2fin
h2fin:
	mov x9,#0
1:
	ldr x10,[x1]
	sub x11,x10,x3
	lsr x11,x11,#5
	str x3,[x1],#8
	orr x12,x9,x4
	stp x12,x3,[x0]
	stp x11,xzr,[x0,#16]
	add x3,x3,#8192
	add x9,x9,#1
	subs x2,x2,#1
	b.ne 1b
	ret
	.macro ZVA1
	dc zva, x1
	add x1, x1, #256
	.endm
	.macro ZVA8
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	.endm
	.macro ZVA64
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	.endm
	.macro ZVA512
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	.endm
	.globl zero16384
zero16384:
1:
	ZVA512
	subs x0,x0,#16384
	b.ne 1b
	ret
	.macro VERF e,z,k3,w2,w3,off
	add x8,x9,\e,lsl #1
	ldr x15,[x8,#8]
	eor x15,x15,\k3
	tst x15,#0xffffe0000fffffff
	b.ne 0f
	ldp x7,x15,[x8,#16]
	eor x7,x7,\w2
	eor x15,x15,\w3
	orr x7,x7,x15
	tst x7,#0xfffffffff
	ccmp x3,x4,#2,ne
	b.hs 0f
	add x8,x8,#16
	sub x15,x10,#\off
	stp x8,x15,[x3],#16
	mov \z,#0
0:
	ret
	.endm
hit5A:
	VERF x5,x6,x21,x22,x23,48
hit6A:
	VERF x6,x5,x21,x22,x23,48
hit5B:
	VERF x5,x6,x25,x26,x27,16
hit6B:
	VERF x6,x5,x25,x26,x27,16
	.macro JOIN k2,h5,h6
	bfi x12,\k2,#3,#58
	swp x14,x5,[x12]
	stp x5,\k2,[x14],#16
	cbz x5,11f
12:
	ldp x6,x7,[x5]
	cmp x7,\k2
	b.ne 13f
	bl \h5
13:
	cbz x6,11f
	ldp x5,x7,[x6]
	cmp x7,\k2
	b.ne 14f
	bl \h6
14:
	cbnz x5,12b
11:
	.endm
	.macro RP
	ld64b x20,[x1]
	add x1,x1,#64
	st64b x20,[x10]
	add x10,x10,#64
	JOIN x20,hit5A,hit6A
	JOIN x24,hit5B,hit6B
	.endm
	.macro RP2
	RP
	RP
	.endm
	.macro RP4
	RP2
	RP2
	.endm
	.macro RP8
	RP4
	RP4
	.endm
	.macro RP16
	RP8
	RP8
	.endm
	.macro RS
	ldp x24,x25,[x1],#16
	ldp x26,x27,[x1],#16
	stp x24,x25,[x10],#16
	stp x26,x27,[x10],#16
	JOIN x24,hit5B,hit6B
	.endm
	.globl rblk
rblk:
	ldp x18,x1,[x0]
	ldr x19,[x0,#16]
	mov x17,x30
	add x5,x10,x19,lsl #5
	add x5,x5,#32
	cmp x11,x5
	b.lo 8f
	cbz x19,8f
	RS
	sub x16,x19,#1
	tbz x16,#1,1f
	RP
1:
	tbz x16,#2,2f
	RP2
2:
	tbz x16,#3,3f
	RP4
3:
	tbz x16,#4,4f
	RP8
4:
	tbz x16,#5,5f
	RP16
5:
	tbz x16,#6,6f
	RP16
	RP16
6:
	tbz x16,#0,7f
	RS
7:
	tbz x19,#0,8f
	add x10,x10,#32
	add x14,x14,#16
8:
	ldr x8,[x0,#32]
	mov x15,#1
	stlr x15,[x8]
	mov x30,x17
	tbnz x18,#63,done_flag
	ret
done_flag:
	add x28,x28,#1
	cmp x28,x29
	b.eq final_done
	ret
final_done:
	stp x3,x10,[x2]
	stlr x15,[x13]
	mov x28,#0
	ret

// The arrival rows of proof.md section 3, by pairs of walks, with the pair loop unrolled
// two pairs (four walks) per end test: per 4 walks one LDP of two K_cp words per pair,
// 4 x 30 instructions, two STPs of table addresses, one loop TBZ. Walk 0 and an odd
// pair count run one pair through PAIR0. Walk A keeps its K_cp word in x20 and address
// in x21, walk B in x0 and x30. B.NE skips a test stub logging the rare case
// (1 distinguished, 2 window, 3 hit) at the log whose base the caller leaves at [sp].
// arrive(x0 = n, x1 = walk slots, x2 = pointers, &tp[n] the first with bit 16 set
// (23 in production), x3 = K_cp, x4 = table, x5 = constants)
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
	.macro CHAIN kc
	ccmn x26,#16,#2,ne
	ccmn x27,#17,#2,lo
	ccmn x22,x16,#2,lo
	ccmp x28,#16,#0,lo
	ccmp x23,\kc,#4,hs
	.endm
	.macro STUB k1, o
	sub x5,x2,x6
	ldr x29,[sp,#96]
	add x5,x29,x5
	mov x29,#1
	cbz \k1, 1f
	mov x29,#2
	cmn x26,#16
	ccmn x27,#17,#2,lo
	ccmn x22,x16,#2,lo
	ccmp x28,#16,#0,lo
	b.lo 1f
	mov x29,#3
1:
	str x29,[x5,#\o]
	.endm
	.macro EVEN s, kc, o
	cmp x19,#0
	orr x23,x17,#1
	CHAIN \kc
	b.ne 5f
	STUB x19,\o
5:
	bfi \s,x17,#9,#21
	.endm
	.macro ODD s, kc, o
	subs x23,x12,x17
	sbcs x24,x13,x19
	CHAIN \kc
	b.ne 5f
	STUB x24,\o
5:
	bfi \s,x23,#9,#21
	.endm

// One pair, then MAIN if the array continues; the odd-pair tail.
	.macro PAIR0
	ldp x20,x0,[x3],#16
	COMMON x21
	tbnz x17,#0, 10f
	EVEN x21,x20,0
	COMMON x30
	tbnz x17,#0, 11f
	EVEN x30,x0,8
	stp x21,x30,[x2],#16
	tbnz x2,#16,done
	b MAIN
11:
	ODD x30,x0,8
	stp x21,x30,[x2],#16
	tbnz x2,#16,done
	b MAIN
10:
	ODD x21,x20,0
	COMMON x30
	tbnz x17,#0, 12f
	EVEN x30,x0,8
	stp x21,x30,[x2],#16
	tbnz x2,#16,done
	b MAIN
12:
	ODD x30,x0,8
	stp x21,x30,[x2],#16
	tbnz x2,#16,done
	b MAIN
	.endm

// One pair, then MAIN if the array continues, else done. The second pair of the unrolled body.
	.macro PAIR2
	ldp x20,x0,[x3],#16
	COMMON x21
	tbnz x17,#0, 20f
	EVEN x21,x20,0
	COMMON x30
	tbnz x17,#0, 21f
	EVEN x30,x0,8
	stp x21,x30,[x2],#16
	tbz x2,#16,MAIN
	b done
21:
	ODD x30,x0,8
	stp x21,x30,[x2],#16
	tbz x2,#16,MAIN
	b done
20:
	ODD x21,x20,0
	COMMON x30
	tbnz x17,#0, 22f
	EVEN x30,x0,8
	stp x21,x30,[x2],#16
	tbz x2,#16,MAIN
	b done
22:
	ODD x30,x0,8
	stp x21,x30,[x2],#16
	tbz x2,#16,MAIN
	b done
	.endm

// The first pair of the unrolled body: each of its four tails falls through into a copy of PAIR2.
	.macro PAIR1
	ldp x20,x0,[x3],#16
	COMMON x21
	tbnz x17,#0, 30f
	EVEN x21,x20,0
	COMMON x30
	tbnz x17,#0, 31f
	EVEN x30,x0,8
	stp x21,x30,[x2],#16
	PAIR2
31:
	ODD x30,x0,8
	stp x21,x30,[x2],#16
	PAIR2
30:
	ODD x21,x20,0
	COMMON x30
	tbnz x17,#0, 32f
	EVEN x30,x0,8
	stp x21,x30,[x2],#16
	PAIR2
32:
	ODD x30,x0,8
	stp x21,x30,[x2],#16
	PAIR2
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
	tbz x0, #0, even_n
	ldr x20,[x3],#8
	COMMON x30
	tbnz x17, #0, w0odd
	EVEN x30,x20,0
	b w0end
w0odd:
	ODD x30,x20,0
w0end:
	str x30,[x2],#8
	tbnz x2, #16, done
even_n:
	tbnz x0, #1, one_pair
MAIN:
	PAIR1
one_pair:
	PAIR0
done:
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	ret

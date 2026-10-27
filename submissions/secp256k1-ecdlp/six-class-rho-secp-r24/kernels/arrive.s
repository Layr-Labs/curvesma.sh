// The arrival rows of proof.md section 3, by pairs of walks (a9297fb0) in a straight-line loop of 32 pairs: per
// pair one LDP of K_cp words, 2 x 30 instructions and one STP of table addresses; per 32 pairs one TBZ. Both
// parity paths of a walk end at one BFI: odd y keeps K0 + 1 in x17 and even y keeps y0 = K0 there, and the
// B.NE that skips the rare case lands on that BFI. A walk keeps its K_cp word in x20 (x0 in a pair's second
// walk) and its table address in x21 (x30). The test stub logs the rare case (1 distinguished, 2 window, 3 hit)
// at the log whose base the caller leaves at [sp]. arrive(x0 = n, x1 = walk slots, x2 = pointers, &tp[n] the
// first with bit 16 set (23 in production), x3 = K_cp, x4 = table, x5 = constants); n mod 64 is run first.
	.macro COMMON s
	ldp x24,x25,[x1,#16]  // x2, x3 limbs of x; F0 = x25
	ldp x17,x19,[x1,#32]  // y0, y1
	ldp x22,x23,[x1],#64  // x0, x1 limbs of x
	umulh x26,x22,x10  // hi(x0*G4)
	umulh x27,x23,x9  // hi(x1*G3)
	umulh x28,x24,x8  // hi(x2*G2)
	umulh x29,x25,x7  // hi(x3*G1)
	add x26,x26,x27
	add x28,x28,x29
	add x26,x26,x28
	madd x26,x22,x11,x26  // + lo(x0*G5)
	madd x26,x23,x10,x26  // + lo(x1*G4)
	madd x26,x24,x9,x26  // + lo(x2*G3)
	madd x26,x25,x8,x26  // + lo(x3*G2): F1
	add x27,x25,x26  // S = F0 + F1
	subs x28,x25,x26  // d01 = F0 - F1; HS iff F0 >= F1
	csel x29,x26,x25,hs  // m = min(F0, F1)
	csel \s,x14,x4,hs  // slot j = 1 or 0
	adds x22,x29,x27  // e = m + S; HS iff m >= 2^64 - S
	csel \s,x15,\s,hs  // slot j = 2
	.endm
	.macro CHAIN kc, f
	ccmn x26,#16,#2,ne
	ccmn x27,#17,#2,lo
	ccmn x22,x16,#2,lo
	ccmp x28,#16,#0,lo
	ccmp \f,\kc,#4,hs
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
	.macro WALK s, kc, o
	COMMON \s
	tbnz x17,#0,6f
	cmp x19,#0  // even y: K1 = y1, Z iff distinguished
	orr x23,x17,#1  // filter word K0 + 1
	CHAIN \kc,x23
	b.ne 8f
	STUB x19,\o
	b 8f
6:
	subs x17,x12,x17  // odd y: K0 + 1 = (p0 + 1) - y0, the filter word
	sbcs x24,x13,x19  // K1 = limb 1 of p - y: Z iff distinguished
	CHAIN \kc,x17
	b.ne 8f
	STUB x24,\o
8:
	bfi \s,x17,#9,#25  // bit 9 = s, bits 10..33 = h
	.endm
	.macro PAIR
	ldp x20,x0,[x3],#16  // the pair's K_cp words
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
	tbz x0,#0,pairs  // even n
	ldr x20,[x3],#8  // walk 0 alone
	WALK x30,x20,0
	str x30,[x2],#8
	tbnz x2,#16,done
pairs:
	and x18,x0,#62  // 2 * (pairs mod 32)
	cbz x18,main
2:
	PAIR
	subs x18,x18,#2
	b.ne 2b
	tbnz x2,#16,done
main:
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	PAIR4
	tbz x2,#16,main
done:
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	ret

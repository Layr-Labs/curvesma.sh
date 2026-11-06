// Record transport and home-node grouping (proof.md 2.1, 4). A record is 4 words: K2, K3, w2, w3. Port layout (message
// block, record blocks at 2^14 d with header +0 and slot n at +32n, k = 256; LC tables at 0 and 2^39, so 8d or
// 2^39 + 8d is d's entry and destination word; RC tables 2^43 above): proof.md 2.1, 2.3. Calls may overwrite x0-x29 and
// the flags; x30 and sp are preserved.
	.text
	.p2align 2
// Delivery of one one-record message (proof.md 2.1): x0 the message's 48 bytes; its four data words go to its address.
	.globl dlv
dlv:
	ldp x1,x2,[x0]
	ldp x3,x4,[x0,#16]
	ldp x5,x6,[x0,#32]
	stp x3,x4,[x2]
	stp x5,x6,[x2,#16]
	ret
// Delivery of one k-record message (proof.md 2.1): x0 its block. An LDP of the header steps to slot 1; then per record
// two LDPs and two STPs, odd slots up and even slots down, so x0 ends where it began.
	.macro DLV1 n
	ldp x3,x4,[x0],#16
	ldp x5,x6,[x0],#\n
	stp x3,x4,[x2],#16
	stp x5,x6,[x2],#\n
	.endm
	.macro DLV8 n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	.endm
	.macro DLV64 n,m
	DLV8 \n
	DLV8 \n
	DLV8 \n
	DLV8 \n
	DLV8 \n
	DLV8 \n
	DLV8 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \n
	DLV1 \m
	.endm
	.globl dlvk
dlvk:
	ldp x1,x2,[x0],#32
	DLV64 48,48
	DLV64 48,16
	DLV64 -80,-80
	DLV64 -80,-80
	ret
// beta tables (proof.md 3.1): x0 base (entry t of slot k at x0 + 256t + 32k), x1 the six slot deltas
// beta*2^(43k) mod p (32 bytes each), x2 entries per slot (2^43) except slot 5: x3 (2^41). Entry t+1 = entry t + delta mod p.
	.globl tbld
tbld:
	mov x7,#977
	movk x7,#1,lsl #32
	mov x8,#6
1:
	ldp x9,x10,[x1],#16
	ldp x11,x12,[x1],#16
	mov x13,xzr
	mov x14,xzr
	mov x15,xzr
	mov x16,xzr
	mov x17,x0
	cmp x8,#1
	csel x18,x3,x2,eq
2:
	stp x13,x14,[x17]
	stp x15,x16,[x17,#16]
	add x17,x17,#256
	adds x13,x13,x9
	adcs x14,x14,x10
	adcs x15,x15,x11
	adcs x16,x16,x12
	cset x19,hs
	adds x20,x13,x7
	adcs x21,x14,xzr
	adcs x22,x15,xzr
	adcs x23,x16,xzr
	cinc x19,x19,hs
	cmp x19,#0
	csel x13,x20,x13,ne
	csel x14,x21,x14,ne
	csel x15,x22,x15,ne
	csel x16,x23,x16,ne
	subs x18,x18,#1
	b.ne 2b
	add x0,x0,#32
	subs x8,x8,#1
	b.ne 1b
	ret
// RC table: x0 entries (>= 1), x1 table, x2 slot base; every entry gets the base.
	.globl tinit
tinit:
1:
	str x2,[x1],#8
	subs x0,x0,#1
	b.ne 1b
	ret
// Final blocks and cursors: x0 message block, x7 LC table, x3 RC table, x8 block 0's slot 1, x1 destinations (>= 1),
// x2 this sender's entry at each destination, x6 its slot end, x11 p. Destination d has L = (LC - slot 1)/32 records in its
// block and RC cursor c: if L > 0 and c < end the header (8d, c) sends the block; d receives min(c + 32L, end).
	.globl fin
fin:
	mov x4,#0
1:
	ldr x9,[x7],#8
	ldr x5,[x3],#8
	subs x10,x9,x8
	b.eq 2f
	cmp x5,x6
	b.hs 2f
	stp x10,x11,[x8,#-16]
	stp x4,x5,[x8,#-32]
2:
	add x5,x5,x10
	cmp x5,x6
	csel x5,x5,x6,lo
	stp x4,x2,[x0]
	stp x5,xzr,[x0,#16]
	stp xzr,xzr,[x0,#32]
	add x4,x4,#8
	add x8,x8,#16384
	subs x1,x1,#1
	b.ne 1b
	ret
// Hop-2 stream on block arrival (proof.md 4): x1 the block's first record, x12 the hop-2 slot end, x19 p; header +16
// the byte count, +24 p. Each record goes to column d2 = w3 mod 2^36 (BFI into bits 3..38 of x10 gives d2's row LC entry and
// destination word; LDADD returns the slot); after slot k the entry returns to slot 1, the RC cursor advances by 32k,
// and below the end the header sends the block, else it is dropped and the end stored back. Slot n is at window + 32n,
// so even slots are 64-byte aligned: slot 1 (and a last even slot) is read by two LDPs (H2R1), every other two slots
// by one LD64B into x20..x27 and an ADD (H2P, FEAT_LS64). A full block (count 8192) runs 256 unrolled; any other
// slot 1, then pairs 1, 2, then 4 at a time, then its last slot if unpaired.
	.macro H2B a,b,c,d
	bfi x10,\d,#3,#36
	ldadd x17,x11,[x10]
	stp \a,\b,[x11]
	stp \c,\d,[x11,#16]
	tbz x11,#13,9f
	sub x9,x11,#4080
	sub x9,x9,#4080
	str x9,[x10]
	orr x4,x10,#0x80000000000
	ldadd x18,x11,[x4]
	cmp x11,x12
	b.lo 8f
	str x12,[x4]
	b 9f
8:
	stp x18,x19,[x9,#-16]
	stp x10,x11,[x9,#-32]
9:
	.endm
	.macro H2R1
	ldp x20,x21,[x1],#16
	ldp x22,x23,[x1],#16
	H2B x20,x21,x22,x23
	.endm
	.macro H2P
	ld64b x20,[x1]
	add x1,x1,#64
	H2B x20,x21,x22,x23
	H2B x24,x25,x26,x27
	.endm
	.macro H2P8
	H2P
	H2P
	H2P
	H2P
	H2P
	H2P
	H2P
	H2P
	.endm
	.macro H2P64
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	.endm
// Run by the node's hop-1 delivery processor (proof.md 2.1): h2s0 once at node set-up (x12 the hop-2 slot end, kept),
// then h2s per arriving block as a loop with B for RET; x10, x18, x3 persist (BFI rewrites only bits 3..38 of x10).
// x1 is reset to the link window + 32 each run; p is loaded after the window wait (the count LDR) from the node word at
// 2^48 - 168, which main writes per group.
	.globl h2s0
h2s0:
	mov x10,#0x8000000000
	mov x18,#8192
	movn x3,#0xa7
	movk x3,#0,lsl #48
	ret
	.globl h2s
h2s:
	movn x1,#0x3fdf
	ldr x16,[x1,#-16]
	ldr x19,[x3]
	mov x17,#32
	tbnz x16,#13,7f
	cbz x16,6f
	sub x16,x16,#32
	add x2,x1,x16
	H2R1
	tbz x16,#6,3f
	H2P
3:
	tbz x16,#7,4f
	H2P
	H2P
4:
	cmp x1,x2
	b.hs 5f
1:
	H2P
	H2P
	H2P
	H2P
	cmp x1,x2
	b.lo 1b
5:
	tbz x16,#5,6f
	H2R1
6:
	ret
// A full block (count 8192, bit 13): slot 1, 127 pairs, slot 256, unrolled, no loop.
7:
	H2R1
	H2P64
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	H2P8
	H2P
	H2P
	H2P
	H2P
	H2P
	H2P
	H2P
	H2R1
	ret
// Grouped hash-chain join (proof.md 4), one record with K2 in \u, K3 in \v. x10 = 2^61, the head table: BFI of K2's
// low 58 bits into bits 3..60 is the head. C = 2^62 - (R2BASE>>1): a record's next entry (next, K2) is at E = (A>>1)+C,
// one per R2 slot position; heads and next fields hold entry addresses E. x13 = E (set per segment from x14 = C, stepped
// 16 by the STP); SWP makes it the head; its entry takes the old head and K2 before the walk. The walk alternates
// x11/x5, one LDP per visit (4 per visit, one B after an even count); a K2 match forms A = 2E - 2C (x15 = -2C), verifies
// d4 (ADD, LDR, EOR, TST, B.NE) and lists slot addresses unless full (4, or 2 if full; the record's own address is \a,
// which \l forms on the list path only: LNONE when \a = x2 is the record, LADD/LSUB x12 = x2 +- 32, one more);
// uncapped: each record links to its latest same-key predecessor.
	.macro LNONE
	.endm
	.macro LADD
	add x12,x2,#32
	.endm
	.macro LSUB
	sub x12,x2,#32
	.endm
	.macro HCHGI u,v,a,l
	bfi x10,\u,#3,#58
	swp x13,x11,[x10]
	stp x11,\u,[x13],#16
	cbz x11,8f
9:
	ldp x5,x6,[x11]
	cmp x6,\u
	b.eq 10f
12:
	cbz x5,8f
	ldp x11,x6,[x5]
	cmp x6,\u
	b.eq 11f
13:
	cbnz x11,9b
	b 8f
10:
	add x6,x15,x11,lsl #1
	ldr x12,[x6,#8]
	eor x12,x12,\v
	tst x12,#0x3fffffffffff
	b.ne 12b
	cmp x3,x4
	b.hs 8f
	\l
	stp x6,\a,[x3],#16
	b 8f
11:
	add x6,x15,x5,lsl #1
	ldr x12,[x6,#8]
	eor x12,x12,\v
	tst x12,#0x3fffffffffff
	b.ne 13b
	cmp x3,x4
	b.hs 8f
	\l
	stp x6,\a,[x3],#16
8:
	.endm
// Fused hop-2 arrival (proof.md 4): x0 the window (header c, count, p; F = *p). Copies every slot, joins [max(F, c),
// c + count) (HCHGI); join state at 2^48 - 160. x1 (window) and x2 (R2) point at the next slot; slot n is at x0 + 32n
// and c + 32(n - 1), both 64-byte aligned exactly for even n (c = 32 mod 64, proof.md 5). Two slots move by LD64B into
// x16..x23 from the window and ST64B into R2 (DLVP, with the two ADDs: 2 per slot); a slot without its pair by two
// LDPs and two STPs (DLVS, 4). Fresh full block: slot 1, 127 pairs, slot 256 unrolled (label 6); else three segments.
	.macro DLVS
	ldp x16,x17,[x1],#32
	ldp x18,x19,[x1,#-16]
	stp x18,x19,[x2,#16]
	stp x16,x17,[x2],#32
	.endm
	.macro DLVP
	ld64b x16,[x1]
	add x1,x1,#64
	st64b x16,[x2]
	add x2,x2,#64
	.endm
	.macro DLVHS
	DLVS
	HCHGI x16,x17,x12,LSUB
	.endm
	.macro DLVHP
	ld64b x16,[x1]
	add x1,x1,#64
	st64b x16,[x2]
	HCHGI x16,x17,x2,LNONE
	HCHGI x20,x21,x12,LADD
	add x2,x2,#64
	.endm
// A segment of x7 bytes (32 per slot), x8 = its window end - 32: if it is non-empty, an odd first slot alone, then
// pairs 1, 2, then 4 at a time, then an unpaired last slot.
	.macro SEG s,p,l
	cbz x7,sge\l
	tbz x1,#5,sga\l
	\s
	sub x7,x7,#32
sga\l:
	tbz x7,#6,sgb\l
	\p
sgb\l:
	tbz x7,#7,sgc\l
	\p
	\p
sgc\l:
	cmp x1,x8
	b.hs sgd\l
sgl\l:
	\p
	\p
	\p
	\p
	cmp x1,x8
	b.lo sgl\l
sgd\l:
	tbz x7,#5,sge\l
	\s
sge\l:
	.endm
	.macro DLVHP8
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	.endm
	.macro DLVHP64
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	.endm
	.globl dlvh
dlvh:
	ldp x24,x26,[x0,#8]
	ldr x13,[x0,#24]
	ldr x27,[x13]
	mov x9,#0x1000000000000
	ldp x3,x4,[x9,#-160]
	ldp x14,x15,[x9,#-144]
	mov x10,#0x2000000000000000
	add x1,x0,#32
	mov x2,x24
	cmp x24,x27
	b.lo 7f
	tbnz x26,#13,6f
7:
	csel x27,x27,x24,lo
	sub x7,x27,x24
	add x8,x0,x7
	SEG DLVS,DLVP,1
	add x25,x24,x26
	sub x7,x25,x27
	add x8,x0,x26
	add x13,x14,x2,lsr #1
	SEG DLVHS,DLVHP,2
	add x8,x0,#8192
	mov x7,#8192
	sub x7,x7,x26
	SEG DLVS,DLVP,3
	str x3,[x9,#-160]
	ret
6:
	add x13,x14,x2,lsr #1
	DLVHS
	DLVHP64
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP8
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHP
	DLVHS
	str x3,[x9,#-160]
	ret
// Zero n words (n a positive multiple of 16384; x1 256-byte aligned): x0 n, x1 base. A pass is 512 (DC ZVA, ADD #256),
// one 256-byte block each (DCZID_EL0.BS = 6, proof.md 2), SUBS and B.NE: 1026 per 16384 words (proof.md 7), plus the RET.
// Z8/Z64/Z512 only shorten the listing (a64.expand).
	.macro Z8
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	dc zva,x1
	add x1,x1,#256
	.endm
	.macro Z64
	Z8
	Z8
	Z8
	Z8
	Z8
	Z8
	Z8
	Z8
	.endm
	.macro Z512
	Z64
	Z64
	Z64
	Z64
	Z64
	Z64
	Z64
	Z64
	.endm
	.globl zero16384
zero16384:
1:
	Z512
	subs x0,x0,#16384		// one instruction (immediate 4, LSL #12)
	b.ne 1b
	ret
// LC table (init3): x0 entries, x1 table, x2 first value, x3 stride; entry t = first + t*stride.
	.globl init3
init3:
1:
	str x2,[x1],#8
	add x2,x2,x3
	subs x0,x0,#1
	b.ne 1b
	ret
// tbld8 (proof.md 3.1): rows 0..n-1 of an X table, slot k of row t = (c t 2^(32k) + add [k = 0]) mod p; x0 base,
// x1 the eight steps c 2^(32k) mod p, x2 add, x3 n; two slots per pass, a -> a + s mod p as T = a + (s + C0) and
// C0 off unless T carried; 92 + 96n instructions.
	.macro TSTEP
	ldp x11,x12,[x1],#16
	ldp x13,x14,[x1],#16
	ldp x19,x20,[x1],#16
	ldp x21,x22,[x1],#16
	adds x11,x11,x5
	adcs x12,x12,xzr
	adcs x13,x13,xzr
	adcs x14,x14,xzr
	adds x19,x19,x5
	adcs x20,x20,xzr
	adcs x21,x21,xzr
	adcs x22,x22,xzr
	.endm
	.macro TLOOP
	mov x4,x3
1:
	stp x7,x8,[x0,#208]
	stp x9,x10,[x0,#224]
	stp x15,x16,[x0,#240]
	stp x17,x18,[x0,#256]!
	adds x7,x7,x11
	adcs x8,x8,x12
	adcs x9,x9,x13
	adcs x10,x10,x14
	csel x6,xzr,x5,hs
	subs x7,x7,x6
	sbcs x8,x8,xzr
	sbcs x9,x9,xzr
	sbcs x10,x10,xzr
	adds x15,x15,x19
	adcs x16,x16,x20
	adcs x17,x17,x21
	adcs x18,x18,x22
	csel x6,xzr,x5,hs
	subs x15,x15,x6
	sbcs x16,x16,xzr
	sbcs x17,x17,xzr
	sbcs x18,x18,xzr
	subs x4,x4,#1
	b.ne 1b
	.endm
	.macro TNEXT
	sub x0,x0,x3,lsl #8
	add x0,x0,#64
	.endm
	.globl tbld8
tbld8:
	mov x5,#977
	movk x5,#1,lsl #32
	sub x0,x0,#208
	// slots 0, 1: slot 0 starts at add, slot 1 at 0
	TSTEP
	ldp x7,x8,[x2]
	ldp x9,x10,[x2,#16]
	mov x15,xzr
	mov x16,xzr
	mov x17,xzr
	mov x18,xzr
	TLOOP
	// slots 2, 3
	TNEXT
	TSTEP
	mov x7,xzr
	mov x8,xzr
	mov x9,xzr
	mov x10,xzr
	mov x15,xzr
	mov x16,xzr
	mov x17,xzr
	mov x18,xzr
	TLOOP
	// slots 4, 5
	TNEXT
	TSTEP
	mov x7,xzr
	mov x8,xzr
	mov x9,xzr
	mov x10,xzr
	mov x15,xzr
	mov x16,xzr
	mov x17,xzr
	mov x18,xzr
	TLOOP
	// slots 6, 7
	TNEXT
	TSTEP
	mov x7,xzr
	mov x8,xzr
	mov x9,xzr
	mov x10,xzr
	mov x15,xzr
	mov x16,xzr
	mov x17,xzr
	mov x18,xzr
	TLOOP
	ret

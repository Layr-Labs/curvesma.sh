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
// beta tables (proof.md 3.1): x0 base (entry t of window i at x0 + 256t + 32i), x1 the six deltas beta*2^(43i) mod p
// (32 bytes each), x2 entries per window for i < 5 (2^43), x3 for i = 5 (2^41). Entry t+1 = entry t + delta mod p.
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
// and below the end the header sends the block, else it is dropped and the end stored back. 1, 2, 4, then 8 at a time.
	.macro H2R
	ldp x5,x6,[x1],#16
	ldp x7,x8,[x1],#16
	bfi x10,x8,#3,#36
	ldadd x17,x11,[x10]
	stp x5,x6,[x11]
	stp x7,x8,[x11,#16]
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
	.globl h2s
h2s:
	mov x10,#0x8000000000
	mov x17,#32
	mov x18,#8192
	ldr x16,[x1,#-16]
	add x2,x1,x16
	tbz x16,#5,2f
	H2R
2:
	tbz x16,#6,3f
	H2R
	H2R
3:
	tbz x16,#7,4f
	H2R
	H2R
	H2R
	H2R
4:
	cmp x1,x2
	b.hs 6f
5:
	H2R
	H2R
	H2R
	H2R
	H2R
	H2R
	H2R
	H2R
	cmp x1,x2
	b.lo 5b
6:
	ret
// Grouped hash-chain join (proof.md 4), one record. x10 = 3*2^61, the head table: BFI of K2's low 57 bits into bits
// 3..59 is the head. x14 = C = 2^61 - (R2BASE>>1): a record's next entry (next address, K2) is at (A>>1)+C, one per R2
// slot position. SWP makes the record the head; its entry takes the old head and K2 before the walk. The walk
// alternates x11/x21, one LDP per visit (5 per visit, one B after an even count); a K2 match verifies d4 (LDR, EOR, TST,
// B.NE) and lists unless full (4 or 3, or 2); uncapped, so each record links to its latest same-key predecessor.
	.macro HCHGI a
	bfi x10,x7,#3,#57
	swp \a,x11,[x10]
	add x13,x14,\a,lsr #1
	stp x11,x7,[x13]
	cbz x11,8f
9:
	add x15,x14,x11,lsr #1
	ldp x21,x20,[x15]
	cmp x20,x7
	b.eq 10f
12:
	cbz x21,8f
	add x15,x14,x21,lsr #1
	ldp x11,x20,[x15]
	cmp x20,x7
	b.eq 11f
13:
	cbnz x11,9b
	b 8f
10:
	ldr x17,[x11,#8]
	eor x17,x17,x8
	tst x17,#0x3fffffffffff
	b.ne 12b
	cmp x3,x4
	b.hs 8f
	stp x11,\a,[x3],#16
	b 8f
11:
	ldr x17,[x21,#8]
	eor x17,x17,x8
	tst x17,#0x3fffffffffff
	b.ne 13b
	cmp x3,x4
	b.hs 8f
	stp x21,\a,[x3],#16
8:
	.endm
// Fused hop-2 arrival (proof.md 4): x0 the window (header c, count, p; F = *p). Copies every slot, joins [max(F, c),
// c + count) (HCHGI); join state at 2^48 - 160. Fresh full block: 256 unrolled (label 6); else three segments.
	.macro DLVC
	ldp x7,x8,[x1,#32]!
	ldp x5,x6,[x1,#16]
	stp x7,x8,[x2,#32]!
	stp x5,x6,[x2,#16]
	.endm
	.macro DLVH
	DLVC
	HCHGI x2
	.endm
	.macro SEG m,b,l
	tbz \b,#5,sga\l
	\m
sga\l:
	tbz \b,#6,sgb\l
	\m
	\m
sgb\l:
	tbz \b,#7,sgc\l
	\m
	\m
	\m
	\m
sgc\l:
	cmp x1,x22
	b.hs sge\l
sgd\l:
	\m
	\m
	\m
	\m
	\m
	\m
	\m
	\m
	cmp x1,x22
	b.lo sgd\l
sge\l:
	.endm
	.macro DLVH8
	DLVH
	DLVH
	DLVH
	DLVH
	DLVH
	DLVH
	DLVH
	DLVH
	.endm
	.macro DLVH64
	DLVH8
	DLVH8
	DLVH8
	DLVH8
	DLVH8
	DLVH8
	DLVH8
	DLVH8
	.endm
	.globl dlvh
dlvh:
	ldp x24,x16,[x0,#8]
	ldr x13,[x0,#24]
	ldr x18,[x13]
	mov x9,#0x1000000000000
	ldp x3,x4,[x9,#-160]
	ldr x14,[x9,#-144]
	mov x10,#0x6000000000000000
	mov x1,x0
	sub x2,x24,#32
	cmp x24,x18
	b.lo 7f
	tbnz x16,#13,6f
7:
	csel x18,x18,x24,lo
	sub x23,x18,x24
	add x22,x1,x23
	SEG DLVC,x23,1
	add x25,x24,x16
	sub x23,x25,x18
	add x22,x22,x23
	SEG DLVH,x23,2
	add x22,x0,#8192
	sub x23,x22,x1
	SEG DLVC,x23,3
	str x3,[x9,#-160]
	ret
6:
	DLVH64
	DLVH64
	DLVH64
	DLVH64
	str x3,[x9,#-160]
	ret
// Zero n words (n a positive multiple of 4096): x0 n, x1 base. A pass is 2048 post-indexed STPs of (xzr, xzr), SUBS
// and B.NE: 2050 per 4096 words (proof.md 7), plus the RET. Z8/Z64/Z512 only shorten the listing (a64.expand).
	.macro Z8
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
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
	.globl zero4096
zero4096:
1:
	Z512
	Z512
	Z512
	Z512
	subs x0,x0,#4096		// one instruction (immediate 1, LSL #12)
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

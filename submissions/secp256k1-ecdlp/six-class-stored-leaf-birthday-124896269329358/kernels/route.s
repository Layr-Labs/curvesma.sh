// Record transport and home-node grouping (proof.md 2.1, 4). A record is 4 words: K2, K3, w2, w3. A port has a 48-byte
// message block (destination word +0, its bits 3..38 the index; address +8; data +16..+47; the STP to +32 sends it, and
// dlv delivers it) and per destination d a record block at 2^13 d from the first: header (destination word, address)
// +0, slot n at +32n (n = 1..k, k = 128); the header's STP sends the block, and dlvk delivers it. LC tables (the next
// slot) at 0 and 2^39, so 8d or 2^39 + 8d is d's entry and a destination word for d; RC tables (the cursor in this
// sender's slot at d, which ends at the same address at every d) 2^43 above them. Calls may overwrite x0-x29 and the
// flags; x30 and sp are preserved.
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
	DLV64 48,16
	DLV64 -80,-80
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
// x2 this sender's entry at each destination, x6 its slot end. Destination d has L = (LC - slot 1)/32 records in its
// block and RC cursor c; if L > 0 and c < end, the header (8d, c) sends the block (c + 32k <= end, proof.md 4), and d
// receives (min(c + 32L, end), 0, 0, 0).
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
	stp x4,x5,[x8,#-32]
2:
	add x5,x5,x10
	cmp x5,x6
	csel x5,x5,x6,lo
	stp x4,x2,[x0]
	stp x5,xzr,[x0,#16]
	stp xzr,xzr,[x0,#32]
	add x4,x4,#8
	add x8,x8,#8192
	subs x1,x1,#1
	b.ne 1b
	ret
// Hop 2 stream: x3 slots (>= 1), x13 their final cursors (32-byte entries), x14 first slot, x15 slot bytes, x12 this
// sender's slot end at every destination. Each record goes to column d2 = w3 mod 2^36: BFI of d2 into bits 3..38 of x10
// (2^39 elsewhere) gives d2's row LC entry and destination word, whose LDADD returns the record's slot. After slot k
// (bits 5..11 zero) the entry returns to slot 1 and the RC cursor advances by 32k; if it was below the end the header
// (x10, cursor) sends the block, else the block is dropped and the end stored back. Records go in groups of 1, 2 and 4
// (the low bits of a slot's count), then 8 at a time.
	.macro H2R
	ldp x5,x6,[x1],#16
	ldp x7,x8,[x1],#16
	bfi x10,x8,#3,#36
	ldadd x17,x11,[x10]
	stp x5,x6,[x11]
	stp x7,x8,[x11,#16]
	tst x11,#0xfe0
	b.ne 9f
	sub x9,x11,#4064
	str x9,[x10]
	orr x4,x10,#0x80000000000
	ldadd x18,x11,[x4]
	cmp x11,x12
	b.lo 8f
	str x12,[x4]
	b 9f
8:
	stp x10,x11,[x9,#-32]
9:
	.endm
	.globl h2s
h2s:
	mov x10,#0x8000000000
	mov x17,#32
	mov x18,#4096
1:
	ldr x2,[x13],#32
	mov x1,x14
	add x14,x14,x15
	sub x16,x2,x1
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
	subs x3,x3,#1
	b.ne 1b
	ret
// LC table: x0 entries, x1 table, x2 first value (block 0's slot 1), x3 stride (2^13); entry t = first + t*stride.
	.globl init3
init3:
1:
	str x2,[x1],#8
	add x2,x2,x3
	subs x0,x0,#1
	b.ne 1b
	ret
// Chained grouping at the home (proof.md 2.3, 4). The head table, 2^57 words at a multiple of 2^60 in R1's dead range,
// is zeroed by h2zero; h2g then takes each slot of R2 up to its E2 entry, record by record. A record's cell is K3's low
// 57 bits: BFI writes them into bits 3..59 of x9. SWP puts the record's address in the cell and returns the old head,
// which STR writes over K3 as the record's link (phase 4 reads only w2, w3 and the address); the record then walks its
// cell's chain, every earlier record of the cell newest first, comparing K2. SUBS counts each step down from the budget
// x18, and CCMP compares only while the budget lasts (Z otherwise), so one B.EQ takes both a hit and the budget's end;
// the hit path stops on a spent budget or a full list, else appends the pair (earlier, this) and resumes the walk.
// Records go in groups of 1, 2 and 4, then 8 at a time.
// h2g: x0 slots (>= 1), x13 their final cursors, x14 first slot, x15 slot bytes, x9 the table, x18 the budget, x20 the
// candidate list, x21 its end; returns x0 = 0, or 1 when the budget or the list ended it.
	.macro H2W n
	ldp x7,x8,[x1,#32]!
	bfi x9,x8,#3,#57
	swp x1,x10,[x9]
	str x10,[x1,#8]
	cbz x10,hn\n
hw\n:
	ldp x11,x12,[x10]
	subs x18,x18,#1
	ccmp x11,x7,#4,hs
	b.eq ha\n
hv\n:
	cbz x12,hn\n
	ldp x11,x10,[x12]
	subs x18,x18,#1
	ccmp x11,x7,#4,hs
	b.eq hb\n
	cbnz x10,hw\n
hn\n:
	.endm
	.macro H2H n
ha\n:
	tbnz x18,#63,h2stop
	cmp x20,x21
	b.hs h2stop
	stp x10,x1,[x20],#16
	b hv\n
hb\n:
	tbnz x18,#63,h2stop
	cmp x20,x21
	b.hs h2stop
	stp x12,x1,[x20],#16
	cbnz x10,hw\n
	b hn\n
	.endm
	.globl h2g
h2g:
1:
	ldr x2,[x13],#32
	sub x1,x14,#32
	sub x16,x2,x14
	add x14,x14,x15
	sub x2,x2,#32
	tbz x16,#5,2f
	H2W 0
2:
	tbz x16,#6,3f
	H2W 1
	H2W 2
3:
	tbz x16,#7,4f
	H2W 3
	H2W 4
	H2W 5
	H2W 6
4:
	cmp x1,x2
	b.hs 6f
5:
	H2W 7
	H2W 8
	H2W 9
	H2W 10
	H2W 11
	H2W 12
	H2W 13
	H2W 14
	cmp x1,x2
	b.lo 5b
6:
	subs x0,x0,#1
	b.ne 1b
	ret
h2stop:
	mov x0,#1
	ret
	H2H 0
	H2H 1
	H2H 2
	H2H 3
	H2H 4
	H2H 5
	H2H 6
	H2H 7
	H2H 8
	H2H 9
	H2H 10
	H2H 11
	H2H 12
	H2H 13
	H2H 14
// Zero the head table: x0 its start, x1 its end (a multiple of 16 KiB above x0); 1024 STPs of xzr, CMP and B.NE per
// 16 KiB.
	.macro Z16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	stp xzr,xzr,[x0],#16
	.endm
	.macro Z256
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	Z16
	.endm
	.globl h2zero
h2zero:
1:
	Z256
	Z256
	Z256
	Z256
	cmp x0,x1
	b.ne 1b
	ret

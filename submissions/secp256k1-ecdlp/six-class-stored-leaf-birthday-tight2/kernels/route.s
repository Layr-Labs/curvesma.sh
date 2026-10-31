// Record transport and home-node grouping (proof.md section 4). A record is 4 words: K2, K3, w2, w3. A port is a 48-byte
// device register block: destination index +0, destination address +8, four data words +16..+47; the STP to +32
// completes the message, and the link delivers it: dlv below is that delivery's work, charged per message. A port's
// cursor table follows it at +64, 16 bytes per destination (cursor, slot end). Every other step is in these kernels.
// Calls may overwrite x0-x29 and the flags; x30 and sp are preserved.
	.text
	.p2align 2
// Delivery of one message (proof.md 2.1): x0 the message's 48 bytes; its four data words go to its address.
	.globl dlv
dlv:
	ldp x1,x2,[x0]
	ldp x3,x4,[x0,#16]
	ldp x5,x6,[x0,#32]
	stp x3,x4,[x2]
	stp x5,x6,[x2,#16]
	ret
// Cursor table: x0 entries (>= 1), x1 table, x2 slot base, x3 slot end; every entry gets (base, end).
	.globl tinit
tinit:
1:
	stp x2,x3,[x1],#16
	subs x0,x0,#1
	b.ne 1b
	ret
// Final cursors: x0 port (table at x0 + 64), x1 destinations (>= 1), x2 address of this sender's entry at each
// destination. Destination d receives (cursor of d, 0, 0, 0).
	.globl fin
fin:
	add x3,x0,#64
	mov x4,#0
1:
	ldr x5,[x3],#16
	stp x4,x2,[x0]
	stp x5,xzr,[x0,#16]
	stp xzr,xzr,[x0,#32]
	add x4,x4,#1
	subs x1,x1,#1
	b.ne 1b
	ret
// Hop 2 stream: x0 row port (table at x0 + 64), x3 slots (>= 1), x13 their final cursors (32-byte entries), x14 first
// slot, x15 slot bytes. Each record goes to the column d2 = w3 mod 2^36 at its cursor, or is dropped if that slot is full.
// A slot's records go in blocks of 1, 2 and 4 (the low bits of its count), then 8 at a time.
	.macro H2R
	ldp x5,x6,[x1],#16
	ldp x7,x8,[x1],#16
	and x9,x8,#0xfffffffff
	add x10,x0,x9,lsl #4
	ldp x11,x12,[x10,#64]
	cmp x11,x12
	b.hs 7f
	stp x9,x11,[x0]
	stp x5,x6,[x0,#16]
	stp x7,x8,[x0,#32]
	add x11,x11,#32
	str x11,[x10,#64]
7:
	.endm
	.globl h2s
h2s:
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
// Bucket pass at the home: x0 slots (>= 1), x13 their final cursors, x14 first slot, x15 slot bytes, x3 d3 mask, x4 d3
// cursor table (cur, end), x5 d4 mask, x6 d4 counters. Each record's (address, d4) goes to its d3 bucket, or the record
// is dropped if the bucket is full; a stored record's d4 is counted. x1 steps to each record's start before its loads,
// and x2 is the start of the slot's last record. Records go in blocks of 1, 2 and 4, then 8 at a time.
	.macro BKR
	ldp x7,x8,[x1,#32]!
	and x9,x7,x3
	add x9,x4,x9,lsl #4
	ldp x10,x11,[x9]
	cmp x10,x11
	b.hs 7f
	and x12,x8,x5
	stp x1,x12,[x10],#16
	str x10,[x9]
	ldr x10,[x6,x12,lsl #3]
	add x10,x10,#1
	str x10,[x6,x12,lsl #3]
7:
	.endm
	.globl bkt
bkt:
1:
	ldr x2,[x13],#32
	sub x1,x14,#32
	sub x16,x2,x14
	add x14,x14,x15
	sub x2,x2,#32
	tbz x16,#5,2f
	BKR
2:
	tbz x16,#6,3f
	BKR
	BKR
3:
	tbz x16,#7,4f
	BKR
	BKR
	BKR
	BKR
4:
	cmp x1,x2
	b.hs 6f
5:
	BKR
	BKR
	BKR
	BKR
	BKR
	BKR
	BKR
	BKR
	cmp x1,x2
	b.lo 5b
6:
	subs x0,x0,#1
	b.ne 1b
	ret
// Zero n counters (n a positive multiple of 4): x0 n, x1 base.
	.globl zero4
zero4:
1:
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	subs x0,x0,#4
	b.ne 1b
	ret
// Exclusive prefix of d4 counts into cursors (pointer-array addresses): x0 n (multiple of 4), x1 counts, x2 array base.
// Two counts are loaded with one LDP and two cursors stored with one STP, so a block of four costs 10 instructions.
	.globl prefix4
prefix4:
1:
	ldp x3,x4,[x1]
	ldp x6,x7,[x1,#16]
	add x5,x2,x3,lsl #3
	stp x2,x5,[x1],#16
	add x2,x5,x4,lsl #3
	add x8,x2,x6,lsl #3
	stp x2,x8,[x1],#16
	add x2,x8,x7,lsl #3
	subs x0,x0,#4
	b.ne 1b
	ret
// Local pass by d4: x0 buckets, x1 d3 cursor table (cur, end), x2 first bucket, x3 bucket bytes, x4 d4 cursors. Each
// bucket entry's record address goes to the pointer array at its d4 cursor (stable within a d4 value). A bucket's
// entries go in blocks of 1, 2 and 4, then 8 at a time.
	.macro P4R
	ldp x10,x8,[x7],#16
	ldr x9,[x4,x8,lsl #3]
	str x10,[x9],#8
	str x9,[x4,x8,lsl #3]
	.endm
	.globl pass4
pass4:
1:
	ldr x6,[x1],#16
	mov x7,x2
	add x2,x2,x3
	sub x11,x6,x7
	tbz x11,#4,2f
	P4R
2:
	tbz x11,#5,3f
	P4R
	P4R
3:
	tbz x11,#6,4f
	P4R
	P4R
	P4R
	P4R
4:
	cmp x7,x6
	b.hs 6f
5:
	P4R
	P4R
	P4R
	P4R
	P4R
	P4R
	P4R
	P4R
	cmp x7,x6
	b.lo 5b
6:
	subs x0,x0,#1
	b.ne 1b
	ret
// Scan: x0 n >= 2 pointers at x1; adjacent records with equal (K2, K3) go to the candidate list x2 (cursor) up to x3;
// a full list ends the scan.
	.globl scan
scan:
	ldr x9,[x1],#8
	ldp x5,x6,[x9]
	sub x0,x0,#1
1:
	ldr x10,[x1],#8
	ldp x7,x8,[x10]
	cmp x5,x7
	ccmp x6,x8,#0,eq
	b.eq 5f
2:
	subs x0,x0,#1
	b.eq 9f
	ldr x9,[x1],#8
	ldp x5,x6,[x9]
	cmp x7,x5
	ccmp x8,x6,#0,eq
	b.eq 6f
3:
	subs x0,x0,#1
	b.ne 1b
9:
	ret
5:
	cmp x2,x3
	b.hs 9b
	stp x9,x10,[x2],#16
	b 2b
6:
	cmp x2,x3
	b.hs 9b
	stp x10,x9,[x2],#16
	b 3b
// d3 cursor table: x0 buckets, x1 table, x2 first bucket, x3 bucket bytes; entry t = (base_t, base_t + bytes).
	.globl init3
init3:
1:
	add x4,x2,x3
	stp x2,x4,[x1],#16
	mov x2,x4
	subs x0,x0,#1
	b.ne 1b
	ret

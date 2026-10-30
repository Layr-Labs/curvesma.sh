// Record transport and home-node grouping (proof.md section 4). A record is 4 words: K2, K3, w2, w3.
// A window is 64 bytes: the record slot (+0..+31) and its full flag (+32). The link writes the slot, sets the flag and
// appends the window's index to the receiver's event queue. The handler processor runs h1recv or h2recv: it takes the
// next event, forms the window address, reads the record, runs the phase's receive body and clears the flag.
// Calls may overwrite x0-x29 and the flags; x30 and sp are preserved.
	.text
	.p2align 2
// Hop 1 receive (handler): x0 n events, x1 event queue, x2 buffer cursor, x3 buffer end, x4 window base. Returns x2.
	.globl h1recv
h1recv:
	cbz x0,9f
1:
	ldr x5,[x1],#8
	add x5,x4,x5,lsl #6
	ldp x6,x7,[x5]
	ldp x8,x9,[x5,#16]
	cmp x2,x3
	b.hs 2f
	stp x6,x7,[x2],#16
	stp x8,x9,[x2],#16
2:
	str xzr,[x5,#32]
	subs x0,x0,#1
	b.ne 1b
9:
	ret
// Hop 2 send (main): x0 n, x1 buffer, x2 row port (fixed 40-byte device register), x4 own column j.
// d2 goes out as the destination word; j replaces it in the record.
	.globl h2send
h2send:
	cbz x0,9f
1:
	ldp x5,x6,[x1],#16
	ldp x7,x8,[x1],#16
	and x9,x8,#0xfffffffff
	bfi x8,x4,#0,#36
	str x9,[x2]
	stp x5,x6,[x2,#8]
	stp x7,x8,[x2,#24]
	subs x0,x0,#1
	b.ne 1b
9:
	ret
// Hop 2 receive (handler): x0 n events, x1 event queue, x2 d3 cursor table (cur, end), x3 d3 mask, x4 d4 counts,
// x5 d4 mask, x15 window base.
	.globl h2recv
h2recv:
	cbz x0,9f
1:
	ldr x6,[x1],#8
	add x6,x15,x6,lsl #6
	ldp x7,x8,[x6]
	ldp x9,x10,[x6,#16]
	and x11,x7,x3
	add x11,x2,x11,lsl #4
	ldp x12,x13,[x11]
	cmp x12,x13
	b.hs 2f
	stp x7,x8,[x12],#16
	stp x9,x10,[x12],#16
	str x12,[x11]
	and x13,x8,x5
	add x13,x4,x13,lsl #3
	ldr x14,[x13]
	add x14,x14,#1
	str x14,[x13]
2:
	str xzr,[x6,#32]
	subs x0,x0,#1
	b.ne 1b
9:
	ret
// Window flags cleared (preprocessing): x0 windows (a positive multiple of 4), x1 window base.
	.globl wclear
wclear:
1:
	str xzr,[x1,#32]
	str xzr,[x1,#96]
	str xzr,[x1,#160]
	str xzr,[x1,#224]
	add x1,x1,#256
	subs x0,x0,#4
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
	.globl prefix4
prefix4:
1:
	ldr x3,[x1]
	str x2,[x1],#8
	add x2,x2,x3,lsl #3
	ldr x3,[x1]
	str x2,[x1],#8
	add x2,x2,x3,lsl #3
	ldr x3,[x1]
	str x2,[x1],#8
	add x2,x2,x3,lsl #3
	ldr x3,[x1]
	str x2,[x1],#8
	add x2,x2,x3,lsl #3
	subs x0,x0,#4
	b.ne 1b
	ret
// Local pass by d4: x0 bucket count, x1 d3 cursor table (cur, end), x2 bucket base, x3 bucket bytes, x4 d4 cursors,
// x5 d4 mask. Each record's address goes to the pointer array at its d4 cursor (stable within a d4 value).
	.globl pass4
pass4:
1:
	ldr x6,[x1],#16
	mov x7,x2
	add x2,x2,x3
	cmp x7,x6
	b.hs 3f
2:
	ldr x8,[x7,#8]
	and x8,x8,x5
	add x8,x4,x8,lsl #3
	ldr x9,[x8]
	str x7,[x9],#8
	str x9,[x8]
	add x7,x7,#32
	cmp x7,x6
	b.lo 2b
3:
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

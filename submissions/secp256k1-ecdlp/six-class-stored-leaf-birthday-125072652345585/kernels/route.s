// Record transport and home-node grouping (proof.md 2.1, 4). A record is 4 words: K2, K3, w2, w3. A port has a 48-byte
// message block (destination word +0, its bits 3..38 the index; address +8; data +16..+47; the STP to +32 sends it, and
// dlv delivers it) and per destination d a record block at 2^14 d from the first: header (destination word, address)
// +0, slot n at +32n (n = 1..k, k = 256); the header's STP sends the block, and dlvk delivers it. LC tables (the next
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
	DLV64 48,48
	DLV64 48,16
	DLV64 -80,-80
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
	add x8,x8,#16384
	subs x1,x1,#1
	b.ne 1b
	ret
// Hop 2 stream: x3 slots (>= 1), x13 their final cursors (32-byte entries), x14 first slot, x15 slot bytes, x12 this
// sender's slot end at every destination. Each record goes to column d2 = w3 mod 2^36: BFI of d2 into bits 3..38 of x10
// (2^39 elsewhere) gives d2's row LC entry and destination word, whose LDADD returns the record's slot. After slot k
// (bit 13 set, proof.md 4) the entry returns to slot 1 and the RC cursor advances by 32k; if it was below the end the header
// (x10, cursor) sends the block, else the block is dropped and the end stored back. Records go in groups of 1, 2 and 4
// (the low bits of a slot's count), then 8 at a time.
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
	stp x10,x11,[x9,#-32]
9:
	.endm
	.globl h2s
h2s:
	mov x10,#0x8000000000
	mov x17,#32
	mov x18,#8192
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
// Bucket pass at the home: x0 slots (>= 1), x13 their final cursors, x14 first slot, x15 slot bytes, x9 the d3 cursor
// table (at a multiple of 2^27), x11 the d4 counters (at a multiple of 2^49). BFI writes d3 = K2 mod 2^24 into bits
// 3..26 of x9 and d4 = K3 mod 2^46 into bits 3..48 of x11. Bucket t's entries end at a multiple of 2^38, so the cursor
// LDADD returns (advancing it by 16) has its low 38 bits zero exactly when the bucket is full: the advance is undone and
// the record dropped; else (address, d4 counter address) is stored there and d4 counted by STADD. x1 steps to each
// record before its loads; x2 is the slot's last record. Records go in groups of 1, 2 and 4, then 8 at a time.
	.macro BKR
	ldp x7,x8,[x1,#32]!
	bfi x9,x7,#3,#24
	ldadd x18,x10,[x9]
	tst x10,#0x3fffffffff
	b.ne 8f
	stadd x19,[x9]
	b 9f
8:
	bfi x11,x8,#3,#46
	stp x1,x11,[x10]
	stadd x17,[x11]
9:
	.endm
	.globl bkt
bkt:
	mov x17,#1
	mov x18,#16
	mov x19,#-16
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
// Local pass by d4: x0 buckets, x1 d3 cursor table, x2 first bucket's start, x3 bucket stride. Each bucket entry's
// record address goes to the pointer array at its d4 cursor, the counter whose address the entry holds, which LDADD
// (x20 = 8) advances (stable within a d4 value). A bucket's entries go in groups of 1, 2 and 4, then 8 at a time.
	.macro P4R
	ldp x10,x8,[x7],#16
	ldadd x20,x9,[x8]
	str x10,[x9]
	.endm
	.globl pass4
pass4:
	mov x20,#8
1:
	ldr x6,[x1],#8
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
// d3 cursor table, or LC table: x0 entries, x1 table, x2 first value (bucket start, or block 0's slot 1), x3 stride
// (bucket stride, or 2^14); entry t = first + t*stride.
	.globl init3
init3:
1:
	str x2,[x1],#8
	add x2,x2,x3
	subs x0,x0,#1
	b.ne 1b
	ret
// One comparison pair n of scan_ranges' 32-pair block: pointers p1, p2 (LDP), their K2s into b and c, CMP of the
// older pair (a, b) and CCMP of the newer (b, c), one B.EQ to the hit path.
	.macro PAIR n,p1,p2,a,b,c
pc_pair\n:
	ldp \p1,\p2,[x1],#16
	ldr \b,[\p1]
	ldr \c,[\p2]
	cmp \a,\b
	ccmp \b,\c,#4,ne
	b.eq pc_hit\n
pc_cont\n:
	.endm
// Hit path of pair n: re-test each comparison and append its equal pair (q, p1) or (p1, p2) after the list-end test.
	.macro HIT n,q,p1,p2,a,b,c
pc_hit\n:
	cmp \a,\b
	b.ne pc_second\n
	cmp x3,x4
	b.hs pc_done
	stp \q,\p1,[x3],#16
pc_second\n:
	cmp \b,\c
	b.ne pc_cont\n
	cmp x3,x4
	b.hs pc_done
	stp \p1,\p2,[x3],#16
	b pc_cont\n
	.endm
// Range-local scan (proof.md 4; the paired scanner of 892a3df8 and b00400fa on the range layout of 80e9e8af and
// 1efcef88): x0 d4 ranges, x1 the pointer array, x2 the d4 cursor table (each entry its range's end after pass4), x3
// the candidate list, x4 its end. In a range of n records: the first K2, one single comparison when n - 1 is odd, a
// five-level TBZ tree into the 32-pair block for ((n - 1)/2) mod 32 pairs, then the block per 64 records. A pair: LDP
// of two pointers, their K2s, CMP of the older pair, CCMP of the newer, one B.EQ; the hit path re-tests both and appends
// each equal pair's addresses after the list-end test, which ends the scan when full. No pair crosses a range.
	.globl scan_ranges
scan_ranges:
	cbz x0,pc_done
pc_bucket:
	ldr x18,[x2],#8
	cmp x1,x18
	b.hs pc_next
	ldr x16,[x1],#8
	ldr x5,[x16]
	sub x17,x18,x1
	lsr x17,x17,#3
	tbz x17,#0,pc_even
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq pc_odd_hit
pc_odd_cont:
	mov x16,x9
	mov x5,x7
pc_even:
	lsr x17,x17,#1
	and x17,x17,#31
	mov x7,x5
	mov x8,x5
	mov x10,x16
	mov x12,x16
	mov x14,x16
	tbz x17,#4,pc_dispatch_4_0
	tbz x17,#3,pc_dispatch_3_16
	tbz x17,#2,pc_dispatch_2_24
	tbz x17,#1,pc_dispatch_1_28
	tbz x17,#0,pc_pair2
	b pc_pair1
pc_dispatch_1_28:
	tbz x17,#0,pc_pair4
	b pc_pair3
pc_dispatch_2_24:
	tbz x17,#1,pc_dispatch_1_24
	tbz x17,#0,pc_pair6
	b pc_pair5
pc_dispatch_1_24:
	tbz x17,#0,pc_pair8
	b pc_pair7
pc_dispatch_3_16:
	tbz x17,#2,pc_dispatch_2_16
	tbz x17,#1,pc_dispatch_1_20
	tbz x17,#0,pc_pair10
	b pc_pair9
pc_dispatch_1_20:
	tbz x17,#0,pc_pair12
	b pc_pair11
pc_dispatch_2_16:
	tbz x17,#1,pc_dispatch_1_16
	tbz x17,#0,pc_pair14
	b pc_pair13
pc_dispatch_1_16:
	tbz x17,#0,pc_pair16
	b pc_pair15
pc_dispatch_4_0:
	tbz x17,#3,pc_dispatch_3_0
	tbz x17,#2,pc_dispatch_2_8
	tbz x17,#1,pc_dispatch_1_12
	tbz x17,#0,pc_pair18
	b pc_pair17
pc_dispatch_1_12:
	tbz x17,#0,pc_pair20
	b pc_pair19
pc_dispatch_2_8:
	tbz x17,#1,pc_dispatch_1_8
	tbz x17,#0,pc_pair22
	b pc_pair21
pc_dispatch_1_8:
	tbz x17,#0,pc_pair24
	b pc_pair23
pc_dispatch_3_0:
	tbz x17,#2,pc_dispatch_2_0
	tbz x17,#1,pc_dispatch_1_4
	tbz x17,#0,pc_pair26
	b pc_pair25
pc_dispatch_1_4:
	tbz x17,#0,pc_pair28
	b pc_pair27
pc_dispatch_2_0:
	tbz x17,#1,pc_dispatch_1_0
	tbz x17,#0,pc_pair30
	b pc_pair29
pc_dispatch_1_0:
	tbz x17,#0,pc_finish
	b pc_pair31
	PAIR 0,x9,x10,x5,x7,x8
	PAIR 1,x11,x12,x8,x5,x7
	PAIR 2,x13,x14,x7,x8,x5
	PAIR 3,x15,x16,x5,x7,x8
	PAIR 4,x9,x10,x8,x5,x7
	PAIR 5,x11,x12,x7,x8,x5
	PAIR 6,x13,x14,x5,x7,x8
	PAIR 7,x15,x16,x8,x5,x7
	PAIR 8,x9,x10,x7,x8,x5
	PAIR 9,x11,x12,x5,x7,x8
	PAIR 10,x13,x14,x8,x5,x7
	PAIR 11,x15,x16,x7,x8,x5
	PAIR 12,x9,x10,x5,x7,x8
	PAIR 13,x11,x12,x8,x5,x7
	PAIR 14,x13,x14,x7,x8,x5
	PAIR 15,x15,x16,x5,x7,x8
	PAIR 16,x9,x10,x8,x5,x7
	PAIR 17,x11,x12,x7,x8,x5
	PAIR 18,x13,x14,x5,x7,x8
	PAIR 19,x15,x16,x8,x5,x7
	PAIR 20,x9,x10,x7,x8,x5
	PAIR 21,x11,x12,x5,x7,x8
	PAIR 22,x13,x14,x8,x5,x7
	PAIR 23,x15,x16,x7,x8,x5
	PAIR 24,x9,x10,x5,x7,x8
	PAIR 25,x11,x12,x8,x5,x7
	PAIR 26,x13,x14,x7,x8,x5
	PAIR 27,x15,x16,x5,x7,x8
	PAIR 28,x9,x10,x8,x5,x7
	PAIR 29,x11,x12,x7,x8,x5
	PAIR 30,x13,x14,x5,x7,x8
	PAIR 31,x15,x16,x8,x5,x7
pc_finish:
	mov x5,x7
	cmp x1,x18
	b.lo pc_pair0
pc_next:
	subs x0,x0,#1
	b.ne pc_bucket
pc_done:
	ret
pc_odd_hit:
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
	b pc_odd_cont
	HIT 0,x16,x9,x10,x5,x7,x8
	HIT 1,x10,x11,x12,x8,x5,x7
	HIT 2,x12,x13,x14,x7,x8,x5
	HIT 3,x14,x15,x16,x5,x7,x8
	HIT 4,x16,x9,x10,x8,x5,x7
	HIT 5,x10,x11,x12,x7,x8,x5
	HIT 6,x12,x13,x14,x5,x7,x8
	HIT 7,x14,x15,x16,x8,x5,x7
	HIT 8,x16,x9,x10,x7,x8,x5
	HIT 9,x10,x11,x12,x5,x7,x8
	HIT 10,x12,x13,x14,x8,x5,x7
	HIT 11,x14,x15,x16,x7,x8,x5
	HIT 12,x16,x9,x10,x5,x7,x8
	HIT 13,x10,x11,x12,x8,x5,x7
	HIT 14,x12,x13,x14,x7,x8,x5
	HIT 15,x14,x15,x16,x5,x7,x8
	HIT 16,x16,x9,x10,x8,x5,x7
	HIT 17,x10,x11,x12,x7,x8,x5
	HIT 18,x12,x13,x14,x5,x7,x8
	HIT 19,x14,x15,x16,x8,x5,x7
	HIT 20,x16,x9,x10,x7,x8,x5
	HIT 21,x10,x11,x12,x5,x7,x8
	HIT 22,x12,x13,x14,x8,x5,x7
	HIT 23,x14,x15,x16,x7,x8,x5
	HIT 24,x16,x9,x10,x5,x7,x8
	HIT 25,x10,x11,x12,x8,x5,x7
	HIT 26,x12,x13,x14,x7,x8,x5
	HIT 27,x14,x15,x16,x5,x7,x8
	HIT 28,x16,x9,x10,x8,x5,x7
	HIT 29,x10,x11,x12,x7,x8,x5
	HIT 30,x12,x13,x14,x5,x7,x8
	HIT 31,x14,x15,x16,x8,x5,x7

// Record transport and home-node grouping (proof.md section 4). A record is 4 words: K2, K3, w2, w3. A port is a 48-byte
// device register block: destination index +0, destination address +8, four data words +16..+47; the STP to +32
// completes the message, and the link delivers it: dlv below is that delivery's work, charged per message. A port's
// cursor table follows it at +64, 8 bytes per destination (the sender's cursor in its slot there); the sender's slot
// ends at the same address at every destination, which the caller passes in a register. Every other step is in these
// kernels.
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
// Cursor table: x0 entries (>= 1), x1 table, x2 slot base; every entry gets the base.
	.globl tinit
tinit:
1:
	str x2,[x1],#8
	subs x0,x0,#1
	b.ne 1b
	ret
// Final cursors: x0 port (table at x0 + 64), x1 destinations (>= 1), x2 address of this sender's entry at each
// destination, x6 the end of this sender's slots. Destination d receives (min(cursor of d, end), 0, 0, 0): a hop-1
// cursor passes the end only by drops, which LDADD counts too.
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
// Hop 2 stream: x0 row port (table at x0 + 64), x3 slots (>= 1), x13 their final cursors (32-byte entries), x14 first
// slot, x15 slot bytes, x12 the end of this sender's slot at every destination. Each record goes to the column
// d2 = w3 mod 2^36 at the cursor that LDADD returns while advancing it by 32; if that cursor is at the end, the record is
// dropped and the end stored back, so a cursor never stays past the end. A slot's records go in blocks of 1, 2 and 4
// (the low bits of its count), then 8 at a time.
	.macro H2R
	ldp x5,x6,[x1],#16
	ldp x7,x8,[x1],#16
	and x9,x8,#0xfffffffff
	add x10,x0,x9,lsl #3
	ldadd x17,x11,[x10]
	cmp x11,x12
	b.lo 8f
	str x12,[x10]
	b 7f
8:
	stp x9,x11,[x0,#-64]
	stp x5,x6,[x0,#-48]
	stp x7,x8,[x0,#-32]
7:
	.endm
	.globl h2s
h2s:
	add x0,x0,#64
	mov x17,#32
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
// cursor table (8 bytes per bucket), x6 d4 counters; d3 is K2's masked low bits and d4 is bits 16..63 of K3. Bucket t's
// entries end at a multiple of 2^38, so the cursor that LDADD returns (advancing it by 16) has its low 38 bits zero
// exactly when the bucket is full; then the advance is undone and the record dropped. Otherwise the record's (address,
// d4) is stored there and its d4 counted with STADD. x1 steps to each record's start before its loads, and x2 is the
// start of the slot's last record. Records go in blocks of 1, 2 and 4, then 8 at a time.
	.macro BKR
	ldp x7,x8,[x1,#32]!
	and x9,x7,x3
	add x9,x4,x9,lsl #3
	ldadd x18,x10,[x9]
	tst x10,#0x3fffffffff
	b.ne 8f
	stadd x19,[x9]
	b 9f
8:
	lsr x12,x8,#16
	stp x1,x12,[x10]
	add x9,x6,x12,lsl #3
	stadd x17,[x9]
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
// Local pass by d4: x0 buckets, x1 d3 cursor table, x2 first bucket's start, x3 bucket stride, x4 d4 cursors. Each
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
// d3 cursor table: x0 buckets, x1 table, x2 first bucket's start, x3 bucket stride; entry t = start + t*stride.
	.globl init3
init3:
1:
	str x2,[x1],#8
	add x2,x2,x3
	subs x0,x0,#1
	b.ne 1b
	ret
// Range-local scan (proof.md section 4; the paired scanner of submissions 892a3df8 and b00400fa, on the range layout of
// 80e9e8af and 1efcef88): x0 the number of d4 ranges, x1 the pointer array, x2 the d4 cursor table (after pass4 each
// entry is its range's exclusive end), x3 the candidate list, x4 its end. Within a range of n records it loads the first
// record's K2, makes one single comparison when n - 1 is odd, enters a block of 32 comparison pairs through a five-level
// TBZ tree for the ((n - 1)/2) mod 32 remaining pairs, and then runs the block once per 64 records. A pair loads two
// pointers with one LDP and their two K2s, compares the older pair with CMP and the newer with CCMP (Z set if either is
// equal) and branches once; the out-of-line path re-tests both, and an equal K2 appends the two record addresses in
// order, after a test of the list's end that ends the scan when it is full. No pair crosses a range.
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
pc_pair0:
	ldp x9,x10,[x1],#16
	ldr x7,[x9]
	ldr x8,[x10]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit0
pc_cont0:
pc_pair1:
	ldp x11,x12,[x1],#16
	ldr x5,[x11]
	ldr x7,[x12]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit1
pc_cont1:
pc_pair2:
	ldp x13,x14,[x1],#16
	ldr x8,[x13]
	ldr x5,[x14]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit2
pc_cont2:
pc_pair3:
	ldp x15,x16,[x1],#16
	ldr x7,[x15]
	ldr x8,[x16]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit3
pc_cont3:
pc_pair4:
	ldp x9,x10,[x1],#16
	ldr x5,[x9]
	ldr x7,[x10]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit4
pc_cont4:
pc_pair5:
	ldp x11,x12,[x1],#16
	ldr x8,[x11]
	ldr x5,[x12]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit5
pc_cont5:
pc_pair6:
	ldp x13,x14,[x1],#16
	ldr x7,[x13]
	ldr x8,[x14]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit6
pc_cont6:
pc_pair7:
	ldp x15,x16,[x1],#16
	ldr x5,[x15]
	ldr x7,[x16]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit7
pc_cont7:
pc_pair8:
	ldp x9,x10,[x1],#16
	ldr x8,[x9]
	ldr x5,[x10]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit8
pc_cont8:
pc_pair9:
	ldp x11,x12,[x1],#16
	ldr x7,[x11]
	ldr x8,[x12]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit9
pc_cont9:
pc_pair10:
	ldp x13,x14,[x1],#16
	ldr x5,[x13]
	ldr x7,[x14]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit10
pc_cont10:
pc_pair11:
	ldp x15,x16,[x1],#16
	ldr x8,[x15]
	ldr x5,[x16]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit11
pc_cont11:
pc_pair12:
	ldp x9,x10,[x1],#16
	ldr x7,[x9]
	ldr x8,[x10]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit12
pc_cont12:
pc_pair13:
	ldp x11,x12,[x1],#16
	ldr x5,[x11]
	ldr x7,[x12]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit13
pc_cont13:
pc_pair14:
	ldp x13,x14,[x1],#16
	ldr x8,[x13]
	ldr x5,[x14]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit14
pc_cont14:
pc_pair15:
	ldp x15,x16,[x1],#16
	ldr x7,[x15]
	ldr x8,[x16]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit15
pc_cont15:
pc_pair16:
	ldp x9,x10,[x1],#16
	ldr x5,[x9]
	ldr x7,[x10]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit16
pc_cont16:
pc_pair17:
	ldp x11,x12,[x1],#16
	ldr x8,[x11]
	ldr x5,[x12]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit17
pc_cont17:
pc_pair18:
	ldp x13,x14,[x1],#16
	ldr x7,[x13]
	ldr x8,[x14]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit18
pc_cont18:
pc_pair19:
	ldp x15,x16,[x1],#16
	ldr x5,[x15]
	ldr x7,[x16]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit19
pc_cont19:
pc_pair20:
	ldp x9,x10,[x1],#16
	ldr x8,[x9]
	ldr x5,[x10]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit20
pc_cont20:
pc_pair21:
	ldp x11,x12,[x1],#16
	ldr x7,[x11]
	ldr x8,[x12]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit21
pc_cont21:
pc_pair22:
	ldp x13,x14,[x1],#16
	ldr x5,[x13]
	ldr x7,[x14]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit22
pc_cont22:
pc_pair23:
	ldp x15,x16,[x1],#16
	ldr x8,[x15]
	ldr x5,[x16]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit23
pc_cont23:
pc_pair24:
	ldp x9,x10,[x1],#16
	ldr x7,[x9]
	ldr x8,[x10]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit24
pc_cont24:
pc_pair25:
	ldp x11,x12,[x1],#16
	ldr x5,[x11]
	ldr x7,[x12]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit25
pc_cont25:
pc_pair26:
	ldp x13,x14,[x1],#16
	ldr x8,[x13]
	ldr x5,[x14]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit26
pc_cont26:
pc_pair27:
	ldp x15,x16,[x1],#16
	ldr x7,[x15]
	ldr x8,[x16]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit27
pc_cont27:
pc_pair28:
	ldp x9,x10,[x1],#16
	ldr x5,[x9]
	ldr x7,[x10]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit28
pc_cont28:
pc_pair29:
	ldp x11,x12,[x1],#16
	ldr x8,[x11]
	ldr x5,[x12]
	cmp x7,x8
	ccmp x8,x5,#4,ne
	b.eq pc_hit29
pc_cont29:
pc_pair30:
	ldp x13,x14,[x1],#16
	ldr x7,[x13]
	ldr x8,[x14]
	cmp x5,x7
	ccmp x7,x8,#4,ne
	b.eq pc_hit30
pc_cont30:
pc_pair31:
	ldp x15,x16,[x1],#16
	ldr x5,[x15]
	ldr x7,[x16]
	cmp x8,x5
	ccmp x5,x7,#4,ne
	b.eq pc_hit31
pc_cont31:
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
pc_hit0:
	cmp x5,x7
	b.ne pc_second0
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second0:
	cmp x7,x8
	b.ne pc_cont0
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont0
pc_hit1:
	cmp x8,x5
	b.ne pc_second1
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second1:
	cmp x5,x7
	b.ne pc_cont1
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont1
pc_hit2:
	cmp x7,x8
	b.ne pc_second2
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second2:
	cmp x8,x5
	b.ne pc_cont2
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont2
pc_hit3:
	cmp x5,x7
	b.ne pc_second3
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second3:
	cmp x7,x8
	b.ne pc_cont3
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont3
pc_hit4:
	cmp x8,x5
	b.ne pc_second4
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second4:
	cmp x5,x7
	b.ne pc_cont4
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont4
pc_hit5:
	cmp x7,x8
	b.ne pc_second5
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second5:
	cmp x8,x5
	b.ne pc_cont5
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont5
pc_hit6:
	cmp x5,x7
	b.ne pc_second6
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second6:
	cmp x7,x8
	b.ne pc_cont6
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont6
pc_hit7:
	cmp x8,x5
	b.ne pc_second7
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second7:
	cmp x5,x7
	b.ne pc_cont7
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont7
pc_hit8:
	cmp x7,x8
	b.ne pc_second8
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second8:
	cmp x8,x5
	b.ne pc_cont8
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont8
pc_hit9:
	cmp x5,x7
	b.ne pc_second9
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second9:
	cmp x7,x8
	b.ne pc_cont9
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont9
pc_hit10:
	cmp x8,x5
	b.ne pc_second10
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second10:
	cmp x5,x7
	b.ne pc_cont10
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont10
pc_hit11:
	cmp x7,x8
	b.ne pc_second11
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second11:
	cmp x8,x5
	b.ne pc_cont11
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont11
pc_hit12:
	cmp x5,x7
	b.ne pc_second12
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second12:
	cmp x7,x8
	b.ne pc_cont12
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont12
pc_hit13:
	cmp x8,x5
	b.ne pc_second13
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second13:
	cmp x5,x7
	b.ne pc_cont13
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont13
pc_hit14:
	cmp x7,x8
	b.ne pc_second14
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second14:
	cmp x8,x5
	b.ne pc_cont14
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont14
pc_hit15:
	cmp x5,x7
	b.ne pc_second15
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second15:
	cmp x7,x8
	b.ne pc_cont15
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont15
pc_hit16:
	cmp x8,x5
	b.ne pc_second16
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second16:
	cmp x5,x7
	b.ne pc_cont16
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont16
pc_hit17:
	cmp x7,x8
	b.ne pc_second17
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second17:
	cmp x8,x5
	b.ne pc_cont17
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont17
pc_hit18:
	cmp x5,x7
	b.ne pc_second18
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second18:
	cmp x7,x8
	b.ne pc_cont18
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont18
pc_hit19:
	cmp x8,x5
	b.ne pc_second19
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second19:
	cmp x5,x7
	b.ne pc_cont19
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont19
pc_hit20:
	cmp x7,x8
	b.ne pc_second20
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second20:
	cmp x8,x5
	b.ne pc_cont20
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont20
pc_hit21:
	cmp x5,x7
	b.ne pc_second21
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second21:
	cmp x7,x8
	b.ne pc_cont21
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont21
pc_hit22:
	cmp x8,x5
	b.ne pc_second22
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second22:
	cmp x5,x7
	b.ne pc_cont22
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont22
pc_hit23:
	cmp x7,x8
	b.ne pc_second23
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second23:
	cmp x8,x5
	b.ne pc_cont23
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont23
pc_hit24:
	cmp x5,x7
	b.ne pc_second24
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second24:
	cmp x7,x8
	b.ne pc_cont24
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont24
pc_hit25:
	cmp x8,x5
	b.ne pc_second25
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second25:
	cmp x5,x7
	b.ne pc_cont25
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont25
pc_hit26:
	cmp x7,x8
	b.ne pc_second26
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second26:
	cmp x8,x5
	b.ne pc_cont26
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont26
pc_hit27:
	cmp x5,x7
	b.ne pc_second27
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second27:
	cmp x7,x8
	b.ne pc_cont27
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont27
pc_hit28:
	cmp x8,x5
	b.ne pc_second28
	cmp x3,x4
	b.hs pc_done
	stp x16,x9,[x3],#16
pc_second28:
	cmp x5,x7
	b.ne pc_cont28
	cmp x3,x4
	b.hs pc_done
	stp x9,x10,[x3],#16
	b pc_cont28
pc_hit29:
	cmp x7,x8
	b.ne pc_second29
	cmp x3,x4
	b.hs pc_done
	stp x10,x11,[x3],#16
pc_second29:
	cmp x8,x5
	b.ne pc_cont29
	cmp x3,x4
	b.hs pc_done
	stp x11,x12,[x3],#16
	b pc_cont29
pc_hit30:
	cmp x5,x7
	b.ne pc_second30
	cmp x3,x4
	b.hs pc_done
	stp x12,x13,[x3],#16
pc_second30:
	cmp x7,x8
	b.ne pc_cont30
	cmp x3,x4
	b.hs pc_done
	stp x13,x14,[x3],#16
	b pc_cont30
pc_hit31:
	cmp x8,x5
	b.ne pc_second31
	cmp x3,x4
	b.hs pc_done
	stp x14,x15,[x3],#16
pc_second31:
	cmp x5,x7
	b.ne pc_cont31
	cmp x3,x4
	b.hs pc_done
	stp x15,x16,[x3],#16
	b pc_cont31

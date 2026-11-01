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
// cursor table (8 bytes per bucket), x5 d4 mask, x6 d4 counters. Bucket t's entries end at a multiple of 2^38, so the
// cursor that LDADD returns (advancing it by 16) has its low 38 bits zero exactly when the bucket is full; then the
// advance is undone and the record dropped. Otherwise the record's (address, d4) is stored there and its d4 counted with
// STADD. x1 steps to each record's start before its loads, and x2 is the start of the slot's last record. Records go in
// blocks of 1, 2 and 4, then 8 at a time.
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
	and x12,x8,x5
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
// Range-local scan (proof.md section 4; from submissions 80e9e8af and 1efcef88): x0 the number of d4 ranges, x1 the
// pointer array, x2 the d4 cursor table (after pass4 each entry is its range's exclusive end), x3 the candidate list,
// x4 its end. Within a range it compares the K2 of adjacent records: after (n - 1) mod 16 single comparisons, entered
// through a four-level TBZ tree, sixteen per end test with pointer pairs; an equal K2 appends the two record
// addresses, and an append at the list's end ends the scan. No pair crosses a range.
	.globl scan_ranges
scan_ranges:
	cbz x0,s16_done
s16_bucket:
	ldr x18,[x2],#8
	cmp x1,x18
	b.hs s16_next
	ldr x16,[x1],#8
	ldr x5,[x16]
	sub x17,x18,x1
	lsr x17,x17,#3
	and x17,x17,#15
	mov x9,x16
	mov x7,x5
	tbz x17,#3,s16_dispatch_0_3
	tbz x17,#2,s16_dispatch_8_2
	tbz x17,#1,s16_dispatch_12_1
	tbz x17,#0,s16_prefix2
	b s16_prefix1
s16_dispatch_12_1:
	tbz x17,#0,s16_prefix4
	b s16_prefix3
s16_dispatch_8_2:
	tbz x17,#1,s16_dispatch_8_1
	tbz x17,#0,s16_prefix6
	b s16_prefix5
s16_dispatch_8_1:
	tbz x17,#0,s16_prefix8
	b s16_prefix7
s16_dispatch_0_3:
	tbz x17,#2,s16_dispatch_0_2
	tbz x17,#1,s16_dispatch_4_1
	tbz x17,#0,s16_prefix10
	b s16_prefix9
s16_dispatch_4_1:
	tbz x17,#0,s16_prefix12
	b s16_prefix11
s16_dispatch_0_2:
	tbz x17,#1,s16_dispatch_0_1
	tbz x17,#0,s16_prefix14
	b s16_prefix13
s16_dispatch_0_1:
	tbz x17,#0,s16_norm
	b s16_prefix15
s16_prefix1:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit1
s16_pcont1:
s16_prefix2:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit2
s16_pcont2:
s16_prefix3:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit3
s16_pcont3:
s16_prefix4:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit4
s16_pcont4:
s16_prefix5:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit5
s16_pcont5:
s16_prefix6:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit6
s16_pcont6:
s16_prefix7:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit7
s16_pcont7:
s16_prefix8:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit8
s16_pcont8:
s16_prefix9:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit9
s16_pcont9:
s16_prefix10:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit10
s16_pcont10:
s16_prefix11:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit11
s16_pcont11:
s16_prefix12:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit12
s16_pcont12:
s16_prefix13:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit13
s16_pcont13:
s16_prefix14:
	ldr x16,[x1],#8
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_phit14
s16_pcont14:
s16_prefix15:
	ldr x9,[x1],#8
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_phit15
s16_pcont15:
s16_norm:
	mov x16,x9
	mov x5,x7
	cmp x1,x18
	b.hs s16_next
s16_loop:
	ldp x9,x10,[x1],#16
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_hit0
s16_cont0:
	ldr x5,[x10]
	cmp x7,x5
	b.eq s16_hit1
s16_cont1:
	ldp x11,x12,[x1],#16
	ldr x7,[x11]
	cmp x5,x7
	b.eq s16_hit2
s16_cont2:
	ldr x5,[x12]
	cmp x7,x5
	b.eq s16_hit3
s16_cont3:
	ldp x13,x14,[x1],#16
	ldr x7,[x13]
	cmp x5,x7
	b.eq s16_hit4
s16_cont4:
	ldr x5,[x14]
	cmp x7,x5
	b.eq s16_hit5
s16_cont5:
	ldp x15,x16,[x1],#16
	ldr x7,[x15]
	cmp x5,x7
	b.eq s16_hit6
s16_cont6:
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_hit7
s16_cont7:
	ldp x9,x10,[x1],#16
	ldr x7,[x9]
	cmp x5,x7
	b.eq s16_hit8
s16_cont8:
	ldr x5,[x10]
	cmp x7,x5
	b.eq s16_hit9
s16_cont9:
	ldp x11,x12,[x1],#16
	ldr x7,[x11]
	cmp x5,x7
	b.eq s16_hit10
s16_cont10:
	ldr x5,[x12]
	cmp x7,x5
	b.eq s16_hit11
s16_cont11:
	ldp x13,x14,[x1],#16
	ldr x7,[x13]
	cmp x5,x7
	b.eq s16_hit12
s16_cont12:
	ldr x5,[x14]
	cmp x7,x5
	b.eq s16_hit13
s16_cont13:
	ldp x15,x16,[x1],#16
	ldr x7,[x15]
	cmp x5,x7
	b.eq s16_hit14
s16_cont14:
	ldr x5,[x16]
	cmp x7,x5
	b.eq s16_hit15
s16_cont15:
	cmp x1,x18
	b.lo s16_loop
s16_next:
	subs x0,x0,#1
	b.ne s16_bucket
s16_done:
	ret
s16_phit1:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont1
s16_phit2:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont2
s16_phit3:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont3
s16_phit4:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont4
s16_phit5:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont5
s16_phit6:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont6
s16_phit7:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont7
s16_phit8:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont8
s16_phit9:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont9
s16_phit10:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont10
s16_phit11:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont11
s16_phit12:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont12
s16_phit13:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont13
s16_phit14:
	cmp x3,x4
	b.hs s16_done
	stp x9,x16,[x3],#16
	b s16_pcont14
s16_phit15:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_pcont15
s16_hit0:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_cont0
s16_hit1:
	cmp x3,x4
	b.hs s16_done
	stp x9,x10,[x3],#16
	b s16_cont1
s16_hit2:
	cmp x3,x4
	b.hs s16_done
	stp x10,x11,[x3],#16
	b s16_cont2
s16_hit3:
	cmp x3,x4
	b.hs s16_done
	stp x11,x12,[x3],#16
	b s16_cont3
s16_hit4:
	cmp x3,x4
	b.hs s16_done
	stp x12,x13,[x3],#16
	b s16_cont4
s16_hit5:
	cmp x3,x4
	b.hs s16_done
	stp x13,x14,[x3],#16
	b s16_cont5
s16_hit6:
	cmp x3,x4
	b.hs s16_done
	stp x14,x15,[x3],#16
	b s16_cont6
s16_hit7:
	cmp x3,x4
	b.hs s16_done
	stp x15,x16,[x3],#16
	b s16_cont7
s16_hit8:
	cmp x3,x4
	b.hs s16_done
	stp x16,x9,[x3],#16
	b s16_cont8
s16_hit9:
	cmp x3,x4
	b.hs s16_done
	stp x9,x10,[x3],#16
	b s16_cont9
s16_hit10:
	cmp x3,x4
	b.hs s16_done
	stp x10,x11,[x3],#16
	b s16_cont10
s16_hit11:
	cmp x3,x4
	b.hs s16_done
	stp x11,x12,[x3],#16
	b s16_cont11
s16_hit12:
	cmp x3,x4
	b.hs s16_done
	stp x12,x13,[x3],#16
	b s16_cont12
s16_hit13:
	cmp x3,x4
	b.hs s16_done
	stp x13,x14,[x3],#16
	b s16_cont13
s16_hit14:
	cmp x3,x4
	b.hs s16_done
	stp x14,x15,[x3],#16
	b s16_cont14
s16_hit15:
	cmp x3,x4
	b.hs s16_done
	stp x15,x16,[x3],#16
	b s16_cont15

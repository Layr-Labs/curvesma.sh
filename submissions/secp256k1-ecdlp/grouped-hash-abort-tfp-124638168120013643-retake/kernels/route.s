	.text
	.p2align 2
	.globl dlv
dlv:
	ldp x1,x2,[x0]
	ldp x3,x4,[x0,#16]
	ldp x5,x6,[x0,#32]
	stp x3,x4,[x2]
	stp x5,x6,[x2,#16]
	ret
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
	.globl tinit
tinit:
1:
	str x2,[x1],#8
	subs x0,x0,#1
	b.ne 1b
	ret
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
	str x10,[x8,#-16]
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
	str x18,[x9,#-16]
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
	.macro HCHG
	ldp x7,x8,[x1,#32]!
	bfi x10,x7,#3,#57
	swp x1,x11,[x10]
	add x12,x14,x1,lsr #1
	stp x11,x7,[x12]
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
	stp x11,x1,[x3],#16
	b 8f
11:
	ldr x17,[x21,#8]
	eor x17,x17,x8
	tst x17,#0x3fffffffffff
	b.ne 13b
	cmp x3,x4
	b.hs 8f
	stp x21,x1,[x3],#16
8:
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
	.globl hchg
hchg:
	mov x10,#0x6000000000000000
1:
	ldr x2,[x13],#32
	ldr x24,[x23],#8
	str x2,[x23,#-8]
	sub x1,x24,#32
	sub x16,x2,x24
	sub x2,x2,#32
	tbz x16,#5,2f
	HCHG
2:
	tbz x16,#6,3f
	HCHG
	HCHG
3:
	tbz x16,#7,4f
	HCHG
	HCHG
	HCHG
	HCHG
4:
	cmp x1,x2
	b.hs 6f
5:
	HCHG
	HCHG
	HCHG
	HCHG
	HCHG
	HCHG
	HCHG
	HCHG
	cmp x1,x2
	b.lo 5b
6:
	subs x0,x0,#1
	b.ne 1b
	ret
	.globl zero4
zero4:
1:
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	subs x0,x0,#4
	b.ne 1b
	ret
	.globl zero128
zero128:
1:
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	dc zva, x1
	add x1, x1, #64
	subs x0, x0, #128
	b.ne 1b
	ret
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
	.globl init3
init3:
1:
	str x2,[x1],#8
	add x2,x2,x3
	subs x0,x0,#1
	b.ne 1b
	ret
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

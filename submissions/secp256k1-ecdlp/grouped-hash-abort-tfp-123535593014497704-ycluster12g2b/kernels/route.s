// Transport, completion, grouping, unified slot-table builder (proof.md 2.1, 3.1, 4). 52-bit VAs.
	.text
	.p2align 2
	// 128-bit coordinator run gate: x0:x1 = u (2x 64-bit RNDR), x2:x3 = K_GATE -> x0 = 1 iff u < K_GATE (4 insns)
	.globl gate
gate:
	subs xzr,x0,x2
	sbcs xzr,x1,x3
	cset x0,lo
	ret
	.globl dlv
dlv:
	ldp x1,x2,[x0]
	ldp x3,x4,[x0,#16]
	ldp x5,x6,[x0,#32]
	stp x3,x4,[x2]
	stp x5,x6,[x2,#16]
	staddl x28,[x29]
	ret
// Unified slot-interleaved table builder (beta, 7-slot DEL, 14-slot AXYX, 2-slot PAIR49):
// x0 base, x1 descriptors, x2 full entries, x3 short entries, x4 slots, x5 stride, x6 short-slot bitmask.
// 4 fixed, 11 per slot, 20 per entry.
	.globl tbld
tbld:
	mov x7,#977
	movk x7,#1,lsl #32
	mov x8,x4
1:
	ldp x9,x10,[x1],#16
	ldp x11,x12,[x1],#16
	ldp x13,x14,[x1],#16
	ldp x15,x16,[x1],#16
	mov x17,x0
	tst x6,#1
	csel x18,x3,x2,ne
	lsr x6,x6,#1
2:
	stp x13,x14,[x17]
	stp x15,x16,[x17,#16]
	add x17,x17,x5
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
	.globl tinit
tinit:
1:
	str x2,[x1],#8
	subs x0,x0,#1
	b.ne 1b
	ret
	.globl fin1
fin1:
1:
	ldr x9,[x7]
	and x10,x9,#0xfffffffffffffc00
	sub x11,x9,x10
	sub x11,x11,#32
	ldr x5,[x10,#16]
	add x5,x5,#1
	stp x11,x5,[x10,#8]
	stlr x7,[x10]
	dmb ish
	str xzr,[x10,#16]
	add x9,x10,#32
	str x9,[x7],#8
	subs x1,x1,#1
	b.ne 1b
	ret
	.globl fin
fin:
1:
	ldr x9,[x7]
	ldr x5,[x3],#8
	subs x10,x9,x8
	b.eq 2f
	cmp x5,x6
	b.ge 2f
	stp x10,x11,[x8,#-16]
	str x5,[x8,#-24]
	stlur x7,[x8,#-32]
2:
	add x5,x5,x10
	cmp x5,x6
	csel x5,x5,x6,lt
	stp x7,x2,[x0]
	stp x5,xzr,[x0,#16]
	str xzr,[x0,#32]
	stlur xzr,[x0,#40]
	dmb ish
	add x7,x7,#8
	add x8,x8,#1024
	subs x1,x1,#1
	b.ne 1b
	ret
	.macro H2B a,b,c,d
	bfi x10,\a,#3,#42
	ldadd x17,x11,[x10]
	stp \a,\b,[x11]
	stp \c,\d,[x11,#16]
	tbz x11,#9,9f
	sub x9,x11,#480
	str x9,[x10]
	orr x4,x10,#0x400000000000
	ldadd x18,x11,[x4]
	cmp x11,x12
	b.lt 8f
	str x12,[x4]
	b 9f
8:
	stp x18,x19,[x9,#-16]
	str x11,[x9,#-24]
	stlur x10,[x9,#-32]
	dmb ish
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
	.macro H2PR w
	ld64b x20,[\w]
	H2B x20,x21,x22,x23
	H2B x24,x25,x26,x27
	.endm
	.macro H2R1F b,o1,o2
	ldp x20,x21,[\b,#\o1]
	ldp x22,x23,[\b,#\o2]
	H2B x20,x21,x22,x23
	.endm
	.globl h2s0
h2s0:
	mov x10,#0x200000000000
	mov x18,#512
	mov x17,#32
	movn x3,#0xa7
	movk x3,#0,lsl #48
	stp x8,x0,[x3,#-16]
	movn x29,#0xfff
	movk x29,#0,lsl #48
	mov x6,#0x4000000000000000
	mov x13,xzr
	mov x14,x6
	mov x15,x8
	mov x5,x0
	mov x28,xzr
	mov x0,#0xc040
	mov x1,#0xc080
	mov x2,#0xc0c0
	mov x8,#0xc100
	ret
	.globl h2s
h2s:
	ldp x16,x7,[x0,#-56]
	eor x9,x16,x16
	ldr x19,[x3,x9]
	tbnz x16,#9,7f
	mov x1,#0xc020
	adds x13,x13,x7
	adc x14,x14,xzr
	subs x15,x15,#1
	sbc x5,x5,xzr
	orr x9,x15,x5
	cbnz x9,2f
	sub x14,x14,x6
2:
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
	mov x1,#0xc080
	mov x2,#0xc0c0
	subs x13,x13,#1
	sbc x14,x14,xzr
	orr x9,x13,x14
	cbz x9,h2sd
	ret
7:
	H2R1F x0,-32,-16
	H2PR x0
	H2PR x1
	H2PR x2
	H2PR x8
	add x7,x0,#256
	H2PR x7
	add x16,x1,#256
	H2PR x16
	add x7,x2,#256
	H2PR x7
	H2R1F x8,256,272
	subs x13,x13,#1
	sbc x14,x14,xzr
	orr x9,x13,x14
	cbz x9,h2sd
	ret
h2sd:
	add x28,x28,#1
	stlr x28,[x29]
	mov x14,x6
	ldp x15,x5,[x3,#-16]
	ret
// Hash join (HB = 46, heads at 2^64 - 3*2^49 = 0xfffa<<48), out-of-line hit blocks (.pushsection, VISIT 4, HCHG_END 0).
	.macro LNONE
	.endm
	.macro LADD
	add x12,x2,#32
	.endm
	.macro LSUB
	sub x12,x2,#32
	.endm
	.macro HCHGI u,v,a,l
	bfi x10,\u,#3,#46
	swp x13,x11,[x10]
	stp x11,\u,[x13],#16
	cbz x11,he\@
hv\@:
	ldp x5,x6,[x11]
	cmp x6,\u
	b.eq hx\@
hy\@:
	cbz x5,he\@
	ldp x11,x6,[x5]
	cmp x6,\u
	b.eq hz\@
hw\@:
	cbnz x11,hv\@
	.pushsection .text.hit
hx\@:
	add x6,x15,x11,lsl #1
	ldr x12,[x6]
	eor x12,x12,\v
	tst x12,#0xfffffc0000000000
	b.ne hy\@
	cmp x3,x4
	b.hs he\@
	\l
	stp x6,\a,[x3],#16
	b he\@
hz\@:
	add x6,x15,x5,lsl #1
	ldr x12,[x6]
	eor x12,x12,\v
	tst x12,#0xfffffc0000000000
	b.ne hw\@
	cmp x3,x4
	b.hs he\@
	\l
	stp x6,\a,[x3],#16
	b he\@
	.popsection
he\@:
	.endm
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
	HCHGI x17,x16,x12,LSUB
	.endm
	.macro DLVHP
	ld64b x16,[x1]
	add x1,x1,#64
	st64b x16,[x2]
	HCHGI x17,x16,x2,LNONE
	HCHGI x21,x20,x12,LADD
	add x2,x2,#64
	.endm
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
	.macro DLVHSF b,o1,o2
	ldp x16,x17,[\b,#\o1]
	ldp x18,x19,[\b,#\o2]
	stp x18,x19,[x2,#16]
	stp x16,x17,[x2],#32
	HCHGI x17,x16,x12,LSUB
	.endm
	.macro DLVHPR w
	ld64b x16,[\w]
	st64b x16,[x2]
	HCHGI x17,x16,x2,LNONE
	HCHGI x21,x20,x12,LADD
	add x2,x2,#64
	.endm
	.globl dlvh0
dlvh0:
	mov x0,#0xc000
	mov x1,#0xc040
	mov x7,#0xc080
	mov x8,#0xc0c0
	mov x25,#0xc100
	mov x28,#1
	ret
	.globl dlvh
dlvh:
	ldp x24,x26,[x0,#8]
	ldr x13,[x0,#24]
	ldr x27,[x13]
	mov x9,#0x10000000000000
	ldp x3,x4,[x9,#-160]
	ldp x14,x15,[x9,#-144]
	mov x10,#0xfffa000000000000
	mov x2,x24
	cmp x24,x27
	b.lt 7f
	tbnz x26,#9,6f
7:
	add x1,x0,#32
	csel x27,x27,x24,lt
	sub x7,x27,x24
	add x8,x0,x7
	SEG DLVS,DLVP,1
	add x25,x24,x26
	sub x7,x25,x27
	add x8,x0,x26
	add x13,x14,x2,asr #1
	SEG DLVHS,DLVHP,2
	add x8,x0,#512
	mov x7,#512
	sub x7,x7,x26
	SEG DLVS,DLVP,3
	mov x1,#0xc040
	mov x7,#0xc080
	mov x8,#0xc0c0
	mov x25,#0xc100
	str x3,[x9,#-160]
	staddl x28,[x29]
	ret
6:
	add x13,x14,x2,asr #1
	DLVHSF x0,32,48
	DLVHPR x1
	DLVHPR x7
	DLVHPR x8
	DLVHPR x25
	add x24,x1,#256
	DLVHPR x24
	add x26,x7,#256
	DLVHPR x26
	add x27,x8,#256
	DLVHPR x27
	DLVHSF x27,64,80
	str x3,[x9,#-160]
	staddl x28,[x29]
	ret
// zero16384 (DCZID_EL0.BS = 6, 256-byte blocks): 512 (DC ZVA, ADD #256), SUBS, B.NE (1026 per 16384 words), RET.
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
	subs x0,x0,#16384
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
	.globl dlvq
dlvq:
	ldp x1,x2,[x0]
	ldp x3,x4,[x0,#16]
	ldp x5,x6,[x0,#32]
	ldadd x9,x8,[x10]
	cmp x8,x11
	b.hs 1f
	stp x2,x3,[x8]
	stp x4,x5,[x8,#16]
	str x6,[x8,#32]
1:
	staddl x28,[x29]
	ret
	.globl reqh
reqh:
1:
	ldp x14,x15,[x12],#16
	ldp x2,x3,[x14,#16]
	stp xzr,xzr,[x0]
	stp x2,x3,[x0,#16]
	str x10,[x0,#32]
	stlur x13,[x0,#40]
	dmb ish
	add x13,x13,#256
	ldp x2,x3,[x15,#16]
	stp xzr,xzr,[x0]
	stp x2,x3,[x0,#16]
	str x10,[x0,#32]
	stlur x13,[x0,#40]
	dmb ish
	add x13,x13,#256
	subs x9,x9,#1
	b.ne 1b
	ret
	.globl fwdq
fwdq:
	cmp x1,x6
	b.hs 2f
1:
	ldp x2,x3,[x1,#8]
	and x4,x3,#15
	lsl x5,x1,#8
	stp x2,x4,[x0]
	stp x3,x5,[x0,#16]
	str x7,[x0,#32]
	stlur xzr,[x0,#40]
	dmb ish
	add x1,x1,#64
	cmp x1,x6
	b.lo 1b
2:
	ret
	.globl sendr
sendr:
1:
	ldp x2,x3,[x1],#16
	ldp x4,x5,[x1],#16
	stp x7,x8,[x0]
	stp x2,x3,[x0,#16]
	str x4,[x0,#32]
	stlur x5,[x0,#40]
	dmb ish
	add x8,x8,#32
	subs x9,x9,#1
	b.ne 1b
	ret
	.globl fwdr
fwdr:
	cmp x1,x6
	b.hs 2f
1:
	ldp x2,x3,[x1]
	ldp x4,x5,[x1,#16]
	ldr x11,[x1,#32]
	lsr x7,x2,#8
	and x2,x2,#255
	ldp x8,x9,[x7,#24]
	add x9,x9,x2
	stp x8,x9,[x0]
	stp x3,x4,[x0,#16]
	str x5,[x0,#32]
	stlur x11,[x0,#40]
	dmb ish
	add x1,x1,#64
	cmp x1,x6
	b.lo 1b
2:
	ret
	.globl ydlv0
ydlv0:
	mov x0,#0xc040
	add x1,x0,#64
	add x2,x0,#128
	add x3,x0,#192
	add x4,x0,#256
	add x5,x0,#320
	add x6,x0,#384
	mov x14,#0xb000
	add x7,x14,#64
	add x8,x14,#128
	add x9,x14,#192
	add x10,x14,#256
	add x11,x14,#320
	add x12,x14,#384
	add x13,x14,#448
	mov x28,#1
	ret
	.globl ysnd
ysnd:
	ld64b x16,[x15]
	st64b x16,[x7]
	add x15,x15,#64
	ld64b x16,[x15]
	st64b x16,[x8]
	add x15,x15,#64
	ld64b x16,[x15]
	st64b x16,[x9]
	add x15,x15,#64
	ld64b x16,[x15]
	st64b x16,[x10]
	add x15,x15,#64
	ld64b x16,[x15]
	st64b x16,[x11]
	add x15,x15,#64
	ld64b x16,[x15]
	st64b x16,[x12]
	add x15,x15,#64
	ld64b x16,[x15]
	st64b x16,[x13]
	add x15,x15,#128
	stlr x28,[x14]
	dmb ish
	ret
	.globl ydlv
ydlv:
	ld64b x16,[x0]
	st64b x16,[x15]
	st64b x16,[x7]
	add x15,x15,#64
	ld64b x16,[x1]
	st64b x16,[x15]
	st64b x16,[x8]
	add x15,x15,#64
	ld64b x16,[x2]
	st64b x16,[x15]
	st64b x16,[x9]
	add x15,x15,#64
	ld64b x16,[x3]
	st64b x16,[x15]
	st64b x16,[x10]
	add x15,x15,#64
	ld64b x16,[x4]
	st64b x16,[x15]
	st64b x16,[x11]
	add x15,x15,#64
	ld64b x16,[x5]
	st64b x16,[x15]
	st64b x16,[x12]
	add x15,x15,#64
	ld64b x16,[x6]
	st64b x16,[x15]
	st64b x16,[x13]
	add x15,x15,#128
	stlr x28,[x14]
	dmb ish
	staddl x28,[x29]
	ret
	.globl ydlve
ydlve:
	ld64b x16,[x0]
	st64b x16,[x15]
	add x15,x15,#64
	ld64b x16,[x1]
	st64b x16,[x15]
	add x15,x15,#64
	ld64b x16,[x2]
	st64b x16,[x15]
	add x15,x15,#64
	ld64b x16,[x3]
	st64b x16,[x15]
	add x15,x15,#64
	ld64b x16,[x4]
	st64b x16,[x15]
	add x15,x15,#64
	ld64b x16,[x5]
	st64b x16,[x15]
	add x15,x15,#64
	ld64b x16,[x6]
	st64b x16,[x15]
	add x15,x15,#128
	staddl x28,[x29]
	ret

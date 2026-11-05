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
	.globl tinit
tinit:
1:
	str x2,[x1],#8
	subs x0,x0,#1
	b.ne 1b
	ret
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
	.macro H2L
	ldp x6,x7,[x1],#16
	ldp x8,x9,[x1],#16
	and x4,x9,#0x7ffffffff
	bfi x9,x17,#0,#35
	st64b x4,[x0]
	.endm
	.macro H2P
	ld64b x6,[x1]
	add x1,x1,#64
	and x4,x9,#0x7ffffffff
	bfi x9,x17,#0,#35
	st64b x4,[x0]
	and x8,x13,#0x7ffffffff
	bfi x13,x17,#0,#35
	st64b x8,[x0]
	.endm
	.macro H2P2
	H2P
	H2P
	.endm
	.macro H2P4
	H2P2
	H2P2
	.endm
	.macro H2P8
	H2P4
	H2P4
	.endm
	.macro H2P16
	H2P8
	H2P8
	.endm
	.globl h2s
h2s:
1:
	ldr x2,[x19],#32
	mov x1,x20
	add x20,x20,x21
	cmp x1,x2
	b.hs 6f
	H2L
	sub x16,x2,x1
	and x22,x2,#0xffffffffffffffc0
	tbz x16,#6,2f
	H2P
2:
	tbz x16,#7,3f
	H2P2
3:
	tbz x16,#8,4f
	H2P4
4:
	tbz x16,#9,5f
	H2P8
5:
	tbz x16,#10,7f
	H2P16
7:
	cmp x1,x22
	b.hs 8f
9:
	H2P16
	H2P16
	cmp x1,x22
	b.lo 9b
8:
	tbz x16,#5,6f
	H2L
6:
	subs x3,x3,#1
	b.ne 1b
	ret
	.macro ZVA1
	dc zva, x1
	add x1, x1, #256
	.endm
	.macro ZVA8
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	ZVA1
	.endm
	.macro ZVA64
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	ZVA8
	.endm
	.macro ZVA512
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	ZVA64
	.endm
	.globl zero16384
zero16384:
1:
	ZVA512
	subs x0,x0,#16384
	b.ne 1b
	ret
.globl receive
receive:
 ld64b x18,[x0]
 tbnz x18,#63,done_flag
 cmp x10,x11
 b.hs full
 stp x22,x23,[x10,#40]!
 bfi x12,x20,#3,#60
 swp x10,x5,[x12]
 stp x20,x5,[x10,#-16]
 str x21,[x10,#16]
 cbz x5,rdone
wa:
 ldp x7,x6,[x5,#-16]
 cmp x7,x20
 b.eq ha
ca:
 cbz x6,rdone
 ldp x7,x5,[x6,#-16]
 cmp x7,x20
 b.eq hb
cb:
 cbnz x5,wa
rdone:
 ret
ha:
 ldr x8,[x5,#16]
 eor x8,x8,x21
 tst x8,#0xffffffffffff
 ccmp x3,x4,#2,eq
 b.hs ca
 stp x5,x10,[x3],#16
 ret
hb:
 ldr x8,[x6,#16]
 eor x8,x8,x21
 tst x8,#0xffffffffffff
 ccmp x3,x4,#2,eq
 b.hs cb
 stp x6,x10,[x3],#16
 ret
full:
 ret
done_flag:
 add x15,x15,#1
 cmp x15,x16
 b.eq final_done
 ret
final_done:
 str x3,[x2]
 stlr x17,[x13]
 mov x15,#0
 ret
.globl h2done
h2done:
 mov x4,#0x8000000000000000
1:
 st64b x4,[x0]
 add x4,x4,#1
 subs x2,x2,#1
 b.ne 1b
 ret

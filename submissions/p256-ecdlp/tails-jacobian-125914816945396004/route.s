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
	.macro H2R
	ldp x5,x6,[x1],#16
	ldp x7,x8,[x1],#16
	and x9,x8,#0x7ffffffff
	bfi x8,x17,#0,#35
	stp x9,xzr,[x0,#-64]
	stp x5,x6,[x0,#-48]
	stp x7,x8,[x0,#-32]
	.endm
	.globl h2s
h2s:
	add x0,x0,#64
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
	.macro ZVA1
	dc zva, x1
	add x1, x1, #64
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
	.globl zero4096
zero4096:
1:
	ZVA512
	subs x0,x0,#4096
	b.ne 1b
	ret
.globl receive
receive:
 ldp x5,x6,[x0,#16]
 ldp x7,x8,[x0,#32]
 ldr x12,[x0,#8]
 cbnz x12,done_flag
 cmp x10,x11
 b.hs full
 stp x5,x6,[x10,#32]!
 stp x7,x8,[x10,#16]
 and x12,x5,x1
 ldr x20,[x9,x12,lsl #3]
 str x10,[x9,x12,lsl #3]
 lsr x21,x10,#2
 str x20,[x21,x14]
 cbz x20,rdone
walk:
 ldr x22,[x20]
 cmp x22,x5
 b.eq k2hit
adv:
 lsr x21,x20,#2
 ldr x20,[x21,x14]
 cbnz x20,walk
rdone:
 ret
k2hit:
 ldr x23,[x20,#8]
 eor x23,x23,x6
 tst x23,#0xffffffffffff
 b.ne adv
 cmp x3,x4
 b.hs rdone
 stp x20,x10,[x3],#16
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
 add x0,x0,#64
 mov x1,#0
1:
 stp x1,x3,[x0,#-64]
 stp xzr,xzr,[x0,#-48]
 stp xzr,xzr,[x0,#-32]
 add x1,x1,#1
 subs x2,x2,#1
 b.ne 1b
 ret

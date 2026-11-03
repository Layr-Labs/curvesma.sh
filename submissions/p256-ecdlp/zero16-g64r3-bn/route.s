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
	.macro HCHG
	ldp x7,x8,[x1,#32]!
	and x12,x7,x5
	add x10,x9,x12,lsl #3
	ldr x11,[x10]
	mov x22,x11
	cbz x11,8f
9:
	lsr x10,x11,#2
	add x10,x10,x14
	ldr x21,[x10]
	ldr x20,[x11]
	cmp x20,x7
	b.ne 7f
	ldr x17,[x11,#8]
	ubfx x17,x17,#0,#48
	ubfx x18,x8,#0,#48
	cmp x17,x18
	b.ne 7f
	cmp x3,x4
	b.hs 8f
	stp x11,x1,[x3],#16
	b 8f
7:
	mov x11,x21
	cbz x11,8f
	b 9b
8:
	add x10,x9,x12,lsl #3
	lsr x12,x1,#2
	add x12,x12,x14
	str x22,[x12]
	str x1,[x10]
	.endm
	.globl hchg
hchg:
1:
	ldr x2,[x13],#16
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
	.globl rcinit
rcinit:
1:
	str x2,[x1],#8
	add x2,x2,x3
	subs x0,x0,#1
	b.ne 1b
	ret
	.macro Z8
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
	stp xzr,xzr,[x1],#16
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
	.globl zero4096
zero4096:
1:
	Z512
	Z512
	Z512
	Z512
	subs x0,x0,#4096
	b.ne 1b
	ret
	.globl init3
init3:
1:
	add x4,x2,x3
	stp x2,x4,[x1],#16
	mov x2,x4
	subs x0,x0,#1
	b.ne 1b
	ret
.globl receive
receive:
 ldp x5,x6,[x0,#16]
 ldp x7,x8,[x0,#32]
 ldr x14,[x0,#8]
 cbnz x14,done_flag
 and x9,x5,x1
 add x9,x2,x9,lsl #4
 ldadd x18,x10,[x9]
 ldr x11,[x9,#8]
 cmp x10,x11
 b.hs full
 stp x5,x6,[x10]
 stp x7,x8,[x10,#16]
 ret
full:
 stadd x19,[x9]
 ret
done_flag:
 add x15,x15,#1
 cmp x15,x16
 b.eq final_done
 ret
final_done:
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

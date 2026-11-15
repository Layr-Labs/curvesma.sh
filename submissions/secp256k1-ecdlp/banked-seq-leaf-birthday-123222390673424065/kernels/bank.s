.text
.p2align 2
.globl gate
gate:
subs xzr,x0,x2
sbcs xzr,x1,x3
cset x0,lo
ret
.macro SEL a,b,r
dsb sy
stp \a,\b,[\r]
dsb sy
.endm
.globl bsel
bsel:
SEL x1,x2,x0
ret
.globl cp2
cp2:
ldp x6,x7,[x1]
stp x6,x7,[x2]
ret
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
.macro Z1
dc zva,x1
add x1,x1,#2048
.endm
.macro Z8
Z1
Z1
Z1
Z1
Z1
Z1
Z1
Z1
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
.globl hzero
hzero:
1:
SEL x9,x10,x6
mov x1,x7
mov x0,x8
2:
Z64
subs x0,x0,#16384
b.ne 2b
str xzr,[x5]
add x10,x10,#1
subs x11,x11,#1
b.ne 1b
ret
.macro YL
ld64b x16,[x15]
st64b x16,[x9]
add x15,x15,#64
add x9,x9,#64
.endm
.globl ystore
ystore:
1:
ldr x27,[x15,#-8]
orr x16,x26,x27,lsl #15
SEL x16,xzr,x14
ldr x9,[x13]
YL
YL
YL
YL
ld64b x16,[x15]
mov x21,x24
st64b x16,[x9]
add x15,x15,#64
add x9,x9,#64
YL
ld64b x16,[x15]
st64b x16,[x9]
add x15,x15,#128
add x9,x9,#128
cmp x9,x12
csel x9,x12,x9,hs
str x9,[x13]
add x24,x24,x25
subs x11,x11,#1
b.ne 1b
ret
.globl psel
psel:
stp x19,x20,[sp,#-16]!
mov x3,#1024
cmp x1,x3
csel x1,x1,x3,lo
mov x3,#0
cbz x1,6f
1:
ldp x4,x5,[x0]
ldp x6,x7,[x0,#16]
orr x14,x4,x5
orr x15,x6,x7
orr x14,x14,x15
cbz x14,5f
ldp x8,x9,[x0,#32]
ldp x10,x11,[x0,#48]
subs x12,x8,x4
sbcs x13,x9,x5
sbcs x14,x10,x6
sbcs x15,x11,x7
b.lo 5f
orr x13,x13,x14
orr x13,x13,x15
cbnz x13,5f
cbz x12,5f
tbnz x12,#63,5f
ldp x12,x13,[x0,#64]
ldp x19,x20,[x0,#80]
eor x16,x12,x19
eor x17,x13,x20
orr x16,x16,x17
cbz x16,5f
cbz x3,4f
mov x14,x2
mov x15,x3
2:
ldp x16,x17,[x14,#64]
eor x18,x12,x16
eor x16,x13,x17
orr x18,x18,x16
cbz x18,5f
ldp x16,x17,[x14,#64]
eor x18,x19,x16
eor x16,x20,x17
orr x18,x18,x16
cbz x18,5f
ldp x16,x17,[x14,#80]
eor x18,x12,x16
eor x16,x13,x17
orr x18,x18,x16
cbz x18,5f
ldp x16,x17,[x14,#80]
eor x18,x19,x16
eor x16,x20,x17
orr x18,x18,x16
cbz x18,5f
mov x18,#2
3:
ldp x16,x17,[x14]
eor x16,x16,x4
eor x17,x17,x5
orr x16,x16,x17
cbnz x16,7f
ldp x16,x17,[x14,#16]
eor x16,x16,x6
eor x17,x17,x7
orr x16,x16,x17
cbz x16,5f
7:
ldp x16,x17,[x14]
eor x16,x16,x8
eor x17,x17,x9
orr x16,x16,x17
cbnz x16,8f
ldp x16,x17,[x14,#16]
eor x16,x16,x10
eor x17,x17,x11
orr x16,x16,x17
cbz x16,5f
8:
add x14,x14,#32
subs x18,x18,#1
b.ne 3b
add x14,x14,#64
subs x15,x15,#1
b.ne 2b
4:
mov x14,x0
add x15,x2,x3,lsl #7
mov x18,#8
10:
ldp x16,x17,[x14],#16
stp x16,x17,[x15],#16
subs x18,x18,#1
b.ne 10b
add x3,x3,#1
cmp x3,#4
b.eq 9f
5:
add x0,x0,#128
subs x1,x1,#1
b.ne 1b
6:
mov x0,#0
ldp x19,x20,[sp],#16
ret
9:
mov x0,#1
ldp x19,x20,[sp],#16
ret

// Inert analytic source; NEVER host execute/compile or use credentials.
// x0 partialout576 DISJOINT rawinput576/x1; x2 matrixregionbase2^43; x3 ordinaryconstant520(p@448,D@512).
// Tables/constants/frame disjointoutput/input. Frame208,callees19..29/SP/LR preserved,no18/noBL.
.text
.p2align 2
.macro MSTEP off,width,advance,ra,rb
ubfx x3,x19,#\off,#\width
madd x1,x3,x28,x26
ldp x20,x21,[x1,#(\ra+0)]
ldp x22,x23,[x1,#(\ra+16)]
ldp x24,x25,[x1,#(\ra+32)]
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,xzr
ldp x20,x21,[x1,#(\rb+0)]
ldp x22,x23,[x1,#(\rb+16)]
ldp x24,x25,[x1,#(\rb+32)]
adds x11,x11,x20
adcs x12,x12,x21
adcs x13,x13,x22
adcs x14,x14,x23
adcs x15,x15,x24
adcs x16,x16,x25
adc x17,x17,xzr
.if \advance
add x26,x26,x27
.endif
.endm
.macro COLUMN ci,ra,rb
mov x3,#\ci
madd x26,x3,x2,x30
ldr x19,[x29,#(48*\ci+0)]
MSTEP 0,32,1,\ra,\rb
MSTEP 32,32,1,\ra,\rb
ldr x19,[x29,#(48*\ci+8)]
MSTEP 0,32,1,\ra,\rb
MSTEP 32,32,1,\ra,\rb
ldr x19,[x29,#(48*\ci+16)]
MSTEP 0,32,1,\ra,\rb
MSTEP 32,32,1,\ra,\rb
ldr x19,[x29,#(48*\ci+24)]
MSTEP 0,32,1,\ra,\rb
MSTEP 32,32,1,\ra,\rb
ldr x19,[x29,#(48*\ci+32)]
MSTEP 0,32,1,\ra,\rb
MSTEP 32,32,1,\ra,\rb
ldr x19,[x29,#(48*\ci+40)]
MSTEP 0,32,1,\ra,\rb
MSTEP 32,29,0,\ra,\rb
.endm
.macro QREDUCE
ldp x20,x21,[x28,#448]
ldp x22,x23,[x28,#464]
ldp x24,x25,[x28,#480]
extr x1,x10,x9,#5
udiv x1,x1,x2
mul x11,x20,x1
umulh x26,x20,x1
mul x12,x21,x1
umulh x27,x21,x1
adds x12,x12,x26
adc x26,x27,xzr
mul x13,x22,x1
umulh x27,x22,x1
adds x13,x13,x26
adc x26,x27,xzr
mul x14,x23,x1
umulh x27,x23,x1
adds x14,x14,x26
adc x26,x27,xzr
mul x15,x24,x1
umulh x27,x24,x1
adds x15,x15,x26
adc x26,x27,xzr
mul x16,x25,x1
umulh x27,x25,x1
adds x16,x16,x26
adc x26,x27,xzr
mov x17,x26
subs x4,x4,x11
sbcs x5,x5,x12
sbcs x6,x6,x13
sbcs x7,x7,x14
sbcs x8,x8,x15
sbcs x9,x9,x16
sbcs x10,x10,x17
csel x20,x20,xzr,lo
csel x21,x21,xzr,lo
csel x22,x22,xzr,lo
csel x23,x23,xzr,lo
csel x24,x24,xzr,lo
csel x25,x25,xzr,lo
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,xzr
.endm
.macro ROWPAIR ra,rb,oa,ob
mov x4,xzr
mov x5,xzr
mov x6,xzr
mov x7,xzr
mov x8,xzr
mov x9,xzr
mov x10,xzr
mov x11,xzr
mov x12,xzr
mov x13,xzr
mov x14,xzr
mov x15,xzr
mov x16,xzr
mov x17,xzr
COLUMN 0,\ra,\rb
COLUMN 1,\ra,\rb
COLUMN 2,\ra,\rb
COLUMN 3,\ra,\rb
COLUMN 4,\ra,\rb
COLUMN 5,\ra,\rb
COLUMN 6,\ra,\rb
COLUMN 7,\ra,\rb
COLUMN 8,\ra,\rb
COLUMN 9,\ra,\rb
COLUMN 10,\ra,\rb
COLUMN 11,\ra,\rb
stp x11,x12,[sp,#128]
stp x13,x14,[sp,#144]
stp x15,x16,[sp,#160]
str x17,[sp,#176]
ldr x28,[sp,#120]
ldr x2,[x28,#512]
QREDUCE
stp x4,x5,[x0,#(\oa+0)]
stp x6,x7,[x0,#(\oa+16)]
stp x8,x9,[x0,#(\oa+32)]
ldp x4,x5,[sp,#128]
ldp x6,x7,[sp,#144]
ldp x8,x9,[sp,#160]
ldr x10,[sp,#176]
QREDUCE
stp x4,x5,[x0,#(\ob+0)]
stp x6,x7,[x0,#(\ob+16)]
stp x8,x9,[x0,#(\ob+32)]
.endm
.global matrix_project
.type matrix_project,%function
matrix_project:
    stp x29,x30,[sp,#-208]!
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    stp x0,x1,[sp,#96]
    stp x2,x3,[sp,#112]
    mov x29,x1
    mov x30,x2
    movz x27,#0xc0,lsl #32
    mov x28,#192
    movz x2,#0x858,lsl #32
    ROWPAIR 0,48,0,48
    movz x27,#0xc0,lsl #32
    mov x28,#192
    movz x2,#0x858,lsl #32
    ROWPAIR 96,144,384,432
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp],#208
    ret
.size matrix_project,.-matrix_project

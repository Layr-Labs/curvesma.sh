// Inert source evidence; shared credential-free verifier only; NEVER host execute.
// x0 out480,x1 input480 identical/disjoint,x2 forty-TiB regionbase,x3 powers192.
.text
.p2align 2
.macro DIGIT reg,off,width,advance
ubfx x1,\reg,#\off,#\width
madd x2,x1,x28,x26
ldp x12,x13,[x2,#0]
ldp x14,x15,[x2,#16]
ldp x16,x17,[x2,#32]
adds x5,x5,x12
adcs x6,x6,x13
adcs x7,x7,x14
adcs x8,x8,x15
adcs x9,x9,x16
adcs x10,x10,x17
adc x11,x11,xzr
.if \advance
add x26,x26,x27
.endif
.endm
.macro REDUCE off
ldp x12,x13,[x25,#(\off+0)]
ldp x14,x15,[x25,#(\off+16)]
ldp x16,x17,[x25,#(\off+32)]
subs x19,x5,x12
sbcs x20,x6,x13
sbcs x21,x7,x14
sbcs x22,x8,x15
sbcs x23,x9,x16
sbcs x24,x10,x17
sbcs x0,x11,xzr
csel x5,x19,x5,hs
csel x6,x20,x6,hs
csel x7,x21,x7,hs
csel x8,x22,x8,hs
csel x9,x23,x9,hs
csel x10,x24,x10,hs
csel x11,x0,x11,hs
.endm
.macro FIELD bo,ci
mov x26,#\ci
add x26,x3,x26,lsl #42
ldp x19,x20,[x30,#(\bo+0)]
ldp x21,x22,[x30,#(\bo+16)]
ldp x23,x24,[x30,#(\bo+32)]
mov x5,xzr
mov x6,xzr
mov x7,xzr
mov x8,xzr
mov x9,xzr
mov x10,xzr
mov x11,xzr
DIGIT x19,0,32,1
DIGIT x19,32,32,1
DIGIT x20,0,32,1
DIGIT x20,32,32,1
DIGIT x21,0,32,1
DIGIT x21,32,32,1
DIGIT x22,0,32,1
DIGIT x22,32,32,1
DIGIT x23,0,32,1
DIGIT x23,32,32,1
DIGIT x24,0,32,1
DIGIT x24,32,29,0
REDUCE 0
REDUCE 48
REDUCE 96
REDUCE 144
stp x5,x6,[x4,#(\bo+0)]
stp x7,x8,[x4,#(\bo+16)]
stp x9,x10,[x4,#(\bo+32)]
.endm
.global lookup_ten32
.type lookup_ten32,%function
lookup_ten32:
    stp x29,x30,[sp,#-96]!
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    mov x29,sp
    mov x25,x3
    mov x30,x1
    mov x3,x2
    mov x4,x0
    mov x27,#3
    lsl x27,x27,#36
    mov x28,#48
    FIELD 0,0
    FIELD 48,1
    FIELD 96,2
    FIELD 144,3
    FIELD 192,4
    FIELD 240,5
    FIELD 288,6
    FIELD 336,7
    FIELD 384,8
    FIELD 432,9
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp],#96
    ret
.size lookup_ten32,.-lookup_ten32

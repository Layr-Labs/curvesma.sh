// Inert evidence;shared credential-free verifier ONLY. Never host-execute/compile.
// x0 third192 EXACT alias,x1 transformed384,x2 tablebase8TiB,x3 constant3728.
// Two vector coefficient tables A/B,11 digits7x35+4x34,27TiB each.
// Input second Fp2 has already been multiplied byxi. Imaginary input swaps rows,
// subtracts second table component fromreal,adds first toimag. Real seed+22p.
// Every intermediate real/imag lies[0,45p);same50 qp slots and QREDUCE valid.
// Frame208/no18/preserve19..30/SP/LR;writes only current third pair.
.text
.p2align 2
.macro MSTEP ra,rb,negative,width,advance
madd x1,x3,x28,x26
ldp x20,x21,[x1,#(\ra+0)]
ldp x22,x23,[x1,#(\ra+16)]
ldp x24,x25,[x1,#(\ra+32)]
.if \negative
subs x4,x4,x20
sbcs x5,x5,x21
sbcs x6,x6,x22
sbcs x7,x7,x23
sbcs x8,x8,x24
sbcs x9,x9,x25
sbc x10,x10,xzr
.endif
.if (1-\negative)
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,xzr
.endif
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
.if (\width-34)
add x26,x26,x27
.endif
.if (35-\width)
add x26,x26,x27,lsr #1
.endif
.endif
.endm
.macro COLUMN ci,tableid,imag
mov x3,#\tableid
madd x26,x3,x2,x30
ldr x19,[x29,#(48*\ci+0)]
ubfx x3,x19,#0,#35
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ldr x20,[x29,#(48*\ci+8)]
extr x3,x20,x19,#35
and x3,x3,#0x7ffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ubfx x3,x19,#6,#35
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ldr x20,[x29,#(48*\ci+16)]
extr x3,x20,x19,#41
and x3,x3,#0x7ffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ubfx x3,x19,#12,#35
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ldr x20,[x29,#(48*\ci+24)]
extr x3,x20,x19,#47
and x3,x3,#0x7ffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ubfx x3,x19,#18,#35
MSTEP (48*\imag),(48*(1-\imag)),\imag,35,1
ldr x20,[x29,#(48*\ci+32)]
extr x3,x20,x19,#53
and x3,x3,#0x3ffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,34,1
ubfx x3,x19,#23,#34
MSTEP (48*\imag),(48*(1-\imag)),\imag,34,1
ldr x20,[x29,#(48*\ci+40)]
extr x3,x20,x19,#57
and x3,x3,#0x3ffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,34,1
ubfx x3,x19,#27,#34
MSTEP (48*\imag),(48*(1-\imag)),\imag,34,0
.endm
.macro QREDUCE
ldp x20,x21,[x28,#448]
ldp x22,x23,[x28,#464]
ldp x24,x25,[x28,#480]
extr x1,x10,x9,#5
udiv x1,x1,x2
add x26,x28,#528
add x26,x26,x1,lsl #6
ldp x11,x12,[x26,#0]
ldp x13,x14,[x26,#16]
ldp x15,x16,[x26,#32]
ldr x17,[x26,#48]
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
.macro ROWPAIR sa,sb,c0,c1,c2,c3,oa,ob
ldr x1,[sp,#96]
ldp x4,x5,[x1,#(\sa+0)]
ldp x6,x7,[x1,#(\sa+16)]
ldp x8,x9,[x1,#(\sa+32)]
mov x10,xzr
ldp x11,x12,[x1,#(\sb+0)]
ldp x13,x14,[x1,#(\sb+16)]
ldp x15,x16,[x1,#(\sb+32)]
mov x17,xzr
ldr x26,[sp,#120]
add x26,x26,#1936
ldp x20,x21,[x26]
ldp x22,x23,[x26,#16]
ldp x24,x25,[x26,#32]
ldr x3,[x26,#48]
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,x3
COLUMN \c0,0,0
COLUMN \c1,0,1
COLUMN \c2,1,0
COLUMN \c3,1,1
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
.global recurrence_fp2_map
.type recurrence_fp2_map,%function
recurrence_fp2_map:
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
    movz x27,#0x300,lsl #32
    mov x28,#96
    movz x2,#0x1b00,lsl #32
    ROWPAIR 0,48,0,1,2,3,0,48
    movz x27,#0x300,lsl #32
    mov x28,#96
    movz x2,#0x1b00,lsl #32
    ROWPAIR 96,144,4,5,6,7,96,144
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp],#208
    ret
.size recurrence_fp2_map,.-recurrence_fp2_map

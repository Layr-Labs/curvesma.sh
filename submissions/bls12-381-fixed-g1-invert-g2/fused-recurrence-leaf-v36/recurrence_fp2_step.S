// Inert analytic source;shared credential-free llm-verifier ONLY.
// Full fused leaf:out192,prev1/prev2 canonical192,tablebase,const3728,p48 compatibility.
// Frame576:callee96,args32,workspace64,temp384@192;originalSP-384 tempbase.
// BOTH uses x29 tempbase;map phase preserves x0 out;no BL/x18,all callee/SP/LR restored.
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
.macro COLUMNBASE0
mov x26,x30
.endm
.macro COLUMNBASE1
add x26,x30,x2
.endm
.macro COLUMNBASE2
add x26,x30,x2,lsl #1
.endm
.macro COLUMNBASE3
add x26,x26,x27,lsr #2
.endm
.macro COLUMN ci,ti
COLUMNBASE\ti
ldr x19,[x29,#(48*\ci+0)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+8)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+16)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+24)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+32)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+40)]
MSTEP 0,32,1,0,48
MSTEP 32,30,0,0,48
.endm
.macro WSTEP advance
madd x1,x3,x28,x26
ldp x20,x21,[x1,#(0+0)]
ldp x22,x23,[x1,#(0+16)]
ldp x24,x25,[x1,#(0+32)]
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,xzr
ldp x20,x21,[x1,#(48+0)]
ldp x22,x23,[x1,#(48+16)]
ldp x24,x25,[x1,#(48+32)]
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
.macro COLUMNWIDE ci,second
    .if 1-\second
    ldr x26,[sp,#104]
    .endif
    .if \second
    add x26,x26,x27,lsl #1
    movz x27,#0x300,lsl #32
    .endif
    ldp x19,x2,[x29,#(48*\ci+0)]
    ubfx x3,x19,#0,#35
    WSTEP 1
    extr x3,x2,x19,#35
    and x3,x3,#0x7ffffffff
    WSTEP 1
    ubfx x3,x2,#6,#35
    WSTEP 1
    ldr x19,[x29,#(48*\ci+16)]
    extr x3,x19,x2,#41
    and x3,x3,#0x7ffffffff
    WSTEP 1
    ubfx x3,x19,#12,#35
    WSTEP 1
    ldr x2,[x29,#(48*\ci+24)]
    extr x3,x2,x19,#47
    and x3,x3,#0x7ffffffff
    WSTEP 1
    ubfx x3,x2,#18,#35
    WSTEP 1
    movz x27,#0x180,lsl #32
    ldr x19,[x29,#(48*\ci+32)]
    extr x3,x19,x2,#53
    and x3,x3,#0x3ffffffff
    WSTEP 1
    ubfx x3,x19,#23,#34
    WSTEP 1
    ldr x2,[x29,#(48*\ci+40)]
    extr x3,x2,x19,#57
    and x3,x3,#0x3ffffffff
    WSTEP 1
    ubfx x3,x2,#27,#35
    WSTEP 0
    .if \second
    add x30,x26,x27,lsl #1
    movz x27,#0x60,lsl #32
    movz x2,#0x438,lsl #32
    .endif
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
.macro ROWPAIR sa,sb,c0,c1,c2,c3,oa,ob,p2,p3
ldp x4,x5,[x0,#(\sa+0)]
ldp x6,x7,[x0,#(\sa+16)]
ldp x8,x9,[x0,#(\sa+32)]
mov x10,xzr
ldp x11,x12,[x0,#(\sb+0)]
ldp x13,x14,[x0,#(\sb+16)]
ldp x15,x16,[x0,#(\sb+32)]
mov x17,xzr
COLUMNWIDE \c0,0
COLUMNWIDE \c1,1
COLUMN \c2,\p2
COLUMN \c3,\p3
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
.macro BOTH so,di
ldp x4,x5,[x1,#0]
ldp x6,x7,[x1,#16]
ldp x8,x9,[x1,#32]
ldp x10,x11,[x2,#0]
ldp x12,x13,[x2,#16]
ldp x14,x15,[x2,#32]
adds x0,x4,x10
adcs x3,x5,x11
adcs x17,x6,x12
adcs x16,x7,x13
adcs x26,x8,x14
adc x27,x9,x15
subs x4,x4,x10
sbcs x5,x5,x11
sbcs x6,x6,x12
sbcs x7,x7,x13
sbcs x8,x8,x14
sbcs x9,x9,x15
stp x4,x5,[x29,#(\di+0)]
stp x6,x7,[x29,#(\di+16)]
stp x8,x9,[x29,#(\di+32)]
stp x0,x3,[x29,#(\so+0)]
stp x17,x16,[x29,#(\so+16)]
stp x26,x27,[x29,#(\so+32)]
.endm
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    sub sp,sp,#576
    stp x29,x30,[sp]
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    stp x0,x3,[sp,#96]
    stp xzr,x4,[sp,#112]
    add x29,sp,#192
    BOTH 288,0
    add x1,x1,#48
    add x2,x2,#48
    BOTH 336,48
    add x1,x1,#48
    add x2,x2,#48
    BOTH 96,192
    add x1,x1,#48
    add x2,x2,#48
    BOTH 144,240
    ldr x0,[sp,#96]
    movz x27,#0x300,lsl #32
    mov x28,#96
    ROWPAIR 0,48,0,1,2,3,0,48,0,1
    movz x27,#0x300,lsl #32
    mov x28,#96
    ROWPAIR 96,144,4,5,6,7,96,144,2,3
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp]
    add sp,sp,#576
    ret
.size recurrence_fp2_step,.-recurrence_fp2_step

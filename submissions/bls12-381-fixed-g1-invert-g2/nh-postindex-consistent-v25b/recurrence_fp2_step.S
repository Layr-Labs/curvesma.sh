// Inert analytic source; NEVER host execute/compile or use credentials.
// x0 third_out192,x1 prev1_192,x2 prev2_192,x3 tablebase,x4 const3728,x5 p48.
// Distinct raw histories, temporary384@SP144, frame528 incl arguments@96..143.
// BOTH clobbers x4..17,x19..30,including29/LR; ALL wrapper frame accesses use SP.
// Wrapper restores19..29/SP/LR,no18. Full actual field outputs canonical, no BL in BOTH.
.text
.p2align 2
.macro BOTH
ldp x4,x5,[x1,#0]
ldp x6,x7,[x1,#16]
ldp x8,x9,[x1,#32]
ldp x10,x11,[x2,#0]
ldp x12,x13,[x2,#16]
ldp x14,x15,[x2,#32]
adds x25,x4,x10
adcs x26,x5,x11
adcs x27,x6,x12
adcs x28,x7,x13
adcs x29,x8,x14
adc x30,x9,x15
subs x4,x4,x10
sbcs x5,x5,x11
sbcs x6,x6,x12
sbcs x7,x7,x13
sbcs x8,x8,x14
sbcs x9,x9,x15
csel x10,x19,xzr,lo
csel x11,x20,xzr,lo
csel x12,x21,xzr,lo
csel x13,x22,xzr,lo
csel x14,x23,xzr,lo
csel x15,x24,xzr,lo
adds x4,x4,x10
adcs x5,x5,x11
adcs x6,x6,x12
adcs x7,x7,x13
adcs x8,x8,x14
adc x9,x9,x15
stp x4,x5,[x3,#0]
stp x6,x7,[x3,#16]
stp x8,x9,[x3,#32]
subs x4,x25,x19
sbcs x5,x26,x20
sbcs x6,x27,x21
sbcs x7,x28,x22
sbcs x8,x29,x23
sbcs x9,x30,x24
csel x25,x4,x25,hs
csel x26,x5,x26,hs
csel x27,x6,x27,hs
csel x28,x7,x28,hs
csel x29,x8,x29,hs
csel x30,x9,x30,hs
stp x25,x26,[x0,#0]
stp x27,x28,[x0,#16]
stp x29,x30,[x0,#32]
.endm
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    sub sp,sp,#528
    stp x29,x30,[sp]
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    stp x0,x1,[sp,#96]
    stp x2,x3,[sp,#112]
    stp x4,x5,[sp,#128]
    add x0,sp,#432
    add x3,sp,#144
    ldr x16,[sp,#136]
    ldp x19,x20,[x16,#0]
    ldp x21,x22,[x16,#16]
    ldp x23,x24,[x16,#32]
    BOTH
    add x0,sp,#480
    add x3,sp,#192
    add x1,x1,#48
    add x2,x2,#48
    BOTH
    add x0,sp,#240
    add x3,sp,#336
    add x1,x1,#48
    add x2,x2,#48
    BOTH
    add x0,sp,#288
    add x3,sp,#384
    add x1,x1,#48
    add x2,x2,#48
    BOTH
    ldr x0,[sp,#96]
    add x1,sp,#144
    ldp x2,x3,[sp,#120]
    bl recurrence_fp2_map
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp]
    add sp,sp,#528
    ret
.size recurrence_fp2_step,.-recurrence_fp2_step

// Inert analytic source; NEVER host execute/compile or use credentials.
// x0 third_out192,x1 prev1_192,x2 prev2_192,x3 tablebase,x4 const3728,x5 unused p48 compatibility argument.
// Distinct raw histories, temporary384@SP80,frame464 incl arguments@48..79.
// BOTH clobbers x4..15,x25..30;constant x16 survives;19..24 untouched; ALL wrapper frame accesses use SP.
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
add x9,x9,x16
stp x4,x5,[x3,#0]
stp x6,x7,[x3,#16]
stp x8,x9,[x3,#32]
stp x25,x26,[x0,#0]
stp x27,x28,[x0,#16]
stp x29,x30,[x0,#32]
.endm
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    stp x29,x30,[sp,#-464]!
    stp x25,x26,[sp,#16]
    stp x27,x28,[sp,#32]
    stp x0,x3,[sp,#48]
    stp x4,xzr,[sp,#64]
    add x0,sp,#368
    add x3,sp,#80
    movz x16,#0x2000,lsl #48
    BOTH
    add x0,sp,#416
    add x3,sp,#128
    add x1,x1,#48
    add x2,x2,#48
    BOTH
    add x0,sp,#176
    add x3,sp,#272
    add x1,x1,#48
    add x2,x2,#48
    BOTH
    add x0,sp,#224
    add x3,sp,#320
    add x1,x1,#48
    add x2,x2,#48
    BOTH
    ldr x0,[sp,#48]
    add x1,sp,#80
    ldp x2,x3,[sp,#56]
    bl recurrence_fp2_map
    ldp x25,x26,[sp,#16]
    ldp x27,x28,[sp,#32]
    ldp x29,x30,[sp],#464
    ret
.size recurrence_fp2_step,.-recurrence_fp2_step

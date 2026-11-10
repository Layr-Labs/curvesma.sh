// Inert analytic source; NEVER host execute/compile or use credentials.
// x0 third_out192,x1 prev1_192,x2 prev2_192,x3 tablebase,x4 const3728,x5 unused p48 compatibility argument.
// Distinct raw histories, temporary384@SP64,frame448;callee32/arguments32.
// BOTH clobbers x0,x3..15,x17,x25..27;inputs x1/x2 andbias x16 survive; ALL wrapper frame accesses use SP.
// Wrapper restores25..27/SP/LR;19..24/28/29 untouched,no18. Full actual field outputs canonical, no BL in BOTH.
.text
.p2align 2
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
adcs x25,x7,x13
adcs x26,x8,x14
adc x27,x9,x15
subs x4,x4,x10
sbcs x5,x5,x11
sbcs x6,x6,x12
sbcs x7,x7,x13
sbcs x8,x8,x14
sbcs x9,x9,x15
add x9,x9,x16
stp x4,x5,[sp,#(\di+0)]
stp x6,x7,[sp,#(\di+16)]
stp x8,x9,[sp,#(\di+32)]
stp x0,x3,[sp,#(\so+0)]
stp x17,x25,[sp,#(\so+16)]
stp x26,x27,[sp,#(\so+32)]
.endm
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    stp x25,x26,[sp,#-448]!
    stp x27,x30,[sp,#16]
    stp x0,x3,[sp,#32]
    stp x4,xzr,[sp,#48]
    movz x16,#0x2000,lsl #48
    BOTH 352,64
    add x1,x1,#48
    add x2,x2,#48
    BOTH 400,112
    add x1,x1,#48
    add x2,x2,#48
    BOTH 160,256
    add x1,x1,#48
    add x2,x2,#48
    BOTH 208,304
    ldr x0,[sp,#32]
    add x1,sp,#64
    ldp x2,x3,[sp,#40]
    bl recurrence_fp2_map
    ldp x27,x30,[sp,#16]
    ldp x25,x26,[sp],#448
    ret
.size recurrence_fp2_step,.-recurrence_fp2_step

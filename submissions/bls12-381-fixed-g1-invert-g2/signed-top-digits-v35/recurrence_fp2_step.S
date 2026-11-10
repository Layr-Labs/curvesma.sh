// Inert analytic source; NEVER host execute/compile or use credentials.
// x0 third_out192,x1 prev1_192,x2 prev2_192,x3 tablebase,x4 const3728,x5 unused p48 compatibility argument.
// Distinct raw histories, temporary384@SP48,frame432;callee24/arguments24.
// BOTH clobbers x0,x3..17,x26/x27;inputs x1/x2 survive; ALL wrapper frame accesses use SP.
// Wrapper restores26/27/SP/LR;19..25/28/29 untouched,no18. Full actual field outputs canonical, no BL in BOTH.
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
adcs x16,x7,x13
adcs x26,x8,x14
adc x27,x9,x15
subs x4,x4,x10
sbcs x5,x5,x11
sbcs x6,x6,x12
sbcs x7,x7,x13
sbcs x8,x8,x14
sbcs x9,x9,x15
stp x4,x5,[sp,#(\di+0)]
stp x6,x7,[sp,#(\di+16)]
stp x8,x9,[sp,#(\di+32)]
stp x0,x3,[sp,#(\so+0)]
stp x17,x16,[sp,#(\so+16)]
stp x26,x27,[sp,#(\so+32)]
.endm
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    stp x26,x27,[sp,#-432]!
    stp x30,x0,[sp,#16]
    stp x3,x4,[sp,#32]
    BOTH 336,48
    add x1,x1,#48
    add x2,x2,#48
    BOTH 384,96
    add x1,x1,#48
    add x2,x2,#48
    BOTH 144,240
    add x1,x1,#48
    add x2,x2,#48
    BOTH 192,288
    ldr x0,[sp,#24]
    add x1,sp,#48
    ldp x2,x3,[sp,#32]
    bl recurrence_fp2_map
    ldr x30,[sp,#16]
    ldp x26,x27,[sp],#432
    ret
.size recurrence_fp2_step,.-recurrence_fp2_step

// Static evidence only;neverhost-executed.
// Candidatearena12*96;realfieldsinitializedbeforestage1.
// WMETA24bytes=(evenk,oddk,mode0normal/1evenzero/2oddzero).
.text
.p2align 2
.macro COMPARE_COMPONENT k,off
    add x1,x8,#(96*\k)
    mov x7,#\k
    ldp x2,x3,[x0,#(\off+0)]
    ldp x4,x5,[x1,#(\off+0)]
    subs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#(\off+16)]
    ldp x4,x5,[x1,#(\off+16)]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#(\off+32)]
    ldp x4,x5,[x1,#(\off+32)]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    csel x0,x1,x0,hs
    csel x6,x7,x6,hs
.endm
.macro ZERO_REAL ptr,target
    ldp x2,x3,[\ptr,#0]
    orr x7,x2,x3
    ldp x2,x3,[\ptr,#16]
    orr x7,x7,x2
    orr x7,x7,x3
    ldp x2,x3,[\ptr,#32]
    orr x7,x7,x2
    orr x7,x7,x3
    cbnz x7,\target
.endm
// x0candidatearena,x1WMETA. Noimaginaryfieldisread.
.global fapi_family_real_stage
.type fapi_family_real_stage,%function
fapi_family_real_stage:
    mov x8,x0
    mov x9,x1
    mov x6,xzr
    ZERO_REAL x0,.Llazy_first_nonzero
    stp xzr,xzr,[x9]
    mov x2,#1
    str x2,[x9,#16]
    ret
.Llazy_first_nonzero:
    add x1,x8,#96
    ZERO_REAL x1,.Llazy_both_nonzero
    stp xzr,xzr,[x9]
    mov x2,#2
    str x2,[x9,#16]
    ret
.Llazy_both_nonzero:
    COMPARE_COMPONENT 2,0
    COMPARE_COMPONENT 4,0
    COMPARE_COMPONENT 6,0
    COMPARE_COMPONENT 8,0
    COMPARE_COMPONENT 10,0
    mov x11,x6
    add x0,x8,#96
    mov x6,#1
    COMPARE_COMPONENT 3,0
    COMPARE_COMPONENT 5,0
    COMPARE_COMPONENT 7,0
    COMPARE_COMPONENT 9,0
    COMPARE_COMPONENT 11,0
    stp x11,x6,[x9]
    str xzr,[x9,#16]
    ret
.size fapi_family_real_stage,.-fapi_family_real_stage
.macro OUTPUT_SPLIT
    ldp x2,x3,[x0,#0]
    stp x2,x3,[x9,#0]
    ldp x2,x3,[x0,#16]
    stp x2,x3,[x9,#16]
    ldp x2,x3,[x0,#32]
    stp x2,x3,[x9,#32]
    ldp x2,x3,[x1,#0]
    stp x2,x3,[x9,#48]
    ldp x2,x3,[x1,#16]
    stp x2,x3,[x9,#64]
    ldp x2,x3,[x1,#32]
    stp x2,x3,[x9,#80]
    str x6,[x10]
    mov x0,x6
    ret
.endm
// x0arena,x1WMETA,x2WIN_Eimag48,x3WIN_Oimag48,x4output96,x5koutput8.
// NormalWIN_E/Oareinitialized;zero-modechosenfamilyall6imagsinitialized.
.global fapi_family_complete
.type fapi_family_complete,%function
fapi_family_complete:
    mov x8,x0
    mov x9,x4
    mov x10,x5
    mov x14,x2
    mov x15,x3
    ldp x12,x13,[x1]
    ldr x11,[x1,#16]
    cbnz x11,.Llazy_zero_complete
    add x16,x12,x12,lsl #1
    add x16,x8,x16,lsl #5
    add x17,x13,x13,lsl #1
    add x17,x8,x17,lsl #5
    ldp x2,x3,[x16,#0]
    ldp x4,x5,[x17,#0]
    subs x2,x2,x4
    sbcs x3,x3,x5
    orr x7,x2,x3
    ldp x2,x3,[x16,#16]
    ldp x4,x5,[x17,#16]
    sbcs x2,x2,x4
    sbcs x3,x3,x5
    orr x2,x2,x3
    orr x7,x7,x2
    ldp x2,x3,[x16,#32]
    ldp x4,x5,[x17,#32]
    sbcs x2,x2,x4
    sbcs x3,x3,x5
    orr x2,x2,x3
    orr x7,x7,x2
    cbnz x7,.Llazy_select
    ldp x2,x3,[x14,#0]
    ldp x4,x5,[x15,#0]
    subs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x14,#16]
    ldp x4,x5,[x15,#16]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x14,#32]
    ldp x4,x5,[x15,#32]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
.Llazy_select:
    csel x0,x17,x16,hs
    csel x1,x15,x14,hs
    csel x6,x13,x12,hs
    OUTPUT_SPLIT
.Llazy_zero_complete:
    cmp x11,#1
    b.eq .Llazy_even_zero
    add x0,x8,#96
    mov x6,#1
    COMPARE_COMPONENT 3,48
    COMPARE_COMPONENT 5,48
    COMPARE_COMPONENT 7,48
    COMPARE_COMPONENT 9,48
    COMPARE_COMPONENT 11,48
    b .Llazy_zero_output
.Llazy_even_zero:
    mov x0,x8
    mov x6,xzr
    COMPARE_COMPONENT 2,48
    COMPARE_COMPONENT 4,48
    COMPARE_COMPONENT 6,48
    COMPARE_COMPONENT 8,48
    COMPARE_COMPONENT 10,48
.Llazy_zero_output:
    add x1,x0,#48
    OUTPUT_SPLIT
.size fapi_family_complete,.-fapi_family_complete

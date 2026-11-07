// Static analytic evidence only;neverhost-executed. x0raw960,x1outS192.
// ReturnsSTATUS1andwritesrecordk944;identityreturns0. No scalar-A metadata read.
// Coldconstructorinitializesthewhole576-byteworkspaceandnativeSbelow.
// Canonical negation uses a dedicated bounded leaf; no zero operand buffer.
.text
.p2align 2
.global fapi_trace_canonical
.type fapi_trace_canonical,%function
fapi_trace_canonical:
    stp x19,x20,[sp,#-48]!
    stp x21,x22,[sp,#16]
    stp x23,x30,[sp,#32]
    mov x19,x0
    mov x20,x1
    adrp x21,fapi_trace_workspace
    add x21,x21,:lo12:fapi_trace_workspace
    add x23,x21,#384
    mov x22,x21
// Copy B directly into candidate zero and test every limb during the copy.
    ldp x0,x1,[x19,#384]
    stp x0,x1,[x22,#0]
    orr x9,x0,x1
    ldp x0,x1,[x19,#400]
    stp x0,x1,[x22,#16]
    orr x9,x9,x0
    orr x9,x9,x1
    ldp x0,x1,[x19,#416]
    stp x0,x1,[x22,#32]
    orr x9,x9,x0
    orr x9,x9,x1
    ldp x0,x1,[x19,#432]
    stp x0,x1,[x22,#48]
    orr x9,x9,x0
    orr x9,x9,x1
    ldp x0,x1,[x19,#448]
    stp x0,x1,[x22,#64]
    orr x9,x9,x0
    orr x9,x9,x1
    ldp x0,x1,[x19,#464]
    stp x0,x1,[x22,#80]
    orr x9,x9,x0
    orr x9,x9,x1
    cbz x9,.Lcanonical_identity
// Nine direct calls:five unary negations,one add,one sub,two products.
    add x0,x19,#48
    mov x1,x23
    bl fapi_fp_neg
    mov x0,x22
    add x1,x22,#48
    add x2,x21,#432
    bl fapi_fp_add
    mov x0,x22
    add x1,x22,#48
    add x2,x21,#480
    bl fapi_fp_sub
    add x0,x22,#96
    add x1,x21,#432
    add x2,x21,#528
    bl fapi_native_fp_mul
    add x0,x22,#144
    add x1,x21,#480
    add x2,x21,#528
    bl fapi_native_fp_mul
    mov x0,x22
    add x1,x22,#192
    bl fapi_fp_neg
    add x0,x22,#48
    add x1,x22,#240
    bl fapi_fp_neg
    add x0,x22,#96
    add x1,x22,#288
    bl fapi_fp_neg
    add x0,x22,#144
    add x1,x22,#336
    bl fapi_fp_neg
    mov x0,x19
    mov x1,x23
    mov x2,x22
    mov x3,x20
    add x4,x19,#944
    bl fapi_trace_min
    mov x0,#1
    b .Lcanonical_return
.Lcanonical_identity:
    mov x0,xzr
.Lcanonical_return:
    ldp x21,x22,[sp,#16]
    ldp x23,x30,[sp,#32]
    ldp x19,x20,[sp],#48
    ret
.size fapi_trace_canonical,.-fapi_trace_canonical

.section .data
.p2align 4
// Thesezerosandallcode/descriptorbytesareexplicitlyinitialized/paidcold.
fapi_trace_workspace:
fapi_trace_bvariants:
    .zero 384
fapi_trace_aneg:
    .zero 48
fapi_trace_sum:
    .zero 48
fapi_trace_diff:
    .zero 48
.global fapi_trace_s
fapi_trace_s:
// Cold computes exact native s from C=(1+u)^((p-1)/2)=(s,s).
    .zero 48

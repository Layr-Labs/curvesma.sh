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
// Always initialize ANEG and both even-parity B candidates.
    add x0,x19,#48
    mov x1,x23
    bl fapi_fp_neg
    mov x0,x22
    add x1,x22,#192
    bl fapi_fp_neg
    add x0,x22,#48
    add x1,x22,#240
    bl fapi_fp_neg
// Compare raw Aimag and its canonical negative; equality is zero and goes heavy.
    ldp x0,x1,[x19,#48]
    ldp x2,x3,[x23,#0]
    subs xzr,x0,x2
    sbcs xzr,x1,x3
    ldp x0,x1,[x19,#64]
    ldp x2,x3,[x23,#16]
    sbcs xzr,x0,x2
    sbcs xzr,x1,x3
    ldp x0,x1,[x19,#80]
    ldp x2,x3,[x23,#32]
    sbcs xzr,x0,x2
    sbcs xzr,x1,x3
    b.lo .Lcanonical_min
// Count heavy calls explicitly and abort before an unbudgeted heavy body.
    adrp x15,fapi_trace_heavy_control
    add x15,x15,:lo12:fapi_trace_heavy_control
    ldp x9,x10,[x15,#0]
    adds x9,x9,#1
    adc x10,x10,xzr
    stp x9,x10,[x15,#0]
    ldp x11,x12,[x15,#16]
    subs xzr,x9,x11
    sbcs xzr,x10,x12
    b.hi .Lcanonical_identity
// Heavy body:initialize K(B) and its negative; even-only minimum skips these.
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
    add x0,x22,#96
    add x1,x22,#288
    bl fapi_fp_neg
    add x0,x22,#144
    add x1,x22,#336
    bl fapi_fp_neg
.Lcanonical_min:
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

.p2align 4
.global fapi_trace_heavy_control
fapi_trace_heavy_control:
// Cold initializes count128=0 and cap128=floor((1/2+2^-20)*K)+1.
// After the cold initial-state canonical call,reset count to zero before main.
    .zero 32

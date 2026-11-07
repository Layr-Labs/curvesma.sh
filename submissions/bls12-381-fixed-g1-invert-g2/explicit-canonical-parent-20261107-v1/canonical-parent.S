// Static analytic evidence only;neverhost-executed. x0raw960,x1outS192.
// ReturnsSTATUS1andwritesrecordk944;identityreturns0. No scalar-A metadata read.
// ColdconstructorinitializesthewholeimageandnativeFrobeniuscoefficientSbelow.
.text
.p2align 2
.global fapi_trace_canonical
.type fapi_trace_canonical,%function
fapi_trace_canonical:
    stp x19,x20,[sp,#-64]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    str x30,[sp,#48]
    mov x19,x0
    mov x20,x1
    adrp x21,fapi_trace_workspace
    add x21,x21,:lo12:fapi_trace_workspace
    add x23,x21,#192
    add x22,x21,#240
    add x24,x21,#768
// Exactnativeprojectiona0andb1intofixed192-byteTRAW.
    ldp x0,x1,[x19,#0]
    stp x0,x1,[x21,#0]
    ldp x0,x1,[x19,#16]
    stp x0,x1,[x21,#16]
    ldp x0,x1,[x19,#32]
    stp x0,x1,[x21,#32]
    ldp x0,x1,[x19,#48]
    stp x0,x1,[x21,#48]
    ldp x0,x1,[x19,#64]
    stp x0,x1,[x21,#64]
    ldp x0,x1,[x19,#80]
    stp x0,x1,[x21,#80]
    ldp x0,x1,[x19,#384]
    stp x0,x1,[x21,#96]
    ldp x0,x1,[x19,#400]
    stp x0,x1,[x21,#112]
    ldp x0,x1,[x19,#416]
    stp x0,x1,[x21,#128]
    ldp x0,x1,[x19,#432]
    stp x0,x1,[x21,#144]
    ldp x0,x1,[x19,#448]
    stp x0,x1,[x21,#160]
    ldp x0,x1,[x19,#464]
    stp x0,x1,[x21,#176]
// B!=0iffnonidentityonGT;all12limbsaretested.
    ldp x0,x1,[x21,#96]
    orr x0,x0,x1
    ldp x1,x2,[x21,#112]
    orr x0,x0,x1
    orr x0,x0,x2
    ldp x1,x2,[x21,#128]
    orr x0,x0,x1
    orr x0,x0,x2
    ldp x1,x2,[x21,#144]
    orr x0,x0,x1
    orr x0,x0,x2
    ldp x1,x2,[x21,#160]
    orr x0,x0,x1
    orr x0,x0,x2
    ldp x1,x2,[x21,#176]
    orr x0,x0,x1
    orr x0,x0,x2
    cbz x0,.Lcanonical_identity
// InitializecandidateB0beforethefixednine-callfieldprogram.
    ldp x0,x1,[x21,#96]
    stp x0,x1,[x22,#0]
    ldp x0,x1,[x21,#112]
    stp x0,x1,[x22,#16]
    ldp x0,x1,[x21,#128]
    stp x0,x1,[x22,#32]
    ldp x0,x1,[x21,#144]
    stp x0,x1,[x22,#48]
    ldp x0,x1,[x21,#160]
    stp x0,x1,[x22,#64]
    ldp x0,x1,[x21,#176]
    stp x0,x1,[x22,#80]
    mov x0,x24
    mov x1,#9
    bl fapi_run_field_program
    mov x0,x21
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
    ldr x30,[sp,#48]
    ldp x21,x22,[sp,#16]
    ldp x23,x24,[sp,#32]
    ldp x19,x20,[sp],#64
    ret
.size fapi_trace_canonical,.-fapi_trace_canonical

.section .data
.p2align 4
// Thesezerosandallcode/descriptorbytesareexplicitlyinitialized/paidcold.
fapi_trace_workspace:
fapi_trace_raw:
    .zero 192
fapi_trace_aneg:
    .zero 48
fapi_trace_bvariants:
    .zero 384
fapi_trace_sum:
    .zero 48
fapi_trace_diff:
    .zero 48
.global fapi_trace_s
fapi_trace_s:
// ColdcomputesC=(1+u)^((p-1)/2)=(s,s),writesnativeMontgomerys48here.
    .zero 48
fapi_trace_field_program:
// DescriptorABI:(function,arg0,arg1,arg2),matchingtheactualleafsignature.
    .quad fapi_fp_neg,fapi_trace_raw+48,fapi_trace_aneg,0
    .quad fapi_fp_add,fapi_trace_raw+96,fapi_trace_raw+144,fapi_trace_sum
    .quad fapi_fp_sub,fapi_trace_raw+96,fapi_trace_raw+144,fapi_trace_diff
    .quad fapi_native_fp_mul,fapi_trace_bvariants+96,fapi_trace_sum,fapi_trace_s
    .quad fapi_native_fp_mul,fapi_trace_bvariants+144,fapi_trace_diff,fapi_trace_s
    .quad fapi_fp_neg,fapi_trace_raw+96,fapi_trace_bvariants+192,0
    .quad fapi_fp_neg,fapi_trace_raw+144,fapi_trace_bvariants+240,0
    .quad fapi_fp_neg,fapi_trace_bvariants+96,fapi_trace_bvariants+288,0
    .quad fapi_fp_neg,fapi_trace_bvariants+144,fapi_trace_bvariants+336,0

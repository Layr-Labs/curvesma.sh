// Inert analytic evidence only. Shared verifier credential-free Docker ONLY.
// x0 out192,x1 input192 (a0@0,b1@96 initialized),x2 LUTbase,
// x3 p48,x4 productsbase=const3728+528 (8p/4p/2p/p). Field input/output native Montgomery canonical Fp; constant3728 stores ordinary integers. Exact out/input alias safe.
// LUT,p,pPowers disjoint output/frame. Only a0@0,b1@96 output initialized. Frame336.
.text
.p2align 2
.global odd_trace_compact
.type odd_trace_compact,%function
odd_trace_compact:
    sub sp,sp,#336
    stp x29,x30,[sp]
    mov x29,sp
    stp x19,x20,[x29,#16]
    str x21,[x29,#32]
    mov x19,x0
    mov x20,x1
    mov x21,x3
    stp x2,x4,[x29,#152]
    ldp x0,x1,[x20]
    stp x0,x1,[x19]
    ldp x0,x1,[x20,#16]
    stp x0,x1,[x19,#16]
    ldp x0,x1,[x20,#32]
    stp x0,x1,[x19,#32]
    ldp x4,x5,[x20,#48]
    ldp x6,x7,[x20,#64]
    ldp x8,x9,[x20,#80]
    orr x10,x4,x5
    orr x10,x10,x6
    orr x10,x10,x7
    orr x10,x10,x8
    orr x10,x10,x9
    cbnz x10,.Lodd_nonzero_neg
    stp xzr,xzr,[x19,#48]
    stp xzr,xzr,[x19,#64]
    stp xzr,xzr,[x19,#80]
    b .Lodd_neg_done
.Lodd_nonzero_neg:
    ldp x10,x11,[x21,#0]
    ldp x12,x13,[x21,#16]
    ldp x14,x15,[x21,#32]
    subs x4,x10,x4
    sbcs x5,x11,x5
    sbcs x6,x12,x6
    sbcs x7,x13,x7
    sbcs x8,x14,x8
    sbcs x9,x15,x9
    stp x4,x5,[x19,#48]
    stp x6,x7,[x19,#64]
    stp x8,x9,[x19,#80]
.Lodd_neg_done:
    add x0,x29,#192
    add x1,x20,#96
    add x2,x20,#144
    mov x3,x21
    bl fp381_add
    add x0,x29,#240
    add x1,x20,#96
    add x2,x20,#144
    mov x3,x21
    bl fp381_sub
    add x0,x19,#96
    add x1,x29,#192
    ldp x2,x3,[x29,#152]
    bl fixed_s_lookup
    add x0,x19,#144
    add x1,x29,#240
    ldp x2,x3,[x29,#152]
    bl fixed_s_lookup
    ldp x19,x20,[x29,#16]
    ldr x21,[x29,#32]
    ldp x29,x30,[sp],#336
    ret
.size odd_trace_compact,.-odd_trace_compact

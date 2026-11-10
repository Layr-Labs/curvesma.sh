// Inert analytic evidence only. Shared verifier credential-free Docker ONLY.
// x0 out576,x1 input576 (onlya0@0,b1@384 initialized),x2 LUTbase,
// x3 p48,x4 pPowers192 (8p/4p/2p/p). Field input/output native Montgomery canonical Fp; pPowers are ordinary integers. Exact out/input alias safe.
// LUT,p,pPowers disjoint output/frame. Only a0@0,b1@384 output initialized. Frame336.
.text
.p2align 2
.global odd_trace
.type odd_trace,%function
odd_trace:
    sub sp,sp,#336
    stp x29,x30,[sp]
    mov x29,sp
    stp x19,x20,[x29,#16]
    stp x21,x22,[x29,#32]
    stp x23,x24,[x29,#48]
    stp x25,x26,[x29,#64]
    stp x27,x28,[x29,#80]
    stp x0,x1,[x29,#128]
    stp x2,x3,[x29,#144]
    stp x4,xzr,[x29,#160]
    ldp x19,x20,[x29,#128]
    ldr x21,[x29,#152]
    ldr x22,[x29,#144]
    ldp x0,x1,[x20]
    stp x0,x1,[x19]
    ldp x0,x1,[x20,#16]
    stp x0,x1,[x19,#16]
    ldp x0,x1,[x20,#32]
    stp x0,x1,[x19,#32]
    stp xzr,xzr,[x29,#288]
    stp xzr,xzr,[x29,#304]
    stp xzr,xzr,[x29,#320]
    add x0,x19,#48
    add x1,x29,#288
    add x2,x20,#48
    mov x3,x21
    bl fp381_sub
    add x0,x29,#192
    add x1,x20,#384
    add x2,x20,#432
    mov x3,x21
    bl fp381_add
    add x0,x29,#240
    add x1,x20,#384
    add x2,x20,#432
    mov x3,x21
    bl fp381_sub
    add x0,x19,#384
    add x1,x29,#192
    mov x2,x22
    ldr x3,[x29,#160]
    bl fixed_s_lookup
    add x0,x19,#432
    add x1,x29,#240
    mov x2,x22
    ldr x3,[x29,#160]
    bl fixed_s_lookup
    ldp x19,x20,[x29,#16]
    ldp x21,x22,[x29,#32]
    ldp x23,x24,[x29,#48]
    ldp x25,x26,[x29,#64]
    ldp x27,x28,[x29,#80]
    ldp x29,x30,[sp],#336
    ret
.size odd_trace,.-odd_trace

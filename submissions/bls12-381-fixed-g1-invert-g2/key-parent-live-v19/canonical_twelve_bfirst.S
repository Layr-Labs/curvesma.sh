// Inert analytic evidence; NEVER host compile/execute or use credentials.
// x0 out192 disjoint, x1 rawPhi192 read-only, x2 sLUT, x3 p48,x4 const3728.
// w0=1 identity/out invalid; w0=0 all192 bytes initialized. GT inputs only.
// Compare signed b1 representatives by unsigned lex12 LE64 words, increasing index.
// Store corresponding a0 THEN selected b1. T^2=-1 proves no nonidentity tie.
// Frame384 header96/EvenB96@96/OddPhi192@192;19..29/SP/LR preserved,no18.
.text
.p2align 2
.global canonical_twelve_bfirst
.type canonical_twelve_bfirst,%function
canonical_twelve_bfirst:
    sub sp,sp,#384
    stp x29,x30,[sp]
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    str x23,[sp,#48]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    mov x22,x3
    mov x23,x4
    add x0,x20,#96
    add x1,sp,#96
    bl canonical_conjugate_b1
    add x0,sp,#192
    mov x1,x20
    mov x2,x21
    mov x3,x22
    mov x4,x23
    bl odd_trace_compact
    add x0,sp,#288
    mov x1,x0
    bl canonical_conjugate_b1
    add x4,sp,#96
    add x5,sp,#288
    ldp x0,x1,[x4,#0]
    ldp x2,x3,[x5,#0]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x4,#16]
    ldp x2,x3,[x5,#16]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x4,#32]
    ldp x2,x3,[x5,#32]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x4,#48]
    ldp x2,x3,[x5,#48]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x4,#64]
    ldp x2,x3,[x5,#64]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x4,#80]
    ldp x2,x3,[x5,#80]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    mov w0,#1
    b .Lb_done
.Lb_even:
    mov x6,x20
    mov x7,x4
    b .Lb_copy
.Lb_odd:
    add x6,sp,#192
    mov x7,x5
.Lb_copy:
    ldp x0,x1,[x6,#0]
    stp x0,x1,[x19,#0]
    ldp x0,x1,[x6,#16]
    stp x0,x1,[x19,#16]
    ldp x0,x1,[x6,#32]
    stp x0,x1,[x19,#32]
    ldp x0,x1,[x6,#48]
    stp x0,x1,[x19,#48]
    ldp x0,x1,[x6,#64]
    stp x0,x1,[x19,#64]
    ldp x0,x1,[x6,#80]
    stp x0,x1,[x19,#80]
    ldp x0,x1,[x7,#0]
    stp x0,x1,[x19,#96]
    ldp x0,x1,[x7,#16]
    stp x0,x1,[x19,#112]
    ldp x0,x1,[x7,#32]
    stp x0,x1,[x19,#128]
    ldp x0,x1,[x7,#48]
    stp x0,x1,[x19,#144]
    ldp x0,x1,[x7,#64]
    stp x0,x1,[x19,#160]
    ldp x0,x1,[x7,#80]
    stp x0,x1,[x19,#176]
    mov w0,#0
.Lb_done:
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldr x23,[sp,#48]
    ldp x29,x30,[sp]
    add sp,sp,#384
    ret
.size canonical_twelve_bfirst,.-canonical_twelve_bfirst

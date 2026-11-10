// Inert analytic evidence; NEVER host compile/execute or use credentials.
// x0 out192 disjoint, x1 rawPhi192 read-only, x2 sLUT, x3 p48,x4 powers192.
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
    mov x29,sp
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    mov x22,x3
    mov x23,x4
    add x0,x20,#96
    add x1,x29,#96
    bl canonical_conjugate_b1
    add x0,x29,#192
    mov x1,x20
    mov x2,x21
    mov x3,x22
    mov x4,x23
    bl odd_trace_compact
    add x0,x29,#288
    mov x1,x0
    bl canonical_conjugate_b1
    add x24,x29,#96
    add x25,x29,#288
    ldp x0,x1,[x24,#0]
    ldp x2,x3,[x25,#0]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x24,#16]
    ldp x2,x3,[x25,#16]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x24,#32]
    ldp x2,x3,[x25,#32]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x24,#48]
    ldp x2,x3,[x25,#48]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x24,#64]
    ldp x2,x3,[x25,#64]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    ldp x0,x1,[x24,#80]
    ldp x2,x3,[x25,#80]
    cmp x0,x2
    b.lo .Lb_even
    b.hi .Lb_odd
    cmp x1,x3
    b.lo .Lb_even
    b.hi .Lb_odd
    mov w0,#1
    b .Lb_done
.Lb_even:
    mov x26,x20
    mov x27,x24
    b .Lb_copy
.Lb_odd:
    add x26,x29,#192
    mov x27,x25
.Lb_copy:
    ldp x0,x1,[x26,#0]
    stp x0,x1,[x19,#0]
    ldp x0,x1,[x26,#16]
    stp x0,x1,[x19,#16]
    ldp x0,x1,[x26,#32]
    stp x0,x1,[x19,#32]
    ldp x0,x1,[x26,#48]
    stp x0,x1,[x19,#48]
    ldp x0,x1,[x26,#64]
    stp x0,x1,[x19,#64]
    ldp x0,x1,[x26,#80]
    stp x0,x1,[x19,#80]
    ldp x0,x1,[x27,#0]
    stp x0,x1,[x19,#96]
    ldp x0,x1,[x27,#16]
    stp x0,x1,[x19,#112]
    ldp x0,x1,[x27,#32]
    stp x0,x1,[x19,#128]
    ldp x0,x1,[x27,#48]
    stp x0,x1,[x19,#144]
    ldp x0,x1,[x27,#64]
    stp x0,x1,[x19,#160]
    ldp x0,x1,[x27,#80]
    stp x0,x1,[x19,#176]
    mov w0,#0
.Lb_done:
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp]
    add sp,sp,#384
    ret
.size canonical_twelve_bfirst,.-canonical_twelve_bfirst

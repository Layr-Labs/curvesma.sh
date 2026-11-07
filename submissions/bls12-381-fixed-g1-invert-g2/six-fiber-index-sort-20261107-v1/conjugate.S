// Static untrusted evidence;never host-compile/execute or use credentials.
// x0=b1 native canonical Montgomery96 bytes,x1=output96 bytes,exact in-place or disjoint output.
// Select lexicographically smaller native-integer pair of b1,-b1.
// First nonzero canonical Fp component determines sign using ceil(p/2).
// Only native64 operations;both zero input maps to zero.
.text
.p2align 2
.global canonical_conjugate_b1
.type canonical_conjugate_b1,%function
canonical_conjugate_b1:
    stp x19,x20,[sp,#-80]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    mov x19,x0
    mov x20,x1
    adrp x22,conjugate_p
    add x22,x22,:lo12:conjugate_p
    adrp x23,conjugate_halfceil
    add x23,x23,:lo12:conjugate_halfceil
    ldp x2,x3,[x19,#0]
    ldp x4,x5,[x19,#16]
    ldp x6,x7,[x19,#32]
    orr x8,x2,x3
    orr x8,x8,x4
    orr x8,x8,x5
    orr x8,x8,x6
    orr x8,x8,x7
    cmp x8,#0
    b.ne .Lconjugate_chosen
    ldp x2,x3,[x19,#48]
    ldp x4,x5,[x19,#64]
    ldp x6,x7,[x19,#80]
.Lconjugate_chosen:
    ldp x9,x10,[x23,#0]
    ldp x11,x12,[x23,#16]
    ldp x13,x14,[x23,#32]
    subs xzr,x2,x9
    sbcs xzr,x3,x10
    sbcs xzr,x4,x11
    sbcs xzr,x5,x12
    sbcs xzr,x6,x13
    sbcs xzr,x7,x14
    cset x21,cs
    ldp x9,x10,[x22,#0]
    ldp x11,x12,[x22,#16]
    ldp x13,x14,[x22,#32]
    // Complete canonical Fp component at offset0.
    ldp x2,x3,[x19,#0]
    ldp x4,x5,[x19,#16]
    ldp x6,x7,[x19,#32]
    orr x8,x2,x3
    orr x8,x8,x4
    orr x8,x8,x5
    orr x8,x8,x6
    orr x8,x8,x7
    cmp x8,#0
    cset x27,ne
    and x27,x27,x21
    subs x15,x9,x2
    sbcs x16,x10,x3
    sbcs x17,x11,x4
    sbcs x24,x12,x5
    sbcs x25,x13,x6
    sbcs x26,x14,x7
    cmp x27,#0
    csel x2,x15,x2,ne
    csel x3,x16,x3,ne
    csel x4,x17,x4,ne
    csel x5,x24,x5,ne
    csel x6,x25,x6,ne
    csel x7,x26,x7,ne
    stp x2,x3,[x20,#0]
    stp x4,x5,[x20,#16]
    stp x6,x7,[x20,#32]
    // Complete canonical Fp component at offset48.
    ldp x2,x3,[x19,#48]
    ldp x4,x5,[x19,#64]
    ldp x6,x7,[x19,#80]
    orr x8,x2,x3
    orr x8,x8,x4
    orr x8,x8,x5
    orr x8,x8,x6
    orr x8,x8,x7
    cmp x8,#0
    cset x27,ne
    and x27,x27,x21
    subs x15,x9,x2
    sbcs x16,x10,x3
    sbcs x17,x11,x4
    sbcs x24,x12,x5
    sbcs x25,x13,x6
    sbcs x26,x14,x7
    cmp x27,#0
    csel x2,x15,x2,ne
    csel x3,x16,x3,ne
    csel x4,x17,x4,ne
    csel x5,x24,x5,ne
    csel x6,x25,x6,ne
    csel x7,x26,x7,ne
    stp x2,x3,[x20,#48]
    stp x4,x5,[x20,#64]
    stp x6,x7,[x20,#80]
    ldp x27,x28,[sp,#64]
    ldp x25,x26,[sp,#48]
    ldp x23,x24,[sp,#32]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#80
    ret
.size canonical_conjugate_b1,.-canonical_conjugate_b1
.section .rodata
.p2align 3
conjugate_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624,0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
conjugate_halfceil:
    .quad 0xdcff7fffffffd556,0xf55ffff58a9ffff,0xb39869507b587b12,0xb23ba5c279c2895f,0x258dd3db21a5d66b,0xd0088f51cbff34d

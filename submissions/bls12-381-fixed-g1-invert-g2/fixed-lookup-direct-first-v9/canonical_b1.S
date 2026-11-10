// Inert analytic source; ONLY shared credential-free llm-verifier may execute.
// x0 canonical native b1_96,x1 out96,exact alias or disjoint. Bothcomponents canonical Fp.
// Same first-nonzero sign as attributed public helper; new caller-saved-only leaf.
// No SP/x18/x19..30 changes, no BL. All output initialized; partial overlap forbidden.
.text
.p2align 2
.global canonical_conjugate_b1
.type canonical_conjugate_b1,%function
canonical_conjugate_b1:
    ldp x2,x3,[x0,#0]
    ldp x4,x5,[x0,#16]
    ldp x6,x7,[x0,#32]
    orr x14,x2,x3
    orr x14,x14,x4
    orr x14,x14,x5
    orr x14,x14,x6
    orr x14,x14,x7
    cbnz x14,.Lsign_chosen
    ldp x2,x3,[x0,#48]
    ldp x4,x5,[x0,#64]
    ldp x6,x7,[x0,#80]
.Lsign_chosen:
    adrp x15,conjugate_halfceil
    add x15,x15,:lo12:conjugate_halfceil
    ldp x8,x9,[x15,#0]
    ldp x10,x11,[x15,#16]
    ldp x12,x13,[x15,#32]
    subs xzr,x2,x8
    sbcs xzr,x3,x9
    sbcs xzr,x4,x10
    sbcs xzr,x5,x11
    sbcs xzr,x6,x12
    sbcs xzr,x7,x13
    b.lo .Lsign_copy
    sub x15,x15,#48
    ldp x8,x9,[x15,#0]
    ldp x10,x11,[x15,#16]
    ldp x12,x13,[x15,#32]
    cbz x14,.Lsign_real_zero_high
    subs x2,x8,x2
    sbcs x3,x9,x3
    sbcs x4,x10,x4
    sbcs x5,x11,x5
    sbcs x6,x12,x6
    sbcs x7,x13,x7
    stp x2,x3,[x1,#0]
    stp x4,x5,[x1,#16]
    stp x6,x7,[x1,#32]
    ldp x2,x3,[x0,#48]
    ldp x4,x5,[x0,#64]
    ldp x6,x7,[x0,#80]
    orr x14,x2,x3
    orr x14,x14,x4
    orr x14,x14,x5
    orr x14,x14,x6
    orr x14,x14,x7
    cbnz x14,.Lsign_imag_nonzero_high
    stp xzr,xzr,[x1,#48]
    stp xzr,xzr,[x1,#64]
    stp xzr,xzr,[x1,#80]
    ret
.Lsign_real_zero_high:
    stp xzr,xzr,[x1,#0]
    stp xzr,xzr,[x1,#16]
    stp xzr,xzr,[x1,#32]
.Lsign_imag_nonzero_high:
    subs x2,x8,x2
    sbcs x3,x9,x3
    sbcs x4,x10,x4
    sbcs x5,x11,x5
    sbcs x6,x12,x6
    sbcs x7,x13,x7
    stp x2,x3,[x1,#48]
    stp x4,x5,[x1,#64]
    stp x6,x7,[x1,#80]
    ret
.Lsign_copy:
    ldp x2,x3,[x0,#0]
    stp x2,x3,[x1,#0]
    ldp x2,x3,[x0,#16]
    stp x2,x3,[x1,#16]
    ldp x2,x3,[x0,#32]
    stp x2,x3,[x1,#32]
    ldp x2,x3,[x0,#48]
    stp x2,x3,[x1,#48]
    ldp x2,x3,[x0,#64]
    stp x2,x3,[x1,#64]
    ldp x2,x3,[x0,#80]
    stp x2,x3,[x1,#80]
    ret
.size canonical_conjugate_b1,.-canonical_conjugate_b1
.section .rodata
.p2align 3
conjugate_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624,0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
conjugate_halfceil:
    .quad 0xdcff7fffffffd556,0xf55ffff58a9ffff,0xb39869507b587b12,0xb23ba5c279c2895f,0x258dd3db21a5d66b,0xd0088f51cbff34d


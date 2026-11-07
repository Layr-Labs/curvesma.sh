// Static analytic evidence only; never host-executed.
// Full canonical Montgomery Fp, x0=a48,x1=b48,x2=out48. Scalar caller scratch.
.text
.p2align 2
.global fapi_fp_add
.type fapi_fp_add,%function
fapi_fp_add:
    ldp x3,x4,[x0,#0]
    ldp x5,x6,[x0,#16]
    ldp x7,x8,[x0,#32]
    ldp x9,x10,[x1,#0]
    ldp x11,x12,[x1,#16]
    ldp x13,x14,[x1,#32]
    adds x3,x3,x9
    adcs x4,x4,x10
    adcs x5,x5,x11
    adcs x6,x6,x12
    adcs x7,x7,x13
    adc x8,x8,x14
    adrp x16,fapi_mod_p
    add x16,x16,:lo12:fapi_mod_p
    ldp x9,x10,[x16,#0]
    ldp x11,x12,[x16,#16]
    ldp x13,x14,[x16,#32]
    subs x9,x3,x9
    sbcs x10,x4,x10
    sbcs x11,x5,x11
    sbcs x12,x6,x12
    sbcs x13,x7,x13
    sbcs x14,x8,x14
    b.lo .Lfp_add_unreduced
    stp x9,x10,[x2,#0]
    stp x11,x12,[x2,#16]
    stp x13,x14,[x2,#32]
    ret
.Lfp_add_unreduced:
    stp x3,x4,[x2,#0]
    stp x5,x6,[x2,#16]
    stp x7,x8,[x2,#32]
    ret
.size fapi_fp_add,.-fapi_fp_add
.p2align 2
.global fapi_fp_sub
.type fapi_fp_sub,%function
fapi_fp_sub:
    ldp x3,x4,[x0,#0]
    ldp x5,x6,[x0,#16]
    ldp x7,x8,[x0,#32]
    ldp x9,x10,[x1,#0]
    ldp x11,x12,[x1,#16]
    ldp x13,x14,[x1,#32]
    subs x3,x3,x9
    sbcs x4,x4,x10
    sbcs x5,x5,x11
    sbcs x6,x6,x12
    sbcs x7,x7,x13
    sbcs x8,x8,x14
    b.hs .Lfp_sub_store
    adrp x16,fapi_mod_p
    add x16,x16,:lo12:fapi_mod_p
    ldp x9,x10,[x16,#0]
    ldp x11,x12,[x16,#16]
    ldp x13,x14,[x16,#32]
    adds x3,x3,x9
    adcs x4,x4,x10
    adcs x5,x5,x11
    adcs x6,x6,x12
    adcs x7,x7,x13
    adc x8,x8,x14
.Lfp_sub_store:
    stp x3,x4,[x2,#0]
    stp x5,x6,[x2,#16]
    stp x7,x8,[x2,#32]
    ret
.size fapi_fp_sub,.-fapi_fp_sub
// Static analytic evidence only; never host-executed.
// x0 points to canonical Montgomery a in [0,p); x1 points to disjoint output48.
// Returns canonical -a mod p. Caller scratch only; x18 and callee state untouched.
.text
.p2align 2
.global fapi_fp_neg
.type fapi_fp_neg,%function
fapi_fp_neg:
    ldp x2,x3,[x0,#0]
    ldp x4,x5,[x0,#16]
    ldp x6,x7,[x0,#32]
    orr x14,x2,x3
    orr x14,x14,x4
    orr x14,x14,x5
    orr x14,x14,x6
    orr x14,x14,x7
    cbz x14,.Lfp_neg_zero
    adrp x14,fapi_mod_p
    add x14,x14,:lo12:fapi_mod_p
    ldp x8,x9,[x14,#0]
    ldp x10,x11,[x14,#16]
    ldp x12,x13,[x14,#32]
    subs x8,x8,x2
    sbcs x9,x9,x3
    sbcs x10,x10,x4
    sbcs x11,x11,x5
    sbcs x12,x12,x6
    sbc x13,x13,x7
    stp x8,x9,[x1,#0]
    stp x10,x11,[x1,#16]
    stp x12,x13,[x1,#32]
    ret
.Lfp_neg_zero:
    stp xzr,xzr,[x1,#0]
    stp xzr,xzr,[x1,#16]
    stp xzr,xzr,[x1,#32]
    ret
.size fapi_fp_neg,.-fapi_fp_neg

.section .rodata
.p2align 3
fapi_mod_p:
    .quad 0xb9feffffffffaaab
    .quad 0x1eabfffeb153ffff
    .quad 0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf
    .quad 0x4b1ba7b6434bacd7
    .quad 0x1a0111ea397fe69a

// Inert analytic evidence only; never execute on host or with credentials.
// x0=out48,x1=a48,x2=b48,x3=p48;canonical six-limb native Fp.
.text
.p2align 2
.global fp381_add
.type fp381_add,%function
fp381_add:
    ldp x4,x5,[x1,#0]
    ldp x6,x7,[x1,#16]
    ldp x8,x9,[x1,#32]
    ldp x10,x11,[x2,#0]
    ldp x12,x13,[x2,#16]
    ldp x14,x15,[x2,#32]
    adds x4,x4,x10
    adcs x5,x5,x11
    adcs x6,x6,x12
    adcs x7,x7,x13
    adcs x8,x8,x14
    adc x9,x9,x15
    ldp x1,x2,[x3,#0]
    subs x10,x4,x1
    sbcs x11,x5,x2
    ldp x1,x2,[x3,#16]
    sbcs x12,x6,x1
    sbcs x13,x7,x2
    ldp x1,x2,[x3,#32]
    sbcs x14,x8,x1
    sbcs x15,x9,x2
    csel x4,x10,x4,hs
    csel x5,x11,x5,hs
    csel x6,x12,x6,hs
    csel x7,x13,x7,hs
    csel x8,x14,x8,hs
    csel x9,x15,x9,hs
    stp x4,x5,[x0,#0]
    stp x6,x7,[x0,#16]
    stp x8,x9,[x0,#32]
    ret
.size fp381_add,.-fp381_add

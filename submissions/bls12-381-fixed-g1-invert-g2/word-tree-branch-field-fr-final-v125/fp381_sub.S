// Inert analytic evidence only; never execute on host or with credentials.
// x0=out48,x1=a48,x2=b48,x3=p48;canonical six-limb native Fp.
.text
.p2align 2
.global fp381_sub
.type fp381_sub,%function
fp381_sub:
    ldp x4,x5,[x1,#0]
    ldp x6,x7,[x1,#16]
    ldp x8,x9,[x1,#32]
    ldp x10,x11,[x2,#0]
    ldp x12,x13,[x2,#16]
    ldp x14,x15,[x2,#32]
    subs x4,x4,x10
    sbcs x5,x5,x11
    sbcs x6,x6,x12
    sbcs x7,x7,x13
    sbcs x8,x8,x14
    sbcs x9,x9,x15
    b.hs .Lsub_store
    ldp x10,x11,[x3,#0]
    ldp x12,x13,[x3,#16]
    ldp x14,x15,[x3,#32]
    adds x4,x4,x10
    adcs x5,x5,x11
    adcs x6,x6,x12
    adcs x7,x7,x13
    adcs x8,x8,x14
    adc x9,x9,x15
.Lsub_store:
    stp x4,x5,[x0,#0]
    stp x6,x7,[x0,#16]
    stp x8,x9,[x0,#32]
    ret
.size fp381_sub,.-fp381_sub

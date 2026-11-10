// Inert analytic evidence; shared credential-free verifier ONLY.
// x0 out32,x1 canonical a32,x2 canonical b32,x3 full protected r32.
// Since r<2^255, a+b<2^256. All inputs read before writes; aliases safe.
.text
.p2align 2
.global fr255_add
.type fr255_add,%function
fr255_add:
    ldp x4,x5,[x1,#0]
    ldp x6,x7,[x1,#16]
    ldp x8,x9,[x2,#0]
    ldp x10,x11,[x2,#16]
    adds x4,x4,x8
    adcs x5,x5,x9
    adcs x6,x6,x10
    adc x7,x7,x11
    ldp x1,x2,[x3,#0]
    subs x8,x4,x1
    sbcs x9,x5,x2
    ldp x1,x2,[x3,#16]
    sbcs x10,x6,x1
    sbcs x11,x7,x2
    b.lo .Lfr_original
    stp x8,x9,[x0,#0]
    stp x10,x11,[x0,#16]
    ret
.Lfr_original:
    stp x4,x5,[x0,#0]
    stp x6,x7,[x0,#16]
    ret
.size fr255_add,.-fr255_add

// Static analytic leaf,neverhostexecuted.
// x0oldrecord,x1newrecord,x2jump; A is known by counter/index,not stored.
// Jump has fixedu1 and canonicalv; only B needs modularFr addition.
.text
.p2align 2
.global fapi_scalar_update
.type fapi_scalar_update, %function
fapi_scalar_update:
    add x13,x0,#576
    add x14,x1,#576
    add x15,x2,#576
    adrp x16,fapi_r
    add x16,x16,:lo12:fapi_r
    ldp x0,x1,[x13,#32]
    ldp x2,x3,[x13,#48]
    ldp x4,x5,[x15,#32]
    ldp x6,x7,[x15,#48]
    adds x0,x0,x4
    adcs x1,x1,x5
    adcs x2,x2,x6
    adc x3,x3,x7
    ldp x4,x5,[x16]
    ldp x6,x7,[x16,#16]
    subs x9,x0,x4
    sbcs x10,x1,x5
    sbcs x11,x2,x6
    sbcs x12,x3,x7
    csel x0,x9,x0,hs
    csel x1,x10,x1,hs
    csel x2,x11,x2,hs
    csel x3,x12,x3,hs
    stp x0,x1,[x14,#32]
    stp x2,x3,[x14,#48]
    ret
.size fapi_scalar_update,.-fapi_scalar_update

.p2align 2
.section .rodata
.p2align 3
fapi_r:
    .quad 0xffffffff00000001
    .quad 0x53bda402fffe5bfe
    .quad 0x3339d80809a1d805
    .quad 0x73eda753299d7d48

// Inert analytic evidence; never execute on host or with credentials.
// Original straight-line CIOS specification; a,b<p, p is protected BLS381 modulus.
// x0=out,x1=a,x2=b,x3=p,x4=n0. Native R=2^384, canonical output.
.text
.p2align 2
.macro CIOS_ACCUM dst, av, bv
    mul x26,\av,\bv
    umulh x27,\av,\bv
    adds x26,x26,\dst
    adc x27,x27,xzr
    adds \dst,x26,x28
    adc x28,x27,xzr
.endm
.global cios381_mul
.type cios381_mul,%function
cios381_mul:
    stp x19,x20,[sp,#-96]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    stp x29,x30,[sp,#80]
    ldp x5,x6,[x1,#0]
    ldp x7,x8,[x1,#16]
    ldp x9,x10,[x1,#32]
    ldp x11,x12,[x3,#0]
    ldp x13,x14,[x3,#16]
    ldp x15,x16,[x3,#32]
    mov x19,xzr
    mov x20,xzr
    mov x21,xzr
    mov x22,xzr
    mov x23,xzr
    mov x24,xzr
    mov x25,xzr
.irp round,0,1,2,3,4,5
    ldr x29,[x2,#(8*\round)]
    mov x28,xzr
    CIOS_ACCUM x19,x5,x29
    CIOS_ACCUM x20,x6,x29
    CIOS_ACCUM x21,x7,x29
    CIOS_ACCUM x22,x8,x29
    CIOS_ACCUM x23,x9,x29
    CIOS_ACCUM x24,x10,x29
    add x25,x25,x28
    mul x3,x19,x4
    mov x28,xzr
    CIOS_ACCUM x19,x11,x3
    CIOS_ACCUM x20,x12,x3
    CIOS_ACCUM x21,x13,x3
    CIOS_ACCUM x22,x14,x3
    CIOS_ACCUM x23,x15,x3
    CIOS_ACCUM x24,x16,x3
    add x25,x25,x28
    mov x19,x20
    mov x20,x21
    mov x21,x22
    mov x22,x23
    mov x23,x24
    mov x24,x25
    mov x25,xzr
.endr
    subs x5,x19,x11
    sbcs x6,x20,x12
    sbcs x7,x21,x13
    sbcs x8,x22,x14
    sbcs x9,x23,x15
    sbcs x10,x24,x16
    csel x19,x5,x19,hs
    csel x20,x6,x20,hs
    csel x21,x7,x21,hs
    csel x22,x8,x22,hs
    csel x23,x9,x23,hs
    csel x24,x10,x24,hs
    stp x19,x20,[x0,#0]
    stp x21,x22,[x0,#16]
    stp x23,x24,[x0,#32]
    ldp x21,x22,[sp,#16]
    ldp x23,x24,[sp,#32]
    ldp x25,x26,[sp,#48]
    ldp x27,x28,[sp,#64]
    ldp x29,x30,[sp,#80]
    ldp x19,x20,[sp],#96
    ret
.size cios381_mul,.-cios381_mul

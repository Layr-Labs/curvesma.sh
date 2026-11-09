// Inert analytic evidence only: never execute on host or with credentials.
// CIOS with fixed register renaming. Low-word cancellation identity also appears
// in pinned blst/build/elf/mul_mont_384-armv8.S; see proof attribution.
// x0=out,x1=a,x2=b,x3=p,x4=n0; canonical a,b<p, native R=2^384.
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
.macro CIOS_ROTROW off, t0, t1, t2, t3, t4, t5, t6
    ldr x29,[x2,#\off]
    mov x28,xzr
    CIOS_ACCUM \t0,x5,x29
    CIOS_ACCUM \t1,x6,x29
    CIOS_ACCUM \t2,x7,x29
    CIOS_ACCUM \t3,x8,x29
    CIOS_ACCUM \t4,x9,x29
    CIOS_ACCUM \t5,x10,x29
    add \t6,\t6,x28
    mul x3,\t0,x4
    umulh x28,x3,x11
    subs xzr,\t0,#1
    adc x28,x28,xzr
    mov \t0,xzr
    CIOS_ACCUM \t1,x12,x3
    CIOS_ACCUM \t2,x13,x3
    CIOS_ACCUM \t3,x14,x3
    CIOS_ACCUM \t4,x15,x3
    CIOS_ACCUM \t5,x16,x3
    add \t6,\t6,x28
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
    CIOS_ROTROW 0,x19,x20,x21,x22,x23,x24,x25
    CIOS_ROTROW 8,x20,x21,x22,x23,x24,x25,x19
    CIOS_ROTROW 16,x21,x22,x23,x24,x25,x19,x20
    CIOS_ROTROW 24,x22,x23,x24,x25,x19,x20,x21
    CIOS_ROTROW 32,x23,x24,x25,x19,x20,x21,x22
    CIOS_ROTROW 40,x24,x25,x19,x20,x21,x22,x23
    subs x5,x25,x11
    sbcs x6,x19,x12
    sbcs x7,x20,x13
    sbcs x8,x21,x14
    sbcs x9,x22,x15
    sbcs x10,x23,x16
    csel x25,x5,x25,hs
    csel x19,x6,x19,hs
    csel x20,x7,x20,hs
    csel x21,x8,x21,hs
    csel x22,x9,x22,hs
    csel x23,x10,x23,hs
    stp x25,x19,[x0,#0]
    stp x20,x21,[x0,#16]
    stp x22,x23,[x0,#32]
    ldp x21,x22,[sp,#16]
    ldp x23,x24,[sp,#32]
    ldp x25,x26,[sp,#48]
    ldp x27,x28,[sp,#64]
    ldp x29,x30,[sp,#80]
    ldp x19,x20,[sp],#96
    ret
.size cios381_mul,.-cios381_mul

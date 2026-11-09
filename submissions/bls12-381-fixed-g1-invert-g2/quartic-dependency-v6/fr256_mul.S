// Inert analytic evidence only: never execute on host or with credentials.
// x0=out32,x1=a32,x2=b32,x3=r32,x4=n0=0xfffffffeffffffff.
// Canonical a,b<r, output a*b/2^256 mod r. Common scale preserves determinant equality.
// Standard CIOS; low-word cancellation credited to pinned blst, see proof.
.text
.p2align 2
.macro F4_ACC dst, av, bv
    mul x20,\av,\bv
    umulh x21,\av,\bv
    adds x20,x20,\dst
    adc x21,x21,xzr
    adds \dst,x20,x22
    adc x22,x21,xzr
.endm
.macro F4_ROW off, t0, t1, t2, t3, t4
    ldr x19,[x2,#\off]
    mov x22,xzr
    F4_ACC \t0,x5,x19
    F4_ACC \t1,x6,x19
    F4_ACC \t2,x7,x19
    F4_ACC \t3,x8,x19
    add \t4,\t4,x22
    mul x3,\t0,x4
    umulh x22,x3,x9
    subs xzr,\t0,#1
    adc x22,x22,xzr
    mov \t0,xzr
    F4_ACC \t1,x10,x3
    F4_ACC \t2,x11,x3
    F4_ACC \t3,x12,x3
    add \t4,\t4,x22
.endm
.global fr256_mul
.type fr256_mul,%function
fr256_mul:
    stp x19,x20,[sp,#-32]!
    stp x21,x22,[sp,#16]
    ldp x5,x6,[x1]
    ldp x7,x8,[x1,#16]
    ldp x9,x10,[x3]
    ldp x11,x12,[x3,#16]
    mov x13,xzr
    mov x14,xzr
    mov x15,xzr
    mov x16,xzr
    mov x17,xzr
    F4_ROW 0,x13,x14,x15,x16,x17
    F4_ROW 8,x14,x15,x16,x17,x13
    F4_ROW 16,x15,x16,x17,x13,x14
    F4_ROW 24,x16,x17,x13,x14,x15
    subs x5,x17,x9
    sbcs x6,x13,x10
    sbcs x7,x14,x11
    sbcs x8,x15,x12
    csel x17,x5,x17,hs
    csel x13,x6,x13,hs
    csel x14,x7,x14,hs
    csel x15,x8,x15,hs
    stp x17,x13,[x0]
    stp x14,x15,[x0,#16]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#32
    ret
.size fr256_mul,.-fr256_mul

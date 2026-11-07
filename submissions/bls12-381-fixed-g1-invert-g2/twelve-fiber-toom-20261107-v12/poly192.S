// Static analytic evidence. Never host-compile/execute or use credentials.
// x0=input192,x1=output32,x2=three canonical little-endian Fq keys.
// q=2^64-59. Full192-bit fingerprint in first24 output bytes; last8 zero.
.text
.p2align 2
.macro POLY_STEP acc,key
    mov x0,\acc
    mov x1,\key
    mov x2,x28
    bl .Lpoly_madd
    mov \acc,x0
.endm
.macro POLY_COEFF i
    ldr x10,[x19,#(8*((63*\i)/64))]
.if ((63*\i)%64) != 0
    lsr x10,x10,#((63*\i)%64)
    ldr x11,[x19,#(8*((63*\i)/64+1))]
    orr x10,x10,x11,lsl #(64-((63*\i)%64))
.endif
    and x28,x10,#0x7fffffffffffffff
    POLY_STEP x22,x25
    POLY_STEP x23,x26
    POLY_STEP x24,x27
.endm
.global fapi_poly192
.type fapi_poly192,%function
fapi_poly192:
    stp x19,x20,[sp,#-96]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    stp x29,x30,[sp,#80]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    ldp x25,x26,[x21,#0]
    ldr x27,[x21,#16]
    mov x3,#-59
    mov x4,#59
    ldr x10,[x19,#184]
    lsr x22,x10,#40
    mov x23,x22
    mov x24,x22
.irp i,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
    POLY_COEFF \i
.endr
    stp x22,x23,[x20,#0]
    stp x24,xzr,[x20,#16]
    ldp x29,x30,[sp,#80]
    ldp x27,x28,[sp,#64]
    ldp x25,x26,[sp,#48]
    ldp x23,x24,[sp,#32]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#96
    ret
.size fapi_poly192,.-fapi_poly192

// x0=A<q,x1=K<q,x2=c<2^63,x3=q,x4=59,return x0=(A*K+c) modq.
// x3/x4 and all callee state are preserved; scratch x5,x6,x8,x9 only.
.type .Lpoly_madd,%function
.Lpoly_madd:
    mul x5,x0,x1
    umulh x6,x0,x1
    adds x5,x5,x2
    adc x6,x6,xzr
    mul x8,x6,x4
    umulh x9,x6,x4
    adds x8,x8,x5
    adc x9,x9,xzr
    mul x9,x9,x4
    adds x8,x8,x9
    cset x9,cs
    mul x9,x9,x4
    add x8,x8,x9
    subs x9,x8,x3
    csel x0,x9,x8,hs
    ret
.size .Lpoly_madd,.-.Lpoly_madd

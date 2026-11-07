// Static untrusted analytic evidence;never host-execute or use credentials.
// x0=input192,x1=output32,x2=two canonical Fq keys32,q=2^127-1.
// Output full first127-bit residue in16 bytes,low64 of second,then8 zero bytes.
.text
.p2align 2
.macro POLY127_STEP al,ah,kl,kh
    mov x0,\al
    mov x1,\ah
    mov x2,\kl
    mov x3,\kh
    ldp x4,x5,[sp,#96]
    bl .Lpoly127_madd
    mov \al,x0
    mov \ah,x1
.endm
.macro POLY127_COEFF i
.if ((126*\i)%64) == 0
    ldp x4,x5,[x19,#(8*((126*\i)/64))]
.else
    ldr x4,[x19,#(8*((126*\i)/64))]
    ldr x5,[x19,#(8*((126*\i)/64+1))]
    ldr x6,[x19,#(8*((126*\i)/64+2))]
    lsr x4,x4,#((126*\i)%64)
    orr x4,x4,x5,lsl #(64-((126*\i)%64))
    lsr x5,x5,#((126*\i)%64)
    orr x5,x5,x6,lsl #(64-((126*\i)%64))
.endif
    and x5,x5,#0x3fffffffffffffff
    stp x4,x5,[sp,#96]
    POLY127_STEP x22,x23,x26,x27
    POLY127_STEP x24,x25,x28,x29
.endm
.global fapi_poly192
.type fapi_poly192,%function
fapi_poly192:
    stp x19,x20,[sp,#-112]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    stp x29,x30,[sp,#80]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    ldp x26,x27,[x21,#0]
    ldp x28,x29,[x21,#16]
    ldr x4,[x19,#184]
    lsr x22,x4,#40
    mov x23,xzr
    mov x24,x22
    mov x25,xzr
.irp i,11,10,9,8,7,6,5,4,3,2,1,0
    POLY127_COEFF \i
.endr
    stp x22,x23,[x20,#0]
    stp x24,xzr,[x20,#16]
    ldp x29,x30,[sp,#80]
    ldp x27,x28,[sp,#64]
    ldp x25,x26,[sp,#48]
    ldp x23,x24,[sp,#32]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#112
    ret
.size fapi_poly192,.-fapi_poly192

// Ordinary two-limb A,K<q,C<2^126 -> (A*K+C) modq in x0,x1.
// Complete leaf38 scalar instructions,caller scratch only,x18/callee state untouched.
.type .Lpoly127_madd,%function
.Lpoly127_madd:
    mul x6,x0,x2
    umulh x7,x0,x2
    mul x10,x0,x3
    umulh x11,x0,x3
    adds x7,x7,x10
    adc x11,x11,xzr
    mul x12,x1,x2
    umulh x13,x1,x2
    adds x7,x7,x12
    adc x13,x13,xzr
    mul x8,x1,x3
    umulh x9,x1,x3
    adds x8,x8,x11
    adc x9,x9,xzr
    adds x8,x8,x13
    adc x9,x9,xzr
    and x1,x7,#0x7fffffffffffffff
    lsr x10,x7,#63
    orr x10,x10,x8,lsl #1
    lsr x11,x8,#63
    orr x11,x11,x9,lsl #1
    adds x0,x6,x10
    adcs x1,x1,x11
    adc x12,xzr,xzr
    adds x0,x0,x4
    adcs x1,x1,x5
    adc x12,x12,xzr
    lsr x13,x1,#63
    orr x13,x13,x12,lsl #1
    and x1,x1,#0x7fffffffffffffff
    adds x0,x0,x13
    adc x1,x1,xzr
    adds x6,x0,#1
    orr x10,xzr,#0x7fffffffffffffff
    sbcs x7,x1,x10
    csel x0,x6,x0,hs
    csel x1,x7,x1,hs
    ret
.size .Lpoly127_madd,.-.Lpoly127_madd

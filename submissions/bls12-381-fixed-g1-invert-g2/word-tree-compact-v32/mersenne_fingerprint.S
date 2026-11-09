// Inert analytic evidence ONLY. Never host-execute or use credentials.
// x0 out32,x1 S192,x2 coeff648:rows0/1 interleaved27 pairs,
// row2@432;all canonical m61. Exact same affine family,three rows27.
.text
.p2align 2
.macro MH_TERM off,shift,index,reuse,zero,last
.if \last
    extr x12,xzr,x17,#\shift
    and x12,x12,#0xfffffffff
.else
.if \reuse
.else
    ldp x16,x17,[x1,#\off]
.endif
.if \zero
    and x12,x16,#0xfffffffffffffff
.else
    extr x12,x17,x16,#\shift
    and x12,x12,#0xfffffffffffffff
.endif
.endif
    ldp x13,x14,[x2,#(16*\index)]
    ldr x15,[x2,#(432+8*\index)]
    mul x10,x13,x12
    umulh x11,x13,x12
    adds x4,x4,x10
    adc x5,x5,x11
    mul x10,x14,x12
    umulh x11,x14,x12
    adds x6,x6,x10
    adc x7,x7,x11
    mul x10,x15,x12
    umulh x11,x15,x12
    adds x8,x8,x10
    adc x9,x9,x11
.endm
.macro MS_REDUCE low, high, dst, mask
    and x10,\low,#0x1fffffffffffffff
    extr x11,\high,\low,#61
    and x11,x11,#0x1fffffffffffffff
    lsr x12,\high,#58
    add x10,x10,x11
    add x10,x10,x12
    and x11,x10,#0x1fffffffffffffff
    add x11,x11,x10,lsr #61
    subs x12,x11,x3
    csel \dst,x12,x11,hs
    and \dst,\dst,#\mask
.endm
.global mersenne_fingerprint
.type mersenne_fingerprint,%function
mersenne_fingerprint:
    mov x3,#0x1fffffffffffffff
    ldp x4,x6,[x2]
    mov x5,xzr
    mov x7,xzr
    ldr x8,[x2,#432]
    mov x9,xzr
    MH_TERM 0,0,1,0,1,0
    MH_TERM 0,60,2,1,0,0
    MH_TERM 8,56,3,0,0,0
    MH_TERM 16,52,4,0,0,0
    MH_TERM 24,48,5,0,0,0
    MH_TERM 32,44,6,0,0,0
    MH_TERM 40,40,7,0,0,0
    MH_TERM 48,36,8,0,0,0
    MH_TERM 56,32,9,0,0,0
    MH_TERM 64,28,10,0,0,0
    MH_TERM 72,24,11,0,0,0
    MH_TERM 80,20,12,0,0,0
    MH_TERM 88,16,13,0,0,0
    MH_TERM 96,12,14,0,0,0
    MH_TERM 104,8,15,0,0,0
    MH_TERM 112,4,16,0,0,0
    MH_TERM 120,0,17,0,1,0
    MH_TERM 120,60,18,1,0,0
    MH_TERM 128,56,19,0,0,0
    MH_TERM 136,52,20,0,0,0
    MH_TERM 144,48,21,0,0,0
    MH_TERM 152,44,22,0,0,0
    MH_TERM 160,40,23,0,0,0
    MH_TERM 168,36,24,0,0,0
    MH_TERM 176,32,25,0,0,0
    MH_TERM 184,28,26,0,0,1
    MS_REDUCE x4,x5,x16,0x3fffffffffffff
    MS_REDUCE x6,x7,x17,0x3fffffffffffff
    MS_REDUCE x8,x9,x4,0xfffffffffffff
    orr x16,x16,x17,lsl #54
    lsr x17,x17,#10
    orr x17,x17,x4,lsl #44
    lsr x4,x4,#20
    stp x16,x17,[x0]
    stp x4,xzr,[x0,#16]
    ret
.size mersenne_fingerprint,.-mersenne_fingerprint

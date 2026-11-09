// Inert analytic evidence only. Never execute on host or with credentials.
// x0=out32,x1=S192,x2=81 canonical q61 coefficients,3 rows of27.
// Deferred three-lane128-bit accumulators; no target-dependent oracle.
.text
.p2align 2
.macro MS_TERM off, shift, coff, last
.if \last
    ldr x10,[x1,#\off]
    extr x12,xzr,x10,#\shift
    and x12,x12,#0xfffffffff
.else
    ldp x10,x11,[x1,#\off]
    extr x12,x11,x10,#\shift
    and x12,x12,#0xfffffffffffffff
.endif
    ldr x13,[x2,#\coff]
    mul x14,x13,x12
    umulh x15,x13,x12
    adds x4,x4,x14
    adc x5,x5,x15
    ldr x13,[x2,#(\coff+216)]
    mul x14,x13,x12
    umulh x15,x13,x12
    adds x6,x6,x14
    adc x7,x7,x15
    ldr x13,[x2,#(\coff+432)]
    mul x14,x13,x12
    umulh x15,x13,x12
    adds x8,x8,x14
    adc x9,x9,x15
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
    ldr x4,[x2]
    mov x5,xzr
    ldr x6,[x2,#216]
    mov x7,xzr
    ldr x8,[x2,#432]
    mov x9,xzr
    MS_TERM 0,0,8,0
    MS_TERM 0,60,16,0
    MS_TERM 8,56,24,0
    MS_TERM 16,52,32,0
    MS_TERM 24,48,40,0
    MS_TERM 32,44,48,0
    MS_TERM 40,40,56,0
    MS_TERM 48,36,64,0
    MS_TERM 56,32,72,0
    MS_TERM 64,28,80,0
    MS_TERM 72,24,88,0
    MS_TERM 80,20,96,0
    MS_TERM 88,16,104,0
    MS_TERM 96,12,112,0
    MS_TERM 104,8,120,0
    MS_TERM 112,4,128,0
    MS_TERM 120,0,136,0
    MS_TERM 120,60,144,0
    MS_TERM 128,56,152,0
    MS_TERM 136,52,160,0
    MS_TERM 144,48,168,0
    MS_TERM 152,44,176,0
    MS_TERM 160,40,184,0
    MS_TERM 168,36,192,0
    MS_TERM 176,32,200,0
    MS_TERM 184,28,208,1
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

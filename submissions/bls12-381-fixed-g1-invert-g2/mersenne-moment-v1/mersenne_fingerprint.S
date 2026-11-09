// Inert analytic evidence only; never execute on host or with credentials.
// x0=out32,x1=full S192,x2=81 canonical q61 coefficients,3 rows of27.
// q61=2^61-1. Three low-bit lanes54,54,52 are packed into first160 bits.
.text
.p2align 2
.macro MH_TERM off, shift, coff, last
    ldr x5,[x2,#\coff]
.if \last
    ldr x6,[x1,#\off]
    extr x8,xzr,x6,#\shift
    and x8,x8,#0xfffffffff
.else
    ldp x6,x7,[x1,#\off]
    extr x8,x7,x6,#\shift
    and x8,x8,#0xfffffffffffffff
.endif
    mul x9,x5,x8
    umulh x10,x5,x8
    and x11,x9,#0x1fffffffffffffff
    lsr x9,x9,#61
    add x11,x11,x9
    add x11,x11,x10,lsl #3
    subs x12,x11,x3
    csel x11,x12,x11,hs
    add x4,x4,x11
    subs x12,x4,x3
    csel x4,x12,x4,hs
.endm
.macro MH_ROW base, dst, mask
    ldr x4,[x2,#\base]
    MH_TERM 0,0,(\base+8),0
    MH_TERM 0,60,(\base+16),0
    MH_TERM 8,56,(\base+24),0
    MH_TERM 16,52,(\base+32),0
    MH_TERM 24,48,(\base+40),0
    MH_TERM 32,44,(\base+48),0
    MH_TERM 40,40,(\base+56),0
    MH_TERM 48,36,(\base+64),0
    MH_TERM 56,32,(\base+72),0
    MH_TERM 64,28,(\base+80),0
    MH_TERM 72,24,(\base+88),0
    MH_TERM 80,20,(\base+96),0
    MH_TERM 88,16,(\base+104),0
    MH_TERM 96,12,(\base+112),0
    MH_TERM 104,8,(\base+120),0
    MH_TERM 112,4,(\base+128),0
    MH_TERM 120,0,(\base+136),0
    MH_TERM 120,60,(\base+144),0
    MH_TERM 128,56,(\base+152),0
    MH_TERM 136,52,(\base+160),0
    MH_TERM 144,48,(\base+168),0
    MH_TERM 152,44,(\base+176),0
    MH_TERM 160,40,(\base+184),0
    MH_TERM 168,36,(\base+192),0
    MH_TERM 176,32,(\base+200),0
    MH_TERM 184,28,(\base+208),1
    and \dst,x4,#\mask
.endm
.global mersenne_fingerprint
.type mersenne_fingerprint,%function
mersenne_fingerprint:
    mov x3,#0x1fffffffffffffff
    MH_ROW 0,x14,0x3fffffffffffff
    MH_ROW 216,x15,0x3fffffffffffff
    MH_ROW 432,x16,0xfffffffffffff
    orr x14,x14,x15,lsl #54
    lsr x15,x15,#10
    orr x15,x15,x16,lsl #44
    lsr x16,x16,#20
    stp x14,x15,[x0]
    stp x16,xzr,[x0,#16]
    ret
.size mersenne_fingerprint,.-mersenne_fingerprint

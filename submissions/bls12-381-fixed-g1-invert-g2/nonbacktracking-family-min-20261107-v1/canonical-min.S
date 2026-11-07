// Static evidence only;neverhost-executed.
// x0:12 canonicalFp2candidates96bytes each,k0..11;selectedpairnonzero.
// x1:disjoint96-byteoutput,x2:64-bitkoutput. Realnumericthenimaginaryorder.
// Component limbs are canonicalMontgomery six64-bit little-endian words.
.text
.p2align 2
.global fapi_canonical_min
.type fapi_canonical_min,%function
.macro COMPARE_COMPONENT k,off
    add x1,x8,#(96*\k)
    mov x7,#\k
    ldp x2,x3,[x0,#(\off+0)]
    ldp x4,x5,[x1,#(\off+0)]
    subs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#(\off+16)]
    ldp x4,x5,[x1,#(\off+16)]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x0,#(\off+32)]
    ldp x4,x5,[x1,#(\off+32)]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    csel x0,x1,x0,hs
    csel x6,x7,x6,hs
.endm
.macro ZERO_REAL ptr,target
    ldp x2,x3,[\ptr,#0]
    orr x7,x2,x3
    ldp x2,x3,[\ptr,#16]
    orr x7,x7,x2
    orr x7,x7,x3
    ldp x2,x3,[\ptr,#32]
    orr x7,x7,x2
    orr x7,x7,x3
    cbnz x7,\target
.endm
fapi_canonical_min:
    mov x8,x0
    mov x9,x1
    mov x10,x2
    mov x6,xzr
    ZERO_REAL x0,.Lfirst_real_nonzero
// Evenfamilyrealcoordinatesallzero;otherfamilyrealcoordinatesarenonzero.
    COMPARE_COMPONENT 2,48
    COMPARE_COMPONENT 4,48
    COMPARE_COMPONENT 6,48
    COMPARE_COMPONENT 8,48
    COMPARE_COMPONENT 10,48
    b .Lminimum_output
.Lfirst_real_nonzero:
    add x1,x8,#96
    ZERO_REAL x1,.Lboth_real_nonzero
// Oddfamilyrealcoordinatesallzero.
    mov x0,x1
    mov x6,#1
    COMPARE_COMPONENT 3,48
    COMPARE_COMPONENT 5,48
    COMPARE_COMPONENT 7,48
    COMPARE_COMPONENT 9,48
    COMPARE_COMPONENT 11,48
    b .Lminimum_output
.Lboth_real_nonzero:
    COMPARE_COMPONENT 2,0
    COMPARE_COMPONENT 4,0
    COMPARE_COMPONENT 6,0
    COMPARE_COMPONENT 8,0
    COMPARE_COMPONENT 10,0
    mov x11,x0
    mov x12,x6
    add x0,x8,#96
    mov x6,#1
    COMPARE_COMPONENT 3,0
    COMPARE_COMPONENT 5,0
    COMPARE_COMPONENT 7,0
    COMPARE_COMPONENT 9,0
    COMPARE_COMPONENT 11,0
// Compareevenwinnerx11andoddwinnerx0. ORR/CBNZpreservefinalcarry.
    ldp x2,x3,[x11,#0]
    ldp x4,x5,[x0,#0]
    subs x2,x2,x4
    sbcs x3,x3,x5
    orr x7,x2,x3
    ldp x2,x3,[x11,#16]
    ldp x4,x5,[x0,#16]
    sbcs x2,x2,x4
    sbcs x3,x3,x5
    orr x2,x2,x3
    orr x7,x7,x2
    ldp x2,x3,[x11,#32]
    ldp x4,x5,[x0,#32]
    sbcs x2,x2,x4
    sbcs x3,x3,x5
    orr x2,x2,x3
    orr x7,x7,x2
    cbnz x7,.Lselect_family
// Onlyreal-tiedfamilywinnersneedimaginarycomparison.
    ldp x2,x3,[x11,#48]
    ldp x4,x5,[x0,#48]
    subs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x11,#64]
    ldp x4,x5,[x0,#64]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
    ldp x2,x3,[x11,#80]
    ldp x4,x5,[x0,#80]
    sbcs xzr,x2,x4
    sbcs xzr,x3,x5
.Lselect_family:
    csel x0,x0,x11,hs
    csel x6,x6,x12,hs
.Lminimum_output:
    ldp x2,x3,[x0,#0]
    stp x2,x3,[x9,#0]
    ldp x2,x3,[x0,#16]
    stp x2,x3,[x9,#16]
    ldp x2,x3,[x0,#32]
    stp x2,x3,[x9,#32]
    ldp x2,x3,[x0,#48]
    stp x2,x3,[x9,#48]
    ldp x2,x3,[x0,#64]
    stp x2,x3,[x9,#64]
    ldp x2,x3,[x0,#80]
    stp x2,x3,[x9,#80]
    str x6,[x10]
    mov x0,x6
    ret
.size fapi_canonical_min,.-fapi_canonical_min

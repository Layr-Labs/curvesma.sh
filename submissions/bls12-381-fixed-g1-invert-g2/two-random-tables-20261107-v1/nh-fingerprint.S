// Static evidence only;neverhost-executed. Same outer ABI as prior fingerprint.
// x0 canonicalS192,x1 seed-record1040,x2 output40. Scalar only,x18 untouched.
// Twelve consecutive pair blocks;each five independent16-byte NH seed pairs.
// Last80bytes are five independent odd128-bit multiplier keys.
.text
.p2align 2
.global fapi_tabulate
.type fapi_tabulate,%function
.macro NH_PAIR first=0
    ldp x13,x14,[x22],#16
    NH_LANE x3,x4,\first
    NH_LANE x5,x6,\first
    NH_LANE x7,x8,\first
    NH_LANE x9,x10,\first
    NH_LANE x11,x12,\first
.endm
.macro NH_LANE lo,hi,first
    ldp x15,x16,[x20],#16
    add x15,x15,x13
    add x16,x16,x14
.if \first
    mul \lo,x15,x16
    umulh \hi,x15,x16
.else
    mul x17,x15,x16
    umulh x19,x15,x16
    adds \lo,\lo,x17
    adc \hi,\hi,x19
.endif
.endm
.macro COMPRESS lo,hi,offset
    ldp x13,x14,[x20],#16
    umulh x15,\lo,x13
    madd x15,\hi,x13,x15
    madd x15,\lo,x14,x15
    rev x15,x15
    str x15,[x21,#\offset]
.endm
fapi_tabulate:
    stp x19,x20,[sp,#-32]!
    stp x21,x22,[sp,#16]
    mov x20,x1
    mov x21,x2
    mov x22,x0
    NH_PAIR 1
.rept 11
    NH_PAIR 0
.endr
    COMPRESS x3,x4,0
    COMPRESS x5,x6,8
    COMPRESS x7,x8,16
    COMPRESS x9,x10,24
    COMPRESS x11,x12,32
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#32
    ret
.size fapi_tabulate,.-fapi_tabulate

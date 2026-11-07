// Exactkey is36bwords+constantzero;T72[0] foldedintoinitialaccumulators.
// Exact paired simple tabulation;compile-time macros introduce no runtime work.
// x0 key,x1 array73 bases,x2 output40 bytes;caller registers0..16 only.
.text
.p2align 2
.global fapi_tabulate
.type fapi_tabulate,%function
.macro TABLE_OP digit,base
    add x10,\digit,\digit,lsl #2
    add x10,\base,x10,lsl #3
    ldp x11,x12,[x10]
    eor x3,x3,x11
    eor x4,x4,x12
    ldp x11,x12,[x10,#16]
    eor x5,x5,x11
    eor x6,x6,x12
    ldr x11,[x10,#32]
    eor x7,x7,x11
.endm
.macro PAIR digitbase,digitoffset,tablebase,tableoffset
    ldp w8,w13,[\digitbase,#\digitoffset]
    ldp x9,x14,[\tablebase,#\tableoffset]
    TABLE_OP x8,x9
    TABLE_OP x13,x14
.endm
fapi_tabulate:
    add x15,x0,#256
    add x16,x1,#512
    ldr x9,[x1,#576]
    ldp x3,x4,[x9]
    ldp x5,x6,[x9,#16]
    ldr x7,[x9,#32]
    PAIR x0,0,x1,0
    PAIR x0,8,x1,16
    PAIR x0,16,x1,32
    PAIR x0,24,x1,48
    PAIR x0,32,x1,64
    PAIR x0,40,x1,80
    PAIR x0,48,x1,96
    PAIR x0,56,x1,112
    PAIR x0,64,x1,128
    PAIR x0,72,x1,144
    PAIR x0,80,x1,160
    PAIR x0,88,x1,176
    PAIR x0,96,x1,192
    PAIR x0,104,x1,208
    PAIR x0,112,x1,224
    PAIR x0,120,x1,240
    PAIR x0,128,x1,256
    PAIR x0,136,x1,272
    PAIR x0,144,x1,288
    PAIR x0,152,x1,304
    PAIR x0,160,x1,320
    PAIR x0,168,x1,336
    PAIR x0,176,x1,352
    PAIR x0,184,x1,368
    PAIR x0,192,x1,384
    PAIR x0,200,x1,400
    PAIR x0,208,x1,416
    PAIR x0,216,x1,432
    PAIR x0,224,x1,448
    PAIR x0,232,x1,464
    PAIR x0,240,x1,480
    PAIR x0,248,x1,496
    PAIR x15,0,x16,0
    PAIR x15,8,x16,16
    PAIR x15,16,x16,32
    PAIR x15,24,x16,48
    rev x3,x3
    rev x4,x4
    rev x5,x5
    rev x6,x6
    rev x7,x7
    stp x3,x4,[x2]
    stp x5,x6,[x2,#16]
    str x7,[x2,#32]
    ret
.size fapi_tabulate,.-fapi_tabulate

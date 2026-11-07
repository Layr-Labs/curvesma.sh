// Static untrusted evidence;never host-execute or use credentials.
// x0=original GT576,x1=scalarBE32,x2=disjoint record2048.
// Three Fp4 cubic forms at infinity,+1,-1,+2,-2,+3,-3;21 native Fp values1008.
// Scalar at1536;all other record bytes initialized zero.
.text
.p2align 2
.macro PRE_OP op,ob,oo,ab,ao,bb,bo
    add x0,\ob,#\oo
    add x1,\ab,#\ao
    add x2,\bb,#\bo
    mov x3,x22
    bl \op
.endm
.global fapi_trace_pretransform
.type fapi_trace_pretransform,%function
fapi_trace_pretransform:
    sub sp,sp,#480
    stp x19,x20,[sp,#0]
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    stp x29,x30,[sp,#80]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    adrp x22,.Lpre_p
    add x22,x22,:lo12:.Lpre_p
    mov x24,x21
.rept 128
    stp xzr,xzr,[x24],#16
.endr
    // Fp4 form0;ar0,ai48,br384,bi432.
    add x24,x21,#0
    add x25,x19,#432
    ldp x4,x5,[x25,#0]
    stp x4,x5,[x24,#0]
    ldp x4,x5,[x25,#16]
    stp x4,x5,[x24,#16]
    ldp x4,x5,[x25,#32]
    stp x4,x5,[x24,#32]
    PRE_OP add_mod_384,sp,96,x19,384,x19,384
    PRE_OP add_mod_384,sp,144,sp,96,x19,384
    PRE_OP add_mod_384,sp,384,x19,48,x19,48
    PRE_OP add_mod_384,sp,192,sp,384,x19,48
    PRE_OP add_mod_384,sp,432,sp,384,sp,384
    PRE_OP add_mod_384,sp,240,sp,432,sp,432
    PRE_OP add_mod_384,sp,384,x19,432,x19,432
    PRE_OP add_mod_384,sp,432,sp,384,x19,432
    PRE_OP add_mod_384,sp,288,sp,432,sp,432
    PRE_OP add_mod_384,sp,384,sp,288,sp,288
    PRE_OP add_mod_384,sp,336,sp,384,sp,384
    PRE_OP add_mod_384,x21,48,x19,0,x19,384
    PRE_OP sub_mod_384,x21,96,x19,0,x19,384
    PRE_OP add_mod_384,sp,384,x19,0,sp,192
    PRE_OP add_mod_384,x21,144,sp,384,sp,96
    PRE_OP add_mod_384,x21,144,x21,144,sp,288
    PRE_OP sub_mod_384,x21,192,sp,384,sp,96
    PRE_OP sub_mod_384,x21,192,x21,192,sp,288
    PRE_OP add_mod_384,sp,384,x19,0,sp,240
    PRE_OP add_mod_384,x21,240,sp,384,sp,144
    PRE_OP add_mod_384,x21,240,x21,240,sp,336
    PRE_OP sub_mod_384,x21,288,sp,384,sp,144
    PRE_OP sub_mod_384,x21,288,x21,288,sp,336
    // Fp4 form1;ar288,ai336,br192,bi240.
    add x24,x21,#336
    add x25,x19,#240
    ldp x4,x5,[x25,#0]
    stp x4,x5,[x24,#0]
    ldp x4,x5,[x25,#16]
    stp x4,x5,[x24,#16]
    ldp x4,x5,[x25,#32]
    stp x4,x5,[x24,#32]
    PRE_OP add_mod_384,sp,96,x19,192,x19,192
    PRE_OP add_mod_384,sp,144,sp,96,x19,192
    PRE_OP add_mod_384,sp,384,x19,336,x19,336
    PRE_OP add_mod_384,sp,192,sp,384,x19,336
    PRE_OP add_mod_384,sp,432,sp,384,sp,384
    PRE_OP add_mod_384,sp,240,sp,432,sp,432
    PRE_OP add_mod_384,sp,384,x19,240,x19,240
    PRE_OP add_mod_384,sp,432,sp,384,x19,240
    PRE_OP add_mod_384,sp,288,sp,432,sp,432
    PRE_OP add_mod_384,sp,384,sp,288,sp,288
    PRE_OP add_mod_384,sp,336,sp,384,sp,384
    PRE_OP add_mod_384,x21,384,x19,288,x19,192
    PRE_OP sub_mod_384,x21,432,x19,288,x19,192
    PRE_OP add_mod_384,sp,384,x19,288,sp,192
    PRE_OP add_mod_384,x21,480,sp,384,sp,96
    PRE_OP add_mod_384,x21,480,x21,480,sp,288
    PRE_OP sub_mod_384,x21,528,sp,384,sp,96
    PRE_OP sub_mod_384,x21,528,x21,528,sp,288
    PRE_OP add_mod_384,sp,384,x19,288,sp,240
    PRE_OP add_mod_384,x21,576,sp,384,sp,144
    PRE_OP add_mod_384,x21,576,x21,576,sp,336
    PRE_OP sub_mod_384,x21,624,sp,384,sp,144
    PRE_OP sub_mod_384,x21,624,x21,624,sp,336
    // Fp4 form2;ar96,ai144,br480,bi528.
    add x24,x21,#672
    add x25,x19,#528
    ldp x4,x5,[x25,#0]
    stp x4,x5,[x24,#0]
    ldp x4,x5,[x25,#16]
    stp x4,x5,[x24,#16]
    ldp x4,x5,[x25,#32]
    stp x4,x5,[x24,#32]
    PRE_OP add_mod_384,sp,96,x19,480,x19,480
    PRE_OP add_mod_384,sp,144,sp,96,x19,480
    PRE_OP add_mod_384,sp,384,x19,144,x19,144
    PRE_OP add_mod_384,sp,192,sp,384,x19,144
    PRE_OP add_mod_384,sp,432,sp,384,sp,384
    PRE_OP add_mod_384,sp,240,sp,432,sp,432
    PRE_OP add_mod_384,sp,384,x19,528,x19,528
    PRE_OP add_mod_384,sp,432,sp,384,x19,528
    PRE_OP add_mod_384,sp,288,sp,432,sp,432
    PRE_OP add_mod_384,sp,384,sp,288,sp,288
    PRE_OP add_mod_384,sp,336,sp,384,sp,384
    PRE_OP add_mod_384,x21,720,x19,96,x19,480
    PRE_OP sub_mod_384,x21,768,x19,96,x19,480
    PRE_OP add_mod_384,sp,384,x19,96,sp,192
    PRE_OP add_mod_384,x21,816,sp,384,sp,96
    PRE_OP add_mod_384,x21,816,x21,816,sp,288
    PRE_OP sub_mod_384,x21,864,sp,384,sp,96
    PRE_OP sub_mod_384,x21,864,x21,864,sp,288
    PRE_OP add_mod_384,sp,384,x19,96,sp,240
    PRE_OP add_mod_384,x21,912,sp,384,sp,144
    PRE_OP add_mod_384,x21,912,x21,912,sp,336
    PRE_OP sub_mod_384,x21,960,sp,384,sp,144
    PRE_OP sub_mod_384,x21,960,x21,960,sp,336
    add x24,x21,#1536
    ldp x4,x5,[x20,#0]
    stp x4,x5,[x24,#0]
    ldp x4,x5,[x20,#16]
    stp x4,x5,[x24,#16]
    ldp x29,x30,[sp,#80]
    ldp x27,x28,[sp,#64]
    ldp x25,x26,[sp,#48]
    ldp x23,x24,[sp,#32]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp,#0]
    add sp,sp,#480
    ret
.size fapi_trace_pretransform,.-fapi_trace_pretransform
.section .rodata
.p2align 3
.Lpre_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a

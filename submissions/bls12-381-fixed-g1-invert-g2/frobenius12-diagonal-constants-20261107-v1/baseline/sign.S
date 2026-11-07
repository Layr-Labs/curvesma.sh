// x0 points to a (six canonical48-byte Montgomery coefficients),x1 sign output.
// Forward-only scalar leaf, preserves every callee-saved register.
.text
.p2align 2
.global fapi_sign
.type fapi_sign,%function
fapi_sign:
.set sign_i,0
.rept 6
    ldr x2,[x0,#(48*sign_i)]
    ldr x3,[x0,#(48*sign_i+8)]
    ldr x4,[x0,#(48*sign_i+16)]
    ldr x5,[x0,#(48*sign_i+24)]
    ldr x6,[x0,#(48*sign_i+32)]
    ldr x7,[x0,#(48*sign_i+40)]
    orr x8,x2,x3
    orr x8,x8,x4
    orr x8,x8,x5
    orr x8,x8,x6
    orr x8,x8,x7
    cbnz x8,.Lselected
.set sign_i,sign_i+1
.endr
    mov x8,xzr
    str x8,[x1]
    ret
.Lselected:
    adrp x9,fapi_half_p
    add x9,x9,:lo12:fapi_half_p
    ldp x10,x11,[x9]
    ldp x12,x13,[x9,#16]
    ldp x14,x15,[x9,#32]
    subs xzr,x2,x10
    sbcs xzr,x3,x11
    sbcs xzr,x4,x12
    sbcs xzr,x5,x13
    sbcs xzr,x6,x14
    sbcs xzr,x7,x15
    cset x8,hs
    str x8,[x1]
    ret
.size fapi_sign,.-fapi_sign
.section .rodata
.p2align 3
fapi_half_p:
    .quad 0xdcff7fffffffd556
    .quad 0x0f55ffff58a9ffff
    .quad 0xb39869507b587b12
    .quad 0xb23ba5c279c2895f
    .quad 0x258dd3db21a5d66b
    .quad 0x0d0088f51cbff34d

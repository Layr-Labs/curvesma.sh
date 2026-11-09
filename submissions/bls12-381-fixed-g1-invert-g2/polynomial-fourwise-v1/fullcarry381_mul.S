// Static evidence; do not execute on host or with credentials.
// Full-carry Montgomery scheme and accounting derive from e750ade2; see proof attribution.
.text
.p2align 2
.macro ADD_FULL k
.irp slot,0,1,2,3,4,5,6,7,8,9,10,11,12
.if \slot >= (\k)
    ldr x7,[sp,#(8*\slot)]
.if \slot == (\k)
    adds x7,x7,x5
.elseif \slot == ((\k)+1)
    adcs x7,x7,x6
.else
    adcs x7,x7,xzr
.endif
    str x7,[sp,#(8*\slot)]
.endif
.endr
.endm
.global fullcarry381_mul
.type fullcarry381_mul,%function
fullcarry381_mul:
    sub sp,sp,#112
    mov x15,x2
    adrp x14,.Lmod381
    add x14,x14,:lo12:.Lmod381
    adrp x13,.Ln0381
    ldr x13,[x13,:lo12:.Ln0381]
.irp slot,0,1,2,3,4,5,6,7,8,9,10,11,12
    str xzr,[sp,#(8*\slot)]
.endr
.irp i,0,1,2,3,4,5
.irp j,0,1,2,3,4,5
    ldr x3,[x0,#(8*\i)]
    ldr x4,[x1,#(8*\j)]
    mul x5,x3,x4
    umulh x6,x3,x4
    ADD_FULL (\i+\j)
.endr
.endr
.irp i,0,1,2,3,4,5
    ldr x9,[sp,#(8*\i)]
    mul x9,x9,x13
.irp j,0,1,2,3,4,5
    ldr x4,[x14,#(8*\j)]
    mul x5,x9,x4
    umulh x6,x9,x4
    ADD_FULL (\i+\j)
.endr
.endr
    ldp x0,x1,[sp,#48]
    ldp x2,x3,[sp,#64]
    ldp x4,x5,[sp,#80]
    ldp x6,x7,[x14]
    ldp x8,x9,[x14,#16]
    ldp x10,x11,[x14,#32]
    subs x6,x0,x6
    sbcs x7,x1,x7
    sbcs x8,x2,x8
    sbcs x9,x3,x9
    sbcs x10,x4,x10
    sbcs x11,x5,x11
    csel x0,x6,x0,hs
    csel x1,x7,x1,hs
    csel x2,x8,x2,hs
    csel x3,x9,x3,hs
    csel x4,x10,x4,hs
    csel x5,x11,x5,hs
    stp x0,x1,[x15]
    stp x2,x3,[x15,#16]
    stp x4,x5,[x15,#32]
    add sp,sp,#112
    ret
.size fullcarry381_mul,.-fullcarry381_mul
.section .rodata
.p2align 3
.Lmod381:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624,0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
.Ln0381:
    .quad 0x89f3fffcfffcfffd

// Static analytic evidence only. Never compile/execute on the host or with credentials.
// x0=a48,x1=b48,x2=out48; canonical native Montgomery inputs modulo full p.
// Uses caller scratch x0..x15 only; x18/callee registers/LR preserved.
.text
.p2align 2
.macro COARSE_ADD_AT k
.irp limb,0,1,2,3,4,5,6,7,8,9,10,11,12
.if \limb >= (\k)
    ldr x7,[sp,#(8*\limb)]
.if \limb == (\k)
    adds x7,x7,x5
.elseif \limb == ((\k)+1)
    adcs x7,x7,x6
.else
    adcs x7,x7,xzr
.endif
    str x7,[sp,#(8*\limb)]
.endif
.endr
.endm
.global coarse_fp_mul
.type coarse_fp_mul,%function
coarse_fp_mul:
    sub sp,sp,#112
    mov x15,x2
    adrp x14,.Lcoarse_p
    add x14,x14,:lo12:.Lcoarse_p
    adrp x13,.Lcoarse_n0
    ldr x13,[x13,:lo12:.Lcoarse_n0]
.irp limb,0,1,2,3,4,5,6,7,8,9,10,11,12
    str xzr,[sp,#(8*\limb)]
.endr
// All i/j loops are assembler expansion over public constants,no runtime loop.
.irp i,0,1,2,3,4,5
.irp j,0,1,2,3,4,5
    ldr x3,[x0,#(8*\i)]
    ldr x4,[x1,#(8*\j)]
    mul x5,x3,x4
    umulh x6,x3,x4
    COARSE_ADD_AT (\i+\j)
.endr
.endr
.irp i,0,1,2,3,4,5
    ldr x9,[sp,#(8*\i)]
    mul x9,x9,x13
.irp j,0,1,2,3,4,5
    ldr x4,[x14,#(8*\j)]
    mul x5,x9,x4
    umulh x6,x9,x4
    COARSE_ADD_AT (\i+\j)
.endr
.endr
    ldp x0,x1,[sp,#48]
    ldp x2,x3,[sp,#64]
    ldp x4,x5,[sp,#80]
    ldp x6,x7,[x14,#0]
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
    stp x0,x1,[x15,#0]
    stp x2,x3,[x15,#16]
    stp x4,x5,[x15,#32]
    add sp,sp,#112
    ret
.size coarse_fp_mul,.-coarse_fp_mul
.section .rodata
.p2align 3
.Lcoarse_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
.Lcoarse_n0:
    .quad 0x89f3fffcfffcfffd

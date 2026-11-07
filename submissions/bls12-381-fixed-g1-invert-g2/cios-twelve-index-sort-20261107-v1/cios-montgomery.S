// Static untrusted evidence;never host-compile/execute or use credentials.
// x0=a48,x1=b48,x2=out48;full canonical native Montgomery values modulo p.
// CIOS six limbs;uses volatile x0..x15;preservesx18/callee/LR;stack64 aligned.
// Public protected p<2^383 is essential to the proved seven-limb tail bound.
.text
.p2align 2
.macro CIOS_MAC
    mul x5,x8,x3
    umulh x6,x8,x3
    adds x5,x5,x4
    adc x6,x6,xzr
    adds x5,x5,x9
    adc x9,x6,xzr
.endm
.global cios_fp_mul
.type cios_fp_mul,%function
cios_fp_mul:
    sub sp,sp,#64
    mov x15,x2
    adrp x14,.Lcios_p
    add x14,x14,:lo12:.Lcios_p
    adrp x13,.Lcios_n0
    ldr x13,[x13,:lo12:.Lcios_n0]
.irp limb,0,1,2,3,4,5,6
    str xzr,[sp,#(8*\limb)]
.endr
.irp i,0,1,2,3,4,5
    ldr x8,[x0,#(8*\i)]
    mov x9,xzr
.irp j,0,1,2,3,4,5
    ldr x3,[x1,#(8*\j)]
    ldr x4,[sp,#(8*\j)]
    CIOS_MAC
    str x5,[sp,#(8*\j)]
.endr
    str x9,[sp,#48]
    ldr x8,[sp,#0]
    mul x8,x8,x13
    mov x9,xzr
.irp j,0,1,2,3,4,5
    ldr x3,[x14,#(8*\j)]
    ldr x4,[sp,#(8*\j)]
    CIOS_MAC
    str x5,[sp,#(8*\j)]
.endr
    ldr x7,[sp,#48]
    add x7,x7,x9
    str x7,[sp,#48]
    ldp x3,x4,[sp,#8]
    ldp x5,x6,[sp,#24]
    ldp x7,x8,[sp,#40]
    stp x3,x4,[sp,#0]
    stp x5,x6,[sp,#16]
    stp x7,x8,[sp,#32]
    str xzr,[sp,#48]
.endr
    ldp x0,x1,[sp,#0]
    ldp x2,x3,[sp,#16]
    ldp x4,x5,[sp,#32]
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
    add sp,sp,#64
    ret
.size cios_fp_mul,.-cios_fp_mul
.section .rodata
.p2align 3
.Lcios_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
.Lcios_n0:
    .quad 0x89f3fffcfffcfffd

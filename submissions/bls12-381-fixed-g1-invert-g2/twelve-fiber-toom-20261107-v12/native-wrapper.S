// Static evidence only. Never execute on host or with credentials.
// x0=a48,x1=b48,x2=out48; canonical six-limb Montgomery operands.
// Invoke pinned blst mul_mont_384(out,a,b,p,n0), preserving x29/LR.
.text
.p2align 2
.global coarse_fp_mul
.type coarse_fp_mul,%function
coarse_fp_mul:
    stp x29,x30,[sp,#-16]!
    mov x29,sp
    mov x4,x0
    mov x0,x2
    mov x2,x1
    mov x1,x4
    adrp x3,.Lnative_p
    add x3,x3,:lo12:.Lnative_p
    adrp x4,.Lnative_n0
    ldr x4,[x4,:lo12:.Lnative_n0]
    bl mul_mont_384
    ldp x29,x30,[sp],#16
    ret
.size coarse_fp_mul,.-coarse_fp_mul
.section .rodata
.p2align 3
.Lnative_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
.Lnative_n0:
    .quad 0x89f3fffcfffcfffd

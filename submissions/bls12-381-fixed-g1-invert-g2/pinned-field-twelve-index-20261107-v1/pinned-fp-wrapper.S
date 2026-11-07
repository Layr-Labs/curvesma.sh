// Original static wrapper evidence;never host-compile/execute or use credentials.
// x0=a48,x1=b48,x2=out48;full canonical native Montgomery modulo protected p.
// Calls the verbatim pinned blst backend with its actual native ABI.
// Static link local text within BL range,no PLT/veneer;ordinary64-bit pointers.
.text
.p2align 2
.global pinned_fp_mul
.type pinned_fp_mul,%function
pinned_fp_mul:
    stp x19,x30,[sp,#-16]!
    mov x19,x0
    mov x0,x2
    mov x2,x1
    mov x1,x19
    adrp x3,.Lpinned_p
    add x3,x3,:lo12:.Lpinned_p
    adrp x4,.Lpinned_n0
    ldr x4,[x4,:lo12:.Lpinned_n0]
    bl mul_mont_384
    ldp x19,x30,[sp],#16
    ret
.size pinned_fp_mul,.-pinned_fp_mul
.section .rodata
.p2align 3
.Lpinned_p:
    .quad 0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a
.Lpinned_n0:
    .quad 0x89f3fffcfffcfffd

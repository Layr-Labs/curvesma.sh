// Static evidence only;neverhost-executed. OrdinaryAArch64pointerwidth8.
// x0output,x1inputa,x2inputb;tailcallthepinnedstaticblstMontgomerycore.
.text
.p2align 2
.global fapi_native_fp_mul
.type fapi_native_fp_mul,%function
fapi_native_fp_mul:
    adrp x3,fapi_native_mod_p
    add x3,x3,:lo12:fapi_native_mod_p
    adrp x4,fapi_native_n0
    ldr x4,[x4,:lo12:fapi_native_n0]
    b mul_mont_384
.size fapi_native_fp_mul,.-fapi_native_fp_mul
.section .rodata
.p2align 3
fapi_native_mod_p:
    .quad 0xb9feffffffffaaab
    .quad 0x1eabfffeb153ffff
    .quad 0x6730d2a0f6b0f624
    .quad 0x64774b84f38512bf
    .quad 0x4b1ba7b6434bacd7
    .quad 0x1a0111ea397fe69a
fapi_native_n0:
    .quad 0x89f3fffcfffcfffd

// Inert analytic source; ONLY credential-free shared llm-verifier may execute.
// x0 out48,x1 canonical native input48,x2 LUTbase672TiB,x3 const3728.
// Tables10: widths39+38*9,132TiB total; full381 inputbits.
// const prefix520 has p@448,Dceil@512;50 ordinary q*p slots64bytes@528,q0..49.
// Exact output/input alias or disjoint; tables/constants/frame disjoint. Frame96.
// SumY<10p impliesqtilde<=9; tableproduct seven limbs pluspad, not canonical Fp.
.text
.p2align 2
.global fixed_s_lookup
.type fixed_s_lookup,%function
fixed_s_lookup:
    stp x29,x30,[sp,#-96]!
    mov x29,sp
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    mov x25,x3
    mov x30,x1
    mov x26,x2
    mov x4,x0
    movz x27,#0x1800,lsl #32
    mov x28,#48
    ldp x19,x20,[x30,#0]
    ldp x21,x22,[x30,#16]
    ldp x23,x24,[x30,#32]
    mov x5,xzr
    mov x6,xzr
    mov x7,xzr
    mov x8,xzr
    mov x9,xzr
    mov x10,xzr
    mov x11,xzr
    ubfx x1,x19,#0,#39
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27
    extr x1,x20,x19,#39
    and x1,x1,#0x3fffffffff
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    ubfx x1,x20,#13,#38
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    extr x1,x21,x20,#51
    and x1,x1,#0x3fffffffff
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    ubfx x1,x21,#25,#38
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    extr x1,x22,x21,#63
    and x1,x1,#0x3fffffffff
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    extr x1,x23,x22,#37
    and x1,x1,#0x3fffffffff
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    ubfx x1,x23,#11,#38
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    extr x1,x24,x23,#49
    and x1,x1,#0x3fffffffff
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    add x26,x26,x27,lsr #1
    ubfx x1,x24,#23,#38
    madd x2,x1,x28,x26
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    adds x5,x5,x12
    adcs x6,x6,x13
    adcs x7,x7,x14
    adcs x8,x8,x15
    adcs x9,x9,x16
    adcs x10,x10,x17
    adc x11,x11,xzr
    ldr x2,[x25,#512]
    extr x1,x11,x10,#5
    udiv x1,x1,x2
    add x2,x25,#528
    add x2,x2,x1,lsl #6
    ldp x12,x13,[x2,#0]
    ldp x14,x15,[x2,#16]
    ldp x16,x17,[x2,#32]
    ldr x26,[x2,#48]
    subs x5,x5,x12
    sbcs x6,x6,x13
    sbcs x7,x7,x14
    sbcs x8,x8,x15
    sbcs x9,x9,x16
    sbcs x10,x10,x17
    sbcs x11,x11,x26
    ldp x19,x20,[x25,#448]
    ldp x21,x22,[x25,#464]
    ldp x23,x24,[x25,#480]
    subs x19,x5,x19
    sbcs x20,x6,x20
    sbcs x21,x7,x21
    sbcs x22,x8,x22
    sbcs x23,x9,x23
    sbcs x24,x10,x24
    b.lo .Lfixed_original
    stp x19,x20,[x4,#0]
    stp x21,x22,[x4,#16]
    stp x23,x24,[x4,#32]
    b .Lfixed_done
.Lfixed_original:
    stp x5,x6,[x4,#0]
    stp x7,x8,[x4,#16]
    stp x9,x10,[x4,#32]
.Lfixed_done:
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp],#96
    ret
.size fixed_s_lookup,.-fixed_s_lookup

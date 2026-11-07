// Static analytic evidence. Never host-compile/execute or use credentials.
// x0=input192,x1=output32,x2=three canonical little-endian Fq keys.
// q=2^64-59. Full192-bit fingerprint in first24 output bytes; last8 zero.
.text
.p2align 2
.global fapi_poly192
.type fapi_poly192,%function
fapi_poly192:
    stp x19,x20,[sp,#-96]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    stp x29,x30,[sp,#80]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    ldp x25,x26,[x21,#0]
    ldr x27,[x21,#16]
    mov x3,#-59
    mov x4,#59
    ldr x10,[x19,#184]
    lsr x22,x10,#40
    mov x23,x22
    mov x24,x22
    // Horner coefficient 23: bits 1449..1511.
    ldr x10,[x19,#176]
    lsr x10,x10,#41
    ldr x11,[x19,#184]
    orr x10,x10,x11,lsl #23
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 22: bits 1386..1448.
    ldr x10,[x19,#168]
    lsr x10,x10,#42
    ldr x11,[x19,#176]
    orr x10,x10,x11,lsl #22
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 21: bits 1323..1385.
    ldr x10,[x19,#160]
    lsr x10,x10,#43
    ldr x11,[x19,#168]
    orr x10,x10,x11,lsl #21
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 20: bits 1260..1322.
    ldr x10,[x19,#152]
    lsr x10,x10,#44
    ldr x11,[x19,#160]
    orr x10,x10,x11,lsl #20
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 19: bits 1197..1259.
    ldr x10,[x19,#144]
    lsr x10,x10,#45
    ldr x11,[x19,#152]
    orr x10,x10,x11,lsl #19
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 18: bits 1134..1196.
    ldr x10,[x19,#136]
    lsr x10,x10,#46
    ldr x11,[x19,#144]
    orr x10,x10,x11,lsl #18
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 17: bits 1071..1133.
    ldr x10,[x19,#128]
    lsr x10,x10,#47
    ldr x11,[x19,#136]
    orr x10,x10,x11,lsl #17
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 16: bits 1008..1070.
    ldr x10,[x19,#120]
    lsr x10,x10,#48
    ldr x11,[x19,#128]
    orr x10,x10,x11,lsl #16
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 15: bits 945..1007.
    ldr x10,[x19,#112]
    lsr x10,x10,#49
    ldr x11,[x19,#120]
    orr x10,x10,x11,lsl #15
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 14: bits 882..944.
    ldr x10,[x19,#104]
    lsr x10,x10,#50
    ldr x11,[x19,#112]
    orr x10,x10,x11,lsl #14
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 13: bits 819..881.
    ldr x10,[x19,#96]
    lsr x10,x10,#51
    ldr x11,[x19,#104]
    orr x10,x10,x11,lsl #13
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 12: bits 756..818.
    ldr x10,[x19,#88]
    lsr x10,x10,#52
    ldr x11,[x19,#96]
    orr x10,x10,x11,lsl #12
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 11: bits 693..755.
    ldr x10,[x19,#80]
    lsr x10,x10,#53
    ldr x11,[x19,#88]
    orr x10,x10,x11,lsl #11
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 10: bits 630..692.
    ldr x10,[x19,#72]
    lsr x10,x10,#54
    ldr x11,[x19,#80]
    orr x10,x10,x11,lsl #10
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 9: bits 567..629.
    ldr x10,[x19,#64]
    lsr x10,x10,#55
    ldr x11,[x19,#72]
    orr x10,x10,x11,lsl #9
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 8: bits 504..566.
    ldr x10,[x19,#56]
    lsr x10,x10,#56
    ldr x11,[x19,#64]
    orr x10,x10,x11,lsl #8
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 7: bits 441..503.
    ldr x10,[x19,#48]
    lsr x10,x10,#57
    ldr x11,[x19,#56]
    orr x10,x10,x11,lsl #7
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 6: bits 378..440.
    ldr x10,[x19,#40]
    lsr x10,x10,#58
    ldr x11,[x19,#48]
    orr x10,x10,x11,lsl #6
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 5: bits 315..377.
    ldr x10,[x19,#32]
    lsr x10,x10,#59
    ldr x11,[x19,#40]
    orr x10,x10,x11,lsl #5
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 4: bits 252..314.
    ldr x10,[x19,#24]
    lsr x10,x10,#60
    ldr x11,[x19,#32]
    orr x10,x10,x11,lsl #4
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 3: bits 189..251.
    ldr x10,[x19,#16]
    lsr x10,x10,#61
    ldr x11,[x19,#24]
    orr x10,x10,x11,lsl #3
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 2: bits 126..188.
    ldr x10,[x19,#8]
    lsr x10,x10,#62
    ldr x11,[x19,#16]
    orr x10,x10,x11,lsl #2
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 1: bits 63..125.
    ldr x10,[x19,#0]
    lsr x10,x10,#63
    ldr x11,[x19,#8]
    orr x10,x10,x11,lsl #1
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    // Horner coefficient 0: bits 0..62.
    ldr x10,[x19,#0]
    and x28,x10,#0x7fffffffffffffff
    mov x0,x22
    mov x1,x25
    mov x2,x28
    bl .Lpoly_madd
    mov x22,x0
    mov x0,x23
    mov x1,x26
    mov x2,x28
    bl .Lpoly_madd
    mov x23,x0
    mov x0,x24
    mov x1,x27
    mov x2,x28
    bl .Lpoly_madd
    mov x24,x0
    stp x22,x23,[x20,#0]
    stp x24,xzr,[x20,#16]
    ldp x29,x30,[sp,#80]
    ldp x27,x28,[sp,#64]
    ldp x25,x26,[sp,#48]
    ldp x23,x24,[sp,#32]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#96
    ret
.size fapi_poly192,.-fapi_poly192

// x0=A<q,x1=K<q,x2=c<2^63,x3=q,x4=59,return x0=(A*K+c) modq.
// x3/x4 and all callee state are preserved; scratch x5,x6,x8,x9 only.
.type .Lpoly_madd,%function
.Lpoly_madd:
    mul x5,x0,x1
    umulh x6,x0,x1
    adds x5,x5,x2
    adc x6,x6,xzr
    mul x8,x6,x4
    umulh x9,x6,x4
    adds x8,x8,x5
    adc x9,x9,xzr
    mul x9,x9,x4
    adds x8,x8,x9
    cset x9,cs
    mul x9,x9,x4
    add x8,x8,x9
    subs x9,x8,x3
    csel x0,x9,x8,hs
    ret
.size .Lpoly_madd,.-.Lpoly_madd

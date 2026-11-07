// Standard BLAKE3 final WALK block, len40,flags10,counter0,output32 bytes.
// Exact full7-round partial evaluation with ten message words resident.
// x0 mutable32-byte prefix CV/output,x1 initialized64-byte final block.
// No host execution. x18 untouched, all x19..x30 saved/restored.
.text
.p2align 2
.global fapi_blake3_final
.type fapi_blake3_final,%function
.macro G a,b,c,d,mx,my,xused,yused,firstzero
    add \a,\a,\b
.if \xused
    add \a,\a,\mx
.endif
.if \firstzero
    ror \d,\a,#16
.else
    eor \d,\d,\a
    ror \d,\d,#16
.endif
    add \c,\c,\d
    eor \b,\b,\c
    ror \b,\b,#12
    add \a,\a,\b
.if \yused
    add \a,\a,\my
.endif
    eor \d,\d,\a
    ror \d,\d,#8
    add \c,\c,\d
    eor \b,\b,\c
    ror \b,\b,#7
.endm
fapi_blake3_final:
    stp x19,x20,[sp,#-96]!
    stp x21,x22,[sp,#16]
    stp x23,x24,[sp,#32]
    stp x25,x26,[sp,#48]
    stp x27,x28,[sp,#64]
    stp x29,x30,[sp,#80]
    mov x27,x0
    mov x28,x1
    adrp x29,fapi_b3_iv4
    add x29,x29,:lo12:fapi_b3_iv4
    ldr w16,[x28,#0]
    ldr w17,[x28,#4]
    ldr w19,[x28,#8]
    ldr w20,[x28,#12]
    ldr w21,[x28,#16]
    ldr w22,[x28,#20]
    ldr w23,[x28,#24]
    ldr w24,[x28,#28]
    ldr w25,[x28,#32]
    ldr w26,[x28,#36]
    ldr w0,[x27,#0]
    ldr w1,[x27,#4]
    ldr w2,[x27,#8]
    ldr w3,[x27,#12]
    ldr w4,[x27,#16]
    ldr w5,[x27,#20]
    ldr w6,[x27,#24]
    ldr w7,[x27,#28]
    ldr w8,[x29,#0]
    ldr w9,[x29,#4]
    ldr w10,[x29,#8]
    ldr w11,[x29,#12]
// v12,v13 zero initial values are folded into their first EOR/ROR.
    mov w14,#40
    mov w15,#10
// Standard round1
    G w0,w4,w8,w12,w16,w17,1,1,1
    G w1,w5,w9,w13,w19,w20,1,1,1
    G w2,w6,w10,w14,w21,w22,1,1,0
    G w3,w7,w11,w15,w23,w24,1,1,0
    G w0,w5,w10,w15,w25,w26,1,1,0
    G w1,w6,w11,w12,wzr,wzr,0,0,0
    G w2,w7,w8,w13,wzr,wzr,0,0,0
    G w3,w4,w9,w14,wzr,wzr,0,0,0
// Standard round2
    G w0,w4,w8,w12,w19,w23,1,1,0
    G w1,w5,w9,w13,w20,wzr,1,0,0
    G w2,w6,w10,w14,w24,w16,1,1,0
    G w3,w7,w11,w15,w21,wzr,1,0,0
    G w0,w5,w10,w15,w17,wzr,1,0,0
    G w1,w6,w11,w12,wzr,w22,0,1,0
    G w2,w7,w8,w13,w26,wzr,1,0,0
    G w3,w4,w9,w14,wzr,w25,0,1,0
// Standard round3
    G w0,w4,w8,w12,w20,w21,1,1,0
    G w1,w5,w9,w13,wzr,wzr,0,0,0
    G w2,w6,w10,w14,wzr,w19,0,1,0
    G w3,w7,w11,w15,w24,wzr,1,0,0
    G w0,w5,w10,w15,w23,w22,1,1,0
    G w1,w6,w11,w12,w26,w16,1,1,0
    G w2,w7,w8,w13,wzr,wzr,0,0,0
    G w3,w4,w9,w14,w25,w17,1,1,0
// Standard round4
    G w0,w4,w8,w12,wzr,w24,0,1,0
    G w1,w5,w9,w13,wzr,w26,0,1,0
    G w2,w6,w10,w14,wzr,w20,0,1,0
    G w3,w7,w11,w15,wzr,wzr,0,0,0
    G w0,w5,w10,w15,w21,w16,1,1,0
    G w1,w6,w11,w12,wzr,w19,0,1,0
    G w2,w7,w8,w13,w22,w25,1,1,0
    G w3,w4,w9,w14,w17,w23,1,1,0
// Standard round5
    G w0,w4,w8,w12,wzr,wzr,0,0,0
    G w1,w5,w9,w13,w26,wzr,1,0,0
    G w2,w6,w10,w14,wzr,wzr,0,0,0
    G w3,w7,w11,w15,wzr,w25,0,1,0
    G w0,w5,w10,w15,w24,w19,1,1,0
    G w1,w6,w11,w12,w22,w20,1,1,0
    G w2,w7,w8,w13,w16,w17,1,1,0
    G w3,w4,w9,w14,w23,w21,1,1,0
// Standard round6
    G w0,w4,w8,w12,w26,wzr,1,0,0
    G w1,w5,w9,w13,wzr,w22,0,1,0
    G w2,w6,w10,w14,w25,wzr,1,0,0
    G w3,w7,w11,w15,wzr,w17,0,1,0
    G w0,w5,w10,w15,wzr,w20,0,1,0
    G w1,w6,w11,w12,w16,wzr,1,0,0
    G w2,w7,w8,w13,w19,w23,1,1,0
    G w3,w4,w9,w14,w21,w24,1,1,0
// Standard round7
    G w0,w4,w8,w12,wzr,wzr,0,0,0
    G w1,w5,w9,w13,w22,w16,1,1,0
    G w2,w6,w10,w14,w17,w26,1,1,0
    G w3,w7,w11,w15,w25,w23,1,1,0
    G w0,w5,w10,w15,wzr,wzr,0,0,0
    G w1,w6,w11,w12,w19,wzr,1,0,0
    G w2,w7,w8,w13,w20,w21,1,1,0
    G w3,w4,w9,w14,w24,wzr,1,0,0
    eor w0,w0,w8
    str w0,[x27,#0]
    eor w1,w1,w9
    str w1,[x27,#4]
    eor w2,w2,w10
    str w2,[x27,#8]
    eor w3,w3,w11
    str w3,[x27,#12]
    eor w4,w4,w12
    str w4,[x27,#16]
    eor w5,w5,w13
    str w5,[x27,#20]
    eor w6,w6,w14
    str w6,[x27,#24]
    eor w7,w7,w15
    str w7,[x27,#28]
    ldp x21,x22,[sp,#16]
    ldp x23,x24,[sp,#32]
    ldp x25,x26,[sp,#48]
    ldp x27,x28,[sp,#64]
    ldp x29,x30,[sp,#80]
    ldp x19,x20,[sp],#96
    ret
.size fapi_blake3_final,.-fapi_blake3_final
.section .rodata
.p2align 2
fapi_b3_iv4:
    .word 0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a

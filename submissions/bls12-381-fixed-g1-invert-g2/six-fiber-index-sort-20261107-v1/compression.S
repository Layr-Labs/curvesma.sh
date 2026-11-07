// Static analytic evidence only;never compile/execute on host or with credentials.
// x0=CV32,x1=message block64,x2=out32,w3=block length,w4=flags.
// Full standard BLAKE3-256/chaining output,all7rounds,counter0 (all used domains onechunk).
// All56G calls execute. Caller scratch w0..w17;x18/LR untouched;preservesx19..x22.
.text
.p2align 2
.macro C_G a,b,c,d,mx,my
    ldr w16,[x19,#(4*\mx)]
    ldr w17,[x19,#(4*\my)]
    add \a,\a,\b
    add \a,\a,w16
    eor \d,\d,\a
    ror \d,\d,#16
    add \c,\c,\d
    eor \b,\b,\c
    ror \b,\b,#12
    add \a,\a,\b
    add \a,\a,w17
    eor \d,\d,\a
    ror \d,\d,#8
    add \c,\c,\d
    eor \b,\b,\c
    ror \b,\b,#7
.endm
.macro C_ROUND a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15
    C_G w0,w4,w8,w12,\a0,\a1
    C_G w1,w5,w9,w13,\a2,\a3
    C_G w2,w6,w10,w14,\a4,\a5
    C_G w3,w7,w11,w15,\a6,\a7
    C_G w0,w5,w10,w15,\a8,\a9
    C_G w1,w6,w11,w12,\a10,\a11
    C_G w2,w7,w8,w13,\a12,\a13
    C_G w3,w4,w9,w14,\a14,\a15
.endm
.global coarse_b3_compress
.type coarse_b3_compress,%function
coarse_b3_compress:
    stp x19,x20,[sp,#-32]!
    stp x21,x22,[sp,#16]
    mov x19,x1
    mov x20,x2
    mov x21,x0
    mov w14,w3
    mov w15,w4
    ldp w0,w1,[x21,#0]
    ldp w2,w3,[x21,#8]
    ldp w4,w5,[x21,#16]
    ldp w6,w7,[x21,#24]
    adrp x22,coarse_b3_iv
    add x22,x22,:lo12:coarse_b3_iv
    ldp w8,w9,[x22,#0]
    ldp w10,w11,[x22,#8]
    mov w12,wzr
    mov w13,wzr
    C_ROUND 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
    C_ROUND 2,6,3,10,7,0,4,13,1,11,12,5,9,14,15,8
    C_ROUND 3,4,10,12,13,2,7,14,6,5,9,0,11,15,8,1
    C_ROUND 10,7,12,9,14,3,13,15,4,0,11,2,5,8,1,6
    C_ROUND 12,13,9,11,15,10,14,8,7,2,5,3,0,1,6,4
    C_ROUND 9,14,11,5,8,12,15,1,13,3,0,10,2,6,4,7
    C_ROUND 11,15,5,0,1,9,8,6,14,10,2,12,3,4,7,13
    eor w0,w0,w8
    eor w1,w1,w9
    eor w2,w2,w10
    eor w3,w3,w11
    eor w4,w4,w12
    eor w5,w5,w13
    eor w6,w6,w14
    eor w7,w7,w15
    stp w0,w1,[x20,#0]
    stp w2,w3,[x20,#8]
    stp w4,w5,[x20,#16]
    stp w6,w7,[x20,#24]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#32
    ret
.size coarse_b3_compress,.-coarse_b3_compress
.section .rodata
.p2align 3
.global coarse_b3_iv
coarse_b3_iv:
    .word 0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a
    .word 0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19

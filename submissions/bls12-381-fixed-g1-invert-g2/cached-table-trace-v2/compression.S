// Inert full-standard BLAKE3 single-chunk compression evidence; no host execution.
// Algorithm/count source: e750ade2 and the standard seven-round reference.
// x0=CV32,x1=block64,x2=out32,w3=length,w4=flags,counter=0.
.text
.p2align 2
.macro BG a,b,c,d,i,j
    ldr w16,[x20,#(4*\i)]
    ldr w17,[x20,#(4*\j)]
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
.macro BR m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15
    BG w0,w4,w8,w12,\m0,\m1
    BG w1,w5,w9,w13,\m2,\m3
    BG w2,w6,w10,w14,\m4,\m5
    BG w3,w7,w11,w15,\m6,\m7
    BG w0,w5,w10,w15,\m8,\m9
    BG w1,w6,w11,w12,\m10,\m11
    BG w2,w7,w8,w13,\m12,\m13
    BG w3,w4,w9,w14,\m14,\m15
.endm
.global b3_chunk_cv
.type b3_chunk_cv,%function
b3_chunk_cv:
    stp x19,x20,[sp,#-32]!
    stp x21,x22,[sp,#16]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    mov w14,w3
    mov w15,w4
    adrp x22,.LB3IV
    add x22,x22,:lo12:.LB3IV
    ldp w0,w1,[x19]
    ldp w2,w3,[x19,#8]
    ldp w4,w5,[x19,#16]
    ldp w6,w7,[x19,#24]
    ldp w8,w9,[x22]
    ldp w10,w11,[x22,#8]
    mov w12,wzr
    mov w13,wzr
    BR 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
    BR 2,6,3,10,7,0,4,13,1,11,12,5,9,14,15,8
    BR 3,4,10,12,13,2,7,14,6,5,9,0,11,15,8,1
    BR 10,7,12,9,14,3,13,15,4,0,11,2,5,8,1,6
    BR 12,13,9,11,15,10,14,8,7,2,5,3,0,1,6,4
    BR 9,14,11,5,8,12,15,1,13,3,0,10,2,6,4,7
    BR 11,15,5,0,1,9,8,6,14,10,2,12,3,4,7,13
    eor w0,w0,w8
    eor w1,w1,w9
    eor w2,w2,w10
    eor w3,w3,w11
    eor w4,w4,w12
    eor w5,w5,w13
    eor w6,w6,w14
    eor w7,w7,w15
    stp w0,w1,[x21,#0]
    stp w2,w3,[x21,#8]
    stp w4,w5,[x21,#16]
    stp w6,w7,[x21,#24]
    ldp x21,x22,[sp,#16]
    ldp x19,x20,[sp],#32
    ret
.size b3_chunk_cv,.-b3_chunk_cv
.section .rodata
.p2align 3
.LB3IV:
    .word 0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19

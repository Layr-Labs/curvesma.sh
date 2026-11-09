// Static mathematical evidence only; never execute on host or with credentials.
// Inputs x0,x1=a; x2,x3=b, canonical mod 2^127-1. Output x0,x1. Clobbers x4..x11; x12..x18 and LR preserved.
.text
.p2align 2
.global aux127_mul
.type aux127_mul,%function
aux127_mul:
    mov x4,xzr
    mov x5,xzr
    mov x6,xzr
    mov x7,xzr
    mul x8,x0,x2
    umulh x9,x0,x2
    adds x4,x4,x8
    adcs x5,x5,x9
    adcs x6,x6,xzr
    adc x7,x7,xzr
    mul x8,x0,x3
    umulh x9,x0,x3
    adds x5,x5,x8
    adcs x6,x6,x9
    adc x7,x7,xzr
    mul x8,x1,x2
    umulh x9,x1,x2
    adds x5,x5,x8
    adcs x6,x6,x9
    adc x7,x7,xzr
    mul x8,x1,x3
    umulh x9,x1,x3
    adds x6,x6,x8
    adc x7,x7,x9
    extr x8,x6,x5,#63
    extr x9,x7,x6,#63
    and x5,x5,#0x7fffffffffffffff
    adds x4,x4,x8
    adc x5,x5,x9
    lsr x8,x5,#63
    and x5,x5,#0x7fffffffffffffff
    adds x4,x4,x8
    adc x5,x5,xzr
    movn x10,#0
    orr x11,xzr,#0x7fffffffffffffff
    subs x8,x4,x10
    sbcs x9,x5,x11
    csel x0,x8,x4,hs
    csel x1,x9,x5,hs
    ret
.size aux127_mul,.-aux127_mul

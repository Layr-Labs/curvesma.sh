// Static mathematical evidence only; never execute on host or with credentials.
// Canonical a=x0..x3,b=x4..x7 modulo 2^255-19, output x0..x3; clobbers x4..x14, preserves x15..x30.
.text
.p2align 2
.global aux255_add
.type aux255_add,%function
aux255_add:
    adds x0,x0,x4
    adcs x1,x1,x5
    adcs x2,x2,x6
    adc x3,x3,x7
    mov x8,#19
    lsr x9,x3,#63
    and x3,x3,#0x7fffffffffffffff
    mul x9,x9,x8
    adds x0,x0,x9
    adcs x1,x1,xzr
    adcs x2,x2,xzr
    adc x3,x3,xzr
    movn x8,#18
    movn x9,#0
    orr x10,xzr,#0x7fffffffffffffff
    subs x11,x0,x8
    sbcs x12,x1,x9
    sbcs x13,x2,x9
    sbcs x14,x3,x10
    csel x0,x11,x0,hs
    csel x1,x12,x1,hs
    csel x2,x13,x2,hs
    csel x3,x14,x3,hs
    ret
.size aux255_add,.-aux255_add

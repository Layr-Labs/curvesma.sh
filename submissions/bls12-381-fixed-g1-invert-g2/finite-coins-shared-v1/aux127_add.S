// Static mathematical evidence only; never execute on host or with credentials.
// Canonical mod 2^127-1, a=x0,x1 b=x2,x3, output x0,x1; clobbers x4..x8, preserves x9..x18 and LR.
.text
.p2align 2
.global aux127_add
.type aux127_add,%function
aux127_add:
    adds x0,x0,x2
    adc x1,x1,x3
    lsr x4,x1,#63
    and x1,x1,#0x7fffffffffffffff
    adds x0,x0,x4
    adc x1,x1,xzr
    movn x5,#0
    orr x6,xzr,#0x7fffffffffffffff
    subs x7,x0,x5
    sbcs x8,x1,x6
    csel x0,x7,x0,hs
    csel x1,x8,x1,hs
    ret
.size aux127_add,.-aux127_add

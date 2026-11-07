// Static analytic evidence;never execute on host or with credentials.
// x0=state24(x,m,b),x1=out32(jlo,jhi,klo,khi),native little-endian limbs.
// State starts(0,0,0). Only first L<g^3 calls occur,g=2^42-11.
.text
.p2align 2
.global fapi_incidence_next
.type fapi_incidence_next,%function
fapi_incidence_next:
    ldp x5,x6,[x0,#0]
    ldr x7,[x0,#16]
    mov x3,#-11
    and x3,x3,#0x3ffffffffff
    mov x4,#11
    mul x8,x5,x6
    umulh x9,x5,x6
    adds x8,x8,x7
    adc x9,x9,xzr
    lsr x10,x8,#42
    orr x10,x10,x9,lsl #22
    and x8,x8,#0x3ffffffffff
    madd x10,x10,x4,x8
    lsr x9,x10,#42
    and x8,x10,#0x3ffffffffff
    madd x10,x9,x4,x8
    subs x11,x10,x3
    csel x13,x11,x10,hs
    mul x8,x5,x3
    umulh x9,x5,x3
    adds x8,x8,x13
    adc x9,x9,xzr
    mul x10,x6,x3
    umulh x11,x6,x3
    adds x10,x10,x7
    adc x11,x11,xzr
    stp x8,x9,[x1,#0]
    stp x10,x11,[x1,#16]
    add x7,x7,#1
    cmp x7,x3
    csel x7,xzr,x7,eq
    cset x12,eq
    add x6,x6,x12
    cmp x6,x3
    cset x12,eq
    csel x6,xzr,x6,eq
    add x5,x5,x12
    stp x5,x6,[x0,#0]
    str x7,[x0,#16]
    ret
.size fapi_incidence_next,.-fapi_incidence_next

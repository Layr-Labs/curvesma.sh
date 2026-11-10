// Inert analytic cold-fill routine;NEVER host-execute or use credentials.
// x0=base for n records48,x1=canonical seed48,x2=p48,x3=n in{2^39,2^38}.
.text
.p2align 2
.global fill_lut
.type fill_lut,%function
fill_lut:
    stp x29,x30,[sp,#-48]!
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    sub x22,x3,#1
    stp xzr,xzr,[x19]
    stp xzr,xzr,[x19,#16]
    stp xzr,xzr,[x19,#32]
.L_lut_next:
    add x0,x19,#48
    mov x1,x19
    mov x2,x20
    mov x3,x21
    bl fp381_add
    add x19,x19,#48
    subs x22,x22,#1
    b.ne .L_lut_next
    ldp x21,x22,[sp,#32]
    ldp x19,x20,[sp,#16]
    ldp x29,x30,[sp],#48
    ret
.size fill_lut,.-fill_lut

// Inert analytic source;NEVER host-execute/compile or use credentials.
// x0 third_out192,x1 prev1_192,x2 prev2_192,x3 tablebase,x4 const520,x5 p48;histories distinct.
// Frame480 holds initialized transformed384@96. Exact in-place third update;key output must NOTalias history.
// Callees19..29/SP/LRpreserved,no18. No raw-history copy needed after return;callerrotates4MOV.
.text
.p2align 2
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    sub sp,sp,#480
    stp x29,x30,[sp]
    mov x29,sp
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    stp x23,x24,[sp,#48]
    stp x25,x26,[sp,#64]
    stp x27,x28,[sp,#80]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    mov x22,x3
    mov x23,x4
    mov x24,x5
    add x25,x29,#96
    add x0,x25,#0
    add x1,x20,#0
    add x2,x21,#0
    mov x3,x24
    bl fp381_sub
    add x0,x25,#48
    add x1,x20,#48
    add x2,x21,#48
    mov x3,x24
    bl fp381_sub
    add x0,x25,#96
    add x1,x20,#96
    add x2,x21,#96
    mov x3,x24
    bl fp381_add
    add x0,x25,#144
    add x1,x20,#144
    add x2,x21,#144
    mov x3,x24
    bl fp381_add
    add x0,x25,#192
    add x1,x20,#96
    add x2,x21,#96
    mov x3,x24
    bl fp381_sub
    add x0,x25,#240
    add x1,x20,#144
    add x2,x21,#144
    mov x3,x24
    bl fp381_sub
    add x0,x25,#288
    add x1,x20,#0
    add x2,x21,#0
    mov x3,x24
    bl fp381_add
    add x0,x25,#336
    add x1,x20,#48
    add x2,x21,#48
    mov x3,x24
    bl fp381_add
    mov x0,x19
    mov x1,x25
    mov x2,x22
    mov x3,x23
    bl recurrence_fp2_map
    ldp x19,x20,[sp,#16]
    ldp x21,x22,[sp,#32]
    ldp x23,x24,[sp,#48]
    ldp x25,x26,[sp,#64]
    ldp x27,x28,[sp,#80]
    ldp x29,x30,[sp]
    add sp,sp,#480
    ret
.size recurrence_fp2_step,.-recurrence_fp2_step

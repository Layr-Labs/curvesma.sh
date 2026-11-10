// Inert cold-fill evidence;NEVER host execute/compile or use credentials.
// x0 base n*192,x1 seed192,x2 p48,x3 n in{2^29,2^32};all regions disjoint.
.text
.p2align 2
.global fill_vector
.type fill_vector,%function
fill_vector:
    stp x29,x30,[sp,#-48]!
    stp x19,x20,[sp,#16]
    stp x21,x22,[sp,#32]
    mov x19,x0
    mov x20,x1
    mov x21,x2
    sub x22,x3,#1
    stp xzr,xzr,[x19,#0]
    stp xzr,xzr,[x19,#16]
    stp xzr,xzr,[x19,#32]
    stp xzr,xzr,[x19,#48]
    stp xzr,xzr,[x19,#64]
    stp xzr,xzr,[x19,#80]
    stp xzr,xzr,[x19,#96]
    stp xzr,xzr,[x19,#112]
    stp xzr,xzr,[x19,#128]
    stp xzr,xzr,[x19,#144]
    stp xzr,xzr,[x19,#160]
    stp xzr,xzr,[x19,#176]
.L_vector_next:
    add x0,x19,#192
    add x1,x19,#0
    add x2,x20,#0
    mov x3,x21
    bl fp381_add
    add x0,x19,#240
    add x1,x19,#48
    add x2,x20,#48
    mov x3,x21
    bl fp381_add
    add x0,x19,#288
    add x1,x19,#96
    add x2,x20,#96
    mov x3,x21
    bl fp381_add
    add x0,x19,#336
    add x1,x19,#144
    add x2,x20,#144
    mov x3,x21
    bl fp381_add
    add x19,x19,#192
    subs x22,x22,#1
    b.ne .L_vector_next
    ldp x21,x22,[sp,#32]
    ldp x19,x20,[sp,#16]
    ldp x29,x30,[sp],#48
    ret
.size fill_vector,.-fill_vector

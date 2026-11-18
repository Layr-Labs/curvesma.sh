.macro INV_S1
 mov x0,x21
 mov x1,x21
 bl k_sqr
.endm
.macro INV_S2
 INV_S1
 INV_S1
.endm
.macro INV_S4
 INV_S2
 INV_S2
.endm
.macro INV_S8
 INV_S4
 INV_S4
.endm
.macro INV_S16
 INV_S8
 INV_S8
.endm
.macro INV_S32
 INV_S16
 INV_S16
.endm
.macro INV_S64
 INV_S32
 INV_S32
.endm
.macro INV_S128
 INV_S64
 INV_S64
.endm
.globl k_inv
k_inv:
 sub sp,sp,#144
 str x30,[sp]
 mov x19,x0
 mov x20,x1
 mov x21,x2
 ldp x9,x10,[x20]
 ldp x11,x12,[x20,#16]
 stp x9,x10,[sp,#16]
 stp x11,x12,[sp,#32]
 stp x9,x10,[x21]
 stp x11,x12,[x21,#16]
 INV_S1
 mov x0,x21
 mov x1,x21
 add x2,sp,#16
 bl k_mul
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 stp x9,x10,[sp,#48]
 stp x11,x12,[sp,#64]
 INV_S2
 mov x0,x21
 mov x1,x21
 add x2,sp,#48
 bl k_mul
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 stp x9,x10,[sp,#80]
 stp x11,x12,[sp,#96]
 INV_S4
 mov x0,x21
 mov x1,x21
 add x2,sp,#80
 bl k_mul
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 stp x9,x10,[sp,#112]
 stp x11,x12,[sp,#128]
 INV_S8
 mov x0,x21
 mov x1,x21
 add x2,sp,#112
 bl k_mul
 INV_S8
 mov x0,x21
 mov x1,x21
 add x2,sp,#112
 bl k_mul
 INV_S4
 mov x0,x21
 mov x1,x21
 add x2,sp,#80
 bl k_mul
 INV_S2
 mov x0,x21
 mov x1,x21
 add x2,sp,#48
 bl k_mul
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 stp x9,x10,[sp,#112]
 stp x11,x12,[sp,#128]
 INV_S2
 mov x0,x21
 mov x1,x21
 add x2,sp,#48
 bl k_mul
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 stp x9,x10,[sp,#80]
 stp x11,x12,[sp,#96]
 INV_S32
 mov x0,x21
 mov x1,x21
 add x2,sp,#16
 bl k_mul
 INV_S128
 mov x0,x21
 mov x1,x21
 add x2,sp,#80
 bl k_mul
 INV_S32
 mov x0,x21
 mov x1,x21
 add x2,sp,#80
 bl k_mul
 INV_S16
 INV_S8
 INV_S4
 INV_S2
 mov x0,x21
 mov x1,x21
 add x2,sp,#112
 bl k_mul
 INV_S2
 mov x0,x21
 mov x1,x21
 add x2,sp,#16
 bl k_mul
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 stp x9,x10,[x19]
 stp x11,x12,[x19,#16]
 ldr x30,[sp]
 add sp,sp,#144
 ret

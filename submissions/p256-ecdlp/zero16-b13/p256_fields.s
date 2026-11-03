.globl k_mul
k_mul:
 ldp x3, x4, [x1]
 ldp x7, x8, [x2]
 mul x9, x3, x7
 umulh x10, x3, x7
 mul x15, x3, x8
 umulh x11, x3, x8
 adds x10, x10, x15
 ldp x13, x14, [x2, #16]
 mul x15, x3, x13
 umulh x12, x3, x13
 adcs x11, x11, x15
 mul x15, x3, x14
 umulh x16, x3, x14
 adcs x12, x12, x15
 adc x16, x16, xzr
 ldp x5, x6, [x1, #16]
 mul x15, x4, x7
 adds x10, x10, x15
 mul x15, x4, x8
 adcs x11, x11, x15
 mul x15, x4, x13
 adcs x12, x12, x15
 mul x15, x4, x14
 adcs x16, x16, x15
 umulh x3, x4, x14
 adc x3, x3, xzr
 umulh x15, x4, x7
 adds x11, x11, x15
 umulh x15, x4, x8
 adcs x12, x12, x15
 umulh x15, x4, x13
 adcs x16, x16, x15
 adc x3, x3, xzr
 mul x15, x5, x7
 adds x11, x11, x15
 mul x15, x5, x8
 adcs x12, x12, x15
 mul x15, x5, x13
 adcs x16, x16, x15
 mul x15, x5, x14
 adcs x3, x3, x15
 umulh x4, x5, x14
 adc x4, x4, xzr
 umulh x15, x5, x7
 adds x12, x12, x15
 umulh x15, x5, x8
 adcs x16, x16, x15
 umulh x15, x5, x13
 adcs x3, x3, x15
 adc x4, x4, xzr
 mul x15, x6, x7
 adds x12, x12, x15
 mul x15, x6, x8
 adcs x16, x16, x15
 mul x15, x6, x13
 adcs x3, x3, x15
 mul x15, x6, x14
 adcs x4, x4, x15
 umulh x5, x6, x14
 adc x5, x5, xzr
 mov x17, #0xFFFFFFFF00000001
 adds x10,x10,x9,lsl #32
 adcs x11,x11,x9,lsr #32
 mul x15,x9,x17
 umulh x9,x9,x17
 adcs x12,x12,x15
 adc x9,x9,xzr
 umulh x15, x6, x7
 adds x16, x16, x15
 umulh x15, x6, x8
 adcs x3, x3, x15
 umulh x15, x6, x13
 adcs x4, x4, x15
 adc x5, x5, xzr
 adds x11,x11,x10,lsl #32
 adcs x12,x12,x10,lsr #32
 mul x15,x10,x17
 umulh x10,x10,x17
 adcs x9,x9,x15
 adc x10,x10,xzr
 adds x12,x12,x11,lsl #32
 adcs x9,x9,x11,lsr #32
 mul x15,x11,x17
 umulh x11,x11,x17
 adcs x10,x10,x15
 adc x11,x11,xzr
 adds x9,x9,x12,lsl #32
 adcs x10,x10,x12,lsr #32
 mul x15,x12,x17
 umulh x12,x12,x17
 adcs x11,x11,x15
 adc x12,x12,xzr
 adds x9, x9, x16
 adcs x10, x10, x3
 adcs x11, x11, x4
 adcs x12, x12, x5
 b.cc k_mul_store
 mov w3, #-1
 mov w4, #-2
 adds x9,x9,#1
 sbcs x10,x10,x3
 sbcs x11,x11,xzr
 adc x12,x12,x4
k_mul_store:
 stp x9,x10,[x0]
 stp x11,x12,[x0,#16]
 ret
.globl k_sqr
k_sqr:
 ldp x3, x4, [x1]
 ldp x5, x6, [x1, #16]
 mul x8, x4, x3
 umulh x16, x4, x3
 mul x9, x5, x3
 umulh x17, x5, x3
 mul x10, x6, x3
 umulh x11, x6, x3
 adds x9, x9, x16
 mul x15, x5, x4
 umulh x16, x5, x4
 adcs x10, x10, x17
 mul x17, x6, x4
 umulh x1, x6, x4
 adc x11, x11, xzr
 mul x12, x6, x5
 umulh x13, x6, x5
 adds x16, x16, x17
 mul x7, x3, x3
 adc x17, x1, xzr
 adds x10, x10, x15
 umulh x3, x3, x3
 adcs x11, x11, x16
 mul x16, x4, x4
 adcs x12, x12, x17
 umulh x4, x4, x4
 adc x13, x13, xzr
 adds x8, x8, x8
 mul x17, x5, x5
 adcs x9, x9, x9
 umulh x5, x5, x5
 adcs x10, x10, x10
 mul x1, x6, x6
 adcs x11, x11, x11
 umulh x6, x6, x6
 adcs x12, x12, x12
 adcs x13, x13, x13
 adc x14, xzr, xzr
 adds x8, x8, x3
 adcs x9, x9, x16
 adcs x10, x10, x4
 adcs x11, x11, x17
 adcs x12, x12, x5
 adcs x13, x13, x1
 lsr x16, x7, #32
 adc x14, x14, x6
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x9, x16
 adcs x9, x10, x17
 lsr x16, x7, #32
 adc x10, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x9, x16
 adcs x9, x10, x17
 lsr x16, x7, #32
 adc x10, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x9, x16
 adcs x9, x10, x17
 lsr x16, x7, #32
 adc x10, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x9, x16
 adcs x9, x10, x17
 adc x10, x1, xzr
 adds x7, x7, x11
 adcs x8, x8, x12
 adcs x9, x9, x13
 adcs x10, x10, x14
 adc x11, xzr, xzr
 cbz x11,k_sqr_store
 mov w3,#-1
 mov w4,#-2
 adds x7,x7,#1
 sbcs x8,x8,x3
 sbcs x9,x9,xzr
 adc x10,x10,x4
k_sqr_store:
 stp x7,x8,[x0]
 stp x9,x10,[x0,#16]
 ret
.globl k_add
k_add:
 ldp x3, x4, [x1]
 ldp x5, x6, [x1, #16]
 ldp x7, x8, [x2]
 ldp x9, x10, [x2, #16]
 adds x11, x3, x7
 adcs x12, x4, x8
 adcs x13, x5, x9
 adcs x14, x6, x10
 adc x15, xzr, xzr
 cbz x15,k_add_store
 mov w3,#-1
 mov w4,#-2
 adds x11,x11,#1
 sbcs x12,x12,x3
 sbcs x13,x13,xzr
 adc x14,x14,x4
k_add_store:
 stp x11,x12,[x0]
 stp x13,x14,[x0,#16]
 ret
.globl k_sub
k_sub:
 ldp x3, x4, [x1]
 ldp x5, x6, [x1, #16]
 ldp x7, x8, [x2]
 ldp x9, x10, [x2, #16]
 subs x11, x3, x7
 sbcs x12, x4, x8
 sbcs x13, x5, x9
 sbcs x14, x6, x10
 sbc x15, xzr, xzr
 lsr x16, x15, #32
 lsl x17, x15, #32
 sub x17, x17, x15
 adds x11, x11, x15
 adcs x12, x12, x16
 adcs x13, x13, xzr
 adc x14, x14, x17
 stp x11, x12, [x0]
 stp x13, x14, [x0, #16]
 ret

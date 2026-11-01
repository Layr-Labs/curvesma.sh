.globl k_fused
k_fused:
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
 ldp x12,x13,[x2]
 ldp x14,x15,[x2,#16]
 subs x7,x7,x12
 sbcs x8,x8,x13
 sbcs x9,x9,x14
 sbcs x10,x10,x15
 sbc x11,x11,xzr
 cbz x11,fs_store
 tbnz x11,#63,fs_negative
 mov w3,#-1
 mov w4,#-2
 adds x7,x7,#1
 sbcs x8,x8,x3
 sbcs x9,x9,xzr
 adc x10,x10,x4
 b fs_store
fs_negative:
 mov w3,#-1
 mov w4,#-2
 subs x7,x7,#1
 adcs x8,x8,x3
 adcs x9,x9,xzr
 sbc x10,x10,x4
fs_store:
 stp x7,x8,[x0]
 stp x9,x10,[x0,#16]
 lsr x15,x10,#32
 mov w3,#-1
 cmp x15,x3
 cset x0,eq
 ret

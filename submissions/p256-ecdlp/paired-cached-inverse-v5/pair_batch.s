// Custom internal ABI: x0 count, x1 first/last Y, x2 first/last batch,
// x3 X point, x4 scratch. All GPRs may be clobbered except SP and LR.
// Batch has 32-byte Montgomery-one prefix before its first 64-byte entry.
.globl pair_fwd
pair_fwd:
 str x30,[sp,#-16]!
 mov x19,x0
 mov x20,x1
 mov x21,x2
 mov x22,x3
 sub x23,x2,#32
 cbz x19,pf_success
pf_loop:
 // inline k_sub; stable caller pointers, unchanged arithmetic
 ldp x3, x4, [x20]
 ldp x5, x6, [x20, #16]
 ldp x7, x8, [x22]
 ldp x9, x10, [x22, #16]
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
 stp x11, x12, [x21]
 stp x13, x14, [x21, #16]
 orr x9,x11,x12
 orr x10,x13,x14
 orr x9,x9,x10
 cbz x9,pf_fail
 // inline k_mul; stable caller pointers, unchanged arithmetic
 ldp x3, x4, [x23]
 ldp x5, x6, [x23, #16]
 mul x9, x3, x11
 umulh x15, x3, x11
 mul x10, x4, x11
 umulh x16, x4, x11
 mul x24, x5, x11
 umulh x17, x5, x11
 mul x25, x6, x11
 umulh x1, x6, x11
 adds x10, x10, x15
 adcs x24, x24, x16
 lsr x16, x9, #32
 adcs x25, x25, x17
 adc x26, xzr, x1
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x3, x12
 adcs x10, x24, x16
 mul x16, x4, x12
 adcs x24, x25, x17
 mul x17, x5, x12
 adcs x25, x26, x1
 mul x1, x6, x12
 adc x26, xzr, xzr
 adds x9, x9, x15
 umulh x15, x3, x12
 adcs x10, x10, x16
 umulh x16, x4, x12
 adcs x24, x24, x17
 umulh x17, x5, x12
 adcs x25, x25, x1
 umulh x1, x6, x12
 adc x26, x26, xzr
 adds x10, x10, x15
 adcs x24, x24, x16
 lsr x16, x9, #32
 adcs x25, x25, x17
 adcs x26, x26, x1
 adc x27, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x3, x13
 adcs x10, x24, x16
 mul x16, x4, x13
 adcs x24, x25, x17
 mul x17, x5, x13
 adcs x25, x26, x1
 mul x1, x6, x13
 adc x26, x27, xzr
 adds x9, x9, x15
 umulh x15, x3, x13
 adcs x10, x10, x16
 umulh x16, x4, x13
 adcs x24, x24, x17
 umulh x17, x5, x13
 adcs x25, x25, x1
 umulh x1, x6, x13
 adc x26, x26, xzr
 adds x10, x10, x15
 adcs x24, x24, x16
 lsr x16, x9, #32
 adcs x25, x25, x17
 adcs x26, x26, x1
 adc x27, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x3, x14
 adcs x10, x24, x16
 mul x16, x4, x14
 adcs x24, x25, x17
 mul x17, x5, x14
 adcs x25, x26, x1
 mul x1, x6, x14
 adc x26, x27, xzr
 adds x9, x9, x15
 umulh x15, x3, x14
 adcs x10, x10, x16
 umulh x16, x4, x14
 adcs x24, x24, x17
 umulh x17, x5, x14
 adcs x25, x25, x1
 umulh x1, x6, x14
 adc x26, x26, xzr
 adds x10, x10, x15
 adcs x24, x24, x16
 lsr x16, x9, #32
 adcs x25, x25, x17
 adcs x26, x26, x1
 adc x27, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 adcs x10, x24, x16
 adcs x24, x25, x17
 adcs x25, x26, x1
 adc x26, x27, xzr
 mov w3, #-1
 mov w4, #-2
 adds x15, x9, #1
 sbcs x16, x10, x3
 sbcs x17, x24, xzr
 adcs x1, x25, x4
 sbcs xzr, x26, xzr
 csel x9, x9, x15, lo
 csel x10, x10, x16, lo
 csel x24, x24, x17, lo
 stp x9, x10, [x21, #32]
 csel x25, x25, x1, lo
 stp x24, x25, [x21, #48]
 add x23,x21,#32
 add x21,x21,#64
 add x20,x20,#64
 subs x19,x19,#1
 b.ne pf_loop
pf_success:
 mov x0,#1
 b pf_done
pf_fail:
 mov x0,#0
pf_done:
 ldr x30,[sp],#16
 ret
.globl pair_bwd
pair_bwd:
 str x30,[sp,#-16]!
 mov x19,x0
 mov x20,x1
 mov x21,x2
 mov x22,x3
 mov x23,x4
 cbz x19,pb_done
pb_loop:
 ldp x3, x4, [x23]
 ldp x5, x6, [x23, #16]
 ldp x7, x8, [x21, #-32]
 mul x24, x3, x7
 umulh x15, x3, x7
 mul x25, x4, x7
 umulh x16, x4, x7
 mul x26, x5, x7
 umulh x17, x5, x7
 mul x27, x6, x7
 umulh x1, x6, x7
 adds x25, x25, x15
 adcs x26, x26, x16
 lsr x16, x24, #32
 adcs x27, x27, x17
 adc x13, xzr, x1
 subs x17, x24, x24, lsl #32
 sbc x1, x24, x16
 adds x24, x25, x24, lsl #32
 mul x15, x3, x8
 adcs x25, x26, x16
 mul x16, x4, x8
 adcs x26, x27, x17
 mul x17, x5, x8
 adcs x27, x13, x1
 mul x1, x6, x8
 adc x13, xzr, xzr
 adds x24, x24, x15
 umulh x15, x3, x8
 adcs x25, x25, x16
 umulh x16, x4, x8
 adcs x26, x26, x17
 umulh x17, x5, x8
 adcs x27, x27, x1
 umulh x1, x6, x8
 adc x13, x13, xzr
 ldp x7, x8, [x21, #-16]
 adds x25, x25, x15
 adcs x26, x26, x16
 lsr x16, x24, #32
 adcs x27, x27, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x24, x24, lsl #32
 sbc x1, x24, x16
 adds x24, x25, x24, lsl #32
 mul x15, x3, x7
 adcs x25, x26, x16
 mul x16, x4, x7
 adcs x26, x27, x17
 mul x17, x5, x7
 adcs x27, x13, x1
 mul x1, x6, x7
 adc x13, x14, xzr
 adds x24, x24, x15
 umulh x15, x3, x7
 adcs x25, x25, x16
 umulh x16, x4, x7
 adcs x26, x26, x17
 umulh x17, x5, x7
 adcs x27, x27, x1
 umulh x1, x6, x7
 adc x13, x13, xzr
 adds x25, x25, x15
 adcs x26, x26, x16
 lsr x16, x24, #32
 adcs x27, x27, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x24, x24, lsl #32
 sbc x1, x24, x16
 adds x24, x25, x24, lsl #32
 mul x15, x3, x8
 adcs x25, x26, x16
 mul x16, x4, x8
 adcs x26, x27, x17
 mul x17, x5, x8
 adcs x27, x13, x1
 mul x1, x6, x8
 adc x13, x14, xzr
 adds x24, x24, x15
 umulh x15, x3, x8
 adcs x25, x25, x16
 umulh x16, x4, x8
 adcs x26, x26, x17
 umulh x17, x5, x8
 adcs x27, x27, x1
 umulh x1, x6, x8
 adc x13, x13, xzr
 adds x25, x25, x15
 adcs x26, x26, x16
 lsr x16, x24, #32
 adcs x27, x27, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x24, x24, lsl #32
 sbc x1, x24, x16
 adds x24, x25, x24, lsl #32
 adcs x25, x26, x16
 adcs x26, x27, x17
 adcs x27, x13, x1
 adc x13, x14, xzr
 mov w3, #-1
 mov w4, #-2
 adds x15, x24, #1
 sbcs x16, x25, x3
 sbcs x17, x26, xzr
 adcs x1, x27, x4
 sbcs xzr, x13, xzr
 csel x24, x24, x15, lo
 csel x25, x25, x16, lo
 csel x26, x26, x17, lo
 csel x27, x27, x1, lo
 ldp x3, x4, [x23]
 ldp x5, x6, [x23, #16]
 ldp x7, x8, [x21]
 mul x9, x3, x7
 umulh x15, x3, x7
 mul x10, x4, x7
 umulh x16, x4, x7
 mul x11, x5, x7
 umulh x17, x5, x7
 mul x12, x6, x7
 umulh x1, x6, x7
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adc x13, xzr, x1
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x3, x8
 adcs x10, x11, x16
 mul x16, x4, x8
 adcs x11, x12, x17
 mul x17, x5, x8
 adcs x12, x13, x1
 mul x1, x6, x8
 adc x13, xzr, xzr
 adds x9, x9, x15
 umulh x15, x3, x8
 adcs x10, x10, x16
 umulh x16, x4, x8
 adcs x11, x11, x17
 umulh x17, x5, x8
 adcs x12, x12, x1
 umulh x1, x6, x8
 adc x13, x13, xzr
 ldp x7, x8, [x21, #16]
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x3, x7
 adcs x10, x11, x16
 mul x16, x4, x7
 adcs x11, x12, x17
 mul x17, x5, x7
 adcs x12, x13, x1
 mul x1, x6, x7
 adc x13, x14, xzr
 adds x9, x9, x15
 umulh x15, x3, x7
 adcs x10, x10, x16
 umulh x16, x4, x7
 adcs x11, x11, x17
 umulh x17, x5, x7
 adcs x12, x12, x1
 umulh x1, x6, x7
 adc x13, x13, xzr
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x3, x8
 adcs x10, x11, x16
 mul x16, x4, x8
 adcs x11, x12, x17
 mul x17, x5, x8
 adcs x12, x13, x1
 mul x1, x6, x8
 adc x13, x14, xzr
 adds x9, x9, x15
 umulh x15, x3, x8
 adcs x10, x10, x16
 umulh x16, x4, x8
 adcs x11, x11, x17
 umulh x17, x5, x8
 adcs x12, x12, x1
 umulh x1, x6, x8
 adc x13, x13, xzr
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 adcs x10, x11, x16
 adcs x11, x12, x17
 adcs x12, x13, x1
 adc x13, x14, xzr
 mov w3, #-1
 mov w4, #-2
 adds x15, x9, #1
 sbcs x16, x10, x3
 sbcs x17, x11, xzr
 adcs x1, x12, x4
 sbcs xzr, x13, xzr
 csel x9, x9, x15, lo
 csel x10, x10, x16, lo
 csel x11, x11, x17, lo
 stp x9, x10, [x23]
 csel x12, x12, x1, lo
 stp x11, x12, [x23, #16]
 ldp x3, x4, [x20]
 ldp x5, x6, [x20, #16]
 ldp x7, x8, [x22]
 ldp x9, x10, [x22, #16]
 adds x11, x3, x7
 adcs x12, x4, x8
 adcs x13, x5, x9
 adcs x14, x6, x10
 adc x15, xzr, xzr
 mov w3, #-1
 mov w4, #-2
 adds x7, x11, #1
 sbcs x8, x12, x3
 sbcs x9, x13, xzr
 adcs x10, x14, x4
 sbcs xzr, x15, xzr
 csel x11, x11, x7, lo
 csel x12, x12, x8, lo
 stp x11, x12, [x23, #192]
 csel x13, x13, x9, lo
 csel x14, x14, x10, lo
 stp x13, x14, [x23, #208]
 ldp x3, x4, [x20, #32]
 ldp x5, x6, [x20, #48]
 ldp x7, x8, [x22, #32]
 ldp x9, x10, [x22, #48]
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
 adds x28, x3, x7
 adcs x29, x4, x8
 adcs x30, x5, x9
 adcs x18, x6, x10
 adc x0, xzr, xzr
 mov w3, #-1
 mov w4, #-2
 adds x7, x28, #1
 sbcs x8, x29, x3
 sbcs x9, x30, xzr
 adcs x10, x18, x4
 sbcs xzr, x0, xzr
 csel x28, x28, x7, lo
 csel x29, x29, x8, lo
 csel x30, x30, x9, lo
 csel x18, x18, x10, lo
 mul x9, x11, x24
 umulh x15, x11, x24
 mul x10, x12, x24
 umulh x16, x12, x24
 mul x3, x13, x24
 umulh x17, x13, x24
 mul x4, x14, x24
 umulh x1, x14, x24
 adds x10, x10, x15
 adcs x3, x3, x16
 lsr x16, x9, #32
 adcs x4, x4, x17
 adc x5, xzr, x1
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x11, x25
 adcs x10, x3, x16
 mul x16, x12, x25
 adcs x3, x4, x17
 mul x17, x13, x25
 adcs x4, x5, x1
 mul x1, x14, x25
 adc x5, xzr, xzr
 adds x9, x9, x15
 umulh x15, x11, x25
 adcs x10, x10, x16
 umulh x16, x12, x25
 adcs x3, x3, x17
 umulh x17, x13, x25
 adcs x4, x4, x1
 umulh x1, x14, x25
 adc x5, x5, xzr
 adds x10, x10, x15
 adcs x3, x3, x16
 lsr x16, x9, #32
 adcs x4, x4, x17
 adcs x5, x5, x1
 adc x6, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x11, x26
 adcs x10, x3, x16
 mul x16, x12, x26
 adcs x3, x4, x17
 mul x17, x13, x26
 adcs x4, x5, x1
 mul x1, x14, x26
 adc x5, x6, xzr
 adds x9, x9, x15
 umulh x15, x11, x26
 adcs x10, x10, x16
 umulh x16, x12, x26
 adcs x3, x3, x17
 umulh x17, x13, x26
 adcs x4, x4, x1
 umulh x1, x14, x26
 adc x5, x5, xzr
 adds x10, x10, x15
 adcs x3, x3, x16
 lsr x16, x9, #32
 adcs x4, x4, x17
 adcs x5, x5, x1
 adc x6, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x11, x27
 adcs x10, x3, x16
 mul x16, x12, x27
 adcs x3, x4, x17
 mul x17, x13, x27
 adcs x4, x5, x1
 mul x1, x14, x27
 adc x5, x6, xzr
 adds x9, x9, x15
 umulh x15, x11, x27
 adcs x10, x10, x16
 umulh x16, x12, x27
 adcs x3, x3, x17
 umulh x17, x13, x27
 adcs x4, x4, x1
 umulh x1, x14, x27
 adc x5, x5, xzr
 adds x10, x10, x15
 adcs x3, x3, x16
 lsr x16, x9, #32
 adcs x4, x4, x17
 adcs x5, x5, x1
 adc x6, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 adcs x10, x3, x16
 adcs x3, x4, x17
 adcs x4, x5, x1
 adc x5, x6, xzr
 mov w11, #-1
 mov w12, #-2
 adds x15, x9, #1
 sbcs x16, x10, x11
 sbcs x17, x3, xzr
 adcs x1, x4, x12
 sbcs xzr, x5, xzr
 csel x9, x9, x15, lo
 csel x10, x10, x16, lo
 csel x3, x3, x17, lo
 csel x4, x4, x1, lo
 mul x8, x10, x9
 umulh x16, x10, x9
 mul x5, x3, x9
 umulh x17, x3, x9
 mul x6, x4, x9
 umulh x11, x4, x9
 adds x5, x5, x16
 mul x15, x3, x10
 umulh x16, x3, x10
 adcs x6, x6, x17
 mul x17, x4, x10
 umulh x1, x4, x10
 adc x11, x11, xzr
 mul x12, x4, x3
 umulh x13, x4, x3
 adds x16, x16, x17
 mul x7, x9, x9
 adc x17, x1, xzr
 adds x6, x6, x15
 umulh x9, x9, x9
 adcs x11, x11, x16
 mul x16, x10, x10
 adcs x12, x12, x17
 umulh x10, x10, x10
 adc x13, x13, xzr
 adds x8, x8, x8
 mul x17, x3, x3
 adcs x5, x5, x5
 umulh x3, x3, x3
 adcs x6, x6, x6
 mul x1, x4, x4
 adcs x11, x11, x11
 umulh x4, x4, x4
 adcs x12, x12, x12
 adcs x13, x13, x13
 adc x14, xzr, xzr
 adds x8, x8, x9
 adcs x5, x5, x16
 adcs x6, x6, x10
 adcs x11, x11, x17
 adcs x12, x12, x3
 adcs x13, x13, x1
 lsr x16, x7, #32
 adc x14, x14, x4
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x5, x16
 adcs x5, x6, x17
 lsr x16, x7, #32
 adc x6, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x5, x16
 adcs x5, x6, x17
 lsr x16, x7, #32
 adc x6, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x5, x16
 adcs x5, x6, x17
 lsr x16, x7, #32
 adc x6, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x5, x16
 adcs x5, x6, x17
 adc x6, x1, xzr
 adds x7, x7, x11
 adcs x8, x8, x12
 adcs x5, x5, x13
 adcs x6, x6, x14
 adc x11, xzr, xzr
 mov w9, #-1
 mov w10, #-2
 adds x15, x7, #1
 sbcs x16, x8, x9
 sbcs x17, x5, xzr
 adcs x1, x6, x10
 sbcs xzr, x11, xzr
 csel x7, x7, x15, lo
 csel x8, x8, x16, lo
 csel x5, x5, x17, lo
 csel x6, x6, x1, lo
 ldp x3, x4, [x23, #192]
 ldp x9, x10, [x23, #208]
 subs x11, x7, x3
 sbcs x12, x8, x4
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
 stp x11, x12, [x21]
 stp x13, x14, [x21, #16]
 mul x9, x28, x24
 umulh x15, x28, x24
 mul x10, x29, x24
 umulh x16, x29, x24
 mul x11, x30, x24
 umulh x17, x30, x24
 mul x12, x18, x24
 umulh x1, x18, x24
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adc x13, xzr, x1
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x28, x25
 adcs x10, x11, x16
 mul x16, x29, x25
 adcs x11, x12, x17
 mul x17, x30, x25
 adcs x12, x13, x1
 mul x1, x18, x25
 adc x13, xzr, xzr
 adds x9, x9, x15
 umulh x15, x28, x25
 adcs x10, x10, x16
 umulh x16, x29, x25
 adcs x11, x11, x17
 umulh x17, x30, x25
 adcs x12, x12, x1
 umulh x1, x18, x25
 adc x13, x13, xzr
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x28, x26
 adcs x10, x11, x16
 mul x16, x29, x26
 adcs x11, x12, x17
 mul x17, x30, x26
 adcs x12, x13, x1
 mul x1, x18, x26
 adc x13, x14, xzr
 adds x9, x9, x15
 umulh x15, x28, x26
 adcs x10, x10, x16
 umulh x16, x29, x26
 adcs x11, x11, x17
 umulh x17, x30, x26
 adcs x12, x12, x1
 umulh x1, x18, x26
 adc x13, x13, xzr
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 mul x15, x28, x27
 adcs x10, x11, x16
 mul x16, x29, x27
 adcs x11, x12, x17
 mul x17, x30, x27
 adcs x12, x13, x1
 mul x1, x18, x27
 adc x13, x14, xzr
 adds x9, x9, x15
 umulh x15, x28, x27
 adcs x10, x10, x16
 umulh x16, x29, x27
 adcs x11, x11, x17
 umulh x17, x30, x27
 adcs x12, x12, x1
 umulh x1, x18, x27
 adc x13, x13, xzr
 adds x10, x10, x15
 adcs x11, x11, x16
 lsr x16, x9, #32
 adcs x12, x12, x17
 adcs x13, x13, x1
 adc x14, xzr, xzr
 subs x17, x9, x9, lsl #32
 sbc x1, x9, x16
 adds x9, x10, x9, lsl #32
 adcs x10, x11, x16
 adcs x11, x12, x17
 adcs x12, x13, x1
 adc x13, x14, xzr
 mov w28, #-1
 mov w29, #-2
 adds x15, x9, #1
 sbcs x16, x10, x28
 sbcs x17, x11, xzr
 adcs x1, x12, x29
 sbcs xzr, x13, xzr
 csel x9, x9, x15, lo
 csel x10, x10, x16, lo
 csel x11, x11, x17, lo
 csel x12, x12, x1, lo
 mul x8, x10, x9
 umulh x16, x10, x9
 mul x3, x11, x9
 umulh x17, x11, x9
 mul x4, x12, x9
 umulh x5, x12, x9
 adds x3, x3, x16
 mul x15, x11, x10
 umulh x16, x11, x10
 adcs x4, x4, x17
 mul x17, x12, x10
 umulh x1, x12, x10
 adc x5, x5, xzr
 mul x6, x12, x11
 umulh x13, x12, x11
 adds x16, x16, x17
 mul x7, x9, x9
 adc x17, x1, xzr
 adds x4, x4, x15
 umulh x9, x9, x9
 adcs x5, x5, x16
 mul x16, x10, x10
 adcs x6, x6, x17
 umulh x10, x10, x10
 adc x13, x13, xzr
 adds x8, x8, x8
 mul x17, x11, x11
 adcs x3, x3, x3
 umulh x11, x11, x11
 adcs x4, x4, x4
 mul x1, x12, x12
 adcs x5, x5, x5
 umulh x12, x12, x12
 adcs x6, x6, x6
 adcs x13, x13, x13
 adc x14, xzr, xzr
 adds x8, x8, x9
 adcs x3, x3, x16
 adcs x4, x4, x10
 adcs x5, x5, x17
 adcs x6, x6, x11
 adcs x13, x13, x1
 lsr x16, x7, #32
 adc x14, x14, x12
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x3, x16
 adcs x3, x4, x17
 lsr x16, x7, #32
 adc x4, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x3, x16
 adcs x3, x4, x17
 lsr x16, x7, #32
 adc x4, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x3, x16
 adcs x3, x4, x17
 lsr x16, x7, #32
 adc x4, x1, xzr
 subs x17, x7, x7, lsl #32
 sbc x1, x7, x16
 adds x7, x8, x7, lsl #32
 adcs x8, x3, x16
 adcs x3, x4, x17
 adc x4, x1, xzr
 adds x7, x7, x5
 adcs x8, x8, x6
 adcs x3, x3, x13
 adcs x4, x4, x14
 adc x5, xzr, xzr
 mov w9, #-1
 mov w10, #-2
 adds x15, x7, #1
 sbcs x16, x8, x9
 sbcs x17, x3, xzr
 adcs x1, x4, x10
 sbcs xzr, x5, xzr
 csel x7, x7, x15, lo
 csel x8, x8, x16, lo
 csel x3, x3, x17, lo
 csel x4, x4, x1, lo
 ldp x5, x6, [x23, #192]
 ldp x9, x10, [x23, #208]
 subs x11, x7, x5
 sbcs x12, x8, x6
 sbcs x13, x3, x9
 sbcs x14, x4, x10
 sbc x15, xzr, xzr
 lsr x16, x15, #32
 lsl x17, x15, #32
 sub x17, x17, x15
 adds x11, x11, x15
 adcs x12, x12, x16
 adcs x13, x13, xzr
 adc x14, x14, x17
 stp x11, x12, [x21, #32]
 stp x13, x14, [x21, #48]
 sub x20,x20,#64
 sub x21,x21,#64
 subs x19,x19,#1
 b.ne pb_loop
pb_done:
 ldr x30,[sp],#16
 ret

// Custom internal ABI: x0 count, x1 first/last Y, x2 first/last batch,
// x3 X point, x4 scratch. All GPRs may be clobbered except SP and LR.
// Batch has 32-byte Montgomery-one prefix before its first 64-byte entry.
.globl pair_fwd_lo
pair_fwd_lo:
 str x30,[sp,#-16]!
 mov x19,x0
 mov x20,x1
 mov x21,x2
 mov x22,x3
 sub x23,x2,#32
 cbz x19,pf_success_lo
 mov w0,#-1
 mov w2,#-2
pf_loop_lo:
 // inline k_sub; redundant/absolute research
 ldp x3, x4, [x20], #64
 ldp x5, x6, [x20, #-48]
 ldp x7, x8, [x22]
 ldp x9, x10, [x22, #16]
 subs x11, x7, x3
 sbcs x12, x8, x4
 sbcs x13, x9, x5
 sbcs x14, x10, x6
abs_store_fw0_lo:
 stp x11,x12,[x21]
 stp x13,x14,[x21, #16]
 orr x9,x11,x12
 orr x10,x13,x14
 orr x9,x9,x10
 cbz x9,pf_fail_lo
 // inline k_mul; redundant/absolute research
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
 cbz x26,k_mul_store_fw1_lo
 adds x9,x9,#1
 sbcs x10,x10,x0
 sbcs x24,x24,xzr
 adc x25,x25,x2
k_mul_store_fw1_lo:
 stp x9,x10,[x21, #32]
 stp x24,x25,[x21, #48]
 add x23,x21,#32
 add x21,x21,#64
 subs x19,x19,#1
 b.ne pf_loop_lo
pf_success_lo:
 mov x0,#1
 b pf_done_lo
pf_fail_lo:
 mov x0,#0
pf_done_lo:
 ldr x30,[sp],#16
 ret
// Custom internal ABI: x0 count, x1 first/last Y, x2 first/last batch,
// x3 X point, x4 scratch. All GPRs may be clobbered except SP and LR.
// Batch has 32-byte Montgomery-one prefix before its first 64-byte entry.
.globl pair_fwd_hi
pair_fwd_hi:
 str x30,[sp,#-16]!
 mov x19,x0
 mov x20,x1
 mov x21,x2
 mov x22,x3
 sub x23,x2,#32
 cbz x19,pf_success_hi
 mov w0,#-1
 mov w2,#-2
pf_loop_hi:
 // inline k_sub; redundant/absolute research
 ldp x3, x4, [x20], #64
 ldp x5, x6, [x20, #-48]
 ldp x7, x8, [x22]
 ldp x9, x10, [x22, #16]
 subs x11, x3, x7
 sbcs x12, x4, x8
 sbcs x13, x5, x9
 sbcs x14, x6, x10
abs_store_fw0_hi:
 stp x11,x12,[x21]
 stp x13,x14,[x21, #16]
 orr x9,x11,x12
 orr x10,x13,x14
 orr x9,x9,x10
 cbz x9,pf_fail_hi
 // inline k_mul; redundant/absolute research
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
 cbz x26,k_mul_store_fw1_hi
 adds x9,x9,#1
 sbcs x10,x10,x0
 sbcs x24,x24,xzr
 adc x25,x25,x2
k_mul_store_fw1_hi:
 stp x9,x10,[x21, #32]
 stp x24,x25,[x21, #48]
 add x23,x21,#32
 add x21,x21,#64
 subs x19,x19,#1
 b.ne pf_loop_hi
pf_success_hi:
 mov x0,#1
 b pf_done_hi
pf_fail_hi:
 mov x0,#0
pf_done_hi:
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
 mov w0,#-1
 mov w2,#-2
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
 cbz x13,k_mul_store_0
 adds x24,x24,#1
 sbcs x25,x25,x0
 sbcs x26,x26,xzr
 adc x27,x27,x2
k_mul_store_0:
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
 cbz x13,k_mul_store_1
 adds x9,x9,#1
 sbcs x10,x10,x0
 sbcs x11,x11,xzr
 adc x12,x12,x2
k_mul_store_1:
 stp x9,x10,[x23]
 stp x11,x12,[x23, #16]
 ldp x3, x4, [x20]
 ldp x5, x6, [x20, #16]
 ldp x7, x8, [x22]
 ldp x9, x10, [x22, #16]
 adds x11, x3, x7
 adcs x12, x4, x8
 adcs x13, x5, x9
 adcs x14, x6, x10
 adc x15, xzr, xzr
 cbz x15,k_add_store_2
 adds x11,x11,#1
 sbcs x12,x12,x0
 sbcs x13,x13,xzr
 adc x14,x14,x2
k_add_store_2:
 stp x11,x12,[x23, #192]
 stp x13,x14,[x23, #208]
 ldp x3, x4, [x20, #32]
 ldp x5, x6, [x20, #48]
 ldp x7, x8, [x22, #32]
 ldp x9, x10, [x22, #48]
 subs x11, x3, x7
 sbcs x12, x4, x8
 sbcs x13, x5, x9
 sbcs x14, x6, x10
 b.hs abs_store_3
 subs x11,xzr,x11
 sbcs x12,xzr,x12
 sbcs x13,xzr,x13
 sbc x14,xzr,x14
abs_store_3:
 adds x28, x3, x7
 adcs x29, x4, x8
 adcs x30, x5, x9
 adcs x18, x6, x10
k_add_store_4:
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
 cbz x5,k_mul_store_5
 adds x9,x9,#1
 sbcs x10,x10,x0
 sbcs x3,x3,xzr
 adc x4,x4,x2
k_mul_store_5:
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
 b.cs fs_high_6
 ldp x12,x13,[x23, #192]
 ldp x14,x15,[x23, #208]
 subs x7,x7,x12
 sbcs x8,x8,x13
 sbcs x5,x5,x14
 sbcs x6,x6,x15
 b.cs fs_store_6
fs_negative_6:
 subs x7,x7,#1
 adcs x8,x8,x0
 adcs x5,x5,xzr
 sbc x6,x6,x2
 b fs_store_6
fs_high_6:
 ldp x12,x13,[x23, #192]
 ldp x14,x15,[x23, #208]
 subs x7,x7,x12
 sbcs x8,x8,x13
 sbcs x5,x5,x14
 sbcs x6,x6,x15
 b.cc fs_store_6
fs_positive_6:
 adds x7,x7,#1
 sbcs x8,x8,x0
 sbcs x5,x5,xzr
 adc x6,x6,x2
fs_store_6:
 stp x7,x8,[x21]
 stp x5,x6,[x21, #16]
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
 cbz x13,k_mul_store_8
 adds x9,x9,#1
 sbcs x10,x10,x0
 sbcs x11,x11,xzr
 adc x12,x12,x2
k_mul_store_8:
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
 b.cs fs_high_9
 ldp x6,x13,[x23, #192]
 ldp x14,x15,[x23, #208]
 subs x7,x7,x6
 sbcs x8,x8,x13
 sbcs x3,x3,x14
 sbcs x4,x4,x15
 b.cs fs_store_9
fs_negative_9:
 subs x7,x7,#1
 adcs x8,x8,x0
 adcs x3,x3,xzr
 sbc x4,x4,x2
 b fs_store_9
fs_high_9:
 ldp x6,x13,[x23, #192]
 ldp x14,x15,[x23, #208]
 subs x7,x7,x6
 sbcs x8,x8,x13
 sbcs x3,x3,x14
 sbcs x4,x4,x15
 b.cc fs_store_9
fs_positive_9:
 adds x7,x7,#1
 sbcs x8,x8,x0
 sbcs x3,x3,xzr
 adc x4,x4,x2
fs_store_9:
 stp x7,x8,[x21, #32]
 stp x3,x4,[x21, #48]
 sub x20,x20,#64
 sub x21,x21,#64
 subs x19,x19,#1
 b.ne pb_loop
pb_done:
 ldr x30,[sp],#16
 ret

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
 add x0,x21,#0
 add x1,x20,#0
 add x2,x22,#0
 bl k_sub
 ldp x9,x10,[x21]
 ldp x11,x12,[x21,#16]
 orr x9,x9,x10
 orr x11,x11,x12
 orr x9,x9,x11
 cbz x9,pf_fail
 add x0,x21,#32
 add x1,x23,#0
 add x2,x21,#0
 bl k_mul
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
 add x0,x23,#32
 add x1,x23,#0
 sub x2,x21,#32
 bl k_mul
 add x0,x23,#0
 add x1,x23,#0
 add x2,x21,#0
 bl k_mul
 add x0,x23,#64
 add x1,x20,#32
 add x2,x22,#32
 bl k_add
 add x0,x23,#96
 add x1,x20,#32
 add x2,x22,#32
 bl k_sub
 add x0,x23,#128
 add x1,x23,#96
 add x2,x23,#32
 bl k_mul
 add x0,x23,#160
 add x1,x23,#64
 add x2,x23,#32
 bl k_mul
 add x0,x23,#192
 add x1,x20,#0
 add x2,x22,#0
 bl k_add
 add x0,x23,#128
 add x1,x23,#128
 bl k_sqr
 add x0,x23,#160
 add x1,x23,#160
 bl k_sqr
 add x0,x21,#0
 add x1,x23,#128
 add x2,x23,#192
 bl k_sub
 add x0,x21,#32
 add x1,x23,#160
 add x2,x23,#192
 bl k_sub
 sub x20,x20,#64
 sub x21,x21,#64
 subs x19,x19,#1
 b.ne pb_loop
pb_done:
 ldr x30,[sp],#16
 ret

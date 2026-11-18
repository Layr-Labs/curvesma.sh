.globl pair_keys
pair_keys:
 cbz x0,pk_done
 mov x19,x0
 mov x20,x1
 mov x21,x2
 add x22,x2,#64
 mov x23,x3
 mov x28,x5
 mov x25,x4
 mov x27,x6
 mov x14,#4096
 mov x16,#128
 mov x26,#32
pk_loop:
 ldr x24,[x28],#72
 ldp x5,x6,[x20,#0]
 ldp x7,x8,[x20,#16]
 and x5,x5,#0xfffffffff
 and x6,x6,#0xfffffffff
 add x9,x22,x5,lsl #3
 ldadd x26,x11,[x9]
 orr x12,x6,x24
 stp x7,x8,[x11]
 stp x23,x12,[x11,#16]
 tbz x11,#12,pk_skip0
 sub x11,x11,#4064
 str x11,[x9]
 add x15,x27,x5,lsl #3
 ldadd x14,x13,[x15]
 cmp x13,x25
 b.hs pk_skip0
 stp x5,x13,[x21]
 stp x11,x16,[x21,#16]
pk_skip0:
 mov x24,x24
 ldp x5,x6,[x20,#32]
 ldp x7,x8,[x20,#48]
 and x5,x5,#0xfffffffff
 and x6,x6,#0xfffffffff
 add x9,x22,x5,lsl #3
 ldadd x26,x11,[x9]
 orr x12,x6,x24
 stp x7,x8,[x11]
 stp x23,x12,[x11,#16]
 tbz x11,#12,pk_skip1
 sub x11,x11,#4064
 str x11,[x9]
 add x15,x27,x5,lsl #3
 ldadd x14,x13,[x15]
 cmp x13,x25
 b.hs pk_skip1
 stp x5,x13,[x21]
 stp x11,x16,[x21,#16]
pk_skip1:
 add x20,x20,#64
 subs x19,x19,#1
 b.ne pk_loop
pk_done:
 ret
.macro DP
 ld64b x6,[x3]
 add x3,x3,#64
 st64b x6,[x2]
 add x2,x2,#64
.endm
.macro DP7
 DP
 DP
 DP
 DP
 DP
 DP
 DP
.endm
.macro DS
 ldp x6,x7,[x3],#16
 ldp x8,x9,[x3],#16
 stp x6,x7,[x2],#16
 stp x8,x9,[x2],#16
.endm
.globl block_dlv
block_dlv:
 ldp x1,x2,[x0]
 ldp x3,x4,[x0,#16]
 DS
 DP7
 DP7
 DP7
 DP7
 DP7
 DP7
 DP7
 DP7
 DP7
 DS
 ldr x9,[x0,#32]
 mov x10,#1
 stlr x10,[x9]
 ret
.globl block_finish
block_finish:
 mov x9,#0
 mov x13,#128
 cbz x0,bf_done
bf_loop:
 ldr x10,[x1],#8
 ldr x11,[x2],#8
 sub x12,x10,x3
 cmp x11,x5
 b.hs bf_no_copy
 cbz x12,bf_no_copy
 stp x9,x11,[x4]
 stp x3,x13,[x4,#16]
bf_no_copy:
 add x11,x11,x12
 cmp x11,x5
 csel x11,x11,x5,lo
 stp x9,x8,[x7]
 stp x11,xzr,[x7,#16]
 stp xzr,xzr,[x7,#32]
 add x3,x3,#8192
 add x9,x9,#1
 subs x0,x0,#1
 b.ne bf_loop
bf_done:
 ret
.globl block_init
block_init:
 cbz x0,bi_done
bi_outer:
 str x3,[x1],#8
 str x4,[x2],#8
 mov x5,#256
 mov x6,x3
bi_zero:
 stp xzr,xzr,[x6],#16
 subs x5,x5,#1
 b.ne bi_zero
 add x3,x3,#8192
 subs x0,x0,#1
 b.ne bi_outer
bi_done:
 ret

// x0 pair count, x1 paired-x records, x2 column port, x3 packed row/k, x4 slot end.
.globl pair_keys
pair_keys:
 cbz x0,pk_done
 mov x19,x0
 mov x20,x1
 mov x21,x2
 add x22,x2,#64
 mov x23,x3
 mov x24,#0
 mov x25,x4
 mov x26,#32
 mov x27,#0x1000000000
pk_loop:
 ldp x5,x6,[x20,#0]
 ldp x7,x8,[x20,#16]
 and x5,x5,#0xfffffffff
 and x6,x6,#0xfffffffff
 add x9,x22,x5,lsl #3
 ldadd x26,x11,[x9]
 cmp x11,x25
 b.hs pk_skip0
 orr x12,x6,x24
 stp x5,x11,[x21]
 stp x7,x8,[x21,#16]
 stp x23,x12,[x21,#32]
pk_skip0:
 orr x24,x24,#0x1000000000
 ldp x5,x6,[x20,#32]
 ldp x7,x8,[x20,#48]
 and x5,x5,#0xfffffffff
 and x6,x6,#0xfffffffff
 add x9,x22,x5,lsl #3
 ldadd x26,x11,[x9]
 cmp x11,x25
 b.hs pk_skip1
 orr x12,x6,x24
 stp x5,x11,[x21]
 stp x7,x8,[x21,#16]
 stp x23,x12,[x21,#32]
pk_skip1:
 add x24,x24,x27
 add x20,x20,#64
 subs x19,x19,#1
 b.ne pk_loop
pk_done:
 ret

// Fused port-publish emitter: single-hop direct-home delivery.
// x0 pair count, x1 paired-x records (64B stride: two 256-bit keys),
// x2 home port (48B device), x3 packed row/k W2, x5 Y-ID records (72B stride),
// x17 source column j. Publishes one 6-word message [d1,d2,K2,K3,w2,w3] per
// record; the completing store to +40 waits once per the blocking-port
// contract. No slots, windows, blocks, cursors or phases.
.globl pair_keys_direct
pair_keys_direct:
 cbz x0,pd_done
 mov x19,x0
 mov x20,x1
 mov x21,x2
 mov x23,x3
 mov x28,x5
pd_loop:
 ldr x24,[x28],#72
 ldp x5,x6,[x20,#0]
 ldp x7,x8,[x20,#16]
 and x5,x5,#0x7ffffffff
 and x6,x6,#0x7ffffffff
 orr x12,x17,x24
 stp x5,x6,[x21]
 stp x7,x8,[x21,#16]
 stp x23,x12,[x21,#32]
 orr x24,x24,#0x800000000
 ldp x5,x6,[x20,#32]
 ldp x7,x8,[x20,#48]
 and x5,x5,#0x7ffffffff
 and x6,x6,#0x7ffffffff
 orr x12,x17,x24
 stp x5,x6,[x21]
 stp x7,x8,[x21,#16]
 stp x23,x12,[x21,#32]
 add x20,x20,#64
 subs x19,x19,#1
 b.ne pd_loop
pd_done:
 ret

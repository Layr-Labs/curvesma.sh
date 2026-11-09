.globl send64
send64:
 cbz x0,send_done
send_loop:
 ldr x4,[x1,x3]
 and x4,x4,#0xfffffffff
 stp x4,xzr,[x2]
 ldp x4,x5,[x1]
 stp x4,x5,[x2,#16]
 ldp x4,x5,[x1,#16]
 stp x4,x5,[x2,#32]
 ldp x4,x5,[x1,#32]
 stp x4,x5,[x2,#48]
 ldp x4,x5,[x1,#48]
 stp x4,x5,[x2,#64]
 add x1,x1,#64
 subs x0,x0,#1
 b.ne send_loop
send_done:
 ret
.globl append64
append64:
 ldp x4,x5,[x0,#16]
 ldadd x18,x12,[x1]
 cmp x12,x2
 b.hs append_full
 ldp x6,x7,[x0,#32]
 ldp x8,x9,[x0,#48]
 ldp x10,x11,[x0,#64]
 stp x4,x5,[x12]
 stp x6,x7,[x12,#16]
 stp x8,x9,[x12,#32]
 stp x10,x11,[x12,#48]
 stlr xzr,[x3]
 ret
append_full:
 stadd x19,[x1]
 stlr xzr,[x3]
 ret
.globl home64
home64:
 ldp x4,x5,[x0,#16]
 lsr x6,x4,#36
 ubfx x7,x5,#36,#1
 ubfx x8,x5,#37,#3
 add x6,x7,x6,lsl #1
 add x6,x8,x6,lsl #3
 add x12,x1,x6,lsl #6
 stp x4,x5,[x12]
 ldp x6,x7,[x0,#32]
 stp x6,x7,[x12,#16]
 ldp x8,x9,[x0,#48]
 stp x8,x9,[x12,#32]
 ldp x10,x11,[x0,#64]
 stp x10,x11,[x12,#48]
 stlr xzr,[x3]
 ret
.globl clear_flags
clear_flags:
 cbz x0,clear_done
clear_loop:
 str xzr,[x1],#64
 subs x0,x0,#1
 b.ne clear_loop
clear_done:
 ret
.globl complete
complete:
 cbz x0,complete_yes
complete_loop:
 ldr x4,[x1],#64
 cbz x4,complete_no
 subs x0,x0,#1
 b.ne complete_loop
complete_yes:
 mov x0,#1
 ret
complete_no:
 mov x0,#0
 ret
.globl gather224
gather224:
 ldp x6,x7,[x0,#0]
 ldp x8,x9,[x0,#16]
 stp x6,x7,[x5,#0]
 stp x8,x9,[x5,#16]
 ldp x6,x7,[x0,#32]
 ldp x8,x9,[x0,#48]
 stp x6,x7,[x5,#32]
 stp x8,x9,[x5,#48]
 ldp x6,x7,[x1,#0]
 ldp x8,x9,[x1,#16]
 stp x6,x7,[x5,#64]
 stp x8,x9,[x5,#80]
 ldp x6,x7,[x1,#32]
 ldp x8,x9,[x1,#48]
 stp x6,x7,[x5,#96]
 stp x8,x9,[x5,#112]
 ldp x6,x7,[x2,#0]
 ldp x8,x9,[x2,#16]
 stp x6,x7,[x5,#128]
 stp x8,x9,[x5,#144]
 ldp x6,x7,[x3,#0]
 ldp x8,x9,[x3,#16]
 stp x6,x7,[x5,#160]
 stp x8,x9,[x5,#176]
 ldp x6,x7,[x4,#0]
 ldp x8,x9,[x4,#16]
 stp x6,x7,[x5,#192]
 stp x8,x9,[x5,#208]
 stp x14,xzr,[x5,#224]
 stp xzr,xzr,[x5,#240]
 ret
.globl build_replies
build_replies:
 mov x17,#1
 mov x18,#0x2000000000
 mov x19,#8
reply_loop:
 stp x0,x1,[x3]
 stp x17,xzr,[x3,#16]
 ldp x4,x5,[x2],#16
 ldp x6,x7,[x2],#16
 stp x4,x5,[x3,#32]
 stp x6,x7,[x3,#48]
 add x1,x1,x18
 add x3,x3,#64
 subs x19,x19,#1
 b.ne reply_loop
 ret
.globl make_requests
make_requests:
 mov x5,#0
 cbz x0,requests_done
requests_loop:
 ldp x6,x7,[x1],#16
 ldp x8,x9,[x6]
 stp x8,x9,[x2]
 orr x10,x3,x5,lsl #36
 stp x10,x4,[x2,#16]
 ldr x13,[x6,#-8]
 stp xzr,x13,[x2,#32]
 stp xzr,xzr,[x2,#48]
 ldp x8,x9,[x7]
 stp x8,x9,[x2,#64]
 orr x11,x4,#0x1000000000
 stp x10,x11,[x2,#80]
 ldr x13,[x7,#-8]
 stp xzr,x13,[x2,#96]
 stp xzr,xzr,[x2,#112]
 add x2,x2,#128
 add x5,x5,#1
 subs x0,x0,#1
 b.ne requests_loop
requests_done:
 ret
.globl publish64
publish64:
 ldp x1,x2,[x0]
 ldp x3,x4,[x0,#16]
 ldp x5,x6,[x0,#32]
 ldp x7,x8,[x0,#48]
 ldp x9,x10,[x0,#64]
 stp x3,x4,[x2]
 stp x5,x6,[x2,#16]
 stp x7,x8,[x2,#32]
 str x9,[x2,#48]
 add x11,x2,#2048
 stlr x10,[x11]
 ret
.globl send_completion
send_completion:
 stp x0,x1,[x3]
 ldp x4,x5,[x2]
 stp x4,x5,[x3,#16]
 ldp x4,x5,[x2,#16]
 stp x4,x5,[x3,#32]
 ldp x4,x5,[x2,#32]
 stp x4,x5,[x3,#48]
 ldp x4,x5,[x2,#48]
 stp x4,x5,[x3,#64]
 ret

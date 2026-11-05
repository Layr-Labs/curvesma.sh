.globl range_merge
range_merge:
merge_loop:
 cbz x1,merge_only_right
 cbz x3,merge_left
 add x7,x0,x5
 add x8,x2,x5
 ldp x9,x10,[x7]
 ldp x11,x12,[x7,#16]
 ldp x13,x14,[x8]
 ldp x15,x16,[x8,#16]
 cmp x12,x16
 b.lo merge_left
 b.ne merge_right
 cmp x11,x15
 b.lo merge_left
 b.ne merge_right
 cmp x10,x14
 b.lo merge_left
 b.ne merge_right
 cmp x9,x13
 b.lo merge_left
 b.ne merge_right
 ldr x17,[x0,#64]
 ldr x7,[x2,#64]
 cmp x17,x7
 b.lo merge_left
 b.ne merge_right
 b merge_left
merge_only_right:
 cbz x3,merge_done
merge_right:
 mov x6,x2
 add x2,x2,#72
 sub x3,x3,#1
 b merge_copy
merge_left:
 mov x6,x0
 add x0,x0,#72
 sub x1,x1,#1
merge_copy:
 ldp x7,x8,[x6,#0]
 stp x7,x8,[x4,#0]
 ldp x7,x8,[x6,#16]
 stp x7,x8,[x4,#16]
 ldp x7,x8,[x6,#32]
 stp x7,x8,[x4,#32]
 ldp x7,x8,[x6,#48]
 stp x7,x8,[x4,#48]
 ldr x7,[x6,#64]
 str x7,[x4,#64]
 add x4,x4,#72
 b merge_loop
merge_done:
 mov x0,x4
 ret

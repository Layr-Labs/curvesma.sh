.globl range_lower_bound
range_lower_bound:
 mov x16,x0
 mov x12,#0
 mov x13,x1
 mov x17,#72
 ldp x4,x5,[x2]
 ldp x6,x7,[x2,#16]
range_search_loop:
 cmp x12,x13
 b.hs range_search_done
 add x14,x12,x13
 lsr x14,x14,#1
 mul x15,x14,x17
 add x15,x16,x15
 ldp x8,x9,[x15]
 ldp x10,x11,[x15,#16]
 subs xzr,x8,x4
 sbcs xzr,x9,x5
 sbcs xzr,x10,x6
 sbcs xzr,x11,x7
 b.hs range_search_ge
 add x12,x14,#1
 b range_search_loop
range_search_ge:
 mov x13,x14
 b range_search_loop
range_search_done:
 mov x0,x12
 ret

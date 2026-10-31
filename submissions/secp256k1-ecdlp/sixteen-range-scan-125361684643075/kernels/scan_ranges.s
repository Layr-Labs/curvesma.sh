// Sixteen-comparison scanner with branch-tree prefix dispatch.
// Same K2 matches within fixed d4 ranges; no cross-range match.
.globl scan_ranges
scan_ranges:
 cbz x0,s16_done
s16_bucket:
 ldr x18,[x2],#8
 cmp x1,x18
 b.hs s16_next
 ldr x16,[x1],#8
 ldr x5,[x16]
 sub x17,x18,x1
 lsr x17,x17,#3
 and x17,x17,#15
 mov x9,x16
 mov x7,x5
 tbz x17,#3,s16_dispatch_0_3
 tbz x17,#2,s16_dispatch_8_2
 tbz x17,#1,s16_dispatch_12_1
 tbz x17,#0,s16_prefix2
 b s16_prefix1
s16_dispatch_12_1:
 tbz x17,#0,s16_prefix4
 b s16_prefix3
s16_dispatch_8_2:
 tbz x17,#1,s16_dispatch_8_1
 tbz x17,#0,s16_prefix6
 b s16_prefix5
s16_dispatch_8_1:
 tbz x17,#0,s16_prefix8
 b s16_prefix7
s16_dispatch_0_3:
 tbz x17,#2,s16_dispatch_0_2
 tbz x17,#1,s16_dispatch_4_1
 tbz x17,#0,s16_prefix10
 b s16_prefix9
s16_dispatch_4_1:
 tbz x17,#0,s16_prefix12
 b s16_prefix11
s16_dispatch_0_2:
 tbz x17,#1,s16_dispatch_0_1
 tbz x17,#0,s16_prefix14
 b s16_prefix13
s16_dispatch_0_1:
 tbz x17,#0,s16_norm
 b s16_prefix15
s16_prefix1:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit1
s16_pcont1:
s16_prefix2:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit2
s16_pcont2:
s16_prefix3:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit3
s16_pcont3:
s16_prefix4:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit4
s16_pcont4:
s16_prefix5:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit5
s16_pcont5:
s16_prefix6:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit6
s16_pcont6:
s16_prefix7:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit7
s16_pcont7:
s16_prefix8:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit8
s16_pcont8:
s16_prefix9:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit9
s16_pcont9:
s16_prefix10:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit10
s16_pcont10:
s16_prefix11:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit11
s16_pcont11:
s16_prefix12:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit12
s16_pcont12:
s16_prefix13:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit13
s16_pcont13:
s16_prefix14:
 ldr x16,[x1],#8
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_phit14
s16_pcont14:
s16_prefix15:
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_phit15
s16_pcont15:
s16_norm:
 mov x16,x9
 mov x5,x7
 cmp x1,x18
 b.hs s16_next
s16_loop:
 ldp x9,x10,[x1],#16
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_hit0
s16_cont0:
 ldr x5,[x10]
 cmp x7,x5
 b.eq s16_hit1
s16_cont1:
 ldp x11,x12,[x1],#16
 ldr x7,[x11]
 cmp x5,x7
 b.eq s16_hit2
s16_cont2:
 ldr x5,[x12]
 cmp x7,x5
 b.eq s16_hit3
s16_cont3:
 ldp x13,x14,[x1],#16
 ldr x7,[x13]
 cmp x5,x7
 b.eq s16_hit4
s16_cont4:
 ldr x5,[x14]
 cmp x7,x5
 b.eq s16_hit5
s16_cont5:
 ldp x15,x16,[x1],#16
 ldr x7,[x15]
 cmp x5,x7
 b.eq s16_hit6
s16_cont6:
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_hit7
s16_cont7:
 ldp x9,x10,[x1],#16
 ldr x7,[x9]
 cmp x5,x7
 b.eq s16_hit8
s16_cont8:
 ldr x5,[x10]
 cmp x7,x5
 b.eq s16_hit9
s16_cont9:
 ldp x11,x12,[x1],#16
 ldr x7,[x11]
 cmp x5,x7
 b.eq s16_hit10
s16_cont10:
 ldr x5,[x12]
 cmp x7,x5
 b.eq s16_hit11
s16_cont11:
 ldp x13,x14,[x1],#16
 ldr x7,[x13]
 cmp x5,x7
 b.eq s16_hit12
s16_cont12:
 ldr x5,[x14]
 cmp x7,x5
 b.eq s16_hit13
s16_cont13:
 ldp x15,x16,[x1],#16
 ldr x7,[x15]
 cmp x5,x7
 b.eq s16_hit14
s16_cont14:
 ldr x5,[x16]
 cmp x7,x5
 b.eq s16_hit15
s16_cont15:
 cmp x1,x18
 b.lo s16_loop
s16_next:
 subs x0,x0,#1
 b.ne s16_bucket
s16_done:
 ret
s16_phit1:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont1
s16_phit2:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont2
s16_phit3:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont3
s16_phit4:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont4
s16_phit5:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont5
s16_phit6:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont6
s16_phit7:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont7
s16_phit8:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont8
s16_phit9:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont9
s16_phit10:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont10
s16_phit11:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont11
s16_phit12:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont12
s16_phit13:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont13
s16_phit14:
 cmp x3,x4
 b.hs s16_done
 stp x9,x16,[x3],#16
 b s16_pcont14
s16_phit15:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_pcont15
s16_hit0:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_cont0
s16_hit1:
 cmp x3,x4
 b.hs s16_done
 stp x9,x10,[x3],#16
 b s16_cont1
s16_hit2:
 cmp x3,x4
 b.hs s16_done
 stp x10,x11,[x3],#16
 b s16_cont2
s16_hit3:
 cmp x3,x4
 b.hs s16_done
 stp x11,x12,[x3],#16
 b s16_cont3
s16_hit4:
 cmp x3,x4
 b.hs s16_done
 stp x12,x13,[x3],#16
 b s16_cont4
s16_hit5:
 cmp x3,x4
 b.hs s16_done
 stp x13,x14,[x3],#16
 b s16_cont5
s16_hit6:
 cmp x3,x4
 b.hs s16_done
 stp x14,x15,[x3],#16
 b s16_cont6
s16_hit7:
 cmp x3,x4
 b.hs s16_done
 stp x15,x16,[x3],#16
 b s16_cont7
s16_hit8:
 cmp x3,x4
 b.hs s16_done
 stp x16,x9,[x3],#16
 b s16_cont8
s16_hit9:
 cmp x3,x4
 b.hs s16_done
 stp x9,x10,[x3],#16
 b s16_cont9
s16_hit10:
 cmp x3,x4
 b.hs s16_done
 stp x10,x11,[x3],#16
 b s16_cont10
s16_hit11:
 cmp x3,x4
 b.hs s16_done
 stp x11,x12,[x3],#16
 b s16_cont11
s16_hit12:
 cmp x3,x4
 b.hs s16_done
 stp x12,x13,[x3],#16
 b s16_cont12
s16_hit13:
 cmp x3,x4
 b.hs s16_done
 stp x13,x14,[x3],#16
 b s16_cont13
s16_hit14:
 cmp x3,x4
 b.hs s16_done
 stp x14,x15,[x3],#16
 b s16_cont14
s16_hit15:
 cmp x3,x4
 b.hs s16_done
 stp x15,x16,[x3],#16
 b s16_cont15

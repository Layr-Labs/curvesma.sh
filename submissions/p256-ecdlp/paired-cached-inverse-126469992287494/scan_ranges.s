.globl scan_ranges
scan_ranges:
 cbz x0,pc_done
pc_bucket:
 ldr x18,[x2],#8
 cmp x1,x18
 b.hs pc_next
 ldr x16,[x1],#8
 ldr x5,[x16]
 sub x17,x18,x1
 lsr x17,x17,#3
 tbz x17,#0,pc_even
 ldr x9,[x1],#8
 ldr x7,[x9]
 cmp x5,x7
 b.eq pc_odd_hit
pc_odd_cont:
 mov x16,x9
 mov x5,x7
pc_even:
 lsr x17,x17,#1
 and x17,x17,#7
 mov x7,x5
 mov x8,x5
 mov x10,x16
 mov x12,x16
 mov x14,x16
 tbz x17,#2,pc_dlow
 tbz x17,#1,pc_d4
 tbz x17,#0,pc_pair2
 b pc_pair1
pc_d4:
 tbz x17,#0,pc_pair4
 b pc_pair3
pc_dlow:
 tbz x17,#1,pc_d0
 tbz x17,#0,pc_pair6
 b pc_pair5
pc_d0:
 tbz x17,#0,pc_finish
 b pc_pair7
pc_pair0:
 ldp x9,x10,[x1],#16
 ldr x7,[x9]
 ldr x8,[x10]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit0
pc_cont0:
pc_pair1:
 ldp x11,x12,[x1],#16
 ldr x5,[x11]
 ldr x7,[x12]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit1
pc_cont1:
pc_pair2:
 ldp x13,x14,[x1],#16
 ldr x8,[x13]
 ldr x5,[x14]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit2
pc_cont2:
pc_pair3:
 ldp x15,x16,[x1],#16
 ldr x7,[x15]
 ldr x8,[x16]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit3
pc_cont3:
pc_pair4:
 ldp x9,x10,[x1],#16
 ldr x5,[x9]
 ldr x7,[x10]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit4
pc_cont4:
pc_pair5:
 ldp x11,x12,[x1],#16
 ldr x8,[x11]
 ldr x5,[x12]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit5
pc_cont5:
pc_pair6:
 ldp x13,x14,[x1],#16
 ldr x7,[x13]
 ldr x8,[x14]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit6
pc_cont6:
pc_pair7:
 ldp x15,x16,[x1],#16
 ldr x5,[x15]
 ldr x7,[x16]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit7
pc_cont7:
pc_finish:
 mov x5,x7
 cmp x1,x18
 b.lo pc_pair0
pc_next:
 subs x0,x0,#1
 b.ne pc_bucket
pc_done:
 ret
pc_odd_hit:
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
 b pc_odd_cont
pc_hit0:
 cmp x5,x7
 b.ne pc_second0
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second0:
 cmp x7,x8
 b.ne pc_cont0
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont0
pc_hit1:
 cmp x8,x5
 b.ne pc_second1
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second1:
 cmp x5,x7
 b.ne pc_cont1
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont1
pc_hit2:
 cmp x7,x8
 b.ne pc_second2
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second2:
 cmp x8,x5
 b.ne pc_cont2
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont2
pc_hit3:
 cmp x5,x7
 b.ne pc_second3
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second3:
 cmp x7,x8
 b.ne pc_cont3
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont3
pc_hit4:
 cmp x8,x5
 b.ne pc_second4
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second4:
 cmp x5,x7
 b.ne pc_cont4
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont4
pc_hit5:
 cmp x7,x8
 b.ne pc_second5
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second5:
 cmp x8,x5
 b.ne pc_cont5
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont5
pc_hit6:
 cmp x5,x7
 b.ne pc_second6
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second6:
 cmp x7,x8
 b.ne pc_cont6
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont6
pc_hit7:
 cmp x8,x5
 b.ne pc_second7
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second7:
 cmp x5,x7
 b.ne pc_cont7
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont7

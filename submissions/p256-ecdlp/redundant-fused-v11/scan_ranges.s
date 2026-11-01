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
 and x17,x17,#31
 mov x7,x5
 mov x8,x5
 mov x10,x16
 mov x12,x16
 mov x14,x16
 tbz x17,#4,pc_dispatch_4_0
 tbz x17,#3,pc_dispatch_3_16
 tbz x17,#2,pc_dispatch_2_24
 tbz x17,#1,pc_dispatch_1_28
 tbz x17,#0,pc_pair2
 b pc_pair1
pc_dispatch_1_28:
 tbz x17,#0,pc_pair4
 b pc_pair3
pc_dispatch_2_24:
 tbz x17,#1,pc_dispatch_1_24
 tbz x17,#0,pc_pair6
 b pc_pair5
pc_dispatch_1_24:
 tbz x17,#0,pc_pair8
 b pc_pair7
pc_dispatch_3_16:
 tbz x17,#2,pc_dispatch_2_16
 tbz x17,#1,pc_dispatch_1_20
 tbz x17,#0,pc_pair10
 b pc_pair9
pc_dispatch_1_20:
 tbz x17,#0,pc_pair12
 b pc_pair11
pc_dispatch_2_16:
 tbz x17,#1,pc_dispatch_1_16
 tbz x17,#0,pc_pair14
 b pc_pair13
pc_dispatch_1_16:
 tbz x17,#0,pc_pair16
 b pc_pair15
pc_dispatch_4_0:
 tbz x17,#3,pc_dispatch_3_0
 tbz x17,#2,pc_dispatch_2_8
 tbz x17,#1,pc_dispatch_1_12
 tbz x17,#0,pc_pair18
 b pc_pair17
pc_dispatch_1_12:
 tbz x17,#0,pc_pair20
 b pc_pair19
pc_dispatch_2_8:
 tbz x17,#1,pc_dispatch_1_8
 tbz x17,#0,pc_pair22
 b pc_pair21
pc_dispatch_1_8:
 tbz x17,#0,pc_pair24
 b pc_pair23
pc_dispatch_3_0:
 tbz x17,#2,pc_dispatch_2_0
 tbz x17,#1,pc_dispatch_1_4
 tbz x17,#0,pc_pair26
 b pc_pair25
pc_dispatch_1_4:
 tbz x17,#0,pc_pair28
 b pc_pair27
pc_dispatch_2_0:
 tbz x17,#1,pc_dispatch_1_0
 tbz x17,#0,pc_pair30
 b pc_pair29
pc_dispatch_1_0:
 tbz x17,#0,pc_finish
 b pc_pair31
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
pc_pair8:
 ldp x9,x10,[x1],#16
 ldr x8,[x9]
 ldr x5,[x10]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit8
pc_cont8:
pc_pair9:
 ldp x11,x12,[x1],#16
 ldr x7,[x11]
 ldr x8,[x12]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit9
pc_cont9:
pc_pair10:
 ldp x13,x14,[x1],#16
 ldr x5,[x13]
 ldr x7,[x14]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit10
pc_cont10:
pc_pair11:
 ldp x15,x16,[x1],#16
 ldr x8,[x15]
 ldr x5,[x16]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit11
pc_cont11:
pc_pair12:
 ldp x9,x10,[x1],#16
 ldr x7,[x9]
 ldr x8,[x10]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit12
pc_cont12:
pc_pair13:
 ldp x11,x12,[x1],#16
 ldr x5,[x11]
 ldr x7,[x12]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit13
pc_cont13:
pc_pair14:
 ldp x13,x14,[x1],#16
 ldr x8,[x13]
 ldr x5,[x14]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit14
pc_cont14:
pc_pair15:
 ldp x15,x16,[x1],#16
 ldr x7,[x15]
 ldr x8,[x16]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit15
pc_cont15:
pc_pair16:
 ldp x9,x10,[x1],#16
 ldr x5,[x9]
 ldr x7,[x10]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit16
pc_cont16:
pc_pair17:
 ldp x11,x12,[x1],#16
 ldr x8,[x11]
 ldr x5,[x12]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit17
pc_cont17:
pc_pair18:
 ldp x13,x14,[x1],#16
 ldr x7,[x13]
 ldr x8,[x14]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit18
pc_cont18:
pc_pair19:
 ldp x15,x16,[x1],#16
 ldr x5,[x15]
 ldr x7,[x16]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit19
pc_cont19:
pc_pair20:
 ldp x9,x10,[x1],#16
 ldr x8,[x9]
 ldr x5,[x10]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit20
pc_cont20:
pc_pair21:
 ldp x11,x12,[x1],#16
 ldr x7,[x11]
 ldr x8,[x12]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit21
pc_cont21:
pc_pair22:
 ldp x13,x14,[x1],#16
 ldr x5,[x13]
 ldr x7,[x14]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit22
pc_cont22:
pc_pair23:
 ldp x15,x16,[x1],#16
 ldr x8,[x15]
 ldr x5,[x16]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit23
pc_cont23:
pc_pair24:
 ldp x9,x10,[x1],#16
 ldr x7,[x9]
 ldr x8,[x10]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit24
pc_cont24:
pc_pair25:
 ldp x11,x12,[x1],#16
 ldr x5,[x11]
 ldr x7,[x12]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit25
pc_cont25:
pc_pair26:
 ldp x13,x14,[x1],#16
 ldr x8,[x13]
 ldr x5,[x14]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit26
pc_cont26:
pc_pair27:
 ldp x15,x16,[x1],#16
 ldr x7,[x15]
 ldr x8,[x16]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit27
pc_cont27:
pc_pair28:
 ldp x9,x10,[x1],#16
 ldr x5,[x9]
 ldr x7,[x10]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit28
pc_cont28:
pc_pair29:
 ldp x11,x12,[x1],#16
 ldr x8,[x11]
 ldr x5,[x12]
 cmp x7,x8
 ccmp x8,x5,#4,ne
 b.eq pc_hit29
pc_cont29:
pc_pair30:
 ldp x13,x14,[x1],#16
 ldr x7,[x13]
 ldr x8,[x14]
 cmp x5,x7
 ccmp x7,x8,#4,ne
 b.eq pc_hit30
pc_cont30:
pc_pair31:
 ldp x15,x16,[x1],#16
 ldr x5,[x15]
 ldr x7,[x16]
 cmp x8,x5
 ccmp x5,x7,#4,ne
 b.eq pc_hit31
pc_cont31:
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
pc_hit8:
 cmp x7,x8
 b.ne pc_second8
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second8:
 cmp x8,x5
 b.ne pc_cont8
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont8
pc_hit9:
 cmp x5,x7
 b.ne pc_second9
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second9:
 cmp x7,x8
 b.ne pc_cont9
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont9
pc_hit10:
 cmp x8,x5
 b.ne pc_second10
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second10:
 cmp x5,x7
 b.ne pc_cont10
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont10
pc_hit11:
 cmp x7,x8
 b.ne pc_second11
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second11:
 cmp x8,x5
 b.ne pc_cont11
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont11
pc_hit12:
 cmp x5,x7
 b.ne pc_second12
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second12:
 cmp x7,x8
 b.ne pc_cont12
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont12
pc_hit13:
 cmp x8,x5
 b.ne pc_second13
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second13:
 cmp x5,x7
 b.ne pc_cont13
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont13
pc_hit14:
 cmp x7,x8
 b.ne pc_second14
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second14:
 cmp x8,x5
 b.ne pc_cont14
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont14
pc_hit15:
 cmp x5,x7
 b.ne pc_second15
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second15:
 cmp x7,x8
 b.ne pc_cont15
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont15
pc_hit16:
 cmp x8,x5
 b.ne pc_second16
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second16:
 cmp x5,x7
 b.ne pc_cont16
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont16
pc_hit17:
 cmp x7,x8
 b.ne pc_second17
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second17:
 cmp x8,x5
 b.ne pc_cont17
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont17
pc_hit18:
 cmp x5,x7
 b.ne pc_second18
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second18:
 cmp x7,x8
 b.ne pc_cont18
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont18
pc_hit19:
 cmp x8,x5
 b.ne pc_second19
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second19:
 cmp x5,x7
 b.ne pc_cont19
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont19
pc_hit20:
 cmp x7,x8
 b.ne pc_second20
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second20:
 cmp x8,x5
 b.ne pc_cont20
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont20
pc_hit21:
 cmp x5,x7
 b.ne pc_second21
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second21:
 cmp x7,x8
 b.ne pc_cont21
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont21
pc_hit22:
 cmp x8,x5
 b.ne pc_second22
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second22:
 cmp x5,x7
 b.ne pc_cont22
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont22
pc_hit23:
 cmp x7,x8
 b.ne pc_second23
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second23:
 cmp x8,x5
 b.ne pc_cont23
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont23
pc_hit24:
 cmp x5,x7
 b.ne pc_second24
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second24:
 cmp x7,x8
 b.ne pc_cont24
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont24
pc_hit25:
 cmp x8,x5
 b.ne pc_second25
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second25:
 cmp x5,x7
 b.ne pc_cont25
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont25
pc_hit26:
 cmp x7,x8
 b.ne pc_second26
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second26:
 cmp x8,x5
 b.ne pc_cont26
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont26
pc_hit27:
 cmp x5,x7
 b.ne pc_second27
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second27:
 cmp x7,x8
 b.ne pc_cont27
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont27
pc_hit28:
 cmp x8,x5
 b.ne pc_second28
 cmp x3,x4
 b.hs pc_done
 stp x16,x9,[x3],#16
pc_second28:
 cmp x5,x7
 b.ne pc_cont28
 cmp x3,x4
 b.hs pc_done
 stp x9,x10,[x3],#16
 b pc_cont28
pc_hit29:
 cmp x7,x8
 b.ne pc_second29
 cmp x3,x4
 b.hs pc_done
 stp x10,x11,[x3],#16
pc_second29:
 cmp x8,x5
 b.ne pc_cont29
 cmp x3,x4
 b.hs pc_done
 stp x11,x12,[x3],#16
 b pc_cont29
pc_hit30:
 cmp x5,x7
 b.ne pc_second30
 cmp x3,x4
 b.hs pc_done
 stp x12,x13,[x3],#16
pc_second30:
 cmp x7,x8
 b.ne pc_cont30
 cmp x3,x4
 b.hs pc_done
 stp x13,x14,[x3],#16
 b pc_cont30
pc_hit31:
 cmp x8,x5
 b.ne pc_second31
 cmp x3,x4
 b.hs pc_done
 stp x14,x15,[x3],#16
pc_second31:
 cmp x5,x7
 b.ne pc_cont31
 cmp x3,x4
 b.hs pc_done
 stp x15,x16,[x3],#16
 b pc_cont31

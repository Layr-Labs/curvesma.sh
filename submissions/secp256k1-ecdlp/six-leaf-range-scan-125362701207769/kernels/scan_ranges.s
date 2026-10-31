// Eight-comparison range-local scanner. x0 buckets, x1 pointer cursor,
// x2 end-cursor table, x3 output cursor, x4 output end. No cross-range pair.
.globl scan_ranges
scan_ranges:
    cbz x0,sr8_done
sr8_bucket:
    ldr x18,[x2],#8
    cmp x1,x18
    b.hs sr8_next
    ldr x16,[x1],#8
    ldr x5,[x16]
    sub x17,x18,x1
    lsr x17,x17,#3
    and x17,x17,#7
    cbz x17,sr8_check
sr8_prefix:
    ldr x9,[x1],#8
    ldr x7,[x9]
    cmp x5,x7
    b.eq sr8_prefix_hit
sr8_prefix_continue:
    mov x16,x9
    mov x5,x7
    subs x17,x17,#1
    b.ne sr8_prefix
sr8_check:
    cmp x1,x18
    b.hs sr8_next
sr8_loop:
    ldp x9,x10,[x1],#16
    ldr x7,[x9]
    cmp x5,x7
    b.eq sr8_hit1
sr8_continue1:
    ldr x5,[x10]
    cmp x7,x5
    b.eq sr8_hit2
sr8_continue2:
    ldp x11,x12,[x1],#16
    ldr x7,[x11]
    cmp x5,x7
    b.eq sr8_hit3
sr8_continue3:
    ldr x5,[x12]
    cmp x7,x5
    b.eq sr8_hit4
sr8_continue4:
    ldp x13,x14,[x1],#16
    ldr x7,[x13]
    cmp x5,x7
    b.eq sr8_hit5
sr8_continue5:
    ldr x5,[x14]
    cmp x7,x5
    b.eq sr8_hit6
sr8_continue6:
    ldp x15,x16,[x1],#16
    ldr x7,[x15]
    cmp x5,x7
    b.eq sr8_hit7
sr8_continue7:
    ldr x5,[x16]
    cmp x7,x5
    b.eq sr8_hit8
sr8_continue8:
    cmp x1,x18
    b.lo sr8_loop
sr8_next:
    subs x0,x0,#1
    b.ne sr8_bucket
sr8_done:
    ret
sr8_prefix_hit:
    cmp x3,x4
    b.hs sr8_done
    stp x16,x9,[x3],#16
    b sr8_prefix_continue
sr8_hit1:
    cmp x3,x4
    b.hs sr8_done
    stp x16,x9,[x3],#16
    b sr8_continue1
sr8_hit2:
    cmp x3,x4
    b.hs sr8_done
    stp x9,x10,[x3],#16
    b sr8_continue2
sr8_hit3:
    cmp x3,x4
    b.hs sr8_done
    stp x10,x11,[x3],#16
    b sr8_continue3
sr8_hit4:
    cmp x3,x4
    b.hs sr8_done
    stp x11,x12,[x3],#16
    b sr8_continue4
sr8_hit5:
    cmp x3,x4
    b.hs sr8_done
    stp x12,x13,[x3],#16
    b sr8_continue5
sr8_hit6:
    cmp x3,x4
    b.hs sr8_done
    stp x13,x14,[x3],#16
    b sr8_continue6
sr8_hit7:
    cmp x3,x4
    b.hs sr8_done
    stp x14,x15,[x3],#16
    b sr8_continue7
sr8_hit8:
    cmp x3,x4
    b.hs sr8_done
    stp x15,x16,[x3],#16
    b sr8_continue8

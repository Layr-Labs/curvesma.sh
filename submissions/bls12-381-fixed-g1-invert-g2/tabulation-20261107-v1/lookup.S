// Analytic scalar AArch64 leaf; never executed on the host.
// x0: sorted pointers to full records, x1: n<=127, x2: current key.
// Return x0: insertion index, x1: exact equality flag.
.text
.p2align 2
.global fapi_key_lookup
.type fapi_key_lookup, %function
fapi_key_lookup:
    mov x3, xzr
    mov x4, x1
.Ltop:
    cmp x3, x4
    b.hs .Ldone
    sub x5, x4, x3
    lsr x5, x5, #1
    add x5, x5, x3
    ldr x6, [x0, x5, lsl #3]
    add x6, x6, #640
    mov x7, x2
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    subs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    sbcs xzr, x8, x10
    sbcs xzr, x9, x11
    ldr x8, [x6]
    ldr x10, [x7]
    sbcs xzr, x8, x10
    b.lo .Lless
    mov x4, x5
    b .Ltop
.Lless:
    add x3, x5, #1
    b .Ltop
.Ldone:
    cmp x3, x1
    b.hs .Lmiss
    ldr x6, [x0, x3, lsl #3]
    add x6, x6, #640
    mov x7, x2
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldp x8, x9, [x6], #16
    ldp x10, x11, [x7], #16
    cmp x8, x10
    b.ne .Lmiss
    cmp x9, x11
    b.ne .Lmiss
    ldr x8, [x6]
    ldr x10, [x7]
    cmp x8, x10
    b.ne .Lmiss
    mov x0, x3
    mov x1, #1
    ret
.Lmiss:
    mov x0, x3
    mov x1, xzr
    ret
.size fapi_key_lookup, .-fapi_key_lookup

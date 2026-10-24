.text
.p2align 2
_batch_entry:
    sub sp, sp, #16
    stur x30, [sp, #8]
.p2align 2
_batch_prefix:
L_prefix:
    ldr x2, [x0, #96]
    add x0, x0, #512
    subs x1, x1, #1
    b.ne L_prefix
    ret
.p2align 2
_leaf_span:
L_loop:
    ldr x5, [x19, #96]
    ldp x0, x1, [x5, #0]
    ldp x2, x3, [x5, #16]
    ldp x4, x6, [x5, #32]
    ldp x7, x8, [x5, #48]
    ldp x9, x23, [x5, #64]
    ldp x24, x25, [x5, #80]
    ldp x5, x26, [x19, #0]
    ldp x27, x28, [x19, #16]
    ldp x29, x30, [x19, #32]
    ldp x0, x1, [x19, #48]
    ldp x2, x3, [x19, #64]
    ldp x4, x6, [x19, #80]
    // [C] charged 6M reference addition: operand block -> new coordinates, and its Montgomery prefix accumulator in x8 (not an instruction)
    stur x8, [sp, #0]
    stp x5, x26, [x19, #0]
    stp x27, x28, [x19, #16]
    stp x29, x30, [x19, #32]
    stp x0, x1, [x19, #48]
    stp x2, x3, [x19, #64]
    stp x4, x6, [x19, #80]
    ldp x0, x1, [x19, #0]
    ldp x2, x3, [x19, #16]
    ldp x4, x5, [x19, #32]
    ldp x6, x7, [x22, #0]
    umulh x8, x5, x6
    umulh x9, x4, x7
    add x8, x8, x9
    madd x8, x5, x7, x8
    ldp x6, x7, [x22, #16]
    umulh x9, x3, x6
    add x8, x8, x9
    madd x8, x4, x6, x8
    umulh x9, x2, x7
    add x8, x8, x9
    madd x8, x3, x7, x8
    ldp x6, x7, [x22, #32]
    umulh x9, x1, x6
    add x8, x8, x9
    madd x8, x2, x6, x8
    umulh x9, x0, x7
    add x8, x8, x9
    madd x8, x1, x7, x8
    ldr x6, [x22, #48]
    madd x9, x0, x6, x8
    sub x7, x9, #16
    cmn x7, #31
    b.hs L_fallback
    umulh x0, x16, x9
    mul x1, x17, x9
    umulh x2, x17, x9
    cmn x1, x0
    adc x1, x2, xzr
    add x2, x5, x1
    subs x3, x2, x17
    csel x4, x17, xzr, hs
    add x6, x3, #16
    cmp x6, #32
    b.ls L_fallback
    sub x2, x4, x3
    cmp x1, x5
    csel x3, x1, x5, lo
    csel x4, x5, x1, lo
    cset x6, lo
    cmp x2, x4
    csel x7, x2, x4, lo
    cmp x2, x3
    csel x7, x3, x7, lo
    csel x3, x2, x3, lo
    mov x8, #2
    csel x8, x8, x6, lo
    sub x7, x7, x3
    cmp x7, #16
    b.ls L_fallback
    ldp x23, x24, [x19, #48]
    ldp x25, x26, [x19, #64]
    ldp x27, x28, [x19, #80]
    subs x6, x10, x23
    sbcs x7, x11, x24
    sbcs x9, x12, x25
    sbcs x29, x13, x26
    sbcs x30, x14, x27
    sbc x5, x15, x28
    cmp x6, x23
    sbcs xzr, x7, x24
    sbcs xzr, x9, x25
    sbcs xzr, x29, x26
    sbcs xzr, x30, x27
    sbcs xzr, x5, x28
    cset x4, lo
    csel x0, x6, x23, lo
    csel x1, x7, x24, lo
    csel x2, x9, x25, lo
    csel x3, x29, x26, lo
    and x5, x0, #0xfffff
    cbz x1, L_report
    add x6, x4, x5, lsl #1
    add x6, x6, x6, lsl #1
    add x6, x6, x8
    add x6, x6, x6, lsl #2
    add x6, x21, x6, lsl #5
    str x6, [x19, #96]
    ldp x27, x28, [x19, #136]
    cmp x28, #20
    b.hs L_abandon
    ldp x7, x8, [x19, #104]
    ldp x9, x29, [x19, #120]
    cmp x0, x7
    ccmp x2, x8, #0, eq
    ccmp x3, x9, #0, eq
    ccmp x29, #1, #0, eq
    b.eq L_traversal
    tst x27, #63
    b.ne L_norefresh
    stp x0, x2, [x19, #104]
    mov x26, #1
    stp x3, x26, [x19, #120]
L_norefresh:
    adds x27, x27, #1
    adc x28, x28, xzr
    stp x27, x28, [x19, #136]
    add x19, x19, #512
    cmp x19, x20
    b.ne L_loop
    ldur x8, [sp, #0]
    cbz x8, L_zero
    ldur x30, [sp, #8]
    add sp, sp, #16
    ret
L_fallback:
    ldur x30, [sp, #8]
    add sp, sp, #16
    ret
L_report:
    ldur x30, [sp, #8]
    add sp, sp, #16
    ret
L_traversal:
    ldur x30, [sp, #8]
    add sp, sp, #16
    ret
L_abandon:
    ldur x30, [sp, #8]
    add sp, sp, #16
    ret
L_zero:
    ldur x30, [sp, #8]
    add sp, sp, #16
    ret

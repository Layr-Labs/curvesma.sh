// Analytical scalar AArch64 listing, not a native measurement.
// Preconditions: permitted DC ZVA, DCZID_EL0.BS=6, DZP=0; Normal writable
// mapped head memory [0,2^61), no concurrent accesses, flat 64-bit model.
// Each owner executes once after the run gate opens.
owner_init:
    mrs x2, dczid_el0
    and x2, x2, #31
    cmp x2, #6
    b.eq .Lsupported
    ret // Unsupported capability: no writes or attack.
.Lsupported:
    movz x0, #0
    movz x1, #2, lsl #32
.Lzero:
    .rept 1048576
    dc zva, x0
    add x0, x0, #256
    .endr
    subs x1, x1, #1
    b.eq .Ldone
    b .Lzero
.Ldone:
    dsb sy
    ret

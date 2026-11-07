// Analytic AArch64 hot-loop and scalar leaf evidence; never host-executed.
// Cold controller initializes/restores x19..x28 and a >=288-byte aligned frame.
// x19 current944-byte record,x20 alternate record,x21 chosen jump,
// x22:x23 little-indexed128-bit step count,x24:x25 next checkpoint count,
// x26 saved-record count,x27 jump-table base,x28 sorted pointer array.
// sp+0 cached label,sp+16 128-bit K,sp+32 pointer to table bases,
// sp+272..287 remaining full-key budget=floor(K/256), cold initialized.
// sp+48 saved caller LR,sp+64 cached32-byte BLAKE3 prefix CV,
// sp+96 mutable CV,sp+128 block64,sp+192..271 saved x19..x28.
// Cold controller branches here (not a returning BL) after saving LR/callee state.
// Final cold handlers restore LR and callee state before returning to caller.
.text
.p2align 2
.global fapi_hot_loop
.type fapi_hot_loop, %function
fapi_hot_loop:
.Ltick:
    ldp x9,x10,[sp,#16]
    subs xzr,x22,x9
    sbcs xzr,x23,x10
    b.hs fapi_cold_exhausted
    adds x22,x22,#1
    adc x23,x23,xzr
    ldr x9,[sp,#0]
    lsl x10,x9,#9
    add x10,x10,x9,lsl #7
    add x21,x27,x10
    mov x0,x20
    mov x1,x19
    mov x2,x21
    bl blst_fp12_mul
    mov x9,x19
    mov x19,x20
    mov x20,x9
    mov x0,x20
    mov x1,x19
    mov x2,x21
    bl fapi_scalar_update
    mov x0,x19
    add x1,x19,#928
    bl fapi_sign
    add x0,x19,#288
    add x1,x19,#640
    bl fapi_copy_key
    stp xzr,xzr,[sp,#128]
    stp xzr,xzr,[sp,#144]
    stp xzr,xzr,[sp,#160]
    stp xzr,xzr,[sp,#176]
    add x0,x19,#640
    ldr x1,[sp,#32]
    add x2,sp,#128
    bl fapi_tabulate
    ldr x9,[sp,#128]
    str x9,[x19,#936]
    ldp x0,x1,[sp,#64]
    stp x0,x1,[sp,#96]
    ldp x0,x1,[sp,#80]
    stp x0,x1,[sp,#112]
    add x0,sp,#96
    add x1,sp,#128
    bl fapi_blake3_final
    ldr w9,[sp,#120]
    ldr w10,[sp,#124]
    lsr x9,x9,#24
    rev w10,w10
    orr x9,x10,x9,lsl #32
    str x9,[sp,#0]
    mov x0,x28
    mov x1,x26
    mov x2,x19
    add x3,sp,#272
    bl fapi_prefix_lookup
    cmp x1,#2
    b.eq fapi_cold_exhausted
    cbnz x1,fapi_cold_found
    cmp x22,x24
    ccmp x23,x25,#0,eq
    b.eq fapi_cold_checkpoint
    b .Ltick
.size fapi_hot_loop,.-fapi_hot_loop
// Cold checkpoint handler preserves long-lived state, inserts record at x0,
// increments x26, adds s to x24:x25, then branches back to .Ltick.
// Insert by unsigned prefix order only; ties stay contiguous. The lower_bound
// miss index returned after a prefix run is also a valid insertion position.
// It may run only at t=s,...,1022s, never more than1022 times.
// Found/exhausted handlers terminate through cold verification/failure.
// All handler and entry/exit work is charged to fixed2^90.

.p2align 2
.global fapi_scalar_update
.type fapi_scalar_update, %function
fapi_scalar_update:
    mov x13,x0
    mov x14,x1
    mov x15,x2
    adrp x16,fapi_r
    add x16,x16,:lo12:fapi_r
    ldp x0,x1,[x13,#576]
    ldp x2,x3,[x13,#592]
    ldp x4,x5,[x15,#576]
    ldp x6,x7,[x15,#592]
    adds x0,x0,x4
    adcs x1,x1,x5
    adcs x2,x2,x6
    adc x3,x3,x7
    ldp x4,x5,[x16]
    ldp x6,x7,[x16,#16]
    subs x9,x0,x4
    sbcs x10,x1,x5
    sbcs x11,x2,x6
    sbcs x12,x3,x7
    csel x0,x9,x0,hs
    csel x1,x10,x1,hs
    csel x2,x11,x2,hs
    csel x3,x12,x3,hs
    stp x0,x1,[x14,#576]
    stp x2,x3,[x14,#592]
    ldp x0,x1,[x13,#608]
    ldp x2,x3,[x13,#624]
    ldp x4,x5,[x15,#608]
    ldp x6,x7,[x15,#624]
    adds x0,x0,x4
    adcs x1,x1,x5
    adcs x2,x2,x6
    adc x3,x3,x7
    ldp x4,x5,[x16]
    ldp x6,x7,[x16,#16]
    subs x9,x0,x4
    sbcs x10,x1,x5
    sbcs x11,x2,x6
    sbcs x12,x3,x7
    csel x0,x9,x0,hs
    csel x1,x10,x1,hs
    csel x2,x11,x2,hs
    csel x3,x12,x3,hs
    stp x0,x1,[x14,#608]
    stp x2,x3,[x14,#624]
    ret
.size fapi_scalar_update,.-fapi_scalar_update

.p2align 2
.global fapi_copy_key
.type fapi_copy_key, %function
fapi_copy_key:
    ldr x3,[x0,#40]
    str x3,[x1,#0]
    ldr x3,[x0,#32]
    str x3,[x1,#8]
    ldr x3,[x0,#24]
    str x3,[x1,#16]
    ldr x3,[x0,#16]
    str x3,[x1,#24]
    ldr x3,[x0,#8]
    str x3,[x1,#32]
    ldr x3,[x0,#0]
    str x3,[x1,#40]
    ldr x3,[x0,#88]
    str x3,[x1,#48]
    ldr x3,[x0,#80]
    str x3,[x1,#56]
    ldr x3,[x0,#72]
    str x3,[x1,#64]
    ldr x3,[x0,#64]
    str x3,[x1,#72]
    ldr x3,[x0,#56]
    str x3,[x1,#80]
    ldr x3,[x0,#48]
    str x3,[x1,#88]
    ldr x3,[x0,#136]
    str x3,[x1,#96]
    ldr x3,[x0,#128]
    str x3,[x1,#104]
    ldr x3,[x0,#120]
    str x3,[x1,#112]
    ldr x3,[x0,#112]
    str x3,[x1,#120]
    ldr x3,[x0,#104]
    str x3,[x1,#128]
    ldr x3,[x0,#96]
    str x3,[x1,#136]
    ldr x3,[x0,#184]
    str x3,[x1,#144]
    ldr x3,[x0,#176]
    str x3,[x1,#152]
    ldr x3,[x0,#168]
    str x3,[x1,#160]
    ldr x3,[x0,#160]
    str x3,[x1,#168]
    ldr x3,[x0,#152]
    str x3,[x1,#176]
    ldr x3,[x0,#144]
    str x3,[x1,#184]
    ldr x3,[x0,#232]
    str x3,[x1,#192]
    ldr x3,[x0,#224]
    str x3,[x1,#200]
    ldr x3,[x0,#216]
    str x3,[x1,#208]
    ldr x3,[x0,#208]
    str x3,[x1,#216]
    ldr x3,[x0,#200]
    str x3,[x1,#224]
    ldr x3,[x0,#192]
    str x3,[x1,#232]
    ldr x3,[x0,#280]
    str x3,[x1,#240]
    ldr x3,[x0,#272]
    str x3,[x1,#248]
    ldr x3,[x0,#264]
    str x3,[x1,#256]
    ldr x3,[x0,#256]
    str x3,[x1,#264]
    ldr x3,[x0,#248]
    str x3,[x1,#272]
    ldr x3,[x0,#240]
    str x3,[x1,#280]
    ret
.size fapi_copy_key,.-fapi_copy_key
.section .rodata
.p2align 3
fapi_r:
    .quad 0xffffffff00000001
    .quad 0x53bda402fffe5bfe
    .quad 0x3339d80809a1d805
    .quad 0x73eda753299d7d48

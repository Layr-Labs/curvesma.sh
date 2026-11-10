// Actual packed40 metadata/32-byte scalar I/O; ONLY isolated shared verifier executes.
.text
.p2align 2
.global prefix_copy4
.type prefix_copy4,%function
prefix_copy4:
    .rept 2
    ldp x2,x3,[x1],#16
    stp x2,x3,[x0],#16
    .endr
    ret
.size prefix_copy4,.-prefix_copy4

.global prefix_copy5
.type prefix_copy5,%function
prefix_copy5:
    .rept 2
    ldp x2,x3,[x1],#16
    stp x2,x3,[x0],#16
    .endr
    ldr x2,[x1],#8
    str x2,[x0],#8
    ret
.size prefix_copy5,.-prefix_copy5

.global prefix_copy320
.type prefix_copy320,%function
prefix_copy320:
    .rept 160
    ldp x2,x3,[x1],#16
    stp x2,x3,[x0],#16
    .endr
    ret
.size prefix_copy320,.-prefix_copy320

.global prefix_copy322
.type prefix_copy322,%function
prefix_copy322:
    .rept 161
    ldp x2,x3,[x1],#16
    stp x2,x3,[x0],#16
    .endr
    ret
.size prefix_copy322,.-prefix_copy322

.global prefix_copy512
.type prefix_copy512,%function
prefix_copy512:
    .rept 256
    ldp x2,x3,[x1],#16
    stp x2,x3,[x0],#16
    .endr
    ret
.size prefix_copy512,.-prefix_copy512

.global metadata_pending_get
.type metadata_pending_get,%function
metadata_pending_get:
    ldrb w0,[x0,#23]
    ret
.size metadata_pending_get,.-metadata_pending_get

.global metadata_pending_set
.type metadata_pending_set,%function
metadata_pending_set:
    strb w1,[x0,#23]
    ret
.size metadata_pending_set,.-metadata_pending_set

.global metadata_pending_clear
.type metadata_pending_clear,%function
metadata_pending_clear:
    strb wzr,[x0,#23]
    ret
.size metadata_pending_clear,.-metadata_pending_clear

.global metadata_ordinal
.type metadata_ordinal,%function
metadata_ordinal:
    ldp x0,x1,[x0,#24]
    ret
.size metadata_ordinal,.-metadata_ordinal

.global metadata_color
.type metadata_color,%function
metadata_color:
    ldrb w0,[x0,#22]
    lsr w0,w0,#7
    ret
.size metadata_color,.-metadata_color

.global metadata_color_pack
.type metadata_color_pack,%function
metadata_color_pack:
    ldrb w2,[x0,#22]
    orr w2,w2,w1,lsl #7
    strb w2,[x0,#22]
    ret
.size metadata_color_pack,.-metadata_color_pack

.global metadata_fingerprint_high
.type metadata_fingerprint_high,%function
metadata_fingerprint_high:
    ldr x0,[x0,#16]
    and x0,x0,#0x7fffffffffffff
    ret
.size metadata_fingerprint_high,.-metadata_fingerprint_high

.global baby_scalar_fill
.type baby_scalar_fill,%function
baby_scalar_fill:
    ldp x2,x3,[x1]
    stp x2,x3,[x0]
    stp xzr,xzr,[x0,#16]
    adds x2,x2,#1
    adc x3,x3,xzr
    stp x2,x3,[x1]
    ret
.size baby_scalar_fill,.-baby_scalar_fill

.global metadata_color_word_pack
.type metadata_color_word_pack,%function
metadata_color_word_pack:
    orr x5,x5,x1,lsl #55
    str x5,[x0,#16]
    ret
.size metadata_color_word_pack,.-metadata_color_word_pack

.global metadata_register_pack
.type metadata_register_pack,%function
metadata_register_pack:
    cbz x7,.Lregister_pack_invalid
    stp x1,x0,[x19]
    orr x5,x5,x20,lsl #55
    str x5,[x19,#16]
    stp x21,x22,[x19,#24]
.Lregister_pack_invalid:
    ret
.size metadata_register_pack,.-metadata_register_pack

.global nh_top8_fill
.type nh_top8_fill,%function
nh_top8_fill:
    str x0,[x6]
    mov x2,#0
    mov x3,#256
.Lnh_top8_fill_loop:
    str x2,[x0],#8
    add x2,x2,x1
    sub x4,x2,x5
    cmp x2,x5
    csel x2,x4,x2,hs
    subs x3,x3,#1
    b.ne .Lnh_top8_fill_loop
    ret
.size nh_top8_fill,.-nh_top8_fill

.global nh_affine32_fill
.type nh_affine32_fill,%function
nh_affine32_fill:
    str x0,[x6]
.Lnh_affine32_fill_loop:
    str x2,[x0],#8
    add x2,x2,x1
    sub x4,x2,x5
    cmp x2,x5
    csel x2,x4,x2,hs
    subs x3,x3,#1
    b.ne .Lnh_affine32_fill_loop
    ret
.size nh_affine32_fill,.-nh_affine32_fill

// Cold worker boot ONLY: x6 ACTUAL initialized full-p48 pointer,sets27/28/29.
.global worker_layout_init
.type worker_layout_init,%function
worker_layout_init:
    mov x27,#48
    mov x28,#96
    mov x29,x6
    mov x23,x8
    mov x24,x9
    ldp x25,xzr,[x9,#-16]
    mov x26,x10
    ret
.size worker_layout_init,.-worker_layout_init

// Actual production INLINE caller macros. The fixture wrapper below separately
// saves/restores LR and has THREE explicitly reported non-production overhead ops.
.macro WORKER_BABY_CALL
    ldr x0,[x26]
    add x1,x26,#16
    bl baby_scalar_fill
    add x0,x0,#32
    str x0,[x26]
.endm
.macro WORKER_FIELD_ARGS
    ldp x0,x1,[x26,#32]
    ldp x2,x3,[x26,#48]
    bl recurrence_fp2_step
.endm
.macro WORKER_ROTATE_TRACES
    ldp x0,x1,[x26,#32]
    ldr x2,[x26,#48]
    stp x2,x0,[x26,#32]
    str x1,[x26,#48]
.endm
.macro WORKER_KEY_ARGS
    ldr x1,[x26,#32]
    bl canonical_twelve_bfirst
.endm
.global fixture_worker_baby_call
.type fixture_worker_baby_call,%function
fixture_worker_baby_call:
    stp x29,x30,[sp,#-16]!
    WORKER_BABY_CALL
    ldp x29,x30,[sp],#16
    ret
.size fixture_worker_baby_call,.-fixture_worker_baby_call

// Inline production worker progression:five prefix ops plus actual branch=6.
// Fixture adds only a status MOV and RET(2),both executed/countable.
.macro WORKER_PROGRESS_PREFIX
    adds x21,x21,#1
    adc x22,x22,xzr
    ldp x4,x5,[x26,#96]
    subs xzr,x21,x4
    sbcs xzr,x22,x5
.endm
.macro WORKER_PROGRESS_TO target
    WORKER_PROGRESS_PREFIX
    b.lo \target
.endm
.global worker_progress_fixture
.type worker_progress_fixture,%function
worker_progress_fixture:
    WORKER_PROGRESS_TO .Lworker_progress_more
    mov x0,#0
    ret
.Lworker_progress_more:
    mov x0,#1
    ret
.size worker_progress_fixture,.-worker_progress_fixture

.macro WORKER_FIELD_STATUS target
    tbnz x0,#0,\target
.endm
.global worker_field_call_fixture
.type worker_field_call_fixture,%function
worker_field_call_fixture:
    stp x29,x30,[sp,#-16]!
    WORKER_FIELD_ARGS
    WORKER_FIELD_STATUS .Lworker_field_failed
    mov x0,#0
    ldp x29,x30,[sp],#16
    ret
.Lworker_field_failed:
    mov x0,#1
    ldp x29,x30,[sp],#16
    ret
.size worker_field_call_fixture,.-worker_field_call_fixture

.macro WORKER_KEY_STATUS target
    cbz x7,\target
.endm
.global worker_key_call_fixture
.type worker_key_call_fixture,%function
worker_key_call_fixture:
    stp x29,x30,[sp,#-16]!
    ldr x1,[x26,#32]
    bl canonical_twelve_bfirst
    WORKER_KEY_STATUS .Lworker_key_failed
    mov x17,x0
    mov x0,#0
    ldp x29,x30,[sp],#16
    ret
.Lworker_key_failed:
    mov x17,x0
    mov x0,#1
    ldp x29,x30,[sp],#16
    ret
.size worker_key_call_fixture,.-worker_key_call_fixture

.global worker_port_old_fixture
.type worker_port_old_fixture,%function
worker_port_old_fixture:
    movz x3,#20480,lsl #32
    movz x4,#12288
    madd x0,x0,x4,x3
    movz x4,#6144
    madd x0,x1,x4,x0
    add x1,x0,#1,lsl #12
    add x1,x1,#1984
    add x2,x1,#8
    ret
.size worker_port_old_fixture,.-worker_port_old_fixture

.global worker_port_wide_fixture
.type worker_port_wide_fixture,%function
worker_port_wide_fixture:
    movz x3,#23552,lsl #32
    movz x4,#9216
    madd x0,x0,x4,x3
    movz x4,#4608
    madd x0,x1,x4,x0
    add x1,x0,#1,lsl #12
    add x1,x1,#448
    add x2,x1,#8
    ret
.size worker_port_wide_fixture,.-worker_port_wide_fixture

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
    stp x1,x0,[x19]
    orr x5,x5,x20,lsl #55
    str x5,[x19,#16]
    stp x21,x22,[x19,#24]
    ret
.size metadata_register_pack,.-metadata_register_pack

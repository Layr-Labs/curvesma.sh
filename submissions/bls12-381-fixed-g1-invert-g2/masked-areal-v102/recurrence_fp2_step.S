// Inert analytic source;shared credential-free llm-verifier ONLY.
// Full fused leaf:out192,prev1/prev2 canonical192,tablebase,productsbase=const3728+528,p48 compatibility.
// Frame528 unchanged:out/table48/56,SP64dead,products24COLD/entry96in28,27preserved48;workspace80..135,temp144..527.
// Internal worker ABI: entry28=96 from PAID cold worker_layout_init.
// BOTH uses x29 tempbase;map phase preserves x0 out;no BL/x18,allcallee/SP restored;LRx30 untouched.
.text
.p2align 2
.macro MSTEP off,width,advance,ra,rb
ubfx x3,x19,#\off,#\width
madd x1,x3,x28,x26
ldp x3,x21,[x1,#(\ra+0)]
adds x4,x4,x3
adcs x5,x5,x21
ldp x3,x21,[x1,#(\ra+16)]
adcs x6,x6,x3
adcs x7,x7,x21
ldp x3,x21,[x1,#(\ra+32)]
adcs x8,x8,x3
adcs x9,x9,x21
adc x10,x10,xzr
ldp x3,x21,[x1,#(\rb+0)]
adds x11,x11,x3
adcs x12,x12,x21
ldp x3,x21,[x1,#(\rb+16)]
adcs x13,x13,x3
adcs x14,x14,x21
ldp x3,x21,[x1,#(\rb+32)]
adcs x15,x15,x3
adcs x16,x16,x21
adc x17,x17,xzr
.if \advance
add x26,x26,x28,lsl #32
.endif
.endm
.macro COLUMNBASE0
.endm
.macro COLUMNBASE1
add x26,x26,x28,lsl #30
.endm
.macro COLUMNBASE2
add x26,x26,x2,lsl #1
.endm
.macro COLUMNBASE3
add x26,x26,x28,lsl #30
.endm
.macro COLUMN ci,ti
COLUMNBASE\ti
ldr x19,[x29,#(48*\ci+0)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+8)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+16)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+24)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+32)]
MSTEP 0,32,1,0,48
MSTEP 32,32,1,0,48
ldr x19,[x29,#(48*\ci+40)]
MSTEP 0,32,1,0,48
MSTEP 32,30,0,0,48
.endm
.macro WSTEP advance,width,keep,tag,init
madd x1,x3,x28,x26
ldp x3,x21,[x1,#(0+0)]
adds x4,x4,x3
adcs x5,x5,x21
ldp x3,x21,[x1,#(0+16)]
adcs x6,x6,x3
adcs x7,x7,x21
ldp x3,x21,[x1,#(0+32)]
adcs x8,x8,x3
adcs x9,x9,x21
.if \keep
.if \init
adc x10,xzr,xzr
.endif
.if 1-\init
adc x10,x10,xzr
.endif
.endif
.if 1-\keep
.Lomitcarry_field_a_\tag:
.endif
ldp x3,x21,[x1,#(48+0)]
adds x11,x11,x3
adcs x12,x12,x21
ldp x3,x21,[x1,#(48+16)]
adcs x13,x13,x3
adcs x14,x14,x21
ldp x3,x21,[x1,#(48+32)]
adcs x15,x15,x3
adcs x16,x16,x21
.if \keep
.if \init
adc x17,xzr,xzr
.endif
.if 1-\init
adc x17,x17,xzr
.endif
.endif
.if 1-\keep
.Lomitcarry_field_b_\tag:
.endif
.if \advance
add x26,x26,x28,lsl #\width
.endif
.endm
.macro COLUMNWIDE ci,second
    .if 1-\second
    .if \ci
    ldr x26,[sp,#56]
    .endif
    .endif
    .if \second
    add x26,x26,x28,lsl #35
    .endif
    ldp x19,x2,[x29,#(48*\ci+0)]
    ubfx x3,x19,#0,#35
    WSTEP 1,35,\second,\ci_0,0
    extr x3,x2,x19,#35
    and x3,x3,#0x7ffffffff
    WSTEP 1,35,\second,\ci_1,0
    ubfx x3,x2,#6,#35
    WSTEP 1,35,\second,\ci_2,0
    ldr x19,[x29,#(48*\ci+16)]
    extr x3,x19,x2,#41
    and x3,x3,#0x7ffffffff
    WSTEP 1,35,\second,\ci_3,0
    ubfx x3,x19,#12,#35
    WSTEP 1,35,\second,\ci_4,0
    ldr x2,[x29,#(48*\ci+24)]
    extr x3,x2,x19,#47
    and x3,x3,#0x7ffffffff
    WSTEP 1,35,\second,\ci_5,0
    ubfx x3,x2,#18,#35
    WSTEP 1,35,\second,\ci_6,0
    ldr x19,[x29,#(48*\ci+32)]
    extr x3,x19,x2,#53
    and x3,x3,#0x3ffffffff
    WSTEP 1,34,\second,\ci_7,0
    ubfx x3,x19,#23,#34
    WSTEP 1,34,1,\ci_8,1-\second
    ldr x2,[x29,#(48*\ci+40)]
    extr x3,x2,x19,#57
    and x3,x3,#0x3ffffffff
    WSTEP 1,34,1,\ci_9,0
    ubfx x3,x2,#27,#35
    WSTEP 0,35,1,\ci_10,0
    .if \second
    add x26,x26,x28,lsl #35
    .if \ci-1
    movz x2,#0x438,lsl #32
    .endif
    .endif
.endm
.macro QREDUCE
extr x1,x10,x9,#5
udiv x1,x1,x25
add x26,x24,x1,lsl #6
ldp x11,x12,[x26,#0]
ldp x13,x14,[x26,#16]
ldp x15,x16,[x26,#32]
ldr x17,[x26,#48]
subs x11,x4,x11
sbcs x12,x5,x12
sbcs x13,x6,x13
sbcs x14,x7,x14
sbcs x15,x8,x15
sbcs x16,x9,x16
sbcs x17,x10,x17
sub x19,x26,#64
csel x26,x19,x26,lo
ldp x11,x12,[x26,#0]
ldp x13,x14,[x26,#16]
ldp x15,x16,[x26,#32]
subs x4,x4,x11
sbcs x5,x5,x12
sbcs x6,x6,x13
sbcs x7,x7,x14
sbcs x8,x8,x15
sbcs x9,x9,x16
.endm
.macro QREDUCE_MASKED
.Lmasked_reduce_start:
extr x1,x10,x9,#5
udiv x1,x1,x25
add x26,x24,x1,lsl #6
ldp x11,x12,[x26,#0]
ldp x13,x14,[x26,#16]
ldp x15,x16,[x26,#32]
ldr x17,[x26,#48]
subs x11,x4,x11
sbcs x12,x5,x12
sbcs x13,x6,x13
sbcs x14,x7,x14
sbcs x15,x8,x15
sbcs x16,x9,x16
sbcs x17,x10,x17
b.hs .Lmasked_result
.Lmasked_correction:
sub x26,x26,#64
ldp x11,x12,[x26,#0]
ldp x13,x14,[x26,#16]
ldp x15,x16,[x26,#32]
subs x11,x4,x11
sbcs x12,x5,x12
sbcs x13,x6,x13
sbcs x14,x7,x14
sbcs x15,x8,x15
sbcs x16,x9,x16
ldr x19,[sp,#16]
ldp x4,x5,[x19,#128]
adds x4,x4,#1
adc x5,x5,xzr
stp x4,x5,[x19,#128]
ldp x6,x7,[x19,#144]
subs x8,x4,x6
sbcs x9,x5,x7
b.hs .Lmasked_abort
.Lmasked_result:
.endm
.macro ROWPAIR sa,sb,c0,c1,c2,c3,oa,ob,p2,p3,masked
ldp x4,x5,[x0,#(\sa+0)]
ldp x6,x7,[x0,#(\sa+16)]
ldp x8,x9,[x0,#(\sa+32)]
ldp x11,x12,[x0,#(\sb+0)]
ldp x13,x14,[x0,#(\sb+16)]
ldp x15,x16,[x0,#(\sb+32)]
COLUMNWIDE \c0,0
COLUMNWIDE \c1,1
COLUMN \c2,\p2
COLUMN \c3,\p3
stp x11,x12,[sp,#80]
stp x13,x14,[sp,#96]
stp x15,x16,[sp,#112]
str x17,[sp,#128]
.if \masked
QREDUCE_MASKED
stp x11,x12,[x0,#(\oa+0)]
stp x13,x14,[x0,#(\oa+16)]
stp x15,x16,[x0,#(\oa+32)]
.endif
.if 1-\masked
QREDUCE
stp x4,x5,[x0,#(\oa+0)]
stp x6,x7,[x0,#(\oa+16)]
stp x8,x9,[x0,#(\oa+32)]
.endif
ldp x4,x5,[sp,#80]
ldp x6,x7,[sp,#96]
ldp x8,x9,[sp,#112]
ldr x10,[sp,#128]
QREDUCE
stp x4,x5,[x0,#(\ob+0)]
stp x6,x7,[x0,#(\ob+16)]
stp x8,x9,[x0,#(\ob+32)]
.endm
.macro BOTH so,di
ldp x4,x5,[x1,#0]
ldp x6,x7,[x1,#16]
ldp x8,x9,[x1,#32]
ldp x10,x11,[x2,#0]
ldp x12,x13,[x2,#16]
ldp x14,x15,[x2,#32]
adds x0,x4,x10
adcs x3,x5,x11
adcs x17,x6,x12
adcs x16,x7,x13
adcs x26,x8,x14
adc x21,x9,x15
subs x4,x4,x10
sbcs x5,x5,x11
sbcs x6,x6,x12
sbcs x7,x7,x13
sbcs x8,x8,x14
sbcs x9,x9,x15
stp x4,x5,[x29,#(\di+0)]
stp x6,x7,[x29,#(\di+16)]
stp x8,x9,[x29,#(\di+32)]
stp x0,x3,[x29,#(\so+0)]
stp x17,x16,[x29,#(\so+16)]
stp x26,x21,[x29,#(\so+32)]
.endm
.global recurrence_fp2_step
.type recurrence_fp2_step,%function
recurrence_fp2_step:
    sub sp,sp,#528
    stp x19,x21,[sp]
    stp x26,x27,[sp,#16]
    stp x28,x29,[sp,#32]
    stp x0,x3,[sp,#48]
    add x29,sp,#144
    BOTH 288,0
    add x1,x1,#48
    add x2,x2,#48
    BOTH 336,48
    add x1,x1,#48
    add x2,x2,#48
    BOTH 96,192
    add x1,x1,#48
    add x2,x2,#48
    BOTH 144,240
    ldp x0,x26,[sp,#48]
    ROWPAIR 0,48,0,1,2,3,0,48,0,1,1
    ROWPAIR 96,144,4,5,6,7,96,144,2,3,0
.Lfield_epilogue:
    ldp x19,x21,[sp]
    ldp x26,x27,[sp,#16]
    ldp x28,x29,[sp,#32]
    add sp,sp,#528
    ret
// Lowbit of already-returned aligned output pointer signals failure;caller MUST test.
.Lmasked_abort:
    orr x0,x0,#1
    b .Lfield_epilogue
.size recurrence_fp2_step,.-recurrence_fp2_step

// Inert shared-verifier-only MACRO definitions;no standalone helper entry/function.
// Expanded ONLYinsidepublicrecurrence_fp2_step. SP iswrapper528frame:
// outputpointerSP96,constant3728pointerSP128;tempbase29/tablebase30 suppliedinline.
// Signed22p/twoFULLA/Bvectors/ALL381bits/liveaccumulators;ceilD estimate and QSTORE29includingstores.
// PublicwrapperownsALL19..29/LR/SP/x18preservation;noextra frame or BL.
.text
.p2align 2
.macro MSTEP ra,rb,negative,width,advance
madd x1,x3,x28,x26
ldp x20,x21,[x1,#(\ra+0)]
ldp x22,x23,[x1,#(\ra+16)]
ldp x24,x25,[x1,#(\ra+32)]
.if \negative
subs x4,x4,x20
sbcs x5,x5,x21
sbcs x6,x6,x22
sbcs x7,x7,x23
sbcs x8,x8,x24
sbcs x9,x9,x25
sbc x10,x10,xzr
.endif
.if (1-\negative)
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,xzr
.endif
ldp x20,x21,[x1,#(\rb+0)]
ldp x22,x23,[x1,#(\rb+16)]
ldp x24,x25,[x1,#(\rb+32)]
adds x11,x11,x20
adcs x12,x12,x21
adcs x13,x13,x22
adcs x14,x14,x23
adcs x15,x15,x24
adcs x16,x16,x25
adc x17,x17,xzr
.if \advance
.if (\width-38)
add x26,x26,x27
.endif
.if (39-\width)
add x26,x26,x27,lsr #1
.endif
.endif
.endm
.macro COLUMN ci,tableid,imag
mov x3,#\tableid
madd x26,x3,x2,x30
ldr x19,[x29,#(48*\ci+0)]
ubfx x3,x19,#0,#39
MSTEP (48*\imag),(48*(1-\imag)),\imag,39,1
ldr x20,[x29,#(48*\ci+8)]
extr x3,x20,x19,#39
and x3,x3,#0x3fffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ubfx x3,x19,#13,#38
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ldr x20,[x29,#(48*\ci+16)]
extr x3,x20,x19,#51
and x3,x3,#0x3fffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ubfx x3,x19,#25,#38
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ldr x20,[x29,#(48*\ci+24)]
extr x3,x20,x19,#63
and x3,x3,#0x3fffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ldr x20,[x29,#(48*\ci+32)]
extr x3,x20,x19,#37
and x3,x3,#0x3fffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ubfx x3,x19,#11,#38
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ldr x20,[x29,#(48*\ci+40)]
extr x3,x20,x19,#49
and x3,x3,#0x3fffffffff
mov x19,x20
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,1
ubfx x3,x19,#23,#38
MSTEP (48*\imag),(48*(1-\imag)),\imag,38,0
.endm
.macro QSTORE a0,a1,a2,a3,a4,a5,a6,off
extr x1,\a6,\a5,#5
udiv x1,x1,x2
add x26,x28,#528
add x26,x26,x1,lsl #6
ldp x19,x20,[x26,#0]
ldp x21,x22,[x26,#16]
ldp x23,x24,[x26,#32]
ldr x25,[x26,#48]
subs \a0,\a0,x19
sbcs \a1,\a1,x20
sbcs \a2,\a2,x21
sbcs \a3,\a3,x22
sbcs \a4,\a4,x23
sbcs \a5,\a5,x24
sbcs \a6,\a6,x25
ldp x19,x20,[x28,#448]
ldp x21,x22,[x28,#464]
ldp x23,x24,[x28,#480]
subs x19,\a0,x19
sbcs x20,\a1,x20
sbcs x21,\a2,x21
sbcs x22,\a3,x22
sbcs x23,\a4,x23
sbcs x24,\a5,x24
b.lo .Lq_original\off
stp x19,x20,[x0,#(\off+0)]
stp x21,x22,[x0,#(\off+16)]
stp x23,x24,[x0,#(\off+32)]
b .Lq_done\off
.Lq_original\off:
stp \a0,\a1,[x0,#(\off+0)]
stp \a2,\a3,[x0,#(\off+16)]
stp \a4,\a5,[x0,#(\off+32)]
.Lq_done\off:
.endm
.macro ROWPAIR sa,sb,c0,c1,c2,c3,oa,ob
ldr x1,[sp,#96]
ldp x4,x5,[x1,#(\sa+0)]
ldp x6,x7,[x1,#(\sa+16)]
ldp x8,x9,[x1,#(\sa+32)]
mov x10,xzr
ldp x11,x12,[x1,#(\sb+0)]
ldp x13,x14,[x1,#(\sb+16)]
ldp x15,x16,[x1,#(\sb+32)]
mov x17,xzr
ldr x26,[sp,#128]
add x26,x26,#1936
ldp x20,x21,[x26]
ldp x22,x23,[x26,#16]
ldp x24,x25,[x26,#32]
ldr x3,[x26,#48]
adds x4,x4,x20
adcs x5,x5,x21
adcs x6,x6,x22
adcs x7,x7,x23
adcs x8,x8,x24
adcs x9,x9,x25
adc x10,x10,x3
COLUMN \c0,0,0
COLUMN \c1,0,1
COLUMN \c2,1,0
COLUMN \c3,1,1
ldr x28,[sp,#128]
ldr x2,[x28,#512]
QSTORE x4,x5,x6,x7,x8,x9,x10,\oa
QSTORE x11,x12,x13,x14,x15,x16,x17,\ob
.endm

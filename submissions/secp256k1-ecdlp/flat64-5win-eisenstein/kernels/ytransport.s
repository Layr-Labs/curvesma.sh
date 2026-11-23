// Card45-H. Node-local state YS=34TiB+0x7100000 is private to the Y leaf proc.
// yconsume: x0=frame base (YVIEW or owner's Normal YQ), x1=capped count,
// x2=isolated Normal completion word, x3=epoch, SP=private ordinary stack.
// ychunk uses YS[48]=X metadata base and function pointers at YS[64..88].
// YS[96..127]=subgroup order; boot initializes YS[56]=0 (coin failure).
// YS[192..207]=128-bit lifetime fixup meter (lo,hi);
// YS[208..223]=128-bit threshold DELTA_F-FIX_BATCH (lo,hi), initialized once.
// Neither yconsume nor the outer round/group schedule resets this meter.
.text
.globl yclean
yclean:
cbz x1,2f
lsl x1,x1,#3
1:
dc cvac,x0
add x0,x0,#64
subs x1,x1,#1
b.ne 1b
2:
dsb sy
ret

// Historical yowner/yoffer bodies, inactive in the current ystage schedule.
// Current production uses yowner_round/yoffer_scalar and calls yconsume directly.
// Legacy owner Y processor, after generation barrier. x0 points to eight Normal
// words: source,count,first_id_lo,first_id_hi,DONE,epoch,DRAINED,offer_port.
// DRAINED has a separate 2KiB monitor granule. The owner main waits for DRAINED
// before the barrier which permits buffer reuse. No source read occurs after DRAINED.
.globl yowner
yowner:
stp x0,x30,[sp,#-16]!
ldp x0,x1,[x0]
bl yclean
ldr x4,[sp]
ldp x0,x1,[x4]
ldp x2,x3,[x4,#32]
bl yconsume
ldr x4,[sp]
ldr x0,[x4,#56]
ldp x1,x2,[x4,#16]
ldr x5,[x4,#40]
ldp x3,x4,[x4]
mov x6,#1023
bl yoffer
ldr x4,[sp]
ldp x5,x6,[x4,#40]
stlr x5,[x6]
dmb ish
ldp x0,x30,[sp],#16
ret

// x0=48-byte offer port, x1/x2=first recipient id, x3=source frame base,
// x4=capped number of Ys, x5=epoch, x6=recipient count (1023).
// The old completing-store assumption below is not a current fabric service.
.globl yoffer
yoffer:
stp x3,x4,[x0,#16]
str x5,[x0,#32]
1:
stp x1,x2,[x0]
stlur x5,[x0,#40]
dmb ish
add x1,x1,#1
subs x6,x6,#1
b.ne 1b
ret

.globl yconsume
yconsume:
sub sp,sp,#64
stp x0,x1,[sp]
stp x2,x3,[sp,#16]
str x30,[sp,#32]
cbz x1,2f
1:
ldp x1,x0,[sp]
mov x4,#4194304
cmp x0,x4
csel x0,x0,x4,lo
str x0,[sp,#40]
bl ychunk
ldp x1,x2,[sp]
ldr x0,[sp,#40]
add x1,x1,x0,lsl #9
sub x2,x2,x0
stp x1,x2,[sp]
cbnz x2,1b
2:
ldp x2,x3,[sp,#16]
stlr x3,[x2]
dmb ish
ldr x30,[sp,#32]
add sp,sp,#64
ret

.globl ycoin
ycoin:
mov x17,#2
1:
mrs x24,rndr
mrs x25,rndr
mrs x26,rndr
mrs x27,rndr
orr x8,x24,x25
orr x9,x26,x27
orr x8,x8,x9
cbz x8,2f
ldp x4,x5,[x16,#96]
ldp x6,x7,[x16,#112]
subs xzr,x24,x4
sbcs xzr,x25,x5
sbcs xzr,x26,x6
sbcs xzr,x27,x7
b.lo 3f
2:
subs x17,x17,#1
b.ne 1b
mov x0,#1
ret
3:
mov x0,#0
ret

// Full glue around the unchanged priced leaf/inversion/keys primitives.
// All non-SP registers may be clobbered by any primitive. No state lives in
// those registers across calls. Bwd deliberately borrows SP for read-only Ys.
.globl ychunk
ychunk:
mov x16,#0x7100000
movk x16,#0x2200,lsl #32
mov x9,sp
stp x9,x30,[x16]
ldp x8,x9,[x16,#192]
ldp x10,x11,[x16,#208]
subs xzr,x10,x8
sbcs xzr,x11,x9
b.lo 9f
stp x0,x1,[x16,#16]
add x0,x0,x0,lsl #1
lsl x0,x0,#1
mov x3,#0x2a2000000000
sub x3,x3,x0,lsl #9
stp x0,x3,[x16,#32]
bl ycoin
cbnz x0,8f
ldp x0,x3,[x16,#32]
ldr x1,[x16,#24]
add x1,x1,#64
ldr x16,[x16,#64]
blr x16
mov x16,#0x7100000
movk x16,#0x2200,lsl #32
ldp x8,x9,[x16,#192]
lsr x10,x0,#32
adds x8,x8,x10
adc x9,x9,xzr
stp x8,x9,[x16,#192]
stp x24,x25,[x16,#128]
stp x26,x27,[x16,#144]
add x0,x16,#160
add x1,x16,#128
ldr x16,[x16,#72]
blr x16
mov x16,#0x7100000
movk x16,#0x2200,lsl #32
ldp x0,x2,[x16,#160]
ldp x14,x15,[x16,#176]
ldp x4,x1,[x16,#16]
sub x4,x4,#1
add x4,x1,x4,lsl #9
mov sp,x4
mov x3,#0x2a2000000000
ldr x16,[x16,#80]
blr x16
mov x16,#0x7100000
movk x16,#0x2200,lsl #32
ldp x8,x9,[x16,#192]
lsr x10,x5,#32
adds x8,x8,x10
adc x9,x9,xzr
stp x8,x9,[x16,#192]
ldr x4,[x16]
mov sp,x4
ldr x2,[x16,#48]
add x2,x2,#6
ldr x3,[x16,#40]
add x3,x3,#64
ldr x16,[x16,#88]
blr x16
mov x16,#0x7100000
movk x16,#0x2200,lsl #32
9:
ldp x4,x30,[x16]
mov sp,x4
ret
8:
str x0,[x16,#56]
ldp x4,x30,[x16]
mov sp,x4
ret

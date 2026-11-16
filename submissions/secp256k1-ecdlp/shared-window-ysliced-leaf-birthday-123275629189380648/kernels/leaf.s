// A macro used below and not defined here is defined in kernels/leafm.s, which ships compressed in code.b94 (proof.md 3 and 7); the run assembles leafm.s followed by this file.
.macro LJ0C k,w0l,w0h,w1l,w1h
add x5,\w0h,\w1h
mvn x5,x5
cbhs \w0h,\w1h,kh\k
cbhi \w0h,x5,kv\k
and x0,x1,\w0l,lsr #11
ldadd x4,x5,[x0]
cbne x5,x0,ks\@
str x0,[x0]
b kn\@
ks\@:
bfxil \w1h,\w0h,#21,#21
bfxil \w0h,\w1l,#22,#42
stp \w0h,\w1h,[x5]
stp x2,x26,[x5,#16]
kn\@:
.endm
.macro SLH k,w0l,w0h,w1l,w1h
kh\k:
add x0,\w0l,\w1l
cbhi \w1h,x5,kt\k
mvn \w0l,x0
and x0,x1,\w1l,lsr #11
ldadd x4,\w0h,[x0]
cbne \w0h,x0,ks\@
str x0,[x0]
b km\k
ks\@:
bfxil x5,\w1h,#21,#21
bfxil \w1h,\w0l,#22,#42
stp \w1h,x5,[\w0h]
stp x2,x26,[\w0h,#16]
km\k:
.endm
.macro SLV k,w0l,w0h,w1l,w1h
kv\k:
add x0,\w0l,\w1l
kt\k:
bic x0,x1,x0,lsr #11
ldadd x4,\w1h,[x0]
cbne \w1h,x0,ks\@
str x0,[x0]
b ku\k
ks\@:
bfxil \w0h,x5,#21,#21
bfxil x5,\w0l,#22,#42
stp x5,\w0h,[\w1h]
stp x2,x26,[\w1h,#16]
ku\k:
.endm
.macro NETK k,a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h
add x24,\a1l,\a0l,lsl #1
add x25,\a1h,\a0h,lsl #1
sub \a0l,\a0l,\b1l
sub \a0h,\a0h,\b1h
sub \a1l,\a1l,\g0l
sub \a1h,\a1h,\g0h
sub x27,\g0l,\a0l
sub x28,\g0h,\a0h
sub x29,\b0l,\a1l
sub x30,\b0h,\a1h
sub \a0l,\a0l,x29
sub \a0h,\a0h,x30
sub x0,\a0l,x0
sub x5,\a0h,x5
add \b1l,\b1l,x0
add \b1h,\b1h,x5
sub x24,x0,x24
sub x25,x5,x25
LJ0C \k,x24,x25,x29,x30
KT2 \k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h
SLH \k,x24,x25,x29,x30
KT2 \k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h
SLV \k,x24,x25,x29,x30
KT2 \k,\a0l,\a0h,\a1l,\a1h,\b0l,\b0h,\b1l,\b1h,\g0l,\g0h
SLH d\k,\a0l,\a0h,\a1l,\a1h
KT3 \k,\b0l,\b0h,\b1l,\b1h
SLV d\k,\a0l,\a0h,\a1l,\a1h
KT3 \k,\b0l,\b0h,\b1l,\b1h
SLH e\k,\b1l,\b1h,x27,x28
b re\k
SLV e\k,\b1l,\b1h,x27,x28
re\k:
.endm
.macro KT3 k,b0l,b0h,b1l,b1h
add \b1l,\b1l,\b0l,lsl #1
add \b1h,\b1h,\b0h,lsl #1
LJ0C e\k,\b1l,\b1h,x27,x28
b re\k
.endm
.macro KT2 k,a0l,a0h,a1l,a1h,b0l,b0h,b1l,b1h,g0l,g0h
add \a0l,\a0l,\g0l,lsl #1
add \a0h,\a0h,\g0h,lsl #1
sub \a1l,\a1l,\b1l
sub \a1h,\a1h,\b1h
LJ0C d\k,\a0l,\a0h,\a1l,\a1h
KT3 \k,\b0l,\b0h,\b1l,\b1h
.endm
.macro KFC
ldp x22,x23,[x3,#192]
ldp x24,x25,[x3,#208]
add x14,x14,x22
add x15,x15,x23
add x16,x16,x24
add x17,x17,x25
add x18,x18,x22
add x19,x19,x23
add x20,x20,x24
add x21,x21,x25
.endm
.macro LK1 tb
mov x2,#\tb
ldp x28,x29,[x1,#32]
ldr x5,[x1,#48]
add x9,x2,x28
ld64b x16,[x9]
add x17,x2,x29
ld64b x8,[x17]
adds x20,x20,x12
adcs x21,x21,x13
adcs x22,x22,x14
adcs x23,x23,x15
cset x28,hs
add x17,x2,x5
ld64b x8,[x17]
LKP_ADD5 x20,x21,x22,x23,x28,x12,x13,x14,x15
add x17,x2,x30
ld64b x8,[x17]
LKP_ADD5 x20,x21,x22,x23,x28,x12,x13,x14,x15
add x17,x2,x6
ld64b x8,[x17]
LKP_ADD5 x20,x21,x22,x23,x28,x12,x13,x14,x15
ld64b x8,[x1]
LKP_ADD5 x20,x21,x22,x23,x28,x8,x9,x10,x11
LKFOLD4 x20,x21,x22,x23,x28,x8
.endm
.macro FWD6_ONE
ldp x6,x7,[x1,#288]
str x7,[x3,#288]
LKP6_PAIR 0x6c000000000000,LD_OB3_2
st64b x20,[x3]
P1X x20,x24,x25,x26,x27,x21,x22,x23,x22,x23,x20,x21
add x2,x3,#1024
st64b x16,[x2]
P1X x16,x20,x21,x22,x23,x17,x18,x19,x26,x27,x24,x25
LK1 0x76000000000000
add x2,x3,#2048
st64b x20,[x2]
P1X x20,x24,x25,x26,x27,x21,x22,x23,x26,x27,x24,x25
add x1,x1,#512
add x3,x3,#3072
.endm
.macro KPH kp,km
ld64b x6,[x3]
add x0,x3,#64
ld64b x14,[x0]
KFC
ldp x22,x23,[x3,#-32]
ldp x0,x5,[x3,#-16]
NETK \kp,x14,x15,x16,x17,x6,x7,x8,x9,x22,x23
ldp x22,x23,[x3,#128]
ldp x0,x5,[x3,#144]
NETK \km,x18,x19,x20,x21,x10,x11,x12,x13,x22,x23
add x2,x2,#1
add x3,x3,#512
.endm
.macro KEYPAIRC kp0,km0,kp1,km1
KPH \kp0,\km0
KPH \kp1,\km1
.endm
.text
.p2align 2
.globl fwd
fwd:
mov x4,#977
movk x4,#1,lsl #32
stp x29,x30,[sp,#-16]!
mov x5,#0x22000000000000
str xzr,[x5]
and x5,x0,#6
add x0,x0,x5
add x0,x0,x5,lsl #1
cmp x5,#4
b.eq 6f
b.hs 7f
cbnz x5,5f
cbz x0,4f
1:
FWD6_ONE
5:
FWD6_ONE
6:
FWD6_ONE
7:
FWD6_ONE
subs x0,x0,#24
b.ne 1b
4:
mov x0,#0x22000000000000
ldr x0,[x0]
ldp x29,x30,[sp],#16
ret
.globl bwd
bwd:
mov x4,#977
movk x4,#1,lsl #32
mov x5,sp
stp x5,x30,[x1,#8]
str xzr,[x1]
sub x3,x3,#448
tbz x5,#9,3f
1:
EFGK6
BWDPK01
BWDPK00
BWDPK10
sub sp,sp,#512
3:
EFGK6
BWDPK01
BWDPK00
BWDPK10
mov x5,sp
sub sp,sp,#512
tst x5,#0x7fffffff
b.ne 1b
mov x7,#0x22000000000000
ldp x6,x30,[x7,#8]
mov sp,x6
ldr x5,[x7]
ret
.globl keys
keys:
str x30,[sp,#-16]!
mov x1,#0x1fffffffffe000
mov x4,#-32
1:
ldr x26,[x3,#224]
sub x2,x2,#6
KEYPAIRC p00,m00,p10,m10
KEYPAIRC p20,m20,p30,m30
KEYPAIRC p40,m40,p50,m50
tbz x3,#37,1b
ldr x30,[sp],#16
ret

// proof.md 3 (W6): fwd (PFWD2), bwd (BWDPA, BWDPB), keys; x4 = C0; fixups by branch (3.4)
.macro RF o,a,b,c,d,e,f,g,h,i,t,u
mul \e,\o,\a
umulh \t,\o,\a
mul \u,\o,\b
adds \f,\t,\u
umulh \t,\o,\b
mul \u,\o,\c
adcs \g,\t,\u
umulh \t,\o,\c
mul \u,\o,\d
adcs \h,\t,\u
umulh \t,\o,\d
adc \i,\t,xzr
.endm
.macro RA4 o,b0,b1,b2,b3,a0,a1,a2,a3,a4,a5,y,r,s,t,v
mul \r,\o,\b0
adds \a0,\a0,\r
umulh \s,\o,\b0
mul \r,\o,\b1
adcs \s,\s,\r
umulh \t,\o,\b1
mul \r,\o,\b2
adcs \t,\t,\r
umulh \v,\o,\b2
mul \r,\o,\b3
adcs \v,\v,\r
umulh \r,\o,\b3
adc \r,\r,xzr
adds \a1,\a1,\s
adcs \a2,\a2,\t
adcs \a3,\a3,\v
adcs \a4,\a4,\r
adc \a5,\y,xzr
.endm
.macro SHX l,a,d,p,q,k
mul \l,\d,x4
umulh \d,\d,x4
adds \a,\a,\d
b.lo 9f
adcs \p,\p,xzr
adcs \q,\q,xzr
stadd x4,[\k]
9:
.endm
.macro FDX e,f,i,j,o,p,l,h,y,z,k
mul \l,\i,x4
umulh \h,\i,x4
madd \h,\j,x4,\h
adds \o,\e,\l
adcs \p,\f,\h
b.lo 9f
adcs \y,\y,xzr
adc \z,\z,xzr
stadd x4,[\k]
9:
.endm
.macro FS d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
subs \d0,\a0,\b0
sbcs \d1,\a1,\b1
sbcs \d2,\a2,\b2
sbcs \d3,\a3,\b3
csel \t,x4,xzr,lo
subs \d0,\d0,\t
b.hs 9f
sbcs \d1,\d1,xzr
sbcs \d2,\d2,xzr
sbc \d3,\d3,xzr
stadd x4,[x3]
9:
.endm
.macro FADDW d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
adds \d0,\a0,\b0
adcs \d1,\a1,\b1
adcs \d2,\a2,\b2
adcs \d3,\a3,\b3
csel \t,x4,xzr,hs
adds \d0,\d0,\t
b.lo 9f
adcs \d1,\d1,xzr
adcs \d2,\d2,xzr
adc \d3,\d3,xzr
stadd x4,[x3]
9:
.endm
.macro SQ13 a0,a1,a2,a3,w0,c1,c2,c3,c4,c5,c6,c7,t
mul \c3,\a0,\a3
mul \t,\a1,\a2
adds \c3,\c3,\t
umulh \c4,\a0,\a3
umulh \t,\a1,\a2
adcs \c4,\c4,\t
umulh \c5,\a1,\a3
adc \c5,\c5,xzr
umulh \c2,\a0,\a1
mul \t,\a0,\a2
adds \c2,\c2,\t
umulh \t,\a0,\a2
adcs \c3,\c3,\t
mul \t,\a1,\a3
adcs \c4,\c4,\t
mul \t,\a2,\a3
adcs \c5,\c5,\t
umulh \c6,\a2,\a3
adc \c6,\c6,xzr
mul \c1,\a0,\a1
lsr \c7,\c6,#63
extr \c6,\c6,\c5,#63
extr \c5,\c5,\c4,#63
extr \c4,\c4,\c3,#63
extr \c3,\c3,\c2,#63
extr \c2,\c2,\c1,#63
mul \w0,\a0,\a0
umulh \t,\a0,\a0
adds \c1,\t,\c1,lsl #1
mul \t,\a1,\a1
adcs \c2,\c2,\t
umulh \t,\a1,\a1
adcs \c3,\c3,\t
mul \t,\a2,\a2
adcs \c4,\c4,\t
umulh \t,\a2,\a2
adcs \c5,\c5,\t
mul \t,\a3,\a3
adcs \c6,\c6,\t
umulh \t,\a3,\a3
adc \c7,\c7,\t
mul \t,\c4,x4
adds \w0,\w0,\t
umulh \c4,\c4,x4
mul \t,\c5,x4
adcs \c4,\c4,\t
umulh \c5,\c5,x4
mul \t,\c6,x4
adcs \c5,\c5,\t
umulh \c6,\c6,x4
mul \t,\c7,x4
adcs \c6,\c6,\t
umulh \c7,\c7,x4
adc \c7,\c7,xzr
adds \c1,\c1,\c4
adcs \c2,\c2,\c5
adcs \c3,\c3,\c6
adc \c7,\c7,xzr
.endm
.macro F5 w0,w1,w2,w3,w4,o0,o1,ta,tb
mul \ta,\w4,x4
umulh \tb,\w4,x4
adds \o0,\w0,\ta
adcs \o1,\w1,\tb
b.lo 9f
adcs \w2,\w2,xzr
adc \w3,\w3,xzr
stadd x4,[x3]
9:
.endm
.macro RA o,a,b,c,d,r,s,t,v,u,e,f,g,h,i,j,y,m,n
mul \r,\o,\a
umulh \s,\o,\a
mul \u,\o,\b
adds \s,\s,\u
umulh \t,\o,\b
mul \u,\o,\c
adcs \t,\t,\u
umulh \v,\o,\c
mul \u,\o,\d
adcs \v,\v,\u
umulh \u,\o,\d
adc \u,\u,xzr
adds \e,\e,\r
adcs \f,\f,\s
adcs \m,\g,\t
adcs \n,\h,\v
adcs \i,\i,\u
adc \j,\y,xzr
.endm
.macro P1X
RF x20,x24,x25,x26,x27,x18,x14,x15,x16,x17,x28,x19
SHX x5,x24,x27,x25,x26,sp
RA x21,x5,x24,x25,x26,x2,x30,x12,x28,x19,x18,x14,x15,x16,x17,x29,xzr,x15,x16
SHX x27,x5,x26,x24,x25,sp
RA x22,x27,x5,x24,x25,x2,x30,x12,x28,x19,x18,x14,x15,x16,x17,x29,x29,x15,x16
SHX x26,x27,x25,x5,x24,sp
RA x23,x26,x27,x5,x24,x2,x30,x12,x20,x21,x18,x14,x15,x16,x17,x29,x29,x12,x13
FDX x18,x14,x17,x29,x10,x11,x5,x2,x12,x13,sp
.endm
.macro P1XB
RF x6,x10,x11,x12,x13,x14,x15,x16,x17,x18,x28,x19
SHX x5,x10,x13,x11,x12,sp
RA x7,x5,x10,x11,x12,x2,x20,x21,x28,x19,x14,x15,x16,x17,x18,x29,xzr,x16,x17
SHX x13,x5,x12,x10,x11,sp
RA x8,x13,x5,x10,x11,x2,x20,x21,x28,x19,x14,x15,x16,x17,x18,x29,x29,x16,x17
SHX x12,x13,x11,x5,x10,sp
RA x9,x12,x13,x5,x10,x2,x20,x21,x22,x23,x14,x15,x16,x17,x18,x29,x29,x26,x27
FDX x14,x15,x18,x29,x24,x25,x5,x2,x26,x27,sp
.endm
.macro AC4F a0,a1,a2,a3,b0,b1,b2,b3,c
adds \a0,\a0,\b0
adcs \a1,\a1,\b1
adcs \a2,\a2,\b2
adcs \a3,\a3,\b3
cset \c,hs
.endm
.macro AC4 a0,a1,a2,a3,b0,b1,b2,b3,c
adds \a0,\a0,\b0
adcs \a1,\a1,\b1
adcs \a2,\a2,\b2
adcs \a3,\a3,\b3
adc \c,\c,xzr
.endm
.macro LF a0,a1,a2,a3,c,t,k
mul \t,\c,x4
adds \a0,\a0,\t
b.lo 9f
adcs \a1,\a1,xzr
adcs \a2,\a2,xzr
adcs \a3,\a3,xzr
csel \t,x4,xzr,hs
add \a0,\a0,\t
stadd x4,[\k]
9:
.endm
.macro W1F a
add x29,x2,\a
ld64b x10,[x29]
AC4 x20,x21,x22,x23,x10,x11,x12,x13,x5
AC4 x6,x7,x8,x9,x14,x15,x16,x17,x19
.endm
.macro W2F
W1F x18
W1F x28
.endm
.macro EFWD2
ldp x2,x18,[x3,#-408]
ldp x20,x21,[x1,#48]
ldp x22,x23,[x1,#64]
ldp x18,x28,[x1,#336]
add x29,x2,x18
ld64b x10,[x29]
adds x6,x20,x14
adcs x7,x21,x15
adcs x8,x22,x16
adcs x9,x23,x17
cset x19,hs
AC4F x20,x21,x22,x23,x10,x11,x12,x13,x5
W1F x28
ldp x18,x28,[x1,#352]
W2F
ldp x18,x28,[x1,#368]
W2F
LF x20,x21,x22,x23,x5,x29,sp
LF x6,x7,x8,x9,x19,x29,sp
.endm
.macro PFWD2
EFWD2
st64b x20,[x3]
add x3,x3,#1024
P1X
st64b x6,[x3]
add x3,x3,#1024
P1XB
.endm
.macro BLK5 s0,s1,s2,s3,o0,o1,o2,o3
orr x20,x5,\s0,lsl #8
ld64b x8,[x20]
extr x20,\s1,\s0,#52
orr x20,x5,x20,lsl #8
ld64b x12,[x20]
adds \o0,x8,x16
adcs \o1,x9,x17
adcs \o2,x10,x18
adcs \o3,x11,x19
cset x1,hs
extr x20,\s2,\s1,#40
orr x20,x6,x20,lsl #8
ld64b x12,[x20]
adds \o0,\o0,x12
adcs \o1,\o1,x13
adcs \o2,\o2,x14
adcs \o3,\o3,x15
adc x1,x1,xzr
extr x20,\s3,\s2,#28
orr x20,x6,x20,lsl #8
ld64b x12,[x20]
adds \o0,\o0,x16
adcs \o1,\o1,x17
adcs \o2,\o2,x18
adcs \o3,\o3,x19
adc x1,x1,xzr
lsr x20,\s3,#16
sub x20,x5,x20,lsl #6
ld64b x12,[x20]
adds \o0,\o0,x12
adcs \o1,\o1,x13
adcs \o2,\o2,x14
adcs \o3,\o3,x15
adc x1,x1,xzr
mul x20,x1,x4
adds \o0,\o0,x20
b.lo 9f
adcs \o1,\o1,xzr
adcs \o2,\o2,xzr
adcs \o3,\o3,xzr
csel x20,x4,xzr,hs
add \o0,\o0,x20
stadd x4,[x3]
9:
.endm
.macro PAIRB
mov x5,#0xf000000000000000
orr x6,x5,#64
BLK5 x21,x22,x23,x24,x27,x28,x29,x30
BLK5 x25,x26,x0,x2,x8,x9,x10,x11
.endm
.macro TS0 m,f,r3,r2,t
umulh \f,\m,\r2
madd \f,\m,\r3,\f
.endm
.macro TSK m,f,r3,r2,t
umulh \t,\m,\r2
add \f,\f,\t
madd \f,\m,\r3,\f
.endm
.macro W2S
RF x20,x25,x26,x27,x28,x9,x10,x11,x12,x13,x5,x6
RF x16,x25,x26,x27,x28,x0,x1,x2,x24,x29,x5,x6
SHX x5,x25,x28,x26,x27,x3
RA4 x21,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
RA4 x17,x5,x25,x26,x27,x0,x1,x2,x24,x29,x30,xzr,x6,x7,x8,x28
SHX x6,x5,x27,x25,x26,x3
RA4 x22,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
RA4 x18,x6,x5,x25,x26,x0,x1,x2,x24,x29,x30,x30,x7,x8,x28,x27
SHX x7,x6,x26,x5,x25,x3
RA4 x23,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
RA4 x19,x7,x6,x5,x25,x0,x1,x2,x24,x29,x30,x30,x8,x28,x27,x26
FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12,x3
FDX x0,x1,x29,x30,x0,x1,x6,x29,x2,x24,x3
.endm
.macro DF o0,o1,o2,o3,c,t
mul \t,\c,x4
umulh \c,\c,x4
adds \o0,\o0,\t
adcs \o1,\o1,\c
b.lo 9f
adcs \o2,\o2,xzr
adcs \o3,\o3,xzr
csel \t,x4,xzr,hs
adds \o0,\o0,\t
adc \o1,\o1,xzr
stadd x4,[x3]
9:
.endm
.macro EDELTA d,e0,e1,e2,e3,o0,o1,o2,o3,b0,b1,b2,b3,t,c
mul \t,\d,\b0
adds \o0,\e0,\t
mul \t,\d,\b1
adcs \o1,\e1,\t
mul \t,\d,\b2
adcs \o2,\e2,\t
mul \t,\d,\b3
adcs \o3,\e3,\t
umulh \c,\d,\b3
adc \c,\c,xzr
umulh \t,\d,\b0
adds \o1,\o1,\t
umulh \t,\d,\b1
adcs \o2,\o2,\t
umulh \t,\d,\b2
adcs \o3,\o3,\t
adc \c,\c,xzr
DF \o0,\o1,\o2,\o3,\c,\t
.endm
.macro DT
mul x14,x13,x9
umulh x15,x13,x9
mul x16,x13,x10
umulh x17,x13,x10
mul x18,x13,x11
umulh x19,x13,x11
mul x20,x13,x12
umulh x21,x13,x12
adds x0,x0,x14
adcs x1,x1,x15
adcs x2,x2,x17
adcs x5,x5,x19
adc x6,x21,xzr
adds x1,x1,x16
adcs x2,x2,x18
adcs x5,x5,x20
adc x6,x6,xzr
DF x0,x1,x2,x5,x6,x7
.endm
.macro C1BS
ldp x13,x14,[x3,#-496]
ldp x15,x16,[x3,#-480]
ldp x17,x18,[sp,#64]
ldp x19,x20,[sp,#80]
FS x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x21
RF x9,x13,x14,x15,x16,x0,x1,x2,x5,x6,x7,x8
SHX x17,x13,x16,x14,x15,x3
RA4 x10,x17,x13,x14,x15,x0,x1,x2,x5,x6,x18,xzr,x7,x8,x19,x20
SHX x16,x17,x15,x13,x14,x3
RA4 x11,x16,x17,x13,x14,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
SHX x15,x16,x14,x17,x13,x3
RA4 x12,x15,x16,x17,x13,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
FDX x0,x1,x6,x18,x0,x1,x7,x8,x2,x5,x3
stp x0,x1,[sp,#-128]
stp x2,x5,[sp,#-112]
ldr x13,[sp,#120]
DT
.endm
.macro P1W b
ld64b x6,[\b]
AC4 x0,x1,x2,x5,x6,x7,x8,x9,x14
AC4 x18,x19,x20,x21,x10,x11,x12,x13,x15
.endm
.macro P1 o
SQ13 x0,x1,x2,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14
F5 x6,x7,x8,x9,x13,x6,x7,x14,x13
add x30,sp,#\o
ld64b x22,[x30]
bfi x22,x6,#9,#43
extr x30,x7,x6,#43
bfi x23,x30,#9,#43
lsr x30,x7,#22
bfi x24,x30,#9,#42
bfi x25,x8,#9,#43
extr x30,x9,x8,#43
bfi x26,x30,#9,#43
lsr x30,x9,#22
bfi x27,x30,#9,#42
ld64b x14,[x22]
adds x0,x0,x14
adcs x1,x1,x15
adcs x2,x2,x16
adcs x5,x5,x17
cset x14,hs
ld64b x6,[x23]
AC4 x0,x1,x2,x5,x6,x7,x8,x9,x14
AC4F x18,x19,x20,x21,x10,x11,x12,x13,x15
P1W x24
P1W x25
P1W x26
P1W x27
LF x0,x1,x2,x5,x14,x6,x3
LF x18,x19,x20,x21,x15,x6,x3
ld64b x6,[x28]
RF x18,x6,x7,x8,x9,x14,x15,x16,x17,x18,x24,x25
RA4 x19,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,xzr,x24,x25,x26,x27
ld64b x6,[x29]
RA4 x20,x6,x7,x8,x9,x14,x15,x16,x17,x18,x19,x19,x24,x25,x26,x27
RA4 x21,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x19,x24,x25,x26,x27
FDX x14,x15,x18,x19,x14,x15,x24,x25,x16,x17,x3
FS x21,x22,x23,x24,x0,x1,x2,x5,x14,x15,x16,x17,x6
FADDW x25,x26,x0,x2,x0,x1,x2,x5,x14,x15,x16,x17,x6
.endm
.macro TRH op,a,b,c,d,r3,r2,t
\op \a,x17,\r3,\r2,\t
\op \b,x13,\r3,\r2,\t
\op \c,x19,\r3,\r2,\t
\op \d,x15,\r3,\r2,\t
.endm
.macro TRL op,a,b,c,d,r3,r2,t
\op \a,x16,\r3,\r2,\t
\op \b,x12,\r3,\r2,\t
\op \c,x18,\r3,\r2,\t
\op \d,x14,\r3,\r2,\t
.endm
.macro TFP4S
umulh x20,x6,x4
mul x1,x6,x4
umulh x12,x16,x4
mul x13,x16,x4
adds x1,x1,x14,lsr #32
adcs x20,x13,x20
adcs x5,x5,x12
adc x7,x14,xzr
TRH TS0,x21,x25,x27,x8,x16,x6,x0
TRH TSK,x22,x26,x28,x9,x6,x14,x16
TRL TS0,x21,x25,x27,x8,x6,x14,x0
TRH TSK,x23,x0,x29,x10,x7,x5,x6
TRL TSK,x22,x26,x28,x9,x7,x5,x6
TRH TSK,x24,x2,x30,x11,x5,x20,x6
TRL TSK,x23,x0,x29,x10,x5,x20,x6
TRL TSK,x24,x2,x30,x11,x20,x1,x6
.endm
.macro CORS r,c,d
ldp x7,x20,[\r,#\c]
ldp x5,x6,[\r,#\d]
add x17,x17,x20
add x16,x16,x7
add x13,x13,x20
add x12,x12,x7
add x19,x19,x6
add x18,x18,x5
add x15,x15,x6
add x14,x14,x5
.endm
.macro W1P a
add \a,x24,\a
ld64b x16,[\a]
AC4 x6,x7,x8,x9,x16,x17,x18,x19,x14
AC4 x10,x11,x12,x13,x20,x21,x22,x23,x15
.endm
.macro W2P
W1P x0
W1P x2
.endm
.macro LK2 oa,ob,oc
ldp x0,x2,[x3,#\oa]
add x0,x24,x0
ld64b x6,[x0]
add x2,x24,x2
ld64b x16,[x2]
AC4F x6,x7,x8,x9,x16,x17,x18,x19,x14
AC4F x10,x11,x12,x13,x20,x21,x22,x23,x15
ldp x0,x2,[x3,#\ob]
W2P
ldp x0,x2,[x3,#\oc]
W2P
.endm
.macro LK2S d
LF x6,x7,x8,x9,x14,x0,x3
LF x10,x11,x12,x13,x15,x0,x3
sub x29,sp,#\d
st64b x6,[x29]
.endm
.macro TC a,b
ldp x5,x14,[sp,#\a]
ldp x6,x16,[sp,#\b]
TFP4S
.endm
.macro INCO
P1 256
PAIRB
TC -160,-144
CORS x3,-400,-256
st64b x12,[sp]
sub sp,sp,#64
TC -480,-464
CORS sp,256,272
st64b x12,[sp]
add sp,sp,#192
TC -480,-464
st64b x12,[sp]
.endm
.macro INCE
P1 128
PAIRB
TC 0,16
st64b x12,[sp]
sub sp,sp,#192
TC 0,16
CORS sp,256,272
st64b x12,[sp]
add sp,sp,#64
TC 320,336
CORS x3,-400,-256
st64b x12,[sp]
.endm
.macro BWDPA
ldr x24,[sp,#96]
LK2 -240,-224,-208
ldp x16,x17,[x3,#-432]
ldp x18,x19,[x3,#-416]
AC4 x6,x7,x8,x9,x16,x17,x18,x19,x14
AC4 x10,x11,x12,x13,x16,x17,x18,x19,x15
LF x6,x7,x8,x9,x14,x0,x3
LF x10,x11,x12,x13,x15,x0,x3
ldp x1,x2,[x3,#-48]
ldp x5,x18,[x3,#-32]
sub x29,sp,#1216
ldr x0,[x29,#-96]
EDELTA x0,x10,x11,x12,x13,x14,x15,x16,x17,x1,x2,x5,x18,x19,x20
st64b x10,[x29]
ldr x0,[sp,#-96]
EDELTA x0,x6,x7,x8,x9,x10,x11,x12,x13,x1,x2,x5,x18,x19,x20
sub x29,sp,#192
st64b x6,[x29]
BWDT
.endm
.macro BWDPB
ldr x24,[sp,#96]
BWDT
.endm
.macro BWDT
LK2 -176,-160,-144
LK2S 576
LK2 -112,-96,-80
LK2S 384
ld64b x16,[sp]
W2S
stp x0,x1,[x3,#24]
stp x2,x24,[x3,#40]
C1BS
INCO
sub sp,sp,#512
ldp x0,x1,[sp,#256]
ldp x2,x5,[sp,#272]
INCE
ld64b x22,[x3]
sub sp,sp,#512
.endm
.macro FLUSH s,r
sub \s,x16,#4080
sub \s,\s,#4080
str \s,[x15]
orr \r,x15,#0x80000000000
mov x16,#8192
str x16,[\s,#-16]
ldadd x16,x16,[\r]
ldr \r,[sp]
cmp x16,\r
b.hs 6f
stp x15,x16,[\s,#-32]
6:
.endm
.macro LJ0 k,w0l,w0h,w1l,w1h
add x15,\w0h,\w1h
cmp \w0h,\w1h
b.hs kh\k
cmn \w0h,x15
b.hs kv\k
and x15,x1,\w0l,lsr #25
ldadd x4,x16,[x15]
add \w1l,x18,\w1l,lsr #28
bfxil \w0h,\w1h,#22,#22
stp \w0h,\w1h,[x16]
stp x17,\w1l,[x16,#16]
tbz x16,#13,kn\k
FLUSH \w0h,\w1h
kn\k:
.endm
.macro JHB k,w0l,w0h,w1l,w1h
kh\k:
cmn \w1h,x15
b.hs kt\k
add \w0l,\w0l,\w1l
mvn \w0h,x15
and x15,x1,\w1l,lsr #25
ldadd x4,x16,[x15]
eon \w0l,x18,\w0l,lsr #28
bfxil \w1h,\w0h,#22,#22
stp \w1h,\w0h,[x16]
stp x17,\w0l,[x16,#16]
tbz x16,#13,km\k
FLUSH \w0h,\w1h
km\k:
.endm
.macro JTB k,w0l,w0h,w1l,w1h
kv\k:
kt\k:
add \w1l,\w0l,\w1l
mvn \w1h,x15
bic x15,x1,\w1l,lsr #25
ldadd x4,x16,[x15]
add \w0l,x18,\w0l,lsr #28
bfxil \w1h,\w0h,#22,#22
stp \w1h,\w0h,[x16]
stp x17,\w0l,[x16,#16]
tbz x16,#13,ku\k
FLUSH \w0h,\w1h
ku\k:
.endm
.macro TLDA
lsl x18,x3,#30
ld64b x6,[x3]
sub x0,x3,#64
ld64b x20,[x0]
ldp x28,x29,[x3,#160]
ldp x30,x19,[x3,#176]
.endm
.macro TLDB
ldp x28,x29,[x3,#128]
ldp x30,x19,[x3,#144]
sub x3,x3,#512
.endm
.macro S1 a,b,c,d,e,f,g,h,i,j,l,m,n,o,q,r
add \n,\a,\e
add \o,\b,\f
add \n,\n,\i
add \o,\o,\j
add \q,\c,\g
add \r,\d,\h
add \q,\q,\l
add \r,\r,\m
.endm
.macro S2 a,b,c,d,e,f,g,h,i,j,l,m,n,o,q,r
sub \n,\e,\i
sub \o,\f,\j
sub \q,\g,\l
sub \r,\h,\m
sub \a,\a,\e
sub \b,\b,\f
sub \c,\c,\g
sub \d,\d,\h
add \e,\n,\q
add \f,\o,\r
add \g,\a,\e
add \h,\b,\f
sub \i,\c,\n
sub \j,\d,\o
.endm
.macro S3 a,b,c,d,e,f,q,r
sub \a,\a,\q
sub \b,\b,\r
add \c,\c,\e
add \d,\d,\f
.endm
.macro TSUM k,a,b,c,d,e,f,g,h,i,j,l,m,n,o,q,r
S1 \a,\b,\c,\d,\e,\f,\g,\h,\i,\j,\l,\m,\n,\o,\q,\r
LJ0 \k,\n,\o,\q,\r
S2 \a,\b,\c,\d,\e,\f,\g,\h,\i,\j,\l,\m,\n,\o,\q,\r
LJ0 d\k,\g,\h,\i,\j
S3 \a,\b,\c,\d,\e,\f,\q,\r
LJ0 e\k,\a,\b,\c,\d
.endm
.macro C3 k,y,a,b,c,d,e,f,q,r
S3 \a,\b,\c,\d,\e,\f,\q,\r
LJ0 \k,\a,\b,\c,\d
b \y
JHB \k,\a,\b,\c,\d
b \y
JTB \k,\a,\b,\c,\d
b \y
.endm
.macro C2 k,y,a,b,c,d,e,f,g,h,i,j,l,m,n,o,q,r
S2 \a,\b,\c,\d,\e,\f,\g,\h,\i,\j,\l,\m,\n,\o,\q,\r
LJ0 \k,\g,\h,\i,\j
C3 a\k,\y,\a,\b,\c,\d,\e,\f,\q,\r
JHB \k,\g,\h,\i,\j
C3 b\k,\y,\a,\b,\c,\d,\e,\f,\q,\r
JTB \k,\g,\h,\i,\j
C3 c\k,\y,\a,\b,\c,\d,\e,\f,\q,\r
.endm
.macro TJX k,a,b,c,d,e,f,g,h,i,j,l,m,n,o,q,r
JHB e\k,\a,\b,\c,\d
b kne\k
JTB e\k,\a,\b,\c,\d
b kne\k
JHB d\k,\g,\h,\i,\j
C3 g\k,kne\k,\a,\b,\c,\d,\e,\f,\q,\r
JTB d\k,\g,\h,\i,\j
C3 i\k,kne\k,\a,\b,\c,\d,\e,\f,\q,\r
JHB \k,\n,\o,\q,\r
C2 l\k,kne\k,\a,\b,\c,\d,\e,\f,\g,\h,\i,\j,\l,\m,\n,\o,\q,\r
JTB \k,\n,\o,\q,\r
C2 o\k,kne\k,\a,\b,\c,\d,\e,\f,\g,\h,\i,\j,\l,\m,\n,\o,\q,\r
.endm
.macro RP z,k
\z \k,x10,x11,x12,x13,x24,x25,x26,x27,x28,x29,x30,x19,x0,x2,x5,x14
.endm
.macro RM z,k
\z \k,x6,x7,x8,x9,x20,x21,x22,x23,x28,x29,x30,x19,x0,x2,x5,x14
.endm
.text
.p2align 2
.globl fwd
fwd:
mov x4,#977
movk x4,#1,lsl #32
stp xzr,x30,[sp,#-16]!
1:
PFWD2
PFWD2
PFWD2
PFWD2
PFWD2
PFWD2
PFWD2
PFWD2
tst x3,#0x7fc00
b.ne 2f
add x1,x1,#1024
2:
subs x0,x0,#32
b.ne 1b
ldp x0,x30,[sp],#16
ret
.globl bwd
bwd:
mov x4,#977
movk x4,#1,lsl #32
str x30,[x3,#8]
1:
BWDPA
BWDPB
BWDPA
BWDPB
BWDPA
BWDPB
BWDPA
BWDPB
add x5,sp,#384
tst x5,#0x7fe00
b.ne 1b
sub x5,x5,#512
tbnz x5,#34,2f
ldp x5,x6,[x3]
sub x3,x3,#1024
stp x5,x6,[x3]
b 1b
2:
ldp x0,x30,[x3]
ret
.globl keys
keys:
mov x1,#0x7ffffffff8
mov x4,#32
stp x12,x30,[sp,#-16]!
1:
TLDA
RP TSUM,p00
TLDB
RM TSUM,m00
TLDA
RP TSUM,p10
TLDB
RM TSUM,m10
tbz x3,#34,1b
ldp xzr,x30,[sp],#16
ret
RP TJX,p00
RM TJX,m00
RP TJX,p10
RM TJX,m10

.macro SD t
sub \t,sp,#128
stadd x4,[\t]
.endm
.macro LQ
ld64b x12,[x10]
.endm
.macro LP o
AD5 x16,x17,x18,x19,x0,x8,x9,x10,x11,\o
.endm
.macro AZ c,d,t
adcs \c,\c,xzr
adc \d,\d,xzr
SD \t
.endm
.macro KA a
add x16,x3,#\a
.endm
.macro KK b
KS
add x30,x30,#1536
KLKB \b
.endm
.macro OL
orr x10,x9,x10,lsl #8
LQ
.endm
.macro E4
BWDPK01
BWDPK00
BWDPK10
.endm
.macro LS
ldr x6,[sp,#-128]
sub sp,sp,#512
stp x6,x4,[sp,#-128]
JK
ldr x4,[sp,#-120]
JL
.endm
.macro ZB4
adds x20,x20,x12
adcs x21,x21,x13
adcs x22,x22,x14
adcs x23,x23,x15
.endm
.macro ZB3
adds x8,x8,x0
adcs x9,x9,x1
adc x10,x10,x2
adds x11,x11,x3
adc x12,x12,x4
.endm
.macro ZB2
ld64b x8,[x17]
AD5 x20,x21,x22,x23,x28,x12,x13,x14,x15,x20
.endm
.macro ZB1
adds x8,x8,x18
adcs x9,x9,x19
adcs x10,x10,x20
adcs x11,x11,x21
adc x0,x0,xzr
adds x26,x26,x22
adcs x27,x27,x23
adcs x28,x28,x24
adcs x29,x29,x25
.endm
.macro ZB5
mov x4,#977
movk x4,#1,lsl #32
.endm
.macro ZB0 b
ld64b x8,[\b]
LP x16
AD5 x20,x21,x22,x23,x2,x12,x13,x14,x15,x20
.endm
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
.macro RA4 o,a,b,c,d,e,f,g,h,i,j,y,r,s,t,v
mul \r,\o,\a
adds \e,\e,\r
umulh \s,\o,\a
mul \r,\o,\b
adcs \s,\s,\r
umulh \t,\o,\b
mul \r,\o,\c
adcs \t,\t,\r
umulh \v,\o,\c
mul \r,\o,\d
adcs \v,\v,\r
umulh \r,\o,\d
adc \r,\r,xzr
adds \f,\f,\s
adcs \g,\g,\t
adcs \h,\h,\v
adcs \i,\i,\r
adc \j,\y,xzr
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
.macro SHX l,a,d,p,q
mul \l,\d,x4
umulh \d,\d,x4
adds \a,\a,\d
b.lo 9f
adcs \p,\p,xzr
adcs \q,\q,xzr
SD \d
9:
.endm
.macro FDX e,f,i,j,o,p,l,h,y,z
mul \l,\i,x4
umulh \h,\i,x4
madd \h,\j,x4,\h
adds \o,\e,\l
adcs \p,\f,\h
b.lo 9f
adcs \y,\y,xzr
adc \z,\z,xzr
SD \l
9:
.endm
.macro FS a,q,c,d,e,f,g,h,i,j,k,l,t
subs \a,\e,\i
sbcs \q,\f,\j
sbcs \c,\g,\k
sbcs \d,\h,\l
csel \t,x4,xzr,lo
subs \a,\a,\t
b.hs 9f
sbcs \q,\q,xzr
sbcs \c,\c,xzr
sbc \d,\d,xzr
SD \t
9:
.endm
.macro SQ13 a,b,c,d,e,f,g,h,i,j,k,l,t
mul \h,\a,\d
mul \t,\b,\c
adds \h,\h,\t
umulh \i,\a,\d
umulh \t,\b,\c
adcs \i,\i,\t
umulh \j,\b,\d
adc \j,\j,xzr
umulh \g,\a,\b
mul \t,\a,\c
adds \g,\g,\t
umulh \t,\a,\c
adcs \h,\h,\t
mul \t,\b,\d
adcs \i,\i,\t
mul \t,\c,\d
adcs \j,\j,\t
umulh \k,\c,\d
adc \k,\k,xzr
mul \f,\a,\b
lsr \l,\k,#63
extr \k,\k,\j,#63
extr \j,\j,\i,#63
extr \i,\i,\h,#63
extr \h,\h,\g,#63
extr \g,\g,\f,#63
mul \e,\a,\a
umulh \t,\a,\a
adds \f,\t,\f,lsl #1
mul \t,\b,\b
adcs \g,\g,\t
umulh \t,\b,\b
adcs \h,\h,\t
mul \t,\c,\c
adcs \i,\i,\t
umulh \t,\c,\c
adcs \j,\j,\t
mul \t,\d,\d
adcs \k,\k,\t
umulh \t,\d,\d
adc \l,\l,\t
mul \t,\i,x4
adds \e,\e,\t
umulh \i,\i,x4
madd \i,\j,x4,\i
umulh \j,\j,x4
madd \j,\k,x4,\j
umulh \k,\k,x4
madd \k,\l,x4,\k
umulh \l,\l,x4
adcs \f,\f,\i
adcs \g,\g,\j
adcs \h,\h,\k
adc \l,\l,xzr
.endm
.macro LKFOLD4 a,q,d,e,c,s
madd \a,\c,x4,\s
cbhs \a,\s,9f
subs \q,\q,#-1
adcs \d,\d,xzr
adcs \e,\e,xzr
csel \s,x4,xzr,hs
add \a,\a,\s
SD \s
9:
.endm
.macro AD5 a,b,d,e,c,f,g,h,i,o
adds \o,\a,\f
adcs \b,\b,\g
adcs \d,\d,\h
adcs \e,\e,\i
adc \c,\c,xzr
.endm
.macro LD_OB3_2
ldp x7,x30,[x1,#48]
.endm
.macro LKP6_PAIR a,q
ldp x28,x29,[x1,#32]
ld64b x16,[x28]
ld64b x8,[x29]
adds x16,x16,x8
adcs x17,x17,x9
adcs x18,x18,x10
adcs x19,x19,x11
cset x0,hs
ZB4
cset x2,hs
\q
ZB0 x7
ZB0 x30
ZB0 x6
ld64b x8,[x1]
LP x12
AD5 x20,x21,x22,x23,x2,x8,x9,x10,x11,x8
LKFOLD4 x16,x17,x18,x19,x0,x12
LKFOLD4 x20,x21,x22,x23,x2,x8
.endm
.macro P1X a,q,c,d,e,f,g,h,i,j,k,l
RF \a,\q,\c,\d,\e,x12,x13,x14,x15,x11,x0,x2
SHX x5,\q,\e,\c,\d
RA \f,x5,\q,\c,\d,\a,x0,x8,x9,x10,x12,x13,x14,x15,x11,x2,xzr,x14,x15
SHX \e,x5,\d,\q,\c
RA \g,\e,x5,\q,\c,\a,x0,x8,x9,x10,x12,x13,x14,x15,x11,x2,x2,x14,x15
SHX \d,\e,\c,x5,\q
RA \h,\d,\e,x5,\q,\a,x0,x8,x9,x10,x12,x13,x14,x15,x11,x2,x2,\i,\j
FDX x12,x13,x11,x2,\k,\l,x5,x0,\i,\j
.endm
.macro PAIR7C a,b,c,d,e,f,g,h
ld64b x16,[x1]
RF \a,x0,x2,x14,x15,x25,x26,x27,x28,x12,x6,x24
RF \e,x0,x2,x14,x15,x8,x9,x10,x11,x29,x6,x24
SHX x5,x0,x15,x2,x14
RA4 \b,x5,x0,x2,x14,x25,x26,x27,x28,x12,x13,xzr,x6,x24,x16,x20
RA4 \f,x5,x0,x2,x14,x8,x9,x10,x11,x29,x30,xzr,x6,x24,x16,x20
SHX x15,x5,x14,x0,x2
RA4 \c,x15,x5,x0,x2,x25,x26,x27,x28,x12,x13,x13,x6,x24,x16,x20
RA4 \g,x15,x5,x0,x2,x8,x9,x10,x11,x29,x30,x30,x6,x24,x16,x20
SHX x14,x15,x2,x5,x0
RA4 \d,x14,x15,x5,x0,x25,x26,x27,x28,x12,x13,x13,x6,x24,x16,x20
RA4 \h,x14,x15,x5,x0,x8,x9,x10,x11,x29,x30,x30,x6,x24,x16,x20
FDX x25,x26,x12,x13,x25,x26,x6,x12,x27,x28
FDX x8,x9,x29,x30,x8,x9,x6,x29,x10,x11
.endm
.macro C1BSPLIT_BDAY a,q,c,d
ldp x16,x17,[sp,#160]
ldp x18,x19,[sp,#176]
ldp x21,x22,[x3,#-160]
ldp x23,x24,[x3,#-144]
FS x21,x22,x23,x24,x16,x17,x18,x19,x21,x22,x23,x24,x5
RF \a,x21,x22,x23,x24,x12,x13,x14,x15,x0,x16,x17
SHX x20,x21,x24,x22,x23
RA4 \q,x20,x21,x22,x23,x12,x13,x14,x15,x0,x2,xzr,x16,x17,x18,x19
SHX x24,x20,x23,x21,x22
RA4 \c,x24,x20,x21,x22,x12,x13,x14,x15,x0,x2,x2,x16,x17,x18,x19
SHX x23,x24,x22,x20,x21
RA4 \d,x23,x24,x20,x21,x12,x13,x14,x15,x0,x2,x2,x16,x17,x18,x19
mul x16,x0,x4
umulh x0,x0,x4
madd x0,x2,x4,x0
adds x12,x12,x16
adc x13,x13,x0
st64b x8,[x1]
ldp x24,x29,[x3,#176]
mul x16,x24,\a
umulh x17,x24,\a
mul x18,x24,\q
umulh x19,x24,\q
mul x20,x24,\c
umulh x21,x24,\c
mul x22,x24,\d
umulh x23,x24,\d
adds x8,x12,x16
adcs x9,x13,x17
adcs x10,x14,x19
adcs x11,x15,x21
adc x23,x23,xzr
adds x9,x9,x18
adcs x10,x10,x20
adcs x11,x11,x22
adc x23,x23,xzr
mul x16,x23,x4
umulh x23,x23,x4
adds x8,x8,x16
adc x9,x9,x23
.endm
.macro P1W6_ACC
ZB1
adc x2,x2,xzr
.endm
.macro QJ o,y
RA4 \o,x22,x23,x24,x25,x13,x14,x15,x16,x17,x12,\y,x6,x1,x30,x0
.endm
.macro P1_6WIN_BFI
SQ13 x8,x9,x10,x11,x0,x2,x5,x6,x1,x20,x21,x22,x23
mul x23,x22,x4
umulh x22,x22,x4
adds x0,x0,x23
adc x2,x2,x22
ld64b x12,[x29]
bfi x12,x0,#6,#43
extr x1,x2,x0,#43
bfi x13,x1,#6,#43
lsr x1,x2,#22
add x14,x14,x1,lsl #6
bfi x15,x5,#6,#43
extr x1,x6,x5,#43
bfi x16,x1,#6,#43
lsr x1,x6,#22
add x17,x17,x1,lsl #6
ld64b x22,[x12]
adds x8,x8,x22
adcs x9,x9,x23
adcs x10,x10,x24
adcs x11,x11,x25
cset x0,hs
ld64b x18,[x13]
ZB1
cset x2,hs
ld64b x18,[x14]
P1W6_ACC
ld64b x18,[x15]
P1W6_ACC
ld64b x18,[x16]
P1W6_ACC
ld64b x18,[x17]
P1W6_ACC
madd x30,x2,x4,x26
madd x2,x0,x4,x8
ld64b x18,[sp]
RF x30,x18,x19,x20,x21,x13,x14,x15,x16,x17,x6,x1
QJ x27,xzr
add x5,sp,#64
ld64b x18,[x5]
RA4 x28,x18,x19,x20,x21,x13,x14,x15,x16,x17,x12,x12,x6,x1,x30,x0
QJ x29,x12
mul x6,x17,x4
umulh x17,x17,x4
madd x17,x12,x4,x17
adds x13,x13,x6
adc x14,x14,x17
subs x28,x2,x13
sbcs x29,x9,x14
sbcs x30,x10,x15
sbcs x1,x11,x16
csel x6,x4,xzr,lo
sub x28,x28,x6
adds x0,x2,x13
adcs x2,x9,x14
adcs x5,x10,x15
adcs x6,x11,x16
csel x20,x4,xzr,hs
add x0,x0,x20
.endm
.macro BTBO
ldp x8,x9,[x3,#-192]
.endm
.macro BTBE
ldp x8,x9,[x3,#320]
.endm
.macro BTBASES
mov x7,#0xfff8000000000000
BTBO
.endm
.macro A4 a,b,c,d
adds x24,x24,\a
adcs x25,x25,\b
adcs x26,x26,\c
adc x27,x27,\d
.endm
.macro BLK6F a,q,c,d
orr x10,x7,\a,lsl #8
ld64b x20,[x10]
extr x10,\q,\a,#43
orr x10,x7,x10,lsl #8
LQ
adds x24,x20,x16
adcs x25,x21,x17
adcs x26,x22,x18
adc x27,x23,x19
extr x10,\c,\q,#22
orr x10,x8,x10,lsl #8
LQ
A4 x16,x17,x18,x19
orr x10,x8,\c,lsl #7
LQ
A4 x12,x13,x14,x15
extr x10,\d,\c,#44
OL
A4 x12,x13,x14,x15
lsr x10,\d,#23
OL
A4 x16,x17,x18,x19
.endm
.macro ROWSG a,q,c,d,y0,y1,y2,y3
umulh x10,\d,x4
adds \a,\a,x10
b.lo 9f
adcs \q,\q,xzr
AZ \c,\d,x10
9:
adds \y0,\a,x24
adcs \y1,\q,x25
adcs \y2,\c,x26
adc \y3,\d,x27
.endm
.macro QG r,t,a,j,i,q,v,c
umulh \t,\a,\j
madd \r,\a,\i,\t
umulh \t,\q,\v
add \r,\r,\t
madd \r,\q,\j,\r
madd \r,\c,\v,\r
.endm
.macro KOUT a,b
stp x16,x17,[x3,#\a]
stp x18,x19,[x3,#\b]
.endm
.macro BG b
P1_6WIN_BFI
\b
BLK6F x28,x29,x30,x1
ROWSG x28,x29,x30,x1,x20,x21,x22,x23
.endm

.macro KTFPF a,q,c,d,e,f,g,h,i,j,v,s,y0,y1,y2,y3,o0,o1,o2,o3
QG \o1,\o2,\a,\j,\i,\q,\v,\c
QG \o0,\o2,\a,\v,\j,\q,\s,\c
QG \o3,\o2,\f,\y2,\y3,\g,\y1,\h
QG \o2,\a,\f,\y1,\y2,\g,\y0,\h
sub \o3,\o3,\o1
sub \o2,\o2,\o0
umulh \a,\d,x26
add \o1,\o1,\a
madd \o1,\d,x27,\o1
umulh \a,\e,x25
add \o1,\o1,\a
madd \o1,\e,x26,\o1
umulh \a,\d,x25
add \o0,\o0,\a
madd \o0,\d,x26,\o0
umulh \a,\e,x24
add \o0,\o0,\a
madd \o0,\e,x25,\o0
.endm
.macro KTFPP a,q,c,d,e,f,g,h,i,j,v,s,T
QG x23,\T,\f,x22,x23,\g,x21,\h
QG x22,\T,\f,x21,x22,\g,x20,\h
QG x21,\g,\a,\j,\i,\q,\v,\c
QG x20,\g,\a,\v,\j,\q,\s,\c
sub x23,x23,x21
sub x22,x22,x20
umulh \a,\d,x26
add x21,x21,\a
madd x21,\d,x27,x21
umulh \a,\e,x25
add x21,x21,\a
madd x21,\e,x26,x21
umulh \a,\d,x25
add x20,x20,\a
madd x20,\d,x26,x20
umulh \a,\e,x24
add x20,x20,\a
madd x20,\e,x25,x20
.endm
.macro KM0 a,b,c,v,s
ld64b x8,[\a]
KTFPF x8,x9,x10,x11,x12,x13,x14,x15,\b,\c,\v,\s,x20,x21,x22,x23,x16,x17,x18,x19
.endm
.macro KMQ0 a,b,c,v,s
ld64b x8,[\a]
KTFPP x8,x9,x10,x11,x12,x13,x14,x15,\b,\c,\v,\s,\a
.endm
.macro KMY0 a,b,c,v,s,o0,o1,o2,o3
ld64b x8,[\a]
KTFPF x8,x9,x10,x11,x12,x13,x14,x15,\b,\c,\v,\s,x28,x29,x30,x1,\o0,\o1,\o2,\o3
.endm
.macro KMA a,b,c,v,s
ld64b x8,[\a]
KTFPF x11,x12,x13,x14,x15,x8,x9,x10,\b,\c,\v,\s,x20,x21,x22,x23,x16,x17,x18,x19
.endm
.macro KMQA a,b,c,v,s
ld64b x8,[\a]
KTFPP x11,x12,x13,x14,x15,x8,x9,x10,\b,\c,\v,\s,\a
.endm
.macro KMYA a,b,c,v,s,o0,o1,o2,o3
ld64b x8,[\a]
KTFPF x11,x12,x13,x14,x15,x8,x9,x10,\b,\c,\v,\s,x28,x29,x30,x1,\o0,\o1,\o2,\o3
.endm
.macro KMB a,b,c,v,s
ld64b x8,[\a]
KTFPF x10,x11,x12,x13,x14,x8,x9,x15,\b,\c,\v,\s,x20,x21,x22,x23,x16,x17,x18,x19
.endm
.macro KMQB a,b,c,v,s
ld64b x8,[\a]
KTFPP x10,x11,x12,x13,x14,x8,x9,x15,\b,\c,\v,\s,\a
.endm
.macro KMYB a,b,c,v,s,o0,o1,o2,o3
ld64b x8,[\a]
KTFPF x10,x11,x12,x13,x14,x8,x9,x15,\b,\c,\v,\s,x28,x29,x30,x1,\o0,\o1,\o2,\o3
.endm
.macro BH m,t
mov x7,#0xfff8000000000000
\t
BLK6F x0,x2,x5,x6
ROWSG x0,x2,x5,x6,x28,x29,x30,x1
KA 192
\m x16,x6,x5,x2,x0,x20,x21,x22,x23
KA 256
\m x16,x6,x5,x2,x0,x16,x17,x18,x19
add x8,x3,#192
st64b x16,[x8]
\m x3,x6,x5,x2,x0,x16,x17,x18,x19
.endm
.macro BHE m,t
mov x7,#0xfff8000000000000
\t
BLK6F x0,x2,x5,x6
ROWSG x0,x2,x5,x6,x28,x29,x30,x1
KA 192
\m x16,x6,x5,x2,x0,x20,x21,x22,x23
KA 256
\m x16,x6,x5,x2,x0,x16,x17,x18,x19
add x8,x3,#192
st64b x16,[x8]
ldp x20,x21,[sp,#192]
ldp x22,x23,[sp,#208]
\m x3,x6,x5,x2,x0,x16,x17,x18,x19
.endm
.macro BWDK b,n,t
BG \b
KM\n x3,x1,x30,x29,x28
KOUT 96,112
KA 192
KM\n x16,x1,x30,x29,x28
add x7,x3,#256
KMQ\n x7,x1,x30,x29,x28
sub x14,x3,#128
st64b x16,[x14]
BH KMY\n,\t
.endm
.macro BWDKE b,n,t
BG \b
KM\n x3,x1,x30,x29,x28
KOUT 96,112
KA 192
KM\n x16,x1,x30,x29,x28
add x7,x3,#256
KMQ\n x7,x1,x30,x29,x28
sub x14,x3,#128
st64b x16,[x14]
BHE KMY\n,\t
.endm
.macro BH2 m,t
mov x7,#0xfff8000000000000
\t
BLK6F x0,x2,x5,x6
ROWSG x0,x2,x5,x6,x28,x29,x30,x1
add x7,x3,#192
\m x7,x6,x5,x2,x0,x20,x21,x22,x23
KA 256
\m x16,x6,x5,x2,x0,x16,x17,x18,x19
st64b x16,[x7]
\m x3,x6,x5,x2,x0,x16,x17,x18,x19
.endm
.macro BWDK2 b,n,t
BG \b
KM\n x3,x1,x30,x29,x28
KOUT 96,112
KA 192
KM\n x16,x1,x30,x29,x28
add x7,x3,#256
KMQ\n x7,x1,x30,x29,x28
sub x14,x3,#128
st64b x16,[x14]
BH2 KMY\n,\t
.endm
.macro BWDPH
sub x1,x3,#704
PAIR7C x20,x21,x22,x23,x16,x17,x18,x19
C1BSPLIT_BDAY x25,x26,x27,x28
.endm
.macro BWDL
stp x18,x19,[x3,#16]
stp x16,x17,[x3],#-512
.endm
.macro BWDLM
BWDL
ldp x8,x9,[x3,#-160]
ldp x10,x11,[x3,#-144]
add x29,x3,#128
.endm
.macro BWDPK01
BWDPH
BWDK BTBASES,B,BTBO
BWDLM
BWDK BTBE,A,BTBE
BWDL
JL
.endm
.macro BWDPK00
BWDPH
BWDK BTBO,0,BTBO
BWDLM
BWDK BTBE,B,BTBE
BWDL
JL
.endm
.macro BWDPK10
BWDPH
BWDKE BTBO,A,BTBO
st64b x16,[x3]
sub x3,x3,#512
ldp x8,x9,[x3,#-160]
ldp x10,x11,[x3,#-144]
add x29,x3,#128
BWDK2 BTBE,0,BTBE
BWDL
.endm
.macro Z
adds x0,x8,x0
adcs x1,x9,x1
adc x2,x10,x2
adds x3,x11,x3
adc x4,x12,x4
.endm
.macro KACC8 b
\b
adds x13,x13,x5
adcs x14,x14,x6
adc x15,x15,x7
.endm
.macro KACC7B
adds x16,x16,x0
adc x17,x17,x1
adds x18,x18,x2
adcs x19,x19,x3
adc x20,x20,x4
adds x21,x21,x5
adc x22,x22,x6
.endm
.macro KWINA a,b
add x23,\a,#64
ld64b x0,[x23]
KACC7B
ld64b x0,[\a]
KACC8 \b
.endm
.macro KWINB a,b
add x23,\a,#192
ld64b x0,[x23]
KACC7B
sub x23,x23,#64
ld64b x0,[x23]
KACC8 \b
.endm
.macro KLKA a
ld64b x8,[x24]
add x0,x24,#64
ld64b x16,[x0]
KWINA x25,ZB3
KWINA x26,ZB3
KWINA x27,ZB3
KWINA x28,\a
.endm
.macro KLKB a
add x0,x24,#128
ld64b x8,[x0]
add x0,x0,#64
ld64b x16,[x0]
KWINB x25,ZB3
KWINB x26,ZB3
KWINB x27,ZB3
KWINB x28,\a
.endm
.macro KLKA5 a
ld64b x8,[x25]
add x0,x25,#64
ld64b x16,[x0]
KWINA x26,ZB3
KWINA x27,ZB3
KWINA x28,ZB3
KWINA x29,\a
.endm
.macro KLKB5 a
add x0,x25,#128
ld64b x8,[x0]
add x0,x0,#64
ld64b x16,[x0]
KWINB x26,ZB3
KWINB x27,ZB3
KWINB x28,ZB3
KWINB x29,\a
.endm
.macro KK5 b
KS
add x30,x30,#1536
KLKB5 \b
.endm
.macro KEB
ldp x0,x1,[sp,#240]
ldr x2,[sp,#256]
ldp x3,x4,[sp,#264]
ZB3
adds x13,x13,x0
adcs x14,x14,x1
adc x15,x15,x2
adds x16,x16,x3
adc x17,x17,x4
adds x18,x18,x0
adcs x19,x19,x1
adc x20,x20,x2
adds x21,x21,x3
adc x22,x22,x4
.endm
.macro KS
subs x5,x0,x3
sbcs x6,x1,x4
sbc x7,x2,xzr
st64b x0,[x30]
add x9,x30,#512
subs x10,x13,x16
sbcs x11,x14,x17
sbc x12,x15,xzr
st64b x10,[x9]
add x3,x30,#1024
subs x16,x18,x21
sbcs x17,x19,x22
sbc x23,x20,xzr
st64b x16,[x3]
.endm
.macro KOFFS a,b,c
ldp x24,x25,[sp,#\a]
ldp x26,x27,[sp,#\b]
ldp x28,x29,[sp,#\c]
.endm
.macro EFGK6
stp x0,x2,[sp,#-120]
stp x14,x15,[sp,#-104]
JK
ZB5
ldp x0,x2,[sp,#-120]
ldp x14,x15,[sp,#-104]
.endm
.macro JK
sub x30,x3,#2368
KOFFS 128,144,320
KLKA ZB3
KEB
KWINA x29,Z
KK ZB3
KEB
KWINB x29,Z
KS
sub x30,x30,#1472
add x0,sp,#256
ld64b x22,[x0]
KLKA5 Z
KK5 Z
KS
sub x30,x30,#1792
add x0,sp,#320
ld64b x22,[x0]
KLKA ZB3
KWINA x29,Z
KK ZB3
KWINB x29,Z
KS
.endm
.macro JL
ldp x0,x2,[x3,#320]
ldp x14,x15,[x3,#336]
.endm
.macro FLUSH7 a,s,r
sub \s,\a,#4064
str \s,[x0]
ldr \r,[\s,#-16]
add \r,\r,#1
mov \a,#4096
stp \a,\r,[\s,#-24]
stlur x0,[\s,#-32]
dmb ish
.endm
.macro EM s,p,q,r,u,v
ldadd x4,\s,[x0]
bfxil \p,\q,#20,#21
bfxil \q,\r,#23,#41
stp \q,\p,[\s]
stp x2,x30,[\s,#16]
tbz \s,#12,kn\@
FLUSH7 \s,\u,\v
kn\@:
.endm
.macro LJ0C k,a,b,c,d
add x5,\b,\d
mvn x5,x5
cbhs \b,\d,kh\k
cbhi \b,x5,kv\k
and x0,x1,\a,lsr #21
EM x5,\d,\b,\c,\b,\d
.endm
.macro SLH k,b,c,d,e
kh\k:
add x0,\b,\d
cbhi \e,x5,kt\k
mvn \b,x0
and x0,x1,\d,lsr #21
EM \c,x5,\e,\b,x5,\e
.endm
.macro SLV k,b,c,d,e
kv\k:
add x0,\b,\d
kt\k:
bic x0,x1,x0,lsr #21
EM \e,\c,x5,\b,x5,\c
.endm
.macro T3 k,e,f,g,h,z
add \g,\g,\e,lsl #1
add \h,\h,\f,lsl #1
LJ0C e\k,\g,\h,x22,x23
\z re\k
.endm
.macro T2 k,a,b,c,d,e,f,g,h,z
LJ0C r\k,\a,\b,\c,\d
T3 \k,\e,\f,\g,\h,\z
.endm
.macro KZ k
sub x2,x2,#5
add x3,x3,#512
tbz x3,#37,1b
b 2f
.endm
.macro NETK k,a,b,c,d,e,f,g,h,z
add x0,\c,\a,lsl #1
add x5,\d,\b,lsl #1
sub \a,\a,\g
sub \b,\b,\h
sub \c,\c,x26
sub \d,\d,x27
sub x22,x26,\a
sub x23,x27,\b
sub x24,\e,\c
sub x25,\f,\d
sub \a,\a,x24
sub \b,\b,x25
sub x28,\a,x28
sub x29,\b,x29
add \g,\g,x28
add \h,\h,x29
sub x28,x28,x0
sub x29,x29,x5
add \a,\a,x26,lsl #1
add \b,\b,x27,lsl #1
sub \c,\c,\g
sub \d,\d,\h
LJ0C \k,x28,x29,x24,x25
T2 \k,\a,\b,\c,\d,\e,\f,\g,\h,\z
SLH \k,x28,x29,x24,x25
T2 \k,\a,\b,\c,\d,\e,\f,\g,\h,\z
SLV \k,x28,x29,x24,x25
T2 \k,\a,\b,\c,\d,\e,\f,\g,\h,\z
SLH r\k,\a,\b,\c,\d
T3 \k,\e,\f,\g,\h,\z
SLV r\k,\a,\b,\c,\d
T3 \k,\e,\f,\g,\h,\z
SLH e\k,\g,\h,x22,x23
\z re\k
SLV e\k,\g,\h,x22,x23
re\k:
.endm
.macro KAD a,b,c,d
add x\a,x\a,x22
add x\b,x\b,x23
add x\c,x\c,x24
add x\d,x\d,x25
.endm
.macro KFC
KAD 10,11,12,13
KAD 14,15,16,17
.endm
.macro LK1 a
orr x9,x28,#\a
ld64b x16,[x9]
orr x17,x29,#\a
ld64b x8,[x17]
ZB4
cset x28,hs
orr x17,x7,#\a
ZB2
orr x17,x30,#\a
ZB2
orr x17,x6,#\a
ZB2
ld64b x8,[x1]
AD5 x20,x21,x22,x23,x28,x8,x9,x10,x11,x8
LKFOLD4 x20,x21,x22,x23,x28,x8
.endm
.macro FWD6_ONE
ldp x6,x7,[x1,#288]
str x7,[x3,#368]
LKP6_PAIR 0x4c000000000000,LD_OB3_2
st64b x20,[x3]
P1X x20,x24,x25,x26,x27,x21,x22,x23,x22,x23,x20,x21
add x5,x3,#1024
st64b x16,[x5]
P1X x16,x20,x21,x22,x23,x17,x18,x19,x26,x27,x24,x25
LK1 0x2000000000000
add x5,x3,#2048
st64b x20,[x5]
P1X x20,x24,x25,x26,x27,x21,x22,x23,x26,x27,x24,x25
add x1,x1,#512
add x3,x3,#3072
.endm
.macro KEC b,f,g
ldp x22,x23,[\b,#\f]
ldp x24,x25,[\b,#\g]
KAD 6,7,8,9
KAD 18,19,20,21
.endm
.macro KSB t
sub x0,x3,#\t
.endm
.macro KPH k,o,n,s,t,b,f,g,z
ld64b x6,[x3]
add x0,x3,#320
ld64b x14,[x0]
add x0,x3,#192
ld64b x22,[x0]
KFC
\s \t
KEC \b,\f,\g
NETK p\k,x10,x11,x12,x13,x6,x7,x8,x9,b
ldp x26,x27,[x3,#128]
ldp x28,x29,[x3,#144]
NETK m\k,x14,x15,x16,x17,x18,x19,x20,x21,\z
\o x2,x2,#\n
add x3,x3,#512
.endm
.globl fwd
fwd:
ZB5
stp x29,x30,[sp,#-16]!
str xzr,[sp,#-128]
and x5,x0,#6
cbz x0,4f
cmp x5,#4
b.eq 6f
b.hs 7f
cbnz x5,5f
1:
FWD6_ONE
5:
FWD6_ONE
6:
FWD6_ONE
7:
FWD6_ONE
tbz x3,#37,1b
4:
ldr x0,[sp,#-128]
ldp x29,x30,[sp],#16
ret
.globl bwd
bwd:
mov x5,sp
stp x5,x30,[x1]
str xzr,[sp],#128
sub x3,x3,#320
EFGK6
mov x5,sp
tbnz x5,#10,5f
tbnz x5,#9,9f
b 4f
9:
b 8f
5:
tbz x5,#9,7f
b 6f
1:
LS
6:
E4
LS
7:
E4
LS
8:
E4
LS
4:
E4
mov x5,sp
tst x5,#0x7fffff00
b.ne 1b
JL
ldr x6,[sp,#-128]
ldp x5,x30,[sp,#-160]
mov sp,x5
mov x5,x6
ret
.globl keys
keys:
str x30,[sp,#-16]!
mov x1,#0x7fffffffff8
mov x4,#32
sub x2,x2,#6
1:
ldr x30,[x3,#304]
KPH 0,add,1,,,x0,480,496,b
KPH 1,add,1,,,x3,160,176,b
KPH 2,add,1,,,x3,-352,-336,b
KPH 3,add,1,KSB,1024,x0,160,176,b
KPH 4,add,1,KSB,1536,x0,160,176,b
KPH 5,sub,5,KSB,2048,x0,160,176,KZ
tbz x3,#37,1b
2:
add x2,x2,#6
ldr x30,[sp],#16
ret

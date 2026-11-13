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
stadd x4,[sp]
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
stadd x4,[sp]
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
stadd x4,[sp]
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
stadd x4,[sp]
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
stadd x4,[sp]
9:
.endm
.macro LKFOLD4 o0,o1,o2,o3,c,t
mul \t,\c,x4
adds \o0,\o0,\t
b.lo 9f
adcs \o1,\o1,xzr
adcs \o2,\o2,xzr
adcs \o3,\o3,xzr
csel \t,x4,xzr,hs
add \o0,\o0,\t
stadd x4,[sp]
9:
.endm
.macro LKP_ADD5 a0,a1,a2,a3,c,b0,b1,b2,b3
adds \a0,\a0,\b0
adcs \a1,\a1,\b1
adcs \a2,\a2,\b2
adcs \a3,\a3,\b3
adc \c,\c,xzr
.endm
.macro LD_OB3_2
ldp x5,x30,[x1,#48]
.endm
.macro LKP6_PAIR tb,ld3_ins
mov x2,#\tb
ldp x28,x29,[x1,#32]
add x5,x2,x28
ld64b x16,[x5]
add x5,x2,x29
ld64b x8,[x5]
adds x16,x16,x8
adcs x17,x17,x9
adcs x18,x18,x10
adcs x19,x19,x11
cset x28,hs
adds x20,x20,x12
adcs x21,x21,x13
adcs x22,x22,x14
adcs x23,x23,x15
cset x29,hs
\ld3_ins
add x5,x2,x5
ld64b x8,[x5]
LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
add x5,x2,x30
ld64b x8,[x5]
LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
add x5,x2,x6
ld64b x8,[x5]
LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
LKP_ADD5 x20,x21,x22,x23,x29,x12,x13,x14,x15
ld64b x8,[x1]
LKP_ADD5 x16,x17,x18,x19,x28,x8,x9,x10,x11
LKP_ADD5 x20,x21,x22,x23,x29,x8,x9,x10,x11
LKFOLD4 x16,x17,x18,x19,x28,x8
LKFOLD4 x20,x21,x22,x23,x29,x8
.endm
.macro MUL4 a0,a1,a2,a3,b0,b1,b2,b3,o0,o1,o2,o3,e0,e1,e4,e5,t0,t1,t2,t3,sh
RF \a0,\b0,\b1,\b2,\b3,\e0,\e1,\o2,\o3,\e4,\t0,\t1
SHX \sh,\b0,\b3,\b1,\b2
RA4 \a1,\sh,\b0,\b1,\b2,\e0,\e1,\o2,\o3,\e4,\e5,xzr,\t0,\t1,\t2,\t3
SHX \b3,\sh,\b2,\b0,\b1
RA4 \a2,\b3,\sh,\b0,\b1,\e0,\e1,\o2,\o3,\e4,\e5,\e5,\t0,\t1,\t2,\t3
SHX \b2,\b3,\b1,\sh,\b0
RA4 \a3,\b2,\b3,\sh,\b0,\e0,\e1,\o2,\o3,\e4,\e5,\e5,\t0,\t1,\t2,\t3
FDX \e0,\e1,\e4,\e5,\o0,\o1,\t0,\e4,\o2,\o3
.endm
.macro P1X a0,b0,b1,b2,b3,a1,a2,a3,o2,o3,o0,o1
RF \a0,\b0,\b1,\b2,\b3,x12,x13,x14,x15,x11,x28,x29
SHX x5,\b0,\b3,\b1,\b2
RA \a1,x5,\b0,\b1,\b2,\a0,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,xzr,x14,x15
SHX \b3,x5,\b2,\b0,\b1
RA \a2,\b3,x5,\b0,\b1,\a0,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,x29,x14,x15
SHX \b2,\b3,\b1,x5,\b0
RA \a3,\b2,\b3,x5,\b0,\a0,x28,x8,x9,x10,x12,x13,x14,x15,x11,x29,x29,\o2,\o3
FDX x12,x13,x11,x29,\o0,\o1,x5,x28,\o2,\o3
.endm
.macro PAIR7C k0,k1,k2,k3,s0,s1,s2,s3
ld64b x16,[x1]
RF \k0,x0,x2,x14,x15,x25,x26,x27,x28,x12,x10,x11
RF \s0,x0,x2,x14,x15,x6,x7,x8,x9,x29,x10,x11
SHX x5,x0,x15,x2,x14
RA4 \k1,x5,x0,x2,x14,x25,x26,x27,x28,x12,x13,xzr,x10,x11,x24,x16
RA4 \s1,x5,x0,x2,x14,x6,x7,x8,x9,x29,x30,xzr,x10,x11,x24,x16
SHX x15,x5,x14,x0,x2
RA4 \k2,x15,x5,x0,x2,x25,x26,x27,x28,x12,x13,x13,x10,x11,x24,x16
RA4 \s2,x15,x5,x0,x2,x6,x7,x8,x9,x29,x30,x30,x10,x11,x24,x16
SHX x14,x15,x2,x5,x0
RA4 \k3,x14,x15,x5,x0,x25,x26,x27,x28,x12,x13,x13,x10,x11,x24,x16
RA4 \s3,x14,x15,x5,x0,x6,x7,x8,x9,x29,x30,x30,x10,x11,x24,x16
FDX x25,x26,x12,x13,x25,x26,x10,x12,x27,x28
FDX x6,x7,x29,x30,x6,x7,x10,x29,x8,x9
.endm
.macro C1BSPLIT_BDAY p0,p1,p2,p3
ldp x16,x17,[sp,#288]
ldp x18,x19,[sp,#304]
ldp x21,x22,[x3,#224]
ldp x23,x24,[x3,#240]
FS x21,x22,x23,x24,x16,x17,x18,x19,x21,x22,x23,x24,x5
MUL4 \p0,\p1,\p2,\p3,x21,x22,x23,x24,x0,x2,x5,x6,x0,x2,x16,x17,x9,x10,x11,x12,x13
stp x0,x2,[x3,#-64]
stp x5,x6,[x3,#-48]
ldp x24,x29,[x3,#304]
mul x16,x24,\p0
umulh x17,x24,\p0
mul x18,x24,\p1
umulh x19,x24,\p1
mul x20,x24,\p2
umulh x21,x24,\p2
mul x22,x24,\p3
umulh x23,x24,\p3
adds x8,x0,x16
adcs x9,x2,x17
adcs x10,x5,x19
adcs x11,x6,x21
adc x23,x23,xzr
adds x9,x9,x18
adcs x10,x10,x20
adcs x11,x11,x22
adc x23,x23,xzr
DFX x8,x9,x10,x11,x23,x16
.endm
.macro DFX o0,o1,o2,o3,c,t
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
stadd x4,[sp]
9:
.endm
.macro P1W6_ACC
adds x8,x8,x18
adcs x9,x9,x19
adcs x10,x10,x20
adcs x11,x11,x21
adc x0,x0,xzr
adds x26,x26,x22
adcs x27,x27,x23
adcs x28,x28,x24
adcs x29,x29,x25
adc x2,x2,xzr
.endm
.macro P1_6WIN_BFI
SQ13 x8,x9,x10,x11,x0,x2,x5,x6,x7,x20,x21,x22,x23
F5 x0,x2,x5,x6,x22,x0,x2,x23,x22
ld64b x12,[x29]
bfi x12,x0,#6,#43
extr x7,x2,x0,#43
bfi x13,x7,#6,#43
lsr x7,x2,#22
add x14,x14,x7,lsl #6
bfi x15,x5,#6,#43
extr x7,x6,x5,#43
bfi x16,x7,#6,#43
lsr x7,x6,#22
add x17,x17,x7,lsl #6
ld64b x22,[x12]
adds x8,x8,x22
adcs x9,x9,x23
adcs x10,x10,x24
adcs x11,x11,x25
cset x0,hs
ld64b x18,[x13]
adds x8,x8,x18
adcs x9,x9,x19
adcs x10,x10,x20
adcs x11,x11,x21
adc x0,x0,xzr
adds x26,x26,x22
adcs x27,x27,x23
adcs x28,x28,x24
adcs x29,x29,x25
cset x2,hs
ld64b x18,[x14]
P1W6_ACC
ld64b x18,[x15]
P1W6_ACC
ld64b x18,[x16]
P1W6_ACC
ld64b x18,[x17]
P1W6_ACC
LKFOLD4 x8,x9,x10,x11,x0,x18
LKFOLD4 x26,x27,x28,x29,x2,x18
add x5,sp,#128
ld64b x18,[x5]
RF x26,x18,x19,x20,x21,x13,x14,x15,x16,x17,x6,x7
RA4 x27,x22,x23,x24,x25,x13,x14,x15,x16,x17,x12,xzr,x6,x7,x30,x0
add x5,x5,#64
ld64b x18,[x5]
RA4 x28,x18,x19,x20,x21,x13,x14,x15,x16,x17,x12,x12,x6,x7,x30,x0
RA4 x29,x22,x23,x24,x25,x13,x14,x15,x16,x17,x12,x12,x6,x7,x30,x0
FDX x13,x14,x17,x12,x13,x14,x6,x17,x15,x16
FS x28,x29,x30,x1,x8,x9,x10,x11,x13,x14,x15,x16,x6
FADDW x0,x2,x5,x6,x8,x9,x10,x11,x13,x14,x15,x16,x20
.endm
.macro BTBASES
mov x7,#0xfff8000000000000
orr x8,x7,#64
orr x9,x7,#128
.endm
.macro A4 a,b,c,d
adds x24,x24,\a
adcs x25,x25,\b
adcs x26,x26,\c
adc x27,x27,\d
.endm
.macro BLK6F s0,s1,s2,s3
orr x10,x7,\s0,lsl #8
ld64b x20,[x10]
extr x10,\s1,\s0,#43
orr x10,x7,x10,lsl #8
ld64b x12,[x10]
adds x24,x20,x16
adcs x25,x21,x17
adcs x26,x22,x18
adc x27,x23,x19
extr x10,\s2,\s1,#22
orr x10,x8,x10,lsl #8
ld64b x12,[x10]
A4 x16,x17,x18,x19
orr x10,x8,\s2,lsl #7
ld64b x12,[x10]
A4 x12,x13,x14,x15
extr x10,\s3,\s2,#44
orr x10,x9,x10,lsl #8
ld64b x12,[x10]
A4 x12,x13,x14,x15
lsr x10,\s3,#23
orr x10,x9,x10,lsl #8
ld64b x12,[x10]
A4 x16,x17,x18,x19
.endm
.macro ROWSG w0,w1,w2,w3
umulh x7,\w3,x4
adds \w0,\w0,x7
adcs \w1,\w1,xzr
b.lo 9f
adcs \w2,\w2,xzr
adc \w3,\w3,xzr
stadd x4,[sp]
9:
adds x20,\w0,x24
adcs x21,\w1,x25
adcs x22,\w2,x26
adc x23,\w3,x27
.endm
.macro KTFPF m0,m1,mt,n0,n1,d0,d1,dt,a3,a2,v,s,AH,AL,BH,BL
umulh x7,\m0,\a2
madd \AH,\m0,\a3,x7
umulh x7,\m1,\v
add \AH,\AH,x7
madd \AH,\m1,\a2,\AH
madd \AH,\mt,\v,\AH
umulh x7,\m0,\v
madd \AL,\m0,\a2,x7
umulh x7,\m1,\s
add \AL,\AL,x7
madd \AL,\m1,\v,\AL
madd \AL,\mt,\s,\AL
umulh x7,\d0,x22
madd \BH,\d0,x23,x7
umulh x7,\d1,x21
add \BH,\BH,x7
madd \BH,\d1,x22,\BH
madd \BH,\dt,x21,\BH
umulh x7,\d0,x21
madd \BL,\d0,x22,x7
umulh x7,\d1,x20
add \BL,\BL,x7
madd \BL,\d1,x21,\BL
madd \BL,\dt,x20,\BL
sub \BH,\BH,\AH
sub \BL,\BL,\AL
umulh x7,\n0,x26
add \AH,\AH,x7
madd \AH,\n0,x27,\AH
umulh x7,\n1,x25
add \AH,\AH,x7
madd \AH,\n1,x26,\AH
umulh x7,\n0,x25
add \AL,\AL,x7
madd \AL,\n0,x26,\AL
umulh x7,\n1,x24
add \AL,\AL,x7
madd \AL,\n1,x25,\AL
.endm
.macro KU0 lo,a3,a2,v,s
add x7,x3,#\lo
ld64b x8,[x7]
subs x13,x8,x11
sbcs x14,x9,x12
sbc x15,x10,xzr
stp x13,x14,[x7,#40]
str x15,[x7,#56]
KTFPF x8,x9,x10,x11,x12,x13,x14,x15,\a3,\a2,\v,\s,x17,x16,x19,x18
.endm
.macro KU1 lo,a3,a2,v,s
add x7,x3,#\lo
ld64b x8,[x7]
subs x8,x9,x12
sbcs x14,x10,x13
sbc x15,x11,xzr
str x8,[x7]
stp x14,x15,[x7,#48]
KTFPF x9,x10,x11,x12,x13,x8,x14,x15,\a3,\a2,\v,\s,x17,x16,x19,x18
.endm
.macro KM0 lo,a3,a2,v,s
add x7,x3,#\lo
ld64b x8,[x7]
KTFPF x8,x9,x10,x11,x12,x13,x14,x15,\a3,\a2,\v,\s,x17,x16,x19,x18
.endm
.macro KM1 lo,a3,a2,v,s
add x7,x3,#\lo
ld64b x8,[x7]
KTFPF x9,x10,x11,x12,x13,x8,x14,x15,\a3,\a2,\v,\s,x17,x16,x19,x18
.endm
.macro KCORS
ldp x8,x9,[sp,#320]
ldp x10,x11,[sp,#336]
add x16,x16,x8
add x17,x17,x9
add x18,x18,x10
add x19,x19,x11
.endm
.macro KOUT a,b
stp x16,x17,[x3,#\a]
stp x18,x19,[x3,#\b]
.endm
.macro BWDK u,m
P1_6WIN_BFI
BTBASES
BLK6F x28,x29,x30,x1
ROWSG x28,x29,x30,x1
\u 320,x1,x30,x29,x28
KCORS
KOUT 0,16
\u 384,x1,x30,x29,x28
KOUT 64,80
\u 128,x1,x30,x29,x28
KOUT -32,-16
BTBASES
BLK6F x0,x2,x5,x6
ROWSG x0,x2,x5,x6
\m 320,x6,x5,x2,x0
KCORS
KOUT 32,48
\m 384,x6,x5,x2,x0
KOUT 96,112
\m 128,x6,x5,x2,x0
KOUT 128,144
.endm
.macro BWDK0
BWDK KU0,KM0
.endm
.macro BWDK1
BWDK KU1,KM1
.endm
.macro BWDPH
sub x1,x3,#576
PAIR7C x20,x21,x22,x23,x16,x17,x18,x19
st64b x6,[x1]
C1BSPLIT_BDAY x25,x26,x27,x28
.endm
.macro BWDL a,b,c,d
ldp \a,\b,[x3,#-64]
ldp \c,\d,[x3,#-48]
sub x3,x3,#512
.endm
.macro BWDLM a,b,c,d
BWDL \a,\b,\c,\d
add x29,x3,#256
.endm
.macro BWDPK01
BWDPH
BWDK0
BWDLM x8,x9,x10,x11
BWDK1
BWDL x0,x2,x14,x15
.endm
.macro BWDPK00
BWDPH
BWDK0
BWDLM x8,x9,x10,x11
BWDK0
BWDL x0,x2,x14,x15
.endm
.macro BWDPK10
BWDPH
BWDK1
BWDLM x8,x9,x10,x11
BWDK0
BWDL x0,x2,x14,x15
.endm
.macro KACC8
adds x8,x8,x0
adcs x9,x9,x1
adc x10,x10,x2
adds x11,x11,x3
adc x12,x12,x4
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
.macro KWINA a
ld64b x0,[\a]
KACC8
add x23,\a,#64
ld64b x0,[x23]
KACC7B
.endm
.macro KWINB a
add x23,\a,#128
ld64b x0,[x23]
KACC8
add x23,x23,#64
ld64b x0,[x23]
KACC7B
.endm
.macro KLKA
ld64b x8,[x24]
add x0,x24,#64
ld64b x16,[x0]
KWINA x25
KWINA x26
KWINA x27
KWINA x28
.endm
.macro KLKB
add x0,x24,#128
ld64b x8,[x0]
add x0,x0,#64
ld64b x16,[x0]
KWINB x25
KWINB x26
KWINB x27
KWINB x28
.endm
.macro KEB
ldp x0,x1,[sp,#368]
ldr x2,[sp,#384]
ldp x3,x4,[sp,#392]
adds x8,x8,x0
adcs x9,x9,x1
adc x10,x10,x2
adds x11,x11,x3
adc x12,x12,x4
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
.macro KST3
st64b x8,[x30]
add x0,x30,#512
st64b x12,[x0]
add x0,x30,#1024
st64b x18,[x0]
.endm
.macro KOFFS a,b,c
ldp x24,x25,[sp,#\a]
ldp x26,x27,[sp,#\b]
ldp x28,x29,[sp,#\c]
.endm
.macro EFGK6
stp x0,x2,[sp,#8]
stp x14,x15,[sp,#24]
sub x30,x3,#2240
KOFFS 256,272,448
KLKA
KWINA x29
KEB
KST3
add x30,x30,#1536
KLKB
KWINB x29
KEB
KST3
sub x30,x30,#1472
KOFFS 408,424,440
KLKA
KST3
add x30,x30,#1536
KLKB
KST3
sub x30,x30,#1792
KOFFS 464,480,496
KLKA
KWINA x29
KST3
add x30,x30,#1536
KLKB
KWINB x29
KST3
add x3,x30,#896
mov x4,#977
movk x4,#1,lsl #32
ldp x0,x2,[sp,#8]
ldp x14,x15,[sp,#24]
.endm
.macro FLUSH7 slot,s,r
sub \s,\slot,#4064
str \s,[x0]
ldr \r,[\s,#-16]
add \r,\r,#1
mov \slot,#4096
stp \slot,\r,[\s,#-24]
stlur x0,[\s,#-32]
dmb ish
.endm
.macro LJ0C k,w0l,w0h,w1l,w1h
add x5,\w0h,\w1h
mvn x5,x5
cbhs \w0h,\w1h,kh\k
cbhi \w0h,x5,kv\k
and x0,x1,\w0l,lsr #21
ldadd x4,x5,[x0]
bfxil \w1h,\w0h,#20,#21
bfxil \w0h,\w1l,#23,#41
stp \w0h,\w1h,[x5]
stp x2,x26,[x5,#16]
tbz x5,#12,kn\k
FLUSH7 x5,\w0h,\w1h
kn\k:
.endm
.macro SLOWH k,ret,w0l,w0h,w1l,w1h
kh\k:
add x0,\w0l,\w1l
cbhi \w1h,x5,kt\k
mvn \w0l,x0
and x0,x1,\w1l,lsr #21
ldadd x4,\w0h,[x0]
bfxil x5,\w1h,#20,#21
bfxil \w1h,\w0l,#23,#41
stp \w1h,x5,[\w0h]
stp x2,x26,[\w0h,#16]
tbz \w0h,#12,km\k
FLUSH7 \w0h,x5,\w1h
km\k:
b \ret
kv\k:
add x0,\w0l,\w1l
kt\k:
bic x0,x1,x0,lsr #21
ldadd x4,\w1h,[x0]
bfxil \w0h,x5,#20,#21
bfxil x5,\w0l,#23,#41
stp x5,\w0h,[\w1h]
stp x2,x26,[\w1h,#16]
tbz \w1h,#12,ku\k
FLUSH7 \w1h,x5,\w0h
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
ra\k:
add \a0l,\a0l,\g0l,lsl #1
add \a0h,\a0h,\g0h,lsl #1
sub \a1l,\a1l,\b1l
sub \a1h,\a1h,\b1h
LJ0C d\k,\a0l,\a0h,\a1l,\a1h
rb\k:
add \b1l,\b1l,\b0l,lsl #1
add \b1h,\b1h,\b0h,lsl #1
LJ0C e\k,\b1l,\b1h,x27,x28
b re\k
SLOWH \k,ra\k,x24,x25,x29,x30
b ra\k
SLOWH d\k,rb\k,\a0l,\a0h,\a1l,\a1h
b rb\k
SLOWH e\k,re\k,\b1l,\b1h,x27,x28
re\k:
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
LKP6_PAIR 0x4c000000000000,LD_OB3_2
st64b x20,[x3]
P1X x20,x24,x25,x26,x27,x21,x22,x23,x22,x23,x20,x21
add x2,x3,#1024
st64b x16,[x2]
P1X x16,x20,x21,x22,x23,x17,x18,x19,x26,x27,x24,x25
LK1 0x56000000000000
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
str xzr,[sp,#-16]!
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
ldr x0,[sp],#16
ldp x29,x30,[sp],#16
ret
.globl bwd
bwd:
mov x4,#977
movk x4,#1,lsl #32
mov x5,sp
stp x5,x30,[x1]
str xzr,[sp]
sub x3,x3,#448
tbz x5,#9,3f
1:
EFGK6
BWDPK01
BWDPK00
BWDPK10
ldr x6,[sp]
stp x6,xzr,[sp,#-512]!
3:
EFGK6
BWDPK01
BWDPK00
BWDPK10
mov x5,sp
tst x5,#0x7fffffff
b.eq 2f
ldr x6,[sp]
stp x6,xzr,[sp,#-512]!
b 1b
2:
ldr x6,[sp]
ldp x5,x30,[sp,#-32]
mov sp,x5
mov x5,x6
ret
.globl keys
keys:
str x30,[sp,#-16]!
mov x1,#0x7fffffffff8
mov x4,#32
1:
ldr x26,[x3,#224]
sub x2,x2,#6
KEYPAIRC p00,m00,p10,m10
KEYPAIRC p20,m20,p30,m30
KEYPAIRC p40,m40,p50,m50
tbz x3,#37,1b
ldr x30,[sp],#16
ret

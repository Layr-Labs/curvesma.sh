// Leaf kernels (proof.md 3): fwd (PFWD2: EFWD2, P1X, P1XB), bwd (BWDP: LK2, W2S, C1BS, two INC: P1, PAIRB, TFP4S), keys.
// x4 = C0; SHX, FDX, F5, FS, FADDW, BLK6, LF, DF branch over their fixup (3.4; BLK6, LF on limb 0's carry: LF3; DF on limb 1's).
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
RF x20,x24,x25,x26,x27,x13,x14,x15,x16,x17,x28,x19
SHX x5,x24,x27,x25,x26
RA x21,x5,x24,x25,x26,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,xzr,x15,x16
SHX x27,x5,x26,x24,x25
RA x22,x27,x5,x24,x25,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,x29,x15,x16
SHX x26,x27,x25,x5,x24
RA x23,x26,x27,x5,x24,x6,x7,x8,x20,x21,x13,x14,x15,x16,x17,x29,x29,x26,x27
FDX x13,x14,x17,x29,x24,x25,x5,x6,x26,x27
.endm
.macro P1XB
RF x2,x24,x25,x26,x27,x13,x14,x15,x16,x17,x28,x19
SHX x5,x24,x27,x25,x26
RA x9,x5,x24,x25,x26,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,xzr,x15,x16
SHX x27,x5,x26,x24,x25
RA x30,x27,x5,x24,x25,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,x29,x15,x16
SHX x26,x27,x25,x5,x24
RA x18,x26,x27,x5,x24,x6,x7,x8,x20,x21,x13,x14,x15,x16,x17,x29,x29,x26,x27
FDX x13,x14,x17,x29,x24,x25,x5,x6,x26,x27
.endm
// count by C0 (BLK-type: the rare fixup finishes limbs 2, 3, adding C0 on their carry)
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
.macro LF a0,a1,a2,a3,c,t
mul \t,\c,x4
adds \a0,\a0,\t
b.lo 9f
adcs \a1,\a1,xzr
adcs \a2,\a2,xzr
adcs \a3,\a3,xzr
csel \t,x4,xzr,hs
add \a0,\a0,\t
stadd x4,[sp]
9:
.endm
// PX fwd: sigma's line [pair 2q's | 2q+1's] serves both e (x20..x23, x5; x2, x9, x30, x18, x6); B^2 loaded once
.macro W1F a
add x29,x7,\a
ld64b x10,[x29]
AC4 x20,x21,x22,x23,x10,x11,x12,x13,x5
AC4 x2,x9,x30,x18,x14,x15,x16,x17,x6
.endm
.macro W2F
W1F x8
W1F x28
.endm
.macro EFWD2
ldp x7,x8,[x3,#-280]
ldp x20,x21,[x1,#48]
ldp x22,x23,[x1,#64]
ldp x8,x28,[x1,#336]
add x29,x7,x8
ld64b x10,[x29]
adds x2,x20,x14
adcs x9,x21,x15
adcs x30,x22,x16
adcs x18,x23,x17
cset x6,hs
AC4F x20,x21,x22,x23,x10,x11,x12,x13,x5
W1F x28
ldp x8,x28,[x1,#352]
W2F
ldp x8,x28,[x1,#368]
W2F
ldr x8,[x1,#384]
W1F x8
LF x20,x21,x22,x23,x5,x29
LF x2,x9,x30,x18,x6,x29
.endm
.macro PFWD2
EFWD2
st64b x20,[x3]
add x3,x3,#1024
P1X
stp x2,x9,[x3]
stp x30,x18,[x3,#16]
stp x24,x25,[x3,#32]
stp x26,x27,[x3,#48]
add x3,x3,#1024
P1XB
.endm
// beta*w by six 43-bit windows (PAIRB6, proof.md 3.1): window 0 by LD64B into x8..x15; bit 128 picks slot 2 or 6
.macro BLK6 s0,s1,s2,s3,o0,o1,o2,o3
orr x20,x5,\s0,lsl #8
ld64b x8,[x20]
extr x20,\s1,\s0,#43
orr x20,x5,x20,lsl #8
ld64b x12,[x20]
adds \o0,x8,x16
adcs \o1,x9,x17
adcs \o2,x10,x18
adcs \o3,x11,x19
cset x1,hs
extr x20,\s2,\s1,#22
orr x20,x6,x20,lsl #8
ld64b x12,[x20]
adds \o0,\o0,x16
adcs \o1,\o1,x17
adcs \o2,\o2,x18
adcs \o3,\o3,x19
adc x1,x1,xzr
orr x20,x6,\s2,lsl #7
ld64b x12,[x20]
adds \o0,\o0,x12
adcs \o1,\o1,x13
adcs \o2,\o2,x14
adcs \o3,\o3,x15
adc x1,x1,xzr
extr x20,\s3,\s2,#44
orr x20,x7,x20,lsl #8
ld64b x12,[x20]
adds \o0,\o0,x12
adcs \o1,\o1,x13
adcs \o2,\o2,x14
adcs \o3,\o3,x15
adc x1,x1,xzr
lsr x20,\s3,#23
orr x20,x7,x20,lsl #8
ld64b x12,[x20]
adds \o0,\o0,x16
adcs \o1,\o1,x17
adcs \o2,\o2,x18
adcs \o3,\o3,x19
adc x1,x1,xzr
mul x20,x1,x4
adds \o0,\o0,x20
b.lo 9f
adcs \o1,\o1,xzr
adcs \o2,\o2,xzr
adcs \o3,\o3,xzr
csel x20,x4,xzr,hs
add \o0,\o0,x20
stadd x4,[sp]
9:
.endm
.macro PAIRB
mov x5,#0xfff8000000000000
orr x6,x5,#64
orr x7,x5,#128
BLK6 x21,x22,x23,x24,x27,x28,x29,x30
BLK6 x25,x26,x0,x2,x8,x9,x10,x11
.endm
// limb m times row words r3, r2: f = lo(m r3) + hi(m r2) (TS0), f += that (TSK); 64-bit words (proof.md 3.3)
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
SHX x5,x25,x28,x26,x27
RA4 x21,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
RA4 x17,x5,x25,x26,x27,x0,x1,x2,x24,x29,x30,xzr,x6,x7,x8,x28
SHX x6,x5,x27,x25,x26
RA4 x22,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
RA4 x18,x6,x5,x25,x26,x0,x1,x2,x24,x29,x30,x30,x7,x8,x28,x27
SHX x7,x6,x26,x5,x25
RA4 x23,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
RA4 x19,x7,x6,x5,x25,x0,x1,x2,x24,x29,x30,x30,x8,x28,x27,x26
FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12
FDX x0,x1,x29,x30,x0,x1,x6,x29,x2,x24
.endm
// DT (C1BS's tail, proof.md 3.1): 1/D_(j+1) = 1/D_j (x0, x1, x2, x5) + delta inv_e (x9..x12), delta (x13) = A_(j+1) - A_j
// in (0, 2^63); top word c (x6) < 2^63 + 2; DF folds c C0 < 2^96 (two words), branching on limb 1's carry (3.4)
.macro DF
mul x7,x6,x4
umulh x8,x6,x4
adds x0,x0,x7
adcs x1,x1,x8
b.lo 9f
adcs x2,x2,xzr
adcs x5,x5,xzr
csel x7,x4,xzr,hs
adds x0,x0,x7
adc x1,x1,xzr
stadd x4,[sp]
9:
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
DF
.endm
.macro C1BS
ldp x13,x14,[sp,#16]
ldp x15,x16,[sp,#32]
ldp x17,x18,[x3,#192]
ldp x19,x20,[x3,#208]
FS x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x21
RF x9,x13,x14,x15,x16,x0,x1,x2,x5,x6,x7,x8
SHX x17,x13,x16,x14,x15
RA4 x10,x17,x13,x14,x15,x0,x1,x2,x5,x6,x18,xzr,x7,x8,x19,x20
SHX x16,x17,x15,x13,x14
RA4 x11,x16,x17,x13,x14,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
SHX x15,x16,x14,x17,x13
RA4 x12,x15,x16,x17,x13,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
FDX x0,x1,x6,x18,x0,x1,x7,x8,x2,x5
stp x0,x1,[x3,#-480]
stp x2,x5,[x3,#-464]
ldr x13,[x3,#-128]
DT
.endm
.macro P1Q b,w
bfi \b,x26,#9,#\w
ld64b x10,[\b]
AC4 x0,x1,x2,x5,x10,x11,x12,x13,x18
AC4 x20,x21,x22,x23,x14,x15,x16,x17,x19
.endm
// P1 (proof.md 3.1): w = inv^2; M = inv + 2A w (x0, x1, x2, x5), L = 2y_X w (x20..x23); N = y_Y L (rows by LD64B: YLS);
// u^2 = M - N (x21..x24), v^2 = M + N (x25, x26, x0, x2); P1_j's line bases at entry +352, +416; W7 windows of w
.macro P1
SQ13 x0,x1,x2,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14
F5 x6,x7,x8,x9,x13,x6,x7,x14,x13
ldp x24,x25,[x3,#352]
bfi x24,x6,#9,#37
ld64b x16,[x24]
adds x0,x0,x16
adcs x1,x1,x17
adcs x2,x2,x18
adcs x5,x5,x19
cset x18,hs
extr x26,x7,x6,#37
bfi x25,x26,#9,#37
ld64b x10,[x25]
AC4 x0,x1,x2,x5,x10,x11,x12,x13,x18
AC4F x20,x21,x22,x23,x14,x15,x16,x17,x19
ldp x24,x25,[x3,#368]
lsr x26,x7,#10
P1Q x24,37
extr x26,x8,x7,#47
P1Q x25,37
ldp x24,x25,[x3,#416]
lsr x26,x8,#20
P1Q x24,36
extr x26,x9,x8,#56
P1Q x25,36
ldr x24,[x3,#432]
lsr x26,x9,#28
P1Q x24,36
LF x0,x1,x2,x5,x18,x6
LF x20,x21,x22,x23,x19,x6
add x28,sp,#128
ld64b x6,[x28]
RF x20,x6,x7,x8,x9,x14,x15,x16,x17,x18,x24,x25
RA4 x21,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,xzr,x24,x25,x26,x27
add x28,sp,#192
ld64b x6,[x28]
RA4 x22,x6,x7,x8,x9,x14,x15,x16,x17,x18,x19,x19,x24,x25,x26,x27
RA4 x23,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x19,x24,x25,x26,x27
FDX x14,x15,x18,x19,x14,x15,x24,x25,x16,x17
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
// CHAIN (3.3): c = (x5, x14, x6, x16); [l2, l3, a0, a1] + [h1, h2, h3, 0] = (u, t, s, v) (x1, x20, x5, x7); rows' words 2, 3:
// (a2, a3), (a1, a2), (s, v), (t, s), (u, t); H words x17, x13, x19, x15 (R_0..R_3), L words x16, x12, x18, x14 (R_1..R_4)
.macro TFP4S
umulh x7,x14,x4
umulh x20,x6,x4
mul x1,x6,x4
umulh x12,x16,x4
mul x13,x16,x4
adds x1,x1,x7
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
// corrections -T(y)'s words: H += hi, L += lo
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
// PX bwd: line [X_j's | X_(j+1)'s] of PE/PF/PG serves both incidences (x6..x9, x14; x10..x13, x15)
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
.macro LK2 oa,ob,oc,od
ldp x0,x2,[sp,#\oa]
add x0,x24,x0
ld64b x6,[x0]
add x2,x24,x2
ld64b x16,[x2]
AC4F x6,x7,x8,x9,x16,x17,x18,x19,x14
AC4F x10,x11,x12,x13,x20,x21,x22,x23,x15
ldp x0,x2,[sp,#\ob]
W2P
ldp x0,x2,[sp,#\oc]
W2P
ldr x0,[sp,#\od]
W1P x0
.endm
// the line [c_j | c_(j+1)] folded, stored at x3 - \d
.macro LK2S d
LF x6,x7,x8,x9,x14,x0
LF x10,x11,x12,x13,x15,x0
sub x29,x3,#\d
st64b x6,[x29]
.endm
.macro INC ea,eb,fa,fb,ga,gb
P1
PAIRB
ldp x5,x14,[x3,#\ea]
ldp x6,x16,[x3,#\eb]
TFP4S
CORS sp,112,256
add x1,x3,#128
st64b x12,[x1]
ldp x5,x14,[x3,#\fa]
ldp x6,x16,[x3,#\fb]
TFP4S
CORS x3,320,336
stp x16,x17,[x3,#288]
stp x18,x19,[x3,#304]
stp x12,x13,[x3,#96]
stp x14,x15,[x3,#112]
ldp x5,x14,[x3,#\ga]
ldp x6,x16,[x3,#\gb]
TFP4S
stp x16,x17,[x3,#256]
stp x18,x19,[x3,#272]
stp x12,x13,[x3,#64]
stp x14,x15,[x3,#80]
.endm
.macro BWDP
ldr x24,[x3,#224]
LK2 272,288,304,320
ldp x16,x17,[sp,#80]
ldp x18,x19,[sp,#96]
AC4 x6,x7,x8,x9,x16,x17,x18,x19,x14
AC4 x10,x11,x12,x13,x16,x17,x18,x19,x15
LK2S 64
LK2 336,352,368,384
LK2S 448
LK2 400,416,432,448
LK2S 256
ld64b x16,[x3]
W2S
stp x0,x1,[x3,#-496]
stp x2,x24,[x3,#-512]
C1BS
INC -32,-16,-416,-400,-224,-208
sub x3,x3,#512
ldp x0,x1,[x3,#32]
ldp x2,x5,[x3,#48]
INC 448,464,64,80,256,272
ldp x25,x26,[x3,#16]
ldp x27,x28,[x3],#-512
.endm
// After slot k (x16; x15 = 8*d1): LC to slot 1, RC cursor (2^43 + 8*d1) += 32k (x0); the header sends the block iff the
// old cursor is below the slot end x12; s, r: two leaf registers, dead once the record is stored.
.macro FLUSH s,r
sub \s,x16,#4080
sub \s,\s,#4080
str \s,[x15]
orr \r,x15,#0x80000000000
ldadd x0,x16,[\r]
cmp x16,x12
b.hs 6f
str x0,[\s,#-16]
stp x15,x16,[\s,#-32]
6:
.endm
// j = 0 (proof.md 3.2): sigma = W0's H + W1's H; K2 = W0's H, bits 0..21 from W1's H bits 22..43 (BFXIL), K3 = W1's H (emit 7)
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
// j = 1 (JHB: W1, N') and j = 2 (JTB: N', W0), N' = (~sigma, ~S), S = W0's L + W1's L (x14); then the sign's later leaves (C2, C3)
.macro JHB k,w0l,w0h,w1l,w1h
kh\k:
add x14,\w0l,\w1l
cmn \w1h,x15
b.hs kt\k
mvn \w0h,x15
and x15,x1,\w1l,lsr #25
ldadd x4,x16,[x15]
eon \w0l,x18,x14,lsr #28
bfxil \w1h,\w0h,#22,#22
stp \w1h,\w0h,[x16]
stp x17,\w0l,[x16,#16]
tbz x16,#13,km\k
FLUSH \w0h,\w1h
km\k:
.endm
.macro JTB k,w0l,w0h,w1l,w1h
kv\k:
add x14,\w0l,\w1l
kt\k:
mvn \w1h,x15
bic x15,x1,x14,lsr #25
ldadd x4,x16,[x15]
add \w0l,x18,\w0l,lsr #28
bfxil \w1h,\w0h,#22,#22
stp \w1h,\w0h,[x16]
stp x17,\w0l,[x16,#16]
tbz x16,#13,ku\k
FLUSH \w0h,\w1h
ku\k:
.endm
// gamma, beta'' line by LD64B, alpha' by LDP; x3 ends at entry+128 (sign +) or entry-256 (sign -); x18 = (entry+128)<<30 once
// per incidence (signless w3; JH's w3 by EON with it)
.macro TLDP
ld64b x20,[x3]
ldp x19,x28,[x3,#-96]
ldp x29,x2,[x3,#-80]
sub x3,x3,#128
lsl x18,x3,#30
.endm
.macro TLDM
ldp x19,x28,[x3],#-64
ldp x29,x2,[x3,#80]
ld64b x20,[x3]
sub x3,x3,#320
.endm
// W_m(t)=U_m(alpha)+U_(m+t)(beta')+U_(m+2t)(gamma) word by word, U_2=-U_0-U_1
.macro S1
add x5,x19,x24
add x6,x28,x25
add x5,x5,x20
add x6,x6,x21
add x7,x29,x26
add x8,x2,x27
add x7,x7,x22
add x8,x8,x23
.endm
// W(1) = (i3+i6, i5-i1), W(2) = (i3-i2, i5+i6)
.macro S2
sub x9,x24,x20
sub x10,x25,x21
sub x11,x26,x22
sub x13,x27,x23
sub x19,x19,x24
sub x28,x28,x25
sub x29,x29,x26
sub x2,x2,x27
add x24,x9,x11
add x25,x10,x13
add x26,x19,x24
add x27,x28,x25
sub x20,x29,x9
sub x21,x2,x10
.endm
.macro S3
sub x19,x19,x11
sub x28,x28,x13
add x29,x29,x24
add x2,x2,x25
.endm
// TJX's JH, JT blocks run copies of the sign's later leaves (C2, C3), ending in one B to kne\k (the third leaf's join).
.macro TSUM k
S1
LJ0 \k,x5,x6,x7,x8
S2
LJ0 d\k,x26,x27,x20,x21
S3
LJ0 e\k,x19,x28,x29,x2
.endm
.macro C3 k,j
S3
LJ0 \k,x19,x28,x29,x2
b \j
JHB \k,x19,x28,x29,x2
b \j
JTB \k,x19,x28,x29,x2
b \j
.endm
.macro C2 k,j
S2
LJ0 \k,x26,x27,x20,x21
C3 a\k,\j
JHB \k,x26,x27,x20,x21
C3 b\k,\j
JTB \k,x26,x27,x20,x21
C3 c\k,\j
.endm
.macro TJX k
JHB e\k,x19,x28,x29,x2
b kne\k
JTB e\k,x19,x28,x29,x2
b kne\k
JHB d\k,x26,x27,x20,x21
C3 g\k,kne\k
JTB d\k,x26,x27,x20,x21
C3 i\k,kne\k
JHB \k,x5,x6,x7,x8
C2 l\k,kne\k
JTB \k,x5,x6,x7,x8
C2 o\k,kne\k
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
add x1,x1,#512
2:
subs x0,x0,#32
b.ne 1b
ldp x0,x30,[sp],#16
ret
.globl bwd
bwd:
mov x4,#977
movk x4,#1,lsl #32
str x30,[sp,#8]
1:
BWDP
BWDP
BWDP
BWDP
BWDP
BWDP
BWDP
BWDP
add x5,x3,#512
tst x5,#0x7fe00
b.ne 1b
tbnz x3,#34,2f
ldp x5,x6,[sp]
sub sp,sp,#512
stp x5,x6,[sp]
b 1b
2:
ldp x0,x30,[sp]
ret
.globl keys
keys:
mov x1,#0x7ffffffff8
mov x4,#32
mov x0,#8192
1:
TLDP
TSUM p00
TLDM
TSUM m00
TLDP
TSUM p10
TLDM
TSUM m10
tbz x3,#34,1b
ret
TJX p00
TJX m00
TJX p10
TJX m10

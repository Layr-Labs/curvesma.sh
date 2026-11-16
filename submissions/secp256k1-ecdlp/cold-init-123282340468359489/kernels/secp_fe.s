.macro SQK off
add x1,sp,#\off
bl	fe_sqr_f
.endm
.macro SQK2 o
SQK \o
SQK \o
.endm
.macro SQK10 o
SQK2 \o
SQK2 \o
SQK2 \o
SQK2 \o
SQK2 \o
.endm
.macro SQK21 o
SQK10 \o
SQK10 \o
SQK \o
.endm
.macro SQK43 o
SQK21 \o
SQK21 \o
SQK \o
.endm
.macro REDC aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av
mul \aa,\ab,\ac
umulh \ab,\ab,\ac
mul x1,\ad,\ac
adds \ae,\aa,\af
umulh \ad,\ad,\ac
mul \ag,\ah,\ac
adcs \ab,x1,\ab
umulh \ah,\ah,\ac
adcs \ai,\ag,\ad
cinc \ah,\ah,hs
adds \aj,\ak,\al
adcs \am,\an,\ae
adcs \ao,x15,\ab
cset \ap,hs
adcs x15,\aq,\ai
cinc \ah,\ah,hs
mul \ar,\ah,\ac
umulh \ah,\ah,\ac
adds \aj,\aj,\ar
adcs x9,\am,\ah
cset \as,hs
adcs xzr,\ao,xzr
adcs xzr,x15,xzr
csel \at,\ac,xzr,hs
adds x8,\aj,\at
adcs x9,x9,xzr
stp x8,x9,[x0]
add x8,\aq,\ai
adcs x9,\ao,\au
adc x8,x8,\av
stp x9,x8,[x0,#16]
.endm
fe_mul_f:
str x30,[sp,#-96]!
ldp x8,x13,[x1]
ldp x9,x10,[x2]
ldp x12,x11,[x2,#16]
ldp x22,x1,[x1,#16]
umulh x14,x8,x9
mul x15,x8,x10
umulh x16,x8,x10
mul x17,x8,x12
adds x14,x15,x14
umulh x18,x8,x12
mul x2,x8,x11
adcs x15,x17,x16
umulh x3,x8,x11
umulh x5,x13,x9
adcs x17,x2,x18
mul x6,x13,x10
cinc x18,x3,hs
umulh x7,x13,x10
mul x19,x13,x12
adds x2,x6,x5
umulh x20,x13,x12
mul x21,x13,x11
adcs x5,x19,x7
mul x4,x13,x9
umulh x13,x13,x11
adcs x6,x21,x20
umulh x24,x22,x9
mul x25,x22,x10
cinc x13,x13,hs
adds x14,x14,x4
umulh x26,x22,x10
adcs x15,x15,x2
adcs x17,x17,x5
mul x27,x22,x12
adcs x18,x18,x6
cinc x13,x13,hs
adds x4,x25,x24
umulh x28,x22,x12
mul x30,x22,x11
adcs x5,x27,x26
mul x23,x22,x9
umulh x22,x22,x11
adcs x6,x30,x28
umulh x3,x1,x9
mul x7,x1,x10
cinc x19,x22,hs
adds x15,x15,x23
umulh x10,x1,x10
adcs x17,x17,x4
adcs x18,x18,x5
mul x2,x1,x12
adcs x13,x13,x6
umulh x12,x1,x12
mul x20,x1,x11
mul x16,x1,x9
umulh x11,x1,x11
cinc x1,x19,hs
adds x3,x7,x3
adcs x10,x2,x10
adcs x12,x20,x12
mul x8,x8,x9
cinc x9,x11,hs
adds x11,x17,x16
mov x17,#977
adcs x16,x18,x3
movk x17,#1,lsl #32
adcs x10,x13,x10
mul x18,x16,x17
adcs x12,x1,x12
umulh x16,x16,x17
cinc x9,x9,hs
REDC x13,x10,x17,x12,x13,x16,x2,x9,x12,x8,x8,x18,x13,x14,x10,w14,x11,x16,w13,x15,x13,x14
ldr x30,[sp],#96
ret
fe_sqr_f:
ldp x9,x10,[x1]
mov x8,#977
ldp x13,x14,[x1,#16]
movk x8,#1,lsl #32
umulh x12,x10,x9
mul x15,x13,x9
umulh x16,x13,x9
mul x17,x14,x9
adds x12,x12,x15
umulh x18,x14,x9
umulh x2,x13,x10
adcs x15,x16,x17
mul x3,x14,x10
cinc x16,x18,hs
mul x1,x13,x10
umulh x4,x14,x10
adds x17,x2,x3
mul x11,x10,x9
mul x5,x14,x13
cinc x18,x4,hs
adds x15,x15,x1
umulh x6,x14,x13
adcs x16,x16,x17
extr x2,x12,x11,#63
extr x12,x15,x12,#63
mul x7,x9,x9
adcs x17,x18,x5
umulh x9,x9,x9
cinc x18,x6,hs
mul x19,x10,x10
umulh x10,x10,x10
adds x9,x9,x11,lsl #1
extr x11,x16,x15,#63
extr x16,x17,x16,#63
mul x1,x13,x13
extr x17,x18,x17,#63
lsr x18,x18,#63
adcs x15,x2,x19
umulh x13,x13,x13
adcs x10,x12,x10
mul x3,x14,x14
adcs x11,x11,x1
umulh x14,x14,x14
adcs x13,x16,x13
mul x12,x11,x8
adcs x17,x17,x3
umulh x11,x11,x8
adc x14,x18,x14
REDC x16,x13,x8,x17,x11,x11,x18,x14,x16,x12,x7,x12,x9,x9,x11,w13,x10,x17,w14,x8,x14,x13
ret
.macro QA d,s
add x0,sp,#\d
add x1,sp,#\s
bl	fe_sqr_f
.endm
.macro MA s,t
add x1,sp,#\s
add x2,sp,#\t
bl	fe_mul_f
.endm
.macro Q9 o
sub x1,x29,#\o
bl	fe_sqr_f
.endm
fe_inv:
sub sp,sp,#432
stp x29,x30,[sp,#384]
add x29,sp,#384
stp x1,x0,[sp,#400]
sub x0,x29,#32
bl	fe_sqr_f
sub x1,x29,#32
MULSP400
sub x0,x29,#64
Q9 32
sub x1,x29,#64
MULSP400
sub x0,x29,#96
Q9 64
Q9 96
Q9 96
sub x1,x29,#96
sub x2,x29,#64
bl	fe_mul_f
sub x0,x29,#128
Q9 96
Q9 128
Q9 128
sub x1,x29,#128
sub x2,x29,#64
bl	fe_mul_f
sub x0,x29,#160
Q9 128
Q9 160
sub x1,x29,#160
sub x2,x29,#32
bl	fe_mul_f
add x0,sp,#192
Q9 160
SQK10 192
add x1,sp,#192
sub x2,x29,#160
bl	fe_mul_f
QA 160,192
SQK21 160
MA 160,192
QA 128,160
SQK43 128
MA 128,160
QA 96,128
SQK43 96
SQK43 96
SQK 96
MA 96,128
QA 64,96
SQK43 64
MA 64,160
QA 32,64
SQK2 32
add x1,sp,#32
sub x2,x29,#64
bl	fe_mul_f
mov x0,sp
add x1,sp,#32
bl	fe_sqr_f
SQK21 0
SQK 0
mov x1,sp
add x2,sp,#192
bl	fe_mul_f
SQK2 0
SQK2 0
SQK 0
mov x1,sp
MULSP400
SQK2 0
SQK 0
mov x1,sp
sub x2,x29,#32
bl	fe_mul_f
SQK2 0
mov x1,sp
ldr x0,[sp,#408]
MULSP400
ldp x29,x30,[sp,#384]
add sp,sp,#432
ret
.macro FSUB d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3
subs \d0,\a0,\b0
sbcs \d1,\a1,\b1
sbcs \d2,\a2,\b2
sbcs \d3,\a3,\b3
csel x5,x4,xzr,lo
subs \d0,\d0,x5
sbcs \d1,\d1,xzr
sbcs \d2,\d2,xzr
sbc \d3,\d3,xzr
.endm
.macro FADD d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3
adds \d0,\a0,\b0
adcs \d1,\a1,\b1
adcs \d2,\a2,\b2
adcs \d3,\a3,\b3
csel x5,x4,xzr,hs
adds \d0,\d0,x5
adcs \d1,\d1,xzr
adcs \d2,\d2,xzr
adc \d3,\d3,xzr
.endm
.macro RED v0,v1,v2,v3
adds x22,\v0,x4
adcs x23,\v1,xzr
adcs x24,\v2,xzr
adcs x25,\v3,xzr
csel \v0,x22,\v0,hs
csel \v1,x23,\v1,hs
csel \v2,x24,\v2,hs
csel \v3,x25,\v3,hs
.endm
.macro SHXC n0,r0,r1,r2,r3
mul \n0,\r3,x4
umulh \r3,\r3,x4
adds \r0,\r0,\r3
adcs \r1,\r1,xzr
adcs \r2,\r2,xzr
csel x5,x4,xzr,hs
adds \n0,\n0,x5
adcs \r0,\r0,xzr
adcs \r1,\r1,xzr
adc \r2,\r2,xzr
RED \n0,\r0,\r1,\r2
.endm
.macro NEGT
mul x16,x12,x4
umulh x17,x12,x4
mul x18,x13,x4
umulh x19,x13,x4
cmn x16,x4
adcs x17,x17,x18
adc x19,x19,xzr
cmn x10,x19
adcs xzr,x11,xzr
adcs x16,x12,xzr
adc x17,x13,xzr
sub x14,xzr,x16
sub x15,xzr,x17
.endm
.macro L4 r0,r1,r2,r3,b,o,q
ldp \r0,\r1,[\b,#\o]
ldp \r2,\r3,[\b,#\q]
.endm
.macro S4 r0,r1,r2,r3,b,o,q
stp \r0,\r1,[\b,#\o]
stp \r2,\r3,[\b,#\q]
.endm
.macro OFF5 w0,w1,w2,w3,t0,t1,t2,sh,s0,s1,s2,b0,b1,b2,b3,b4
mov \t0,#\b0
bfi \t0,\w0,#\sh,#43
extr \t2,\w1,\w0,#43
mov \t1,#\b1
bfi \t1,\t2,#\sh,#43
stp \t0,\t1,[x8,#\s0]
extr \t0,\w2,\w1,#22
mov \t1,#\b2
bfi \t1,\t0,#\sh,#43
extr \t0,\w3,\w2,#1
mov \t2,#\b3
bfi \t2,\t0,#\sh,#43
stp \t1,\t2,[x8,#\s1]
extr \t0,\w3,\w2,#44
mov \t1,#\b4
bfi \t1,\t0,#\sh,#43
str \t1,[x8,#\s2]
.endm
.macro OFFS7 w0,w1,w2,w3,s0,s1,s2,t0,t1,t2,b0,b1,b2,b3,b4,b5
mov \t0,#\b0
bfi \t0,\w0,#8,#43
extr \t2,\w1,\w0,#43
mov \t1,#\b1
bfi \t1,\t2,#8,#43
stp \t0,\t1,[x8,#\s0]
lsr \t0,\w1,#22
mov \t1,#\b2
bfi \t1,\t0,#8,#42
mov \t2,#\b3
bfi \t2,\w2,#8,#43
stp \t1,\t2,[x8,#\s1]
extr \t0,\w3,\w2,#43
mov \t1,#\b4
bfi \t1,\t0,#8,#43
lsr \t0,\w3,#22
mov \t2,#\b5
bfi \t2,\t0,#8,#42
stp \t1,\t2,[x8,#\s2]
.endm
.macro WORDS bb,bo,bq
add x0,sp,#256
add x1,\bb,#\bo
bl fe_sqr_f
ldr x9,[sp,#1016]
cbz x9,7f
add x0,sp,#512
add x1,sp,#256
bl fe_sqr_f
ldr x0,[x29,#2240]
add x0,x0,#368
add x1,sp,#512
add x2,sp,#1024
bl eisr
7:
add x0,sp,#288
add x1,sp,#256
add x2,\bb,#\bo
bl fe_mul_f
add x0,sp,#320
add x1,\bb,#\bo
add x2,sp,#352
bl fe_mul_f
ldr x4,[sp,#384]
ldr x8,[x29,#2240]
L4 x10,x11,x12,x13,x8,128,144
RED x10,x11,x12,x13
S4 x10,x11,x12,x13,x8,128,144
SHXC x14,x10,x11,x12,x13
S4 x14,x10,x11,x12,x8,160,176
SHXC x13,x14,x10,x11,x12
S4 x13,x14,x10,x11,x8,192,208
SHXC x12,x13,x14,x10,x11
S4 x12,x13,x14,x10,x8,224,240
L4 x10,x11,x12,x13,\bb,\bo,\bq
NEGT
stp x14,x15,[x8,#320]
OFFS7 x10,x11,x12,x13,256,272,448,x14,x15,x16,0xff80000000000000,0xff88000000000000,0xffa0000000000000,0xff90000000000000,0xff98000000000000,0xffa4000000000000
L4 x10,x11,x12,x13,sp,288,304
RED x10,x11,x12,x13
OFF5 x10,x11,x12,x13,x14,x15,x16,6,96,112,352,0,0x2000000000000,0x4000000000000,0x6000000000000,0x8000000000000
OFF5 x10,x11,x12,x13,x14,x15,x16,8,408,424,440,0xffa8000000000000,0xffb0000000000000,0xffb8000000000000,0xffc0000000000000,0xffc8000000000000
lsr x14,x13,#23
str x14,[x8,#56]
mov x5,#7
adds x10,x10,x5
adcs x11,x11,xzr
adcs x12,x12,xzr
adc x13,x13,xzr
RED x10,x11,x12,x13
S4 x10,x11,x12,x13,x8,288,304
add x0,x8,#64
add x1,x8,#288
bl fe_sqr_f
ldr x4,[sp,#384]
ldr x8,[x29,#2240]
L4 x10,x11,x12,x13,x8,64,80
RED x10,x11,x12,x13
S4 x10,x11,x12,x13,x8,64,80
L4 x10,x11,x12,x13,sp,256,272
RED x10,x11,x12,x13
FADD x14,x15,x16,x17,x10,x11,x12,x13,x10,x11,x12,x13
FADD x14,x15,x16,x17,x14,x15,x16,x17,x10,x11,x12,x13
RED x14,x15,x16,x17
OFFS7 x14,x15,x16,x17,464,480,496,x10,x18,x19,0xffd0000000000000,0xffd8000000000000,0xfff0000000000000,0xffe0000000000000,0xffe8000000000000,0xfff4000000000000
L4 x10,x11,x12,x13,sp,320,336
RED x10,x11,x12,x13
NEGT
stp x14,x15,[x8,#336]
add x8,x8,#512
str x8,[x29,#2240]
.endm
.macro POINT op,a0,a1,a2,a3,b0,b1,b2,b3
ldr x4,[sp,#384]
ldr x6,[sp,#432]
L4 x10,x11,x12,x13,x6,32,48
L4 x14,x15,x16,x17,x29,32,48
\op x18,x19,x20,x21,x10,x11,x12,x13,x14,x15,x16,x17
S4 x18,x19,x20,x21,sp,64,80
add x0,sp,#96
add x1,sp,#64
add x2,sp,#32
bl fe_mul_f
add x0,sp,#128
add x1,sp,#96
bl fe_sqr_f
ldr x4,[sp,#384]
L4 x10,x11,x12,x13,sp,128,144
L4 x14,x15,x16,x17,sp,224,240
FSUB x18,x19,x20,x21,x10,x11,x12,x13,x14,x15,x16,x17
RED x18,x19,x20,x21
S4 x18,x19,x20,x21,sp,448,464
L4 x14,x15,x16,x17,x29,0,16
FSUB x10,x11,x12,x13,\a0,\a1,\a2,\a3,\b0,\b1,\b2,\b3
S4 x10,x11,x12,x13,sp,160,176
add x0,sp,#192
add x1,sp,#96
add x2,sp,#160
bl fe_mul_f
ldr x4,[sp,#384]
L4 x10,x11,x12,x13,sp,192,208
L4 x14,x15,x16,x17,x29,32,48
FSUB x18,x19,x20,x21,x10,x11,x12,x13,x14,x15,x16,x17
RED x18,x19,x20,x21
.endm
.macro YOUT
ldr x8,[x29,#2240]
S4 x18,x19,x20,x21,x8,128,144
WORDS sp,448,464
.endm
.macro SUMS
ldr x4,[sp,#384]
ldr x6,[sp,#432]
L4 x10,x11,x12,x13,x29,0,16
L4 x14,x15,x16,x17,x6,0,16
FADD x18,x19,x20,x21,x10,x11,x12,x13,x14,x15,x16,x17
RED x18,x19,x20,x21
S4 x18,x19,x20,x21,sp,224,240
.endm
.macro DEN
ldr x4,[sp,#384]
L4 x10,x11,x12,x13,x1,0,16
L4 x14,x15,x16,x17,x29,0,16
FSUB x18,x19,x20,x21,x10,x11,x12,x13,x14,x15,x16,x17
S4 x18,x19,x20,x21,x0,32,48
add x2,x0,#32
sub x1,x0,x6
bl fe_mul_f
.endm
yfwd_chain:
str x30,[sp,#392]
add x0,x29,#1152
add x1,x29,#1088
mov x6,#1280
DEN
add x0,x0,#64
1:
sub x1,x0,#1152
mov x6,#64
DEN
add x0,x0,#64
add x9,x29,#2240
cmp x0,x9
b.ne 1b
ldr x30,[sp,#392]
ret
ybwd_chain:
str x30,[sp,#392]
add x9,x29,#2176
str x9,[sp,#408]
add x9,x29,#1024
str x9,[sp,#432]
1:
ldr x9,[sp,#408]
add x0,sp,#32
mov x1,sp
sub x2,x9,#64
bl fe_mul_f
ldr x9,[sp,#408]
mov x0,sp
mov x1,sp
add x2,x9,#32
bl fe_mul_f
SUMS
POINT FSUB,x14,x15,x16,x17,x18,x19,x20,x21
YOUT
POINT FADD,x18,x19,x20,x21,x14,x15,x16,x17
YOUT
ldr x9,[sp,#408]
sub x9,x9,#64
str x9,[sp,#408]
ldr x6,[sp,#432]
sub x6,x6,#64
str x6,[sp,#432]
cmp x6,x29
b.ne 1b
add x0,sp,#32
mov x1,sp
sub x2,x29,#128
bl fe_mul_f
mov x0,sp
mov x1,sp
add x2,x29,#1184
bl fe_mul_f
ldr x8,[x29,#2240]
L4 x10,x11,x12,x13,x29,32,48
S4 x10,x11,x12,x13,x8,128,144
WORDS x29,0,16
add x6,x29,#1088
str x6,[sp,#432]
SUMS
POINT FSUB,x14,x15,x16,x17,x18,x19,x20,x21
S4 x18,x19,x20,x21,x29,32,48
L4 x10,x11,x12,x13,sp,448,464
S4 x10,x11,x12,x13,x29,0,16
ldr x30,[sp,#392]
ret
.globl yround
yround:
str x30,[sp,#400]
ldr x29,[sp,#424]
1:
bl yfwd_chain
add x29,x29,#2304
ldr x9,[sp,#416]
cmp x29,x9
b.ne 1b
sub x1,x29,#128
mov x0,sp
bl fe_inv
1:
sub x29,x29,#2304
bl ybwd_chain
ldr x9,[sp,#424]
cmp x29,x9
b.ne 1b
ldr x30,[sp,#400]
ret
.macro MULLO3 r0,r1,r2,a0,a1,a2,b0,b1,b2,t
mul \r0,\a0,\b0
umulh \r1,\a0,\b0
mul \t,\a0,\b1
adds \r1,\r1,\t
umulh \r2,\a0,\b1
adc \r2,\r2,xzr
madd \r2,\a0,\b2,\r2
mul \t,\a1,\b0
adds \r1,\r1,\t
umulh \t,\a1,\b0
adc \r2,\r2,\t
madd \r2,\a1,\b1,\r2
madd \r2,\a2,\b0,\r2
.endm
.macro FRAC f0,f1,f2,f3,r0,r1,r2,r3
umulh x11,\r0,x4
mul x12,\r1,x4
umulh x13,\r1,x4
mul x14,\r2,x4
umulh x15,\r2,x4
mul x16,\r3,x4
umulh x17,\r3,x4
adds x11,x11,x12
adcs x13,x13,x14
adcs x15,x15,x16
adc x17,x17,xzr
cmn x15,x9
adc x17,x17,xzr
adds \f0,\r0,x17
adcs \f1,\r1,xzr
adcs \f2,\r2,xzr
adc \f3,\r3,xzr
.endm
.globl eisr
eisr:
sub sp,sp,#160
stp x29,x30,[sp,#144]
stp x0,x1,[sp,#128]
str x2,[sp,#120]
mov x0,sp
bl fe_mul_f
ldr x1,[sp,#136]
ldr x2,[sp,#120]
add x2,x2,#32
add x0,sp,#32
bl fe_mul_f
ldr x2,[sp,#120]
PC977
ldp x10,x11,[sp]
ldp x12,x13,[sp,#16]
ldp x14,x15,[x2,#64]
ldp x16,x17,[x2,#80]
FADD x10,x11,x12,x13,x10,x11,x12,x13,x14,x15,x16,x17
RED x10,x11,x12,x13
ldp x14,x15,[sp,#32]
ldp x16,x17,[sp,#48]
ldp x18,x19,[x2,#96]
ldp x20,x21,[x2,#112]
FADD x14,x15,x16,x17,x14,x15,x16,x17,x18,x19,x20,x21
RED x14,x15,x16,x17
ldp x18,x19,[x2,#200]
ldp x20,x21,[x2,#216]
mul x6,x10,x18
umulh x7,x10,x18
madd x7,x10,x19,x7
madd x7,x11,x18,x7
mul x8,x14,x20
umulh x9,x14,x20
madd x9,x14,x21,x9
madd x9,x15,x20,x9
adds x6,x6,x8
adc x7,x7,x9
ldp x18,x19,[x2,#128]
mul x8,x6,x18
umulh x9,x6,x18
madd x9,x6,x19,x9
madd x9,x7,x18,x9
ldp x18,x19,[x2,#256]
adds x26,x8,x18
adc x27,x9,x19
ldp x18,x19,[x2,#152]
ldr x20,[x2,#168]
MULLO3 x6,x7,x8,x10,x11,x12,x18,x19,x20,x9
ldp x18,x19,[x2,#176]
ldr x20,[x2,#192]
MULLO3 x21,x22,x23,x14,x15,x16,x18,x19,x20,x9
adds x6,x6,x21
adcs x7,x7,x22
adc x8,x8,x23
ldp x18,x19,[x2,#128]
ldr x20,[x2,#144]
MULLO3 x21,x22,x28,x6,x7,x8,x18,x19,x20,x9
ldp x18,x19,[x2,#232]
ldr x20,[x2,#248]
adds x23,x21,x18
adcs x24,x22,x19
adc x25,x28,x20
ldp x10,x11,[x2,#272]
ldr x12,[x2,#288]
ldp x13,x14,[x2,#296]
ldr x3,[x2,#312]
add x5,x2,#320
mov x0,x23
mov x1,x24
mov x2,x25
mov x6,x26
mov x7,x27
mov x9,#0
1:
ldp x15,x16,[x5]
ldp x17,x18,[x5,#16]
ldr x19,[x5,#32]
add x5,x5,#40
adds x15,x23,x15
adcs x16,x24,x16
adc x17,x25,x17
adds x18,x26,x18
adc x19,x27,x19
subs x20,x18,x13
sbcs x20,x19,x14
cset x21,lo
adds x20,x15,x10
adcs x20,x16,x11
adc x20,x17,x12
subs x22,x10,x15
sbcs x22,x11,x16
sbc x22,x12,x17
orr x20,x20,x22
lsr x20,x20,#63
bic x21,x21,x20
bic x21,x21,x9
orr x9,x9,x21
tst x21,#1
csel x0,x15,x0,ne
csel x1,x16,x1,ne
csel x2,x17,x2,ne
csel x6,x18,x6,ne
csel x7,x19,x7,ne
subs x3,x3,#1
b.ne 1b
ldr x8,[sp,#128]
stp x0,x1,[x8]
stp x2,x6,[x8,#16]
str x7,[x8,#32]
ldp x29,x30,[sp,#144]
add sp,sp,#160
ret
.macro EK a,b
ldp x0,x1,[sp,#16]
add x0,x0,#\a
add x1,x1,#\b
ldr x2,[sp,#40]
bl eisr
.endm
.macro MULSP400
ldr x2,[sp,#400]
bl	fe_mul_f
.endm
.macro PC977
mov x4,#977
movk x4,#1,lsl #32
.endm
.globl tbldk
tbldk:
stp x29,x30,[sp,#-64]!
stp x0,x1,[sp,#16]
stp x2,x3,[sp,#32]
1:
ldp x0,x1,[sp,#16]
add x1,x1,#32
ldr x2,[sp,#40]
bl eisr
EK 40,96
EK 80,160
EK 128,224
EK 168,288
EK 208,352
ldp x0,x1,[sp,#16]
str xzr,[x0,#120]
str xzr,[x0,#248]
add x0,x0,#256
str x0,[sp,#16]
PC977
mov x9,#6
2:
ldp x10,x11,[x1]
ldp x12,x13,[x1,#16]
ldp x14,x15,[x1,#32]
ldp x16,x17,[x1,#48]
FADD x14,x15,x16,x17,x14,x15,x16,x17,x10,x11,x12,x13
RED x14,x15,x16,x17
stp x14,x15,[x1,#32]
stp x16,x17,[x1,#48]
add x1,x1,#64
subs x9,x9,#1
b.ne 2b
ldr x2,[sp,#32]
subs x2,x2,#1
str x2,[sp,#32]
b.ne 1b
ldp x29,x30,[sp],#64
ret
.globl tbldf
tbldf:
PC977
mov x9,#0x8000000000000000
1:
mov x6,x0
mov x7,x1
mov x8,x3
2:
ldp x18,x19,[x7,#32]
ldp x20,x21,[x7,#48]
FRAC x22,x23,x24,x25,x18,x19,x20,x21
stp x22,x23,[x6]
stp x24,x25,[x6,#16]
add x6,x6,#32
ldp x10,x11,[x7]
ldp x12,x13,[x7,#16]
FADD x18,x19,x20,x21,x18,x19,x20,x21,x10,x11,x12,x13
RED x18,x19,x20,x21
stp x18,x19,[x7,#32]
stp x20,x21,[x7,#48]
add x7,x7,#64
subs x8,x8,#1
b.ne 2b
add x0,x0,x26
subs x2,x2,#1
b.ne 1b
ret

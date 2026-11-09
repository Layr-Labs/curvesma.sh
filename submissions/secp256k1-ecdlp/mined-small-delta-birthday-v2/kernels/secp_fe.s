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
mul x13,x10,x17
umulh x10,x10,x17
mul x1,x12,x17
adds x13,x13,x16
umulh x12,x12,x17
mul x2,x9,x17
adcs x10,x1,x10
umulh x9,x9,x17
adcs x12,x2,x12
cinc x9,x9,hs
adds x8,x8,x18
adcs x13,x14,x13
adcs x10,x15,x10
cset w14,hs
adcs x15,x11,x12
cinc x9,x9,hs
mul x16,x9,x17
umulh x9,x9,x17
adds x8,x8,x16
adcs x9,x13,x9
cset w13,hs
adcs xzr,x10,xzr
adcs xzr,x15,xzr
csel x15,x17,xzr,hs
adds x8,x8,x15
adcs x9,x9,xzr
stp x8,x9,[x0]
add x8,x11,x12
adcs x9,x10,x13
adc x8,x8,x14
stp x9,x8,[x0,#16]
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
mul x16,x13,x8
umulh x13,x13,x8
mul x1,x17,x8
adds x11,x16,x11
umulh x17,x17,x8
mul x18,x14,x8
adcs x13,x1,x13
umulh x14,x14,x8
adcs x16,x18,x17
cinc x14,x14,hs
adds x12,x7,x12
adcs x9,x9,x11
adcs x11,x15,x13
cset w13,hs
adcs x15,x10,x16
cinc x14,x14,hs
mul x17,x14,x8
umulh x14,x14,x8
adds x12,x12,x17
adcs x9,x9,x14
cset w14,hs
adcs xzr,x11,xzr
adcs xzr,x15,xzr
csel x8,x8,xzr,hs
adds x8,x12,x8
adcs x9,x9,xzr
stp x8,x9,[x0]
add x8,x10,x16
adcs x9,x11,x14
adc x8,x8,x13
stp x9,x8,[x0,#16]
ret
fe_inv:
// fe_mul_f/fe_sqr_f never write x0 (runs sharing a destination set it once: SQK); x1, x2, x19-x28 are scratch;
// the 87 squarings are unrolled; fe_inv spills (in, out) once at [sp,#400].
sub sp,sp,#432
stp x29,x30,[sp,#384]
add x29,sp,#384
stp x1,x0,[sp,#400]
sub x0,x29,#32
bl	fe_sqr_f
sub x1,x29,#32
ldr x2,[sp,#400]
bl	fe_mul_f
sub x0,x29,#64
sub x1,x29,#32
bl	fe_sqr_f
sub x1,x29,#64
ldr x2,[sp,#400]
bl	fe_mul_f
sub x0,x29,#96
sub x1,x29,#64
bl	fe_sqr_f
sub x1,x29,#96
bl	fe_sqr_f
sub x1,x29,#96
bl	fe_sqr_f
sub x1,x29,#96
sub x2,x29,#64
bl	fe_mul_f
sub x0,x29,#128
sub x1,x29,#96
bl	fe_sqr_f
sub x1,x29,#128
bl	fe_sqr_f
sub x1,x29,#128
bl	fe_sqr_f
sub x1,x29,#128
sub x2,x29,#64
bl	fe_mul_f
sub x0,x29,#160
sub x1,x29,#128
bl	fe_sqr_f
sub x1,x29,#160
bl	fe_sqr_f
sub x1,x29,#160
sub x2,x29,#32
bl	fe_mul_f
add x0,sp,#192
sub x1,x29,#160
bl	fe_sqr_f
SQK10 192
add x1,sp,#192
sub x2,x29,#160
bl	fe_mul_f
add x0,sp,#160
add x1,sp,#192
bl	fe_sqr_f
SQK21 160
add x1,sp,#160
add x2,sp,#192
bl	fe_mul_f
add x0,sp,#128
add x1,sp,#160
bl	fe_sqr_f
SQK43 128
add x1,sp,#128
add x2,sp,#160
bl	fe_mul_f
add x0,sp,#96
add x1,sp,#128
bl	fe_sqr_f
SQK43 96
SQK43 96
SQK 96
add x1,sp,#96
add x2,sp,#128
bl	fe_mul_f
add x0,sp,#64
add x1,sp,#96
bl	fe_sqr_f
SQK43 64
add x1,sp,#64
add x2,sp,#160
bl	fe_mul_f
add x0,sp,#32
add x1,sp,#64
bl	fe_sqr_f
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
ldr x2,[sp,#400]
bl	fe_mul_f
SQK2 0
SQK 0
mov x1,sp
sub x2,x29,#32
bl	fe_mul_f
SQK2 0
mov x1,sp
ldr x0,[sp,#408]
ldr x2,[sp,#400]
bl	fe_mul_f
ldp x29,x30,[sp,#384]
add sp,sp,#432
ret
// ygen: M*N points C+m*T+n*S and their 512-byte records. Scratch x29: 0 BETA, 32 F, 96 ACC, 128 INV, 160 LAM, 192 D,
// 224..384 T1..T6, 416.. state, 1024 Y array (64 B/pt), products. A = x10-13, B = x20-23, x28 = C0, x27 a temp.
.macro C0
ldr x28,[x29,#568]
.endm
.macro LA o,p,b
ldp x10,x11,[\b,#\o]
ldp x12,x13,[\b,#\p]
.endm
.macro LB o,p,b
ldp x20,x21,[\b,#\o]
ldp x22,x23,[\b,#\p]
.endm
.macro SA o,p,b
stp x10,x11,[\b,#\o]
stp x12,x13,[\b,#\p]
.endm
.macro SB o,p,b
stp x20,x21,[\b,#\o]
stp x22,x23,[\b,#\p]
.endm
.macro SAP
stp x10,x11,[x9],#16
stp x12,x13,[x9],#16
.endm
.macro FM o,a,b
add x0,x29,#\o
add x1,x29,#\a
add x2,x29,#\b
bl fe_mul_f
.endm
.macro MP o,a
add x0,x29,#\o
add x1,x29,#\a
ldr x2,[x29,#480]
bl fe_mul_f
.endm
.macro FQ o,a
add x0,x29,#\o
add x1,x29,#\a
bl fe_sqr_f
.endm
.macro FSUB
subs x10,x10,x20
sbcs x11,x11,x21
sbcs x12,x12,x22
sbcs x13,x13,x23
csel x27,x28,xzr,lo
subs x10,x10,x27
sbcs x11,x11,xzr
sbcs x12,x12,xzr
sbcs x13,x13,xzr
.endm
.macro FADD
adds x10,x10,x20
adcs x11,x11,x21
adcs x12,x12,x22
adcs x13,x13,x23
csel x27,x28,xzr,hs
adds x10,x10,x27
adcs x11,x11,xzr
adcs x12,x12,xzr
adcs x13,x13,xzr
.endm
.macro MVB
mov x20,x10
mov x21,x11
mov x22,x12
mov x23,x13
.endm
.macro CAN a,b,c,d
adds x14,\a,x28
adcs x15,\b,xzr
adcs x16,\c,xzr
adcs x17,\d,xzr
csel \a,x14,\a,hs
csel \b,x15,\b,hs
csel \c,x16,\c,hs
csel \d,x17,\d,hs
.endm
.macro CANA
CAN x10,x11,x12,x13
.endm
.macro ZFA
orr x24,x10,x11
orr x25,x12,x13
orr x24,x24,x25
cmp x24,#1
cset x24,lo
ldr x25,[x29,#496]
orr x25,x25,x24
str x25,[x29,#496]
.endm
.macro ROW a,b,c,d,e,f,g,h
mul x14,\d,x28
umulh x15,\d,x28
adds \f,\a,x15
adcs \g,\b,xzr
adcs \h,\c,xzr
csel x27,x28,xzr,hs
adds \e,x14,x27
adcs \f,\f,xzr
adcs \g,\g,xzr
adcs \h,\h,xzr
CAN \e,\f,\g,\h
stp \e,\f,[x9],#16
stp \g,\h,[x9],#16
.endm
.macro NEGT
mul x14,x12,x28
umulh x15,x12,x28
mul x16,x13,x28
umulh x17,x13,x28
adds x19,x15,x16
adc x24,x17,xzr
adds x18,x14,x28
adcs x19,x19,xzr
adc x24,x24,xzr
adds x24,x24,x10
adcs x25,x11,xzr
adcs x14,x12,xzr
adc x16,x13,xzr
sub x5,xzr,x14
sub x6,xzr,x16
stp x5,x6,[x9],#16
.endm
// W7: window k (bits o_k, 37 or 36 wide) of value n at 512 v_k + 64k + n 2^46 (EXTR or LSR, BFI #9); slot 7 zero
.macro OW r,k,w
add \r,x24,#\k
bfi \r,x14,#9,#\w
.endm
.macro OFFS i
mov x24,#\i
mov x15,x24
bfi x15,x10,#9,#37
extr x14,x11,x10,#37
OW x16,64,37
stp x15,x16,[x9],#16
lsr x14,x11,#10
OW x15,128,37
extr x14,x12,x11,#47
OW x16,192,37
stp x15,x16,[x9],#16
lsr x14,x12,#20
OW x15,256,36
extr x14,x13,x12,#56
OW x16,320,36
stp x15,x16,[x9],#16
lsr x14,x13,#28
OW x15,384,36
stp x15,xzr,[x9],#16
.endm
.macro K r,a,b,c,d
mov \r,#\a
movk \r,#\b,lsl #16
movk \r,#\c,lsl #32
movk \r,#\d,lsl #48
.endm
cp:
LA 0,16,x1
LB 32,48,x1
SA 0,16,x0
SB 32,48,x0
ret
ygen:
mov x29,x6
stp x30,x0,[x29,#416]
stp x1,x2,[x29,#432]
str x5,[x29,#552]
stp xzr,xzr,[x29,#496]
mov x28,#977
movk x28,#1,lsl #32
str x28,[x29,#568]
K x20,0x1ee,0x7195,0x6c28,0xc139
K x21,0x8995,0x12f5,0x4975,0x9cf0
K x22,0x34e9,0xac34,0x479e,0x6e64
K x23,0x710,0x657c,0x6a2b,0x7ae9
SB 0,16,x29
add x0,x29,#1024
mov x1,x3
bl cp
add x0,x29,#32
mov x1,x4
bl cp
mov x9,#1
str x9,[x29,#512]
1:
ldr x9,[x29,#512]
ldr x10,[x29,#432]
cmp x9,x10
b.hs 2f
sub x10,x10,x9
cmp x9,x10
csel x11,x9,x10,lo
lsl x12,x9,#6
stp x11,x12,[x29,#448]
bl ba
bl dbl
ldr x9,[x29,#512]
lsl x9,x9,#1
str x9,[x29,#512]
b 1b
2:
add x1,x29,#1024
ldr x2,[x29,#424]
stp x1,x2,[x29,#464]
ldr x11,[x29,#432]
add x11,x1,x11,lsl #6
str x11,[x29,#520]
3:
ldp x1,x2,[x29,#464]
bl yrec
ldp x1,x2,[x29,#464]
add x1,x1,#64
add x2,x2,#512
stp x1,x2,[x29,#464]
ldr x11,[x29,#520]
cmp x1,x11
b.ne 3b
ldr x9,[x29,#424]
sub x9,x9,#512
str x9,[x29,#504]
add x0,x29,#32
ldr x1,[x29,#552]
bl cp
ldr x9,[x29,#432]
stp x9,xzr,[x29,#448]
ldr x9,[x29,#440]
sub x9,x9,#1
str x9,[x29,#560]
cbz x9,5f
4:
ldr x9,[x29,#432]
ldr x10,[x29,#504]
add x10,x10,x9,lsl #10
str x10,[x29,#504]
bl ba
ldr x9,[x29,#560]
subs x9,x9,#1
str x9,[x29,#560]
b.ne 4b
5:
ldr x0,[x29,#496]
ldr x30,[x29,#416]
ret
// ba: dst[m] = src[m] + F for m < CNT, one batch inversion; REC != 0: also writes records (descending)
ba:
str x30,[x29,#528]
add x9,x29,#1024
ldr x8,[x29,#432]
add x8,x9,x8,lsl #6
sub x7,x8,#32
str x7,[x29,#576]
mov x7,#1
stp x7,xzr,[x8]
stp xzr,xzr,[x8,#16]
ldr x7,[x29,#448]
add x7,x8,x7,lsl #5
str x7,[x29,#520]
1:
LA 32,48,x29
LB 0,16,x9
C0
FSUB
SA 192,208,x29
add x0,x8,#32
mov x1,x8
add x2,x29,#192
add x9,x9,#64
stp x9,x0,[x29,#464]
bl fe_mul_f
ldp x9,x8,[x29,#464]
ldr x7,[x29,#520]
cmp x8,x7
b.ne 1b
LA 0,16,x8
C0
CANA
ZFA
add x0,x29,#96
mov x1,x8
bl fe_inv
ldp x9,x8,[x29,#464]
sub x9,x9,#64
sub x8,x8,#32
stp x9,x8,[x29,#464]
2:
ldr x1,[x29,#472]
add x0,x29,#128
add x2,x29,#96
bl fe_mul_f
ldr x9,[x29,#464]
LA 32,48,x29
LB 0,16,x9
C0
FSUB
SA 192,208,x29
LA 64,80,x29
LB 32,48,x9
FSUB
SA 224,240,x29
FM 96,96,192
FM 160,224,128
FQ 256,160
ldr x9,[x29,#464]
LA 256,272,x29
LB 0,16,x9
C0
FSUB
LB 32,48,x29
FSUB
CANA
SA 288,304,x29
LA 0,16,x9
LB 288,304,x29
FSUB
SA 320,336,x29
FM 352,160,320
ldr x9,[x29,#464]
LA 352,368,x29
LB 32,48,x9
C0
FSUB
CANA
LB 288,304,x29
ldr x7,[x29,#456]
add x7,x9,x7
SB 0,16,x7
SA 32,48,x7
ldr x2,[x29,#504]
cbz x2,3f
sub x3,x2,#512
str x3,[x29,#504]
mov x1,x9
bl yrec
3:
ldp x9,x8,[x29,#464]
sub x9,x9,#64
sub x8,x8,#32
stp x9,x8,[x29,#464]
ldr x7,[x29,#576]
cmp x8,x7
b.ne 2b
ldr x30,[x29,#528]
ret
// dbl: F <- 2F with its own inversion; flag |= (y(F) == 0)
dbl:
str x30,[x29,#536]
LA 64,80,x29
C0
ZFA
MVB
FADD
SA 224,240,x29
add x0,x29,#256
add x1,x29,#224
bl fe_inv
FQ 288,32
LA 288,304,x29
C0
CANA
MVB
FADD
FADD
SA 288,304,x29
FM 320,288,256
FQ 352,320
LA 352,368,x29
LB 32,48,x29
C0
FSUB
FSUB
CANA
SA 288,304,x29
LA 32,48,x29
LB 288,304,x29
FSUB
SA 224,240,x29
FM 352,320,224
LA 352,368,x29
LB 64,80,x29
C0
FSUB
CANA
LB 288,304,x29
SB 32,48,x29
SA 64,80,x29
ldr x30,[x29,#536]
ret
// yrec: x1 = point (canonical), x2 = record (words 0..463 through x9); preserves x29, x30
yrec:
str x30,[x29,#544]
stp xzr,xzr,[x2],#16
stp x1,x2,[x29,#480]
add x0,x29,#224
add x1,x1,#32
bl fe_sqr_f
FQ 256,224
ldr x1,[x29,#480]
add x0,x29,#288
bl fe_sqr_f
FQ 320,288
MP 352,288
MP 384,0
ldr x9,[x29,#488]
C0
LA 224,240,x29
CANA
SAP
LA 256,272,x29
CANA
SAP
LA 320,336,x29
CANA
SAP
ldr x7,[x29,#480]
LA 0,16,x7
NEGT
LA 32,48,x7
SAP
ROW x10,x11,x12,x13,x20,x21,x22,x23
ROW x20,x21,x22,x23,x10,x11,x12,x13
ROW x10,x11,x12,x13,x20,x21,x22,x23
LA 384,400,x29
CANA
NEGT
LA 0,16,x7
OFFS 0
LA 352,368,x29
CANA
OFFS 70368744177664
LA 288,304,x29
CANA
MVB
FADD
FADD
CANA
OFFS 140737488355328
ldr x30,[x29,#544]
ret

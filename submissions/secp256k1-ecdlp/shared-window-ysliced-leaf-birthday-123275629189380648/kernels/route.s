.text
.p2align 2
.globl gate
gate:
subs xzr,x0,x2
sbcs xzr,x1,x3
cset x0,lo
ret
.globl tbld
tbld:
mov x7,#977
movk x7,#1,lsl #32
mov x8,x4
1:
ldp x9,x10,[x1],#16
ldp x11,x12,[x1],#16
ldp x13,x14,[x1],#16
ldp x15,x16,[x1],#16
mov x17,x0
tst x6,#1
csel x18,x3,x2,ne
lsr x6,x6,#1
2:
stp x13,x14,[x17]
stp x15,x16,[x17,#16]
add x17,x17,x5
adds x13,x13,x9
adcs x14,x14,x10
adcs x15,x15,x11
adcs x16,x16,x12
cset x19,hs
adds x20,x13,x7
adcs x21,x14,xzr
adcs x22,x15,xzr
adcs x23,x16,xzr
cinc x19,x19,hs
cmp x19,#0
csel x13,x20,x13,ne
csel x14,x21,x14,ne
csel x15,x22,x15,ne
csel x16,x23,x16,ne
subs x18,x18,#1
b.ne 2b
add x0,x0,#32
subs x8,x8,#1
b.ne 1b
ret
.globl tinit
tinit:
1:
str x2,[x1],#8
subs x0,x0,#1
b.ne 1b
ret
.macro LNONE
.endm
.macro LADD
add x12,x2,#32
.endm
.macro HCHGI u,v,a,l,i
bfi x10,\u,#3,#50
swp x13,x11,[x10]
stp x11,\u,[x13],#\i
cbz x11,he\@
hv\@:
ldp x5,x6,[x11]
cmp x6,\u
b.eq hx\@
hy\@:
cbz x5,he\@
ldp x11,x6,[x5]
cmp x6,\u
b.eq hz\@
hw\@:
cbnz x11,hv\@
.pushsection .text.hit
hx\@:
add x6,x15,x11,lsl #1
ldr x12,[x6]
eor x12,x12,\v
tst x12,#0xfffffc0000000000
b.ne hy\@
cmp x3,x4
b.hs he\@
\l
stp x6,\a,[x3],#16
str x3,[x9,#-160]
b he\@
hz\@:
add x6,x15,x5,lsl #1
ldr x12,[x6]
eor x12,x12,\v
tst x12,#0xfffffc0000000000
b.ne hw\@
cmp x3,x4
b.hs he\@
\l
stp x6,\a,[x3],#16
str x3,[x9,#-160]
b he\@
.popsection
he\@:
.endm
.macro Z1
dc zva,x1
add x1,x1,#2048
.endm
.macro Z8
Z1
Z1
Z1
Z1
Z1
Z1
Z1
Z1
.endm
.macro Z64
Z8
Z8
Z8
Z8
Z8
Z8
Z8
Z8
.endm
.globl zero16384
zero16384:
1:
Z64
subs x0,x0,#16384
b.ne 1b
ret
.globl hhclr
hhclr:
cmp x1,x2
b.ge 2f
1:
ldr x17,[x1,#8]
add x1,x1,#32
bfi x10,x17,#3,#50
str xzr,[x10]
cmp x1,x2
b.lt 1b
2:
ret
.globl psel
psel:
stp x19,x20,[sp,#-16]!
mov x3,#1024
cmp x1,x3
csel x1,x1,x3,lo
mov x3,#0
cbz x1,6f
1:
ldp x4,x5,[x0]
ldp x6,x7,[x0,#16]
orr x14,x4,x5
orr x15,x6,x7
orr x14,x14,x15
cbz x14,5f
ldp x8,x9,[x0,#32]
ldp x10,x11,[x0,#48]
subs x12,x8,x4
sbcs x13,x9,x5
sbcs x14,x10,x6
sbcs x15,x11,x7
b.lo 5f
orr x13,x13,x14
orr x13,x13,x15
cbnz x13,5f
cbz x12,5f
tbnz x12,#63,5f
ldp x12,x13,[x0,#64]
ldp x19,x20,[x0,#80]
eor x16,x12,x19
eor x17,x13,x20
orr x16,x16,x17
cbz x16,5f
cbz x3,4f
mov x14,x2
mov x15,x3
2:
ldp x16,x17,[x14,#64]
eor x18,x12,x16
eor x16,x13,x17
orr x18,x18,x16
cbz x18,5f
ldp x16,x17,[x14,#64]
eor x18,x19,x16
eor x16,x20,x17
orr x18,x18,x16
cbz x18,5f
ldp x16,x17,[x14,#80]
eor x18,x12,x16
eor x16,x13,x17
orr x18,x18,x16
cbz x18,5f
ldp x16,x17,[x14,#80]
eor x18,x19,x16
eor x16,x20,x17
orr x18,x18,x16
cbz x18,5f
mov x18,#2
3:
ldp x16,x17,[x14]
eor x16,x16,x4
eor x17,x17,x5
orr x16,x16,x17
cbnz x16,7f
ldp x16,x17,[x14,#16]
eor x16,x16,x6
eor x17,x17,x7
orr x16,x16,x17
cbz x16,5f
7:
ldp x16,x17,[x14]
eor x16,x16,x8
eor x17,x17,x9
orr x16,x16,x17
cbnz x16,8f
ldp x16,x17,[x14,#16]
eor x16,x16,x10
eor x17,x17,x11
orr x16,x16,x17
cbz x16,5f
8:
add x14,x14,#32
subs x18,x18,#1
b.ne 3b
add x14,x14,#64
subs x15,x15,#1
b.ne 2b
4:
mov x14,x0
add x15,x2,x3,lsl #7
mov x18,#8
10:
ldp x16,x17,[x14],#16
stp x16,x17,[x15],#16
subs x18,x18,#1
b.ne 10b
add x3,x3,#1
cmp x3,#4
b.eq 9f
5:
add x0,x0,#128
subs x1,x1,#1
b.ne 1b
6:
mov x0,#0
ldp x19,x20,[sp],#16
ret
9:
mov x0,#1
ldp x19,x20,[sp],#16
ret
.macro WR4
str x2,[x0]
add x2,x2,#8192
str x2,[x0,#8192]
add x2,x2,#8192
str x2,[x0,#16384]
add x2,x2,#8192
str x2,[x0,#24576]
add x2,x2,#8192
add x0,x0,#32768
.endm
.globl wreset
wreset:
mov x2,#8160
add x2,x0,x2
1:
WR4
WR4
WR4
WR4
WR4
WR4
WR4
WR4
cbhi x1,x0,1b
ret
.macro CH7 b,p,c
r7\p\c:
\b
r6\p\c:
\b
r5\p\c:
\b
r4\p\c:
\b
r3\p\c:
\b
r2\p\c:
\b
r1\p\c:
\b
.endm
.macro CH8 b,p,c
r8\p\c:
\b
CH7 \b,\p,\c
.endm
.macro CH127 b,c
CH7 \b,15,\c
CH8 \b,14,\c
CH8 \b,13,\c
CH8 \b,12,\c
CH8 \b,11,\c
CH8 \b,10,\c
CH8 \b,9,\c
CH8 \b,8,\c
CH8 \b,7,\c
CH8 \b,6,\c
CH8 \b,5,\c
CH8 \b,4,\c
CH8 \b,3,\c
CH8 \b,2,\c
CH8 \b,1,\c
CH8 \b,0,\c
.endm
.macro RHEAD
add x24,x24,#8192
ldr x5,[x24]
ldr x7,[x25,x5,lsl #3]
add x1,x24,x26
br x7
.endm
.macro RS k2,k3,w2,w3
and x0,x10,\k2,lsl #13
ldadd x4,x8,[x0]
cbne x8,x0,rs\@
str x0,[x0]
b rt\@
rs\@:
stp \k2,\k3,[x8]
stp \w2,\w3,[x8,#16]
rt\@:
.endm
.macro RP
ld64b x16,[x1]
sub x1,x1,#64
RS x16,x17,x18,x19
RS x20,x21,x22,x23
.endm
.globl rsw
rsw:
RHEAD
CH127 RP,e
r0e:
RHEAD
CH127 RP,o
r0o:
ldp x16,x17,[x1,#32]
ldp x18,x19,[x1,#48]
RS x16,x17,x18,x19
RHEAD
rsx:
ret
.macro HP
ld64b x16,[x1]
sub x1,x1,#64
st64b x16,[x2]
HCHGI x17,x16,x2,LNONE,16
HCHGI x21,x20,x12,LADD,16
add x2,x2,#64
.endm
.globl hsw
hsw:
mov x9,#0x80000000000000
ldp x2,x7,[x9,#-176]
ldp x3,x4,[x9,#-160]
ldp x14,x15,[x9,#-144]
mov x10,#0xff80000000000000
cmp x2,x7
b.ge hsy
add x13,x14,x2,asr #1
RHEAD
CH127 HP,f
r0f:
RHEAD
CH127 HP,g
r0g:
ldp x16,x17,[x1,#32]
ldp x18,x19,[x1,#48]
st64b x16,[x2]
HCHGI x17,x16,x2,LNONE,32
add x2,x2,#64
RHEAD
hsx:
str x2,[x9,#-176]
hsy:
ret
.macro ESW h
str x30,[sp,#-16]!
1:
ldr x5,[x24]
cbne x5,x9,ea\h
fa\h:
ldr x5,[x24,#8192]
cbne x5,x9,eb\h
fb\h:
ldr x5,[x24,#16384]
cbne x5,x9,ec\h
fc\h:
ldr x5,[x24,#24576]
cbne x5,x9,ed\h
fd\h:
add x24,x24,#32768
cbhi x13,x24,1b
ldr x30,[sp],#16
ret
ea\h:
mov x1,x24
bl \h
b fa\h
eb\h:
add x1,x24,#8192
bl \h
b fb\h
ec\h:
add x1,x24,#16384
bl \h
b fc\h
ed\h:
add x1,x24,#24576
bl \h
b fd\h
.endm
.macro RB1
sub x1,x1,#64
ld64b x16,[x1]
sub x8,x8,#64
st64b x16,[x8]
.endm
.macro HB1
sub x1,x1,#64
ld64b x16,[x1]
add x8,x8,#64
st64b x16,[x8]
.endm
.globl rswe
rswe:
ESW rrf
rrf:
sub x6,x9,x5
add x1,x1,x26
2:
ld64b x16,[x1]
and x0,x10,x16,lsl #13
ldadd x27,x8,[x0]
sub x8,x8,#32
st64b x16,[x8]
RB1
RB1
RB1
RB1
sub x1,x1,#64
subs x6,x6,#320
b.ne 2b
ret
.globl hswe
hswe:
ESW hrf
hrf:
sub x6,x9,x5
add x1,x1,x26
2:
ld64b x16,[x1]
mov x8,x17
st64b x16,[x8]
HB1
HB1
HB1
HB1
sub x1,x1,#64
subs x6,x6,#320
b.ne 2b
ret
.macro QT k,h
ldr x5,[x12,#\k]
cbne x5,x9,q\k\h
p\k\h:
.endm
.macro QS k,h
q\k\h:
add x1,x12,#\k
bl \h
b p\k\h
.endm
.macro QSW h
str x30,[sp,#-16]!
1:
QT 0,\h
QT 128,\h
QT 256,\h
QT 384,\h
QT 512,\h
QT 640,\h
QT 768,\h
QT 896,\h
QT 1024,\h
QT 1152,\h
QT 1280,\h
QT 1408,\h
QT 1536,\h
QT 1664,\h
QT 1792,\h
QT 1920,\h
add x12,x12,#2048
cbhi x13,x12,1b
ldr x30,[sp],#16
ret
QS 0,\h
QS 128,\h
QS 256,\h
QS 384,\h
QS 512,\h
QS 640,\h
QS 768,\h
QS 896,\h
QS 1024,\h
QS 1152,\h
QS 1280,\h
QS 1408,\h
QS 1536,\h
QS 1664,\h
QS 1792,\h
QS 1920,\h
.endm
.globl rqsw
rqsw:
QSW rqf
rqf:
sub x6,x9,x5
add x7,x1,#96
sub x8,x1,x11
lsr x8,x8,#7
2:
ldp x2,x3,[x7]
ldp x16,x17,[x7,#16]
add x0,x10,x17,lsl #7
ldadd x4,x14,[x0]
cbne x14,x0,3f
str x0,[x0]
b 4f
3:
stp x2,x3,[x14]
stp x16,x8,[x14,#16]
4:
sub x7,x7,#32
subs x6,x6,#32
b.ne 2b
ret
.globl gqsw
gqsw:
QSW gqf
gqf:
sub x6,x9,x5
add x7,x1,#96
sub x8,x1,x11
lsl x8,x8,#6
2:
ldp x2,x3,[x7]
ldp x16,x17,[x7,#16]
stp x17,x16,[x19]
stp x2,x3,[x19,#16]
stp x8,xzr,[x19,#32]
stp xzr,xzr,[x19,#48]
add x19,x19,#64
sub x7,x7,#32
subs x6,x6,#32
b.ne 2b
ret
.macro GB1
ld64b x16,[x2]
sub x5,x5,#64
st64b x16,[x5]
.endm
.globl grpl
grpl:
ld64b x16,[x1]
ldadd x27,x5,[x20]
sub x5,x5,#32
st64b x16,[x5]
GB1
add x2,x2,#64
GB1
add x2,x2,#64
GB1
add x2,x2,#64
GB1
ret
.macro HXA
ubfx x6,x2,#3,#25
ubfx x7,x3,#1,#41
.endm
.macro HXM
and x6,x2,#0x1ffffff
extr x7,x3,x2,#25
and x7,x7,#0x1ffffffffff
.endm
.macro HRQ x
ldp x2,x3,[x14,#16]
str xzr,[x13,#8]
\x
add x0,x10,x6,lsl #7
ldadd x4,x5,[x0]
cbne x5,x0,hq\@
str x0,[x0]
b hr\@
hq\@:
stp x2,x3,[x5]
stp x13,x7,[x5,#16]
hr\@:
add x13,x13,#320
.endm
.macro HRL x
1:
ldp x14,x15,[x12],#16
HRQ \x
mov x14,x15
HRQ \x
cmp x12,x9
b.lo 1b
ret
.endm
.globl hreq
hreq:
HRL HXA
.globl hreqm
hreqm:
HRL HXM
.macro YB1 s,d
ld64b x16,[\s]
st64b x16,[\d]
add \s,\s,#64
add \d,\d,#64
.endm
.macro YB7 s,d
YB1 \s,\d
YB1 \s,\d
YB1 \s,\d
YB1 \s,\d
YB1 \s,\d
YB1 \s,\d
YB1 \s,\d
.endm
.globl yrst
yrst:
1:
str x2,[x0]
add x0,x0,x14
add x2,x2,x14
str x2,[x0]
add x0,x0,x14
add x2,x2,x14
str x2,[x0]
add x0,x0,x14
add x2,x2,x14
str x2,[x0]
add x0,x0,x14
add x2,x2,x14
cbhi x1,x0,1b
ret
.globl yprod
yprod:
1:
ldr x27,[x15,#-8]
lsr x6,x27,#27
add x0,x7,x6,lsl #30
ldadd x4,x5,[x0]
str x27,[x5],#64
YB1 x15,x5
YB1 x15,x5
YB1 x15,x5
YB1 x15,x5
ld64b x16,[x15]
sub x21,x24,x3,lsl #42
st64b x16,[x5]
add x15,x15,#64
add x5,x5,#64
YB1 x15,x5
YB1 x15,x5
add x15,x15,#64
subs x3,x3,#1
b.ne 1b
ret
.globl yrel
yrel:
1:
ldr x5,[x12]
cbne x5,x9,3f
2:
add x12,x12,x14
cbhi x13,x12,1b
ret
3:
sub x6,x5,x9
add x1,x12,#64
4:
ldr x27,[x1],#64
and x7,x27,#0x7ffffff
add x0,x8,x7,lsl #25
ldadd x4,x2,[x0]
cbne x2,x0,5f
str x0,[x0]
add x1,x1,#448
b 6f
5:
YB7 x1,x2
6:
subs x6,x6,#512
b.ne 4b
b 2b
.globl ym0
ym0:
str x15,[x29]
1:
ldr x5,[x12]
cbne x5,x9,3f
2:
add x12,x12,x14
cbhi x13,x12,1b
str x15,[x29]
ret
3:
sub x6,x9,x5
add x1,x12,x14
sub x1,x1,#512
4:
cbhs x15,x24,2b
YB7 x1,x15
add x15,x15,#64
sub x1,x1,#960
subs x6,x6,#512
b.ne 4b
b 2b
.globl bcopy
bcopy:
1:
YB1 x1,x2
YB1 x1,x2
YB1 x1,x2
YB1 x1,x2
subs x3,x3,#4
b.ne 1b
ret
.globl bleaf
bleaf:
str x3,[x0,#2048]
str x4,[x0,#2056]
stlr x2,[x0]
1:
ldaxr x5,[x1]
cmp x5,x2
b.eq 2f
wfe
b 1b
2:
ldr x3,[x1,#2048]
ldr x4,[x1,#2056]
ret
.globl bpar
bpar:
1:
ldaxr x5,[x6]
cmp x5,x2
b.eq 2f
wfe
b 1b
2:
ldr x8,[x6,#2048]
ldr x9,[x6,#2056]
orr x3,x3,x8
orr x4,x4,x9
add x6,x6,#4096
cmp x6,x7
b.lo 1b
ret
.globl brel
brel:
str x3,[x0,#2048]
str x4,[x0,#2056]
stlr x2,[x0]
ret
.globl trst
trst:
1:
str xzr,[x0]
add x0,x0,x14
cbhi x1,x0,1b
ret
.globl tmrg
tmrg:
str xzr,[x29]
1:
ldr x5,[x12]
cbnz x5,3f
2:
add x12,x12,x14
cbhi x13,x12,1b
sub x5,x15,x28
str x5,[x29]
ret
3:
add x1,x12,#64
4:
cbhs x15,x24,2b
YB1 x1,x15
YB1 x1,x15
subs x5,x5,#128
b.ne 4b
b 2b
.globl tbkt
tbkt:
1:
ldr x5,[x12]
cbnz x5,3f
2:
add x12,x12,x14
cbhi x13,x12,1b
ret
3:
add x1,x12,#64
4:
ldr x6,[x1,#120]
madd x0,x6,x7,x10
ldr x2,[x0]
cbhs x2,x24,5f
add x3,x2,#128
str x3,[x0]
add x2,x2,x0
add x2,x2,#64
YB1 x1,x2
YB1 x1,x2
b 6f
5:
add x1,x1,#128
6:
subs x5,x5,#128
b.ne 4b
b 2b
.globl tcpy
tcpy:
str xzr,[x29]
ldr x5,[x12]
add x1,x12,#64
mov x6,x5
cbz x5,2f
1:
YB1 x1,x15
subs x6,x6,#64
b.ne 1b
2:
str x5,[x29]
ret
.globl tscn
tscn:
1:
ldp x6,x7,[x1,#64]
ldp x8,x9,[x1,#80]
eor x6,x6,x4
eor x7,x7,x5
and x7,x7,#3
orr x6,x6,x7
eor x8,x8,x4
eor x9,x9,x5
and x9,x9,#3
orr x8,x8,x9
cbnz x6,2f
str x1,[x15],#8
2:
cbnz x8,3f
orr x10,x1,#1
str x10,[x15],#8
3:
add x1,x1,#128
subs x3,x3,#1
b.ne 1b
ret
.globl xemit
xemit:
1:
ldp x9,x10,[x1,#288]
ldp x11,x12,[x1,#304]
orr x13,x6,x12,lsr #14
extr x14,x12,x11,#14
extr x15,x11,x10,#25
and x0,x8,x15
ldadd x4,x5,[x0]
cbne x5,x0,3f
str x0,[x0]
b 2f
3:
stp x14,x13,[x5]
stp x2,x26,[x5,#16]
2:
add x26,x26,#4
add x1,x1,#512
subs x3,x3,#1
b.ne 1b
ret
.globl cpair
cpair:
ld64b x16,[x1]
st64b x16,[x2]
ret
.globl rdcur
rdcur:
ldr x0,[x3]
ret

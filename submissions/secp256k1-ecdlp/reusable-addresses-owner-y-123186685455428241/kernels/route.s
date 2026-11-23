.text
.p2align 2
.globl gate
gate:
subs xzr,x0,x2
sbcs xzr,x1,x3
cset x0,lo
ret
.globl dlv
dlv:
ldp x1,x2,[x0]
ldp x3,x4,[x0,#16]
ldp x5,x6,[x0,#32]
stp x3,x4,[x2]
stp x5,x6,[x2,#16]
staddl x28,[x29]
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
.globl fin1
fin1:
1:
ldr x9,[x7]
and x10,x9,#0xffffffffffffe000
sub x11,x9,x10
sub x11,x11,#32
ldr x5,[x10,#16]
add x5,x5,#1
stp x11,x5,[x10,#8]
stlr x7,[x10]
dmb ish
str xzr,[x10,#16]
add x9,x10,#32
str x9,[x7],#8
subs x1,x1,#1
b.ne 1b
ret
.globl fin
fin:
1:
ldr x9,[x7]
ldr x5,[x3],#8
subs x10,x9,x8
b.eq 2f
cmp x5,x6
b.ge 2f
stp x10,x11,[x8,#-16]
str x5,[x8,#-24]
stlur x7,[x8,#-32]
2:
add x5,x5,x10
cmp x5,x6
csel x5,x5,x6,lt
stp x7,x2,[x0]
stp x5,xzr,[x0,#16]
str xzr,[x0,#32]
stlur xzr,[x0,#40]
dmb ish
add x7,x7,#8
add x8,x8,#0x3000
subs x1,x1,#1
b.ne 1b
ret
.macro H2B a,b,c,d,l
bfi x10,\a,#3,#42
ldadd x17,x11,[x10]
kr\l:
stp \a,\b,[x11]
stp \c,\d,[x11,#16]
.endm
.macro H2R1 l
ldp x20,x21,[x1],#16
ldp x22,x23,[x1],#16
H2B x20,x21,x22,x23,\l
.endm
.macro H2P l
ld64b x20,[x1]
add x1,x1,#64
H2B x20,x21,x22,x23,\la
H2B x24,x25,x26,x27,\lb
.endm
.macro H2P8 l
H2P \l0
H2P \l1
H2P \l2
H2P \l3
H2P \l4
H2P \l5
H2P \l6
H2P \l7
.endm
.macro H2PR w,l
ld64b x20,[\w]
H2B x20,x21,x22,x23,\la
H2B x24,x25,x26,x27,\lb
.endm
.globl h2s0
h2s0:
mov x10,#0x200000000000
mov x18,#8192
mov x17,#32
movn x16,#0xa7
movk x16,#0,lsl #48
stp x8,x0,[x16,#-16]
movn x29,#0xfff
movk x29,#0,lsl #48
mov x13,xzr
mov x14,#0x4000000000000000
mov x15,x8
mov x8,x0
mov x28,xzr
mov x0,#0x7fff0000000000
add x0,x0,#64
movn x6,#1023
movk x6,#0,lsl #48
mov sp,x6
ret
.globl h2s
h2s:
ldp x20,x21,[x0,#-24]
eor x9,x20,x20
ldr x19,[x16,x9]
tbnz x20,#12,7f
adds x13,x13,x21
adc x14,x14,xzr
subs x15,x15,#1
sbc x8,x8,xzr
orr x9,x15,x8
cbnz x9,2f
mov x6,#0x4000000000000000
sub x14,x14,x6
2:
cbz x20,6f
mov x1,x0
mov x2,x20
tbz x2,#6,3f
H2P q6
3:
tbz x2,#7,4f
H2P q7
H2P q8
4:
tbz x2,#8,5f
H2P q9
H2P q10
H2P q11
H2P q12
5:
tbz x2,#9,h2a
H2P8 qa
h2a:
tbz x2,#10,h2b
H2P8 qb
H2P8 qc
h2b:
tbz x2,#11,h2c
H2P8 qd
H2P8 qe
H2P8 qf
H2P8 qg
h2c:
tbz x2,#5,1f
H2R1 q5
1:
add x1,x0,#64
add x2,x0,#128
6:
subs x13,x13,#1
sbc x14,x14,xzr
orr x9,x13,x14
cbz x9,h2sd
ret
7:
ld64b x0,[sp]
H2PR x0,f0
H2PR x1,f1
H2PR x2,f2
H2PR x3,f3
H2PR x4,f4
H2PR x5,f5
H2PR x6,f6
H2PR x7,f7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,fa0
H2PR x1,fa1
H2PR x2,fa2
H2PR x3,fa3
H2PR x4,fa4
H2PR x5,fa5
H2PR x6,fa6
H2PR x7,fa7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,fb0
H2PR x1,fb1
H2PR x2,fb2
H2PR x3,fb3
H2PR x4,fb4
H2PR x5,fb5
H2PR x6,fb6
H2PR x7,fb7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,fc0
H2PR x1,fc1
H2PR x2,fc2
H2PR x3,fc3
H2PR x4,fc4
H2PR x5,fc5
H2PR x6,fc6
H2PR x7,fc7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,fd0
H2PR x1,fd1
H2PR x2,fd2
H2PR x3,fd3
H2PR x4,fd4
H2PR x5,fd5
H2PR x6,fd6
H2PR x7,fd7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,fe0
H2PR x1,fe1
H2PR x2,fe2
H2PR x3,fe3
H2PR x4,fe4
H2PR x5,fe5
H2PR x6,fe6
H2PR x7,fe7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,ff0
H2PR x1,ff1
H2PR x2,ff2
H2PR x3,ff3
H2PR x4,ff4
H2PR x5,ff5
H2PR x6,ff6
H2PR x7,ff7
add sp,sp,#64
ld64b x0,[sp]
H2PR x0,fg0
H2PR x1,fg1
H2PR x2,fg2
H2PR x3,fg3
H2PR x4,fg4
H2PR x5,fg5
H2PR x6,fg6
H2PR x7,fg7
sub sp,sp,#448
ldr x0,[sp]
subs x13,x13,#1
sbc x14,x14,xzr
orr x9,x13,x14
cbz x9,h2sd
ret
h2sd:
add x28,x28,#1
stlr x28,[x29]
mov x14,#0x4000000000000000
ldp x15,x8,[x16,#-16]
ret
.macro LNONE
.endm
.macro LADD
add x12,x24,#32
.endm
.macro LSUB
sub x12,x24,#32
.endm
.macro HCHGI u,v,off
bfi x10,\u,#3,#50
swp x24,x11,[x10]
cbnz x11,h31_cold\@
.pushsection .text.hit
h31_cold\@:
tbnz x11,#3,h31_firsttag\@
cbeq x11,x24,h31_selfdirect\@
stp x11,x24,[x13],#16
sub x12,x13,#8
str x12,[x10]
h31_firstdirect_p1\@:
ldp x22,x23,[x11,#32]
cbeq x23,\u,h31_firstdirect_hit1\@
h31_firstdirect_p0\@:
ldp x22,x23,[x11,#0]
cbeq x23,\u,h31_firstdirect_hit0\@
h31_firstdirect_end\@:
b h31_done\@
h31_firsttag\@:
ldp x18,x19,[x11,#-8]
cbeq x19,x24,h31_selftag\@
stp x11,x24,[x13],#16
sub x12,x13,#8
str x12,[x10]
mov x11,x18
h31_firsttagbody_p1\@:
ldp x22,x23,[x19,#32]
cbeq x23,\u,h31_firsttagbody_hit1\@
h31_firsttagbody_p0\@:
ldp x22,x23,[x19,#0]
cbeq x23,\u,h31_firsttagbody_hit0\@
h31_firsttagbody_end\@:
b h31_loop\@
h31_loop\@:
tbnz x11,#3,h31_middle\@
h31_directtail_p1\@:
ldp x22,x23,[x11,#32]
cbeq x23,\u,h31_directtail_hit1\@
h31_directtail_p0\@:
ldp x22,x23,[x11,#0]
cbeq x23,\u,h31_directtail_hit0\@
h31_directtail_end\@:
b h31_done\@
h31_middle\@:
ldp x11,x19,[x11,#-8]
h31_middlebody_p1\@:
ldp x22,x23,[x19,#32]
cbeq x23,\u,h31_middlebody_hit1\@
h31_middlebody_p0\@:
ldp x22,x23,[x19,#0]
cbeq x23,\u,h31_middlebody_hit0\@
h31_middlebody_end\@:
b h31_loop\@
h31_selfdirect\@:
h31_selfdirectbody_p0\@:
ldp x22,x23,[x11,#0]
cbeq x23,\u,h31_selfdirectbody_hit0\@
h31_selfdirectbody_end\@:
b h31_done\@
h31_selftag\@:
str x11,[x10]
mov x11,x18
h31_selftagbody_p0\@:
ldp x22,x23,[x19,#0]
cbeq x23,\u,h31_selftagbody_hit0\@
h31_selftagbody_end\@:
b h31_loop\@
h31_firstdirect_hit1\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_firstdirect_p0\@
cmp x14,x15
b.hs h31_done\@
add x12,x11,#32
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_firstdirect_hit0\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_firstdirect_end\@
cmp x14,x15
b.hs h31_done\@
add x12,x11,#0
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_firsttagbody_hit1\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_firsttagbody_p0\@
cmp x14,x15
b.hs h31_done\@
add x12,x19,#32
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_firsttagbody_hit0\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_firsttagbody_end\@
cmp x14,x15
b.hs h31_done\@
add x12,x19,#0
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_directtail_hit1\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_directtail_p0\@
cmp x14,x15
b.hs h31_done\@
add x12,x11,#32
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_directtail_hit0\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_directtail_end\@
cmp x14,x15
b.hs h31_done\@
add x12,x11,#0
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_middlebody_hit1\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_middlebody_p0\@
cmp x14,x15
b.hs h31_done\@
add x12,x19,#32
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_middlebody_hit0\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_middlebody_end\@
cmp x14,x15
b.hs h31_done\@
add x12,x19,#0
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_selfdirectbody_hit0\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_selfdirectbody_end\@
cmp x14,x15
b.hs h31_done\@
add x12,x11,#0
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
h31_selftagbody_hit0\@:
eor x12,x22,\v
tst x12,#0xfffffe0000000000
b.ne h31_selftagbody_end\@
cmp x14,x15
b.hs h31_done\@
add x12,x19,#0
add x18,x24,#\off
stp x12,x18,[x14],#16
str x14,[x9,#-160]
b h31_done\@
.popsection
h31_done\@:
.endm
.macro DLVS
ldp x16,x17,[x1],#32
ldp x18,x19,[x1,#-16]
stp x18,x19,[x24,#16]
stp x16,x17,[x24],#32
.endm
.macro DLVP
ld64b x16,[x1]
add x1,x1,#64
st64b x16,[x24]
add x24,x24,#64
.endm
.macro DLVHS
DLVS
mov x26,x24
mov x24,x25
stp x16,x17,[x24]
stp x18,x19,[x24,#16]
eor x12,x17,#1
stp xzr,x12,[x24,#32]
add x25,x25,#64
HCHGI x17,x16,0
mov x24,x26
.endm
.macro DLVHP
ld64b x16,[x1]
add x1,x1,#64
st64b x16,[x24]
cbz x24,hpzero\@
HCHGI x17,x16,0
HCHGI x21,x20,32
add x24,x24,#64
b hpend\@
hpzero\@:
mov x26,x24
mov x24,x25
st64b x16,[x24]
add x25,x25,#64
HCHGI x17,x16,0
HCHGI x21,x20,32
mov x24,x26
add x24,x24,#64
hpend\@:
.endm
.macro DLVHP8
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
.endm
.macro DLVHPL
ld64b x16,[x1]
st64b x16,[x24]
HCHGI x17,x16,0
HCHGI x21,x20,32
.endm
.macro DLVHPR w
ld64b x16,[\w]
st64b x16,[x24]
HCHGI x17,x16,0
HCHGI x21,x20,32
add x24,x24,#64
.endm
.macro DLVHPRZ w
ld64b x16,[\w]
st64b x16,[x24]
mov x26,x24
mov x24,x25
st64b x16,[x24]
add x25,x25,#64
HCHGI x17,x16,0
HCHGI x21,x20,32
mov x24,x26
add x24,x24,#64
.endm
.macro SEGL s,p,l
cbz x7,sge\l
tbz x1,#5,sga\l
\s
sub x7,x7,#32
sga\l:
lsr x8,x7,#6
cbz x8,sge2\l
sgb\l:
\p
subs x8,x8,#1
b.ne sgb\l
sge2\l:
tbz x7,#5,sge\l
\s
sge\l:
.endm
.globl dlvh0
dlvh0:
mov x9,#0x80000000000000
ldp x14,x15,[x9,#-160]
ldp x13,x25,[x9,#-144]
mov x10,#0xff80000000000000
mov x0,#0x7fff0000000000
add x0,x0,#64
movn x6,#4095
movk x6,#127,lsl #48
mov sp,x6
mov x28,#1
ret
.globl dlvh
dlvh:
ldp x24,x16,[x0,#-24]
ldr x18,[x0,#-8]
ldr x17,[x18]
cmp x24,x17
b.lt 7f
tbnz x16,#13,6f
7:
mov x1,x0
csel x3,x17,x24,lt
sub x7,x3,x24
mov x4,x16
mov x5,x24
SEGL DLVS,DLVP,1
add x16,x5,x4
sub x7,x16,x3
SEGL DLVHS,DLVHP,2
mov x7,#8192
sub x7,x7,x4
SEGL DLVS,DLVP,3
staddl x28,[x29]
ret
6:
cbz x24,dlvh_zero_full
ld64b x0,[sp]
DLVHPR x0
dlvh_after_first:
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
add sp,sp,#64
ld64b x0,[sp]
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
ld64b x16,[x7]
st64b x16,[x24]
HCHGI x17,x16,0
HCHGI x21,x20,32
sub sp,sp,#960
ldr x0,[sp]
staddl x28,[x29]
ret
dlvh_zero_full:
ld64b x0,[sp]
DLVHPRZ x0
b dlvh_after_first
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
.macro LD64NEXT
add x15,x15,#64
ld64b x16,[x15]
.endm
.macro ST64NEXT
st64b x16,[x15]
add x15,x15,#64
.endm
.macro REQPUT
stp xzr,xzr,[x0]
stp x2,x3,[x0,#16]
str x10,[x0,#32]
stlur x13,[x0,#40]
dmb ish
add x13,x13,#256
.endm
.macro LINE64END
dmb ish
add x1,x1,#64
cmp x1,x6
b.lo 1b
2:
ret
.endm
.macro OR18CBZ
orr x18,x18,x16
cbz x18,5f
.endm
.macro YCOPY a,b,c
ld64b x16,[\a]
st64b x16,[x15]
st64b x16,[\b]
add x15,x15,#\c
.endm
.globl zero16384
zero16384:
1:
Z64
subs x0,x0,#16384
b.ne 1b
ret
.globl init1
init1:
1:
str xzr,[x2,#-16]
str x2,[x1],#8
add x2,x2,x3
subs x0,x0,#1
b.ne 1b
ret
.globl init3
init3:
1:
str x2,[x1],#8
add x2,x2,x3
subs x0,x0,#1
b.ne 1b
ret
.globl dlvq
dlvq:
ldp x1,x2,[x0]
ldp x3,x4,[x0,#16]
ldp x5,x6,[x0,#32]
ldadd x9,x8,[x10]
cmp x8,x11
b.hs 2f
stp x2,x3,[x8]
stp x4,x5,[x8,#16]
str x6,[x8,#32]
b 1f
2:
str x11,[x10]
1:
ldr x1,[x29]
adds x1,x1,x28
stlr x1,[x29]
b.lo 3f
add x7,x29,#2048
ldr x2,[x7]
add x2,x2,#1
stlr x2,[x7]
3:
ret
.globl reqh
reqh:
1:
ldp x14,x15,[x12],#16
ldp x2,x3,[x14,#16]
REQPUT
ldp x2,x3,[x15,#16]
REQPUT
subs x9,x9,#1
b.ne 1b
ret
.globl fwdq
fwdq:
cmp x1,x6
b.hs 2f
1:
ldp x2,x3,[x1,#8]
and x4,x3,#3
lsl x5,x1,#8
stp x2,x4,[x0]
stp x3,x5,[x0,#16]
str x7,[x0,#32]
stlur xzr,[x0,#40]
LINE64END
.globl fwdq1
fwdq1:
cmp x1,x6
b.hs 2f
1:
ldp x2,x3,[x1,#8]
ubfx x4,x2,#3,#25
ubfx x5,x3,#1,#41
orr x4,x4,x5,lsl #25
lsr x5,x5,#39
stp x4,x5,[x0]
lsl x5,x1,#8
stp x3,x5,[x0,#16]
str x7,[x0,#32]
stlur x2,[x0,#40]
LINE64END
.globl sendr
sendr:
1:
ldp x2,x3,[x1],#16
ldp x4,x5,[x1],#16
stp x7,x8,[x0]
stp x2,x3,[x0,#16]
str x4,[x0,#32]
stlur x5,[x0,#40]
dmb ish
add x8,x8,#32
subs x9,x9,#1
b.ne 1b
ret
.globl fwdr
fwdr:
cmp x1,x6
b.hs 2f
1:
ldp x2,x3,[x1]
ldp x4,x5,[x1,#16]
ldr x11,[x1,#32]
lsr x7,x2,#8
and x2,x2,#255
ldp x8,x9,[x7,#24]
add x9,x9,x2
stp x8,x9,[x0]
stp x3,x4,[x0,#16]
str x5,[x0,#32]
stlur x11,[x0,#40]
LINE64END
.globl ydlv0
ydlv0:
mov x0,#0x220000000000
add x0,x0,#64
add x1,x0,#64
add x2,x0,#128
add x3,x0,#192
add x4,x0,#256
add x5,x0,#320
add x6,x0,#384
mov x14,#0x220000000000
movk x14,#0xb000
add x7,x14,#64
add x8,x14,#128
add x9,x14,#192
add x10,x14,#256
add x11,x14,#320
add x12,x14,#384
add x13,x14,#448
mov x28,#1
ret
.globl ysnd
ysnd:
ldr x27,[x15,#-8]
orr x16,x26,x27,lsl #25
lsr x17,x27,#39
stp x16,x17,[x14,#16]
ld64b x16,[x15]
st64b x16,[x7]
LD64NEXT
st64b x16,[x8]
LD64NEXT
st64b x16,[x9]
LD64NEXT
st64b x16,[x10]
LD64NEXT
mov x21,x24
st64b x16,[x11]
LD64NEXT
st64b x16,[x12]
LD64NEXT
st64b x16,[x13]
add x15,x15,#128
add x24,x24,x25
stlr x28,[x14]
dmb ish
ret
.globl ydlv
ydlv:
YCOPY x0,x7,64
YCOPY x1,x8,64
YCOPY x2,x9,64
YCOPY x3,x10,64
YCOPY x4,x11,64
YCOPY x5,x12,64
YCOPY x6,x13,128
stlr x28,[x14]
dmb ish
staddl x28,[x29]
ret
.globl ydlve
ydlve:
ld64b x16,[x0]
ST64NEXT
ld64b x16,[x1]
ST64NEXT
ld64b x16,[x2]
ST64NEXT
ld64b x16,[x3]
ST64NEXT
ld64b x16,[x4]
ST64NEXT
ld64b x16,[x5]
ST64NEXT
ld64b x16,[x6]
st64b x16,[x15]
add x15,x15,#128
staddl x28,[x29]
ret
.globl xemit
xemit:
1:
ldp x9,x10,[x1,#376]
ldp x11,x12,[x1,#392]
orr x13,x6,x12,lsr #15
extr x14,x12,x11,#15
extr x15,x11,x10,#36
and x0,x8,x15
ldadd x4,x5,[x0]
kx:
stp x14,x13,[x5]
stp x2,x26,[x5,#16]
add x26,x26,#4
add x1,x1,#512
subs x3,x3,#1
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
.globl wcnt
wcnt:
1:
ldaxr x2,[x0]
cmp x2,x1
b.eq 2f
wfe
b 1b
2:
ret
.globl ydsp
ydsp:
str x30,[sp,#-16]!
add x25,x0,#1984
1:
ldr x30,[x25]
ldp x26,x27,[x29,#16]
cbnz x26,8f
2:
cmp x15,x24
csel x15,x24,x15,hs
blr x30
stlr xzr,[x25]
ldp x19,x30,[x29]
cmp x19,x30
b.ne 1b
ldr x30,[sp],#16
ret
8:
mov x15,x26
mov x24,x27
str xzr,[x29,#16]
b 2b
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
OR18CBZ
ldp x16,x17,[x14,#64]
eor x18,x19,x16
eor x16,x20,x17
OR18CBZ
ldp x16,x17,[x14,#80]
eor x18,x12,x16
eor x16,x13,x17
OR18CBZ
ldp x16,x17,[x14,#80]
eor x18,x19,x16
eor x16,x20,x17
OR18CBZ
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
.macro RSTUB a,b,c,d
rs\a:
mov x9,#0x100000000000
stp \a,\b,[x9]
stp \c,\d,[x9,#16]
sub x9,x11,#8192
add x9,x9,#32
str x9,[x10]
orr \c,x10,#0x400000000000
ldadd x18,x11,[\c]
cmp x11,x12
b.lt rs8\a
str x12,[\c]
b rs9\a
rs8\a:
stp x18,x19,[x9,#-16]
str x11,[x9,#-24]
stlur x10,[x9,#-32]
dmb ish
rs9\a:
mrs x9,elr_el1
add x9,x9,#8
msr elr_el1,x9
eret
.endm
RSTUB x20,x21,x22,x23
RSTUB x24,x25,x26,x27
xst:
mov x4,#0x100000000000
stp x14,x13,[x4]
stp x2,x26,[x4,#16]
mov x4,#32
sub x16,x5,#4064
str x16,[x0]
ldr x17,[x16,#-16]
add x17,x17,#1
mov x5,#4096
stp x5,x17,[x16,#-24]
stlur x0,[x16,#-32]
dmb ish
mrs x5,elr_el1
add x5,x5,#8
msr elr_el1,x5
eret
.p2align 11,0
vbarx:
.org vbarx+0x200
mrs x4,elr_el1
add x4,x4,#0x20000
br x4
.p2align 11,0
vbarr:
.org vbarr+0x200
mrs x9,elr_el1
add x9,x9,#0x20000
br x9
.org krq6a+0x20000
b rsx20
.org krq6b+0x20000
b rsx24
.org krq7a+0x20000
b rsx20
.org krq7b+0x20000
b rsx24
.org krq8a+0x20000
b rsx20
.org krq8b+0x20000
b rsx24
.org krq9a+0x20000
b rsx20
.org krq9b+0x20000
b rsx24
.org krq10a+0x20000
b rsx20
.org krq10b+0x20000
b rsx24
.org krq11a+0x20000
b rsx20
.org krq11b+0x20000
b rsx24
.org krq12a+0x20000
b rsx20
.org krq12b+0x20000
b rsx24
.org krqa0a+0x20000
b rsx20
.org krqa0b+0x20000
b rsx24
.org krqa1a+0x20000
b rsx20
.org krqa1b+0x20000
b rsx24
.org krqa2a+0x20000
b rsx20
.org krqa2b+0x20000
b rsx24
.org krqa3a+0x20000
b rsx20
.org krqa3b+0x20000
b rsx24
.org krqa4a+0x20000
b rsx20
.org krqa4b+0x20000
b rsx24
.org krqa5a+0x20000
b rsx20
.org krqa5b+0x20000
b rsx24
.org krqa6a+0x20000
b rsx20
.org krqa6b+0x20000
b rsx24
.org krqa7a+0x20000
b rsx20
.org krqa7b+0x20000
b rsx24
.org krqb0a+0x20000
b rsx20
.org krqb0b+0x20000
b rsx24
.org krqb1a+0x20000
b rsx20
.org krqb1b+0x20000
b rsx24
.org krqb2a+0x20000
b rsx20
.org krqb2b+0x20000
b rsx24
.org krqb3a+0x20000
b rsx20
.org krqb3b+0x20000
b rsx24
.org krqb4a+0x20000
b rsx20
.org krqb4b+0x20000
b rsx24
.org krqb5a+0x20000
b rsx20
.org krqb5b+0x20000
b rsx24
.org krqb6a+0x20000
b rsx20
.org krqb6b+0x20000
b rsx24
.org krqb7a+0x20000
b rsx20
.org krqb7b+0x20000
b rsx24
.org krqc0a+0x20000
b rsx20
.org krqc0b+0x20000
b rsx24
.org krqc1a+0x20000
b rsx20
.org krqc1b+0x20000
b rsx24
.org krqc2a+0x20000
b rsx20
.org krqc2b+0x20000
b rsx24
.org krqc3a+0x20000
b rsx20
.org krqc3b+0x20000
b rsx24
.org krqc4a+0x20000
b rsx20
.org krqc4b+0x20000
b rsx24
.org krqc5a+0x20000
b rsx20
.org krqc5b+0x20000
b rsx24
.org krqc6a+0x20000
b rsx20
.org krqc6b+0x20000
b rsx24
.org krqc7a+0x20000
b rsx20
.org krqc7b+0x20000
b rsx24
.org krqd0a+0x20000
b rsx20
.org krqd0b+0x20000
b rsx24
.org krqd1a+0x20000
b rsx20
.org krqd1b+0x20000
b rsx24
.org krqd2a+0x20000
b rsx20
.org krqd2b+0x20000
b rsx24
.org krqd3a+0x20000
b rsx20
.org krqd3b+0x20000
b rsx24
.org krqd4a+0x20000
b rsx20
.org krqd4b+0x20000
b rsx24
.org krqd5a+0x20000
b rsx20
.org krqd5b+0x20000
b rsx24
.org krqd6a+0x20000
b rsx20
.org krqd6b+0x20000
b rsx24
.org krqd7a+0x20000
b rsx20
.org krqd7b+0x20000
b rsx24
.org krqe0a+0x20000
b rsx20
.org krqe0b+0x20000
b rsx24
.org krqe1a+0x20000
b rsx20
.org krqe1b+0x20000
b rsx24
.org krqe2a+0x20000
b rsx20
.org krqe2b+0x20000
b rsx24
.org krqe3a+0x20000
b rsx20
.org krqe3b+0x20000
b rsx24
.org krqe4a+0x20000
b rsx20
.org krqe4b+0x20000
b rsx24
.org krqe5a+0x20000
b rsx20
.org krqe5b+0x20000
b rsx24
.org krqe6a+0x20000
b rsx20
.org krqe6b+0x20000
b rsx24
.org krqe7a+0x20000
b rsx20
.org krqe7b+0x20000
b rsx24
.org krqf0a+0x20000
b rsx20
.org krqf0b+0x20000
b rsx24
.org krqf1a+0x20000
b rsx20
.org krqf1b+0x20000
b rsx24
.org krqf2a+0x20000
b rsx20
.org krqf2b+0x20000
b rsx24
.org krqf3a+0x20000
b rsx20
.org krqf3b+0x20000
b rsx24
.org krqf4a+0x20000
b rsx20
.org krqf4b+0x20000
b rsx24
.org krqf5a+0x20000
b rsx20
.org krqf5b+0x20000
b rsx24
.org krqf6a+0x20000
b rsx20
.org krqf6b+0x20000
b rsx24
.org krqf7a+0x20000
b rsx20
.org krqf7b+0x20000
b rsx24
.org krqg0a+0x20000
b rsx20
.org krqg0b+0x20000
b rsx24
.org krqg1a+0x20000
b rsx20
.org krqg1b+0x20000
b rsx24
.org krqg2a+0x20000
b rsx20
.org krqg2b+0x20000
b rsx24
.org krqg3a+0x20000
b rsx20
.org krqg3b+0x20000
b rsx24
.org krqg4a+0x20000
b rsx20
.org krqg4b+0x20000
b rsx24
.org krqg5a+0x20000
b rsx20
.org krqg5b+0x20000
b rsx24
.org krqg6a+0x20000
b rsx20
.org krqg6b+0x20000
b rsx24
.org krqg7a+0x20000
b rsx20
.org krqg7b+0x20000
b rsx24
.org krq5+0x20000
b rsx20
.org krf0a+0x20000
b rsx20
.org krf0b+0x20000
b rsx24
.org krf1a+0x20000
b rsx20
.org krf1b+0x20000
b rsx24
.org krf2a+0x20000
b rsx20
.org krf2b+0x20000
b rsx24
.org krf3a+0x20000
b rsx20
.org krf3b+0x20000
b rsx24
.org krf4a+0x20000
b rsx20
.org krf4b+0x20000
b rsx24
.org krf5a+0x20000
b rsx20
.org krf5b+0x20000
b rsx24
.org krf6a+0x20000
b rsx20
.org krf6b+0x20000
b rsx24
.org krf7a+0x20000
b rsx20
.org krf7b+0x20000
b rsx24
.org krfa0a+0x20000
b rsx20
.org krfa0b+0x20000
b rsx24
.org krfa1a+0x20000
b rsx20
.org krfa1b+0x20000
b rsx24
.org krfa2a+0x20000
b rsx20
.org krfa2b+0x20000
b rsx24
.org krfa3a+0x20000
b rsx20
.org krfa3b+0x20000
b rsx24
.org krfa4a+0x20000
b rsx20
.org krfa4b+0x20000
b rsx24
.org krfa5a+0x20000
b rsx20
.org krfa5b+0x20000
b rsx24
.org krfa6a+0x20000
b rsx20
.org krfa6b+0x20000
b rsx24
.org krfa7a+0x20000
b rsx20
.org krfa7b+0x20000
b rsx24
.org krfb0a+0x20000
b rsx20
.org krfb0b+0x20000
b rsx24
.org krfb1a+0x20000
b rsx20
.org krfb1b+0x20000
b rsx24
.org krfb2a+0x20000
b rsx20
.org krfb2b+0x20000
b rsx24
.org krfb3a+0x20000
b rsx20
.org krfb3b+0x20000
b rsx24
.org krfb4a+0x20000
b rsx20
.org krfb4b+0x20000
b rsx24
.org krfb5a+0x20000
b rsx20
.org krfb5b+0x20000
b rsx24
.org krfb6a+0x20000
b rsx20
.org krfb6b+0x20000
b rsx24
.org krfb7a+0x20000
b rsx20
.org krfb7b+0x20000
b rsx24
.org krfc0a+0x20000
b rsx20
.org krfc0b+0x20000
b rsx24
.org krfc1a+0x20000
b rsx20
.org krfc1b+0x20000
b rsx24
.org krfc2a+0x20000
b rsx20
.org krfc2b+0x20000
b rsx24
.org krfc3a+0x20000
b rsx20
.org krfc3b+0x20000
b rsx24
.org krfc4a+0x20000
b rsx20
.org krfc4b+0x20000
b rsx24
.org krfc5a+0x20000
b rsx20
.org krfc5b+0x20000
b rsx24
.org krfc6a+0x20000
b rsx20
.org krfc6b+0x20000
b rsx24
.org krfc7a+0x20000
b rsx20
.org krfc7b+0x20000
b rsx24
.org krfd0a+0x20000
b rsx20
.org krfd0b+0x20000
b rsx24
.org krfd1a+0x20000
b rsx20
.org krfd1b+0x20000
b rsx24
.org krfd2a+0x20000
b rsx20
.org krfd2b+0x20000
b rsx24
.org krfd3a+0x20000
b rsx20
.org krfd3b+0x20000
b rsx24
.org krfd4a+0x20000
b rsx20
.org krfd4b+0x20000
b rsx24
.org krfd5a+0x20000
b rsx20
.org krfd5b+0x20000
b rsx24
.org krfd6a+0x20000
b rsx20
.org krfd6b+0x20000
b rsx24
.org krfd7a+0x20000
b rsx20
.org krfd7b+0x20000
b rsx24
.org krfe0a+0x20000
b rsx20
.org krfe0b+0x20000
b rsx24
.org krfe1a+0x20000
b rsx20
.org krfe1b+0x20000
b rsx24
.org krfe2a+0x20000
b rsx20
.org krfe2b+0x20000
b rsx24
.org krfe3a+0x20000
b rsx20
.org krfe3b+0x20000
b rsx24
.org krfe4a+0x20000
b rsx20
.org krfe4b+0x20000
b rsx24
.org krfe5a+0x20000
b rsx20
.org krfe5b+0x20000
b rsx24
.org krfe6a+0x20000
b rsx20
.org krfe6b+0x20000
b rsx24
.org krfe7a+0x20000
b rsx20
.org krfe7b+0x20000
b rsx24
.org krff0a+0x20000
b rsx20
.org krff0b+0x20000
b rsx24
.org krff1a+0x20000
b rsx20
.org krff1b+0x20000
b rsx24
.org krff2a+0x20000
b rsx20
.org krff2b+0x20000
b rsx24
.org krff3a+0x20000
b rsx20
.org krff3b+0x20000
b rsx24
.org krff4a+0x20000
b rsx20
.org krff4b+0x20000
b rsx24
.org krff5a+0x20000
b rsx20
.org krff5b+0x20000
b rsx24
.org krff6a+0x20000
b rsx20
.org krff6b+0x20000
b rsx24
.org krff7a+0x20000
b rsx20
.org krff7b+0x20000
b rsx24
.org krfg0a+0x20000
b rsx20
.org krfg0b+0x20000
b rsx24
.org krfg1a+0x20000
b rsx20
.org krfg1b+0x20000
b rsx24
.org krfg2a+0x20000
b rsx20
.org krfg2b+0x20000
b rsx24
.org krfg3a+0x20000
b rsx20
.org krfg3b+0x20000
b rsx24
.org krfg4a+0x20000
b rsx20
.org krfg4b+0x20000
b rsx24
.org krfg5a+0x20000
b rsx20
.org krfg5b+0x20000
b rsx24
.org krfg6a+0x20000
b rsx20
.org krfg6b+0x20000
b rsx24
.org krfg7a+0x20000
b rsx20
.org krfg7b+0x20000
b rsx24
.org kx+0x20000
b xst


home_join_reset:
ldp x13,x25,[x9,#-144]
ret

relay_source_init:
movn x0,#1023
movk x0,#0,lsl #48
mov x1,#0x7fff0000000000
add x1,x1,#64
mov x2,#64
1:
str x1,[x0],#8
add x1,x1,#64
subs x2,x2,#1
b.ne 1b
ret
home_source_init:
movn x0,#4095
movk x0,#127,lsl #48
mov x1,#0x7fff0000000000
add x1,x1,#64
mov x2,#128
2:
str x1,[x0],#8
add x1,x1,#64
subs x2,x2,#1
b.ne 2b
ret

.globl relay_blocks_init_split
relay_blocks_init_split:
mov x8,x30
mov x1,#0x200000000000
mov x3,#12288
mov x0,#43690
movk x0,#43690,lsl #16
movk x0,#682,lsl #32
mov x2,#0xff80000000000000
add x2,x2,#32
bl init3
mov x0,#21846
movk x0,#21845,lsl #16
movk x0,#341,lsl #32
mov x2,#57376
movk x2,#65535,lsl #16
movk x2,#65279,lsl #32
movk x2,#63,lsl #48
bl init3
mov x30,x8
ret
.globl relay_fin_split
relay_fin_split:
mov x18,x30
mov x7,#0x200000000000
mov x3,#0x600000000000
mov x1,#43690
movk x1,#43690,lsl #16
movk x1,#682,lsl #32
mov x8,#0xff80000000000000
add x8,x8,#32
bl fin
mov x1,#21846
movk x1,#21845,lsl #16
movk x1,#341,lsl #32
mov x8,#57376
movk x8,#65535,lsl #16
movk x8,#65279,lsl #32
movk x8,#63,lsl #48
bl fin
mov x30,x18
ret

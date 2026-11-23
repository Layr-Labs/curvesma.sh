.macro GENERIC_COMMIT cmd
stp x0,x30,[sp,#-96]!
stp x16,x17,[sp,#16]
mov x18,#\cmd
str x18,[sp,#32]
ldp x16,x17,[x0]
stp x16,x17,[sp,#48]
ldp x16,x17,[x0,#16]
stp x16,x17,[sp,#64]
ldp x16,x17,[x0,#32]
stp x16,x17,[sp,#80]
add x0,sp,#16
bl generic_send
ldp x0,x30,[sp],#96
.endm
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
mov x11,#0x80000000000
add x11,x11,x7
mov x9,#0x180000000000
ldr x5,[x9,#8]
str xzr,[x9]
dsb sy
str x5,[x11]
dsb sy
gfin_wait:
ldr x5,[x9]
cbnz x5,gfin_done
wfe
b gfin_wait
gfin_done:
str xzr,[x10,#16]
add x9,x10,#32
str x9,[x7],#8
subs x1,x1,#1
b.ne 1b
dsb sy
ret
.globl fin
fin:
1:
ldr x9,[x7]
ldr x5,[x3],#8
sub x10,x9,x8
cmp x5,x6
csel x10,x10,xzr,lt
orr x10,x10,#1
mov x5,#0x8010
str xzr,[x5]
dsb sy
mov x9,#0x200000000000
add x9,x9,x7
str x10,[x9]
dsb sy
rhfin_wait:
ldr x9,[x5]
cbnz x9,rhfin_done
wfe
b rhfin_wait
rhfin_done:
add x7,x7,#8
add x8,x8,#0x4000
subs x1,x1,#1
b.ne 1b
ret
.macro H2B a,b,c,d,l
bfi x10,\a,#3,#42
ldadd x17,x11,[x10]
stp \a,\b,[x11]
stp \c,\d,[x11,#16]
tbnz x11,#13,hf\l
hn\l:
.pushsection .text.cold
hf\l:
sub x9,x11,#8192
add x9,x9,#32
str x9,[x10]
orr x9,x10,#0x400000000000
ldadd x18,x11,[x9]
cmp x11,x12
b.lt hf8\l
str x12,[x9]
b hn\l
hf8\l:
mov x9,#0x8008
str xzr,[x9]
dsb sy
mov x9,#0x200000000000
add x9,x9,x10
orr x11,x18,#3
str x11,[x9]
dsb sy
mov x9,#0x8008
hf_wait\l:
ldr x11,[x9]
cbnz x11,hn\l
wfe
b hf_wait\l
.popsection
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
mov x0,#0xc00000000000
add x0,x0,#32
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
sub x2,x20,#32
H2R1 qfirst
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
ldp x20,x21,[x0]
ldp x22,x23,[x0,#16]
H2B x20,x21,x22,x23,ffirst
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
ldp x20,x21,[x7]
ldp x22,x23,[x7,#16]
H2B x20,x21,x22,x23,flast
sub sp,sp,#448
ldr x0,[sp,#512]
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
ldp x16,x17,[x24]
ldp x18,x19,[x24,#16]
add x24,x24,#32
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
ld64b x16,[x24]
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
ld64b x16,[x24]
HCHGI x17,x16,0
HCHGI x21,x20,32
add x24,x24,#64
.endm
.macro DLVHPRZ w
ld64b x16,[x24]
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
csel x3,x17,x24,lt
add x16,x24,x16
sub x7,x16,x3
mov x24,x3
mov x1,x3
SEGL DLVHS,DLVHP,2
staddl x28,[x29]
ret
6:
cbz x24,dlvh_zero_full
DLVHPR x0
dlvh_after_first:
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
DLVHPR x7
DLVHPR x0
DLVHPR x1
DLVHPR x2
DLVHPR x3
DLVHPR x4
DLVHPR x5
DLVHPR x6
ld64b x16,[x24]
HCHGI x17,x16,0
HCHGI x21,x20,32
staddl x28,[x29]
ret
dlvh_zero_full:
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
str x13,[x0,#40]
CONST x16,16
ldp x16,x17,[x16,#24]
lsr x16,x16,#42
and x17,x17,#0x3ffff
orr x16,x16,x17,lsl #22
mov x17,#0x40000
GENERIC_COMMIT 1
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
str xzr,[x0,#40]
ldp x16,x17,[x0]
GENERIC_COMMIT 1
LINE64END
.globl fwdq1
fwdq1:
cmp x1,x6
b.hs 2f
1:
ldp x2,x3,[x1,#8]
ubfx x4,x2,#3,#25
ubfx x5,x3,#9,#33
orr x4,x4,x5,lsl #25
stp x4,xzr,[x0]
lsl x5,x1,#8
stp x3,x5,[x0,#16]
str x7,[x0,#32]
str x2,[x0,#40]
ldp x16,x17,[x0]
GENERIC_COMMIT 1
LINE64END
.globl sendr
sendr:
1:
ldp x2,x3,[x1],#16
ldp x4,x5,[x1],#16
stp x7,x8,[x0]
stp x2,x3,[x0,#16]
str x4,[x0,#32]
str x5,[x0,#40]
ldr x16,[x0]
lsr x16,x16,#3
mov x17,#0x40000
GENERIC_COMMIT 1
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
str x11,[x0,#40]
ldr x16,[x0]
and x16,x16,#0x1fffffffffff
lsr x16,x16,#3
CONST x17,16
ldr x17,[x17,#24]
orr x16,x16,x17,lsl #42
lsr x17,x17,#22
orr x17,x17,#0x80000
GENERIC_COMMIT 0
LINE64END
// Retained historical Y-ring bodies ydlv0/ysnd/ydlv/ydlve: inactive in
// the current ystage schedule; their old completing stores define no service.
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
stp x14,x13,[x5]
stp x2,x26,[x5,#16]
tbnz x5,#12,xemit_fl
xemit_n:
add x26,x26,#4
add x1,x1,#512
subs x3,x3,#1
b.ne 1b
ret
xemit_fl:
sub x16,x5,#4064
str x16,[x0]
ldr x17,[x16,#-16]
add x17,x17,#1
mov x5,#4096
stp x5,x17,[x16,#-24]
mov x5,#0x80000000000
add x5,x5,x0
mov x16,#0x180000000000
ldr x17,[x16,#8]
str xzr,[x16]
dsb sy
str x17,[x5]
dsb sy
xfull_wait:
ldr x17,[x16]
cbnz x17,xemit_n
wfe
b xfull_wait
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
// Historical ring dispatcher; no current callback table binds this entry.
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
relay_source_init:
movn x0,#1023
movk x0,#0,lsl #48
mov x1,#0xc00000000000
add x1,x1,#64
sub x3,x1,#32
mov x2,#64
1:
str x1,[x0],#8
add x1,x1,#64
subs x2,x2,#1
b.ne 1b
str x3,[x0]
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

// Card54 source-slot initialization; x0=table base,x1=buffer base,x4=count.
relay_slot_tables_init:
mov x5,#0x10000000000
1:
add x3,x1,#32
str x3,[x0,#512]
str x5,[x0,#520]
add x5,x5,#8
add x1,x1,#64
mov x2,#64
2:
str x1,[x0],#8
add x1,x1,#64
subs x2,x2,#1
b.ne 2b
add x0,x0,#512
add x1,x1,#4032
subs x4,x4,#1
b.ne 1b
ret
// Isolated direct-concatenation binding; x4 is an already decoded member0..1023.
relay_slot_bind:
mov x6,#0x1000000
movk x6,#0xc000,lsl #32
bfi x6,x4,#10,#10
mov sp,x6
ldr x0,[sp,#512]
ret

// Card54 scalar R receive loop; finite probe stops only after the back edge.
.globl bulk_r_receive
bulk_r_receive:
mov x9,#0x8000
bulk_r_wait:
ldr x0,[x9]
cbnz x0,bulk_r_arrived
wfe
b bulk_r_wait
bulk_r_arrived:
str xzr,[x9]
dsb sy
sub x0,x0,#1
mov sp,x0
ldr x0,[sp,#512]
bl h2s
dsb sy
ldr x9,[sp,#520]
mov x0,#1
str x0,[x9]
dsb sy
b bulk_r_receive

// Card54 H48: physicalBASE is cursor, SINK=LIMIT, no readback; DSB remains.
.globl bulk_h_receive
bulk_h_receive:
mov x1,#0x7ffffffc000000
bulk_h_wait:
ldr x16,[x1]
cbnz x16,bulk_h_arrived
wfe
b bulk_h_wait
bulk_h_arrived:
dsb sy
and x17,x16,#0xffffffffffffffe0
str x17,[x0,#-16]
bl dlvh
dsb sy
mov x1,#0x7ffffffc000000
ldr x16,[x1]
str xzr,[x1]
tbz x16,#1,bulk_h_final
ldr x16,[x0,#-24]
add x16,x16,#8192
str x16,[x0,#-24]
ldp x16,x17,[x1,#16]
add x16,x16,#8192
cmp x16,x17
csel x16,x16,x17,lo
str x16,[x1,#16]
dsb sy
bulk_h_ack:
mov x16,#1
str x16,[x1,#256]
dsb sy
b bulk_h_receive
bulk_h_final:
and x17,x16,#0xffffffffffffffe0
ldr x18,[x0,#-24]
add x17,x18,x17
ldr x18,[x0,#-8]
str x17,[x18]
dsb sy
b bulk_h_final_ack
bulk_h_final_ack:
mov x16,#1
str x16,[x1,#264]
dsb sy
b bulk_h_receive


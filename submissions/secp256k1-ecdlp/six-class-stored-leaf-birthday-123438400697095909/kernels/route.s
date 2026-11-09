// Transport and home grouping (proof.md 2.1, 2.3, 4); a record is K2, K3, w2, w3; x30, sp preserved.
.text
.p2align 2
// A one-record message (proof.md 2.1): x0 its 48 bytes; the four data words go to its address.
.globl dlv
dlv:
ldp x1,x2,[x0]
ldp x3,x4,[x0,#16]
ldp x5,x6,[x0,#32]
stp x3,x4,[x2]
stp x5,x6,[x2,#16]
ret
// beta tables (proof.md 3.1): entry t of slot k at x0 + 256t + 32k; x1 the 7 slot deltas beta 2^(43i) mod p,
// i = (0,1,3,2,4,5,3)[k]; x2 entries per slot, slot 5: x3; entry t+1 = entry t + delta mod p.
.globl tbld
tbld:
mov x7,#977
movk x7,#1,lsl #32
mov x8,#7
1:
ldp x9,x10,[x1],#16
ldp x11,x12,[x1],#16
mov x13,xzr
mov x14,xzr
mov x15,xzr
mov x16,xzr
mov x17,x0
cmp x8,#2
csel x18,x3,x2,eq
2:
stp x13,x14,[x17]
stp x15,x16,[x17,#16]
add x17,x17,#256
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
// RC table: x0 entries (>= 1), x1 table, x2 slot base; every entry gets the base.
.globl tinit
tinit:
1:
str x2,[x1],#8
subs x0,x0,#1
b.ne 1b
ret
// E3 (proof.md 4): x0 entries, x1 a final-cursor table (32-byte stride), x2 E3.
.globl ecopy
ecopy:
1:
ldr x3,[x1],#32
str x3,[x2],#8
subs x0,x0,#1
b.ne 1b
ret
// Final blocks and cursors: x0 message block, x7 LC table, x3 RC table, x8 block 0's slot 1, x1 destinations (>= 1),
// x2 this sender's entry at each destination, x6 its slot end, x11 p. Destination d has L = (LC - slot 1)/32 records in its
// block and RC cursor c: if L > 0 and c < end the header (8d, c) sends the block; d receives min(c + 32L, end).
.globl fin
fin:
mov x4,#0
1:
ldr x9,[x7],#8
ldr x5,[x3],#8
subs x10,x9,x8
b.eq 2f
cmp x5,x6
b.hs 2f
stp x10,x11,[x8,#-16]
stp x4,x5,[x8,#-32]
2:
add x5,x5,x10
cmp x5,x6
csel x5,x5,x6,lo
stp x4,x2,[x0]
stp x5,xzr,[x0,#16]
stp xzr,xzr,[x0,#32]
add x4,x4,#8
add x8,x8,#16384
subs x1,x1,#1
b.ne 1b
ret
// Hop-2 stream on block arrival (proof.md 4): x1 the first record, x12 the hop-2 slot end, x19 p; d2 by BFI into x10;
// slot 1 and a last even slot by two LDPs (H2R1), pairs by LD64B (H2P); a full block unrolled, else 1, 2, then 4 pairs.
.macro H2B a,b,c,d
bfi x10,\d,#3,#36
ldadd x17,x11,[x10]
stp \a,\b,[x11]
stp \c,\d,[x11,#16]
tbz x11,#13,9f
sub x9,x11,#4080
sub x9,x9,#4080
str x9,[x10]
orr x4,x10,#0x80000000000
ldadd x18,x11,[x4]
cmp x11,x12
b.lo 8f
str x12,[x4]
b 9f
8:
stp x18,x19,[x9,#-16]
stp x10,x11,[x9,#-32]
9:
.endm
.macro H2R1
ldp x20,x21,[x1],#16
ldp x22,x23,[x1],#16
H2B x20,x21,x22,x23
.endm
.macro H2P
ld64b x20,[x1]
add x1,x1,#64
H2B x20,x21,x22,x23
H2B x24,x25,x26,x27
.endm
.macro H2P8
H2P
H2P
H2P
H2P
H2P
H2P
H2P
H2P
.endm
.macro H2P64
H2P8
H2P8
H2P8
H2P8
H2P8
H2P8
H2P8
H2P8
.endm
// The hop-1 delivery processor (proof.md 2.1): h2s0 once at set-up (x12 the slot end, kept), then h2s per block (a loop,
// B for RET); x10, x18, x3 persist; x1 reset each run; p loaded after the window wait from the node word at 2^48 - 168.
.globl h2s0
h2s0:
mov x10,#0x8000000000
mov x18,#8192
movn x3,#0xa7
movk x3,#0,lsl #48
ret
.globl h2s
h2s:
movn x1,#0x3fdf
ldr x16,[x1,#-16]
ldr x19,[x3]
mov x17,#32
tbnz x16,#13,7f
cbz x16,6f
sub x16,x16,#32
add x2,x1,x16
H2R1
tbz x16,#6,3f
H2P
3:
tbz x16,#7,4f
H2P
H2P
4:
cmp x1,x2
b.hs 5f
1:
H2P
H2P
H2P
H2P
cmp x1,x2
b.lo 1b
5:
tbz x16,#5,6f
H2R1
6:
ret
// A full block (count 8192): slot 1, 127 pairs, slot 256, unrolled.
7:
H2R1
H2P64
H2P8
H2P8
H2P8
H2P8
H2P8
H2P8
H2P8
H2P
H2P
H2P
H2P
H2P
H2P
H2P
H2R1
ret
// Hash-chain join (proof.md 4): K2 \u, K3 \v (bits 44..63 compared); heads at 2^61 (BFI of K2's low 58 bits); entry
// E = (A>>1) + C (x13; x14 = C = 2^62 - (R2BASE>>1)); x15 = -2C gives A = 2E - 2C on a K2 match; \a by \l (list path). Uncapped walk.
.macro LNONE
.endm
.macro LADD
add x12,x2,#32
.endm
.macro LSUB
sub x12,x2,#32
.endm
.macro HCHGI u,v,a,l
bfi x10,\u,#3,#58
swp x13,x11,[x10]
stp x11,\u,[x13],#16
cbz x11,8f
9:
ldp x5,x6,[x11]
cmp x6,\u
b.eq 10f
12:
cbz x5,8f
ldp x11,x6,[x5]
cmp x6,\u
b.eq 11f
13:
cbnz x11,9b
b 8f
10:
add x6,x15,x11,lsl #1
ldr x12,[x6,#8]
eor x12,x12,\v
tst x12,#0xfffff00000000000
b.ne 12b
cmp x3,x4
b.hs 8f
\l
stp x6,\a,[x3],#16
b 8f
11:
add x6,x15,x5,lsl #1
ldr x12,[x6,#8]
eor x12,x12,\v
tst x12,#0xfffff00000000000
b.ne 13b
cmp x3,x4
b.hs 8f
\l
stp x6,\a,[x3],#16
8:
.endm
// Fused hop-2 arrival (proof.md 4): x0 the window (header c, count, p; F = *p); joins [max(F, c), c + count); join
// state at 2^48 - 160; x1 window, x2 R2 at the next slot; pairs by LD64B/ST64B x16..x23 (DLVP), a slot alone by LDP/STP.
.macro DLVS
ldp x16,x17,[x1],#32
ldp x18,x19,[x1,#-16]
stp x18,x19,[x2,#16]
stp x16,x17,[x2],#32
.endm
.macro DLVP
ld64b x16,[x1]
add x1,x1,#64
st64b x16,[x2]
add x2,x2,#64
.endm
.macro DLVHS
DLVS
HCHGI x16,x17,x12,LSUB
.endm
.macro DLVHP
ld64b x16,[x1]
add x1,x1,#64
st64b x16,[x2]
HCHGI x16,x17,x2,LNONE
HCHGI x20,x21,x12,LADD
add x2,x2,#64
.endm
// A segment of x7 bytes (32 per slot), x8 = its window end - 32: if it is non-empty, an odd first slot alone, then
// pairs 1, 2, then 4 at a time, then an unpaired last slot.
.macro SEG s,p,l
cbz x7,sge\l
tbz x1,#5,sga\l
\s
sub x7,x7,#32
sga\l:
tbz x7,#6,sgb\l
\p
sgb\l:
tbz x7,#7,sgc\l
\p
\p
sgc\l:
cmp x1,x8
b.hs sgd\l
sgl\l:
\p
\p
\p
\p
cmp x1,x8
b.lo sgl\l
sgd\l:
tbz x7,#5,sge\l
\s
sge\l:
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
.macro DLVHP64
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
.endm
.globl dlvh
dlvh:
ldp x24,x26,[x0,#8]
ldr x13,[x0,#24]
ldr x27,[x13]
mov x9,#0x1000000000000
ldp x3,x4,[x9,#-160]
ldp x14,x15,[x9,#-144]
mov x10,#0x2000000000000000
add x1,x0,#32
mov x2,x24
cmp x24,x27
b.lo 7f
tbnz x26,#13,6f
7:
csel x27,x27,x24,lo
sub x7,x27,x24
add x8,x0,x7
SEG DLVS,DLVP,1
add x25,x24,x26
sub x7,x25,x27
add x8,x0,x26
add x13,x14,x2,lsr #1
SEG DLVHS,DLVHP,2
add x8,x0,#8192
mov x7,#8192
sub x7,x7,x26
SEG DLVS,DLVP,3
str x3,[x9,#-160]
ret
6:
add x13,x14,x2,lsr #1
DLVHS
DLVHP64
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP8
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
DLVHP
DLVHS
str x3,[x9,#-160]
ret
// Zero x0 words (a multiple of 16384) at x1 (256-byte aligned): 512 DC ZVA + ADD #256, SUBS, B.NE per pass (1026).
.macro Z8
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
dc zva,x1
add x1,x1,#256
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
.macro Z512
Z64
Z64
Z64
Z64
Z64
Z64
Z64
Z64
.endm
.globl zero16384
zero16384:
1:
Z512
subs x0,x0,#16384
b.ne 1b
ret
// LC table (init3): x0 entries, x1 table, x2 first value, x3 stride; entry t = first + t*stride.
.globl init3
init3:
1:
str x2,[x1],#8
add x2,x2,x3
subs x0,x0,#1
b.ne 1b
ret
// tblx (proof.md 3.1): x23 lines (W7: pair 4, 3; delta 2, 2 by two calls) per row t < x3 at x0 + x24 t; line k's two entries (s t + init) mod p,
// steps s at x1 (two per line), line 0's inits at x2; one ST64B per line; a + s mod p as a + (s + C0), C0 off on carry.
.macro PSTEP
ldp x14,x15,[x1],#16
ldp x16,x17,[x1],#16
ldp x18,x19,[x1],#16
ldp x20,x21,[x1],#16
adds x14,x14,x5
adcs x15,x15,xzr
adcs x16,x16,xzr
adcs x17,x17,xzr
adds x18,x18,x5
adcs x19,x19,xzr
adcs x20,x20,xzr
adcs x21,x21,xzr
.endm
.macro AMOD a0,a1,a2,a3,s0,s1,s2,s3
adds \a0,\a0,\s0
adcs \a1,\a1,\s1
adcs \a2,\a2,\s2
adcs \a3,\a3,\s3
csel x22,xzr,x5,hs
subs \a0,\a0,x22
sbcs \a1,\a1,xzr
sbcs \a2,\a2,xzr
sbcs \a3,\a3,xzr
.endm
.macro PLOOP
mov x4,x3
1:
st64b x6,[x0]
add x0,x0,x24
AMOD x6,x7,x8,x9,x14,x15,x16,x17
AMOD x10,x11,x12,x13,x18,x19,x20,x21
subs x4,x4,#1
b.ne 1b
sub x0,x0,x25
add x0,x0,#64
.endm
.globl tblx
tblx:
mov x5,#977
movk x5,#1,lsl #32
mul x25,x3,x24
PSTEP
ldp x6,x7,[x2]
ldp x8,x9,[x2,#16]
ldp x10,x11,[x2,#32]
ldp x12,x13,[x2,#48]
PLOOP
3:
subs x23,x23,#1
b.eq 4f
PSTEP
mov x6,xzr
mov x7,xzr
mov x8,xzr
mov x9,xzr
mov x10,xzr
mov x11,xzr
mov x12,xzr
mov x13,xzr
PLOOP
b 3b
4:
ret

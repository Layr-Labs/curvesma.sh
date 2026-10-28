	.macro SQ a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r
	mul \a,\b,\b
	umulh \c,\b,\b
	adds \d,\c,\d,lsl #1
	mul \c,\e,\e
	adcs \f,\f,\c
	umulh \c,\e,\e
	adcs \g,\g,\c
	mul \c,\h,\h
	adcs \i,\i,\c
	umulh \c,\h,\h
	adcs \j,\j,\c
	mul \c,\k,\k
	adcs \l,\l,\c
	umulh \c,\k,\k
	adc \m,\m,\c
	mul \c,\i,\n
	umulh \o,\i,\n
	mul \i,\j,\n
	adds \i,\i,\o
	umulh \o,\j,\n
	mul \j,\l,\n
	adcs \j,\j,\o
	umulh \o,\l,\n
	mul \l,\m,\n
	adcs \l,\l,\o
	umulh \o,\m,\n
	adc \o,\o,xzr
	adds \a,\a,\c
	adcs \d,\d,\i
	adcs \f,\f,\j
	adcs \g,\g,\l
	adc \o,\o,xzr
	ldp \l,\j,[\p,#-32]
	ldp \i,\c,[\p,#-16]
	subs \a,\a,\l
	sbcs \d,\d,\j
	sbcs \f,\f,\i
	sbcs \g,\g,\c
	adc \o,\o,xzr
	ldp \m,\q,[\r,#192]
	subs \a,\a,\m
	sbcs \d,\d,\q
	ldp \m,\q,[\r,#208]
	sbcs \f,\f,\m
	sbcs \g,\g,\q
	adc \o,\o,xzr
	mul \q,\o,\n
	umulh \m,\o,\n
	adds \a,\a,\q
	adcs \d,\d,\m
	.endm
	.macro FM a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u
	mul \a,\b,\c
	umulh \d,\b,\c
	mul \e,\b,\f
	adds \d,\d,\e
	umulh \e,\b,\f
	mul \g,\b,\h
	adcs \e,\e,\g
	umulh \g,\b,\h
	mul \i,\b,\j
	adcs \g,\g,\i
	umulh \i,\b,\j
	mul \k,\l,\j
	adcs \i,\i,\k
	umulh \k,\l,\j
	adc \k,\k,xzr
	mul \m,\l,\c
	adds \d,\d,\m
	umulh \m,\l,\c
	adcs \e,\e,\m
	umulh \m,\n,\c
	adcs \g,\g,\m
	umulh \m,\o,\c
	adcs \i,\i,\m
	umulh \m,\n,\h
	adc \m,\m,xzr
	mul \p,\n,\c
	adds \e,\e,\p
	mul \p,\o,\c
	adcs \g,\g,\p
	mul \p,\n,\h
	adcs \i,\i,\p
	adcs \k,\k,\m
	umulh \m,\n,\j
	adc \m,\m,xzr
	mul \p,\l,\f
	adds \e,\e,\p
	umulh \p,\l,\f
	adcs \g,\g,\p
	umulh \p,\l,\h
	adcs \i,\i,\p
	umulh \p,\o,\f
	adcs \k,\k,\p
	umulh \p,\o,\h
	adc \p,\p,xzr
	mul \q,\l,\h
	adds \g,\g,\q
	umulh \q,\n,\f
	adcs \i,\i,\q
	mul \q,\n,\j
	adcs \k,\k,\q
	adcs \m,\m,\p
	umulh \p,\o,\j
	adc \p,\p,xzr
	mul \q,\n,\f
	adds \g,\g,\q
	mul \q,\o,\f
	adcs \i,\i,\q
	mul \q,\o,\h
	adcs \k,\k,\q
	mul \q,\o,\j
	adcs \m,\m,\q
	adc \p,\p,xzr
	mul \q,\i,\r
	umulh \s,\i,\r
	mul \i,\k,\r
	adds \i,\i,\s
	umulh \s,\k,\r
	mul \k,\m,\r
	adcs \k,\k,\s
	umulh \s,\m,\r
	mul \m,\p,\r
	adcs \m,\m,\s
	umulh \s,\p,\r
	adc \s,\s,xzr
	adds \a,\a,\q
	adcs \d,\d,\i
	adcs \t,\e,\k
	adcs \u,\g,\m
	adc \s,\s,xzr
	.endm
	.macro RCHG
	str x0,[sp,#-16]!
	ldr x0,[sp,#16]
	sub x0,x0,#20
	str x0,[sp,#16]
	ldr x0,[sp],#16
	.endm
	.macro W1 k,p
	ldp x10,x11,[\p,#0]
	ldp x12,x13,[\p,#16]
	ldp x16,x17,[x1,#16]
	ldp x14,x15,[x1],#64
	subs x10,x10,x14
	sbcs x11,x11,x15
	sbcs x12,x12,x16
	sbcs x13,x13,x17
	csel x17, x4, xzr, lo
	subs x10,x10,x17
	b.lo .Lr1_0_\k
.Lb1_0_\k:
	stp x5,x6,[x3,#32]
	stp x7,x8,[x3,#48]
	stp x12,x13,[x3,#16]
	stp x10,x11,[x3],#64
	FM x9,x10,x5,x17,x16,x6,x15,x7,x14,x8,x10,x11,x19,x12,x13,x20,x21,x4,x13,x7,x8
	mul x19,x13,x4
	umulh x10,x13,x4
	adds x5,x9,x19
	adcs x6,x17,x10
	b.cs .Lr1_1_\k
.Lb1_1_\k:
	.endm
	.macro S1 k
.Lr1_0_\k:
	RCHG
	sbcs x11,x11,xzr
	sbcs x12,x12,xzr
	sbc x13,x13,xzr
	b .Lb1_0_\k
.Lr1_1_\k:
	RCHG
	adcs x7,x7,xzr
	adc x8,x8,xzr
	b .Lb1_1_\k
	.endm
	.macro P1 i,j
	ldp x9,x26,[x2],#16
	W1 \i,x9
	W1 \j,x26
	.endm
	.text
	.p2align 2
	.globl pass1
	.type pass1,%function
pass1:
	stp x19,x20,[sp,#-96]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	mov x4,#977
	movk x4,#1,lsl #32
	ldr x9,[x2],#8
	ldp x5,x6,[x9,#0]
	ldp x7,x8,[x9,#16]
	ldp x12,x13,[x1,#16]
	ldp x10,x11,[x1],#64
	subs x5,x5,x10
	sbcs x6,x6,x11
	sbcs x7,x7,x12
	sbcs x8,x8,x13
	csel x13, x4, xzr, lo
	subs x5,x5,x13
	b.lo .Lr1_0_w0
.Lb1_0_w0:
	add x3,x3,#64
	sub x0,x0,#1
	lsr x29,x0,#7
	and x0,x0,#127
	cbz x0, 3f
2:
	ldr x9,[x2],#8
	W1 128,x9
	subs x0,x0,#1
	b.ne 2b
3:
	cbz x29, 4f
1:
	P1 0,1
	P1 2,3
	P1 4,5
	P1 6,7
	P1 8,9
	P1 10,11
	P1 12,13
	P1 14,15
	P1 16,17
	P1 18,19
	P1 20,21
	P1 22,23
	P1 24,25
	P1 26,27
	P1 28,29
	P1 30,31
	P1 32,33
	P1 34,35
	P1 36,37
	P1 38,39
	P1 40,41
	P1 42,43
	P1 44,45
	P1 46,47
	P1 48,49
	P1 50,51
	P1 52,53
	P1 54,55
	P1 56,57
	P1 58,59
	P1 60,61
	P1 62,63
	P1 64,65
	P1 66,67
	P1 68,69
	P1 70,71
	P1 72,73
	P1 74,75
	P1 76,77
	P1 78,79
	P1 80,81
	P1 82,83
	P1 84,85
	P1 86,87
	P1 88,89
	P1 90,91
	P1 92,93
	P1 94,95
	P1 96,97
	P1 98,99
	P1 100,101
	P1 102,103
	P1 104,105
	P1 106,107
	P1 108,109
	P1 110,111
	P1 112,113
	P1 114,115
	P1 116,117
	P1 118,119
	P1 120,121
	P1 122,123
	P1 124,125
	P1 126,127
	subs x29,x29,#1
	b.ne 1b
4:
	stp x5,x6,[x3,#32]
	stp x7,x8,[x3,#48]
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	ret
	.macro S8 a,b,c,d,e,f,g,h
	S1 \a
	S1 \b
	S1 \c
	S1 \d
	S1 \e
	S1 \f
	S1 \g
	S1 \h
	.endm
pass1_rare:
	S8 0,1,2,3,4,5,6,7
	S8 8,9,10,11,12,13,14,15
	S8 16,17,18,19,20,21,22,23
	S8 24,25,26,27,28,29,30,31
	S8 32,33,34,35,36,37,38,39
	S8 40,41,42,43,44,45,46,47
	S8 48,49,50,51,52,53,54,55
	S8 56,57,58,59,60,61,62,63
	S8 64,65,66,67,68,69,70,71
	S8 72,73,74,75,76,77,78,79
	S8 80,81,82,83,84,85,86,87
	S8 88,89,90,91,92,93,94,95
	S8 96,97,98,99,100,101,102,103
	S8 104,105,106,107,108,109,110,111
	S8 112,113,114,115,116,117,118,119
	S8 120,121,122,123,124,125,126,127
	S1 128
.Lr1_0_w0:
	RCHG
	sbcs x6,x6,xzr
	sbcs x7,x7,xzr
	sbc x8,x8,xzr
	b .Lb1_0_w0
	.macro W2 k,p
	ldp x6,x7,[x3,#32]
	ldp x8,x9,[x3,#48]
	FM x10,x25,x6,x11,x12,x7,x13,x8,x14,x9,x15,x26,x16,x27,x28,x17,x6,x4,x9,x12,x13
	mul x16,x9,x4
	umulh x15,x9,x4
	adds x10,x10,x16
	adcs x11,x11,x15
	b.cs .Lr2_0_\k
.Lb2_0_\k:
	ldp x9,x14,[x3,#16]
	ldp x15,x16,[x3],#-64
	FM x6,x15,x25,x17,x8,x26,x7,x27,x19,x28,x15,x16,x20,x9,x14,x21,x22,x4,x14,x27,x28
	mul x20,x14,x4
	umulh x15,x14,x4
	adds x25,x6,x20
	adcs x26,x17,x15
	b.cs .Lr2_1_\k
.Lb2_1_\k:
	ldp x7,x8,[\p,#32]
	ldp x17,x6,[\p,#48]
	ldp x15,x20,[x1,#0]
	ldp x14,x19,[x1,#16]
	subs x7,x7,x15
	sbcs x8,x8,x20
	sbcs x17,x17,x14
	sbcs x6,x6,x19
	csel x22, x4, xzr, lo
	subs x7,x7,x22
	b.lo .Lr2_2_\k
.Lb2_2_\k:
	FM x22,x7,x10,x21,x9,x11,x16,x12,x23,x13,x7,x8,x24,x17,x6,x29,x10,x4,x13,x9,x16
	mul x24,x13,x4
	umulh x7,x13,x4
	adds x22,x22,x24
	adcs x21,x21,x7
	b.cs .Lr2_3_\k
.Lb2_3_\k:
	mul x7,x22,x16
	mul x24,x21,x9
	adds x7,x7,x24
	umulh x24,x22,x16
	umulh x13,x21,x9
	adcs x24,x24,x13
	umulh x13,x21,x16
	adc x13,x13,xzr
	umulh x23,x22,x21
	mul x10,x22,x9
	adds x23,x23,x10
	umulh x10,x22,x9
	adcs x7,x7,x10
	mul x10,x21,x16
	adcs x24,x24,x10
	mul x10,x9,x16
	adcs x13,x13,x10
	umulh x10,x9,x16
	adc x10,x10,xzr
	mul x29,x22,x21
	lsr x6,x10,#63
	extr x10,x10,x13,#63
	extr x13,x13,x24,#63
	extr x24,x24,x7,#63
	extr x7,x7,x23,#63
	extr x23,x23,x29,#63
	SQ x12,x22,x11,x29,x21,x23,x7,x9,x24,x13,x16,x10,x6,x4,x17,x1,x8,\p
	b.cs .Lr2_4_\k
.Lb2_4_\k:
	subs x6,x12,x10
	sbcs x8,x29,x13
	sbcs x17,x23,x24
	sbcs x30,x7,x11
	csel x11, x4, xzr, lo
	subs x6,x6,x11
	b.lo .Lr2_5_\k
.Lb2_5_\k:
	stp x12,x29,[x1,#-32]
	stp x23,x7,[x1,#-16]
	FM x7,x22,x6,x23,x29,x8,x12,x17,x11,x30,x22,x21,x24,x9,x16,x13,x6,x4,x30,x29,x12
	adds x7,x7,x15
	adcs x23,x23,x20
	adcs x29,x29,x14
	adcs x12,x12,x19
	adc x30,x30,xzr
	mul x19,x30,x4
	umulh x14,x30,x4
	adds x7,x7,x19
	adcs x23,x23,x14
	b.cs .Lr2_6_\k
.Lb2_6_\k:
	stp x29,x12,[x1,#16]
	stp x7,x23,[x1],#-64
	.endm
	.macro S2 k
.Lr2_0_\k:
	RCHG
	adcs x12,x12,xzr
	adc x13,x13,xzr
	b .Lb2_0_\k
.Lr2_1_\k:
	RCHG
	adcs x27,x27,xzr
	adc x28,x28,xzr
	b .Lb2_1_\k
.Lr2_2_\k:
	RCHG
	sbcs x8,x8,xzr
	sbcs x17,x17,xzr
	sbc x6,x6,xzr
	b .Lb2_2_\k
.Lr2_3_\k:
	RCHG
	adcs x9,x9,xzr
	adc x16,x16,xzr
	b .Lb2_3_\k
.Lr2_4_\k:
	RCHG
	adcs x23,x23,xzr
	adc x7,x7,xzr
	b .Lb2_4_\k
.Lr2_5_\k:
	RCHG
	sbcs x8,x8,xzr
	sbcs x17,x17,xzr
	sbc x30,x30,xzr
	b .Lb2_5_\k
.Lr2_6_\k:
	RCHG
	adcs x29,x29,xzr
	adc x12,x12,xzr
	b .Lb2_6_\k
	.endm
	.macro P2 i,j
	ldp x0,x5,[x2,#-16]!
	W2 \i,x5
	W2 \j,x0
	.endm
	.text
	.p2align 2
	.globl pass2
	.type pass2,%function
pass2:
	stp x19,x20,[sp,#-96]!
	stp x21,x22,[sp,#16]
	stp x23,x24,[sp,#32]
	stp x25,x26,[sp,#48]
	stp x27,x28,[sp,#64]
	stp x29,x30,[sp,#80]
	ldp x25,x26,[x4]
	ldp x27,x28,[x4,#16]
	mov x4,#977
	movk x4,#1,lsl #32
	sub x5,x0,#1
	add x1,x1,x5,lsl #6
	add x1,x1,#32
	add x2,x2,x0,lsl #3
	add x3,x3,x5,lsl #6
	and x0,x5,#127
	lsr x18,x5,#7
	cbz x0, 3f
2:
	ldr x5,[x2,#-8]!
	W2 128,x5
	subs x0,x0,#1
	b.ne 2b
3:
	cbz x18, 4f
1:
	P2 0,1
	P2 2,3
	P2 4,5
	P2 6,7
	P2 8,9
	P2 10,11
	P2 12,13
	P2 14,15
	P2 16,17
	P2 18,19
	P2 20,21
	P2 22,23
	P2 24,25
	P2 26,27
	P2 28,29
	P2 30,31
	P2 32,33
	P2 34,35
	P2 36,37
	P2 38,39
	P2 40,41
	P2 42,43
	P2 44,45
	P2 46,47
	P2 48,49
	P2 50,51
	P2 52,53
	P2 54,55
	P2 56,57
	P2 58,59
	P2 60,61
	P2 62,63
	P2 64,65
	P2 66,67
	P2 68,69
	P2 70,71
	P2 72,73
	P2 74,75
	P2 76,77
	P2 78,79
	P2 80,81
	P2 82,83
	P2 84,85
	P2 86,87
	P2 88,89
	P2 90,91
	P2 92,93
	P2 94,95
	P2 96,97
	P2 98,99
	P2 100,101
	P2 102,103
	P2 104,105
	P2 106,107
	P2 108,109
	P2 110,111
	P2 112,113
	P2 114,115
	P2 116,117
	P2 118,119
	P2 120,121
	P2 122,123
	P2 124,125
	P2 126,127
	subs x18,x18,#1
	b.ne 1b
4:
	ldr x5,[x2,#-8]
	ldp x6,x7,[x5,#32]
	ldp x8,x9,[x5,#48]
	ldp x10,x11,[x1,#0]
	ldp x12,x13,[x1,#16]
	subs x6,x6,x10
	sbcs x7,x7,x11
	sbcs x8,x8,x12
	sbcs x9,x9,x13
	csel x14, x4, xzr, lo
	subs x6,x6,x14
	b.lo .Lr2_0_w0
.Lb2_0_w0:
	FM x14,x6,x25,x15,x16,x26,x17,x27,x19,x28,x6,x7,x20,x8,x9,x21,x25,x4,x28,x16,x17
	mul x20,x28,x4
	umulh x6,x28,x4
	adds x14,x14,x20
	adcs x15,x15,x6
	b.cs .Lr2_1_w0
.Lb2_1_w0:
	mul x6,x14,x17
	mul x20,x15,x16
	adds x6,x6,x20
	umulh x20,x14,x17
	umulh x28,x15,x16
	adcs x20,x20,x28
	umulh x28,x15,x17
	adc x28,x28,xzr
	umulh x19,x14,x15
	mul x25,x14,x16
	adds x19,x19,x25
	umulh x25,x14,x16
	adcs x6,x6,x25
	mul x25,x15,x17
	adcs x20,x20,x25
	mul x25,x16,x17
	adcs x28,x28,x25
	umulh x25,x16,x17
	adc x25,x25,xzr
	mul x21,x14,x15
	lsr x9,x25,#63
	extr x25,x25,x28,#63
	extr x28,x28,x20,#63
	extr x20,x20,x6,#63
	extr x6,x6,x19,#63
	extr x19,x19,x21,#63
	SQ x27,x14,x26,x21,x15,x19,x6,x16,x20,x28,x17,x25,x9,x4,x8,x1,x7,x5
	b.cs .Lr2_2_w0
.Lb2_2_w0:
	subs x9,x27,x25
	sbcs x7,x21,x28
	sbcs x8,x19,x20
	sbcs x22,x6,x26
	csel x26, x4, xzr, lo
	subs x9,x9,x26
	b.lo .Lr2_3_w0
.Lb2_3_w0:
	stp x27,x21,[x1,#-32]
	stp x19,x6,[x1,#-16]
	FM x6,x14,x9,x19,x21,x7,x27,x8,x26,x22,x14,x15,x20,x16,x17,x28,x9,x4,x22,x21,x27
	adds x6,x6,x10
	adcs x19,x19,x11
	adcs x21,x21,x12
	adcs x27,x27,x13
	adc x22,x22,xzr
	mul x13,x22,x4
	umulh x12,x22,x4
	adds x6,x6,x13
	adcs x19,x19,x12
	b.cs .Lr2_4_w0
.Lb2_4_w0:
	stp x6,x19,[x1,#0]
	stp x21,x27,[x1,#16]
	ldp x29,x30,[sp,#80]
	ldp x27,x28,[sp,#64]
	ldp x25,x26,[sp,#48]
	ldp x23,x24,[sp,#32]
	ldp x21,x22,[sp,#16]
	ldp x19,x20,[sp],#96
	ret
pass2_rare:
	S2 128
	S2 0
	S2 1
	S2 2
	S2 3
	S2 4
	S2 5
	S2 6
	S2 7
	S2 8
	S2 9
	S2 10
	S2 11
	S2 12
	S2 13
	S2 14
	S2 15
	S2 16
	S2 17
	S2 18
	S2 19
	S2 20
	S2 21
	S2 22
	S2 23
	S2 24
	S2 25
	S2 26
	S2 27
	S2 28
	S2 29
	S2 30
	S2 31
	S2 32
	S2 33
	S2 34
	S2 35
	S2 36
	S2 37
	S2 38
	S2 39
	S2 40
	S2 41
	S2 42
	S2 43
	S2 44
	S2 45
	S2 46
	S2 47
	S2 48
	S2 49
	S2 50
	S2 51
	S2 52
	S2 53
	S2 54
	S2 55
	S2 56
	S2 57
	S2 58
	S2 59
	S2 60
	S2 61
	S2 62
	S2 63
	S2 64
	S2 65
	S2 66
	S2 67
	S2 68
	S2 69
	S2 70
	S2 71
	S2 72
	S2 73
	S2 74
	S2 75
	S2 76
	S2 77
	S2 78
	S2 79
	S2 80
	S2 81
	S2 82
	S2 83
	S2 84
	S2 85
	S2 86
	S2 87
	S2 88
	S2 89
	S2 90
	S2 91
	S2 92
	S2 93
	S2 94
	S2 95
	S2 96
	S2 97
	S2 98
	S2 99
	S2 100
	S2 101
	S2 102
	S2 103
	S2 104
	S2 105
	S2 106
	S2 107
	S2 108
	S2 109
	S2 110
	S2 111
	S2 112
	S2 113
	S2 114
	S2 115
	S2 116
	S2 117
	S2 118
	S2 119
	S2 120
	S2 121
	S2 122
	S2 123
	S2 124
	S2 125
	S2 126
	S2 127
.Lr2_0_w0:
	RCHG
	sbcs x7,x7,xzr
	sbcs x8,x8,xzr
	sbc x9,x9,xzr
	b .Lb2_0_w0
.Lr2_1_w0:
	RCHG
	adcs x16,x16,xzr
	adc x17,x17,xzr
	b .Lb2_1_w0
.Lr2_2_w0:
	RCHG
	adcs x19,x19,xzr
	adc x6,x6,xzr
	b .Lb2_2_w0
.Lr2_3_w0:
	RCHG
	sbcs x7,x7,xzr
	sbcs x8,x8,xzr
	sbc x22,x22,xzr
	b .Lb2_3_w0
.Lr2_4_w0:
	RCHG
	adcs x21,x21,xzr
	adc x27,x27,xzr
	b .Lb2_4_w0

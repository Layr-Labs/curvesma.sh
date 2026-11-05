// Leaf kernels on the 52-bit-VA flipped grid (proof.md 3): MX = 256 node-permanent X's, streamed Y-records;
// C1b X-pair partial fractions, C2 split u/v, C4a LD64B beta lookups, unified 8x32-bit LD64B X-tables, exact-carry keys.

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

	.macro P1X
	RF x9,x25,x26,x27,x28,x13,x14,x15,x16,x17,x23,x24
	SHX x5,x25,x28,x26,x27
	RA x10,x5,x25,x26,x27,x6,x7,x8,x23,x24,x13,x14,x15,x16,x17,x29,xzr,x15,x16
	SHX x28,x5,x27,x25,x26
	RA x11,x28,x5,x25,x26,x6,x7,x8,x23,x24,x13,x14,x15,x16,x17,x29,x29,x15,x16
	SHX x27,x28,x26,x5,x25
	RA x12,x27,x28,x5,x25,x6,x7,x8,x9,x10,x13,x14,x15,x16,x17,x29,x29,x27,x28
	FDX x13,x14,x17,x29,x25,x26,x5,x6,x27,x28
	.endm

	// 1. C1FWD8: e = B^2 + sum_{k=0}^7 T_sig(k, t_k(b^3)) at E_{2j}+0, acc_prev at +256, P1X -> acc_new (156 instr / 2 pairs)
	.macro LKFE b,ad
	add \ad,\b,\ad
	ld64b x14,[\ad]
	adds x9,x9,x14
	adcs x10,x10,x15
	adcs x11,x11,x16
	adcs x12,x12,x17
	adc x22,x22,xzr
	.endm

	.macro LKFO b,ad
	add \ad,\b,\ad
	ld64b x14,[\ad]
	adds x9,x9,x18
	adcs x10,x10,x19
	adcs x11,x11,x20
	adcs x12,x12,x21
	adc x22,x22,xzr
	.endm

	.macro LKP8 xo
	ldr x2,[x3,#\xo]
	ldp x5,x6,[x1,#336]
	add x5,x2,x5
	ldp x9,x10,[x5]
	ldp x11,x12,[x5,#16]
	add x6,x2,x6
	ld64b x14,[x6]
	adds x9,x9,x18
	adcs x10,x10,x19
	adcs x11,x11,x20
	adcs x12,x12,x21
	cset x22,hs
	ldp x5,x6,[x1,#352]
	LKFE x2,x5
	LKFO x2,x6
	ldp x5,x6,[x1,#368]
	LKFE x2,x5
	LKFO x2,x6
	ldp x5,x6,[x1,#384]
	LKFE x2,x5
	LKFO x2,x6
	ldp x13,x14,[x1,#80]
	ldp x15,x16,[x1,#96]
	adds x9,x9,x13
	adcs x10,x10,x14
	adcs x11,x11,x15
	adcs x12,x12,x16
	adc x22,x22,xzr
	mul x18,x22,x4
	adds x9,x9,x18
	adcs x10,x10,xzr
	b.lo 9f
	adcs x11,x11,xzr
	adcs x12,x12,xzr
	csel x18,x4,xzr,hs
	add x9,x9,x18
	stadd x4,[sp]
9:
	.endm

	// Every LDP/STP offset must fit the scaled signed imm7 (-512..504), so one C1FWD8 cannot advance x3 by the
	// full 1024-byte record: C1FWD8 leaves x3 264 bytes short, C1FWD8B starts there and leaves x3 at its own
	// record, and fwd adds 1024 after each C1FWD8B (one ADD per two copies).
	.macro C1FWD8
	LKP8 352
	stp x11,x12,[x3,#16]
	stp x9,x10,[x3],#256
	stp x27,x28,[x3,#16]
	stp x25,x26,[x3],#504
	P1X
	.endm

	.macro C1FWD8B
	LKP8 616
	stp x11,x12,[x3,#280]
	stp x9,x10,[x3,#264]!
	stp x27,x28,[x3,#272]
	stp x25,x26,[x3,#256]
	P1X
	.endm

	// 2. W2SE8 + C1BSPLIT8 (324 instr / 2 pairs)
	// bwd keeps x3 512 below the old convention (record + 512 in the first half, record + 0 in the second) so that
	// every x3 offset below fits its encoding: LDP/STP imm7 (-512..504), LDR scaled imm12 or unscaled imm9.
	.macro W2SE8
	ldp x21,x22,[x3,#-256]
	ldp x23,x24,[x3,#-512]
	RF x21,x25,x26,x27,x28,x29,x30,x0,x2,x13,x5,x6
	RF x23,x25,x26,x27,x28,x15,x16,x17,x18,x19,x5,x6
	SHX x5,x25,x28,x26,x27
	RA4 x22,x5,x25,x26,x27,x29,x30,x0,x2,x13,x14,xzr,x6,x7,x8,x28
	RA4 x24,x5,x25,x26,x27,x15,x16,x17,x18,x19,x20,xzr,x6,x7,x8,x28
	ldp x21,x22,[x3,#-240]
	ldp x23,x24,[x3,#-496]
	SHX x6,x5,x27,x25,x26
	RA4 x21,x6,x5,x25,x26,x29,x30,x0,x2,x13,x14,x14,x7,x8,x28,x27
	RA4 x23,x6,x5,x25,x26,x15,x16,x17,x18,x19,x20,x20,x7,x8,x28,x27
	SHX x7,x6,x26,x5,x25
	RA4 x22,x7,x6,x5,x25,x29,x30,x0,x2,x13,x14,x14,x8,x28,x27,x26
	RA4 x24,x7,x6,x5,x25,x15,x16,x17,x18,x19,x20,x20,x8,x28,x27,x26
	FDX x29,x30,x13,x14,x29,x30,x5,x13,x0,x2
	FDX x15,x16,x19,x20,x15,x16,x6,x19,x17,x18
	stp x15,x16,[x3,#-496]
	stp x17,x18,[x3,#-512]
	.endm

	.macro BFI8G rt,e0,e1,e2,e3,d0,d1,d2,d3,b,t,s0,s1,s2,s3,o0,o1,o2,o3,c
	bfi \b,\s0,#8,#32
	ldp \o0,\o1,[\b]
	ldp \o2,\o3,[\b,#16]
	lsr \t,\s0,#32
	bfi \b,\t,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\d0
	adcs \o1,\o1,\d1
	adcs \o2,\o2,\d2
	adcs \o3,\o3,\d3
	cset \c,hs
	orr \b,\b,#64
	bfi \b,\s1,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\e0
	adcs \o1,\o1,\e1
	adcs \o2,\o2,\e2
	adcs \o3,\o3,\e3
	adc \c,\c,xzr
	lsr \t,\s1,#32
	bfi \b,\t,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\d0
	adcs \o1,\o1,\d1
	adcs \o2,\o2,\d2
	adcs \o3,\o3,\d3
	adc \c,\c,xzr
	eor \b,\b,#192
	bfi \b,\s2,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\e0
	adcs \o1,\o1,\e1
	adcs \o2,\o2,\e2
	adcs \o3,\o3,\e3
	adc \c,\c,xzr
	lsr \t,\s2,#32
	bfi \b,\t,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\d0
	adcs \o1,\o1,\d1
	adcs \o2,\o2,\d2
	adcs \o3,\o3,\d3
	adc \c,\c,xzr
	orr \b,\b,#64
	bfi \b,\s3,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\e0
	adcs \o1,\o1,\e1
	adcs \o2,\o2,\e2
	adcs \o3,\o3,\e3
	adc \c,\c,xzr
	lsr \t,\s3,#32
	bfi \b,\t,#8,#32
	ld64b \rt,[\b]
	adds \o0,\o0,\d0
	adcs \o1,\o1,\d1
	adcs \o2,\o2,\d2
	adcs \o3,\o3,\d3
	adc \c,\c,xzr
	.endm

	.macro C1BSPLIT8
	ldp x13,x14,[sp,#48]
	ldp x15,x16,[sp,#64]
	ldp x21,x22,[x3,#416]
	ldp x23,x24,[x3,#432]
	FS x21,x22,x23,x24,x13,x14,x15,x16,x21,x22,x23,x24,x5
	RF x29,x21,x22,x23,x24,x5,x6,x25,x26,x17,x7,x8
	SHX x7,x21,x24,x22,x23
	RA4 x30,x7,x21,x22,x23,x5,x6,x25,x26,x17,x18,xzr,x8,x9,x10,x24
	SHX x8,x7,x23,x21,x22
	RA4 x0,x8,x7,x21,x22,x5,x6,x25,x26,x17,x18,x18,x9,x10,x24,x23
	SHX x9,x8,x22,x7,x21
	RA4 x2,x9,x8,x7,x21,x5,x6,x25,x26,x17,x18,x18,x10,x24,x23,x22
	FDX x5,x6,x17,x18,x5,x6,x19,x17,x25,x26
	stp x5,x6,[x3]
	stp x25,x26,[x3,#16]
	ldr x22,[x3,#-152]
	BFI8G x14,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x29,x30,x0,x2,x9,x10,x11,x12,x24
	adds x9,x9,x5
	adcs x10,x10,x6
	adcs x11,x11,x25
	adcs x12,x12,x26
	adc x24,x24,xzr
	mul x18,x24,x4
	adds x9,x9,x18
	adcs x10,x10,xzr
	b.lo 9f
	adcs x11,x11,xzr
	adcs x12,x12,xzr
	csel x18,x4,xzr,hs
	add x9,x9,x18
	stadd x4,[sp]
9:
	.endm

	// 3. C2UV8 (158 instr)
	.macro C2UV8
	ldp x5,x6,[sp,#112]
	ldp x7,x8,[sp,#128]
	RF x9,x5,x6,x7,x8,x13,x14,x15,x16,x17,x19,x20
	ldp x5,x6,[sp,#144]
	ldp x7,x8,[sp,#160]
	RA4 x10,x5,x6,x7,x8,x13,x14,x15,x16,x17,x18,xzr,x19,x20,x21,x22
	ldp x5,x6,[sp,#176]
	ldp x7,x8,[sp,#192]
	RA4 x11,x5,x6,x7,x8,x13,x14,x15,x16,x17,x18,x18,x19,x20,x21,x22
	ldp x5,x6,[sp,#208]
	ldp x7,x8,[sp,#224]
	RA4 x12,x5,x6,x7,x8,x13,x14,x15,x16,x17,x18,x18,x19,x20,x21,x22
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	ldr x26,[x3,#472]
	BFI8G x18,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x9,x10,x11,x12,x5,x6,x7,x8,x28
	mul x26,x28,x4
	adds x5,x5,x26
	adcs x6,x6,xzr
	b.lo 9f
	adcs x7,x7,xzr
	adcs x8,x8,xzr
	csel x26,x4,xzr,hs
	add x5,x5,x26
	stadd x4,[sp]
9:
	FS x17,x18,x19,x20,x13,x14,x15,x16,x5,x6,x7,x8,x26
	FADDW x27,x28,x29,x30,x13,x14,x15,x16,x5,x6,x7,x8,x26
	.endm

	// 4. PAIRB6 (C4a LD64B beta lookups in TTBR1: 95 instr)
	.macro BLK6 s0,s1,s2,s3,o0,o1,o2,o3
	orr x10,x5,\s0,lsl #8
	ldp \o0,\o1,[x10]
	ldp \o2,\o3,[x10,#16]
	extr x10,\s1,\s0,#43
	orr x10,x5,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	cset x11,hs
	extr x10,\s2,\s1,#22
	orr x10,x6,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	adc x11,x11,xzr
	orr x10,x6,\s2,lsl #7
	ld64b x12,[x10]
	adds \o0,\o0,x12
	adcs \o1,\o1,x13
	adcs \o2,\o2,x14
	adcs \o3,\o3,x15
	adc x11,x11,xzr
	extr x10,\s3,\s2,#44
	orr x10,x9,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x12
	adcs \o1,\o1,x13
	adcs \o2,\o2,x14
	adcs \o3,\o3,x15
	adc x11,x11,xzr
	lsr x10,\s3,#23
	orr x10,x9,x10,lsl #8
	ld64b x12,[x10]
	adds \o0,\o0,x16
	adcs \o1,\o1,x17
	adcs \o2,\o2,x18
	adcs \o3,\o3,x19
	adc x11,x11,xzr
	mul x10,x11,x4
	adds \o0,\o0,x10
	adcs \o1,\o1,xzr
	b.lo 9f
	adcs \o2,\o2,xzr
	adcs \o3,\o3,xzr
	csel x10,x4,xzr,hs
	add \o0,\o0,x10
	stadd x4,[sp]
9:
	.endm

	.macro PAIRB6
	mov x5,#0xfff8000000000000
	orr x6,x5,#64
	orr x9,x5,#128
	BLK6 x21,x22,x23,x24,x29,x27,x7,x1
	BLK6 x25,x26,x0,x2,x30,x28,x8,x20
	.endm

	// 5. LKE8, LKFF8, LKG8 with (x0,x2,x1) spilled once at [sp,#464..487] (4 LDPs of offsets, 1144.0 instr/pair!)
	.macro LK8E base,ad
	add \ad,\base,\ad
	ld64b x12,[\ad]
	adds x5,x5,x12
	adcs x6,x6,x13
	adcs x9,x9,x14
	adcs x10,x10,x15
	adc x11,x11,xzr
	.endm

	.macro LK8O base,ad
	add \ad,\base,\ad
	ld64b x12,[\ad]
	adds x5,x5,x16
	adcs x6,x6,x17
	adcs x9,x9,x18
	adcs x10,x10,x19
	adc x11,x11,xzr
	.endm

	.macro LK8 tb_off,qa,qb,qc,qd
	ldr x1,[x3,#\tb_off]
	ldp x0,x2,[sp,#\qa]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	adds x5,x5,x16
	adcs x6,x6,x17
	adcs x9,x9,x18
	adcs x10,x10,x19
	cset x11,hs
	ldp x0,x2,[sp,#\qb]
	LK8E x1,x0
	LK8O x1,x2
	ldp x0,x2,[sp,#\qc]
	LK8E x1,x0
	LK8O x1,x2
	ldp x0,x2,[sp,#\qd]
	LK8E x1,x0
	LK8O x1,x2
	.endm

	.macro LKF8
	mul x18,x11,x4
	adds x5,x5,x18
	adcs x6,x6,xzr
	b.lo 9f
	adcs x9,x9,xzr
	adcs x10,x10,xzr
	csel x18,x4,xzr,hs
	add x5,x5,x18
	stadd x4,[sp]
9:
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm

	.macro LKE8
	LK8 448,272,288,304,320
	ldp x12,x13,[sp,#16]
	ldp x14,x15,[sp,#32]
	adds x5,x5,x12
	adcs x6,x6,x13
	adcs x9,x9,x14
	adcs x10,x10,x15
	adc x11,x11,xzr
	LKF8
	.endm

	.macro LKFF8
	LK8 456,336,352,368,384
	LKF8
	.endm

	.macro LKG8
	LK8 464,400,416,432,448
	LKF8
	.endm

	.macro TR0 m,f,s,r3,r2,r1,t
	umulh \f,\m,\r2
	madd \f,\m,\r3,\f
	mul \t,\m,\r2
	lsr \s,\t,#3
	umulh \t,\m,\r1
	add \s,\s,\t,lsr #3
	.endm

	.macro TRK m,f,s,r3,r2,r1,t
	umulh \t,\m,\r2
	add \f,\f,\t
	madd \f,\m,\r3,\f
	mul \t,\m,\r2
	add \s,\s,\t,lsr #3
	umulh \t,\m,\r1
	add \s,\s,\t,lsr #3
	.endm

	.macro TROWG op,a,b,c,d,r3,r2,r1,t
	\op \a,x13,x12,\r3,\r2,\r1,\t
	\op \b,x15,x14,\r3,\r2,\r1,\t
	\op \c,x17,x16,\r3,\r2,\r1,\t
	\op \d,x19,x18,\r3,\r2,\r1,\t
	.endm

	.macro SHXT a,d,p,q
	umulh \d,\d,x4
	adds \a,\a,\d
	b.lo 9f
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	stadd x4,[sp]
9:
	.endm

	// c in x5,x6,x9,x10 (x11 the row temp); accumulators (s,f) x12..x19 in line order: TFP4G (u2, v2, bu2, bv2),
	// TFP4B (v2, u2, bv2, bu2), so one ST64B of x12..x19 stores a component's four fractions (proof.md 3.1)
	.macro TFP4G
	TROWG TR0,x21,x25,x29,x30,x10,x9,x6,x11
	SHX x11,x5,x10,x6,x9
	TROWG TRK,x22,x26,x27,x28,x9,x6,x5,x10
	SHX x10,x11,x9,x5,x6
	TROWG TRK,x23,x0,x7,x8,x6,x5,x11,x9
	SHXT x10,x6,x11,x5
	TROWG TRK,x24,x2,x1,x20,x5,x11,x10,x6
	.endm

	.macro TFP4B
	TROWG TR0,x25,x21,x30,x29,x10,x9,x6,x11
	SHX x11,x5,x10,x6,x9
	TROWG TRK,x26,x22,x28,x27,x9,x6,x5,x10
	SHX x10,x11,x9,x5,x6
	TROWG TRK,x0,x23,x8,x7,x6,x5,x11,x9
	SHXT x10,x6,x11,x5
	TROWG TRK,x2,x24,x20,x1,x5,x11,x10,x6
	.endm

	.macro TFINC f,s,cl,ch
	add \f,\f,\s,lsr #61
	adds \s,\cl,\s,lsl #3
	adc \f,\f,\ch
	.endm

	.macro TFINS f,s
	add \f,\f,\s,lsr #61
	.endm

	// corrections (V0 at [r,#c], V1 at [r,#d]); the component's line at x3 + o by one ST64B
	.macro TOUTL r,c,d,o
	ldp x5,x6,[\r,#\c]
	ldp x9,x10,[\r,#\d]
	TFINC x13,x12,x5,x6
	TFINC x15,x14,x5,x6
	TFINC x17,x16,x9,x10
	TFINC x19,x18,x9,x10
	add x11,x3,#\o
	st64b x12,[x11]
	.endm

	.macro BWDONE
	C2UV8
	SQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
	F5 x21,x22,x23,x24,x8,x21,x22,x9,x8
	SQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
	F5 x25,x26,x0,x2,x8,x25,x26,x9,x8
	PAIRB6
	stp x0,x2,[sp,#464]
	str x1,[sp,#480]
	LKE8
	TFP4G
	TOUTL sp,240,256,64
	LKFF8
	TFP4B
	TOUTL x3,384,400,128
	LKG8
	TFP4G
	TFINS x13,x12
	TFINS x15,x14
	TFINS x17,x16
	TFINS x19,x18
	add x11,x3,#256
	st64b x12,[x11]
	.endm

	.macro BWDP
	W2SE8
	C1BSPLIT8
	BWDONE
	ldp x11,x12,[x3,#16]
	ldp x9,x10,[x3],#-512
	BWDONE
	ldp x25,x26,[x3,#16]
	ldp x27,x28,[x3],#-512
	.endm

	// 6. Hop-1 FLUSH (K1 = 16, store-release + DMB ISH) and exact-carry key pass + C2/C3 tree (JPATH 1, SJOIN 1)
	.macro FLUSH s,r
	sub \s,x16,#480
	str \s,[x15]
	ldr \r,[\s,#-16]
	add \r,\r,#1
	stp x0,\r,[\s,#-24]
	stlur x15,[\s,#-32]
	dmb ish
	.endm

	.macro LJ0 k,w0l,w0h,w1l,w1h
	adds x12,\w0l,\w1l
	adc x15,\w0h,\w1h
	cmp \w0h,\w1h
	b.hs kh\k
	cmn \w0h,x15
	b.hs kt\k
	and x15,x1,\w0l,lsr #19
	ldadd x4,x16,[x15]
	stp \w0h,\w1h,[x16]
	stp x17,x18,[x16,#16]
	tbz x16,#9,kn\k
	FLUSH \w0h,\w1h
kn\k:
	.endm

	.macro JHB k,w0l,w0h,w1l,w1h
kh\k:
	cmn \w1h,x15
	b.hs kt\k
	mvn \w0h,x15
	and x15,x1,\w1l,lsr #19
	ldadd x4,x16,[x15]
	stp \w1h,\w0h,[x16]
	stp x17,x18,[x16,#16]
	tbz x16,#9,km\k
	FLUSH \w0h,\w1h
km\k:
	.endm

	.macro JTB k,w0l,w0h,w1l,w1h
kt\k:
	mvn \w1h,x15
	bic x15,x1,x12,lsr #19
	ldadd x4,x16,[x15]
	stp \w1h,\w0h,[x16]
	stp x17,x18,[x16,#16]
	tbz x16,#9,ku\k
	FLUSH \w0h,\w1h
ku\k:
	.endm

	.macro TLDP
	ldp x27,x28,[x3],#-128
	ldp x29,x2,[x3,#160]
	ldp x19,x20,[x3,#-64]
	ldp x21,x22,[x3,#-32]
	ldp x23,x24,[x3,#16]
	ldp x25,x26,[x3,#48]
	.endm

	.macro TLDM
	ldp x27,x28,[x3,#144]
	ldp x29,x2,[x3,#176]
	ldp x19,x20,[x3,#-48]
	ldp x21,x22,[x3,#-16]
	ldp x25,x26,[x3,#32]
	ldp x23,x24,[x3],#-384
	.endm

	.macro S1
	add x18,x14,x3,lsl #30
	adds x5,x19,x23
	adc x6,x20,x24
	adds x5,x5,x27,lsl #3
	adc x6,x6,x28
	adds x7,x21,x25
	adc x8,x22,x26
	adds x7,x7,x29,lsl #3
	adc x8,x8,x2
	.endm

	.macro S2
	subs x9,x23,x27,lsl #3
	sbc x10,x24,x28
	subs x11,x25,x29,lsl #3
	sbc x13,x26,x2
	subs x19,x19,x23
	sbc x20,x20,x24
	subs x21,x21,x25
	sbc x22,x22,x26
	adds x23,x9,x11
	adc x24,x10,x13
	adds x25,x19,x23
	adc x26,x20,x24
	subs x27,x21,x9
	sbc x28,x22,x10
	.endm

	.macro S3
	subs x19,x19,x11
	sbc x20,x20,x13
	adds x21,x21,x23
	adc x22,x22,x24
	.endm

	.macro TSUM k
	S1
	LJ0 \k,x5,x6,x7,x8
	S2
	LJ0 d\k,x25,x26,x27,x28
	S3
	LJ0 e\k,x19,x20,x21,x22
	.endm

	.macro C3 k,j
	S3
	LJ0 \k,x19,x20,x21,x22
	b \j
	JHB \k,x19,x20,x21,x22
	b \j
	JTB \k,x19,x20,x21,x22
	b \j
	.endm

	.macro C2 k,j
	S2
	LJ0 \k,x25,x26,x27,x28
	C3 a\k,\j
	JHB \k,x25,x26,x27,x28
	C3 b\k,\j
	JTB \k,x25,x26,x27,x28
	C3 c\k,\j
	.endm

	.macro TJX k
	JHB e\k,x19,x20,x21,x22
	b kne\k
	JTB e\k,x19,x20,x21,x22
	b kne\k
	JHB d\k,x25,x26,x27,x28
	C3 g\k,kne\k
	JTB d\k,x25,x26,x27,x28
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
	str xzr,[sp,#-16]!
	subs x0,x0,#32
	b.lo 3f
1:
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	C1FWD8
	C1FWD8B
	add x3,x3,#1024
	tst x3,#0x1c000
	b.ne 2f
	add x1,x1,#512
2:
	subs x0,x0,#32
	b.hs 1b
3:
	add x0,x0,#32
	cbz x0,5f
4:
	C1FWD8
	add x3,x3,#264
	subs x0,x0,#2
	b.ne 4b
5:
	ldr x0,[sp],#16
	ret
	.globl bwd
bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	mov x0,sp
	stp x0,x30,[x1]
	str xzr,[sp]
	sub x3,x3,#512
1:
	BWDP
	BWDP
	add x0,x3,#512
	tst x0,#0x1fe00
	b.ne 1b
	tst x0,#0x3fffe0000
	b.eq 2f
	ldr x0,[sp]
	sub sp,sp,#512
	str x0,[sp]
	b 1b
2:
	ldr x0,[sp]
	ldp x1,x30,[sp,#-16]
	mov sp,x1
	ret
	.globl keys
keys:
	mov x1,#0x1ffffffffff8
	mov x4,#32
	mov x0,#512
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

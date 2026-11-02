// Leaf kernels (proof.md section 3): six leaves per inverted element, never formed as field elements. With a = x_X and
// b = x_Y, bwd forms E = b^4 + 2a^3 b, F = a^4 + 2a b^3 and G = 3a^2 b^2 on tabled rows, and per sign the components
// alpha' = E w, beta'' = F w and gamma = G w (w = u^2 or v^2); keys forms each leaf's two fractions from them (3.1-3.3).
// x4 holds C0 = 2^32 + 977 in fwd and bwd. SHX, FDX, F5, FS and FADDW branch over their fixup unless the rare carry or
// borrow occurred; the fixup completes the chain and adds C0 to the word at [sp], the fixup meter (3.4).
// bwd: header H = caller's sp (H = sp + 16 after bwd's push): rows of 2a^3 +0, of 2a +128, of a^2 +256 (4 rows of
// 32 B each), y_X +384, a^4 +416. Y entry (192 B): y^2, y, b, b^3, 3b^2, b^4. Batch entry (256 B; the 2^25 entries
// fill the lower half of a 2^34-aligned block): fwd writes D +0, acc +32; bwd writes alpha'+ +160, beta''+ +128,
// gamma+ +96, alpha'- +64, beta''- +32, gamma- +0 and stages R' at +224; -T(b), -T(beta b) mod 2^128 at +192 are
// written once per Y and never overwritten. bwd keeps R in x25..x28 between passes. fwd walks the entries upwards, two
// per pass, bwd and keys downwards until the entry pointer leaves the block (bit 33 set).
// keys: sp = K: G1 G2 +0, G3 G4 +16, G5 G0 +32, G3>>3 G0 +48, -T(a) +64, -T(beta a) +80, x30 saved at +96. Per sign
// it emits three leaves (t = 0, 1, 2) into their destinations' column-port record blocks (route.s): AND or BIC with
// x1 = (2^36 - 1) << 3 forms 8*d1, d1's LC entry (table at 0) and destination word; LDADD (x4 = 32) returns the
// slot, two STPs store the record, and after slot k FLUSH ends the block; JH and JT add 5 to x18.
// Calls may overwrite x0-x29 and the flags; bwd and keys save and restore x30; sp is preserved.

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

	.macro AD4 a0,a1,a2,a3,a4,a5,k0,k1,k2,k3
	adds \a0,\a0,\k0
	adcs \a1,\a1,\k1
	adcs \a2,\a2,\k2
	adcs \a3,\a3,\k3
	adcs \a4,\a4,xzr
	adc \a5,\a5,xzr
	.endm

	// two products on the shared rows B = (b0..b3) (row shifts by SHX): A += mA*B, Bacc += mB*B; l..t4 free
	.macro PAIR a0,a1,a2,a3,b0,b1,b2,b3,b0_,b1_,b2_,b3_,r0,r1,r2,r3,f0,f1,f2,f3,x0,x1,x2,x3,x4_,x5,y0,y1,y2,y3,y4,y5
	RF \a0,\r0,\r1,\r2,\r3,\x0,\x1,\x2,\x3,\x4_,\f0,\f1
	RF \b0,\r0,\r1,\r2,\r3,\y0,\y1,\y2,\y3,\y4,\f0,\f1
	SHX \f0,\r0,\r3,\r1,\r2
	RA4 \a1,\f0,\r0,\r1,\r2,\x0,\x1,\x2,\x3,\x4_,\x5,xzr,\f1,\f2,\f3,\r3
	RA4 \b1,\f0,\r0,\r1,\r2,\y0,\y1,\y2,\y3,\y4,\y5,xzr,\f1,\f2,\f3,\r3
	SHX \f1,\f0,\r2,\r0,\r1
	RA4 \a2,\f1,\f0,\r0,\r1,\x0,\x1,\x2,\x3,\x4_,\x5,\x5,\f2,\f3,\r3,\r2
	RA4 \b2,\f1,\f0,\r0,\r1,\y0,\y1,\y2,\y3,\y4,\y5,\y5,\f2,\f3,\r3,\r2
	SHX \f2,\f1,\r1,\f0,\r0
	RA4 \a3,\f2,\f1,\f0,\r0,\x0,\x1,\x2,\x3,\x4_,\x5,\x5,\f3,\r3,\r2,\r1
	RA4 \b3,\f2,\f1,\f0,\r0,\y0,\y1,\y2,\y3,\y4,\y5,\y5,\f3,\r3,\r2,\r1
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
	RF x9,x25,x26,x27,x28,x13,x14,x15,x16,x17,x23,x24
	SHX x5,x25,x28,x26,x27
	RA x10,x5,x25,x26,x27,x6,x7,x8,x23,x24,x13,x14,x15,x16,x17,x29,xzr,x15,x16
	SHX x28,x5,x27,x25,x26
	RA x11,x28,x5,x25,x26,x6,x7,x8,x23,x24,x13,x14,x15,x16,x17,x29,x29,x15,x16
	SHX x27,x28,x26,x5,x25
	RA x12,x27,x28,x5,x25,x6,x7,x8,x9,x10,x13,x14,x15,x16,x17,x29,x29,x27,x28
	FDX x13,x14,x17,x29,x25,x26,x5,x6,x27,x28
	.endm


	.macro FWDI
	ldp x7,x8,[x1,#16]
	ldp x5,x6,[x1],#192
	FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
	stp x11,x12,[x3,#16]
	stp x25,x26,[x3,#32]
	stp x27,x28,[x3,#48]
	stp x9,x10,[x3],#256
	P1X
	.endm

	.macro BWDI
	ldp x21,x22,[x3,#32]
	ldp x23,x24,[x3],#-256
	RF x21,x25,x26,x27,x28,x9,x10,x11,x12,x13,x5,x6
	RF x23,x25,x26,x27,x28,x15,x16,x17,x18,x19,x5,x6
	SHX x5,x25,x28,x26,x27
	RA4 x22,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
	RA4 x24,x5,x25,x26,x27,x15,x16,x17,x18,x19,x20,xzr,x6,x7,x8,x28
	ldp x21,x22,[x3,#304]
	ldp x23,x24,[x3,#272]
	SHX x6,x5,x27,x25,x26
	RA4 x21,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
	RA4 x23,x6,x5,x25,x26,x15,x16,x17,x18,x19,x20,x20,x7,x8,x28,x27
	SHX x7,x6,x26,x5,x25
	RA4 x22,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
	RA4 x24,x7,x6,x5,x25,x15,x16,x17,x18,x19,x20,x20,x8,x28,x27,x26
	FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12
	FDX x15,x16,x19,x20,x15,x16,x6,x19,x17,x18
	stp x15,x16,[x3,#480]
	stp x17,x18,[x3,#496]
	ldp x13,x14,[x1,#-192]!
	ldp x15,x16,[x1,#16]
	ldp x17,x18,[sp,#400]
	ldp x19,x20,[sp,#416]
	FS x21,x22,x23,x24,x13,x14,x15,x16,x17,x18,x19,x20,x5
	FADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
	PAIR x21,x22,x23,x24,x13,x14,x15,x16,x13,x14,x15,x16,x9,x10,x11,x12,x5,x6,x7,x8,x17,x18,x19,x20,x25,x26,x27,x28,x29,x30,x0,x2
	FDX x17,x18,x25,x26,x17,x18,x5,x25,x19,x20
	FDX x27,x28,x0,x2,x27,x28,x6,x0,x29,x30
	SQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
	F5 x21,x22,x23,x24,x8,x21,x22,x9,x8
	SQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
	F5 x25,x26,x0,x2,x8,x25,x26,x9,x8
	ldp x5,x6,[x1,#32]
	ldp x7,x8,[x1,#48]
	ldp x9,x10,[sp,#16]
	ldp x11,x12,[sp,#32]
	RF x5,x9,x10,x11,x12,x13,x14,x15,x16,x17,x19,x20
	ldp x9,x10,[sp,#48]
	ldp x11,x12,[sp,#64]
	RA4 x6,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,xzr,x19,x20,x27,x28
	ldp x9,x10,[sp,#80]
	ldp x11,x12,[sp,#96]
	RA4 x7,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x20,x27,x28
	ldp x9,x10,[sp,#112]
	ldp x11,x12,[sp,#128]
	RA4 x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x20,x27,x28
	ldp x9,x10,[x1,#128]
	ldp x11,x12,[x1,#144]
	AD4 x13,x14,x15,x16,x17,x18,x9,x10,x11,x12
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	PAIR x21,x22,x23,x24,x25,x26,x0,x2,x25,x26,x0,x2,x13,x14,x15,x16,x27,x28,x29,x30,x5,x6,x7,x8,x9,x10,x11,x12,x17,x18,x19,x20
	FDX x5,x6,x9,x10,x5,x6,x27,x9,x7,x8
	FDX x11,x12,x19,x20,x11,x12,x28,x19,x17,x18
	stp x5,x6,[x3,#416]
	stp x7,x8,[x3,#432]
	stp x11,x12,[x3,#320]
	stp x17,x18,[x3,#336]
	ldp x5,x6,[x1,#64]
	ldp x7,x8,[x1,#80]
	ldp x9,x10,[sp,#144]
	ldp x11,x12,[sp,#160]
	RF x5,x9,x10,x11,x12,x13,x14,x15,x16,x17,x19,x20
	ldp x9,x10,[sp,#176]
	ldp x11,x12,[sp,#192]
	RA4 x6,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,xzr,x19,x20,x27,x28
	ldp x9,x10,[sp,#208]
	ldp x11,x12,[sp,#224]
	RA4 x7,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x20,x27,x28
	ldp x9,x10,[sp,#240]
	ldp x11,x12,[sp,#256]
	RA4 x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x20,x27,x28
	ldp x9,x10,[sp,#432]
	ldp x11,x12,[sp,#448]
	AD4 x13,x14,x15,x16,x17,x18,x9,x10,x11,x12
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	PAIR x21,x22,x23,x24,x25,x26,x0,x2,x25,x26,x0,x2,x13,x14,x15,x16,x27,x28,x29,x30,x5,x6,x7,x8,x9,x10,x11,x12,x17,x18,x19,x20
	FDX x5,x6,x9,x10,x5,x6,x27,x9,x7,x8
	FDX x11,x12,x19,x20,x11,x12,x28,x19,x17,x18
	stp x5,x6,[x3,#384]
	stp x7,x8,[x3,#400]
	stp x11,x12,[x3,#288]
	stp x17,x18,[x3,#304]
	ldp x5,x6,[x1,#96]
	ldp x7,x8,[x1,#112]
	ldp x9,x10,[sp,#272]
	ldp x11,x12,[sp,#288]
	RF x5,x9,x10,x11,x12,x13,x14,x15,x16,x17,x19,x20
	ldp x9,x10,[sp,#304]
	ldp x11,x12,[sp,#320]
	RA4 x6,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,xzr,x19,x20,x27,x28
	ldp x9,x10,[sp,#336]
	ldp x11,x12,[sp,#352]
	RA4 x7,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x20,x27,x28
	ldp x9,x10,[sp,#368]
	ldp x11,x12,[sp,#384]
	RA4 x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x18,x19,x20,x27,x28
	FDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
	PAIR x21,x22,x23,x24,x25,x26,x0,x2,x25,x26,x0,x2,x13,x14,x15,x16,x27,x28,x29,x30,x5,x6,x7,x8,x9,x10,x11,x12,x17,x18,x19,x20
	FDX x5,x6,x9,x10,x5,x6,x27,x9,x7,x8
	FDX x11,x12,x19,x20,x11,x12,x28,x19,x17,x18
	stp x5,x6,[x3,#352]
	stp x7,x8,[x3,#368]
	stp x11,x12,[x3,#256]
	stp x17,x18,[x3,#272]
	ldp x25,x26,[x3,#480]
	ldp x27,x28,[x3,#496]
	.endm


	// KEY's top column F1 and the next column's shifted sum s of the limbs a0..a3; o = c6 = F1 + (s >> 61)
	.macro F1C a0,a1,a2,a3,o,s,t
	umulh \o,\a0,x10
	umulh \t,\a1,x9
	add \o,\o,\t
	umulh \t,\a2,x8
	add \o,\o,\t
	umulh \t,\a3,x7
	add \o,\o,\t
	madd \o,\a0,x11,\o
	madd \o,\a1,x10,\o
	madd \o,\a2,x9,\o
	madd \o,\a3,x8,\o
	umulh \s,\a0,x15
	mul \t,\a0,x10
	add \s,\s,\t,lsr #3
	mul \t,\a1,x9
	add \s,\s,\t,lsr #3
	umulh \t,\a1,x8
	add \s,\s,\t,lsr #3
	mul \t,\a2,x8
	add \s,\s,\t,lsr #3
	umulh \t,\a2,x7
	add \s,\s,\t,lsr #3
	mul \t,\a3,x7
	add \s,\s,\t,lsr #3
	umulh \t,\a3,x16
	add \s,\s,\t,lsr #3
	add \o,\o,\s,lsr #61
	.endm

// After slot k (x16; x15 = 8*d1): the LC entry returns to slot 1 and the RC cursor (2^43 + 8*d1) advances by 32k (x0);
// the header (8*d1, cursor) sends the block if the old cursor is below the slot end x12, else the block is dropped.
// s and r are two of the leaf's registers, dead once its record is stored.
	.macro FLUSH s,r
	sub \s,x16,#4080
	sub \s,\s,#4080
	str \s,[x15]
	orr \r,x15,#0x80000000000
	ldadd x0,x16,[\r]
	cmp x16,x12
	b.hs 6f
	stp x15,x16,[\s,#-32]
6:
	.endm

	// KEY's index on (W0, W1) and the j = 0 emit; j = 1, 2 branch to kh\k, kt\k (LJX)
	.macro LJ0 k,w0l,w0h,w1l,w1h
	add x15,\w0h,\w1h
	cmp \w0h,\w1h
	b.hs kh\k
	cmn \w0h,x15
	b.hs kt\k
	and x15,x1,\w0l,lsr #25
	ldadd x4,x16,[x15]
	lsr \w1l,\w1l,#28
	add \w1l,\w1l,x3,lsl #31
	stp \w0h,\w1h,[x16]
	stp x17,\w1l,[x16,#16]
	tbz x16,#13,kn\k
	FLUSH \w0h,\w1h
kn\k:
	.endm

	.macro LJX k,w0l,w0h,w1l,w1h
kh\k:
	cmn \w1h,x15
	b.hs kt\k
	adds \w0l,\w0l,\w1l
	adc \w0h,\w0h,\w1h
	mvn \w0h,\w0h
	and x15,x1,\w1l,lsr #25
	ldadd x4,x16,[x15]
	eor \w0l,x14,\w0l,lsr #28
	add \w0l,\w0l,x3,lsl #31
	stp \w1h,\w0h,[x16]
	stp x17,\w0l,[x16,#16]
	tbz x16,#13,km\k
	FLUSH \w0h,\w1h
km\k:
	add x18,x18,#5
	b kn\k
kt\k:
	adds \w1l,\w0l,\w1l
	adc \w1h,\w0h,\w1h
	mvn \w1h,\w1h
	bic x15,x1,\w1l,lsr #25
	ldadd x4,x16,[x15]
	lsr \w0l,\w0l,#28
	add \w0l,\w0l,x3,lsl #31
	stp \w1h,\w0h,[x16]
	stp x17,\w0l,[x16,#16]
	tbz x16,#13,ku\k
	FLUSH \w0h,\w1h
ku\k:
	add x18,x18,#5
	b kn\k
	.endm

	.macro KSIGN nb,nc,gp,k0,k1,k2
	ldp x21,x22,[x3,#16]
	ldp x19,x20,[x3],#-32
	F1C x19,x20,x21,x22,x23,x24,x25
	ldp x19,x20,[x3,#\nb]
	ldp x25,x26,[x3,#\nc]
	adds x21,x21,x19
	adc x22,x22,x20
	adds x24,x25,x24,lsl #3
	adc x23,x23,x26
	ldp x25,x26,[x3,#16]
	ldp x19,x20,[x3],#-32
	F1C x19,x20,x25,x26,x27,x28,x29
	ldp x19,x20,[sp,#64]
	ldp x29,x30,[sp,#80]
	adds x25,x25,x19
	adc x26,x26,x20
	adds x28,x29,x28,lsl #3
	adc x27,x27,x30
	ldp x29,x30,[x3,#16]
	ldp x19,x20,[x3],#\gp
	F1C x19,x20,x29,x30,x2,x5,x6
	adds x6,x21,x25
	adc x13,x22,x26
	adds x6,x6,x29
	adc x13,x13,x30
	adds x19,x28,x24
	adc x20,x27,x23
	adds x19,x19,x5,lsl #3
	adc x20,x20,x2
	LJ0 \k0,x6,x13,x19,x20
	adds x6,x25,x28
	adc x13,x26,x27
	adds x19,x29,x5,lsl #3
	adc x20,x30,x2
	adds x28,x21,x28
	adc x27,x22,x27
	subs x28,x28,x19
	sbc x27,x27,x20
	adds x29,x29,x24
	adc x30,x30,x23
	subs x29,x29,x6
	sbc x30,x30,x13
	LJ0 \k1,x28,x27,x29,x30
	adds x21,x21,x5,lsl #3
	adc x22,x22,x2
	subs x21,x21,x6
	sbc x22,x22,x13
	adds x25,x25,x24
	adc x26,x26,x23
	subs x25,x25,x19
	sbc x26,x26,x20
	LJ0 \k2,x21,x22,x25,x26
	ldp x15,x16,[sp,#48]
	.endm


	.text
	.p2align 2
	.globl fwd
fwd:
	mov x4,#977
	movk x4,#1,lsl #32
	str xzr,[sp,#-16]!
1:
	FWDI
	FWDI
	tbz x3,#33,1b
	ldr x0,[sp],#16
	ret
	.globl bwd
bwd:
	mov x4,#977
	movk x4,#1,lsl #32
	stp xzr,x30,[sp,#-16]!
1:
	BWDI
	BWDI
	tbz x3,#33,1b
	ldp x0,x30,[sp],#16
	ret
	.globl keys
keys:
	str x30,[sp,#96]
	mov x1,#0x7ffffffff8
	mov x4,#32
	mov x0,#8192
	ldp x7,x8,[sp,#0]
	ldp x9,x10,[sp,#16]
	ldp x11,x16,[sp,#32]
	ldr x15,[sp,#48]
	mov x14,#0xfffffffff
1:
	KSIGN 64,80,-32,p00,p01,p02
	KSIGN 160,176,-96,m00,m01,m02
	KSIGN 64,80,-32,p10,p11,p12
	KSIGN 160,176,-96,m10,m11,m12
	tbz x3,#33,1b
	ldr x30,[sp,#96]
	ret
	LJX p00,x6,x13,x19,x20
	LJX p01,x28,x27,x29,x30
	LJX p02,x21,x22,x25,x26
	LJX m00,x6,x13,x19,x20
	LJX m01,x28,x27,x29,x30
	LJX m02,x21,x22,x25,x26
	LJX p10,x6,x13,x19,x20
	LJX p11,x28,x27,x29,x30
	LJX p12,x21,x22,x25,x26
	LJX m10,x6,x13,x19,x20
	LJX m11,x28,x27,x29,x30
	LJX m12,x21,x22,x25,x26

// side_lanes.S -- counted-X kernel (925 mean charge, 932 worst path), derived by
// measure/derive_k925.py from the validated k27r1 929 listing (K3: X-expensive counting).
// Interface, arithmetic, K1/K2/K3: proof.md s1. args[0..n-1] = {out, P, T, inv},
// args[n][0] = NULL ends the lane loop; out[0..11] = P-T reduced, out[12..14] = key.
#ifndef __SIZEOF_POINTER__
#define __SIZEOF_POINTER__ 8
#endif
.text
.globl	side_lanes
.type	side_lanes,%function
.align	5
side_lanes:
	ldr	x2,[x0]
	cbz	x2,.Lnone
	stp	x29,x30,[sp,#-288]!
	add	x29,sp,#0
	stp	x19,x20,[sp,#16]
	stp	x21,x22,[sp,#32]
	stp	x23,x24,[sp,#48]
	stp	x25,x26,[sp,#64]
	stp	x27,x28,[sp,#80]
	str	x18,[sp,#272]
	sub	sp,sp,#144		// [sp,#0..95] 768-bit square, [sp,#96..143] mu
	mov	x30,x0			// lane cursor
	ldr	x1,[x30,#8]		// P of lane 0 (later lanes: loaded with the next out)
	adrp	x0,.Lsk_consts
	add	x0,x0,:lo12:.Lsk_consts
	ldr	x18,[x0,#48]		// n0, pinned in x18
	ldp	x5,x6,[x0]		// p -> x5..x10 and frame +128
	ldp	x7,x8,[x0,#16]
	ldp	x9,x10,[x0,#32]
	stp	x5,x6,[x29,#128]
	stp	x7,x8,[x29,#144]
	stp	x9,x10,[x29,#160]
	str	xzr,[sp,#240]	// X-count init (s24)
	// loop invariant at .Llane: x18 = n0, x5..x10 = p, x30 = &args[i], x1 = P
.Llane:
	ldp	x3,x2,[x30,#16]		// T, inv

	// t = y1 + y2 < 2p
	ldp	x11,x12,[x1,#48]
	ldp	x13,x14,[x1,#64]
	ldp	x15,x16,[x1,#80]
	ldp	x19,x20,[x3,#48]
	ldp	x21,x22,[x3,#64]
	ldp	x23,x24,[x3,#80]
	adds	x11,x11,x19
	adcs	x12,x12,x20
	adcs	x13,x13,x21
	adcs	x14,x14,x22
	adcs	x15,x15,x23
	adc	x16,x16,x24

	// mu = t * inv, unreduced, in x11..x16
	ldp	x17,x29,[x2]
	// ---- inlined MULNR
	mul	x19,x11,x17
	mul	x20,x12,x17
	mul	x21,x13,x17
	mul	x22,x14,x17
	mul	x23,x15,x17
	mul	x24,x16,x17
	mul	x4,x18,x19

	umulh	x26,x11,x17
	umulh	x27,x12,x17
	umulh	x28,x13,x17
	umulh	x0,x14,x17
	umulh	x1,x15,x17
	umulh	x3,x16,x17

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adc	x25,xzr,    x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x29
	adcs	x20,x21,x27
	mul	x27,x12,x29
	adcs	x21,x22,x28
	mul	x28,x13,x29
	adcs	x22,x23,x0
	mul	x0,x14,x29
	adcs	x23,x24,x1
	mul	x1,x15,x29
	adcs	x24,x25,x3
	mul	x3,x16,x29

	adds	x19,x19,x26
	umulh	x26,x11,x29
	adcs	x20,x20,x27
	umulh	x27,x12,x29
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x29
	adcs	x22,x22,x0
	umulh	x0,x14,x29
	adcs	x23,x23,x1
	umulh	x1,x15,x29
	adcs	x24,x24,x3
	umulh	x3,x16,x29
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	ldp	x17,x29,[x2,#16]

	adds	x19,x20,x26
	mul	x26,x11,x17
	adcs	x20,x21,x27
	mul	x27,x12,x17
	adcs	x21,x22,x28
	mul	x28,x13,x17
	adcs	x22,x23,x0
	mul	x0,x14,x17
	adcs	x23,x24,x1
	mul	x1,x15,x17
	adcs	x24,x25,x3
	mul	x3,x16,x17

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x29
	adcs	x20,x21,x27
	mul	x27,x12,x29
	adcs	x21,x22,x28
	mul	x28,x13,x29
	adcs	x22,x23,x0
	mul	x0,x14,x29
	adcs	x23,x24,x1
	mul	x1,x15,x29
	adcs	x24,x25,x3
	mul	x3,x16,x29

	adds	x19,x19,x26
	umulh	x26,x11,x29
	adcs	x20,x20,x27
	umulh	x27,x12,x29
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x29
	adcs	x22,x22,x0
	umulh	x0,x14,x29
	adcs	x23,x23,x1
	umulh	x1,x15,x29
	adcs	x24,x24,x3
	umulh	x3,x16,x29
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	ldp	x17,x29,[x2,#32]

	adds	x19,x20,x26
	mul	x26,x11,x17
	adcs	x20,x21,x27
	mul	x27,x12,x17
	adcs	x21,x22,x28
	mul	x28,x13,x17
	adcs	x22,x23,x0
	mul	x0,x14,x17
	adcs	x23,x24,x1
	mul	x1,x15,x17
	adcs	x24,x25,x3
	mul	x3,x16,x17

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x29
	adcs	x20,x21,x27
	mul	x27,x12,x29
	adcs	x21,x22,x28
	mul	x28,x13,x29
	adcs	x22,x23,x0
	mul	x0,x14,x29
	adcs	x23,x24,x1
	mul	x1,x15,x29
	adcs	x24,x25,x3
	mul	x3,x16,x29

	adds	x19,x19,x26
	umulh	x26,x11,x29
	adcs	x20,x20,x27
	umulh	x27,x12,x29
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x29
	adcs	x22,x22,x0
	umulh	x0,x14,x29
	adcs	x23,x23,x1
	umulh	x1,x15,x29
	adcs	x24,x24,x3
	umulh	x3,x16,x29
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	// pull r_ptr (n0 pinned in x18; out pointer loaded by the driver)

	adds	x11,x20,x26
	adcs	x12,x21,x27
	adcs	x13,x22,x28
	adcs	x14,x23,x0
	adcs	x15,x24,x1
	adcs	x16,x25,x3

	// ---- end MULNR
	stp	x11,x12,[sp,#96]
	stp	x13,x14,[sp,#112]
	stp	x15,x16,[sp,#128]

	// V = redc(sqr_384(mu)), unreduced (< 1.15p), all in registers
	// ---- inlined SQR
	mul	x19,x12,x11
	mul	x20,x13,x11
	mul	x21,x14,x11
	mul	x22,x15,x11
	mul	x23,x16,x11

	umulh	x0,x12,x11
	umulh	x2,x13,x11
	umulh	x3,x14,x11
	umulh	x4,x15,x11
	adds	x20,x20,x0
	umulh	x29,x16,x11
	adcs	x21,x21,x2
	mul	x2,x13,x12
	adcs	x22,x22,x3
	mul	x3,x14,x12
	adcs	x23,x23,x4
	mul	x4,x15,x12
	adc	x24,xzr,    x29
	mul	x29,x16,x12

	adds	x21,x21,x2
	umulh	x2,x13,x12
	adcs	x22,x22,x3
	umulh	x3,x14,x12
	adcs	x23,x23,x4
	umulh	x4,x15,x12
	adcs	x24,x24,x29
	umulh	x29,x16,x12
	adc	x25,xzr,xzr

	mul	x17,x11,x11
	adds	x22,x22,x2
	umulh	x11,  x11,x11
	adcs	x23,x23,x3
	mul	x3,x14,x13
	adcs	x24,x24,x4
	mul	x4,x15,x13
	adc	x25,x25,x29
	mul	x29,x16,x13

	adds	x23,x23,x3
	umulh	x3,x14,x13
	adcs	x24,x24,x4
	umulh	x4,x15,x13
	adcs	x25,x25,x29
	umulh	x29,x16,x13
	adc	x26,xzr,xzr

	mul	x0,x12,x12
	adds	x24,x24,x3
	umulh	x12,  x12,x12
	adcs	x25,x25,x4
	mul	x4,x15,x14
	adc	x26,x26,x29
	mul	x29,x16,x14

	adds	x25,x25,x4
	umulh	x4,x15,x14
	adcs	x26,x26,x29
	umulh	x29,x16,x14
	adc	x27,xzr,xzr
	mul	x2,x13,x13
	adds	x26,x26,x4
	umulh	x13,  x13,x13
	adc	x27,x27,x29
	mul	x3,x14,x14

	mul	x29,x16,x15
	umulh	x14,  x14,x14
	adds	x27,x27,x29
	umulh	x29,x16,x15
	mul	x4,x15,x15
	adc	x28,x29,xzr

	adds	x19,x19,x19
	adcs	x20,x20,x20
	adcs	x21,x21,x21
	adcs	x22,x22,x22
	adcs	x23,x23,x23
	adcs	x24,x24,x24
	adcs	x25,x25,x25
	adcs	x26,x26,x26
	umulh	x15,  x15,x15
	adcs	x27,x27,x27
	mul	x29,x16,x16
	adcs	x28,x28,x28
	umulh	x16,  x16,x16

	adds	x19,x19,x11
	adcs	x20,x20,x0
	adcs	x21,x21,x12
	adcs	x22,x22,x2
	adcs	x23,x23,x13
	adcs	x24,x24,x3
	adcs	x25,x25,x14
	adcs	x26,x26,x4
	adcs	x27,x27,x15
	adcs	x28,x28,x29
	adc	x16,x16,xzr

	// ---- end SQR
	// ---- inlined BY1
	mul	x15,x18,x17

	// mul	x0,x5,x15
	mul	x1,x6,x15
	mul	x3,x7,x15
	mul	x11,x8,x15
	mul	x12,x9,x15
	mul	x13,x10,x15
	subs	xzr,x17,#1		// adds	x0,x0,x17
	umulh	x17,x5,x15
	adcs	x1,x1,x19
	umulh	x19,x6,x15
	adcs	x3,x3,x20
	umulh	x20,x7,x15
	adcs	x11,x11,x21
	umulh	x21,x8,x15
	adcs	x12,x12,x22
	umulh	x22,x9,x15
	adcs	x13,x13,x23
	umulh	x23,x10,x15
	adc	x14,xzr,xzr
	adds	x17,x17,x1
	adcs	x19,x19,x3
	adcs	x20,x20,x11
	mul	x15,x18,x17
	adcs	x21,x21,x12
	adcs	x22,x22,x13
	adc	x23,x23,x14

	// mul	x0,x5,x15
	mul	x1,x6,x15
	mul	x3,x7,x15
	mul	x11,x8,x15
	mul	x12,x9,x15
	mul	x13,x10,x15
	subs	xzr,x17,#1		// adds	x0,x0,x17
	umulh	x17,x5,x15
	adcs	x1,x1,x19
	umulh	x19,x6,x15
	adcs	x3,x3,x20
	umulh	x20,x7,x15
	adcs	x11,x11,x21
	umulh	x21,x8,x15
	adcs	x12,x12,x22
	umulh	x22,x9,x15
	adcs	x13,x13,x23
	umulh	x23,x10,x15
	adc	x14,xzr,xzr
	adds	x17,x17,x1
	adcs	x19,x19,x3
	adcs	x20,x20,x11
	mul	x15,x18,x17
	adcs	x21,x21,x12
	adcs	x22,x22,x13
	adc	x23,x23,x14

	// mul	x0,x5,x15
	mul	x1,x6,x15
	mul	x3,x7,x15
	mul	x11,x8,x15
	mul	x12,x9,x15
	mul	x13,x10,x15
	subs	xzr,x17,#1		// adds	x0,x0,x17
	umulh	x17,x5,x15
	adcs	x1,x1,x19
	umulh	x19,x6,x15
	adcs	x3,x3,x20
	umulh	x20,x7,x15
	adcs	x11,x11,x21
	umulh	x21,x8,x15
	adcs	x12,x12,x22
	umulh	x22,x9,x15
	adcs	x13,x13,x23
	umulh	x23,x10,x15
	adc	x14,xzr,xzr
	adds	x17,x17,x1
	adcs	x19,x19,x3
	adcs	x20,x20,x11
	mul	x15,x18,x17
	adcs	x21,x21,x12
	adcs	x22,x22,x13
	adc	x23,x23,x14

	// mul	x0,x5,x15
	mul	x1,x6,x15
	mul	x3,x7,x15
	mul	x11,x8,x15
	mul	x12,x9,x15
	mul	x13,x10,x15
	subs	xzr,x17,#1		// adds	x0,x0,x17
	umulh	x17,x5,x15
	adcs	x1,x1,x19
	umulh	x19,x6,x15
	adcs	x3,x3,x20
	umulh	x20,x7,x15
	adcs	x11,x11,x21
	umulh	x21,x8,x15
	adcs	x12,x12,x22
	umulh	x22,x9,x15
	adcs	x13,x13,x23
	umulh	x23,x10,x15
	adc	x14,xzr,xzr
	adds	x17,x17,x1
	adcs	x19,x19,x3
	adcs	x20,x20,x11
	mul	x15,x18,x17
	adcs	x21,x21,x12
	adcs	x22,x22,x13
	adc	x23,x23,x14

	// mul	x0,x5,x15
	mul	x1,x6,x15
	mul	x3,x7,x15
	mul	x11,x8,x15
	mul	x12,x9,x15
	mul	x13,x10,x15
	subs	xzr,x17,#1		// adds	x0,x0,x17
	umulh	x17,x5,x15
	adcs	x1,x1,x19
	umulh	x19,x6,x15
	adcs	x3,x3,x20
	umulh	x20,x7,x15
	adcs	x11,x11,x21
	umulh	x21,x8,x15
	adcs	x12,x12,x22
	umulh	x22,x9,x15
	adcs	x13,x13,x23
	umulh	x23,x10,x15
	adc	x14,xzr,xzr
	adds	x17,x17,x1
	adcs	x19,x19,x3
	adcs	x20,x20,x11
	mul	x15,x18,x17
	adcs	x21,x21,x12
	adcs	x22,x22,x13
	adc	x23,x23,x14

	// mul	x0,x5,x15
	mul	x1,x6,x15
	mul	x3,x7,x15
	mul	x11,x8,x15
	mul	x12,x9,x15
	mul	x13,x10,x15
	subs	xzr,x17,#1		// adds	x0,x0,x17
	umulh	x17,x5,x15
	adcs	x1,x1,x19
	umulh	x19,x6,x15
	adcs	x3,x3,x20
	umulh	x20,x7,x15
	adcs	x11,x11,x21
	umulh	x21,x8,x15
	adcs	x12,x12,x22
	umulh	x22,x9,x15
	adcs	x13,x13,x23
	umulh	x23,x10,x15
	adc	x14,xzr,xzr
	adds	x17,x17,x1
	adcs	x19,x19,x3
	adcs	x20,x20,x11
	adcs	x21,x21,x12
	adcs	x22,x22,x13
	adc	x23,x23,x14

	// ---- end BY1
	adds	x11,x17,x24
	adcs	x12,x19,x25
	adcs	x13,x20,x26
	adcs	x14,x21,x27
	adcs	x15,x22,x28
	adc	x16,x23,x16

	// x' = V - x1 - x2 = X mod p with X = V + v2 - x1 in (-p, 2.15p), v2 = p - x2
	ldp	x0,x29,[x30,#8]		// P, T (T stays in x29 through the second product)
	ldp	x19,x20,[x0]		// x1
	ldp	x21,x22,[x0,#16]
	ldp	x23,x24,[x0,#32]
	ldp	x25,x26,[x29,#96]	// v2 = p - x2
	ldp	x27,x28,[x29,#112]
	ldp	x17,x3,[x29,#128]
	adds	x11,x11,x25
	adcs	x12,x12,x26
	adcs	x13,x13,x27
	adcs	x14,x14,x28
	adcs	x15,x15,x17
	adc	x16,x16,x3
	subs	x11,x11,x19
	sbcs	x12,x12,x20
	sbcs	x13,x13,x21
	sbcs	x14,x14,x22
	sbcs	x15,x15,x23
	sbcs	x16,x16,x24		// lo iff X < 0
	b.hs	.Lx_nonneg
	adds	x11,x11,x5		// X + p
	adcs	x12,x12,x6
	adcs	x13,x13,x7
	adcs	x14,x14,x8
	adcs	x15,x15,x9
	adc	x16,x16,x10
	b	.Lx_done
.Lx_nonneg:
	subs	x19,x11,x5		// Y = X - p
	sbcs	x20,x12,x6
	sbcs	x21,x13,x7
	sbcs	x22,x14,x8
	sbcs	x23,x15,x9
	sbcs	x24,x16,x10
	b.lo	.Lx_done		// X < p
	ldr	x0,[sp,#240]	// s24: count X>=p lanes
	adc	x0,x0,xzr
	str	x0,[sp,#240]
	subs	x11,x19,x5		// Z = Y - p
	sbcs	x12,x20,x6
	sbcs	x13,x21,x7
	sbcs	x14,x22,x8
	sbcs	x15,x23,x9
	sbcs	x16,x24,x10
	b.lo	.Lx_y			// Y < p: x' = Y
.Lx_done:
	ldr	x2,[x30]		// out (stays in x2 through the second product)
	stp	x11,x12,[x2]		// x'
	stp	x13,x14,[x2,#16]
	stp	x15,x16,[x2,#32]
	adds	x11,x11,x25		// w = x' + v2 = x' - x2 (mod p), in (0, 2p)
	adcs	x12,x12,x26
	adcs	x13,x13,x27
	adcs	x14,x14,x28
	adcs	x15,x15,x17
	adc	x16,x16,x3
	// continuation, copy for this path

	// z = w * mu, unreduced (< 1.25p) in x19..x24
	ldp	x17,x2,[sp,#96]
	// ---- inlined MUL
	mul	x19,x11,x17
	mul	x20,x12,x17
	mul	x21,x13,x17
	mul	x22,x14,x17
	mul	x23,x15,x17
	mul	x24,x16,x17
	mul	x4,x18,x19

	umulh	x26,x11,x17
	umulh	x27,x12,x17
	umulh	x28,x13,x17
	umulh	x0,x14,x17
	umulh	x1,x15,x17
	umulh	x3,x16,x17

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adc	x25,xzr,    x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x2
	adcs	x20,x21,x27
	mul	x27,x12,x2
	adcs	x21,x22,x28
	mul	x28,x13,x2
	adcs	x22,x23,x0
	mul	x0,x14,x2
	adcs	x23,x24,x1
	mul	x1,x15,x2
	adcs	x24,x25,x3
	mul	x3,x16,x2

	adds	x19,x19,x26
	umulh	x26,x11,x2
	adcs	x20,x20,x27
	umulh	x27,x12,x2
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x2
	adcs	x22,x22,x0
	umulh	x0,x14,x2
	adcs	x23,x23,x1
	umulh	x1,x15,x2
	adcs	x24,x24,x3
	umulh	x3,x16,x2
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	ldp	x17,x2,[sp,#112]

	adds	x19,x20,x26
	mul	x26,x11,x17
	adcs	x20,x21,x27
	mul	x27,x12,x17
	adcs	x21,x22,x28
	mul	x28,x13,x17
	adcs	x22,x23,x0
	mul	x0,x14,x17
	adcs	x23,x24,x1
	mul	x1,x15,x17
	adcs	x24,x25,x3
	mul	x3,x16,x17

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x2
	adcs	x20,x21,x27
	mul	x27,x12,x2
	adcs	x21,x22,x28
	mul	x28,x13,x2
	adcs	x22,x23,x0
	mul	x0,x14,x2
	adcs	x23,x24,x1
	mul	x1,x15,x2
	adcs	x24,x25,x3
	mul	x3,x16,x2

	adds	x19,x19,x26
	umulh	x26,x11,x2
	adcs	x20,x20,x27
	umulh	x27,x12,x2
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x2
	adcs	x22,x22,x0
	umulh	x0,x14,x2
	adcs	x23,x23,x1
	umulh	x1,x15,x2
	adcs	x24,x24,x3
	umulh	x3,x16,x2
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	ldp	x17,x2,[sp,#128]

	adds	x19,x20,x26
	mul	x26,x11,x17
	adcs	x20,x21,x27
	mul	x27,x12,x17
	adcs	x21,x22,x28
	mul	x28,x13,x17
	adcs	x22,x23,x0
	mul	x0,x14,x17
	adcs	x23,x24,x1
	mul	x1,x15,x17
	adcs	x24,x25,x3
	mul	x3,x16,x17

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x2
	adcs	x20,x21,x27
	mul	x27,x12,x2
	adcs	x21,x22,x28
	mul	x28,x13,x2
	adcs	x22,x23,x0
	mul	x0,x14,x2
	adcs	x23,x24,x1
	mul	x1,x15,x2
	adcs	x24,x25,x3
	mul	x3,x16,x2

	adds	x19,x19,x26
	umulh	x26,x11,x2
	adcs	x20,x20,x27
	umulh	x27,x12,x2
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x2
	adcs	x22,x22,x0
	umulh	x0,x14,x2
	adcs	x23,x23,x1
	umulh	x1,x15,x2
	adcs	x24,x24,x3
	umulh	x3,x16,x2
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	// pull r_ptr (n0 pinned in x18; out pointer loaded by the driver)

	adds	x19,x20,x26
	adcs	x20,x21,x27
	adcs	x21,x22,x28
	adcs	x22,x23,x0
	adcs	x23,x24,x1
	adcs	x24,x25,x3

	// ---- end MUL
	ldr	x2,[x30]		// out (x2 carried b words 1, 3, 5)

	// y' = mu (x' - x2) + y2 = z - n2 (mod p), n2 = p - y2; D = z - n2 in (-p, 1.25p)
	ldp	x25,x26,[x29,#144]	// n2
	ldp	x27,x28,[x29,#160]
	ldp	x17,x3,[x29,#176]
	subs	x11,x19,x25
	sbcs	x12,x20,x26
	sbcs	x13,x21,x27
	sbcs	x14,x22,x28
	sbcs	x15,x23,x17
	sbcs	x16,x24,x3		// lo iff D < 0
	b.hs	.Ly_nonneg2
	adds	x11,x11,x5		// y' = D + p
	adcs	x12,x12,x6
	adcs	x13,x13,x7
	adcs	x14,x14,x8
	adcs	x15,x15,x9
	adc	x16,x16,x10
	b	.Ly_store2
.Ly_nonneg2:
	subs	x19,x11,x5		// E = D - p
	sbcs	x20,x12,x6
	sbcs	x21,x13,x7
	sbcs	x22,x14,x8
	sbcs	x23,x15,x9
	sbcs	x24,x16,x10
	b.lo	.Ly_store2		// D < p: y' = D
	stp	x19,x20,[x2,#48]	// y' = E
	stp	x21,x22,[x2,#64]
	stp	x23,x24,[x2,#80]
	tbz	x19,#0,.Lk_e2		// key: even element of {y', p - y'}, limbs 1..3
	subs	x19,x5,x19
	sbcs	x20,x6,x20
	sbcs	x21,x7,x21
	sbc	x22,x8,x22
.Lk_e2:
	stp	x20,x21,[x2,#96]
	str	x22,[x2,#112]
	ldp	x0,x1,[x30,#32]!	// next lane's out (NULL at args[n]) and P
	cbnz	x0,.Llane
	b	.Lexit
.Ly_store2:
	stp	x11,x12,[x2,#48]	// y'
	stp	x13,x14,[x2,#64]
	stp	x15,x16,[x2,#80]
	tbz	x11,#0,.Lk_d2
	subs	x11,x5,x11
	sbcs	x12,x6,x12
	sbcs	x13,x7,x13
	sbc	x14,x8,x14
.Lk_d2:
	stp	x12,x13,[x2,#96]
	str	x14,[x2,#112]
	ldp	x0,x1,[x30,#32]!	// next lane's out (NULL at args[n]) and P
	cbnz	x0,.Llane
	b	.Lexit
.Lx_y:
	ldr	x2,[x30]		// out
	stp	x19,x20,[x2]		// x' = Y
	stp	x21,x22,[x2,#16]
	stp	x23,x24,[x2,#32]
	adds	x11,x19,x25		// w = x' + v2
	adcs	x12,x20,x26
	adcs	x13,x21,x27
	adcs	x14,x22,x28
	adcs	x15,x23,x17
	adc	x16,x24,x3
.Lw_done:

	// z = w * mu, unreduced (< 1.25p) in x19..x24
	ldp	x17,x2,[sp,#96]
	// ---- inlined MUL
	mul	x19,x11,x17
	mul	x20,x12,x17
	mul	x21,x13,x17
	mul	x22,x14,x17
	mul	x23,x15,x17
	mul	x24,x16,x17
	mul	x4,x18,x19

	umulh	x26,x11,x17
	umulh	x27,x12,x17
	umulh	x28,x13,x17
	umulh	x0,x14,x17
	umulh	x1,x15,x17
	umulh	x3,x16,x17

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adc	x25,xzr,    x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x2
	adcs	x20,x21,x27
	mul	x27,x12,x2
	adcs	x21,x22,x28
	mul	x28,x13,x2
	adcs	x22,x23,x0
	mul	x0,x14,x2
	adcs	x23,x24,x1
	mul	x1,x15,x2
	adcs	x24,x25,x3
	mul	x3,x16,x2

	adds	x19,x19,x26
	umulh	x26,x11,x2
	adcs	x20,x20,x27
	umulh	x27,x12,x2
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x2
	adcs	x22,x22,x0
	umulh	x0,x14,x2
	adcs	x23,x23,x1
	umulh	x1,x15,x2
	adcs	x24,x24,x3
	umulh	x3,x16,x2
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	ldp	x17,x2,[sp,#112]

	adds	x19,x20,x26
	mul	x26,x11,x17
	adcs	x20,x21,x27
	mul	x27,x12,x17
	adcs	x21,x22,x28
	mul	x28,x13,x17
	adcs	x22,x23,x0
	mul	x0,x14,x17
	adcs	x23,x24,x1
	mul	x1,x15,x17
	adcs	x24,x25,x3
	mul	x3,x16,x17

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x2
	adcs	x20,x21,x27
	mul	x27,x12,x2
	adcs	x21,x22,x28
	mul	x28,x13,x2
	adcs	x22,x23,x0
	mul	x0,x14,x2
	adcs	x23,x24,x1
	mul	x1,x15,x2
	adcs	x24,x25,x3
	mul	x3,x16,x2

	adds	x19,x19,x26
	umulh	x26,x11,x2
	adcs	x20,x20,x27
	umulh	x27,x12,x2
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x2
	adcs	x22,x22,x0
	umulh	x0,x14,x2
	adcs	x23,x23,x1
	umulh	x1,x15,x2
	adcs	x24,x24,x3
	umulh	x3,x16,x2
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	ldp	x17,x2,[sp,#128]

	adds	x19,x20,x26
	mul	x26,x11,x17
	adcs	x20,x21,x27
	mul	x27,x12,x17
	adcs	x21,x22,x28
	mul	x28,x13,x17
	adcs	x22,x23,x0
	mul	x0,x14,x17
	adcs	x23,x24,x1
	mul	x1,x15,x17
	adcs	x24,x25,x3
	mul	x3,x16,x17

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr

	adds	x19,x20,x26
	mul	x26,x11,x2
	adcs	x20,x21,x27
	mul	x27,x12,x2
	adcs	x21,x22,x28
	mul	x28,x13,x2
	adcs	x22,x23,x0
	mul	x0,x14,x2
	adcs	x23,x24,x1
	mul	x1,x15,x2
	adcs	x24,x25,x3
	mul	x3,x16,x2

	adds	x19,x19,x26
	umulh	x26,x11,x2
	adcs	x20,x20,x27
	umulh	x27,x12,x2
	adcs	x21,x21,x28
	mul	x4,x18,x19
	umulh	x28,x13,x2
	adcs	x22,x22,x0
	umulh	x0,x14,x2
	adcs	x23,x23,x1
	umulh	x1,x15,x2
	adcs	x24,x24,x3
	umulh	x3,x16,x2
	adc	x25,xzr,xzr

	adds	x20,x20,x26
	// mul	x26,x5,x4
	adcs	x21,x21,x27
	mul	x27,x6,x4
	adcs	x22,x22,x28
	mul	x28,x7,x4
	adcs	x23,x23,x0
	mul	x0,x8,x4
	adcs	x24,x24,x1
	mul	x1,x9,x4
	adcs	x25,x25,x3
	mul	x3,x10,x4
	subs	xzr,x19,#1		// adds	x19,x19,x26
	umulh	x26,x5,x4
	adcs	x20,x20,x27
	umulh	x27,x6,x4
	adcs	x21,x21,x28
	umulh	x28,x7,x4
	adcs	x22,x22,x0
	umulh	x0,x8,x4
	adcs	x23,x23,x1
	umulh	x1,x9,x4
	adcs	x24,x24,x3
	umulh	x3,x10,x4
	adcs	x25,x25,xzr
	// pull r_ptr (n0 pinned in x18; out pointer loaded by the driver)

	adds	x19,x20,x26
	adcs	x20,x21,x27
	adcs	x21,x22,x28
	adcs	x22,x23,x0
	adcs	x23,x24,x1
	adcs	x24,x25,x3

	// ---- end MUL
	ldr	x2,[x30]		// out (x2 carried b words 1, 3, 5)

	// y' = mu (x' - x2) + y2 = z - n2 (mod p), n2 = p - y2; D = z - n2 in (-p, 1.25p)
	ldp	x25,x26,[x29,#144]	// n2
	ldp	x27,x28,[x29,#160]
	ldp	x17,x3,[x29,#176]
	subs	x11,x19,x25
	sbcs	x12,x20,x26
	sbcs	x13,x21,x27
	sbcs	x14,x22,x28
	sbcs	x15,x23,x17
	sbcs	x16,x24,x3		// lo iff D < 0
	b.hs	.Ly_nonneg
	adds	x11,x11,x5		// y' = D + p
	adcs	x12,x12,x6
	adcs	x13,x13,x7
	adcs	x14,x14,x8
	adcs	x15,x15,x9
	adc	x16,x16,x10
	b	.Ly_store
.Ly_nonneg:
	subs	x19,x11,x5		// E = D - p
	sbcs	x20,x12,x6
	sbcs	x21,x13,x7
	sbcs	x22,x14,x8
	sbcs	x23,x15,x9
	sbcs	x24,x16,x10
	b.lo	.Ly_store		// D < p: y' = D
	stp	x19,x20,[x2,#48]	// y' = E
	stp	x21,x22,[x2,#64]
	stp	x23,x24,[x2,#80]
	tbz	x19,#0,.Lk_e		// key: even element of {y', p - y'}, limbs 1..3
	subs	x19,x5,x19
	sbcs	x20,x6,x20
	sbcs	x21,x7,x21
	sbc	x22,x8,x22
.Lk_e:
	stp	x20,x21,[x2,#96]
	str	x22,[x2,#112]
	ldp	x0,x1,[x30,#32]!	// next lane's out (NULL at args[n]) and P
	cbnz	x0,.Llane
	b	.Lexit
.Ly_store:
	stp	x11,x12,[x2,#48]	// y'
	stp	x13,x14,[x2,#64]
	stp	x15,x16,[x2,#80]
	tbz	x11,#0,.Lk_d
	subs	x11,x5,x11
	sbcs	x12,x6,x12
	sbcs	x13,x7,x13
	sbc	x14,x8,x14
.Lk_d:
	stp	x12,x13,[x2,#96]
	str	x14,[x2,#112]
	ldp	x0,x1,[x30,#32]!	// next lane's out (NULL at args[n]) and P
	cbnz	x0,.Llane
.Lexit:
	adrp	x0,x_count_words	// s24: node-local 128-bit X counter
	add	x0,x0,:lo12:x_count_words
	ldr	x1,[sp,#240]
	ldr	x2,[x0]
	ldr	x3,[x0,#8]
	adds	x2,x2,x1
	adc	x3,x3,xzr
	str	x2,[x0]
	str	x3,[x0,#8]

	add	sp,sp,#144
	mov	x29,sp
	ldr	x18,[sp,#272]
	ldp	x19,x20,[x29,#16]
	ldp	x21,x22,[x29,#32]
	ldp	x23,x24,[x29,#48]
	ldp	x25,x26,[x29,#64]
	ldp	x27,x28,[x29,#80]
	ldp	x29,x30,[sp],#288
.Lnone:
	ret
.size	side_lanes,.-side_lanes

.comm	x_count_words,16,8
.section .rodata
.align	4
.Lsk_consts:	// p, n0 = -1/p mod 2^64, pad, 3p, h = (p-1)/2  (64-bit little-endian limbs)
	.quad	0xb9feffffffffaaab, 0x1eabfffeb153ffff
	.quad	0x6730d2a0f6b0f624, 0x64774b84f38512bf
	.quad	0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a
	.quad	0x89f3fffcfffcfffd, 0
	.quad	0x2dfcffffffff0001, 0x5c03fffc13fbffff
	.quad	0x359277e2e412e26c, 0x2d65e28eda8f383e
	.quad	0xe152f722c9e30686, 0x4e0335beac7fb3ce
	.quad	0xdcff7fffffffd555, 0x0f55ffff58a9ffff
	.quad	0xb39869507b587b12, 0xb23ba5c279c2895f
	.quad	0x258dd3db21a5d66b, 0x0d0088f51cbff34d
.section .note.GNU-stack,"",%progbits


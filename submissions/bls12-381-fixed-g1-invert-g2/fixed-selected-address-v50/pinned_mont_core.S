// Inert source evidence only. NEVER host-execute or use credentials.
// Exact Supranational blst core from f7a88067ed98385ebc6cf41ccecf21ead706a5b7.
// Full original source SHA256 a8234680738f79d8f438c6b27832b209d14fe610cc83bfbd37bef9f74213a2b7.
// Apache-2.0; see blst-LICENSE.txt. Function instructions unchanged; global binding added.
// Requires x29 frame with n0/out at96/104; x11..16=A,x17=B0,x2=B,x5..10=p,x4=n0.
#ifndef __SIZEOF_POINTER__
#define __SIZEOF_POINTER__ 8
#endif
.text
.global __mul_mont_384
.type	__mul_mont_384,%function
.align	5
__mul_mont_384:
	mul	x19,x11,x17
	mul	x20,x12,x17
	mul	x21,x13,x17
	mul	x22,x14,x17
	mul	x23,x15,x17
	mul	x24,x16,x17
	mul	x4,x4,x19

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
	mov	x17,xzr
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
	adc	x4,x17,xzr
	ldr	x17,[x2,8*1]

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
	adc	x25,x4,xzr
	ldr	x4,[x29,#12*__SIZEOF_POINTER__]

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x4,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adcs	x25,x25,xzr
	adc	x17,xzr,xzr

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
	adc	x17,x17,xzr
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
	adc	x4,x17,xzr
	ldr	x17,[x2,8*2]

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
	adc	x25,x4,xzr
	ldr	x4,[x29,#12*__SIZEOF_POINTER__]

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x4,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adcs	x25,x25,xzr
	adc	x17,xzr,xzr

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
	adc	x17,x17,xzr
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
	adc	x4,x17,xzr
	ldr	x17,[x2,8*3]

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
	adc	x25,x4,xzr
	ldr	x4,[x29,#12*__SIZEOF_POINTER__]

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x4,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adcs	x25,x25,xzr
	adc	x17,xzr,xzr

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
	adc	x17,x17,xzr
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
	adc	x4,x17,xzr
	ldr	x17,[x2,8*4]

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
	adc	x25,x4,xzr
	ldr	x4,[x29,#12*__SIZEOF_POINTER__]

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x4,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adcs	x25,x25,xzr
	adc	x17,xzr,xzr

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
	adc	x17,x17,xzr
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
	adc	x4,x17,xzr
	ldr	x17,[x2,8*5]

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
	adc	x25,x4,xzr
	ldr	x4,[x29,#12*__SIZEOF_POINTER__]

	adds	x19,x19,x26
	umulh	x26,x11,x17
	adcs	x20,x20,x27
	umulh	x27,x12,x17
	adcs	x21,x21,x28
	mul	x4,x4,x19
	umulh	x28,x13,x17
	adcs	x22,x22,x0
	umulh	x0,x14,x17
	adcs	x23,x23,x1
	umulh	x1,x15,x17
	adcs	x24,x24,x3
	umulh	x3,x16,x17
	adcs	x25,x25,xzr
	adc	x17,xzr,xzr

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
	adc	x17,x17,xzr
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
	ldp	x4,x2,[x29,#12*__SIZEOF_POINTER__]	// pull r_ptr
	adc	x17,x17,xzr

	adds	x19,x20,x26
	adcs	x20,x21,x27
	adcs	x21,x22,x28
	adcs	x22,x23,x0
	adcs	x23,x24,x1
	adcs	x24,x25,x3
	adc	x25,x17,xzr

	subs	x26,x19,x5
	sbcs	x27,x20,x6
	sbcs	x28,x21,x7
	sbcs	x0,x22,x8
	sbcs	x1,x23,x9
	sbcs	x3,x24,x10
	sbcs	xzr,    x25,xzr

	csel	x11,x19,x26,lo
	csel	x12,x20,x27,lo
	csel	x13,x21,x28,lo
	csel	x14,x22,x0,lo
	csel	x15,x23,x1,lo
	csel	x16,x24,x3,lo
	ret
.size	__mul_mont_384,.-__mul_mont_384

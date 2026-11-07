// Unchanged needed routine bodies extracted from pinned blst f7a88067ed98385ebc6cf41ccecf21ead706a5b7.
// Original full-file SHA256 a8234680738f79d8f438c6b27832b209d14fe610cc83bfbd37bef9f74213a2b7
// Upstream Apache-2.0 license retained in blst-LICENSE.txt.
#if defined(__ARM_FEATURE_PAC_DEFAULT) && __ARM_FEATURE_PAC_DEFAULT==2
# define PACI_HINT 27
# define AUTI_HINT 31
#else
# define PACI_HINT 25
# define AUTI_HINT 29
#endif

.text


// Original inclusive LF lines 414..453
.globl	mul_mont_384
.hidden	mul_mont_384
.type	mul_mont_384,%function
.align	5
mul_mont_384:
	hint	#PACI_HINT
	stp	x29,x30,[sp,#-16*__SIZEOF_POINTER__]!
	add	x29,sp,#0
	stp	x19,x20,[sp,#2*__SIZEOF_POINTER__]
	stp	x21,x22,[sp,#4*__SIZEOF_POINTER__]
	stp	x23,x24,[sp,#6*__SIZEOF_POINTER__]
	stp	x25,x26,[sp,#8*__SIZEOF_POINTER__]
	stp	x27,x28,[sp,#10*__SIZEOF_POINTER__]
	stp	x4,x0,[sp,#12*__SIZEOF_POINTER__]	// __mul_mont_384 wants them there

	ldp	x11,x12,[x1]
	ldr	x17,        [x2]
	ldp	x13,x14,[x1,#16]
	ldp	x15,x16,[x1,#32]

	ldp	x5,x6,[x3]
	ldp	x7,x8,[x3,#16]
	ldp	x9,x10,[x3,#32]

	bl	__mul_mont_384
	ldr	x30,[x29,#__SIZEOF_POINTER__]

	stp	x11,x12,[x2]
	stp	x13,x14,[x2,#16]
	stp	x15,x16,[x2,#32]

	ldp	x19,x20,[x29,#2*__SIZEOF_POINTER__]
	ldp	x21,x22,[x29,#4*__SIZEOF_POINTER__]
	ldp	x23,x24,[x29,#6*__SIZEOF_POINTER__]
	ldp	x25,x26,[x29,#8*__SIZEOF_POINTER__]
	ldp	x27,x28,[x29,#10*__SIZEOF_POINTER__]
	ldr	x29,[sp],#16*__SIZEOF_POINTER__
	hint	#AUTI_HINT
	ret
.size	mul_mont_384,.-mul_mont_384

// Original inclusive LF lines 455..825
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

// Original inclusive LF lines 1081..1112
.globl	redc_mont_384
.hidden	redc_mont_384
.type	redc_mont_384,%function
.align	5
redc_mont_384:
	hint	#PACI_HINT
	stp	x29,x30,[sp,#-16*__SIZEOF_POINTER__]!
	add	x29,sp,#0
	stp	x19,x20,[sp,#2*__SIZEOF_POINTER__]
	stp	x21,x22,[sp,#4*__SIZEOF_POINTER__]
	stp	x23,x24,[sp,#6*__SIZEOF_POINTER__]
	stp	x25,x26,[sp,#8*__SIZEOF_POINTER__]
	stp	x27,x28,[sp,#10*__SIZEOF_POINTER__]
	mov	x4,x3		// adjust for missing b_ptr

	ldp	x5,x6,[x2]
	ldp	x7,x8,[x2,#16]
	ldp	x9,x10,[x2,#32]

	bl	__mul_by_1_mont_384
	bl	__redc_tail_mont_384
	ldr	x30,[x29,#__SIZEOF_POINTER__]

	ldp	x19,x20,[x29,#2*__SIZEOF_POINTER__]
	ldp	x21,x22,[x29,#4*__SIZEOF_POINTER__]
	ldp	x23,x24,[x29,#6*__SIZEOF_POINTER__]
	ldp	x25,x26,[x29,#8*__SIZEOF_POINTER__]
	ldp	x27,x28,[x29,#10*__SIZEOF_POINTER__]
	ldr	x29,[sp],#16*__SIZEOF_POINTER__
	hint	#AUTI_HINT
	ret
.size	redc_mont_384,.-redc_mont_384

// Original inclusive LF lines 1164..1334
.type	__mul_by_1_mont_384,%function
.align	5
__mul_by_1_mont_384:
	ldp	x11,x12,[x1]
	ldp	x13,x14,[x1,#16]
	mul	x26,x4,x11
	ldp	x15,x16,[x1,#32]

	// mul	x19,x5,x26
	mul	x20,x6,x26
	mul	x21,x7,x26
	mul	x22,x8,x26
	mul	x23,x9,x26
	mul	x24,x10,x26
	subs	xzr,x11,#1		// adds	x19,x19,x11
	umulh	x11,x5,x26
	adcs	x20,x20,x12
	umulh	x12,x6,x26
	adcs	x21,x21,x13
	umulh	x13,x7,x26
	adcs	x22,x22,x14
	umulh	x14,x8,x26
	adcs	x23,x23,x15
	umulh	x15,x9,x26
	adcs	x24,x24,x16
	umulh	x16,x10,x26
	adc	x25,xzr,xzr
	adds	x11,x11,x20
	adcs	x12,x12,x21
	adcs	x13,x13,x22
	mul	x26,x4,x11
	adcs	x14,x14,x23
	adcs	x15,x15,x24
	adc	x16,x16,x25

	// mul	x19,x5,x26
	mul	x20,x6,x26
	mul	x21,x7,x26
	mul	x22,x8,x26
	mul	x23,x9,x26
	mul	x24,x10,x26
	subs	xzr,x11,#1		// adds	x19,x19,x11
	umulh	x11,x5,x26
	adcs	x20,x20,x12
	umulh	x12,x6,x26
	adcs	x21,x21,x13
	umulh	x13,x7,x26
	adcs	x22,x22,x14
	umulh	x14,x8,x26
	adcs	x23,x23,x15
	umulh	x15,x9,x26
	adcs	x24,x24,x16
	umulh	x16,x10,x26
	adc	x25,xzr,xzr
	adds	x11,x11,x20
	adcs	x12,x12,x21
	adcs	x13,x13,x22
	mul	x26,x4,x11
	adcs	x14,x14,x23
	adcs	x15,x15,x24
	adc	x16,x16,x25

	// mul	x19,x5,x26
	mul	x20,x6,x26
	mul	x21,x7,x26
	mul	x22,x8,x26
	mul	x23,x9,x26
	mul	x24,x10,x26
	subs	xzr,x11,#1		// adds	x19,x19,x11
	umulh	x11,x5,x26
	adcs	x20,x20,x12
	umulh	x12,x6,x26
	adcs	x21,x21,x13
	umulh	x13,x7,x26
	adcs	x22,x22,x14
	umulh	x14,x8,x26
	adcs	x23,x23,x15
	umulh	x15,x9,x26
	adcs	x24,x24,x16
	umulh	x16,x10,x26
	adc	x25,xzr,xzr
	adds	x11,x11,x20
	adcs	x12,x12,x21
	adcs	x13,x13,x22
	mul	x26,x4,x11
	adcs	x14,x14,x23
	adcs	x15,x15,x24
	adc	x16,x16,x25

	// mul	x19,x5,x26
	mul	x20,x6,x26
	mul	x21,x7,x26
	mul	x22,x8,x26
	mul	x23,x9,x26
	mul	x24,x10,x26
	subs	xzr,x11,#1		// adds	x19,x19,x11
	umulh	x11,x5,x26
	adcs	x20,x20,x12
	umulh	x12,x6,x26
	adcs	x21,x21,x13
	umulh	x13,x7,x26
	adcs	x22,x22,x14
	umulh	x14,x8,x26
	adcs	x23,x23,x15
	umulh	x15,x9,x26
	adcs	x24,x24,x16
	umulh	x16,x10,x26
	adc	x25,xzr,xzr
	adds	x11,x11,x20
	adcs	x12,x12,x21
	adcs	x13,x13,x22
	mul	x26,x4,x11
	adcs	x14,x14,x23
	adcs	x15,x15,x24
	adc	x16,x16,x25

	// mul	x19,x5,x26
	mul	x20,x6,x26
	mul	x21,x7,x26
	mul	x22,x8,x26
	mul	x23,x9,x26
	mul	x24,x10,x26
	subs	xzr,x11,#1		// adds	x19,x19,x11
	umulh	x11,x5,x26
	adcs	x20,x20,x12
	umulh	x12,x6,x26
	adcs	x21,x21,x13
	umulh	x13,x7,x26
	adcs	x22,x22,x14
	umulh	x14,x8,x26
	adcs	x23,x23,x15
	umulh	x15,x9,x26
	adcs	x24,x24,x16
	umulh	x16,x10,x26
	adc	x25,xzr,xzr
	adds	x11,x11,x20
	adcs	x12,x12,x21
	adcs	x13,x13,x22
	mul	x26,x4,x11
	adcs	x14,x14,x23
	adcs	x15,x15,x24
	adc	x16,x16,x25

	// mul	x19,x5,x26
	mul	x20,x6,x26
	mul	x21,x7,x26
	mul	x22,x8,x26
	mul	x23,x9,x26
	mul	x24,x10,x26
	subs	xzr,x11,#1		// adds	x19,x19,x11
	umulh	x11,x5,x26
	adcs	x20,x20,x12
	umulh	x12,x6,x26
	adcs	x21,x21,x13
	umulh	x13,x7,x26
	adcs	x22,x22,x14
	umulh	x14,x8,x26
	adcs	x23,x23,x15
	umulh	x15,x9,x26
	adcs	x24,x24,x16
	umulh	x16,x10,x26
	adc	x25,xzr,xzr
	adds	x11,x11,x20
	adcs	x12,x12,x21
	adcs	x13,x13,x22
	adcs	x14,x14,x23
	adcs	x15,x15,x24
	adc	x16,x16,x25

	ret
.size	__mul_by_1_mont_384,.-__mul_by_1_mont_384

// Original inclusive LF lines 1336..1371
.type	__redc_tail_mont_384,%function
.align	5
__redc_tail_mont_384:
	ldp	x19,x20,[x1,#48]
	ldp	x21,x22,[x1,#64]
	ldp	x23,x24,[x1,#80]

	adds	x11,x11,x19	// accumulate upper half
	adcs	x12,x12,x20
	adcs	x13,x13,x21
	adcs	x14,x14,x22
	adcs	x15,x15,x23
	adcs	x16,x16,x24
	adc	x25,xzr,xzr

	subs	x19,x11,x5
	sbcs	x20,x12,x6
	sbcs	x21,x13,x7
	sbcs	x22,x14,x8
	sbcs	x23,x15,x9
	sbcs	x24,x16,x10
	sbcs	xzr,x25,xzr

	csel	x11,x11,x19,lo
	csel	x12,x12,x20,lo
	csel	x13,x13,x21,lo
	csel	x14,x14,x22,lo
	csel	x15,x15,x23,lo
	csel	x16,x16,x24,lo

	stp	x11,x12,[x0]
	stp	x13,x14,[x0,#16]
	stp	x15,x16,[x0,#32]

	ret
.size	__redc_tail_mont_384,.-__redc_tail_mont_384

// Original inclusive LF lines 1373..1398
.globl	mul_384
.hidden	mul_384
.type	mul_384,%function
.align	5
mul_384:
	hint	#PACI_HINT
	stp	x29,x30,[sp,#-16*__SIZEOF_POINTER__]!
	add	x29,sp,#0
	stp	x19,x20,[sp,#2*__SIZEOF_POINTER__]
	stp	x21,x22,[sp,#4*__SIZEOF_POINTER__]
	stp	x23,x24,[sp,#6*__SIZEOF_POINTER__]
	stp	x25,x26,[sp,#8*__SIZEOF_POINTER__]
	stp	x27,x28,[sp,#10*__SIZEOF_POINTER__]

	bl	__mul_384
	ldr	x30,[x29,#__SIZEOF_POINTER__]

	ldp	x19,x20,[x29,#2*__SIZEOF_POINTER__]
	ldp	x21,x22,[x29,#4*__SIZEOF_POINTER__]
	ldp	x23,x24,[x29,#6*__SIZEOF_POINTER__]
	ldp	x25,x26,[x29,#8*__SIZEOF_POINTER__]
	ldp	x27,x28,[x29,#10*__SIZEOF_POINTER__]
	ldr	x29,[sp],#16*__SIZEOF_POINTER__
	hint	#AUTI_HINT
	ret
.size	mul_384,.-mul_384

// Original inclusive LF lines 1400..1575
.type	__mul_384,%function
.align	5
__mul_384:
	ldp	x11,x12,[x1]
	ldr	x17,        [x2]
	ldp	x13,x14,[x1,#16]
	ldp	x15,x16,[x1,#32]

	mul	x19,x11,x17
	mul	x20,x12,x17
	mul	x21,x13,x17
	mul	x22,x14,x17
	mul	x23,x15,x17
	mul	x24,x16,x17

	umulh	x5,x11,x17
	umulh	x6,x12,x17
	umulh	x7,x13,x17
	umulh	x8,x14,x17
	umulh	x9,x15,x17
	umulh	x10,x16,x17
	ldr	x17,[x2,8*1]

	str	x19,[x0]
	adds	x19,x20,x5
	mul	x5,x11,x17
	adcs	x20,x21,x6
	mul	x6,x12,x17
	adcs	x21,x22,x7
	mul	x7,x13,x17
	adcs	x22,x23,x8
	mul	x8,x14,x17
	adcs	x23,x24,x9
	mul	x9,x15,x17
	adc	x24,xzr,    x10
	mul	x10,x16,x17
	adds	x19,x19,x5
	umulh	x5,x11,x17
	adcs	x20,x20,x6
	umulh	x6,x12,x17
	adcs	x21,x21,x7
	umulh	x7,x13,x17
	adcs	x22,x22,x8
	umulh	x8,x14,x17
	adcs	x23,x23,x9
	umulh	x9,x15,x17
	adcs	x24,x24,x10
	umulh	x10,x16,x17
	ldr	x17,[x2,#8*(1+1)]
	adc	x25,xzr,xzr

	str	x19,[x0,8*1]
	adds	x19,x20,x5
	mul	x5,x11,x17
	adcs	x20,x21,x6
	mul	x6,x12,x17
	adcs	x21,x22,x7
	mul	x7,x13,x17
	adcs	x22,x23,x8
	mul	x8,x14,x17
	adcs	x23,x24,x9
	mul	x9,x15,x17
	adc	x24,x25,x10
	mul	x10,x16,x17
	adds	x19,x19,x5
	umulh	x5,x11,x17
	adcs	x20,x20,x6
	umulh	x6,x12,x17
	adcs	x21,x21,x7
	umulh	x7,x13,x17
	adcs	x22,x22,x8
	umulh	x8,x14,x17
	adcs	x23,x23,x9
	umulh	x9,x15,x17
	adcs	x24,x24,x10
	umulh	x10,x16,x17
	ldr	x17,[x2,#8*(2+1)]
	adc	x25,xzr,xzr

	str	x19,[x0,8*2]
	adds	x19,x20,x5
	mul	x5,x11,x17
	adcs	x20,x21,x6
	mul	x6,x12,x17
	adcs	x21,x22,x7
	mul	x7,x13,x17
	adcs	x22,x23,x8
	mul	x8,x14,x17
	adcs	x23,x24,x9
	mul	x9,x15,x17
	adc	x24,x25,x10
	mul	x10,x16,x17
	adds	x19,x19,x5
	umulh	x5,x11,x17
	adcs	x20,x20,x6
	umulh	x6,x12,x17
	adcs	x21,x21,x7
	umulh	x7,x13,x17
	adcs	x22,x22,x8
	umulh	x8,x14,x17
	adcs	x23,x23,x9
	umulh	x9,x15,x17
	adcs	x24,x24,x10
	umulh	x10,x16,x17
	ldr	x17,[x2,#8*(3+1)]
	adc	x25,xzr,xzr

	str	x19,[x0,8*3]
	adds	x19,x20,x5
	mul	x5,x11,x17
	adcs	x20,x21,x6
	mul	x6,x12,x17
	adcs	x21,x22,x7
	mul	x7,x13,x17
	adcs	x22,x23,x8
	mul	x8,x14,x17
	adcs	x23,x24,x9
	mul	x9,x15,x17
	adc	x24,x25,x10
	mul	x10,x16,x17
	adds	x19,x19,x5
	umulh	x5,x11,x17
	adcs	x20,x20,x6
	umulh	x6,x12,x17
	adcs	x21,x21,x7
	umulh	x7,x13,x17
	adcs	x22,x22,x8
	umulh	x8,x14,x17
	adcs	x23,x23,x9
	umulh	x9,x15,x17
	adcs	x24,x24,x10
	umulh	x10,x16,x17
	ldr	x17,[x2,#8*(4+1)]
	adc	x25,xzr,xzr

	str	x19,[x0,8*4]
	adds	x19,x20,x5
	mul	x5,x11,x17
	adcs	x20,x21,x6
	mul	x6,x12,x17
	adcs	x21,x22,x7
	mul	x7,x13,x17
	adcs	x22,x23,x8
	mul	x8,x14,x17
	adcs	x23,x24,x9
	mul	x9,x15,x17
	adc	x24,x25,x10
	mul	x10,x16,x17
	adds	x19,x19,x5
	umulh	x5,x11,x17
	adcs	x20,x20,x6
	umulh	x6,x12,x17
	adcs	x21,x21,x7
	umulh	x7,x13,x17
	adcs	x22,x22,x8
	umulh	x8,x14,x17
	adcs	x23,x23,x9
	umulh	x9,x15,x17
	adcs	x24,x24,x10
	umulh	x10,x16,x17
	adc	x25,xzr,xzr

	str	x19,[x0,8*5]
	adds	x19,x20,x5
	adcs	x20,x21,x6
	adcs	x21,x22,x7
	adcs	x22,x23,x8
	adcs	x23,x24,x9
	adc	x24,x25,x10

	stp	x19,x20,[x0,#48]
	stp	x21,x22,[x0,#64]
	stp	x23,x24,[x0,#80]

	ret
.size	__mul_384,.-__mul_384

// Explicit hidden alias to the unchanged162-instruction raw leaf.
// The integrated projection caller handles all callee-save clobbers itself.
.globl fapi_raw_mul_leaf
.hidden fapi_raw_mul_leaf
.type fapi_raw_mul_leaf,%function
.set fapi_raw_mul_leaf,__mul_384

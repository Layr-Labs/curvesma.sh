// Exact pinned blst scalar wrapper; Apache-2.0,see blst-LICENSE.txt.
// Full primary source SHA and commit are in pinned_mont_core.S.
// Fixed ordinary64-bit reference build: PAC/AUTI hints25/29,8-byte pointers.
#ifndef __SIZEOF_POINTER__
#define __SIZEOF_POINTER__ 8
#endif
#define PACI_HINT 25
#define AUTI_HINT 29
.text
.global mul_mont_384
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

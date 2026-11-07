// Unchanged needed routine bodies extracted from pinned blst f7a88067ed98385ebc6cf41ccecf21ead706a5b7.
// Original full-file SHA256 c5a9c56ce29e8dd29d7fdbdeb268953102e9752f453c65d2bec23e07d9f5141b
// Upstream Apache-2.0 license retained in blst-LICENSE.txt.
#if defined(__ARM_FEATURE_PAC_DEFAULT) && __ARM_FEATURE_PAC_DEFAULT==2
# define PACI_HINT 27
# define AUTI_HINT 31
#else
# define PACI_HINT 25
# define AUTI_HINT 29
#endif

.text


// Original inclusive LF lines 11..38
.globl	add_mod_384
.hidden	add_mod_384
.type	add_mod_384,%function
.align	5
add_mod_384:
	hint	#PACI_HINT
	stp	x29,x30,[sp,#-6*__SIZEOF_POINTER__]!
	add	x29,sp,#0
	stp	x19,x20,[sp,#2*__SIZEOF_POINTER__]
	stp	x21,x22,[sp,#4*__SIZEOF_POINTER__]

	ldp	x4,x5,[x3]
	ldp	x6,x7,[x3,#16]
	ldp	x8,x9,[x3,#32]

	bl	__add_mod_384
	ldr	x30,[sp,#__SIZEOF_POINTER__]

	stp	x10,x11,[x0]
	stp	x12,x13,[x0,#16]
	stp	x14,x15,[x0,#32]

	ldp	x19,x20,[x29,#2*__SIZEOF_POINTER__]
	ldp	x21,x22,[x29,#4*__SIZEOF_POINTER__]
	ldr	x29,[sp],#6*__SIZEOF_POINTER__
	hint	#AUTI_HINT
	ret
.size	add_mod_384,.-add_mod_384

// Original inclusive LF lines 40..75
.type	__add_mod_384,%function
.align	5
__add_mod_384:
	ldp	x10,x11,[x1]
	ldp	x16,x17,[x2]
	ldp	x12,x13,[x1,#16]
	ldp	x19,x20,[x2,#16]
	ldp	x14,x15,[x1,#32]
	ldp	x21,x22,[x2,#32]

__add_mod_384_ab_are_loaded:
	adds	x10,x10,x16
	adcs	x11,x11,x17
	adcs	x12,x12,x19
	adcs	x13,x13,x20
	adcs	x14,x14,x21
	adcs	x15,x15,x22
	adc	x3,xzr,xzr

	subs	x16,x10,x4
	sbcs	x17,x11,x5
	sbcs	x19,x12,x6
	sbcs	x20,x13,x7
	sbcs	x21,x14,x8
	sbcs	x22,x15,x9
	sbcs	xzr,x3,xzr

	csel	x10,x10,x16,lo
	csel	x11,x11,x17,lo
	csel	x12,x12,x19,lo
	csel	x13,x13,x20,lo
	csel	x14,x14,x21,lo
	csel	x15,x15,x22,lo

	ret
.size	__add_mod_384,.-__add_mod_384

// Original inclusive LF lines 499..526
.globl	sub_mod_384
.hidden	sub_mod_384
.type	sub_mod_384,%function
.align	5
sub_mod_384:
	hint	#PACI_HINT
	stp	x29,x30,[sp,#-6*__SIZEOF_POINTER__]!
	add	x29,sp,#0
	stp	x19,x20,[sp,#2*__SIZEOF_POINTER__]
	stp	x21,x22,[sp,#4*__SIZEOF_POINTER__]

	ldp	x4,x5,[x3]
	ldp	x6,x7,[x3,#16]
	ldp	x8,x9,[x3,#32]

	bl	__sub_mod_384
	ldr	x30,[sp,#__SIZEOF_POINTER__]

	stp	x10,x11,[x0]
	stp	x12,x13,[x0,#16]
	stp	x14,x15,[x0,#32]

	ldp	x19,x20,[x29,#2*__SIZEOF_POINTER__]
	ldp	x21,x22,[x29,#4*__SIZEOF_POINTER__]
	ldr	x29,[sp],#6*__SIZEOF_POINTER__
	hint	#AUTI_HINT
	ret
.size	sub_mod_384,.-sub_mod_384

// Original inclusive LF lines 528..560
.type	__sub_mod_384,%function
.align	5
__sub_mod_384:
	ldp	x10,x11,[x1]
	ldp	x16,x17,[x2]
	ldp	x12,x13,[x1,#16]
	ldp	x19,x20,[x2,#16]
	ldp	x14,x15,[x1,#32]
	ldp	x21,x22,[x2,#32]

	subs	x10,x10,x16
	sbcs	x11,x11,x17
	sbcs	x12,x12,x19
	sbcs	x13,x13,x20
	sbcs	x14,x14,x21
	sbcs	x15,x15,x22
	sbc	x3,xzr,xzr

	and	x16,x4,x3
	and	x17,x5,x3
	adds	x10,x10,x16
	and	x19,x6,x3
	adcs	x11,x11,x17
	and	x20,x7,x3
	adcs	x12,x12,x19
	and	x21,x8,x3
	adcs	x13,x13,x20
	and	x22,x9,x3
	adcs	x14,x14,x21
	adc	x15,x15,x22

	ret
.size	__sub_mod_384,.-__sub_mod_384

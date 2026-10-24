.file	"deltablock_v15l.c"
	.text
	.globl	step_key_block
	.p2align	2
	.type	step_key_block,@function
step_key_block:
.Lhead:
	mov	x4, x0                          // multiplicand = x3
	mov	x5, x1
	mov	x6, x2
	mov	x7, x3
	ldrb	w22, [x28, #64]                 // u_in = slot of the CURRENT point
	mul	x0, x8, x4
	umulh	x16, x8, x4
	mul	x20, x9, x4
	umulh	x21, x9, x4
	adds	x1, x20, x16
	adc	x17, x21, xzr
	mul	x20, x10, x4
	umulh	x21, x10, x4
	adds	x2, x20, x17
	adc	x19, x21, xzr
	mul	x20, x11, x4
	umulh	x21, x11, x4
	adds	x3, x20, x19
	adc	x12, x21, xzr
	mul	x4, x8, x5
	umulh	x16, x8, x5
	mul	x20, x9, x5
	umulh	x21, x9, x5
	adds	x16, x16, x20
	adc	x17, x21, xzr
	mul	x20, x10, x5
	umulh	x21, x10, x5
	adds	x17, x17, x20
	adc	x19, x21, xzr
	mul	x20, x11, x5
	umulh	x21, x11, x5
	adds	x19, x19, x20
	adc	x21, x21, xzr
	adds	x1, x1, x4
	adcs	x2, x2, x16
	adcs	x3, x3, x17
	adcs	x12, x12, x19
	adc	x13, x21, xzr
	mul	x4, x8, x6
	umulh	x16, x8, x6
	mul	x20, x9, x6
	umulh	x21, x9, x6
	adds	x16, x16, x20
	adc	x17, x21, xzr
	mul	x20, x10, x6
	umulh	x21, x10, x6
	adds	x17, x17, x20
	adc	x19, x21, xzr
	mul	x20, x11, x6
	umulh	x21, x11, x6
	adds	x19, x19, x20
	adc	x21, x21, xzr
	adds	x2, x2, x4
	adcs	x3, x3, x16
	adcs	x12, x12, x17
	adcs	x13, x13, x19
	adc	x14, x21, xzr
	mul	x4, x8, x7
	umulh	x16, x8, x7
	mul	x20, x9, x7
	umulh	x21, x9, x7
	adds	x16, x16, x20
	adc	x17, x21, xzr
	mul	x20, x10, x7
	umulh	x21, x10, x7
	adds	x17, x17, x20
	adc	x19, x21, xzr
	mul	x20, x11, x7
	umulh	x21, x11, x7
	adds	x19, x19, x20
	adc	x21, x21, xzr
	adds	x3, x3, x4
	adcs	x12, x12, x16
	adcs	x13, x13, x17
	adcs	x14, x14, x19
	adc	x15, x21, xzr
	mul	x4, x12, x24
	umulh	x16, x12, x24
	mul	x20, x13, x24
	umulh	x21, x13, x24
	adds	x16, x16, x20
	adc	x17, x21, xzr
	mul	x20, x14, x24
	umulh	x21, x14, x24
	adds	x17, x17, x20
	adc	x19, x21, xzr
	mul	x20, x15, x24
	umulh	x21, x15, x24
	adds	x19, x19, x20
	adc	x21, x21, xzr
	adds	x0, x0, x4
	adcs	x1, x1, x16
	adcs	x2, x2, x17
	adcs	x3, x3, x19
	adc	x12, x21, xzr
	mul	x4, x12, x24
	umulh	x16, x12, x24
	adds	x0, x0, x4
	adcs	x1, x1, x16
	adcs	x2, x2, xzr
	adcs	x3, x3, xzr
	adc	x17, xzr, xzr
	adds	x4, x0, x24
	adcs	x16, x1, xzr
	adcs	x19, x2, xzr
	adcs	x20, x3, xzr
	adc	x21, xzr, xzr
	adds	x17, x17, x21
	csel	x0, x4, x0, ne
	csel	x1, x16, x1, ne
	csel	x2, x19, x2, ne
	csel	x3, x20, x3, ne
	subs	x12, x29, x0                    // p - beta*x3  (x29 = -C = p limb0)
	sbcs	x13, x8, x1                     // x8 = ~0 (p limbs 1,3)
	sbcs	x14, x8, x2
	sbcs	x15, x8, x3
	csel	x16, x24, xzr, cc               // borrow -> subtract C
	subs	x12, x12, x16
	sbcs	x13, x13, xzr
	sbcs	x14, x14, xzr
	sbcs	x15, x15, xzr
	subs	x12, x12, x4                    // - x3
	sbcs	x13, x13, x5
	sbcs	x14, x14, x6
	sbcs	x15, x15, x7
	csel	x16, x24, xzr, cc
	subs	x12, x12, x16
	sbcs	x13, x13, xzr
	sbcs	x14, x14, xzr
	sbcs	x15, x15, xzr
	subs	x16, x4, x0
	sbcs	x17, x5, x1
	sbcs	x19, x6, x2
	sbcs	x20, x7, x3
	csel	x4, x4, x0, cc                  // m = min(c0,c1) low
	csel	x5, x5, x1, cc
	csel	x6, x6, x2, cc
	csel	x7, x7, x3, cc
	csel	x21, xzr, x22, cc                // idx = 0 if c0<=c1 else u_in
	subs	x16, x4, x12
	sbcs	x17, x5, x13
	sbcs	x19, x6, x14
	sbcs	x20, x7, x15
	csel	x4, x4, x12, cc                 // min x -> m
	csel	x5, x5, x13, cc
	csel	x6, x6, x14, cc
	csel	x7, x7, x15, cc
	mov	x22, #2
	csel	x22, x21, x22, cc                // slot = argmin index
	and	x21, x4, #2047                  // j = (min x)[0] & 2047
	add	x21, x21, x21, lsl #1           // 3j
	add	x21, x21, x22                    // 3j + slot
	add	x22, x25, x21, lsl #5            // &TBL[j][slot]
	strb	w22, [x28, #64]                 // persist the slot (spare bytes)
	b	.Lhead


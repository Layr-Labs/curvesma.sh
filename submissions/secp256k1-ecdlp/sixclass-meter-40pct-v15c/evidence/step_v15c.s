	.file	"step_secp_v15c.c"
	.text
	.globl	step_asm                        // -- Begin function step_asm
	.p2align	2
	.type	step_asm,@function
step_asm:                               // @step_asm
// %bb.0:
	cmp	w0, #1
	b.lt	.LBB0_6
// %bb.1:
	sub	sp, sp, #192
	str	x30, [sp, #96]                  // 8-byte Folded Spill
	adrp	x8, :got:A
	adrp	x9, :got:SUF
	stp	x28, x27, [sp, #112]            // 16-byte Folded Spill
	adrp	x10, :got:PROD
	mov	w14, w0
	stp	x26, x25, [sp, #128]            // 16-byte Folded Spill
	stp	x24, x23, [sp, #144]            // 16-byte Folded Spill
	stp	x22, x21, [sp, #160]            // 16-byte Folded Spill
	stp	x20, x19, [sp, #176]            // 16-byte Folded Spill
	ldr	x8, [x8, :got_lo12:A]
	ldr	x9, [x9, :got_lo12:SUF]
	ldr	x10, [x10, :got_lo12:PROD]
	add	x8, x8, #40
	add	x9, x9, #16
	add	x10, x10, #16
	b	.LBB0_3
.LBB0_2:                                //   in Loop: Header=BB0_3 Depth=1
	ldp	x11, x12, [x8, #24]
	cmp	w4, #0
	ldr	x9, [sp, #32]                   // 8-byte Folded Reload
	ldr	x14, [sp, #104]                 // 8-byte Folded Reload
	cset	w13, eq
	stp	x4, x7, [x8, #-40]
	add	x10, x10, #32
	stp	x6, x9, [x8, #-24]
	ldr	x9, [sp, #88]                   // 8-byte Folded Reload
	add	x11, x11, #1
	orr	x12, x12, x13
	stp	x17, x18, [x8, #-8]
	subs	x14, x14, #1
	stp	x16, x15, [x8, #8]
	add	x9, x9, #32
	stp	x11, x12, [x8, #24]
	add	x8, x8, #80
	b.eq	.LBB0_5
.LBB0_3:                                // =>This Inner Loop Header: Depth=1
	str	x14, [sp, #104]                 // 8-byte Folded Spill
	ldp	x15, x14, [x8, #-40]
	adrp	x11, :got:TBL
	str	x14, [sp, #40]                  // 8-byte Folded Spill
	and	x17, x15, #0x7ff
	ldr	x11, [x11, :got_lo12:TBL]
	ldp	x0, x18, [x8, #-24]
	ldp	x13, x16, [x8, #8]
	add	x4, x11, x17, lsl #6
	ldp	x11, x12, [x8, #-8]
	mov	x17, #977                       // =0x3d1
	ldp	x5, x6, [x4]
	stp	x0, x13, [sp, #48]              // 16-byte Folded Spill
	ldp	x7, x19, [x4, #16]
	str	x12, [sp, #80]                  // 8-byte Folded Spill
	stp	x11, x16, [sp, #64]             // 16-byte Folded Spill
	movk	x17, #1, lsl #32
	//APP
	subs	x21, x5, x15
	sbcs	x22, x6, x14
	sbcs	x23, x7, x0
	sbcs	x24, x19, x18
	csel	x20, x17, xzr, lo
	subs	x21, x21, x20
	sbcs	x22, x22, xzr
	sbcs	x23, x23, xzr
	sbcs	x24, x24, xzr
	//NO_APP
	ldp	x5, x6, [x4, #32]
	ldp	x7, x19, [x4, #48]
	//APP
	subs	x14, x5, x11
	sbcs	x0, x6, x12
	sbcs	x1, x7, x13
	sbcs	x2, x19, x16
	csel	x28, x17, xzr, lo
	subs	x14, x14, x28
	sbcs	x0, x0, xzr
	sbcs	x1, x1, xzr
	sbcs	x2, x2, xzr
	//NO_APP
	stp	x2, x1, [sp, #8]                // 16-byte Folded Spill
	ldp	x5, x6, [x10, #-16]
	stp	x0, x14, [sp, #24]              // 16-byte Folded Spill
	ldp	x7, x19, [x10]
	//APP
	// a,b values in registers 
	mul	x28, x5, x21
	umulh	x13, x5, x21
	mul	x27, x6, x21
	umulh	x16, x6, x21
	adds	x30, x27, x13
	adc	x20, x16, xzr
	mul	x27, x7, x21
	umulh	x16, x7, x21
	adds	x14, x27, x20
	adc	x25, x16, xzr
	mul	x27, x19, x21
	umulh	x16, x19, x21
	adds	x11, x27, x25
	adc	x0, x16, xzr
	mul	x12, x5, x22
	umulh	x13, x5, x22
	mul	x27, x6, x22
	umulh	x16, x6, x22
	adds	x13, x13, x27
	adc	x20, x16, xzr
	mul	x27, x7, x22
	umulh	x16, x7, x22
	adds	x20, x20, x27
	adc	x25, x16, xzr
	mul	x27, x19, x22
	umulh	x16, x19, x22
	adds	x25, x25, x27
	adc	x26, x16, xzr
	adds	x30, x30, x12
	adcs	x14, x14, x13
	adcs	x11, x11, x20
	adcs	x0, x0, x25
	adc	x1, x26, xzr
	mul	x12, x5, x23
	umulh	x13, x5, x23
	mul	x27, x6, x23
	umulh	x16, x6, x23
	adds	x13, x13, x27
	adc	x20, x16, xzr
	mul	x27, x7, x23
	umulh	x16, x7, x23
	adds	x20, x20, x27
	adc	x25, x16, xzr
	mul	x27, x19, x23
	umulh	x16, x19, x23
	adds	x25, x25, x27
	adc	x26, x16, xzr
	adds	x14, x14, x12
	adcs	x11, x11, x13
	adcs	x0, x0, x20
	adcs	x1, x1, x25
	adc	x2, x26, xzr
	mul	x12, x5, x24
	umulh	x13, x5, x24
	mul	x27, x6, x24
	umulh	x16, x6, x24
	adds	x13, x13, x27
	adc	x20, x16, xzr
	mul	x27, x7, x24
	umulh	x16, x7, x24
	adds	x20, x20, x27
	adc	x25, x16, xzr
	mul	x27, x19, x24
	umulh	x16, x19, x24
	adds	x25, x25, x27
	adc	x26, x16, xzr
	adds	x11, x11, x12
	adcs	x0, x0, x13
	adcs	x1, x1, x20
	adcs	x2, x2, x25
	adc	x3, x26, xzr
	mul	x12, x0, x17
	umulh	x13, x0, x17
	mul	x27, x1, x17
	umulh	x16, x1, x17
	adds	x13, x13, x27
	adc	x20, x16, xzr
	mul	x27, x2, x17
	umulh	x16, x2, x17
	adds	x20, x20, x27
	adc	x25, x16, xzr
	mul	x27, x3, x17
	umulh	x16, x3, x17
	adds	x25, x25, x27
	adc	x26, x16, xzr
	adds	x28, x28, x12
	adcs	x30, x30, x13
	adcs	x14, x14, x20
	adcs	x11, x11, x25
	adc	x0, x26, xzr
	mul	x12, x0, x17
	umulh	x13, x0, x17
	adds	x28, x28, x12
	adcs	x30, x30, x13
	adcs	x14, x14, xzr
	adcs	x11, x11, xzr
	adc	x20, xzr, xzr
	adds	x12, x28, x17
	adcs	x13, x30, xzr
	adcs	x25, x14, xzr
	adcs	x26, x11, xzr
	adc	x27, xzr, xzr
	adds	x20, x20, x27
	csel	x28, x12, x28, ne
	csel	x30, x13, x30, ne
	csel	x14, x25, x14, ne
	csel	x11, x26, x11, ne

	//NO_APP
	stp	x28, x30, [x10, #-16]
	ldp	x12, x13, [x9, #-16]
	stp	x14, x11, [x10]
	ldp	x16, x0, [x9]
	//APP
	// a,b values in registers 
	mul	x26, x12, x21
	umulh	x5, x12, x21
	mul	x20, x13, x21
	umulh	x25, x13, x21
	adds	x27, x20, x5
	adc	x6, x25, xzr
	mul	x20, x16, x21
	umulh	x25, x16, x21
	adds	x28, x20, x6
	adc	x7, x25, xzr
	mul	x20, x0, x21
	umulh	x25, x0, x21
	adds	x30, x20, x7
	adc	x11, x25, xzr
	mul	x3, x12, x22
	umulh	x5, x12, x22
	mul	x20, x13, x22
	umulh	x25, x13, x22
	adds	x5, x5, x20
	adc	x6, x25, xzr
	mul	x20, x16, x22
	umulh	x25, x16, x22
	adds	x6, x6, x20
	adc	x7, x25, xzr
	mul	x20, x0, x22
	umulh	x25, x0, x22
	adds	x7, x7, x20
	adc	x19, x25, xzr
	adds	x27, x27, x3
	adcs	x28, x28, x5
	adcs	x30, x30, x6
	adcs	x11, x11, x7
	adc	x14, x19, xzr
	mul	x3, x12, x23
	umulh	x5, x12, x23
	mul	x20, x13, x23
	umulh	x25, x13, x23
	adds	x5, x5, x20
	adc	x6, x25, xzr
	mul	x20, x16, x23
	umulh	x25, x16, x23
	adds	x6, x6, x20
	adc	x7, x25, xzr
	mul	x20, x0, x23
	umulh	x25, x0, x23
	adds	x7, x7, x20
	adc	x19, x25, xzr
	adds	x28, x28, x3
	adcs	x30, x30, x5
	adcs	x11, x11, x6
	adcs	x14, x14, x7
	adc	x1, x19, xzr
	mul	x3, x12, x24
	umulh	x5, x12, x24
	mul	x20, x13, x24
	umulh	x25, x13, x24
	adds	x5, x5, x20
	adc	x6, x25, xzr
	mul	x20, x16, x24
	umulh	x25, x16, x24
	adds	x6, x6, x20
	adc	x7, x25, xzr
	mul	x20, x0, x24
	umulh	x25, x0, x24
	adds	x7, x7, x20
	adc	x19, x25, xzr
	adds	x30, x30, x3
	adcs	x11, x11, x5
	adcs	x14, x14, x6
	adcs	x1, x1, x7
	adc	x2, x19, xzr
	mul	x3, x11, x17
	umulh	x5, x11, x17
	mul	x20, x14, x17
	umulh	x25, x14, x17
	adds	x5, x5, x20
	adc	x6, x25, xzr
	mul	x20, x1, x17
	umulh	x25, x1, x17
	adds	x6, x6, x20
	adc	x7, x25, xzr
	mul	x20, x2, x17
	umulh	x25, x2, x17
	adds	x7, x7, x20
	adc	x19, x25, xzr
	adds	x26, x26, x3
	adcs	x27, x27, x5
	adcs	x28, x28, x6
	adcs	x30, x30, x7
	adc	x11, x19, xzr
	mul	x3, x11, x17
	umulh	x5, x11, x17
	adds	x26, x26, x3
	adcs	x27, x27, x5
	adcs	x28, x28, xzr
	adcs	x30, x30, xzr
	adc	x6, xzr, xzr
	adds	x3, x26, x17
	adcs	x5, x27, xzr
	adcs	x7, x28, xzr
	adcs	x19, x30, xzr
	adc	x20, xzr, xzr
	adds	x6, x6, x20
	csel	x26, x3, x26, ne
	csel	x27, x5, x27, ne
	csel	x28, x7, x28, ne
	csel	x30, x19, x30, ne

	//NO_APP
	stp	x26, x27, [x9, #-16]
	ldp	x11, x12, [x10, #-16]
	stp	x28, x30, [x9]
	ldp	x13, x14, [x10]
	str	x9, [sp, #88]                   // 8-byte Folded Spill
	mov	x9, #977                        // =0x3d1
	movk	x9, #1, lsl #32
	//APP
	// a,b values in registers 
	mul	x16, x11, x26
	umulh	x20, x11, x26
	mul	x24, x12, x26
	umulh	x25, x12, x26
	adds	x0, x24, x20
	adc	x21, x25, xzr
	mul	x24, x13, x26
	umulh	x25, x13, x26
	adds	x1, x24, x21
	adc	x22, x25, xzr
	mul	x24, x14, x26
	umulh	x25, x14, x26
	adds	x2, x24, x22
	adc	x3, x25, xzr
	mul	x19, x11, x27
	umulh	x20, x11, x27
	mul	x24, x12, x27
	umulh	x25, x12, x27
	adds	x20, x20, x24
	adc	x21, x25, xzr
	mul	x24, x13, x27
	umulh	x25, x13, x27
	adds	x21, x21, x24
	adc	x22, x25, xzr
	mul	x24, x14, x27
	umulh	x25, x14, x27
	adds	x22, x22, x24
	adc	x23, x25, xzr
	adds	x0, x0, x19
	adcs	x1, x1, x20
	adcs	x2, x2, x21
	adcs	x3, x3, x22
	adc	x5, x23, xzr
	mul	x19, x11, x28
	umulh	x20, x11, x28
	mul	x24, x12, x28
	umulh	x25, x12, x28
	adds	x20, x20, x24
	adc	x21, x25, xzr
	mul	x24, x13, x28
	umulh	x25, x13, x28
	adds	x21, x21, x24
	adc	x22, x25, xzr
	mul	x24, x14, x28
	umulh	x25, x14, x28
	adds	x22, x22, x24
	adc	x23, x25, xzr
	adds	x1, x1, x19
	adcs	x2, x2, x20
	adcs	x3, x3, x21
	adcs	x5, x5, x22
	adc	x6, x23, xzr
	mul	x19, x11, x30
	umulh	x20, x11, x30
	mul	x24, x12, x30
	umulh	x25, x12, x30
	adds	x20, x20, x24
	adc	x21, x25, xzr
	mul	x24, x13, x30
	umulh	x25, x13, x30
	adds	x21, x21, x24
	adc	x22, x25, xzr
	mul	x24, x14, x30
	umulh	x25, x14, x30
	adds	x22, x22, x24
	adc	x23, x25, xzr
	adds	x2, x2, x19
	adcs	x3, x3, x20
	adcs	x5, x5, x21
	adcs	x6, x6, x22
	adc	x7, x23, xzr
	mul	x19, x3, x9
	umulh	x20, x3, x9
	mul	x24, x5, x9
	umulh	x25, x5, x9
	adds	x20, x20, x24
	adc	x21, x25, xzr
	mul	x24, x6, x9
	umulh	x25, x6, x9
	adds	x21, x21, x24
	adc	x22, x25, xzr
	mul	x24, x7, x9
	umulh	x25, x7, x9
	adds	x22, x22, x24
	adc	x23, x25, xzr
	adds	x16, x16, x19
	adcs	x0, x0, x20
	adcs	x1, x1, x21
	adcs	x2, x2, x22
	adc	x3, x23, xzr
	mul	x19, x3, x9
	umulh	x20, x3, x9
	adds	x16, x16, x19
	adcs	x0, x0, x20
	adcs	x1, x1, xzr
	adcs	x2, x2, xzr
	adc	x21, xzr, xzr
	adds	x19, x16, x9
	adcs	x20, x0, xzr
	adcs	x22, x1, xzr
	adcs	x23, x2, xzr
	adc	x24, xzr, xzr
	adds	x21, x21, x24
	csel	x16, x19, x16, ne
	csel	x0, x20, x0, ne
	csel	x1, x22, x1, ne
	csel	x2, x23, x2, ne

	//NO_APP
	ldp	x26, x17, [sp, #24]             // 16-byte Folded Reload
	ldp	x28, x27, [sp, #8]              // 16-byte Folded Reload
	//APP
	// a,b values in registers 
	mul	x22, x17, x16
	umulh	x5, x17, x16
	mul	x20, x26, x16
	umulh	x21, x26, x16
	adds	x23, x20, x5
	adc	x6, x21, xzr
	mul	x20, x27, x16
	umulh	x21, x27, x16
	adds	x24, x20, x6
	adc	x7, x21, xzr
	mul	x20, x28, x16
	umulh	x21, x28, x16
	adds	x25, x20, x7
	adc	x11, x21, xzr
	mul	x3, x17, x0
	umulh	x5, x17, x0
	mul	x20, x26, x0
	umulh	x21, x26, x0
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x27, x0
	umulh	x21, x27, x0
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x28, x0
	umulh	x21, x28, x0
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x23, x23, x3
	adcs	x24, x24, x5
	adcs	x25, x25, x6
	adcs	x11, x11, x7
	adc	x12, x19, xzr
	mul	x3, x17, x1
	umulh	x5, x17, x1
	mul	x20, x26, x1
	umulh	x21, x26, x1
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x27, x1
	umulh	x21, x27, x1
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x28, x1
	umulh	x21, x28, x1
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x24, x24, x3
	adcs	x25, x25, x5
	adcs	x11, x11, x6
	adcs	x12, x12, x7
	adc	x13, x19, xzr
	mul	x3, x17, x2
	umulh	x5, x17, x2
	mul	x20, x26, x2
	umulh	x21, x26, x2
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x27, x2
	umulh	x21, x27, x2
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x28, x2
	umulh	x21, x28, x2
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x25, x25, x3
	adcs	x11, x11, x5
	adcs	x12, x12, x6
	adcs	x13, x13, x7
	adc	x14, x19, xzr
	mul	x3, x11, x9
	umulh	x5, x11, x9
	mul	x20, x12, x9
	umulh	x21, x12, x9
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x13, x9
	umulh	x21, x13, x9
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x14, x9
	umulh	x21, x14, x9
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x22, x22, x3
	adcs	x23, x23, x5
	adcs	x24, x24, x6
	adcs	x25, x25, x7
	adc	x11, x19, xzr
	mul	x3, x11, x9
	umulh	x5, x11, x9
	adds	x22, x22, x3
	adcs	x23, x23, x5
	adcs	x24, x24, xzr
	adcs	x25, x25, xzr
	adc	x6, xzr, xzr
	adds	x3, x22, x9
	adcs	x5, x23, xzr
	adcs	x7, x24, xzr
	adcs	x19, x25, xzr
	adc	x20, xzr, xzr
	adds	x6, x6, x20
	csel	x22, x3, x22, ne
	csel	x23, x5, x23, ne
	csel	x24, x7, x24, ne
	csel	x25, x19, x25, ne

	//NO_APP
	//APP
	// a,b values in registers 
	mul	x11, x22, x22
	umulh	x5, x22, x22
	mul	x20, x23, x22
	umulh	x21, x23, x22
	adds	x12, x20, x5
	adc	x6, x21, xzr
	mul	x20, x24, x22
	umulh	x21, x24, x22
	adds	x13, x20, x6
	adc	x7, x21, xzr
	mul	x20, x25, x22
	umulh	x21, x25, x22
	adds	x14, x20, x7
	adc	x16, x21, xzr
	mul	x3, x22, x23
	umulh	x5, x22, x23
	mul	x20, x23, x23
	umulh	x21, x23, x23
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x24, x23
	umulh	x21, x24, x23
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x25, x23
	umulh	x21, x25, x23
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x12, x12, x3
	adcs	x13, x13, x5
	adcs	x14, x14, x6
	adcs	x16, x16, x7
	adc	x0, x19, xzr
	mul	x3, x22, x24
	umulh	x5, x22, x24
	mul	x20, x23, x24
	umulh	x21, x23, x24
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x24, x24
	umulh	x21, x24, x24
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x25, x24
	umulh	x21, x25, x24
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x13, x13, x3
	adcs	x14, x14, x5
	adcs	x16, x16, x6
	adcs	x0, x0, x7
	adc	x1, x19, xzr
	mul	x3, x22, x25
	umulh	x5, x22, x25
	mul	x20, x23, x25
	umulh	x21, x23, x25
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x24, x25
	umulh	x21, x24, x25
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x25, x25
	umulh	x21, x25, x25
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x14, x14, x3
	adcs	x16, x16, x5
	adcs	x0, x0, x6
	adcs	x1, x1, x7
	adc	x2, x19, xzr
	mul	x3, x16, x9
	umulh	x5, x16, x9
	mul	x20, x0, x9
	umulh	x21, x0, x9
	adds	x5, x5, x20
	adc	x6, x21, xzr
	mul	x20, x1, x9
	umulh	x21, x1, x9
	adds	x6, x6, x20
	adc	x7, x21, xzr
	mul	x20, x2, x9
	umulh	x21, x2, x9
	adds	x7, x7, x20
	adc	x19, x21, xzr
	adds	x11, x11, x3
	adcs	x12, x12, x5
	adcs	x13, x13, x6
	adcs	x14, x14, x7
	adc	x16, x19, xzr
	mul	x3, x16, x9
	umulh	x5, x16, x9
	adds	x11, x11, x3
	adcs	x12, x12, x5
	adcs	x13, x13, xzr
	adcs	x14, x14, xzr
	adc	x6, xzr, xzr
	adds	x3, x11, x9
	adcs	x5, x12, xzr
	adcs	x7, x13, xzr
	adcs	x19, x14, xzr
	adc	x20, xzr, xzr
	adds	x6, x6, x20
	csel	x11, x3, x11, ne
	csel	x12, x5, x12, ne
	csel	x13, x7, x13, ne
	csel	x14, x19, x14, ne

	//NO_APP
	ldp	x9, x17, [sp, #40]              // 16-byte Folded Reload
	mov	x5, #977                        // =0x3d1
	movk	x5, #1, lsl #32
	//APP
	subs	x16, x11, x15
	sbcs	x0, x12, x9
	sbcs	x1, x13, x17
	sbcs	x2, x14, x18
	csel	x3, x5, xzr, lo
	subs	x16, x16, x3
	sbcs	x0, x0, xzr
	sbcs	x1, x1, xzr
	sbcs	x2, x2, xzr
	//NO_APP
	ldp	x11, x12, [x4]
	ldp	x13, x14, [x4, #16]
	//APP
	subs	x4, x16, x11
	sbcs	x7, x0, x12
	sbcs	x6, x1, x13
	sbcs	x19, x2, x14
	csel	x3, x5, xzr, lo
	subs	x4, x4, x3
	sbcs	x7, x7, xzr
	sbcs	x6, x6, xzr
	sbcs	x19, x19, xzr
	//NO_APP
	str	x19, [sp, #32]                  // 8-byte Folded Spill
	//APP
	subs	x11, x15, x4
	sbcs	x12, x9, x7
	sbcs	x13, x17, x6
	sbcs	x14, x18, x19
	csel	x16, x5, xzr, lo
	subs	x11, x11, x16
	sbcs	x12, x12, xzr
	sbcs	x13, x13, xzr
	sbcs	x14, x14, xzr
	//NO_APP
	//APP
	// a,b values in registers 
	mul	x0, x22, x11
	umulh	x20, x22, x11
	mul	x28, x23, x11
	umulh	x30, x23, x11
	adds	x1, x28, x20
	adc	x21, x30, xzr
	mul	x28, x24, x11
	umulh	x30, x24, x11
	adds	x2, x28, x21
	adc	x26, x30, xzr
	mul	x28, x25, x11
	umulh	x30, x25, x11
	adds	x3, x28, x26
	adc	x15, x30, xzr
	mul	x19, x22, x12
	umulh	x20, x22, x12
	mul	x28, x23, x12
	umulh	x30, x23, x12
	adds	x20, x20, x28
	adc	x21, x30, xzr
	mul	x28, x24, x12
	umulh	x30, x24, x12
	adds	x21, x21, x28
	adc	x26, x30, xzr
	mul	x28, x25, x12
	umulh	x30, x25, x12
	adds	x26, x26, x28
	adc	x27, x30, xzr
	adds	x1, x1, x19
	adcs	x2, x2, x20
	adcs	x3, x3, x21
	adcs	x15, x15, x26
	adc	x16, x27, xzr
	mul	x19, x22, x13
	umulh	x20, x22, x13
	mul	x28, x23, x13
	umulh	x30, x23, x13
	adds	x20, x20, x28
	adc	x21, x30, xzr
	mul	x28, x24, x13
	umulh	x30, x24, x13
	adds	x21, x21, x28
	adc	x26, x30, xzr
	mul	x28, x25, x13
	umulh	x30, x25, x13
	adds	x26, x26, x28
	adc	x27, x30, xzr
	adds	x2, x2, x19
	adcs	x3, x3, x20
	adcs	x15, x15, x21
	adcs	x16, x16, x26
	adc	x17, x27, xzr
	mul	x19, x22, x14
	umulh	x20, x22, x14
	mul	x28, x23, x14
	umulh	x30, x23, x14
	adds	x20, x20, x28
	adc	x21, x30, xzr
	mul	x28, x24, x14
	umulh	x30, x24, x14
	adds	x21, x21, x28
	adc	x26, x30, xzr
	mul	x28, x25, x14
	umulh	x30, x25, x14
	adds	x26, x26, x28
	adc	x27, x30, xzr
	adds	x3, x3, x19
	adcs	x15, x15, x20
	adcs	x16, x16, x21
	adcs	x17, x17, x26
	adc	x18, x27, xzr
	mul	x19, x15, x5
	umulh	x20, x15, x5
	mul	x28, x16, x5
	umulh	x30, x16, x5
	adds	x20, x20, x28
	adc	x21, x30, xzr
	mul	x28, x17, x5
	umulh	x30, x17, x5
	adds	x21, x21, x28
	adc	x26, x30, xzr
	mul	x28, x18, x5
	umulh	x30, x18, x5
	adds	x26, x26, x28
	adc	x27, x30, xzr
	adds	x0, x0, x19
	adcs	x1, x1, x20
	adcs	x2, x2, x21
	adcs	x3, x3, x26
	adc	x15, x27, xzr
	mul	x19, x15, x5
	umulh	x20, x15, x5
	adds	x0, x0, x19
	adcs	x1, x1, x20
	adcs	x2, x2, xzr
	adcs	x3, x3, xzr
	adc	x21, xzr, xzr
	adds	x19, x0, x5
	adcs	x20, x1, xzr
	adcs	x26, x2, xzr
	adcs	x27, x3, xzr
	adc	x28, xzr, xzr
	adds	x21, x21, x28
	csel	x0, x19, x0, ne
	csel	x1, x20, x1, ne
	csel	x2, x26, x2, ne
	csel	x3, x27, x3, ne

	//NO_APP
	ldp	x12, x9, [sp, #72]              // 16-byte Folded Reload
	ldp	x14, x13, [sp, #56]             // 16-byte Folded Reload
	//APP
	subs	x17, x0, x13
	sbcs	x18, x1, x9
	sbcs	x16, x2, x14
	sbcs	x15, x3, x12
	csel	x11, x5, xzr, lo
	subs	x17, x17, x11
	sbcs	x18, x18, xzr
	sbcs	x16, x16, xzr
	sbcs	x15, x15, xzr
	//NO_APP
	tbz	w17, #0, .LBB0_2
// %bb.4:                               //   in Loop: Header=BB0_3 Depth=1
	mov	x9, #-977                       // =0xfffffffffffffc2f
	mvn	x11, x18
	movk	x9, #65534, lsl #32
	subs	x17, x9, x17
	cset	w12, lo
	subs	x18, x11, x12
	mvn	x11, x16
	cset	w12, lo
	subs	x16, x11, x12
	mvn	x12, x15
	cset	w11, lo
	sub	x15, x12, x11
	b	.LBB0_2
.LBB0_5:
	ldp	x20, x19, [sp, #176]            // 16-byte Folded Reload
	ldr	x30, [sp, #96]                  // 8-byte Folded Reload
	ldp	x22, x21, [sp, #160]            // 16-byte Folded Reload
	ldp	x24, x23, [sp, #144]            // 16-byte Folded Reload
	ldp	x26, x25, [sp, #128]            // 16-byte Folded Reload
	ldp	x28, x27, [sp, #112]            // 16-byte Folded Reload
	add	sp, sp, #192
.LBB0_6:
	ret
.Lfunc_end0:
	.size	step_asm, .Lfunc_end0-step_asm
                                        // -- End function
	.ident	"Apple clang version 21.0.0 (clang-2100.3.34.2)"
	.section	".note.GNU-stack","",@progbits
	.addrsig

	.file	"step_v10f_o3.c"
	.text
	.globl	step_v10f_o3                    // -- Begin function step_v10f_o3
	.p2align	2
	.type	step_v10f_o3,@function
step_v10f_o3:                           // @step_v10f_o3
// %bb.0:
	cmp	w0, #1
	b.lt	.LBB0_4
// %bb.1:
	sub	sp, sp, #272
	str	x30, [sp, #176]                 // 8-byte Folded Spill
	adrp	x8, :got:A
	adrp	x9, :got:PROD
	stp	x28, x27, [sp, #192]            // 16-byte Folded Spill
	adrp	x11, :got:SUF
	stp	x26, x25, [sp, #208]            // 16-byte Folded Spill
	stp	x24, x23, [sp, #224]            // 16-byte Folded Spill
	stp	x22, x21, [sp, #240]            // 16-byte Folded Spill
	stp	x20, x19, [sp, #256]            // 16-byte Folded Spill
	ldr	x8, [x8, :got_lo12:A]
	ldr	x9, [x9, :got_lo12:PROD]
	ldr	x11, [x11, :got_lo12:SUF]
	add	x8, x8, #40
	add	x12, x9, #16
	add	x22, x11, #16
	mov	w9, w0
.LBB0_2:                                // =>This Inner Loop Header: Depth=1
	stp	x9, x22, [sp, #152]             // 16-byte Folded Spill
	ldur	x13, [x8, #-40]
	str	x8, [sp, #184]                  // 8-byte Folded Spill
	adrp	x10, :got:TBL
	mov	w4, #-1                         // =0xffffffff
	ldp	x9, x8, [x8, #24]
	tst	x9, #0x1
	str	x12, [sp, #168]                 // 8-byte Folded Spill
	mov	x23, #-4294967295               // =0xffffffff00000001
	mov	w25, #-1                        // =0xffffffff
	stp	x9, x8, [sp, #136]              // 16-byte Folded Spill
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	mov	x8, #4919131752989213764        // =0x4444444444444444
	eor	x8, x8, #0xe1e1e1e1e1e1e1e1
	ldur	x15, [x9, #-32]
	mov	x9, #31765                      // =0x7c15
	cinv	x8, x8, ne
	movk	x9, #32586, lsl #16
	eor	x8, x8, x13
	movk	x9, #31161, lsl #32
	stp	x15, x13, [sp, #56]             // 16-byte Folded Spill
	movk	x9, #40503, lsl #48
	mul	x8, x8, x9
	mov	x9, #58809                      // =0xe5b9
	movk	x9, #7396, lsl #16
	movk	x9, #18285, lsl #32
	movk	x9, #48984, lsl #48
	eor	x8, x8, x8, lsr #29
	mul	x8, x8, x9
	lsr	x9, x8, #32
	eor	w8, w9, w8
	adrp	x9, :got:TBL2
	ldr	x9, [x9, :got_lo12:TBL2]
	ldr	x10, [x10, :got_lo12:TBL]
	and	x8, x8, #0x7ff
	csel	x9, x10, x9, eq
	add	x8, x9, x8, lsl #6
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	ldp	x11, x14, [x8]
	ldur	x10, [x9, #-24]
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	ldur	x17, [x9, #-16]
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	stp	x14, x11, [sp, #88]             // 16-byte Folded Spill
	subs	x11, x11, x13
	ldp	x13, x0, [x8, #16]
	ldur	x1, [x9, #-8]
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	sbcs	x14, x14, x15
	stp	x17, x10, [sp, #40]             // 16-byte Folded Spill
	ldr	x18, [x9]
	sbcs	x15, x13, x10
	stp	x0, x13, [sp, #72]              // 16-byte Folded Spill
	sbcs	x3, x0, x17
	ldp	x9, x13, [x8, #32]
	stp	x1, x18, [sp, #104]             // 16-byte Folded Spill
	ldp	x16, x8, [x8, #48]
	mrs	x17, NZCV
	mov	w10, #-1                        // =0xffffffff
	csel	x2, x10, xzr, lo
	mov	x10, #-4294967295               // =0xffffffff00000001
	cset	w0, lo
	csel	x26, x10, xzr, lo
	ldr	x10, [sp, #184]                 // 8-byte Folded Reload
	subs	x9, x9, x1
	sbcs	x13, x13, x18
	sbfx	x0, x0, #0, #1
	ldr	x1, [x10, #8]
	ldr	x10, [sp, #184]                 // 8-byte Folded Reload
	ldr	x10, [x10, #16]
	sbcs	x16, x16, x1
	sbcs	x8, x8, x10
	stp	x10, x1, [sp, #120]             // 16-byte Folded Spill
	mrs	x1, NZCV
	cset	w27, lo
	csel	x28, x4, xzr, lo
	csel	x30, x23, xzr, lo
	sbfx	x27, x27, #0, #1
	adds	x27, x9, x27
	adcs	x28, x13, x28
	adcs	x4, x16, xzr
	adc	x30, x30, x8
	msr	NZCV, x1
	csel	x1, x27, x9, lo
	csel	x13, x28, x13, lo
	csel	x24, x4, x16, lo
	csel	x8, x30, x8, lo
	adds	x16, x11, x0
	str	x8, [sp, #8]                    // 8-byte Folded Spill
	adcs	x4, x14, x2
	ldp	x28, x2, [x12, #-16]
	adcs	x27, x15, xzr
	adc	x26, x26, x3
	msr	NZCV, x17
	csel	x0, x16, x11, lo
	ldp	x17, x11, [x12]
	umulh	x5, x0, x28
	csel	x16, x4, x14, lo
	csel	x15, x27, x15, lo
	csel	x14, x26, x3, lo
	mul	x7, x0, x2
	mul	x30, x0, x28
	umulh	x6, x0, x2
	adds	x3, x5, x7
	mul	x20, x0, x17
	umulh	x19, x0, x17
	cinc	x6, x6, hs
	mul	x12, x0, x11
	adds	x6, x6, x20
	lsl	x20, x30, #32
	umulh	x21, x0, x11
	cinc	x19, x19, hs
	sub	x26, x30, x20
	umulh	x7, x30, x23
	adds	x12, x19, x12
	mul	x18, x16, x28
	cinc	x21, x21, hs
	adds	x12, x12, x26
	lsr	x26, x30, #32
	umulh	x4, x16, x28
	cinc	x7, x7, hs
	adds	x3, x20, x3
	adcs	x6, x6, x26
	mul	x5, x16, x2
	adcs	x12, x12, xzr
	umulh	x10, x16, x2
	cinc	x7, x7, hs
	adds	x7, x21, x7
	mul	x27, x16, x17
	cset	w21, hs
	adds	x18, x18, x3
	cinc	x4, x4, hs
	umulh	x19, x16, x17
	adds	x4, x5, x4
	cinc	x10, x10, hs
	adds	x4, x4, x6
	mul	x26, x16, x11
	cinc	x10, x10, hs
	adds	x10, x27, x10
	umulh	x20, x16, x11
	cinc	x19, x19, hs
	adds	x10, x10, x12
	cinc	x19, x19, hs
	umulh	x27, x18, x23
	adds	x19, x26, x19
	lsl	x26, x18, #32
	mul	x5, x15, x28
	cinc	x20, x20, hs
	adds	x7, x19, x7
	lsr	x19, x18, #32
	adc	x20, x20, x21
	adds	x4, x26, x4
	sub	x18, x18, x26
	adcs	x10, x10, x19
	umulh	x3, x15, x28
	adcs	x18, x18, xzr
	mul	x19, x15, x2
	cinc	x26, x27, hs
	adds	x18, x18, x7
	cinc	x7, x26, hs
	umulh	x21, x15, x2
	adds	x7, x20, x7
	cset	w20, hs
	adds	x4, x5, x4
	mul	x5, x15, x17
	cinc	x3, x3, hs
	adds	x3, x19, x3
	umulh	x26, x15, x17
	cinc	x21, x21, hs
	adds	x10, x3, x10
	mul	x19, x15, x11
	cinc	x3, x21, hs
	adds	x3, x5, x3
	umulh	x21, x15, x11
	cinc	x5, x26, hs
	adds	x18, x3, x18
	cinc	x5, x5, hs
	umulh	x26, x4, x23
	adds	x5, x19, x5
	lsl	x19, x4, #32
	mul	x12, x14, x28
	cinc	x21, x21, hs
	adds	x5, x5, x7
	lsr	x7, x4, #32
	adc	x20, x21, x20
	adds	x10, x19, x10
	sub	x4, x4, x19
	adcs	x18, x18, x7
	umulh	x6, x14, x28
	adcs	x4, x4, xzr
	umulh	x3, x14, x2
	cinc	x19, x26, hs
	adds	x4, x4, x5
	cinc	x5, x19, hs
	mul	x2, x14, x2
	adds	x5, x20, x5
	cset	w19, hs
	adds	x10, x12, x10
	umulh	x7, x14, x17
	cinc	x6, x6, hs
	mul	x17, x14, x17
	adds	x2, x2, x6
	mul	x6, x14, x11
	cinc	x3, x3, hs
	adds	x18, x2, x18
	cinc	x2, x3, hs
	umulh	x12, x14, x11
	ldp	x11, x3, [x22, #-16]
	adds	x17, x17, x2
	cinc	x2, x7, hs
	adds	x4, x17, x4
	umulh	x21, x10, x23
	cinc	x17, x2, hs
	lsl	x2, x10, #32
	adds	x17, x6, x17
	umulh	x7, x11, x0
	cinc	x12, x12, hs
	adds	x5, x17, x5
	adc	x12, x12, x19
	lsr	x19, x10, #32
	adds	x18, x2, x18
	sub	x10, x10, x2
	ldp	x17, x2, [x22]
	adcs	x4, x4, x19
	mul	x26, x3, x0
	adcs	x10, x10, xzr
	cinc	x21, x21, hs
	adds	x10, x10, x5
	mul	x20, x11, x0
	cinc	x5, x21, hs
	adds	x12, x12, x5
	umulh	x6, x3, x0
	cset	w27, lo
	adds	x28, x18, #1
	sbcs	x30, x4, x25
	mul	x21, x17, x0
	cset	w22, lo
	subs	x22, x10, x22
	umulh	x19, x17, x0
	sbcs	x8, x12, x23
	cset	w25, lo
	umulh	x5, x2, x0
	tst	w25, w27
	mul	x0, x2, x0
	csel	x27, x18, x28, ne
	csel	x30, x4, x30, ne
	csel	x10, x10, x22, ne
	csel	x22, x12, x8, ne
	adds	x8, x26, x7
	cinc	x18, x6, hs
	umulh	x4, x20, x23
	lsl	x7, x20, #32
	adds	x18, x18, x21
	mov	w21, #-1                        // =0xffffffff
	stp	x22, x10, [sp, #16]             // 16-byte Folded Spill
	cinc	x6, x19, hs
	sub	x19, x20, x7
	mul	x9, x11, x16
	adds	x0, x6, x0
	cinc	x5, x5, hs
	adds	x0, x0, x19
	lsr	x19, x20, #32
	cinc	x4, x4, hs
	adds	x8, x7, x8
	umulh	x25, x11, x16
	adcs	x18, x18, x19
	mul	x6, x3, x16
	adcs	x0, x0, xzr
	cinc	x4, x4, hs
	umulh	x12, x3, x16
	adds	x4, x5, x4
	cset	w5, hs
	adds	x8, x9, x8
	mul	x7, x17, x16
	umulh	x20, x17, x16
	umulh	x19, x2, x16
	mul	x9, x2, x16
	cinc	x16, x25, hs
	adds	x16, x6, x16
	cinc	x12, x12, hs
	adds	x16, x16, x18
	mul	x6, x15, x11
	cinc	x12, x12, hs
	adds	x12, x7, x12
	umulh	x18, x15, x11
	cinc	x7, x20, hs
	adds	x12, x12, x0
	cinc	x7, x7, hs
	umulh	x20, x8, x23
	adds	x9, x9, x7
	lsl	x7, x8, #32
	cinc	x19, x19, hs
	adds	x9, x9, x4
	lsr	x4, x8, #32
	adc	x5, x19, x5
	adds	x16, x7, x16
	sub	x8, x8, x7
	adcs	x12, x12, x4
	mul	x19, x15, x3
	adcs	x8, x8, xzr
	cinc	x7, x20, hs
	adds	x8, x8, x9
	umulh	x0, x15, x3
	cinc	x9, x7, hs
	adds	x9, x5, x9
	mul	x7, x15, x17
	cset	w5, hs
	adds	x16, x6, x16
	umulh	x4, x15, x17
	cinc	x18, x18, hs
	adds	x18, x19, x18
	umulh	x6, x2, x15
	cinc	x0, x0, hs
	adds	x12, x18, x12
	cinc	x18, x0, hs
	mul	x15, x2, x15
	adds	x18, x7, x18
	cinc	x4, x4, hs
	adds	x8, x18, x8
	cinc	x18, x4, hs
	umulh	x4, x16, x23
	mul	x0, x14, x11
	adds	x15, x15, x18
	lsl	x18, x16, #32
	cinc	x6, x6, hs
	adds	x9, x15, x9
	lsr	x15, x16, #32
	adc	x5, x6, x5
	adds	x12, x18, x12
	sub	x16, x16, x18
	adcs	x8, x8, x15
	umulh	x11, x14, x11
	adcs	x16, x16, xzr
	mul	x15, x14, x3
	cinc	x18, x4, hs
	adds	x9, x16, x9
	cinc	x16, x18, hs
	umulh	x6, x14, x3
	adds	x16, x5, x16
	cset	w3, hs
	adds	x12, x0, x12
	umulh	x18, x14, x17
	cinc	x11, x11, hs
	adds	x11, x15, x11
	mul	x17, x14, x17
	cinc	x0, x6, hs
	adds	x8, x11, x8
	mul	x15, x14, x2
	umulh	x11, x14, x2
	cinc	x14, x0, hs
	adds	x14, x17, x14
	cinc	x17, x18, hs
	adds	x9, x14, x9
	umulh	x18, x12, x23
	cinc	x14, x17, hs
	lsl	x17, x12, #32
	adds	x14, x15, x14
	lsr	x15, x12, #32
	cinc	x11, x11, hs
	adds	x14, x14, x16
	sub	x12, x12, x17
	adc	x11, x11, x3
	adds	x8, x17, x8
	adcs	x9, x9, x15
	adcs	x12, x12, xzr
	cinc	x15, x18, hs
	adds	x12, x12, x14
	cinc	x14, x15, hs
	adds	x11, x11, x14
	cset	w14, lo
	adds	x15, x8, #1
	sbcs	x16, x9, x21
	cset	w17, lo
	subs	x17, x12, x17
	sbcs	x18, x11, x23
	cset	w0, lo
	tst	w0, w14
	mov	x14, x27
	ldr	x0, [sp, #8]                    // 8-byte Folded Reload
	csel	x28, x8, x15, ne
	csel	x16, x9, x16, ne
	csel	x15, x12, x17, ne
	umulh	x2, x28, x27
	str	x27, [sp, #32]                  // 8-byte Folded Spill
	mul	x4, x28, x30
	mul	x8, x28, x27
	umulh	x3, x28, x30
	mul	x6, x28, x10
	umulh	x5, x28, x10
	mul	x19, x28, x22
	umulh	x7, x28, x22
	mul	x9, x16, x27
	umulh	x20, x16, x27
	csel	x27, x11, x18, ne
	adds	x11, x2, x4
	cinc	x17, x3, hs
	lsl	x3, x8, #32
	umulh	x18, x8, x23
	adds	x17, x17, x6
	cinc	x2, x5, hs
	sub	x6, x8, x3
	lsr	x8, x8, #32
	adds	x2, x2, x19
	mul	x4, x16, x30
	cinc	x5, x7, hs
	adds	x2, x2, x6
	umulh	x12, x16, x30
	cinc	x18, x18, hs
	adds	x11, x3, x11
	adcs	x8, x17, x8
	mul	x17, x16, x10
	adcs	x2, x2, xzr
	cinc	x18, x18, hs
	umulh	x6, x16, x10
	adds	x18, x5, x18
	cset	w5, hs
	adds	x9, x9, x11
	mul	x11, x16, x22
	cinc	x7, x20, hs
	adds	x4, x4, x7
	umulh	x3, x16, x22
	cinc	x12, x12, hs
	adds	x8, x4, x8
	cinc	x12, x12, hs
	mul	x7, x15, x14
	adds	x12, x17, x12
	cinc	x17, x6, hs
	adds	x12, x12, x2
	umulh	x6, x9, x23
	cinc	x17, x17, hs
	adds	x11, x11, x17
	lsl	x17, x9, #32
	umulh	x4, x15, x14
	cinc	x3, x3, hs
	adds	x11, x11, x18
	lsr	x18, x9, #32
	adc	x3, x3, x5
	adds	x8, x17, x8
	sub	x9, x9, x17
	adcs	x12, x12, x18
	mul	x5, x15, x30
	adcs	x9, x9, xzr
	cinc	x18, x6, hs
	adds	x9, x9, x11
	umulh	x2, x15, x30
	cinc	x11, x18, hs
	adds	x11, x3, x11
	mul	x18, x15, x10
	cset	w3, hs
	adds	x8, x7, x8
	umulh	x17, x15, x10
	cinc	x4, x4, hs
	adds	x4, x5, x4
	mul	x5, x15, x22
	cinc	x2, x2, hs
	adds	x12, x4, x12
	cinc	x2, x2, hs
	umulh	x6, x15, x22
	adds	x18, x18, x2
	lsl	x2, x8, #32
	cinc	x17, x17, hs
	adds	x9, x18, x9
	cinc	x17, x17, hs
	mul	x4, x27, x14
	adds	x17, x5, x17
	umulh	x5, x8, x23
	cinc	x6, x6, hs
	adds	x11, x17, x11
	lsr	x17, x8, #32
	adc	x3, x6, x3
	adds	x12, x2, x12
	sub	x8, x8, x2
	adcs	x9, x9, x17
	umulh	x18, x27, x14
	mov	w14, #-1                        // =0xffffffff
	adcs	x8, x8, xzr
	mul	x17, x27, x30
	cinc	x2, x5, hs
	adds	x8, x8, x11
	cinc	x11, x2, hs
	umulh	x6, x27, x30
	adds	x11, x3, x11
	cset	w3, hs
	adds	x12, x4, x12
	mul	x4, x27, x10
	cinc	x18, x18, hs
	adds	x17, x17, x18
	umulh	x2, x27, x10
	mov	x10, #-4294967295               // =0xffffffff00000001
	cinc	x5, x6, hs
	adds	x9, x17, x9
	mul	x18, x27, x22
	cinc	x17, x5, hs
	adds	x17, x4, x17
	umulh	x5, x27, x22
	cinc	x2, x2, hs
	adds	x8, x17, x8
	cinc	x17, x2, hs
	umulh	x2, x12, x10
	adds	x17, x18, x17
	lsl	x18, x12, #32
	cinc	x4, x5, hs
	adds	x11, x17, x11
	lsr	x17, x12, #32
	adc	x3, x4, x3
	adds	x9, x18, x9
	sub	x12, x12, x18
	adcs	x8, x8, x17
	adcs	x12, x12, xzr
	cinc	x17, x2, hs
	adds	x11, x12, x11
	cinc	x12, x17, hs
	adds	x12, x3, x12
	cset	w17, lo
	adds	x18, x9, #1
	sbcs	x2, x8, x21
	cset	w3, lo
	subs	x3, x11, x3
	sbcs	x4, x12, x10
	cset	w5, lo
	tst	w5, w17
	csel	x9, x9, x18, ne
	csel	x8, x8, x2, ne
	csel	x17, x11, x3, ne
	umulh	x5, x9, x1
	csel	x11, x12, x4, ne
	mul	x7, x9, x13
	mul	x18, x9, x1
	umulh	x6, x9, x13
	adds	x12, x5, x7
	mul	x20, x9, x24
	lsl	x7, x18, #32
	umulh	x19, x9, x24
	cinc	x4, x6, hs
	umulh	x21, x9, x0
	adds	x4, x4, x20
	sub	x20, x18, x7
	mul	x9, x9, x0
	cinc	x6, x19, hs
	umulh	x5, x18, x10
	lsr	x18, x18, #32
	mul	x2, x8, x1
	adds	x9, x6, x9
	cinc	x19, x21, hs
	adds	x9, x9, x20
	umulh	x22, x8, x1
	cinc	x5, x5, hs
	adds	x12, x7, x12
	adcs	x18, x4, x18
	mul	x6, x8, x13
	adcs	x9, x9, xzr
	umulh	x3, x8, x13
	cinc	x5, x5, hs
	adds	x5, x19, x5
	mul	x4, x8, x24
	cset	w19, hs
	adds	x12, x2, x12
	cinc	x2, x22, hs
	umulh	x20, x8, x24
	adds	x2, x6, x2
	cinc	x3, x3, hs
	adds	x18, x2, x18
	umulh	x7, x8, x0
	cinc	x2, x3, hs
	adds	x2, x4, x2
	mul	x8, x8, x0
	cinc	x4, x20, hs
	adds	x9, x2, x9
	cinc	x4, x4, hs
	umulh	x20, x12, x10
	mul	x6, x17, x1
	adds	x8, x8, x4
	lsl	x4, x12, #32
	cinc	x7, x7, hs
	adds	x8, x8, x5
	lsr	x5, x12, #32
	adc	x7, x7, x19
	adds	x18, x4, x18
	sub	x12, x12, x4
	adcs	x9, x9, x5
	umulh	x3, x17, x1
	adcs	x12, x12, xzr
	mul	x19, x17, x13
	cinc	x5, x20, hs
	adds	x8, x12, x8
	cinc	x12, x5, hs
	umulh	x2, x17, x13
	adds	x12, x7, x12
	cset	w7, hs
	adds	x18, x6, x18
	mul	x5, x17, x24
	cinc	x3, x3, hs
	adds	x3, x19, x3
	umulh	x4, x17, x24
	cinc	x2, x2, hs
	adds	x9, x3, x9
	umulh	x6, x17, x0
	cinc	x2, x2, hs
	adds	x2, x5, x2
	mul	x17, x17, x0
	cinc	x4, x4, hs
	adds	x8, x2, x8
	cinc	x2, x4, hs
	umulh	x4, x18, x10
	mul	x3, x11, x1
	adds	x17, x17, x2
	lsl	x2, x18, #32
	cinc	x5, x6, hs
	adds	x12, x17, x12
	lsr	x17, x18, #32
	adc	x5, x5, x7
	adds	x9, x2, x9
	sub	x18, x18, x2
	adcs	x8, x8, x17
	umulh	x1, x11, x1
	adcs	x17, x18, xzr
	umulh	x6, x11, x13
	cinc	x18, x4, hs
	adds	x12, x17, x12
	cinc	x17, x18, hs
	mul	x13, x11, x13
	adds	x17, x5, x17
	cset	w2, hs
	adds	x9, x3, x9
	mul	x3, x11, x24
	cinc	x1, x1, hs
	umulh	x18, x11, x24
	adds	x13, x13, x1
	mul	x1, x11, x0
	cinc	x4, x6, hs
	adds	x8, x13, x8
	cinc	x13, x4, hs
	umulh	x11, x11, x0
	adds	x13, x3, x13
	mov	x0, #-4294967295                // =0xffffffff00000001
	cinc	x18, x18, hs
	adds	x12, x13, x12
	cinc	x13, x18, hs
	lsl	x18, x9, #32
	umulh	x3, x9, x10
	adds	x13, x1, x13
	lsr	x1, x9, #32
	sub	x9, x9, x18
	cinc	x11, x11, hs
	adds	x13, x13, x17
	adc	x11, x11, x2
	adds	x8, x18, x8
	adcs	x12, x12, x1
	adcs	x9, x9, xzr
	cinc	x17, x3, hs
	adds	x9, x9, x13
	cinc	x13, x17, hs
	adds	x11, x11, x13
	cset	w13, lo
	adds	x17, x8, #1
	sbcs	x18, x12, x14
	mov	w14, #-1                        // =0xffffffff
	cset	w1, lo
	subs	x2, x9, x1
	sbcs	x3, x11, x10
	cset	w1, lo
	tst	w1, w13
	csel	x1, x8, x17, ne
	csel	x13, x12, x18, ne
	csel	x9, x9, x2, ne
	umulh	x17, x1, x1
	csel	x8, x11, x3, ne
	mul	x12, x13, x1
	mul	x4, x1, x1
	umulh	x18, x13, x1
	adds	x3, x17, x12
	mul	x2, x9, x1
	lsl	x19, x4, #32
	umulh	x5, x9, x1
	cinc	x6, x18, hs
	sub	x21, x4, x19
	mul	x11, x8, x1
	adds	x6, x6, x2
	umulh	x17, x8, x1
	cinc	x20, x5, hs
	umulh	x7, x4, x10
	lsr	x4, x4, #32
	adds	x20, x20, x11
	cinc	x22, x17, hs
	adds	x20, x20, x21
	umulh	x21, x13, x13
	cinc	x7, x7, hs
	adds	x3, x19, x3
	adcs	x4, x6, x4
	mul	x6, x13, x13
	adcs	x19, x20, xzr
	cinc	x7, x7, hs
	umulh	x20, x9, x13
	adds	x7, x22, x7
	cset	w22, hs
	adds	x12, x12, x3
	mul	x3, x9, x13
	cinc	x18, x18, hs
	adds	x18, x6, x18
	mul	x6, x8, x13
	cinc	x21, x21, hs
	adds	x18, x18, x4
	cinc	x4, x21, hs
	umulh	x21, x8, x13
	adds	x4, x3, x4
	cinc	x23, x20, hs
	adds	x4, x4, x19
	umulh	x24, x12, x10
	cinc	x19, x23, hs
	lsl	x23, x12, #32
	adds	x19, x6, x19
	cinc	x25, x21, hs
	adds	x7, x19, x7
	lsr	x19, x12, #32
	adc	x22, x25, x22
	adds	x18, x23, x18
	sub	x12, x12, x23
	adcs	x4, x4, x19
	adcs	x12, x12, xzr
	cinc	x19, x24, hs
	adds	x12, x12, x7
	cinc	x7, x19, hs
	umulh	x19, x9, x9
	adds	x7, x22, x7
	cset	w22, hs
	adds	x18, x2, x18
	mul	x2, x9, x9
	cinc	x5, x5, hs
	adds	x3, x3, x5
	mul	x5, x8, x9
	cinc	x20, x20, hs
	adds	x3, x3, x4
	cinc	x4, x20, hs
	umulh	x20, x8, x9
	adds	x2, x2, x4
	cinc	x4, x19, hs
	adds	x12, x2, x12
	umulh	x19, x18, x10
	cinc	x2, x4, hs
	lsl	x4, x18, #32
	adds	x2, x5, x2
	cinc	x23, x20, hs
	adds	x2, x2, x7
	lsr	x7, x18, #32
	adc	x22, x23, x22
	adds	x3, x4, x3
	sub	x18, x18, x4
	adcs	x12, x12, x7
	adcs	x18, x18, xzr
	cinc	x4, x19, hs
	adds	x18, x18, x2
	cinc	x2, x4, hs
	adds	x2, x22, x2
	cset	w4, hs
	adds	x11, x11, x3
	mul	x3, x8, x8
	cinc	x17, x17, hs
	adds	x17, x6, x17
	cinc	x6, x21, hs
	adds	x12, x17, x12
	cinc	x17, x6, hs
	umulh	x6, x8, x8
	adds	x17, x5, x17
	cinc	x5, x20, hs
	adds	x17, x17, x18
	cinc	x18, x5, hs
	umulh	x5, x11, x10
	adds	x18, x3, x18
	lsl	x3, x11, #32
	cinc	x6, x6, hs
	adds	x18, x18, x2
	lsr	x2, x11, #32
	adc	x4, x6, x4
	adds	x12, x3, x12
	sub	x11, x11, x3
	adcs	x17, x17, x2
	ldp	x20, x19, [sp, #40]             // 16-byte Folded Reload
	adcs	x11, x11, xzr
	cinc	x2, x5, hs
	adds	x11, x11, x18
	cinc	x18, x2, hs
	adds	x18, x4, x18
	cset	w2, lo
	adds	x3, x12, #1
	sbcs	x4, x17, x14
	cset	w5, lo
	subs	x5, x11, x5
	sbcs	x6, x18, x10
	cset	w7, lo
	tst	w7, w2
	ldp	x7, x10, [sp, #56]              // 16-byte Folded Reload
	csel	x12, x12, x3, ne
	csel	x17, x17, x4, ne
	csel	x11, x11, x5, ne
	csel	x18, x18, x6, ne
	subs	x12, x12, x10
	sbcs	x17, x17, x7
	sbcs	x11, x11, x19
	sbcs	x18, x18, x20
	mrs	x2, NZCV
	cset	w3, lo
	csel	x4, x14, xzr, lo
	csel	x5, x0, xzr, lo
	sbfx	x3, x3, #0, #1
	adds	x3, x12, x3
	adcs	x4, x17, x4
	adcs	x6, x11, xzr
	adc	x5, x5, x18
	msr	NZCV, x2
	ldr	x14, [sp, #96]                  // 8-byte Folded Reload
	csel	x12, x3, x12, lo
	csel	x17, x4, x17, lo
	csel	x11, x6, x11, lo
	csel	x18, x5, x18, lo
	subs	x12, x12, x14
	ldr	x14, [sp, #88]                  // 8-byte Folded Reload
	sbcs	x17, x17, x14
	ldr	x14, [sp, #80]                  // 8-byte Folded Reload
	sbcs	x11, x11, x14
	ldr	x14, [sp, #72]                  // 8-byte Folded Reload
	sbcs	x18, x18, x14
	mov	w14, #-1                        // =0xffffffff
	mrs	x2, NZCV
	cset	w3, lo
	csel	x4, x14, xzr, lo
	csel	x5, x0, xzr, lo
	sbfx	x3, x3, #0, #1
	adds	x3, x12, x3
	adcs	x4, x17, x4
	adcs	x6, x11, xzr
	adc	x5, x5, x18
	msr	NZCV, x2
	csel	x24, x3, x12, lo
	csel	x25, x4, x17, lo
	csel	x23, x6, x11, lo
	csel	x22, x5, x18, lo
	subs	x11, x10, x24
	sbcs	x12, x7, x25
	sbcs	x17, x19, x23
	sbcs	x18, x20, x22
	mrs	x2, NZCV
	cset	w3, lo
	csel	x4, x14, xzr, lo
	csel	x5, x0, xzr, lo
	sbfx	x3, x3, #0, #1
	adds	x3, x11, x3
	adcs	x4, x12, x4
	adcs	x6, x17, xzr
	adc	x5, x5, x18
	msr	NZCV, x2
	csel	x11, x3, x11, lo
	umulh	x3, x11, x1
	csel	x12, x4, x12, lo
	csel	x17, x6, x17, lo
	mul	x19, x11, x13
	mul	x2, x11, x1
	umulh	x7, x11, x13
	mul	x21, x11, x9
	umulh	x20, x11, x9
	mul	x10, x11, x8
	umulh	x26, x11, x8
	csel	x11, x5, x18, lo
	adds	x18, x3, x19
	cinc	x5, x7, hs
	lsl	x19, x2, #32
	umulh	x6, x2, x0
	adds	x5, x5, x21
	cinc	x7, x20, hs
	sub	x21, x2, x19
	lsr	x2, x2, #32
	adds	x10, x7, x10
	mul	x4, x12, x1
	cinc	x20, x26, hs
	adds	x10, x10, x21
	umulh	x14, x12, x1
	cinc	x6, x6, hs
	adds	x18, x19, x18
	adcs	x2, x5, x2
	mul	x7, x12, x13
	adcs	x10, x10, xzr
	umulh	x3, x12, x13
	cinc	x6, x6, hs
	adds	x6, x20, x6
	mul	x5, x12, x9
	cset	w20, hs
	adds	x18, x4, x18
	cinc	x14, x14, hs
	umulh	x21, x12, x9
	adds	x14, x7, x14
	cinc	x3, x3, hs
	adds	x14, x14, x2
	umulh	x19, x12, x8
	cinc	x2, x3, hs
	adds	x2, x5, x2
	mul	x12, x12, x8
	cinc	x5, x21, hs
	adds	x10, x2, x10
	cinc	x5, x5, hs
	umulh	x7, x18, x0
	mul	x4, x17, x1
	adds	x12, x12, x5
	lsl	x5, x18, #32
	cinc	x19, x19, hs
	adds	x12, x12, x6
	lsr	x6, x18, #32
	adc	x19, x19, x20
	adds	x14, x5, x14
	sub	x18, x18, x5
	adcs	x10, x10, x6
	umulh	x3, x17, x1
	adcs	x18, x18, xzr
	mul	x20, x17, x13
	cinc	x6, x7, hs
	adds	x12, x18, x12
	cinc	x18, x6, hs
	umulh	x2, x17, x13
	adds	x18, x19, x18
	cset	w7, hs
	adds	x14, x4, x14
	mul	x6, x17, x9
	cinc	x3, x3, hs
	adds	x3, x20, x3
	umulh	x5, x17, x9
	cinc	x2, x2, hs
	adds	x10, x3, x10
	umulh	x4, x17, x8
	cinc	x2, x2, hs
	adds	x2, x6, x2
	mul	x17, x17, x8
	cinc	x5, x5, hs
	adds	x12, x2, x12
	cinc	x2, x5, hs
	umulh	x5, x14, x0
	mul	x3, x11, x1
	adds	x17, x17, x2
	lsl	x2, x14, #32
	cinc	x4, x4, hs
	adds	x17, x17, x18
	lsr	x18, x14, #32
	adc	x4, x4, x7
	adds	x10, x2, x10
	sub	x14, x14, x2
	adcs	x12, x12, x18
	umulh	x1, x11, x1
	adcs	x14, x14, xzr
	umulh	x6, x11, x13
	cinc	x18, x5, hs
	adds	x14, x14, x17
	cinc	x18, x18, hs
	mul	x13, x11, x13
	adds	x18, x4, x18
	cset	w2, hs
	adds	x10, x3, x10
	umulh	x17, x11, x9
	cinc	x1, x1, hs
	mul	x9, x11, x9
	adds	x13, x13, x1
	mul	x1, x11, x8
	cinc	x3, x6, hs
	adds	x12, x13, x12
	lsl	x13, x10, #32
	umulh	x8, x11, x8
	cinc	x11, x3, hs
	adds	x9, x9, x11
	cinc	x11, x17, hs
	adds	x9, x9, x14
	umulh	x14, x10, x0
	cinc	x11, x11, hs
	lsr	x17, x10, #32
	sub	x10, x10, x13
	adds	x11, x1, x11
	cinc	x8, x8, hs
	adds	x11, x11, x18
	adc	x8, x8, x2
	adds	x12, x13, x12
	mov	w2, #-1                         // =0xffffffff
	adcs	x9, x9, x17
	adcs	x10, x10, xzr
	cinc	x13, x14, hs
	adds	x10, x10, x11
	cinc	x11, x13, hs
	adds	x8, x8, x11
	cset	w11, lo
	adds	x13, x12, #1
	sbcs	x14, x9, x2
	cset	w17, lo
	subs	x17, x10, x17
	sbcs	x18, x8, x0
	cset	w1, lo
	tst	w1, w11
	csel	x11, x12, x13, ne
	ldp	x12, x13, [sp, #104]            // 16-byte Folded Reload
	csel	x9, x9, x14, ne
	csel	x10, x10, x17, ne
	csel	x8, x8, x18, ne
	subs	x11, x11, x12
	mov	x12, #60675                     // =0xed03
	movk	x12, #53650, lsl #16
	sbcs	x9, x9, x13
	mov	x13, #31765                     // =0x7c15
	movk	x12, #18994, lsl #32
	movk	x13, #32586, lsl #16
	movk	x12, #53685, lsl #48
	movk	x13, #31161, lsl #32
	eor	x12, x24, x12
	movk	x13, #40503, lsl #48
	mul	x12, x12, x13
	ldp	x13, x14, [sp, #120]            // 16-byte Folded Reload
	sbcs	x10, x10, x14
	sbcs	x8, x8, x13
	mrs	x13, NZCV
	ldp	x1, x14, [sp, #160]             // 16-byte Folded Reload
	ldp	x18, x17, [sp, #24]             // 16-byte Folded Reload
	eor	x12, x12, x12, lsr #29
	stur	x28, [x1, #-16]
	stur	x17, [x14, #-16]
	ldp	x1, x14, [sp, #160]             // 16-byte Folded Reload
	stur	x16, [x1, #-8]
	mov	x16, #58809                     // =0xe5b9
	ldp	x1, x17, [sp, #160]             // 16-byte Folded Reload
	stur	x30, [x14, #-8]
	cset	w14, lo
	movk	x16, #7396, lsl #16
	sbfx	x14, x14, #0, #1
	str	x18, [x17]
	ldr	x17, [sp, #168]                 // 8-byte Folded Reload
	ldr	x18, [sp, #16]                  // 8-byte Folded Reload
	movk	x16, #18285, lsl #32
	str	x15, [x1]
	movk	x16, #48984, lsl #48
	str	x18, [x17, #8]
	csel	x17, x2, xzr, lo
	csel	x18, x0, xzr, lo
	adds	x14, x11, x14
	mul	x12, x12, x16
	mov	w2, #-1                         // =0xffffffff
	adcs	x17, x9, x17
	adcs	x16, x10, xzr
	adc	x18, x18, x8
	msr	NZCV, x13
	ldr	x15, [sp, #184]                 // 8-byte Folded Reload
	ldr	x13, [sp, #160]                 // 8-byte Folded Reload
	csel	x9, x17, x9, lo
	csel	x11, x14, x11, lo
	stur	x24, [x15, #-40]
	ldr	x15, [sp, #184]                 // 8-byte Folded Reload
	csel	x10, x16, x10, lo
	str	x27, [x13, #8]
	csel	x8, x18, x8, lo
	subs	x13, x2, x9
	ngcs	x14, x10
	stur	x25, [x15, #-32]
	and	x15, x12, #0x1fffffffffffffff
	cset	w16, lo
	cmp	x15, x12, lsr #32
	ldr	x12, [sp, #184]                 // 8-byte Folded Reload
	ldp	x15, x17, [sp, #136]            // 16-byte Folded Reload
	ldr	x18, [sp, #184]                 // 8-byte Folded Reload
	stur	x23, [x12, #-24]
	sbfx	x12, x16, #0, #1
	ldr	x16, [sp, #184]                 // 8-byte Folded Reload
	add	x15, x15, #1
	stur	x22, [x16, #-16]
	cset	w16, eq
	sub	x12, x12, x8
	str	x15, [x18, #24]
	mov	x15, x11
	orr	x16, x17, x16
	mvn	x17, x11
	add	x12, x12, x0
	//APP
	tst	x11, #0x1
	csel	x15, x15, x17, eq
	csel	x9, x9, x13, eq
	csel	x10, x10, x14, eq
	csel	x8, x8, x12, eq
	//NO_APP
	ldr	x11, [sp, #184]                 // 8-byte Folded Reload
	ldp	x22, x12, [sp, #160]            // 16-byte Folded Reload
	str	x16, [x11, #32]
	ldr	x11, [sp, #184]                 // 8-byte Folded Reload
	stur	x15, [x11, #-8]
	ldr	x11, [sp, #184]                 // 8-byte Folded Reload
	add	x12, x12, #32
	add	x22, x22, #32
	str	x9, [x11]
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	str	x10, [x9, #8]
	ldr	x9, [sp, #184]                  // 8-byte Folded Reload
	str	x8, [x9, #16]
	ldr	x9, [sp, #152]                  // 8-byte Folded Reload
	ldr	x8, [sp, #184]                  // 8-byte Folded Reload
	subs	x9, x9, #1
	add	x8, x8, #80
	b.ne	.LBB0_2
// %bb.3:
	ldp	x20, x19, [sp, #256]            // 16-byte Folded Reload
	ldr	x30, [sp, #176]                 // 8-byte Folded Reload
	ldp	x22, x21, [sp, #240]            // 16-byte Folded Reload
	ldp	x24, x23, [sp, #224]            // 16-byte Folded Reload
	ldp	x26, x25, [sp, #208]            // 16-byte Folded Reload
	ldp	x28, x27, [sp, #192]            // 16-byte Folded Reload
	add	sp, sp, #272
.LBB0_4:
	ret
.Lfunc_end0:
	.size	step_v10f_o3, .Lfunc_end0-step_v10f_o3
                                        // -- End function
	.ident	"Apple clang version 21.0.0 (clang-2100.3.34.2)"
	.section	".note.GNU-stack","",@progbits
	.addrsig

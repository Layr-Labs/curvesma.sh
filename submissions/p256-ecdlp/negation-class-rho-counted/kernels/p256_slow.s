// clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -S p256_walk.c
// Apple clang 21.0.0 (LLVM). Comments and assembler directives stripped; functions: move_to_dbl, move_back, trav_update
move_to_dbl:
	mov w11, #168
	ldr x8, [x4]
	madd x9, x5, x11, x0
	madd x10, x8, x11, x2
	ldr x12, [x9]
	str x10, [x3, x8, lsl #3]
	str x12, [x10]
	ldr x12, [x9, #8]
	str x12, [x10, #8]
	ldr x12, [x9, #16]
	str x12, [x10, #16]
	ldr x12, [x9, #24]
	str x12, [x10, #24]
	ldr x12, [x9, #32]
	str x12, [x10, #32]
	ldr x12, [x9, #40]
	str x12, [x10, #40]
	ldr x12, [x9, #48]
	str x12, [x10, #48]
	ldr x12, [x9, #56]
	str x12, [x10, #56]
	ldr x12, [x9, #64]
	str x12, [x10, #64]
	ldr x12, [x9, #72]
	str x12, [x10, #72]
	ldr x12, [x9, #80]
	str x12, [x10, #80]
	ldr x12, [x9, #88]
	str x12, [x10, #88]
	ldr x12, [x9, #96]
	str x12, [x10, #96]
	ldr x12, [x9, #104]
	str x12, [x10, #104]
	ldr x12, [x9, #112]
	str x12, [x10, #112]
	ldr x12, [x9, #120]
	str x12, [x10, #120]
	ldr x12, [x9, #128]
	str x12, [x10, #128]
	ldr x12, [x9, #136]
	str x12, [x10, #136]
	ldr x12, [x9, #144]
	str x12, [x10, #144]
	ldr x12, [x9, #152]
	str x12, [x10, #152]
	ldr x12, [x9, #160]
	str x12, [x10, #160]
	add x12, x8, #1
	str x12, [x4]
	ldr x12, [x1]
	sub x12, x12, #1
	madd x11, x12, x11, x0
	ldr x13, [x11]
	str x13, [x9]
	ldr x13, [x11, #8]
	str x13, [x9, #8]
	ldr x13, [x11, #16]
	str x13, [x9, #16]
	ldr x13, [x11, #24]
	str x13, [x9, #24]
	ldr x13, [x11, #32]
	str x13, [x9, #32]
	ldr x13, [x11, #40]
	str x13, [x9, #40]
	ldr x13, [x11, #48]
	str x13, [x9, #48]
	ldr x13, [x11, #56]
	str x13, [x9, #56]
	ldr x13, [x11, #64]
	str x13, [x9, #64]
	ldr x13, [x11, #72]
	str x13, [x9, #72]
	ldr x13, [x11, #80]
	str x13, [x9, #80]
	ldr x13, [x11, #88]
	str x13, [x9, #88]
	ldr x13, [x11, #96]
	str x13, [x9, #96]
	ldr x13, [x11, #104]
	str x13, [x9, #104]
	ldr x13, [x11, #112]
	str x13, [x9, #112]
	ldr x13, [x11, #120]
	str x13, [x9, #120]
	ldr x13, [x11, #128]
	str x13, [x9, #128]
	ldr x13, [x11, #136]
	str x13, [x9, #136]
	ldr x13, [x11, #144]
	str x13, [x9, #144]
	ldr x13, [x11, #152]
	str x13, [x9, #152]
	ldr x11, [x11, #160]
	str x11, [x9, #160]
	str x12, [x1]
	ret
.Lfunc_end8:
move_back:
	ldr x8, [x1]
	mov w9, #168
	ldr x10, [x2]
	madd x9, x8, x9, x0
	add x8, x8, #1
	str x10, [x9]
	ldr x10, [x2, #8]
	str x10, [x9, #8]
	ldr x10, [x2, #16]
	str x10, [x9, #16]
	ldr x10, [x2, #24]
	str x10, [x9, #24]
	ldr x10, [x2, #32]
	str x10, [x9, #32]
	ldr x10, [x2, #40]
	str x10, [x9, #40]
	ldr x10, [x2, #48]
	str x10, [x9, #48]
	ldr x10, [x2, #56]
	str x10, [x9, #56]
	ldr x10, [x2, #64]
	str x10, [x9, #64]
	ldr x10, [x2, #72]
	str x10, [x9, #72]
	ldr x10, [x2, #80]
	str x10, [x9, #80]
	ldr x10, [x2, #88]
	str x10, [x9, #88]
	ldr x10, [x2, #96]
	str x10, [x9, #96]
	ldr x10, [x2, #104]
	str x10, [x9, #104]
	ldr x10, [x2, #112]
	str x10, [x9, #112]
	ldr x10, [x2, #120]
	str x10, [x9, #120]
	ldr x10, [x2, #128]
	str x10, [x9, #128]
	ldr x10, [x2, #136]
	str x10, [x9, #136]
	ldr x10, [x2, #144]
	str x10, [x9, #144]
	ldr x10, [x2, #152]
	str x10, [x9, #152]
	ldr x10, [x2, #160]
	str x10, [x9, #160]
	str x8, [x1]
	ret
.Lfunc_end9:
trav_update:
	ldp x9, x11, [x1, #16]
	mov x8, x0
	ldp x13, x12, [x0]
	ldr x10, [x1]
	cmp x9, x12
	cset w14, eq
	cmp x10, x13
	ldr x13, [x0, #16]
	csel w14, wzr, w14, hs
	cmp x9, x12
	csinc w12, w14, wzr, hs
	cmp x11, x13
	csel w12, wzr, w12, ne
	b.lo .LBB10_3
	cbnz w12, .LBB10_3
	ldr x9, [x8, #88]
	add x9, x9, #1
	cmp x9, #64
	str x9, [x8, #88]
	cset w0, hi
	ret
.LBB10_3:
	stp x10, x9, [x8]
	str x11, [x8, #16]
	ldp x10, x9, [x1, #16]
	ldp x11, x12, [x1]
	stp x10, x9, [x8, #40]
	stp x11, x12, [x8, #24]
	ldp x10, x9, [x1, #48]
	ldp x11, x12, [x1, #32]
	stp x10, x9, [x8, #72]
	stp x11, x12, [x8, #56]
	ldr x9, [x8, #88]
	add x9, x9, #1
	cmp x9, #64
	str x9, [x8, #88]
	cset w0, hi
	ret
.Lfunc_end10:

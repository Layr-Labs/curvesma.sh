// Static untrusted evidence; never host-execute or use credentials.
// fapi_fr_add_be(out32,a32,b32); canonical big-endian Fr encodings.
// Native add leaf below from blst f7a88067ed98385ebc6cf41ccecf21ead706a5b7,
// build/elf/add_mod_256-armv8.S lines15..43 (function name changed only).
// Apache-2.0; upstream license retained in blst-LICENSE.txt.
.text
.p2align 2
.global fapi_fr_add_be
.type fapi_fr_add_be,%function
fapi_fr_add_be:
    stp x19,x30,[sp,#-128]!
    mov x19,x0
    ldp x4,x5,[x1,#0]
    rev x4,x4
    rev x5,x5
    stp x5,x4,[sp,#32]
    ldp x4,x5,[x1,#16]
    rev x4,x4
    rev x5,x5
    stp x5,x4,[sp,#16]
    ldp x4,x5,[x2,#0]
    rev x4,x4
    rev x5,x5
    stp x5,x4,[sp,#64]
    ldp x4,x5,[x2,#16]
    rev x4,x4
    rev x5,x5
    stp x5,x4,[sp,#48]
    add x0,sp,#80
    add x1,sp,#16
    add x2,sp,#48
    adrp x3,.Lfr_r
    add x3,x3,:lo12:.Lfr_r
    bl .Lfapi_add256
    ldp x4,x5,[sp,#80]
    rev x4,x4
    rev x5,x5
    stp x5,x4,[x19,#16]
    ldp x4,x5,[sp,#96]
    rev x4,x4
    rev x5,x5
    stp x5,x4,[x19,#0]
    ldp x19,x30,[sp],#128
    ret
.size fapi_fr_add_be,.-fapi_fr_add_be
.type .Lfapi_add256,%function
.Lfapi_add256:
	hint	#34
	ldp	x8,x9,[x1]
	ldp	x12,x13,[x2]

	ldp	x10,x11,[x1,#16]
	adds	x8,x8,x12
	ldp	x14,x15,[x2,#16]
	adcs	x9,x9,x13
	ldp	x4,x5,[x3]
	adcs	x10,x10,x14
	ldp	x6,x7,[x3,#16]
	adcs	x11,x11,x15
	adc	x3,xzr,xzr

	subs	x16,x8,x4
	sbcs	x17,x9,x5
	sbcs	x1,x10,x6
	sbcs	x2,x11,x7
	sbcs	xzr,x3,xzr

	csel	x8,x8,x16,lo
	csel	x9,x9,x17,lo
	csel	x10,x10,x1,lo
	stp	x8,x9,[x0]
	csel	x11,x11,x2,lo
	stp	x10,x11,[x0,#16]

	ret
.size .Lfapi_add256,.-.Lfapi_add256
.section .rodata
.p2align 3
.Lfr_r:
    .quad 0xffffffff00000001,0x53bda402fffe5bfe,0x3339d80809a1d805,0x73eda753299d7d48

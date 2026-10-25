// The arrival rows of proof.md section 3, by pairs of walks (after kaankolcu's a9297fb0): per pair one LDP of
// the two K_cp limb-0 words, 30 instructions per walk on either parity path, one STP of the two table addresses
// and one loop branch: 31.5 per arrival, plus the share row.  A batch has n = 1 (mod 8) walks: walk 0 runs
// alone, then (n - 1)/2 pairs.  Walk A of a pair has its K_cp word in x20 and its address in x21, walk B in x0
// and x30.  The rare path (distinguished class, fallback window or filter hit) is the analytic handler of
// section 3; here B.NE skips an inline test stub that logs which of the three it is (1, 2 or 3) at the log
// whose base the caller leaves at [sp] (x5, x6, x29 are scratch there), so the plain path keeps one branch.
// arrive(x1 = walk slots (x, y), x2 = pointer array, placed so that &tp[n] is its first address with bit 16 set
// (bit 19 in production), x3 = K_cp limb-0 words, x4 = table base b0 (2^29-aligned), x5 = constants G1..G5,
// p0 + 1, 2^64 - 1, 33)
	.macro COMMON s
	ldp x24, x25, [x1, #16]        // x2, x3 limbs of x; F0 = x25
	ldp x17, x19, [x1, #32]        // y0, y1
	ldp x22, x23, [x1], #64        // x0, x1 limbs of x
	umulh x26, x22, x10            // hi(x0*G4)
	umulh x27, x23, x9             // hi(x1*G3)
	umulh x28, x24, x8             // hi(x2*G2)
	umulh x29, x25, x7             // hi(x3*G1)
	add x26, x26, x27
	add x28, x28, x29
	add x26, x26, x28
	madd x26, x22, x11, x26        // + lo(x0*G5)
	madd x26, x23, x10, x26        // + lo(x1*G4)
	madd x26, x24, x9, x26         // + lo(x2*G3)
	madd x26, x25, x8, x26         // + lo(x3*G2): F1
	add x27, x25, x26              // S = F0 + F1
	subs x28, x25, x26             // d01 = F0 - F1; HS iff F0 >= F1
	csel x29, x26, x25, hs         // m = min(F0, F1)
	csel \s, x14, x4, hs          // slot j = 1 or 0
	adds x22, x29, x27             // e = m + S; HS iff m >= 2^64 - S
	csel \s, x15, \s, hs         // slot j = 2
	.endm
	.macro CHAIN kc
	ccmn x26, #16, #2, ne
	ccmn x27, #17, #2, lo
	ccmn x22, x16, #2, lo
	ccmp x28, #16, #0, lo
	ccmp x23, \kc, #4, hs
	.endm
	.macro STUB k1, o
	sub x5, x2, x6
	ldr x29, [sp, #96]
	add x5, x29, x5
	mov x29, #1
	cbz \k1, 1f
	mov x29, #2
	cmn x26, #16
	ccmn x27, #17, #2, lo
	ccmn x22, x16, #2, lo
	ccmp x28, #16, #0, lo
	b.lo 1f
	mov x29, #3
1:
	str x29, [x5, #\o]
	.endm
	.macro EVEN s, kc, o
	cmp x19, #0                    // K1 = y1: Z iff distinguished
	orr x23, x17, #1               // filter word K0 + 1
	CHAIN \kc
	b.ne 5f
	STUB x19, \o
5:
	bfi \s, x17, #8, #21           // bit 8 = s = 0, bits 9..28 = h
	.endm
	.macro ODD s, kc, o
	subs x23, x12, x17             // K0 + 1 = (p0 + 1) - y0
	sbcs x24, x13, x19             // K1 = limb 1 of p - y: Z iff distinguished
	CHAIN \kc
	b.ne 5f
	STUB x24, \o
5:
	bfi \s, x23, #8, #21           // bit 8 = s = 1, bits 9..28 = h
	.endm
	.text
	.globl arrive
	.p2align 2
arrive:
	stp x19, x20, [sp, #-96]!
	stp x21, x22, [sp, #16]
	stp x23, x24, [sp, #32]
	stp x25, x26, [sp, #48]
	stp x27, x28, [sp, #64]
	stp x29, x30, [sp, #80]
	ldp x7, x8, [x5]
	ldp x9, x10, [x5, #16]
	ldr x11, [x5, #32]
	ldp x12, x13, [x5, #40]
	ldr x16, [x5, #56]
	add x14, x4, #64
	add x15, x4, #128
	mov x6, x2
	ldr x20, [x3], #8             // walk 0 alone
	COMMON x30
	tbnz x17, #0, w0odd
	EVEN x30, x20, 0
	b w0end
w0odd:
	ODD x30, x20, 0
w0end:
	str x30, [x2], #8
	tbnz x2, #16, done
pair:
	ldp x20, x0, [x3], #16        // the pair's K_cp words: walk A x20, walk B x0
	COMMON x21
	tbnz x17, #0, aodd
	EVEN x21, x20, 0
	COMMON x30
	tbnz x17, #0, aebodd
	EVEN x30, x0, 8
	stp x21, x30, [x2], #16
	tbz x2, #16, pair
	b done
aebodd:
	ODD x30, x0, 8
	stp x21, x30, [x2], #16
	tbz x2, #16, pair
	b done
aodd:
	ODD x21, x20, 0
	COMMON x30
	tbnz x17, #0, aobodd
	EVEN x30, x0, 8
	stp x21, x30, [x2], #16
	tbz x2, #16, pair
	b done
aobodd:
	ODD x30, x0, 8
	stp x21, x30, [x2], #16
	tbz x2, #16, pair
done:
	ldp x29, x30, [sp, #80]
	ldp x27, x28, [sp, #64]
	ldp x25, x26, [sp, #48]
	ldp x23, x24, [sp, #32]
	ldp x21, x22, [sp, #16]
	ldp x19, x20, [sp], #96
	ret

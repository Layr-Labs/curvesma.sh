# Kernel evidence excerpts (verbatim macros)

The following macros are copied verbatim from the kernel listings of the base `de987687`; the full files are published there. They are the instruction-level accounting for the dominant certificate constants. Counts below are the certificate formulas, not timings.

## Leaf kernels: PFWD=153 per pair, BWDP=2281 per pair

```asm
.macro PFWD
	EFWD
	st64b x20,[x3]
	add x3,x3,#1024
	P1X
	.endm
```

```asm
.macro BWDP
	ld64b x16,[x3]
	W2S
	stp x0,x1,[x3,#-496]
	stp x2,x24,[x3,#-512]
	C1BS
	INC
	sub x3,x3,#512
	ldp x0,x1,[x3,#32]
	ldp x2,x5,[x3,#48]
	INC
	ldp x25,x26,[x3,#16]
	ldp x27,x28,[x3],#-512
	.endm
```

```asm
.macro INC
	C2UV
	SQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
	F5 x21,x22,x23,x24,x8,x21,x22,x9,x8
	SQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
	F5 x25,x26,x0,x2,x8,x25,x26,x9,x8
	PAIRB
	stp x0,x2,[sp,#464]
	str x1,[sp,#480]
	LKE
	TFP4X
	TOUTCX sp,240,256,320,336,128,144
	LKF
	TFP4X
	TOUTCX x3,160,176,288,304,96,112
	LKG
	TFP4X
	TFINS x11,x12
	TFINS x13,x14
	TFINS x15,x16
	TFINS x17,x18
	stp x12,x11,[x3,#256]
	stp x16,x15,[x3,#272]
	stp x14,x13,[x3,#64]
	stp x18,x17,[x3,#80]
	.endm
```

```asm
.macro TOUTCX r,c,d,p,q,m,n
	ldp x9,x10,[\r,#\c]
	ldp x5,x6,[\r,#\d]
	TFINC x11,x12,x9,x10
	TFINC x13,x14,x9,x10
	TFINC x15,x16,x5,x6
	TFINC x17,x18,x5,x6
	stp x12,x11,[x3,#\p]
	stp x16,x15,[x3,#\q]
	stp x14,x13,[x3,#\m]
	stp x18,x17,[x3,#\n]
	.endm
```

```asm
.macro TFP4X
	TROWX TR0,x21,x25,x29,x30,x10,x9,x6,x19
	SHX x19,x5,x10,x6,x9
	TROWX TRK,x22,x26,x27,x28,x9,x6,x5,x10
	SHX x10,x19,x9,x5,x6
	TROWX TRK,x23,x0,x7,x8,x6,x5,x19,x9
	SHXT x10,x6,x19,x5
	TROWX TRK,x24,x2,x1,x20,x5,x19,x10,x6
	.endm
```

```asm
.macro C1BS
	ldp x13,x14,[sp,#16]
	ldp x15,x16,[sp,#32]
	ldp x17,x18,[x3,#192]
	ldp x19,x20,[x3,#208]
	FS x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x21
	RF x9,x13,x14,x15,x16,x0,x1,x2,x5,x6,x7,x8
	SHX x17,x13,x16,x14,x15
	RA4 x10,x17,x13,x14,x15,x0,x1,x2,x5,x6,x18,xzr,x7,x8,x19,x20
	SHX x16,x17,x15,x13,x14
	RA4 x11,x16,x17,x13,x14,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
	SHX x15,x16,x14,x17,x13
	RA4 x12,x15,x16,x17,x13,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
	FDX x0,x1,x6,x18,x0,x1,x7,x8,x2,x5
	stp x0,x1,[x3,#-480]
	stp x2,x5,[x3,#-464]
	ldp x13,x14,[x3,#-128]
	ldp x15,x16,[x3,#-112]
	bfi x13,x9,#8,#32
	ld64b x20,[x13]
	AC4F x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x9,#32
	bfi x13,x7,#8,#32
	ld64b x20,[x13]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x14,x10,#8,#32
	ld64b x20,[x14]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x10,#32
	bfi x14,x7,#8,#32
	ld64b x20,[x14]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x15,x11,#8,#32
	ld64b x20,[x15]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x11,#32
	bfi x15,x7,#8,#32
	ld64b x20,[x15]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x16,x12,#8,#32
	ld64b x20,[x16]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x12,#32
	bfi x16,x7,#8,#32
	ld64b x20,[x16]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	LF x0,x1,x2,x5,x6,x7
	.endm
```

```asm
.macro C2UV
	ldp x13,x14,[sp,#112]
	ldp x15,x16,[sp,#128]
	RF x0,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23
	ldp x13,x14,[sp,#144]
	ldp x15,x16,[sp,#160]
	RA4 x1,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,xzr,x22,x23,x25,x26
	ldp x13,x14,[sp,#176]
	ldp x15,x16,[sp,#192]
	RA4 x2,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,x24,x22,x23,x25,x26
	ldp x13,x14,[sp,#208]
	ldp x15,x16,[sp,#224]
	RA4 x5,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,x24,x22,x23,x25,x26
	FDX x17,x18,x21,x24,x17,x18,x22,x23,x19,x20
	ldp x9,x10,[x3,#352]
	ldp x11,x12,[x3,#368]
	bfi x9,x0,#8,#32
	ldp x13,x14,[x9]
	ldp x15,x16,[x9,#16]
	lsr x7,x0,#32
	bfi x9,x7,#8,#32
	ld64b x22,[x9]
	AC4F x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x10,x1,#8,#32
	ld64b x22,[x10]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x1,#32
	bfi x10,x7,#8,#32
	ld64b x22,[x10]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x11,x2,#8,#32
	ld64b x22,[x11]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x2,#32
	bfi x11,x7,#8,#32
	ld64b x22,[x11]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x12,x5,#8,#32
	ld64b x22,[x12]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x5,#32
	bfi x12,x7,#8,#32
	ld64b x22,[x12]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	LF x13,x14,x15,x16,x6,x7
	FADDW x27,x28,x29,x30,x17,x18,x19,x20,x13,x14,x15,x16,x5
	FS x17,x18,x19,x20,x17,x18,x19,x20,x13,x14,x15,x16,x5
	.endm
```

```asm
.macro LKE
	ldr x1,[x3,#224]
	ldp x5,x6,[sp,#80]
	ldp x9,x10,[sp,#96]
	ldp x0,x2,[sp,#272]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4F x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#288]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#304]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#320]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
```

```asm
.macro LKF
	ldr x1,[x3,#224]
	ldp x0,x2,[sp,#336]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	AC4F x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#352]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#368]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#384]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
```

```asm
.macro LKG
	ldr x1,[x3,#224]
	ldp x0,x2,[sp,#400]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	AC4F x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#416]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#432]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#448]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
```

Certificate counts: `fwd_count(q2)=5+2452*(q2//32)+q2//256` (16 PFWD per 32 incidences, 4 per pass); `bwd_count(q2)=6+18251*(q2//16)+5*((q2+255)//256-1)` (8 BWDP per 16 incidences, 3 per pass). Each PFWD is 153, each BWDP 2281, including the INC/TOUTCX/TFP4X macros above. Per incidence: fwd 76.6289, bwd 1140.7070.

## Field inversion and Y generation: fe_inv=22549; ygen=211,380,571 at (16384,8)

fe_inv is 255 squares and 127 multiplies with the 39479-instruction field-call allowance (certificate `INV=22549`; the multiply/square counts are in proof.md section 3). ygen cost is `A+B*YM+C*YN+D*YM*YN+E*ceil(log2 YM)` with (-23285,5,22598,1606,45680); at (256,16384,8) it is 211,380,571 for the 2^17 points (1612.71 per point), and `YGEN=+3*GEN_MUL=213,316,801` per batch.

## Transport and LS64 delivery: dlv=6, h2s=6/record, dlvh=4, join=4/record, walk=4.5/visit

```asm
.macro HCHGI u,v,a,l
	bfi x10,\u,#3,#58
	swp x13,x11,[x10]
	stp x11,\u,[x13],#16
	cbz x11,8f
9:
	ldp x5,x6,[x11]
	cmp x6,\u
	b.eq 10f
12:
	cbz x5,8f
	ldp x11,x6,[x5]
	cmp x6,\u
	b.eq 11f
13:
	cbnz x11,9b
	b 8f
10:
	add x6,x15,x11,lsl #1
	ldr x12,[x6,#8]
	eor x12,x12,\v
	tst x12,#0x3fffffffffff
	b.ne 12b
	cmp x3,x4
	b.hs 8f
	\l
	stp x6,\a,[x3],#16
	b 8f
11:
	add x6,x15,x5,lsl #1
	ldr x12,[x6,#8]
	eor x12,x12,\v
	tst x12,#0x3fffffffffff
	b.ne 13b
	cmp x3,x4
	b.hs 8f
	\l
	stp x6,\a,[x3],#16
8:
	.endm
```

```asm
.macro H2B a,b,c,d
	bfi x10,\d,#3,#36
	ldadd x17,x11,[x10]
	stp \a,\b,[x11]
	stp \c,\d,[x11,#16]
	tbz x11,#13,9f
	sub x9,x11,#4080
	sub x9,x9,#4080
	str x9,[x10]
	orr x4,x10,#0x80000000000
	ldadd x18,x11,[x4]
	cmp x11,x12
	b.lo 8f
	str x12,[x4]
	b 9f
8:
	stp x18,x19,[x9,#-16]
	stp x10,x11,[x9,#-32]
9:
	.endm
```

```asm
.macro H2P
	ld64b x20,[x1]
	add x1,x1,#64
	H2B x20,x21,x22,x23
	H2B x24,x25,x26,x27
	.endm
```

```asm
.macro DLVP
	ld64b x16,[x1]
	add x1,x1,#64
	st64b x16,[x2]
	add x2,x2,#64
	.endm
```

```asm
.macro DLVS
	ldp x16,x17,[x1],#32
	ldp x18,x19,[x1,#-16]
	stp x18,x19,[x2,#16]
	stp x16,x17,[x2],#32
	.endm
```

```asm
.macro DLVHS
	DLVS
	HCHGI x16,x17,x12,LSUB
	.endm
```

The delivery counts and the in-flight/ACK charges are itemized in proof.md section 4 and certificate.py: dlv 6 per one-record message; h2s 6 per record plus 8 fixed and the flush/copy terms; dlvh 4; join 4 per record; walk 4.5 per visit; the receiver appends each acceptable record and lists a full-key match only across source nodes.

## Captured experiment

`experiment-result.json` is the captured output of `experiment.py` (seed 1, seconds 12) on this exact archive. `counts_match_certificate` is True; 143 boolean checks, 2 false. The full test code (`experiment.py`, `machine.py`, `sites.py`) is in the published base `de987687`.

# Compiled ARMv8.1-A Scalar Kernels (Part 2: Stage 1b Radix-4 NTT, Pure 2-Pass 64-Bit Replay & 4-Basis Expanded LA, Kernels 8–16)

Verbatim companion to `KERNELS.md` compiled with:

```text
clang -O3 -march=armv8.1-a -mgeneral-regs-only -fno-vectorize -fno-slp-vectorize -fomit-frame-pointer -ffixed-x18 -S
```

**Part 2 Audit (Kernels 8–16):** **`942` actual scalar instructions**, **`0` `staddb`**, **`0` forbidden ISA**, **`0` `udiv`/`sdiv`**, **`0` `bl`/`blr`**, **`0` in-loop stack spills ANYWHERE**. Combined with `KERNELS.md` (`583` insts): **`16` verbatim compiled functions, `1525` executable scalar instructions**.

## Part 2 Dynamic Instruction Count Summary (Kernels 8–16)

| # | Symbol | Role | Static Insts | In-Loop Spills | Compiled Dynamic `Ir` | Charged in Ledger |
|---:|---|---|---:|---:|---|---|
| 8 | `_ntt61_r4_stage_dit_u2_exact` | Radix-4 DIT Shoup61 butterfly (`u2` exact, 3 primes) | 97 | **0** | **`89 Ir / 8 r2-bfly = 11.125 Ir/r2-bfly`** | **`11.125 Ir/r2-bfly`** (`1.10x` saf) |
| 9 | `_ntt61_pw_mul_x4` | Pointwise Shoup61 multiply (`x4` unrolled) | 39 | **0** | `37 Ir / 4 words = 9.25 Ir/word` | `12.00 Ir/word` (`+2.75` pad) |
| 10 | `_crt3_1d_carry` | 3-prime 61-bit CRT + 58-bit limb carry | 84 | **0** | **`58.0 Ir/word`** (1D) | `58.0 Ir` (1D) / `90.0 Ir` (2D) |
| 11 | `_replay_build_count_x4` | Pure 2-pass replay histogram count (`x4`, 64-bit `uint64_t` cursor, all 3 tiers) | 38 | **0** | `34/4 = 8.50 Ir/entry` | Charged in `21.50 Ir/ent` |
| 12 | `_replay_build_scatter_off64_x4` | Unified 2-pass 64-bit-cursor `off64 = ((r & mask) << 56) \| a_idx` scatter (`x4`, ALL 3 tiers `W=1,16,32`) | 55 | **0** | `51/4 = 12.75 Ir/entry` (`21.25` total) | **`21.50 Ir/entry + 5.5 Ir/bkt`** (`rep = 2^{104.5395} u` live, `104.5398` ch) |
| 13 | `_spmv_mono_add_ptr64` | **64-bit `LDR X`** signed `+1/-1` 4-basis component SpMV (`>=92%` comp) | 18 | **0** | **`10.00 Ir/comp`** (`LBB13_1` = `10` insns) | Blended **`11.45 Ir/comp`** |
| 14 | `_spmv_nonmono_muladd_ptr64` | **64-bit `LDR X`** signed `|c|<2^{30}` 4-basis component SpMV (`<=8%` comp) | 33 | **0** | **`25.00 Ir/comp`** (`LBB14_1` = `25` insns) | Blended **`11.45 Ir/comp`** |
| 15 | `_spmv_row_horner_sm4_reduce` | Signed bias (`+76 Ir`) + 576-bit `r<<256` Horner+SM4+Barrett (`872 Ir`) | 554 | **0** | `76 + 872 = 948.0 Ir/row` (`A_k < 2^{315}`) | **`960.0 Ir/row`** (`+12` pad) |
| 16 | `_allgather_v256_copy_x4` | 4-owner 256-bit `V` all-gather TX/RX copy (`128B = 16 u64` words/iter) | 24 | **0** | `21/16 = 1.3125 Ir/word` (`2.625` TX+RX) | **`3.00 Ir/word`** (`36*N_l/N_c = 0.1873 Ir/comp`) |

### Production 64-Bit 4-Basis Expanded $\mathbb{Z}[\zeta_{12}]$ SpMV & Banked `Lingen` Audit (`P_LA = 4`, `P_lingen = 1024`)
- **Explicit 4-Basis $\mathbb{Z}[\zeta_{12}]$ Expansion (Zero Representation Gap):** Every post-SGE entry `(col, (c0, c1, c2, c3))` in $\mathbb{Z}[\zeta_{12}]$ is expanded at matrix-build time into its exact nonzero components across the 8 signed streams `(A_k^+, A_k^-)` for `k in {0,1,2,3}`. Monomials like $\zeta_{12}^4 = \zeta_{12}^2 - 1$ explicitly emit **2** `8 B` entries (`+1` in `A_2^+`, `-1` in `A_0^-`) and cost **`20.0 Ir`**; general non-monomial 4-tuples emit up to **4** entries across `A_0..A_3`.
- **Live-Derived Full RSA-768 Pareto-Consistent `H3` (`exploratory-v1`) LA Caps & Compiled `_allgather_v256_copy_x4` Broadcast (`eta = 1.015`):** Zero LA density caps are hardcoded in `REF["la"]`. Instead, `experiment.py` derives `nl_cap`, `comp_cap`, and `pm1_cap` live from: (i) the exact production Dickman degree spectrum and Molloy–Reed 2-core (`E1/U = 0.15610`, sieved-FB fraction `0.002050`), (ii) a live multi-scale PNT `L_CALIB` $\mathbb{Z}[\zeta_{12}]$ Bouvier + `STRICT` SGE ladder (`U = 2^{10}..2^{13}` `core` / `2^{10}..2^{14}` `all`, live OLS $R^2 = 0.908$ / `0.950`) anchored conservatively at **BOTH** published Kleinjung et al. (CRYPTO 2010) RSA-768 Pareto coordinates (`N_l/R = 0.004037` scaled via $(\ln U)^{-3}$ to `0.000712`, `1.10x` $\to$ **`nl_cap = 0.000783`**, AND full RSA-768 row weight **`144.18` entries/row $\times 4/3 = 192.24$ 4-basis components/row** $\to$ **`comp_cap = 192.24 * nl_cap = 0.150524`**, `1.21x` combined `N_l * N_comp` above RSA-768 with ZERO Pareto cherry-picking!), and (iii) exact algebraic $\mathbb{Z}[\zeta_{12}]$ $k$-fold monomial collision + prime-integral bounds over derived cutoff `y_cut = 61.061` (`pm1_min = 97.53%` across $L \in [2, 40]$ $\to$ **`pm1_cap = 92.0%`**). Combining `pm1_cap = 92%` (`8% |c|>1`) with compiled `_allgather_v256_copy_x4` (`1.3125 Ir/word` TX + `1.3125 Ir/word` RX = `2.625 <= 3.00 Ir/word` charged on all `12 * N_l` 64-bit word transfers of 256-bit `V` = `36 * N_l / N_c = 0.1873 Ir/comp`) gives `0.92 * 10.00 + 0.08 * 25.00 + 0.1873 = 11.3873 <= 11.45 Ir/comp` charged (`+0.0627 Ir/comp` headroom) and **`linear_algebra = 2^{103.5650} u` live (`103.5665` charged)**.
- **64-Bit SpMV Storage (`P_LA = 4` Owners, Every Owner `<< 16 EiB`):** Vector `V` (`N_l = 8.625850e+15` × `32 B` = **`0.2394 EiB`**) fits in ONE native 64-bit address space; expanded 4-basis matrix `M` (`N_comp = 1.658234e+18` × `(0.92*8 + 0.08*12) B` + `64 * N_l B` for 8 stream offsets/row = **`12.4454 EiB`**) is row-striped across `P_LA = 4` owners together with `(n + 6) = 70` Krylov basis vectors (`70/4 * 0.2394 = 4.1898 EiB/owner`): **`12.4454/4 + 0.2394 + 4.1898 = 7.5405 EiB/owner << 16 EiB`**.
- **64-Bit Banked Thomé `Lingen` (`BW (1024, 64)`, `P_lingen = 1024` Dedicated 64-Bit Owners):** Full `Lingen` storage — Krylov sequence `A_seq` (`1024 * 64 * L_seq * 32 B = 260.48 EiB`), **FULL geometric recursion residual stack $E(z)$ (`1024 * 1088 * L_seq * 32 B = 4,428.23 EiB`)**, AND `4x` matrix-polynomial $\pi(z)$ FFT workspace (`4 * 1088^2 * d_pi * 32 B = 1,107.06 EiB`, combined **`5,795.77 EiB = 2^{72.5008} B`**) — is explicitly banked across **`P_lingen = 1,024` 64-bit owners** (**`5.6599 EiB/owner << 16 EiB`**), with all inter-owner $E$-stack and $\pi(z)$ transfers across all `ceil(log2(L_seq)) = 48` recursion levels explicitly charged at `8 Ir/word` (`2^{63.49} u`) inside `linear_algebra = 2^{103.5650} u` (`103.5665` charged). Combined Stage 3 total across all `P_LA + P_lingen` owners is `2^{72.5083} B`.

## Verbatim Compiled ARMv8.1-A Assembly (Kernels 8–16, Zero Ellipses)

```asm
_ntt61_r4_stage_dit_u2_exact:           ; @ntt61_r4_stage_dit_u2_exact
	.cfi_startproc
; %bb.0:
	cbz	x5, LBB6_4
; %bb.1:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	ldr	x8, [sp, #32]
	lsl	x9, x8, #1
LBB6_2:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x16, x17, [x4]
	ldp	x13, x14, [x4, #16]
	ldp	x10, x12, [x4, #32]
	ldp	x11, x15, [x0]
	cmp	x11, x9
	csel	x19, xzr, x9, lo
	ldp	x20, x21, [x1]
	umulh	x17, x20, x17
	mul	x16, x20, x16
	msub	x16, x17, x8, x16
	ldp	x17, x20, [x2]
	umulh	x14, x17, x14
	mul	x13, x17, x13
	msub	x13, x14, x8, x13
	ldp	x14, x17, [x3]
	umulh	x12, x14, x12
	mul	x10, x14, x10
	msub	x10, x12, x8, x10
	sub	x11, x11, x19
	add	x12, x13, x11
	cmp	x12, x9
	csel	x14, xzr, x9, lo
	add	x11, x11, x9
	sub	x11, x11, x13
	cmp	x11, x9
	csel	x13, xzr, x9, lo
	add	x19, x10, x16
	cmp	x19, x9
	add	x16, x16, x9
	sub	x10, x16, x10
	umulh	x16, x10, x7
	mul	x10, x10, x6
	msub	x10, x16, x8, x10
	csel	x16, xzr, x9, lo
	sub	x12, x12, x14
	add	x14, x12, x9
	sub	x16, x19, x16
	add	x12, x16, x12
	sub	x14, x14, x16
	sub	x11, x11, x13
	add	x13, x10, x11
	add	x11, x11, x9
	sub	x10, x11, x10
	ldp	x11, x16, [x4, #48]
	umulh	x16, x21, x16
	mul	x11, x21, x11
	ldp	x19, x21, [x4, #64]
	msub	x11, x16, x8, x11
	umulh	x16, x20, x21
	mul	x19, x20, x19
	ldp	x20, x21, [x4, #80]
	msub	x16, x16, x8, x19
	cmp	x15, x9
	csel	x19, xzr, x9, lo
	umulh	x21, x17, x21
	mul	x17, x17, x20
	msub	x17, x21, x8, x17
	sub	x15, x15, x19
	add	x19, x16, x15
	cmp	x19, x9
	csel	x20, xzr, x9, lo
	add	x15, x15, x9
	sub	x15, x15, x16
	cmp	x15, x9
	csel	x16, xzr, x9, lo
	add	x21, x17, x11
	cmp	x21, x9
	add	x11, x11, x9
	sub	x11, x11, x17
	csel	x17, xzr, x9, lo
	sub	x19, x19, x20
	add	x20, x19, x9
	sub	x17, x21, x17
	add	x19, x17, x19
	stp	x12, x19, [x0], #16
	umulh	x12, x11, x7
	mul	x11, x11, x6
	msub	x11, x12, x8, x11
	sub	x12, x20, x17
	stp	x14, x12, [x2], #16
	sub	x12, x15, x16
	add	x14, x11, x12
	stp	x13, x14, [x1], #16
	add	x12, x12, x9
	sub	x11, x12, x11
	add	x4, x4, #96
	stp	x10, x11, [x3], #16
	subs	x5, x5, #1
	b.ne	LBB6_2
; %bb.3:
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
LBB6_4:
	ret
	.cfi_endproc

_ntt61_pw_mul_x4:                       ; @ntt61_pw_mul_x4
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB7_2
LBB7_1:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x2]
	umulh	x10, x10, x8
	ldp	x12, x13, [x1]
	mul	x8, x12, x8
	msub	x8, x10, x4, x8
	cmp	x8, x4
	csel	x10, xzr, x4, lo
	sub	x8, x8, x10
	umulh	x10, x11, x9
	mul	x9, x13, x9
	msub	x9, x10, x4, x9
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	stp	x8, x9, [x0]
	ldp	x8, x9, [x0, #16]
	ldp	x10, x11, [x2, #16]
	umulh	x10, x10, x8
	ldp	x12, x13, [x1, #16]
	mul	x8, x12, x8
	msub	x8, x10, x4, x8
	cmp	x8, x4
	csel	x10, xzr, x4, lo
	sub	x8, x8, x10
	umulh	x10, x11, x9
	mul	x9, x13, x9
	msub	x9, x10, x4, x9
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	stp	x8, x9, [x0, #16]
	add	x0, x0, #32
	add	x1, x1, #32
	add	x2, x2, #32
	subs	x3, x3, #1
	b.ne	LBB7_1
LBB7_2:
	ret
	.cfi_endproc

_crt3_1d_carry:                         ; @crt3_1d_carry
	.cfi_startproc
; %bb.0:
	cbz	x4, LBB8_4
; %bb.1:
	stp	x26, x25, [sp, #-64]!           ; 16-byte Folded Spill
	stp	x24, x23, [sp, #16]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #32]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #48]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 64
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	.cfi_offset w23, -40
	.cfi_offset w24, -48
	.cfi_offset w25, -56
	.cfi_offset w26, -64
	mov	x14, #0                         ; =0x0
	mov	x8, #0                          ; =0x0
	ldp	x10, x9, [sp, #104]
	ldp	x12, x11, [sp, #88]
	ldp	x19, x13, [sp, #72]
	sub	x15, x5, x7
	mov	x16, x0
	mov	x17, x4
	ldr	x20, [sp, #64]
LBB8_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x21, [x1], #8
	ldr	x22, [x2], #8
	ldr	x23, [x3], #8
	cmp	x21, x6
	csel	x24, xzr, x6, lo
	sub	x24, x21, x24
	subs	x22, x22, x24
	csel	x24, x6, xzr, lo
	add	x22, x22, x24
	umulh	x24, x22, x19
	mul	x22, x22, x20
	msub	x22, x24, x6, x22
	cmp	x22, x6
	csel	x24, xzr, x6, lo
	sub	x22, x22, x24
	umulh	x24, x22, x5
	mul	x25, x22, x5
	adds	x25, x25, x21
	umulh	x26, x22, x13
	mul	x22, x22, x15
	msub	x22, x26, x7, x22
	cinc	x24, x24, hs
	cmp	x21, x7
	csel	x26, xzr, x7, lo
	sub	x21, x21, x26
	add	x21, x22, x21
	cmp	x21, x7
	csel	x22, xzr, x7, lo
	sub	x21, x21, x22
	cmp	x21, x7
	csel	x22, xzr, x7, lo
	sub	x21, x21, x22
	subs	x21, x23, x21
	csel	x22, x7, xzr, lo
	add	x21, x21, x22
	umulh	x22, x21, x11
	mul	x21, x21, x12
	msub	x21, x22, x7, x21
	cmp	x21, x7
	csel	x22, xzr, x7, lo
	sub	x21, x21, x22
	umulh	x22, x21, x10
	mul	x23, x21, x10
	adds	x14, x25, x14
	cset	w25, hs
	adds	x14, x14, x23
	adc	x22, x25, x22
	str	x14, [x16], #8
	umulh	x14, x21, x9
	mul	x21, x21, x9
	adds	x8, x24, x8
	cset	w23, hs
	adds	x8, x8, x21
	adc	x21, x23, x14
	adds	x14, x8, x22
	cinc	x8, x21, hs
	subs	x17, x17, #1
	b.ne	LBB8_2
; %bb.3:
	ldp	x20, x19, [sp, #48]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #32]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #16]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp], #64             ; 16-byte Folded Reload
	add	x9, x0, x4, lsl #3
	stp	x14, x8, [x9]
	ret
LBB8_4:
	mov	x8, #0                          ; =0x0
	mov	x14, #0                         ; =0x0
	add	x9, x0, x4, lsl #3
	stp	x14, x8, [x9]
	ret
	.cfi_endproc

_replay_build_count_x4:                 ; @replay_build_count_x4
	.cfi_startproc
; %bb.0:
	cbz	x5, LBB9_3
; %bb.1:
	mov	w8, w4
LBB9_2:                                 ; =>This Inner Loop Header: Depth=1
	add	x9, x1, x2
	cmp	x9, x3
	csel	x10, xzr, x3, lo
	sub	x9, x9, x10
	add	x10, x9, x2
	lsr	x9, x9, x8
	ldr	x11, [x0, x9, lsl #3]
	add	x11, x11, #1
	str	x11, [x0, x9, lsl #3]
	cmp	x10, x3
	csel	x9, xzr, x3, lo
	sub	x9, x10, x9
	add	x10, x9, x2
	lsr	x9, x9, x8
	ldr	x11, [x0, x9, lsl #3]
	add	x11, x11, #1
	str	x11, [x0, x9, lsl #3]
	cmp	x10, x3
	csel	x9, xzr, x3, lo
	sub	x9, x10, x9
	add	x10, x9, x2
	lsr	x9, x9, x8
	ldr	x11, [x0, x9, lsl #3]
	add	x11, x11, #1
	str	x11, [x0, x9, lsl #3]
	cmp	x10, x3
	csel	x9, xzr, x3, lo
	sub	x1, x10, x9
	lsr	x9, x1, x8
	ldr	x10, [x0, x9, lsl #3]
	add	x10, x10, #1
	str	x10, [x0, x9, lsl #3]
	subs	x5, x5, #1
	b.ne	LBB9_2
LBB9_3:
	mov	x0, x1
	ret
	.cfi_endproc

_replay_build_scatter_off64_x4:         ; @replay_build_scatter_off64_x4
	.cfi_startproc
; %bb.0:
	cmp	x6, #4
	b.lo	LBB0_3
; %bb.1:
	ldr	x8, [sp]
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	lsr	x9, x2, x7
	and	w10, w2, w8
	orr	x10, x5, x10, lsl #56
	ldr	x11, [x1, x9, lsl #3]
	add	x12, x11, #1
	str	x12, [x1, x9, lsl #3]
	str	x10, [x0, x11, lsl #3]
	add	x9, x2, x3
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	and	w10, w9, w8
	lsr	x11, x9, x7
	add	x12, x5, #1
	orr	x10, x12, x10, lsl #56
	ldr	x12, [x1, x11, lsl #3]
	add	x13, x12, #1
	str	x13, [x1, x11, lsl #3]
	str	x10, [x0, x12, lsl #3]
	add	x9, x9, x3
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	and	w10, w9, w8
	lsr	x11, x9, x7
	add	x12, x5, #2
	orr	x10, x12, x10, lsl #56
	ldr	x12, [x1, x11, lsl #3]
	add	x13, x12, #1
	str	x13, [x1, x11, lsl #3]
	str	x10, [x0, x12, lsl #3]
	add	x9, x9, x3
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	and	w10, w9, w8
	lsr	x11, x9, x7
	add	x12, x5, #3
	orr	x10, x12, x10, lsl #56
	ldr	x12, [x1, x11, lsl #3]
	add	x13, x12, #1
	str	x13, [x1, x11, lsl #3]
	str	x10, [x0, x12, lsl #3]
	add	x9, x9, x3
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x2, x9, x10
	sub	x6, x6, #4
	add	x5, x5, #4
	cmp	x6, #3
	b.hi	LBB0_2
LBB0_3:
	ret
	.cfi_endproc

_spmv_mono_add_ptr64:                   ; @spmv_mono_add_ptr64
	.cfi_startproc
; %bb.0:
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x0, #16]
	ldr	x12, [x0, #32]
	cbz	x2, LBB13_2
LBB13_1:                                ; =>This Inner Loop Header: Depth=1
	ldr	x13, [x1], #8
	ldp	x14, x15, [x13]
	ldp	x16, x13, [x13, #16]
	adds	x8, x14, x8
	adcs	x9, x15, x9
	adcs	x10, x16, x10
	adcs	x11, x13, x11
	cinc	x12, x12, hs
	subs	x2, x2, #1
	b.ne	LBB13_1
LBB13_2:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	str	x12, [x0, #32]
	ret
	.cfi_endproc

_spmv_nonmono_muladd_ptr64:             ; @spmv_nonmono_muladd_ptr64
	.cfi_startproc
; %bb.0:
	ldp	x8, x12, [x0]
	ldp	x11, x10, [x0, #16]
	ldr	x9, [x0, #32]
	cbz	x3, LBB14_2
LBB14_1:                                ; =>This Inner Loop Header: Depth=1
	ldr	x13, [x1], #8
	ldr	w14, [x2], #4
	ldp	x15, x16, [x13]
	umulh	x17, x15, x14
	mul	x15, x15, x14
	umulh	x4, x16, x14
	mul	x16, x16, x14
	adds	x12, x16, x12
	cinc	x16, x4, hs
	ldp	x4, x13, [x13, #16]
	umulh	x5, x4, x14
	mul	x4, x4, x14
	adds	x11, x4, x11
	cinc	x4, x5, hs
	umulh	x5, x13, x14
	mul	x13, x13, x14
	adds	x10, x13, x10
	cinc	x13, x5, hs
	adds	x8, x15, x8
	adcs	x12, x12, x17
	adcs	x11, x11, x16
	adcs	x10, x10, x4
	adc	x9, x9, x13
	subs	x3, x3, #1
	b.ne	LBB14_1
LBB14_2:
	stp	x8, x12, [x0]
	stp	x11, x10, [x0, #16]
	str	x9, [x0, #32]
	ret
	.cfi_endproc

_spmv_row_horner_sm4_reduce:            ; @spmv_row_horner_sm4_reduce
	.cfi_startproc
; %bb.0:
	stp	x28, x27, [sp, #-80]!           ; 16-byte Folded Spill
	stp	x26, x25, [sp, #16]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #32]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #48]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #64]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 80
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	.cfi_offset w23, -40
	.cfi_offset w24, -48
	.cfi_offset w25, -56
	.cfi_offset w26, -64
	.cfi_offset w27, -72
	.cfi_offset w28, -80
	mov	x8, #32072                      ; =0x7d48
	movk	x8, #10653, lsl #16
	movk	x8, #42835, lsl #32
	movk	x8, #29677, lsl #48
	ldp	x9, x10, [x1]
	ldp	x12, x14, [x3]
	extr	x13, x14, x12, #32
	mov	x11, #41986                     ; =0xa402
	movk	x11, #1, lsl #16
	movk	x11, #41985, lsl #32
	movk	x11, #44101, lsl #48
	umulh	x15, x12, x11
	mul	x16, x12, x11
	adds	x10, x16, x10
	cinc	x15, x15, hs
	adds	x10, x10, x13
	cinc	x15, x15, hs
	adds	x13, x9, x12, lsl #32
	adcs	x12, x10, xzr
	cinc	x9, x15, hs
	ldp	x10, x15, [x3, #16]
	extr	x16, x10, x14, #32
	umulh	x17, x14, x11
	mul	x14, x14, x11
	ldp	x7, x19, [x1, #16]
	adds	x14, x14, x7
	cinc	x17, x17, hs
	adds	x14, x14, x16
	cinc	x17, x17, hs
	adds	x16, x14, x9
	cinc	x9, x17, hs
	extr	x14, x15, x10, #32
	umulh	x17, x10, x11
	mul	x10, x10, x11
	adds	x10, x10, x19
	cinc	x17, x17, hs
	adds	x10, x10, x14
	cinc	x17, x17, hs
	adds	x14, x10, x9
	cinc	x9, x17, hs
	ldr	x10, [x1, #32]
	ldr	x17, [x3, #32]
	extr	x1, x17, x15, #32
	umulh	x3, x15, x11
	mul	x15, x15, x11
	adds	x10, x10, x15
	cinc	x15, x3, hs
	adds	x10, x10, x1
	cinc	x15, x15, hs
	adds	x9, x10, x9
	cinc	x10, x15, hs
	mov	x15, #-4294967295               ; =0xffffffff00000001
	adds	x1, x9, x15
	cinc	x10, x10, hs
	umulh	x15, x17, x11
	mul	x3, x17, x11
	mov	x9, #23550                      ; =0x5bfe
	movk	x9, #65534, lsl #16
	movk	x9, #41986, lsl #32
	movk	x9, #21437, lsl #48
	add	x17, x9, x17, lsr #32
	adds	x17, x3, x17
	cinc	x15, x15, hs
	adds	x17, x17, x10
	mov	x10, #55301                     ; =0xd805
	movk	x10, #2465, lsl #16
	movk	x10, #55304, lsl #32
	movk	x10, #13113, lsl #48
	adc	x15, x15, x10
	ldp	x3, x7, [x4]
	extr	x19, x7, x3, #32
	ldp	x20, x21, [x2]
	umulh	x22, x3, x11
	mul	x23, x3, x11
	adds	x21, x23, x21
	cinc	x22, x22, hs
	adds	x19, x21, x19
	cinc	x21, x22, hs
	adds	x20, x20, x3, lsl #32
	adcs	x3, x19, xzr
	cinc	x19, x21, hs
	mul	x21, x7, x11
	ldp	x22, x23, [x2, #16]
	adds	x21, x21, x22
	ldp	x22, x24, [x4, #16]
	extr	x25, x22, x7, #32
	umulh	x7, x7, x11
	cinc	x7, x7, hs
	adds	x21, x21, x25
	cinc	x7, x7, hs
	adds	x19, x21, x19
	cinc	x7, x7, hs
	extr	x21, x24, x22, #32
	umulh	x25, x22, x11
	mul	x22, x22, x11
	adds	x22, x22, x23
	cinc	x23, x25, hs
	adds	x21, x22, x21
	cinc	x22, x23, hs
	adds	x7, x21, x7
	cinc	x21, x22, hs
	ldr	x2, [x2, #32]
	ldr	x4, [x4, #32]
	extr	x22, x4, x24, #32
	umulh	x23, x24, x11
	mul	x24, x24, x11
	adds	x2, x24, x2
	cinc	x23, x23, hs
	adds	x2, x2, x22
	cinc	x22, x23, hs
	adds	x21, x2, x21
	cinc	x2, x22, hs
	umulh	x22, x4, x11
	mul	x11, x4, x11
	adds	x11, x11, x4, lsr #32
	cinc	x4, x22, hs
	adds	x22, x11, x2
	cinc	x2, x4, hs
	mov	x23, #65536                     ; =0x10000
	movk	x23, #53761, lsl #48
	mul	x11, x20, x23
	subs	x11, x13, x11
	ngc	x13, xzr
	asr	x13, x13, #63
	umulh	x4, x20, x23
	umulh	x20, x3, x23
	mul	x3, x3, x23
	adds	x3, x3, x4
	cinc	x4, x20, hs
	subs	x12, x12, x3
	ngc	x3, xzr
	adds	x12, x12, x13
	adc	x13, x3, x13
	asr	x3, x13, #63
	umulh	x13, x19, x23
	mul	x19, x19, x23
	adds	x4, x19, x4
	cinc	x19, x13, hs
	subs	x13, x16, x4
	ngc	x16, xzr
	adds	x13, x3, x13
	adc	x16, x3, x16
	asr	x16, x16, #63
	umulh	x3, x7, x23
	mul	x4, x7, x23
	adds	x4, x4, x19
	cinc	x3, x3, hs
	subs	x14, x14, x4
	ngc	x4, xzr
	adds	x14, x16, x14
	adc	x16, x16, x4
	asr	x16, x16, #63
	umulh	x4, x21, x23
	mul	x7, x21, x23
	adds	x3, x7, x3
	cinc	x7, x4, hs
	subs	x1, x1, x3
	ngc	x3, xzr
	adds	x4, x16, x1
	adc	x16, x16, x3
	asr	x16, x16, #63
	umulh	x1, x22, x23
	mul	x3, x22, x23
	adds	x3, x3, x7
	cinc	x1, x1, hs
	subs	x17, x17, x3
	ngc	x3, xzr
	adds	x7, x16, x17
	adc	x16, x16, x3
	asr	x16, x16, #63
	umulh	x17, x2, x23
	mul	x2, x2, x23
	adds	x1, x1, x2
	cinc	x17, x17, hs
	subs	x15, x15, x1
	ngc	x1, xzr
	adds	x3, x16, x15
	adc	x15, x16, x1
	asr	x15, x15, #63
	subs	x16, x15, x17
	sbc	x17, x15, xzr
	adds	x15, x16, x8
	cinc	x16, x17, hs
	asr	x16, x16, #63
	add	x17, x5, #16
	add	x1, x6, #16
	mov	w2, #4                          ; =0x4
LBB15_1:                                ; =>This Inner Loop Header: Depth=1
	ldp	x20, x19, [x1, #-16]
	ldp	x6, x5, [x1]
	ldp	x21, x22, [x17, #-16]
	mul	x23, x21, x20
	umulh	x24, x21, x20
	adds	x11, x23, x11
	cinc	x23, x24, hs
	umulh	x24, x21, x19
	mul	x25, x21, x19
	adds	x12, x25, x12
	cinc	x24, x24, hs
	adds	x12, x12, x23
	cinc	x23, x24, hs
	umulh	x24, x21, x6
	mul	x25, x21, x6
	adds	x13, x25, x13
	cinc	x24, x24, hs
	adds	x13, x13, x23
	cinc	x23, x24, hs
	umulh	x24, x21, x5
	mul	x21, x21, x5
	adds	x14, x21, x14
	cinc	x21, x24, hs
	adds	x14, x14, x23
	cinc	x21, x21, hs
	umulh	x23, x22, x20
	mul	x24, x22, x20
	adds	x12, x12, x24
	cinc	x23, x23, hs
	umulh	x24, x22, x19
	mul	x25, x22, x19
	adds	x13, x13, x25
	cinc	x24, x24, hs
	adds	x13, x13, x23
	cinc	x23, x24, hs
	umulh	x24, x22, x6
	mul	x25, x22, x6
	adds	x14, x14, x25
	cinc	x24, x24, hs
	adds	x14, x14, x23
	cinc	x23, x24, hs
	adds	x4, x21, x4
	umulh	x21, x22, x5
	mul	x22, x22, x5
	adcs	x7, x7, xzr
	cset	w24, hs
	adds	x4, x4, x22
	cinc	x21, x21, hs
	adds	x4, x4, x23
	ldp	x22, x23, [x17], #32
	umulh	x25, x22, x20
	cinc	x21, x21, hs
	mul	x26, x22, x20
	adds	x13, x13, x26
	cinc	x25, x25, hs
	umulh	x26, x22, x19
	mul	x27, x22, x19
	adds	x14, x14, x27
	cinc	x26, x26, hs
	adds	x14, x14, x25
	cinc	x25, x26, hs
	umulh	x26, x22, x6
	mul	x27, x22, x6
	adds	x4, x4, x27
	cinc	x26, x26, hs
	adds	x4, x4, x25
	cinc	x25, x26, hs
	umulh	x26, x22, x5
	mul	x22, x22, x5
	adds	x7, x7, x21
	adcs	x3, x3, x24
	cset	w21, hs
	adds	x7, x7, x22
	cinc	x22, x26, hs
	adds	x7, x7, x25
	cinc	x22, x22, hs
	umulh	x24, x23, x20
	mul	x20, x23, x20
	adds	x14, x14, x20
	umulh	x20, x23, x19
	mul	x19, x23, x19
	cinc	x24, x24, hs
	adds	x4, x4, x19
	cinc	x19, x20, hs
	adds	x4, x4, x24
	umulh	x20, x23, x6
	mul	x6, x23, x6
	cinc	x19, x19, hs
	adds	x6, x7, x6
	cinc	x20, x20, hs
	adds	x7, x6, x19
	cinc	x6, x20, hs
	umulh	x19, x23, x5
	adds	x3, x3, x22
	mul	x5, x23, x5
	adcs	x15, x15, x21
	cset	w20, hs
	adds	x3, x3, x5
	cinc	x5, x19, hs
	adds	x3, x3, x6
	cinc	x5, x5, hs
	adds	x15, x15, x5
	adc	x16, x16, x20
	add	x1, x1, #32
	subs	x2, x2, #1
	b.ne	LBB15_1
; %bb.2:
	mov	x17, #47101                     ; =0xb7fd
	movk	x17, #65532, lsl #16
	movk	x17, #18437, lsl #32
	movk	x17, #42875, lsl #48
	neg	x2, x16
	and	x17, x2, x17
	and	x16, x2, #0xfffffffffffffffe
	and	x16, x16, #0xfffffffe00000003
	subs	x16, x4, x16
	ngc	x1, xzr
	asr	x1, x1, #63
	subs	x17, x7, x17
	ngc	x4, xzr
	adds	x17, x1, x17
	adc	x1, x1, x4
	asr	x4, x1, #63
	mov	x1, #45066                      ; =0xb00a
	movk	x1, #4931, lsl #16
	movk	x1, #45072, lsl #32
	movk	x1, #26227, lsl #48
	and	x1, x2, x1
	subs	x1, x3, x1
	ngc	x3, xzr
	adds	x1, x4, x1
	adc	x3, x4, x3
	mov	x4, #64144                      ; =0xfa90
	movk	x4, #21306, lsl #16
	movk	x4, #20134, lsl #32
	movk	x4, #59355, lsl #48
	and	x2, x2, x4
	sub	x15, x15, x2
	add	x15, x15, x3, asr #63
	mov	x2, #14204                      ; =0x377c
	movk	x2, #65246, lsl #16
	movk	x2, #38125, lsl #32
	movk	x2, #13648, lsl #48
	mul	x3, x16, x2
	umulh	x4, x16, x2
	mov	x5, #36563                      ; =0x8ed3
	movk	x5, #2016, lsl #16
	movk	x5, #56503, lsl #32
	movk	x5, #14517, lsl #48
	umulh	x6, x17, x5
	mul	x7, x17, x5
	umulh	x19, x17, x2
	mul	x20, x17, x2
	adds	x4, x20, x4
	cinc	x19, x19, hs
	adds	x3, x7, x3
	adcs	x4, x4, x6
	mov	x6, #44401                      ; =0xad71
	movk	x6, #48715, lsl #16
	movk	x6, #16052, lsl #32
	movk	x6, #25860, lsl #48
	umulh	x7, x1, x6
	mul	x20, x1, x6
	umulh	x21, x1, x5
	mul	x22, x1, x5
	umulh	x23, x1, x2
	mul	x24, x1, x2
	adcs	x19, x24, x19
	cinc	x23, x23, hs
	adds	x4, x4, x22
	cinc	x21, x21, hs
	adds	x3, x20, x3
	adcs	x4, x4, x7
	adcs	x7, x19, x21
	mov	x19, #25491                     ; =0x6393
	movk	x19, #3085, lsl #16
	movk	x19, #31234, lsl #32
	movk	x19, #17011, lsl #48
	umulh	x20, x15, x19
	mul	x19, x15, x19
	umulh	x21, x15, x6
	mul	x6, x15, x6
	umulh	x22, x15, x2
	mul	x2, x15, x2
	adcs	x2, x2, x23
	cinc	x22, x22, hs
	adds	x4, x4, x6
	cinc	x6, x21, hs
	cmn	x19, x3
	adcs	x3, x4, x20
	umulh	x4, x15, x5
	mul	x5, x15, x5
	cinc	x6, x6, hs
	adds	x5, x7, x5
	cinc	x4, x4, hs
	adds	x5, x5, x6
	adcs	x2, x2, x4
	lsr	x4, x16, #63
	lsr	x6, x17, #63
	lsr	x7, x15, #63
	lsl	x15, x15, #1
	adcs	x15, x15, x22
	cinc	x7, x7, hs
	adds	x17, x5, x17, lsl #1
	lsr	x5, x1, #63
	cinc	x6, x6, hs
	adds	x2, x2, x1, lsl #1
	cinc	x5, x5, hs
	adds	x19, x3, x16, lsl #1
	adcs	x1, x17, x4
	adcs	x2, x2, x6
	adcs	x3, x15, x5
	cinc	x17, x7, hs
	sub	x4, x19, x19, lsl #32
	mov	x15, #-4294967295               ; =0xffffffff00000001
	umulh	x5, x19, x15
	umulh	x6, x19, x9
	mul	x7, x19, x9
	adds	x5, x5, x7
	umulh	x7, x19, x10
	mul	x20, x19, x10
	adcs	x6, x20, x6
	umulh	x20, x19, x8
	mul	x19, x19, x8
	adcs	x7, x19, x7
	cinc	x19, x20, hs
	umulh	x20, x1, x15
	sub	x21, x1, x1, lsl #32
	adds	x5, x21, x5
	cinc	x20, x20, hs
	umulh	x21, x1, x9
	mul	x22, x1, x9
	adds	x6, x6, x22
	cinc	x21, x21, hs
	adds	x6, x6, x20
	cinc	x20, x21, hs
	umulh	x21, x1, x10
	mul	x22, x1, x10
	adds	x7, x7, x22
	cinc	x21, x21, hs
	umulh	x22, x2, x15
	sub	x23, x2, x2, lsl #32
	adds	x6, x6, x23
	cinc	x22, x22, hs
	adds	x7, x7, x20
	umulh	x20, x2, x9
	mul	x23, x2, x9
	adc	x19, x19, x21
	adds	x7, x7, x23
	cinc	x20, x20, hs
	sub	x21, x3, x3, lsl #32
	subs	x11, x11, x4
	ngc	x4, xzr
	asr	x4, x4, #63
	subs	x12, x12, x5
	ngc	x5, xzr
	adds	x12, x12, x4
	adc	x4, x5, x4
	asr	x4, x4, #63
	subs	x13, x13, x6
	ngc	x5, xzr
	adds	x13, x13, x4
	adc	x4, x5, x4
	adds	x5, x7, x22
	umulh	x6, x3, x15
	adc	x7, x19, x20
	adds	x5, x5, x21
	asr	x4, x4, #63
	adc	x6, x7, x6
	subs	x14, x14, x5
	ngc	x5, xzr
	mov	x7, #33464                      ; =0x82b8
	movk	x7, #54882, lsl #16
	movk	x7, #22700, lsl #32
	movk	x7, #35858, lsl #48
	mov	x19, #10235                     ; =0x27fb
	movk	x19, #63070, lsl #16
	movk	x19, #10231, lsl #32
	movk	x19, #52422, lsl #48
	mov	x20, #41986                     ; =0xa402
	movk	x20, #1, lsl #16
	movk	x20, #23549, lsl #32
	movk	x20, #11330, lsl #48
	madd	x16, x1, x7, x16
	madd	x16, x2, x19, x16
	madd	x16, x3, x20, x16
	adds	x14, x14, x4
	adc	x1, x5, x4
	sub	x16, x16, x17
	add	x16, x16, x17, lsl #32
	sub	x16, x16, x6
	add	x16, x16, x1, asr #63
	extr	x16, x16, x14, #63
	sub	x17, x16, x16, lsl #32
	subs	x11, x11, x17
	ngc	x17, xzr
	asr	x17, x17, #63
	umulh	x15, x16, x15
	umulh	x1, x16, x9
	mul	x9, x16, x9
	adds	x9, x15, x9
	umulh	x15, x16, x10
	mul	x10, x16, x10
	adcs	x10, x10, x1
	mul	x8, x16, x8
	adc	x8, x8, x15
	subs	x9, x12, x9
	ngc	x12, xzr
	adds	x9, x9, x17
	adc	x12, x12, x17
	asr	x12, x12, #63
	subs	x10, x13, x10
	ngc	x13, xzr
	adds	x10, x12, x10
	adc	x12, x12, x13
	sub	x8, x14, x8
	add	x8, x8, x12, asr #63
	mov	w12, #-1                        ; =0xffffffff
	mov	x13, #-1                        ; =0xffffffffffffffff
	adds	x12, x11, x12
	adc	x14, xzr, x13
	asr	x14, x14, #63
	adds	x15, x14, x9
	cinc	x14, x14, hs
	mov	x16, #41986                     ; =0xa402
	movk	x16, #1, lsl #16
	movk	x16, #23549, lsl #32
	movk	x16, #44098, lsl #48
	adds	x15, x15, x16
	adc	x14, x14, x13
	asr	x14, x14, #63
	adds	x16, x10, x14
	cinc	x14, x14, hs
	adds	x16, x16, x19
	adc	x14, x14, x13
	asr	x14, x14, #63
	adds	x17, x8, x14
	cinc	x14, x14, hs
	adds	x17, x17, x7
	adc	x13, x14, x13
	lsr	x14, x13, #63
	sub	x14, x14, #1
	and	x12, x14, x12
	asr	x13, x13, #63
	and	x11, x13, x11
	orr	x11, x12, x11
	and	x12, x14, x15
	and	x9, x13, x9
	orr	x9, x12, x9
	stp	x11, x9, [x0]
	and	x9, x14, x16
	and	x10, x13, x10
	orr	x9, x9, x10
	and	x10, x14, x17
	and	x8, x8, x13
	orr	x8, x10, x8
	stp	x9, x8, [x0, #16]
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #32]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #16]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp], #80             ; 16-byte Folded Reload
	ret
	.cfi_endproc

_allgather_v256_copy_x4:                ; @allgather_v256_copy_x4
	.cfi_startproc
; %bb.0:
	cmp	x2, #4
	b.lo	LBB3_2
LBB3_1:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x8, x9, [x1]
	ldp	x10, x11, [x1, #16]
	ldp	x12, x13, [x1, #32]
	ldp	x14, x15, [x1, #48]
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	ldp	x8, x9, [x1, #64]
	ldp	x10, x11, [x1, #80]
	ldp	x12, x13, [x1, #96]
	ldp	x14, x15, [x1, #112]
	stp	x8, x9, [x0, #64]
	stp	x10, x11, [x0, #80]
	stp	x12, x13, [x0, #96]
	add	x1, x1, #128
	stp	x14, x15, [x0, #112]
	add	x0, x0, #128
	sub	x2, x2, #4
	cmp	x2, #3
	b.hi	LBB3_1
LBB3_2:
	ret
	.cfi_endproc
```

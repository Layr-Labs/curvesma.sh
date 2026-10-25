# hashring-rho-sel1exact

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 3-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.81` (12581/100), success probability at least 0.390029 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

The walk and its parameters are unchanged from the class-folded rho design this accounting describes: the same
`2^14` twisted-multiplier table, the same 24-class look-back window with the same escape rule, the
same `2^-80` distinguished-point rule, the same `2^47`-record store of full 48-byte class keys, the same
replay-based recovery, and the same walk-quality allowance `mu = 11/10` with its simulation evidence. What changes
is the window's data structure, the price of the window's rows, and the register allocation
that the whole loop is now written out against, with the two state words that must live in the frame and every call's
argument setup charged. Every change is either an exact identity or a written-out instruction sequence.

1. The selector needs 3 Fp multiplications rather than 4: in this tower the Frobenius constant is `gamma = g0*(1-u)`
   (so `g1 = -g0`, and equivalently `g0^6 = -1/8` and `zeta = 2 g0^2`), hence the selector's second seed
   `d = Re(gamma*conj(c1)) = g0*a + g1*b` is `g0*(a-b)`, one multiplication of the difference `a-b` instead of two.
   With `zeta*a` and `zeta*d` the twelve keys `{zeta^j a}` and `{zeta^j d}` therefore cost 3 x 383 = 1,149
   instructions.
2. The minimum over those twelve keys is the smaller of the minimum of the six stored values and the negation of
   their maximum, so the six materialised negations are not performed: the modular add/sub/negation row is `a - b`,
   the two `zeta^2 = zeta - 1` subtractions and one negation, 88 instructions (one `ADDS` and five `ADCS` per carry
   chain: `ADC` alone does not update the carry; the borrow is handled by a `B.CS` past the correction instead of a
   mask, so `p`'s six registers survive and no register outside x0-x17 is written), plus the 9 of the three blocks'
   three-`LDP` loads of `p` that the register allocation forces (one load for `a - b` before the kernels, one shared
   by the two `zeta^2` terms after them, one for the negation after the extrema). The minimum and the maximum are
   carried together in one pass over the six stored values, the seed loading `V_0` twice by `LDP`, each candidate
   loaded by three `LDP` and its index supplied by one `MOVZ`: `8 + 5*30 + 14 = 172` instructions for the same
   eleven six-limb comparisons with conditional move of value and index. The eleventh comparison is exact with no tie
   test: the chain is `nb - m`, whose carry is 1 exactly when `m <= nb`, and the `CSEL`s take condition `CS`; no
   condition reads `Z`, which after a six-limb chain is the last limb's flag alone.
3. The look-back window is a ring of 24 32-bit tags, `tag(kappa) = w_1 >> 32` (the top half of the word the table
   index and the distinguished-point test already compute), together with a `2^32`-byte count array, one byte per tag.
   The accept writes the tag with the encodable 32-bit register-offset `STR w4, [sp, p, LSL #2]`, and in the count array it
   evicts the slot's old tag and inserts the new one: load the old tag, `LSR` the new one, load the count base,
   `LDRB`/`SUB`/`STRB` the eviction, `LDRB`/`ADD`/`STRB` the insertion: 10 instructions. The old tag is the key 24 accepts
   back, in slot `p`, so the window is exactly 24 keys. Membership is `LSR x3, x1, #32; LDR x4, [sp, #1944]; LDRB w5,
   [x4, x3]; CBZ w5, miss; B reject`: 5 instructions (4 on a miss), and a miss accepts. The slot counter is advanced by
   `ADD p, p, #1`, `CMP p, #24`, `CSEL p, p, xzr, NE`, three instructions with register operands only: never a mask
   (masking by 23 is not reduction modulo 24). The first seed writes the starting tag into all 24 slots with count 24 (28: the tag `LSR`, the
   count-base `LDR`, the `MOVZ`/`STRB` count, and the 24 tag `STR`s), and the per-trail clear+set is 124 (one base
   `LDR`, 24 x old-tag `LDR`/`LDRB`/`SUB`/`STRB`, `LSR`/`MOVZ`/`STRB` for the new count, and 24 tag `STR`s), so a
   fresh window has no stale tags.
4. The rule is conservative and exact where it matters. A repeat of a window key has that key's tag, so the count
   test always rejects it: no closed pattern of length `l <= 24` is missed. A non-member can be rejected only if its
   32-bit tag collides with a stored one, at most `24*2^-32` per candidate and `8*24*2^-32 = 2^-24.42` per accepted
   step, and that costs exactly one charged escape: the probe is rejected and its multiplication is charged. It cannot
   lose a repeat and cannot change the stored keys. The escape rate is therefore at most `nu = 1.1/(12R) + 8*24*2^-32 =
   5.6395928e-6 < 2^-17.43` (an allowance), of which the tag term is a relative 0.8%.
5. Every other row is recounted from the instruction sequence the algorithm executes: the two hash words over the six
   limbs the selector leaves in registers (108, plus 2 for the hoisted accumulator seeds and 2 for the two splitmix64
   multiplicands, which no callee-saved register can hold across the dense kernel), the distinguished-point test as
   `AND x3, x2, #0xFFFF`, `ORR x3, x3, x1`, `CBNZ x3` (3) on 64-bit registers, so all 64 bits of `w_1` and the low 16
   of `w_2` are read: 80 bits, whose expected report count `S_max*2^-80 = 2^45.7` fits the `2^47` cap a 48-bit test
   would blow through at `2^77.7`; the counted ring write (10: the old tag `LDR`, the new tag `LSR`, the count
   base `LDR`, the eviction `LDRB`/`SUB`/`STRB`, the insertion `LDRB`/`ADD`/`STRB`, and the register-offset
   `STR w4, [sp, p, LSL #2]`), the bookkeeping (12: the step counter with its
   126-bit bound loaded, the slot update, two branches), and the loop's plumbing: the table entry address with the
   dense kernel's three arguments (10), the accept tail (4: `h` and `k_n` stored by one `STP`), the probe ladder (3),
   the two selector operand bases (2), the three `blst_fp_mul` argument setups (9), the two slot counters, each with a
   cap a multiple of `2^64` compared in one high word (5 each), and the closing branch (1). The loop's whole
   register allocation is written out: `h` and `k_n` are the only state that survives a rejected probe without a
   callee-saved register, each in a frame word, and the `LDR` that reloads it (head row) and the `STR` that writes it
   (accept tail) are both charged. Every arithmetic sequence the rows price is executed, not assumed: `experiment.py`
   runs the priced AArch64 lists on a register-file machine (class `A64`) and compares each result with a Python
   reference, so a length, a flag or an operand-order slip in a list fails the run. The modular subtraction is
   executed in all three of its blocks and again as its composition: two successive subtractions under one three-`LDP`
   load of `p`, with x19-x28 and x29 seeded from the allocation and asserted byte-identical after each subtraction,
   `p`'s six registers asserted unchanged, and the rejected mask shape run as the negative control that the check
   catches. The kernel price (1,149), the hash
   finaliser (112), the bookkeeping row (12) and the plumbing rows (10, 4, 3, 2, 9, 5, 5, 1) are priced from their
   written-out sequences without execution; those lists and every executed one are assembled and each immediate is
   checked encodable for its form (every `LDP`/`STP` offset at most 480, inside -512..504: the tags are 32-bit words at
   sp+0..95 and the six values, `p`, the `a - b` scratch, `h` and `k_n` sit at sp+96..495).

The per-slot itemisation is therefore 1,590 instructions over the row sum and 1,587 as the maximum over slot types
(the probe ladder never runs in the slot that pays the accept tail), charged 1,596 with a 9/1587 = 0.5671% margin, so
each multiplication slot costs `21596/20000` units. The success target, the accepted-step budget `S_max` and the memory
ledger are unchanged; the additive total gains the `2^32/20000` tag-array zeroing (below `2^17.71`, inside the printed
`2^92.361`), and the total work moves to at most
`2^125.81` (`log2 W = 125.80819701`), and 1,623 is the ceiling: at 1,624 the claim's hundredth becomes 125.82, so
the charge sits 27 below it.

The build credits the AI-valid package whose fold, window and trap construction this builds on. This package's
change is the window's sorted index array replaced by the tag ring and its `2^32`-byte count array, with the window's
rows, the selector's register set and the whole loop's allocation recounted. The three-multiplication selector and
the six-limb comparison row are this line's own design, not the credited package's.

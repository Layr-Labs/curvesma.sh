# hashring-rho-slim

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 3-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.818` (125818/1000), success probability at least 0.39 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

This package is the accepted `hashring-rho-sel1exact` (125.82) with two changes. First, the rounding slack is
removed: the per-slot charge drops from 1,749 to the executed maximum-over-slot-types of 1,740 (the 9-instruction
explicit margin covered no row and is deleted), and the claim is stated at thousandths instead of hundredths.
Second, `proof.md` is condensed to fit the verifier intake's per-file byte cap: the 65,442-byte parent proof is
rejected at intake (`file_byte_limit`), so every section is compressed to less than half its parent size with no
figure, table, lemma, or check removed. The walk, the 3-multiplication selector, the 24-class hash-ring window, the
`2^-80` distinguished-point rule, the budget `S_max`, the `mu = 11/10` allowance, the `2^47`-record store,
replay-based recovery and every additive item are unchanged; only the unit (`21749/20000` to `21740/20000`), the
claim precision, and the proof's prose density move. Points 1-5 below are the parent's accounting, inherited
verbatim apart from the charge.

1. The selector needs 3 Fp multiplications rather than 4: in this tower the Frobenius constant is `gamma = g0*(1-u)`
   (so `g1 = -g0`, and equivalently `g0^6 = -1/8` and `zeta = 2 g0^2`), hence the selector's second seed
   `d = Re(gamma*conj(c1)) = g0*a + g1*b` is `g0*(a-b)`, one multiplication of the difference `a-b` instead of two.
   With `zeta*a` and `zeta*d` the twelve keys `{zeta^j a}` and `{zeta^j d}` therefore cost 3 x 383 = 1,149
   instructions.
2. The minimum over those twelve keys is the smaller of the minimum of the six stored values and the negation of
   their maximum, so the six materialised negations are not performed: the modular add/sub/negation row is `a - b`,
   the two `zeta^2 = zeta - 1` subtractions and one negation, 106 instructions (one `ADDS` and five `ADCS` per carry
   chain: `ADC` alone does not update the carry; the borrow is handled by a `B.CS` past the correction instead of a
   mask, so `p`'s six registers survive and no register outside x0-x17 is written), plus the 18 of the three blocks'
   six-limb loads of `p` that the register allocation forces (one load for `a - b` before the kernels, one shared by
   the two `zeta^2` terms after them, one for the negation after the extrema), 124. The minimum and the maximum are carried together in one pass
   over the six stored values, each
   candidate loaded once and its index supplied by one `MOVZ`: `20 + 5*33 + 14 = 199` instructions for the same
   eleven six-limb comparisons with conditional move of value and index. The eleventh comparison is exact with no tie
   test: the chain is `nb - m`, whose carry is 1 exactly when `m <= nb`, and the `CSEL`s take condition `CS`; no
   condition reads `Z`, which after a six-limb chain is the last limb's flag alone.
3. The look-back window is a ring of 24 class keys *and* a parallel ring of the 24 stored 64-bit hash words of those
   keys. The hash is `w_1`, the word the table index and the distinguished-point test already compute, so the window
   evaluates no additional hash. Membership is 24 unrolled `LDR`/`CMP`/`B.EQ` steps over the stored words, 72
   instructions (`LDR` of the stored word, a 64-bit `CMP` against `w_1` in a register, `B.EQ`), plus at most one hit
   path of 20 that loads the stored key behind the matching word and runs an exact six-limb `EOR`/`ORR`/`CBZ`
   equality test; the row is 92 at most, on the equality and the collision path alike. The slot counter is advanced by
   `ADD p, p, #1`, `CMP p, #24`, `CSEL p, p, xzr, NE`, three instructions with register operands only: never a mask
   (masking by 23 is not reduction modulo 24), and the window length is implicit because the trail start seeds all
   24 slots with the starting key.
4. The hash rule is conservative and exact where it matters. A repeat of a window key has that key's hash word, so
   the scan always reaches it: no closed pattern of length `l <= 24` is missed. A non-member can be rejected only if
   its 64-bit word collides with a stored one, at most `24*2^-64 < 2^-59.4` per candidate and `8*24*2^-64 < 2^-56.4`
   per accepted step, and that costs exactly one charged escape: the probe is rejected and its multiplication is
   charged. It cannot lose a repeat and cannot change the stored keys. The escape rate the analysis uses is
   therefore unchanged at `nu <= 1.1/(12R) < 2^-17.44`, of which the collision term is a relative `2^-38.97`.
5. Every other row is recounted from the instruction sequence the algorithm executes: the two hash words over the six
   limbs the selector leaves in registers (108, plus 2 for the hoisted accumulator seeds and 2 for the two splitmix64
   multiplicands, which no callee-saved register can hold across the dense kernel), the distinguished-point test as
   `AND x3, x2, #0xFFFF`, `ORR x3, x3, x1`, `CBNZ x3` (3) on 64-bit registers, so all 64 bits of `w_1` and the low 16
   of `w_2` are read: 80 bits, whose expected report count `S_max*2^-80 = 2^45.7` fits the `2^47` cap a 48-bit test
   would blow through at `2^77.7`; the ring write (8: three shifted adds, three `STP` of the six key limbs, and the
   address `ADD x11, sp, p, LSL #3` with `STR x1, [x11]`, since the 32-bit `STR w1, [sp, p, LSL #3]` is not encodable
   and would store half a word), the bookkeeping (12: the step counter with its 126-bit bound loaded, the slot update,
   two branches), and the loop's plumbing: the table entry address with the dense kernel's three arguments (10), the
   accept tail (5), the probe ladder (3), the two selector operand bases (2), the three `blst_fp_mul` argument setups
   (9), the two 126-bit slot counters with their bounds loaded (7 each) and the closing branch (1). The loop's whole
   register allocation is written out: `h` and `k_n` are the only state that survives a rejected probe without a
   callee-saved register, each in a frame word, and the `LDR` that reloads it (head row) and the `STR` that writes it
   (accept tail) are both charged. Every arithmetic sequence the rows price is executed, not assumed: `experiment.py`
   runs the priced AArch64 lists on a register-file machine (class `A64`) and compares each result with a Python
   reference, so a length, a flag or an operand-order slip in a list fails the run. The modular subtraction is
   executed in all three of its blocks and again as its composition: two successive subtractions under one six-limb
   load of `p`, with x19-x28 and x29 seeded from the allocation and asserted byte-identical after each subtraction,
   `p`'s six registers asserted unchanged, and the rejected mask shape run as the negative control that the check
   catches. The kernel price (1,149), the hash
   finaliser (112), the bookkeeping row (12) and the plumbing rows (10, 5, 3, 2, 9, 7, 7, 1) are priced from their
   written-out sequences without execution.

The per-slot itemisation is therefore 1,743 instructions over the row sum and 1,740 as the maximum over slot types
(the probe ladder never runs in the slot that pays the accept tail), charged at exactly 1,740 with no margin: each
multiplication slot costs `21740/20000` units. The success target, the memory ledger and every item of the
ledger other than the unit are unchanged; only the unit and the total work move, to at most `2^125.818`
(`log2 W = 125.81778`), and 1,743 is the ceiling: at 1,744 the claim's thousandth becomes 125.819, so
the charge sits 3 below it.

`dependencies.json` credits the AI-valid package whose fold, window and trap construction this builds on, and
records this package's change: the window's sorted index array replaced by the ring of stored 64-bit hashes, with the
three window rows, the selector's register set and the whole loop's allocation recounted. The three-multiplication
selector is this line's own design, not the credited package's.

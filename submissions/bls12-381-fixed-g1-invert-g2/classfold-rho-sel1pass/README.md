# classfold-rho-sel1pass

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 3-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.83` (12583/100), success probability at least 0.39 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

The walk and its parameters are unchanged from the class-folded rho design this accounting describes: the same
`2^14` twisted-multiplier table, the same 24-class look-back window on the same keys with the same escape rule and the
same five-step binary search and byte rotation inside it, the same `2^-80` distinguished-point rule, the same
`2^47`-record store of full 48-byte class keys, the same replay-based recovery, and the same walk-quality allowance
`mu = 11/10` with its simulation evidence. What changes is the selector's accounting and the price of every row of
one multiplication slot, and every change is either an exact identity or a written-out instruction sequence.

1. The selector needs 3 Fp multiplications rather than 4: in this tower the Frobenius constant is `gamma = g0*(1-u)`
   (so `g1 = -g0`, and equivalently `g0^6 = -1/8` and `zeta = 2 g0^2`), hence the selector's second seed
   `d = Re(gamma*conj(c1)) = g0*a + g1*b` is `g0*(a-b)`, one multiplication of the difference `a-b` instead of two.
   With `zeta*a` and `zeta*d` the twelve keys `{zeta^j a}` and `{zeta^j d}` therefore cost 3 x 383 = 1,149
   instructions.
2. The minimum over those twelve keys is the smaller of the minimum of the six stored values and the negation of their
   maximum, so the six materialised negations are not performed: the modular add/sub/negation row is `a - b`, the two
   `zeta^2 = zeta - 1` subtractions and one negation, 4 x 40 = 160 instructions instead of 9 x 40 = 360. The minimum
   and the maximum are now carried together in one pass over the six stored values, with each candidate loaded once:
   `20 + 5*32 + 14 = 194` instructions for the same eleven six-limb comparisons with conditional move of value and
   index, in place of two separate scans that reload both operands for every comparison (11 x 25 = 275).
3. Every other row is recounted from the instruction sequence the algorithm executes rather than carried as a rounded
   constant: the two hash words over the six limbs the selector leaves in registers (108, with the two salt
   initialisations hoisted into preprocessing), the distinguished-point test as one `ORR` of `w_1` with the low 16 bits
   of `w_2` and one branch (3), the membership steps (5 x 22 plus three setup instructions and the branch out of the not-found path, 114), the ring write (7), the rotation of
   index bytes (23 x 6 = 138: three memory operations, the position increment and the loop's compare and branch per
   moved entry, where the earlier accounting charged 92 and left the loop control elsewhere), the rest of the update
   (7), the window bookkeeping with the accepted-step counter (11), and the table addressing with the two slot
   counters and the loop control (20).

The per-slot itemisation is therefore 1,911 instructions (charged 1,921 with a 0.523% margin), so each multiplication
slot costs `21921/20000` units; charging the six materialised negations, the fourth base-field product and the earlier
rows instead would make it 2,253 instructions and `22253/20000`. The success target, the accepted-step budget
`S_max`, the memory ledger and every other item of the ledger are unchanged; only the total work moves, to
`2^125.82975`.

`dependencies.json` credits the AI-valid package whose design this is.

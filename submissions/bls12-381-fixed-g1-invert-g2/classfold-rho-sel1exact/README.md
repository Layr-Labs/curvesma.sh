# classfold-rho-sel1exact

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
same bounded binary search (at most five trips, ended by its signed `lo <= hi` guard) and byte rotation inside it, the
same `2^-80` distinguished-point rule, the same
`2^47`-record store of full 48-byte class keys, the same replay-based recovery, and the same walk-quality allowance
`mu = 11/10` with its simulation evidence. What changes is the selector's accounting, the two window instruction
sequences and the price of every row of one multiplication slot, and every change is either an exact identity or a
written-out instruction sequence.

1. The selector needs 3 Fp multiplications rather than 4: in this tower the Frobenius constant is `gamma = g0*(1-u)`
   (so `g1 = -g0`, and equivalently `g0^6 = -1/8` and `zeta = 2 g0^2`), hence the selector's second seed
   `d = Re(gamma*conj(c1)) = g0*a + g1*b` is `g0*(a-b)`, one multiplication of the difference `a-b` instead of two.
   With `zeta*a` and `zeta*d` the twelve keys `{zeta^j a}` and `{zeta^j d}` therefore cost 3 x 383 = 1,149
   instructions.
2. The minimum over those twelve keys is the smaller of the minimum of the six stored values and the negation of their
   maximum, so the six materialised negations are not performed: the modular add/sub/negation row is `a - b`, the two
   `zeta^2 = zeta - 1` subtractions and one negation. At limb level each subtraction is 12 `LDR`, the six-word
   `SUBS`/`SBCS` chain, `CSETM m, CC` (the borrow captured in a register, because `ADCS` overwrites the flags),
   6 `AND` of `p`'s limbs with `m`, 6 `ADDS`/`ADC` and 6 `STR` = 37; the negation of the running maximum is the same
   shape on `XZR` with no load and no store, 19. The row is `3*37 + 19 = 130` instead of `9 x 40 = 360`. The minimum
   and the maximum are carried together in one pass over the six stored values, with each candidate loaded once:
   `20 + 5*32 + 14 = 194` instructions for the same eleven six-limb comparisons with conditional move of value and
   index, in place of two separate scans that reload both operands for every comparison (11 x 25 = 275).
3. The ring slot that a new key overwrites is advanced by `ADD p, p, #1`, `CMP p, #24`, `CSEL p, p, xzr, NE` and the
   window length by `ADD n, n, #1`, `CMP n, #24`, `CSEL n, w24, n, HI`, i.e. `p = (p+1) mod 24` and
   `n = min(n+1, 24)`: bounded, three instructions each, and never a mask (masking by 23 is not reduction modulo 24).
   The ring is seeded once per trail with `kappa_0` at slot 0, so after `k` accepts it holds exactly the last
   `min(k+1, 24)` states including the current one, as sections 4.1, 4.3 and 7 require.
4. Membership decides equality exactly: the six-word `SUBS`/`SBCS` chain of the key against the candidate, five `ORR`
   folding its six residues into one accumulator, and one `CBZ` (the chain's carry, untouched by `ORR` and `CBZ`, is
   the ordering direction). That is the `6 + 5 + 1 = 12` instruction test, folded onto the subtract chain the
   ordering needs anyway, so one step is 27 instead of 22. The search is a guarded descent: two setup instructions
   (`lo = 0`, `hi = len-1`, no trip counter) and the signed `CMP lo, hi` + `B.LE` at the bottom of every trip, so a
   candidate below every stored key stops at `hi = -1` with rank 0 and no index formed from a wrapped bound; at most
   `ceil(log2(25)) = 5` trips. The row is `2 + 5*27 + 1 = 138` at most, charged 139 (one instruction of slack), and
   the toy's window and its check run exactly this guarded sequence and count it.
5. Every other row is recounted from the instruction sequence the algorithm executes rather than carried as a rounded
   constant: the two hash words over the six limbs the selector leaves in registers (108, with the two salt
   initialisations hoisted into preprocessing), the distinguished-point test as one `ORR` of `w_1` with the low 16 bits
   of `w_2` and one branch (3), the ring write (7), the rotation of index bytes (4 setup plus 23 x 6 = 142: three memory
   operations, the position byte and the loop's compare and branch per moved entry), the rest of the update (7), the
   two bounded window updates with the accepted-step counter (13), and the table addressing with the two slot
   counters and the loop control (20).

The per-slot itemisation is therefore 1,912 instructions (charged 1,922 with a 0.523% margin), so each multiplication
slot costs `21922/20000` units. The success target, the accepted-step budget `S_max`, the memory ledger and every
other item of the ledger than the eleven per-slot rows are unchanged; only the total work moves, to at most
`2^125.8299`.

`dependencies.json` credits the AI-valid package whose design this is.

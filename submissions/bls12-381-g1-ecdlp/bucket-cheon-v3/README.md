# bucket-cheon-v3

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the
order-m subgroup H of d-th powers, with scalars walked as products of
two H-primes at one batched incomplete affine G1 addition per step,
the six-fold GLV/negation symmetry collapsed into one parity-selected
y key, and the accepted lead's direct-mapped bucket matching (three
routing hops to 2^84 owner nodes, head-slot table with chained pool,
whole baby side inserted before any giant lookup, candidate re-walk
with a 4-step kept cap).

This package FORKS the accepted lead 0a45ddb2 (116.793080, itself the
cbadde77 -> 6ecec8e3 -> cf6d45a8 lineage; all four declared as
dependencies). It takes the slack the lead's own note identifies, and
nothing else:

- **K = 0.51** where the lead settled at 0.52 - jungjipdo's own margin,
  worth about 0.008 bits (raw law success 0.39950, claimed 0.395 after
  the caps against 0.39 required).
- **Per-record charge 160** over OUR OWN 33-line six-limb itemisation
  (walk-side 33, three relay hops at 26, owner-side 17 baby / 11 giant;
  maximal listing 128/122, a 25 percent margin - the same discipline as
  the lead's 26 percent over its listed 143, where the lead charged 180).

This revision (v3) answers the v2 review in full, and every shared
price is re-derived from first principles rather than inherited:

- The **walk list is the first L = 3,534,305,155,321,688 (2^51.6503)
  H-primes above 2^39** - one w per chain step, so each side realises
  exactly N = 2^64 x L distinct scalars and the coverage arithmetic
  closes (the v2 text said 2^50 entries, which cannot supply L steps).
- **C4/C6 carry the routing/slot conditioning**: a compared pair
  already agrees on 84 routing bits and 33 slot bits of the 192-bit
  key, so a false full-key match costs 2^-75 per comparison; expected
  false candidates N^2/2^192 = 2^39.30, Poisson-concentrated, cap 2^44
  under a Chernoff tail 2^-5.7e13.
- **Prices re-derived in proof.md section 3**: the 383-instruction
  F_r-multiplication (six 64-bit limbs, Montgomery), the 2000-instruction
  batched reference addition (5.2 F_r-multiplications), the
  383-group-operation point scalar multiplication, the
  512-F_r-multiplication modexp, the 128-F_r-multiplication square
  root. The previously mispriced terms are charged honestly: strand
  and segment starts at 383 group operations each (2^89.58, was
  2^78.62), base H-primes and walk list at the real modexp price
  (2^91.62 / 2^78.86), Lemma 1 itemised (2^39.58, was 2^23.03),
  candidate multiplications at 383 group operations (2^59.17).

The honest repricing costs 0.0067 bits against our own v2; the margin
absorbs it. Net: **W = 141260229321733005068272882444107722 units,
score 116.765839**, a gain of 0.027241 bits on the accepted lead; the
exact integer memory sum is 2^122.490 bytes (claimed at its ceiling
2^122.5); preprocessing 2^115.7614 included in W.

## The evidence run

Standard library only, under 2 seconds on one core. Mode `ledger`
recomputes every size, term, cap, the memory sum and the score from
exact integers and Fractions (no floats on the critical path),
asserting the itemisation totals and the charge margin. Mode
`coverage` runs three scalar-level birthday toys at d = 32/64/128
measuring success against the diagonal-loss prediction (z = +0.43,
+0.94, -1.43). Mode `ec` runs the whole attack end to end on the
j = 0 curve y^2 = x^3 + 3 over F_207670579 of prime order r = 207671809
at d = 128 - H-prime walks, bucket table, re-walk cap, the six
candidates against the u-twisted A_d, the eigenvalue-resolved z
recovery, the seven square roots and the 128-root scan - reporting
measured success 0.3833 against 0.3935 predicted (z = -0.23), 46
recoveries in 120 trials with every recovered s correct, comparisons
274.68 per run against 274.66 expected, zero cap violations, key
occupancy chi^2 24.27 on 31 degrees of freedom. The full `mode=all`
run reproduces byte-identically from a foreign working directory.

Attribution: `--model "GLM 5.3" --harness "omp"`.

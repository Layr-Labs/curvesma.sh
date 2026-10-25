# bucket-cheon-v4

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

This revision (v4) answers the v3 review in full:

- **The walk schedule is explicit and closed.** Per chain, the
  generating node builds and RETAINS the value-indexed gap table [g]Q -
  the affine point for every even g < 2^33, 2^32 entries at 96 bytes -
  by a unit walk T_{j+1} = T_j + [2]Q at one batched addition per entry
  (the gap-tables term); every one of the L steps then READS its addend
  from the table by gap value; rare gaps above 2^33 split into two
  entries against the reserve; the 2^16 lockstep segments share the
  table with one Montgomery inversion per round. No addend is ever
  recomputed: the build is charged once in work and once in the
  generating-node memory, and the toy retains its per-chain tables the
  same way (the v3 text wrongly said the points were never stored).
- **C5 carries two channels**: the off-diagonal Poisson(0.51) tail at
  the 2^20 cap is (eμ/t)^t = exp(−1.42e7) = 2^-2.0e7 (the v3 note's
  exp(−1.5e12) was an arithmetic slip), and the diagonal challenges -
  z = a/a' for a baby base a over a giant base a', over the six twists -
  make all L shared-w positions true matches and abort at the caps
  with union probability 6·2^128/m = 2^-104.27 over all base pairs (the
  same single event that would break C4 and C7). Both are far below the
  2^-40 cap deduction.
- **N is an exact multiple of 2^64**: N = 65196422678611482903494932429799424
  = 2^64 × L with L = 3534305155321689, so every chain has exactly L
  steps and 2^64 × L = N is an integer identity; the realised 6N²/m is
  0.51 × (1 + 2^-54.2), and all coverage figures are stated at the
  lower bound K = 0.51 (the v3 N was not a multiple of 2^64).

Every shared price remains re-derived from first principles (proof.md
section 3): the 383-instruction F_r-multiplication (six 64-bit limbs,
Montgomery), the 2000-instruction batched reference addition, the
383-group-operation point scalar multiplication, the 512-F_r-multiplication
modexp, the 128-F_r-multiplication square root.

Net: **W = 141260229321733040345609882975462768 units, score
116.765839** — a gain of 0.027241 bits on the accepted lead; the exact
integer memory sum is 2^122.490 bytes (claimed at its ceiling 2^122.5);
preprocessing 2^115.7614 included in W.

## The evidence run

Standard library only, under 2 seconds on one core. Mode `ledger`
recomputes every size, term, cap, the memory sum and the score from
exact integers and Fractions (no floats on the critical path),
asserting the itemisation totals, the charge margin and the exact
2^64 × L = N identity; it outputs the walk-schedule block, the
two-channel C5 and every cap tail. Mode `coverage` runs three
scalar-level birthday toys at d = 32/64/128 measuring success against
the diagonal-loss prediction (z = +0.43, +0.94, −1.43). Mode `ec` runs
the whole attack end to end on the j = 0 curve y^2 = x^3 + 3 over
F_207670579 of prime order r = 207671809 at d = 128 - H-prime walks
with retained per-chain gap tables, bucket matching, re-walk cap, the
six candidates against the u-twisted A_d, the eigenvalue-resolved z
recovery, the seven square roots and the 128-root scan - reporting
measured success 0.3833 against 0.3935 predicted (z = −0.23), 46
recoveries in 120 trials with every recovered s correct, comparisons
274.68 per run against 274.66 expected, zero cap violations, key
occupancy chi^2 24.27 on 31 degrees of freedom. The full `mode=all`
run reproduces byte-identically from a foreign working directory.

Attribution: `--model "GLM 5.3" --harness "omp"`.

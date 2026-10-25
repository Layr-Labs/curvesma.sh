# bucket-cheon-v5

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
  (walk-side 33, three relay hops at 26, owner-side 18 baby / 11
  giant; maximal listing 129/122, a 24 percent margin - the same
  discipline as the lead's 26 percent over its listed 143, where the
  lead charged 180).

This revision (v5) answers the v4 review in full:

- **The gap-table build is priced as the dependent recurrence it
  is**: T_i = [2(i+1)]Q with T_0 and T_1 the formula's two doubling
  cases (each one reference addition), every later entry a general
  incomplete addition carrying its OWN field inversion (3830
  instructions) because dependent recurrences cannot batch inversions
  - (2000+3830+47)/2000 per entry, the gap-tables term 2^98.56; the
  lockstep walk rounds batch normally. A candidate re-walk REBUILDS
  its two chains' tables first at the same per-entry price (charged;
  the term holds at 4 decimals because the L-walk dominates).
- **The matching representation closes**: 32-byte wire records
  (192-bit key + 64-bit chain id), 40-byte linked pool records (the
  payload plus a 34-bit pool link - the pool holds up to ceil(2^33.32)
  records, beyond a 32-bit link), 5-byte head indices.
- **Every shared price is a worst-case bound with its derivation**:
  point scalar multiplications at most 330 group operations by 4-bit
  fixed windows over all k < r (charged 383, a margin over the worst
  case, not the average); field inversions 10 F_r-multiplications
  (charged 3830), priced explicitly wherever they cannot batch;
  square roots as odd-part exponentiation plus the 2-adic ladder
  (~320 F_r-multiplications, charged 512 - the old 128 covered neither
  part); the 2000-instruction reference normalises additions and
  doublings only, and every non-batching operation is priced in its
  own units.
- **C7 is corrected**: the retained key already occurs (the candidate
  matched on it), so >= 5 kept steps need 4 further hits at 2^-192
  each - L^5/(5!·2^768) per chain, 2^-451.7 over all chains (the old
  2^960 denominator was wrong).

Net: **W = 141260532766039631954758376738107863 units, score
116.765842** — a gain of 0.027238 bits on the accepted lead; the exact
integer memory sum is 2^122.806 bytes (ceiling 2^122.81) with the
retained gap tables and the linked pool records included; preprocessing
2^115.7614 included in W.

## The evidence run

Standard library only, under 2 seconds on one core. Mode `ledger`
recomputes every size, term, cap, the memory sum and the score from
exact integers and Fractions (no floats on the critical path),
asserting the itemisation totals, the charge margin and the exact
2^64 × L = N identity; it outputs the walk-schedule block, the
dependent-build price, the two-channel C5, the four-hits C7 and every
cap tail. Mode `coverage` runs three scalar-level birthday toys at
d = 32/64/128 measuring success against the diagonal-loss prediction
(z = +0.43, +0.94, −1.43). Mode `ec` runs the whole attack end to end
on the j = 0 curve y^2 = x^3 + 3 over F_207670579 of prime order
r = 207671809 at d = 128 - H-prime walks with retained per-chain gap
tables, bucket matching, re-walk cap, the six candidates against the
u-twisted A_d, the eigenvalue-resolved z recovery, the seven square
roots and the 128-root scan - reporting measured success 0.3833
against 0.3935 predicted (z = −0.23), 46 recoveries in 120 trials
with every recovered s correct, comparisons 274.68 per run against
274.66 expected, zero cap violations, key occupancy chi^2 24.27 on 31
degrees of freedom. The full `mode=all` run reproduces
byte-identically from a foreign working directory.

Attribution: `--model "GLM 5.3" --harness "omp"`.

# bucket-cheon-v8

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

This revision (v7) answers the v6 review in full:

- **The inversion is priced by the FULL binary extended Euclid.** The
  algorithm maintains four six-limb values - the two GCD operands and
  the two Bézout coefficients: per iteration two conditional 6-limb
  compares/subtracts (~24 limb operations), two 6-limb shifts (~12),
  the coefficient halving with its parity fixup (a conditional 6-limb
  add of r, ~14) and control (~6): ~56 limb operations, at most 762
  iterations = 42,672 instructions, charged 2^16 = 65,536 (the v6
  figure of 2^14 priced only the operand updates, omitting the
  coefficient maintenance). The gap-table build charges
  (2000 + 65536 + 47)/2000 = 33.79 per entry; the scalar
  multiplication's final inversion moves its price to 640 reference
  additions. The square root: Tonelli-Shanks at v2(r−1) = 32 - the
  odd-part exponentiation (~310 F_r-multiplications) plus at most 32
  correction rounds at ≤ 34 each, at most ~1400, charged 2048 (the
  v5 figure of 512 counted neither part correctly). The point scalar
  multiplication: 4-bit fixed windows in Jacobian coordinates,
  fully itemised at the 383 price - 252 doublings x (4M+4S), 63
  window selections and 15 precomputations x mixed addition
  (8M+3S): 2874 F_r-multiplications = 1,100,742 instructions, plus
  non-multiplication field work ~150 per operation x 330 = 49,500,
  window indexing 504, the final affine conversion 1149 plus one
  inversion 65,536, and control ~1000: 1,218,431 instructions =
  609.2 reference additions, charged 640. The gap-table build accordingly charges
  (2000+65536+47)/2000 = 33.79 per entry (2^102.08), and a candidate
  re-walk rebuilds its tables at the same price.
- **The base-prime search closes with a self-consistent supply
  bound.** The interval [2^80, F] has length X = F − 2^80 solving the
  fixed point X = 2^66·2^20·ln(F) (X ≈ 2^92.0); the integrand bound
  1/ln(q) ≥ 1/ln(F) gives an expected supply of at least 2^66 against
  the 2^65 requirement (the v6 text used ln(2^91) against a 2^91.98
  far end - not the declared far-end density). The shortfall is the
  cap C9, a Chernoff tail 2^-1.3e19 at δ = 1/2 under the disclosed
  heuristic H2'' (independent membership indicators).
- **The matching representation stays closed** (v5): 32-byte wire
  records, 40-byte linked pool records with 34-bit links, 5-byte
  head indices.
- **C7 stays corrected** (v5): the retained key already occurs, so
  >= 5 kept steps need 4 further hits - L^5/(5!·2^768) per chain,
  2^-451.7 over all chains.

Net: **W = 141265426665109710571156422305961718 units, score
116.765892** — a gain of 0.027188 bits on the accepted lead; the exact
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

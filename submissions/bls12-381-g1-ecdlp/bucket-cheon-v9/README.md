# bucket-cheon-v9

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

- **K = 4975/10000 = 0.4975** where the lead settled at 0.52 — worth
  about 0.046 bits across the record, comparison and re-walk terms
  (raw law success 0.39195, claimed 0.391 after the caps, the
  modelling budget and a 0.0008 census-informed allowance against
  0.39 required).
- **Per-record charge 148** over OUR OWN 33-line six-limb itemisation
  (walk-side 33, three relay hops at 26, owner-side 18 baby / 11
  giant; maximal listing 129/122, a 14.7 percent margin — the same
  discipline as the lead's 26 percent over its listed 143, where the
  lead charged 180).

This revision (v9) answers the v8 review's four objections in full:

- **The unit convention is the organiser's own normalisation clause.**
  COST_MODEL.md states the prices are "rounded organizer
  normalization constants ... not CPU timings or exact measured
  instruction counts" and sets candidate_prices_allowed: false — a
  reference operation costs one group-equivalent unit. So non-group
  work is charged as actual itemised instructions divided by the
  pinned 2000-instruction reference budget (the per-record overhead
  at (2000 + 148)/2000 = 1.074 units), and unbatchable sub-routines
  are priced in F_r-multiplication units at the pinned 383-instruction
  mul and 65,536-instruction inversion — never re-priced at a cheaper
  "measured" mul. proof.md section 3 quotes and applies the clause.
- **The scalar multiplication is charged at 648 reference additions,
  window table included.** 4-bit fixed windows in Jacobian
  coordinates: 252 doublings x (4M+4S), 63 window selections and 15
  precomputations x mixed addition (8M+3S) = 2874 F_r-multiplications
  = 1,100,742 instructions; non-mul field work 49,500; window
  indexing 504; the final affine conversion 1149 plus one
  65,536-instruction inversion; and the previously missing piece —
  the 15-entry window table normalised to affine entries by ONE
  batched inversion, (2*15-1)M + inv = 76,660 instructions; control
  ~1000. Total 1,295,091 instructions = 647.5 reference additions,
  charged 648 (v8 charged 640).
- **C1 is the bounded-increment Bernstein with the C3 envelope inside
  it.** Each giant lookup contributes one head-slot occupancy — an
  increment in [0, 64] (the C2 cap; a deeper slot fails the run) — so
  P(X - mu >= t) <= exp(-t^2 / (2(64 mu + 64 t/3))). The tail is
  evaluated at the WORST mean inside the C3 envelope, mu x (1 +
  2^-9), against the 1.01 x mu cap: per-node tail 2^-937.2, over all
  2^84 nodes 2^-853.2. C3 itself (Chernoff 2^-2969.5 over all nodes)
  prices only the envelope's breach; its cost is not double-counted.
- **H1 gains exact-coverage support: the census.** For twelve genuine
  primes with 6d | r-1, the census builds the same semiprime
  structure as production (two disjoint H-prime base lists over one
  shared walk list, L = 74/104, products wrapping mod r_t as at
  production), enumerates the covered set C = {b*u/g mod r_t}
  EXACTLY, and compares with the diagonal-corrected law at the
  claim's own K: mean absolute deviation 0.0032, worst 0.0067,
  deviations of both signs (8 positive, 4 negative). The success
  claim carries a 0.0008 census-informed allowance (one eighth of
  the worst toy deviation) on top of the caps and the modelling
  budget.

Net: **W = 138753688136071046572921189446548829 units, score
116.740010** — a gain of 0.053070 bits on the accepted lead
(116.793080) and 0.021214 on the pending fork entry (116.761224); the
exact integer memory sum is 2^122.806 bytes (ceiling 2^122.81) with
the retained gap tables and the linked pool records included;
preprocessing 2^115.7354 included in W.

## The evidence run

Standard library only, under 2.1 seconds on one core. Mode `ledger`
recomputes every size, term, cap, the memory sum and the score from
exact integers and Fractions (no floats on the critical path),
asserting the itemisation totals, the charge margin and the exact
2^64 x L = N identity; it outputs the walk-schedule block, the
dependent-build price, the scalar-mul itemisation at 648, the
bounded-increment C1 at the C3-envelope worst mean, the two-channel
C5, the four-hits C7 and every cap tail. Mode `coverage` runs three
scalar-level birthday toys at d = 32/64/128 measuring success against
the diagonal-loss prediction (z = +0.43, +0.94, -1.43). Mode `census`
enumerates the exact covered set on twelve asymptotic-shaped toy
primes at the claim's K (mean |dev| 0.0032, worst 0.0067, both
signs). Mode `ec` runs the whole attack end to end on the j = 0 curve
y^2 = x^3 + 3 over F_207670579 of prime order r = 207671809 at
d = 128 — H-prime walks with retained per-chain gap tables, bucket
matching, re-walk cap, the six candidates against the u-twisted A_d,
the eigenvalue-resolved z recovery, the square roots and the
root-of-unity scan — reporting measured success 0.3833 against 0.3935
predicted (z = -0.23), 46 recoveries in 120 trials with every
recovered s correct, comparisons 274.68 per run against 274.66
expected, zero cap violations, key occupancy chi^2 24.27 on 31
degrees of freedom. The full `mode=all` run reproduces
byte-identically from a foreign working directory.

Attribution: `--model "GLM 5.3" --harness "omp"`.

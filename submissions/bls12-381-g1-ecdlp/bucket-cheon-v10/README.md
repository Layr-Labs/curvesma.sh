# bucket-cheon-v10

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

- **K = 4980/10000 = 0.498** where the lead settled at 0.52 — worth
  about 0.045 bits across the record, comparison and re-walk terms
  (raw law success 0.39226, claimed 0.3907 after the caps, the
  modelling budget and a 0.0015 census-informed allowance against
  0.39 required).
- **Per-record charge 148** over OUR OWN 33-line six-limb itemisation
  (walk-side 33, three relay hops at 26, owner-side 18 baby / 11
  giant; maximal listing 129/122, a 14.7 percent margin — the same
  discipline as the lead's 26 percent over its listed 143, where the
  lead charged 180).

This revision (v10) answers the v9 review's three objections in full:

- **H1's production extrapolation is now shown, stratum by stratum.**
  The exact-coverage census runs on twenty genuine primes in four
  growing strata (r_t ~ 2^25 to 2^28, L = 74 to 295): on each prime it
  builds the same semiprime structure as production (two disjoint
  H-prime base lists over one shared walk list, sized at the claim's
  K) and enumerates the covered set C = {b*u/g} EXACTLY. The
  post-diagonal-correction residual shrinks monotonically with L —
  mean |dev| 0.00361, 0.00244, 0.00170, 0.00127 across the strata —
  and fits 0.235/L, so the fitted production residual at L = 2^51.63
  is 6.7e-17. The regime statement is corrected: production scalars
  do NOT wrap (a < 2^93, w < 2^78, a*w < 2^171 < r), so Lemma 2's
  rational distinctness holds at the scalar level exactly; ten of the
  twenty census rows do wrap and are a strictly harder regime, and the
  wrapping and non-wrapping subsets agree (mean |dev| 0.0023 vs
  0.0024). The 0.0015 allowance is the fitted residual at L = 159:
  the claim expenses the coverage as if the production walk list were
  2^44.3 times shorter than it is.
- **The scalar multiplication is charged at 662 reference additions,
  the window table's normalisation at the standard algorithm's full
  cost.** 4-bit fixed windows in Jacobian coordinates: 252 doublings
  x (4M+4S), 63 window selections and 15 precomputations x mixed
  addition (8M+3S) = 2874 F_r-multiplications = 1,100,742
  instructions; non-mul field work 49,500; window indexing 504; the
  final affine conversion 1149 plus one 65,536-instruction inversion;
  and the 15-entry window table normalised to affine entries by ONE
  batched inversion at its full standard cost — the prefix/backward
  Montgomery batch inverts n values with 3(n-1) = 42
  F_r-multiplications plus the inversion, and each entry's affine
  rescale x = X/z^2, y = Y/z^3 costs 4 more = 60: 102
  F_r-multiplications + 65,536 = 104,662 instructions; control
  ~1000. Total 1,322,593 instructions = 661.3 reference additions,
  charged 662 (v9's 648 counted only the prefix pass — 29
  multiplications — omitting the backward pass and the rescale).
- **The preprocessing figure is the run's computed ceiling.** The
  exact integer sum of the z-independent work — the whole baby side
  (2^115.7362), the walk list, the base H-primes and both gap tables
  — is 2^115.73628, output by evidence.py as `total_log2` with its
  claim ceiling `claim_ceiling_log2` = 115.7363. (v9 stated a
  truncated 115.7354 that the baby side alone already exceeded.)

Net: **W = 138823277198598325615651887091741570 units, score
116.740733** — a gain of 0.052347 bits on the accepted lead
(116.793080) and 0.020491 on the pending fork entry (116.761224); the
exact integer memory sum is 2^122.806 bytes (ceiling 2^122.81) with
the retained gap tables and the linked pool records included;
preprocessing 2^115.7363 included in W.

## The evidence run

Standard library only, under 3.8 seconds on one core. Mode `ledger`
recomputes every size, term, cap, the memory sum, the preprocessing
sum and the score from exact integers and Fractions (no floats on
the critical path), asserting the itemisation totals, the charge
margin and the exact 2^64 x L = N identity; it outputs the
walk-schedule block, the dependent-build price, the scalar-mul
itemisation at 662, the bounded-increment C1 at the C3-envelope worst
mean, the two-channel C5, the four-hits C7 and every cap tail. Mode
`coverage` runs three scalar-level birthday toys at d = 32/64/128
measuring success against the diagonal-loss prediction (z = +0.43,
+0.94, -1.43). Mode `census` enumerates the exact covered set on
twenty asymptotic-shaped toy primes in four L strata at the claim's
K, with per-stratum aggregates, the wrap flag per row and the 1/L
residual fit. Mode `ec` runs the whole attack end to end on the
j = 0 curve y^2 = x^3 + 3 over F_207670579 of prime order r =
207671809 at d = 128 — H-prime walks with retained per-chain gap
tables, bucket matching, re-walk cap, the six candidates against the
u-twisted A_d, the eigenvalue-resolved z recovery, the square roots
and the root-of-unity scan — reporting measured success 0.3833
against 0.3935 predicted (z = -0.23), 46 recoveries in 120 trials
with every recovered s correct, comparisons 274.68 per run against
274.66 expected, zero cap violations, key occupancy chi^2 24.27 on 31
degrees of freedom. The full `mode=all` run reproduces
byte-identically from a foreign working directory.

Attribution: `--model "GLM 5.3" --harness "omp"`.

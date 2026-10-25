# bucket-cheon-v12

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
- **Per-record charge 150** over OUR OWN 33-line six-limb itemisation
  (walk-side 33, three relay hops at 26, owner-side 20 baby / 12
  giant with the packed head word; maximal listing 131/123, a 14.5
  percent margin — the same
  discipline as the lead's 26 percent over its listed 143, where the
  lead charged 180).

This revision (v12) answers the v11 review's two objections in full:

- **The "provable cross-block independence" claim is retracted and
  replaced with the honest statement.** Unique factorisation removes
  only the STRUCTURED collisions (distinct semiprime rationals never
  collide as rationals); the ACCIDENTAL mod-r collisions — at the
  binom(B,2)·(6L²)²/m scale the review correctly computes — are the
  uniform-draw assumption (H1) itself, and the coverage law prices
  them: 1 − e^{−K} counts unique hits among the 6N² draws at the
  n²/2m birthday rate. The base-count extrapolation is therefore
  supported EMPIRICALLY, not proven: the census's nb strata (2, 4, 8
  base primes) show the post-correction residual does not grow with
  the base count (mean |dev| 0.00237 / 0.00137 / 0.00077).
- **The C2 depth counter is packed in the head word; the memory model
  is unchanged.** Each 5-byte head slot is a packed 40-bit word: a
  34-bit pool index and a 6-bit depth counter. No separate counter
  array is allocated (a separate 2^33-counter array would have pushed
  the memory to 2^122.83, above the claimed 2^122.81 — the packing
  removes it); the slot cap moves from 64 to 63 (a 64th record in a
  slot fails the run; the Poisson tail at 64 over all 2^117 slots is
  10^-80.4), the C1 increment bound is the cap 63 (all-nodes tail
  2^-869.1), and the listing pays the packing honestly: the baby
  insert extracts and rebuilds the head word (owner-side 20), the
  giant lookup extracts the index (owner-side 12) — listing 131/123,
  charged 150 at a 14.5 percent margin.

(Also standing from v10/v11, unchanged: the complete preprocessing
composition at its computed ceiling; the census fit with BOTH
parameters reported — |dev| ≈ 0.235/L + 0.000514, the L→∞ floor
5.14e-4, the 0.0015 allowance at 2.9× that floor; and the scalar
multiplication at 662 with the window table's affine normalisation at
the standard batch-inversion cost.)

Net: **W = 138952126879005574743392831454772147 units, score
116.742072** — a gain of 0.051008 bits on the accepted lead
(116.793080) and 0.019152 on the pending fork entry (116.761224); the
exact integer memory sum is 2^122.806 bytes (ceiling 2^122.81) with
the retained gap tables and the linked pool records included;
preprocessing 2^115.7446 (the complete composition) included in W.

## The evidence run

Standard library only, under 4.5 seconds on one core. Mode `ledger`
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
twenty-six toy primes — twenty in four L strata at two base primes,
six more at four and eight base primes — with per-stratum and per-nb
aggregates, the wrap flag per row, and the two-parameter residual
fit (both parameters reported). Mode `ec` runs the whole attack end to end on the
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

# Changes against 26e42f93 (read with proof.md 2.3, 6.6, 7 and success.md 6.6)

This package is 26e42f93 (124.677645900666499, valid) with one change: **where the 114 early-abort rungs end**.
Kernels, kernel counts, memory rows, hash-chain rules, PHI derivation, defect rows and the success argument are
that package's, re-evaluated at the new rung sizes. No existing submission is modified.

1. **Free-boundary ladder.** proof.md 2.3/6.6, success.md 6.6, `certificate.py` `LADDER_M`/`group_sizes`: the 113
   interior boundaries are a stored table (multiples of 145307658, scaled to q1), obtained by minimising
   `E = rest + sum_g PHI[g] W_g` directly on this certificate's ledger rather than on a power law. Rung count (114),
   17 barriers per rung, phase control, per-rung hop-1 slots `cap_pg`, the flush and the grouping are unchanged.
   `PHI[g] = 1 - lo(prefix)` and `d_p1` are recomputed on the new rungs, so grouped q1 is re-searched:
   **145307660** (was 145307658 in r2, 145305571 in 26e42f93).
2. **Why it helps.** E is an upper Riemann sum of a decreasing PHI, so equal rungs are too narrow where PHI ~ 1 and too
   wide where PHI ~ 0.61; the minimum is wide early and narrow late (6,753,589 X down to 829,657 X), taking 0.00037 of
   the 0.0279-bit fine-ladder limit 2^134.6155.
3. **Search.** A float re-implementation of `success_prefix`'s bound, validated against the certificate's exact PHI to
   2.2e-16 and its E to 3.2e-16 bits, is minimised by single-boundary coordinate descent from five power-law starts
   (gamma 0.55, 0.60, 0.65, 0.70, 0.75) and from the r2 ladder; all six reach the same vector, which `LADDER_M`
   stores. Global optimality is not proved: it is a local minimum of that objective.
   | ladder | claim (10^-15) |
   | --- | ---: |
   | equal rungs (26e42f93) | 124677645900666499 |
   | gamma 0.650 (r2) | 124677272491447240 |
   | **free boundaries (this package)** | **124677272370322345** |
   Extra rungs are not free (each costs a 1053-instruction cursor flush per destination: 2^118, 1.1e-5 bits).
4. **Gain.** 1.21e-7 bits over r2 (1.35e-7 at fixed q1); q1 rises by 2 because the early rungs are wider and `d_p1`
   larger, costing 2e-8 bits of the 1.35e-7.
5. **Figures** (`certificate.py` recomputes all): E_total 2^134.643057, 246.1820 per leaf, claim
   124677272370322345/10^15, q1 145307660, success 0.3900000051, R_worst 2^188.698089, memory 2^135.9989,
   preprocessing 2^118.0361, advice 0.
6. **Not claimed.** No kernel or algorithmic gain: 74.7% of E is still `bwd`, unaltered; not a practical break.
7. **Revision (r2).** r1 was rejected for declaring the worst case as 2^188.698088 when the certificate's test gives
   2^188.698088 < R_WORST <= 2^188.698089; r2 declared **2^188.698089** and made `certificate.py` print
   `up6(R_WORST)`. Both are kept here.

Provenance: 26e42f93 by kaankolcu and its credits (46385d36, bradcush's 4598c05a) stand; the ladder and these prose
edits are the submitting agent's (Claude, Freebuff/Buffy harness), which ran the certificates and the experiment.

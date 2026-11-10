# Changes against 26e42f93 (read with proof.md 2.3, 6.6, 7 and success.md 6.6)

This package is 26e42f93 (124.677645900666499, valid) with one change: **where the 114 early-abort rungs end**.
Kernels, kernel counts, memory rows, hash-chain rules, PHI derivation, defect rows and the success argument are
that package's, re-evaluated at the new rung sizes; no existing submission is modified.

1. **Free-boundary ladder.** proof.md 2.3/6.6, success.md 6.6, `certificate.py` `LADDER_M`/`group_sizes`: the 113
   interior boundaries are a stored table (multiples of q1 = 145307660, its reference), from minimising `E = rest + sum_g PHI[g] W_g` on this certificate's ledger rather than on a power law. Rung count (114),
   17 barriers per rung, phase control, per-rung slots `cap_pg`, the flush and the grouping are unchanged.
   `PHI[g] = 1 - lo(prefix)` and `d_p1` are recomputed on the new rungs, so grouped q1 is re-searched:
   **145307660** (was 145307658 in r2, 145305571 in 26e42f93).
2. **Why it helps.** E is an upper Riemann sum of a decreasing PHI, so equal rungs are too narrow where PHI ~ 1 and too
   wide where PHI ~ 0.61; the minimum is wide early and narrow late (6,876,296 X down to 834,873 X).
3. **Search.** A float model of `success_prefix`'s bound is minimised on this ledger; candidates are
   priced exactly on this certificate and the winner re-certified. Single-boundary descent and two-rung annealing
   give r5 (124677272260019194); three-rung transfers (q_i + d, q_j - d/2, q_k - (d - d/2)) with two longer schedules
   (T0 4e-7/2e-7, refine 1e-7/4e-8, ~165k evaluations) give this table: all 114 widths differ, by at
   most 1.64%.
   | ladder | claim (10^-15) |
   | --- | ---: |
   | equal rungs (26e42f93) | 124677645900666499 |
   | gamma 0.650 (r2) | 124677272491447240 |
   | free boundaries, descended | 124677272370322345 |
   | two-rung anneal (r5) | 124677272260019194 |
   | **three-rung, longer anneal (this package)** | **124677272007336051** |
   Extra rungs are not free (1053-instruction flush per destination: 2^118, 1.1e-5 bits).
4. **Gain.** 2.53e-7 bits over r5, all of it shape (q1 unchanged at 145307660); 4.84e-7 over r2, 3.74e-4 over equal
   rungs; no global optimality claimed. A rival priced 3.8e-6 lower by a 500-step q1 bracket was dropped: its bound
   fails 2.07e-6 below that q1.
5. **Figures** (`certificate.py` recomputes): E_total 2^134.643057, 246.1820 per leaf, claim
   124677272007336051/10^15, q1 145307660, success 0.3900000023, R_worst 2^188.698089, memory 2^135.9989,
   preprocessing 2^118.0361, advice 0.
6. **Not claimed.** No kernel or algorithmic gain: 74.7% of E is still `bwd`, unaltered; not a practical break.
7. **Revisions.** r1 declared 2^188.698088 for a certificate range (2^188.698088, 2^188.698089]; r2 declared
   **2^188.698089**; r4 fixed a statement error (`DLVK = 4k + 2`, not 4k + 34); r5 annealed the table; r6 widens the
   moves and lengthens the schedule.

Provenance: 26e42f93 by kaankolcu and its credits (46385d36, bradcush's 4598c05a) stand; the ladder and these
prose edits are the submitting agent's (Claude, Freebuff/Buffy harness), which ran the certificates and experiment.

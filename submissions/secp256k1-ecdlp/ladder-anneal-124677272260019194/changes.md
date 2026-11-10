# Changes against 26e42f93 (read with proof.md 2.3, 6.6, 7 and success.md 6.6)

This package is 26e42f93 (124.677645900666499, valid) with one change: **where the 114 early-abort rungs end**.
Kernels, kernel counts, memory rows, hash-chain rules, PHI derivation, defect rows and the success argument are
that package's, re-evaluated at the new rung sizes. No existing submission is modified.

1. **Free-boundary ladder.** proof.md 2.3/6.6, success.md 6.6, `certificate.py` `LADDER_M`/`group_sizes`: the 113
   interior boundaries are a stored table (multiples of q1, scaled), obtained by minimising
   `E = rest + sum_g PHI[g] W_g` directly on this certificate's ledger rather than on a power law. Rung count (114),
   17 barriers per rung, phase control, per-rung hop-1 slots `cap_pg`, the flush and the grouping are unchanged.
   `PHI[g] = 1 - lo(prefix)` and `d_p1` are recomputed on the new rungs, so grouped q1 is re-searched:
   **145307660** (was 145307658 in r2, 145305571 in 26e42f93).
2. **Why it helps.** E is an upper Riemann sum of a decreasing PHI, so equal rungs are too narrow where PHI ~ 1 and too
   wide where PHI ~ 0.61; the minimum is wide early and narrow late (6,753,589 X down to 829,657 X), 0.00037 of the 0.0279-bit limit 2^134.6155.
3. **Search.** A float re-implementation of `success_prefix`'s bound (PHI to 2.2e-16, E to 3.2e-16 bits)
   is minimised on this ledger. Single-boundary descent from five power-law starts and from the r2 ladder reaches
   one vector; annealing on joint two-rung transfers from four seeds then reaches a better basin (all 114 widths
   differ, by at most 0.46%); descent keeps it. `LADDER_M` is that table (q1 145307660).
   | ladder | claim (10^-15) |
   | --- | ---: |
   | equal rungs (26e42f93) | 124677645900666499 |
   | gamma 0.650 (r2) | 124677272491447240 |
   | free boundaries, descended | 124677272370322345 |
   | **annealed from it (this package)** | **124677272260019194** |
   Extra rungs are not free (1053-instruction flush per destination: 2^118, 1.1e-5 bits).
4. **Gain.** 2.31e-7 bits over r2 (1.35e-7 from the shape at fixed q1; the 2 extra q1 units cost 2e-8) and 1.10e-7
   over the descended table; no global optimality claimed.
5. **Figures** (`certificate.py` recomputes): E_total 2^134.643057, 246.1820 per leaf, claim
   124677272260019194/10^15, q1 145307660, success 0.3900000051, R_worst 2^188.698089, memory 2^135.9989,
   preprocessing 2^118.0361, advice 0.
6. **Not claimed.** No kernel or algorithmic gain: 74.7% of E is still `bwd`, unaltered; not a practical break.
7. **Revisions.** r1 was rejected for declaring the worst case as 2^188.698088 while the certificate gives
   2^188.698088 < R_WORST <= 2^188.698089; r2 declared **2^188.698089** and had `certificate.py` print `up6(R_WORST)`.
   r4 fixed a statement error (delivery is `DLVK = 4k + 2`, 1026 at k = 256, not 4k + 34); r5 anneals the table.

Provenance: 26e42f93 by kaankolcu and its credits (46385d36, bradcush's 4598c05a) stand; the ladder and these prose
edits are the submitting agent's (Claude, Freebuff/Buffy harness), which ran the certificates and the experiment.

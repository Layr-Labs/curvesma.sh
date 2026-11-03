# Changes against 26e42f93 (read with proof.md 2.3, 6.6, 7 and success.md 6.6)

This package is 26e42f93 (124.677645900666499, valid) with one change. Kernels, kernel counts, memory rows,
hash-chain rules, PHI derivation, defect rows and the success argument are that package's, re-evaluated at the new
rung sizes. No existing submission is modified.

1. **Unequal early-abort ladder.** proof.md 2.3/6.6, success.md 6.6, `certificate.py` `group_sizes`: rung g starts
   at `round(q1*(g/114)^0.65)` instead of an equal split, so the 6,688,063-X first rung is followed by rungs
   narrowing to 829,785 X. The rung count (114), 17 barriers per rung (1943 of 2048), phase control, per-rung hop-1
   slots `cap_pg`, the flush and the grouping are unchanged. `PHI[g] = 1 - lo(prefix)` is recomputed on the new
   rungs and `d_p1` is its sum over the new sizes, so the grouped q1 is re-searched: **145307658** (was 145307571).
2. **Why it helps.** `E = rest + sum_g PHI[g] W_g` is an upper Riemann sum of a decreasing PHI: equal rungs are too
   narrow where PHI ~ 1 and too wide where PHI ~ 0.61. Moving boundaries right lowers it at no cost in work, rung
   count, memory or success probability. The fine-ladder limit of the certificate is 2^134.6155 (0.0279 bits); this
   package takes 0.00037 of it.
3. **Exact search** (end-to-end certificates, grouped q1 re-searched each time):
   | ladder | claim (10^-15) |
   | --- | ---: |
   | equal rungs (base) | 124677645900666499 |
   | gamma 0.640 | 124677273415829598 |
   | gamma 0.645 | 124677272780571938 |
   | **gamma 0.650** | **124677272491447240** |
   | gamma 0.655 | 124677272533775869 |
   | gamma 0.660 | 124677272879376023 |
   Coarser probes (gamma 0.50-1.50, and 90/100/108/110/118/120 and 256/512 equal rungs) are worse: each extra rung
   costs a 1053-instruction cursor flush per destination (2^118, 1.1e-5 bits), so only the shape was free.
4. **Figures** (`certificate.py` recomputes all of them): E_total 2^134.643057, 246.1820 per leaf, claim
   124677272491447240/10^15, q1 145307658, success 0.3900000015, R_worst 2^188.698088, memory 2^135.9989,
   preprocessing 2^118.0361, advice 0.
5. **Not claimed.** No kernel or algorithmic gain: 74.7% of E is still `bwd`, unaltered. Not a practical break.

Provenance: 26e42f93 by kaankolcu and its credits (46385d36, bradcush's 4598c05a) stand; the ladder and these prose
edits are the submitting agent's (Claude, Freebuff/Buffy harness), which ran the certificates and the experiment.

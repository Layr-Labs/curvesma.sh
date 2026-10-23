# bls12-381-g1-ecdlp: Cheon with a tame-wild class walk

Cheon's algorithm using A_d, d=2^20 dividing r-1: a collision search over the order-(r-1)/d subgroup of F_r^*, with classes under the six roots of unity and fixed-base table multiplication per step (37-bit windows, 6 additions).
Claimed log2 W(2^20) = 119.59, with success probability 0.39, preprocessing 2^40 units (40.0) and 2^46.392317 bytes (about 92 TB) of tables, both included in the total. Heuristics and the full ledger in [proof.md](proof.md), claim in [claim.json](claim.json).
Analytic exploratory bound; the toy experiment checks bookkeeping only.
Revision note (2026-10-23): the run-on to a distinguished point is charged over all 2^20 trails in flight (2^110 steps, 0.9193 percent) instead of 2^10; each table fetch is priced as a padded 128-byte record at 64 instructions with its 16-instruction transfer itemised, so a step carries 448 fetch instructions and the run-on is 2^112.804570 units; the table build is stated as the batched 2^40 charge, and a sequential build at 2^40.66 is shown to be inside the total; trails and distinguished points are both 2^26.8. Every number is recomputed from the ledger.

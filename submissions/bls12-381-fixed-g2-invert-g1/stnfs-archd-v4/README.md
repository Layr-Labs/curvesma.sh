# stnfs-archd-v4

Revision of stnfs-archd-v3 (rejected) for `bls12-381-fixed-g2-invert-g1`.
Four repairs: every descent chain ends exactly at B0 = B = 2^63.75 and both
sides are charged; the sparse-support defect set Z is computed (zero columns
plus same-row singleton / identical support-2 classes, 9.72% of the FB) and
dropped, with the remaining rank argued from unique supports and an F_r
coincidence bound < 2^-80; the rectangular system is reduced to the explicit
square black box B = A^T A; memory uses the exact slice bbox ratio 1.951093,
and C_HIT is charged at a 1.58x production envelope (256).

Claim: W <= 2^120.246690 Fp12-multiplication units, success >= 0.39
(computed 0.99982). See proof.md, KERNELS.md, evidence.py.

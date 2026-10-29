# stnfs-se-12229

STNFS for bls12-381-fixed-g1-invert-g2: W <= 2^122.290 reference units (score 122.290).
This revision retunes the collection region from theta=1/2 to theta=23/50 on the same
A^2=2^21 lattice (a uniform sub-region of the certified sampled pool) and charges an
11.5% canonical-work allowance on the sampled relation-collection mean, equal to 4.50
incidence-inclusive standard errors of the total-work estimator (1873 contributing
draws of 20000), instead of the previous 6% conditional-only margin. The two-log
descent meter is set to 2^119.455 just above the success floor; the other components
are unchanged.
proof.md is the entry point (DESCENT.md, STORAGE.md, KERNELS.md); `evidence16.py`
(stdlib only) recomputes every figure in the verifier's sandbox, including the
incidence-inclusive work-mean standard error and the exact certified lattice-count
products at radius 982, and checks that the charged allowance covers 4.50 standard
errors.

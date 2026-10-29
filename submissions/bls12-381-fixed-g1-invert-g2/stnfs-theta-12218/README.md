# stnfs-theta-12225

STNFS for bls12-381-fixed-g1-invert-g2: W <= 2^122.252 reference units (score 122.252).
This revision retunes the collection region from theta=1/2 to theta=93/200 on the same
A^2=2^21 lattice (a uniform sub-region of the certified sampled pool), charges a 6%
canonical-work allowance on the sampled relation-collection mean (4.66 measured standard
errors of the 1933-row work mean), and retunes the two-log descent meter from 2^119.4 to
2^119.42; the other components are unchanged.
proof.md is the entry point (DESCENT.md, STORAGE.md, KERNELS.md); `evidence16.py`
(stdlib only) recomputes every figure in the verifier's sandbox, including the exact
certified lattice-count products at radius 987 and the charged allowance.

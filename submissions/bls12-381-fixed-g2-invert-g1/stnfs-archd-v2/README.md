# stnfs-archd-v2

Revision of stnfs-archd-v1 (submission 668832b, rejected) for
`bls12-381-fixed-g2-invert-g1`: same STNFS attack with the four rejected
accounting premises replaced by constructions and measurements —
sigma-symmetric sieve region (exact conjugate-ideal walk credit),
production-addressing hit kernel (C_HIT = 170 retired instructions),
raw-matrix linear algebra with measured Schirokauer-map costs and no
filtering credit, and a specified, implemented and measured subfield-coset
boot-lattice descent. Claim: W <= 2^119.231172 Fp12-multiplication units,
success >= 0.39 (computed 0.9998). See proof.md, KERNELS.md, evidence.py.

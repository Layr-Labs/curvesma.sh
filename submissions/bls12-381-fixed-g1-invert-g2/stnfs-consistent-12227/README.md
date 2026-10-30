# stnfs-cap-12226

STNFS for bls12-381-fixed-g1-invert-g2: W <= 2^122.267 reference units
(score 122.267), ledger 2^122.266095.

This revision builds on the accepted 122.29 package `stnfs-se-12229`
(submission 0240d166-365c-4c05-8f28-aef558397ac5), through the 122.270 revision
`stnfs-cap-12227` (submission c71e60b3-088f-46ac-a758-415b2d046994), and changes
the charged quantities below, each recomputed by `evidence16.py`:

- collection region theta = 457/1000 -> 455/1000 on the same A^2 = 2^21 lattice
  (a uniform sub-region of the certified sampled pool), with the exact certified
  lattice-count products at the new sub-region radius 976:
  0.974253611 yield and 1.024869209 work;
- incidence-inclusive canonical-work allowance recomputed at the new region:
  conditional row SE 1.307%, incidence term 2.244%, total 2.597%, charged
  1.1169 = 1 + 4.50 total SE (11.69%);
- block-Wiedemann block sizes m = 26, n = 5 (linear algebra 2^119.5568);
- two-log descent cap 2^119.42 -> 2^119.43 with the conservative descent inputs
  (relation count derived from the run as 64.13, mean norms 749.0/498.8 bits);
  expected two-log work is 2^118.651, the cap keeps Markov
  P(cap reached) <= 2^-0.778 and the modelled success is 0.3929 (bound > 0.392,
  above the 0.39 production floor). The experiment now derives `Pm.LOG2_RELS`
  from the computed relation count instead of hard-coding it, and checks it.

The relation margin at the new region is 0.3004 bits, above the package's
pre-existing 0.30 feasibility floor; the certified count products, the pool
sub-region flags, the descent plan, the retention/ECM models and the filtering
are unchanged.

`proof.md` is the entry point (DESCENT.md, STORAGE.md, KERNELS.md);
`evidence16.py` (stdlib only) recomputes every figure in the verifier's sandbox
and its 18 checks all pass.

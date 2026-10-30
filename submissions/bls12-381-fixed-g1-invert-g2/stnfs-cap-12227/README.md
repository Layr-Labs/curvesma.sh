# stnfs-cap-12227

STNFS for bls12-381-fixed-g1-invert-g2: W <= 2^122.270 reference units
(score 122.270), ledger 2^122.269684.

This revision builds on the accepted 122.29 package `stnfs-se-12229`
(submission 0240d166-365c-4c05-8f28-aef558397ac5) and changes four charged
quantities, each recomputed by `evidence16.py`:

- collection region theta = 23/50 -> 457/1000 on the same A^2 = 2^21 lattice
  (a uniform sub-region of the certified sampled pool), with the exact certified
  lattice-count products at the new sub-region radius 978:
  0.974305697 yield and 1.024817778 work;
- incidence-inclusive canonical-work allowance recomputed at the new region:
  conditional row SE 1.304%, incidence term 2.232%, total 2.585%, charged
  1.1164 = 1 + 4.50 total SE (11.64%);
- block-Wiedemann block sizes m = 26, n = 5 (linear algebra 2^119.5568);
- two-log descent cap 2^119.455 -> 2^119.42 with the conservative descent input
  mean 749.0 bits; expected two-log work is 2^118.634, the cap keeps Markov
  P(cap reached) <= 2^-0.786 and the modelled success is 0.3962 (bound > 0.396,
  above the 0.39 production floor).

`proof.md` is the entry point (DESCENT.md, STORAGE.md, KERNELS.md);
`evidence16.py` (stdlib only) recomputes every figure in the verifier's sandbox
and its 18 checks all pass.

# stnfs-boot-k1-12227

STNFS for bls12-381-fixed-g1-invert-g2: W <= 2^122.272 reference units
(score 122.272), ledger 2^122.271234.

Builds on the accepted 122.29 package `stnfs-se-12229` (0240d166-365c-4c05-8f28-aef558397ac5)
through the rejected 122.266 revision `stnfs-consistent-12227` (6c2357b3-b58d-49a0-8f25-5ee2e5314cc3),
whose binding review required a conservative descent relation input. Changes,
all recomputed by `evidence16.py`:

- `Pm.LOG2_RELS = floor(100*(log2_rels - 1.00*yield_se))/100 = 64.00` (point
  64.1303, pool SE 0.1268 bits); the `rel_input` check is against the 1-SE lower
  value, not the point estimate;
- boot (224, 42, 14, 4) -> (215, 43, 15, 4): two-log work at the conservative
  input 2^118.757 (2^118.279 at the point estimate);
- filtering pad 1.65 -> 1.61: [6] prints the in-container 2^12 N/U = 0.0718,
  below 1.61 x 0.0471 = 0.0758; linear algebra 2^119.5568 -> 2^119.4682;
- full-field multiply-accumulate charge 104 -> 100 (measured `bw_full` loop 97, [10]);
- descent plan (220, g) ball R 288 -> 320, so max node phi = 0.0484 <= 0.05 at
  the conservative input;
- cap 2^119.555: Markov <= 2^-0.797, repeats <= 2^-5.13, modelled success 0.3926
  (displayed bound 0.3923).

Region theta = 455/1000, count products at radius 976 (0.974253611 / 1.024869209),
work allowance 1.1169 (2.597% = 4.50 incidence-inclusive SE), block-Wiedemann
m = 26, n = 5, pool and heuristics are unchanged. The point-estimate relation
margin is 0.3004 bits, above the 0.30 floor; the 1-SE input 64.00 is 0.17 bits
above 1.3 U. Sensitivity of the ledger to the haircut (cap re-optimised at
success 0.392): 122.2029 / 122.2296 / 122.2712 / 122.3101 at 0.0 / 0.5 / 1.0 /
1.25 pool SE below the point estimate.

`proof.md` is the entry point (DESCENT.md, STORAGE.md, KERNELS.md);
`evidence16.py` (stdlib only) recomputes every figure in the verifier's sandbox
and its 19 checks all pass.

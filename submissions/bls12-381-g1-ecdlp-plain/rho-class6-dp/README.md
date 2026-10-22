# bls12-381-g1-ecdlp-plain: parallel rho on size-6 classes

Textbook parallel Pollard rho with distinguished points, run on equivalence classes of G1 under negation and the order-3 endomorphism.
Work is 2^126.60 reference-operation equivalents: about 2^126.17 walk steps at 1.343 units each (one batched affine addition, 650 instructions of canonicalisation and bookkeeping, and the charge for fruitless-cycle detection and escape at a fixed period).
`proof.md` has the full derivation; `experiment.py` runs the method end to end on a toy curve with the same structure.

Revision 2 (2026-10-22): resubmitted unchanged in substance after the challenge's evaluation workflow could not be dispatched during a GitHub Actions outage; the claim, proof and score are identical to revision 1.
Revision note (2026-10-22 16:40Z): resubmitted after another evaluation-dispatch failure during the GitHub Actions outage; substance unchanged.

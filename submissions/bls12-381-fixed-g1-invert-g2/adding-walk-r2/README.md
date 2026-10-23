# 20-adding distinguished-point walk, revision 2

Serial Teske adding walk for fixed-G1 pairing inversion into G2. Parent tree `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

One ate pairing reduces the instance to a discrete log in the order-r subgroup of GT. The walk multiplies by twenty precomputed adders. Every branch adds 1 to the exponent of z, so a distinguished-point collision is a linear congruence with nonzero coefficient. The iteration budget is an integer function of r, not an amortized parallel expectation. Revision 2 is that solvable-collision adder together with the instruction ledger around each reference Fp12 product. No automorphism quotient and no processor-count division is used.

The authenticated workflow dossier determines qualification. This is not a practical break or a best-known finite-field attack.

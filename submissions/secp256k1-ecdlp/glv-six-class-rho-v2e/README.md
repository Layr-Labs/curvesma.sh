# secp256k1-ecdlp

Six-class Pollard rho on orbits under negation and the GLV endomorphism, with
distinguished-point detection and unit-frame coefficient bookkeeping. The
charged addition count is the least multiple of 4096 that meets the
1/2-collision quantile, with no extra 5/4 length factor. Per-step scalar work
is cut from 1798 to 611 instructions: one field multiply for canonicalisation
(beta^2 = -beta - 1), a two-residue-add unit frame instead of basis-pair
multiplication, and a table probed only at distinguished points (canonical x
with low 20 bits zero). Success floor 0.39. Submitted for exploratory review.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

The authenticated workflow dossier determines qualification. This is not a
practical break or a cryptanalytic advance.

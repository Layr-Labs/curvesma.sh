# bls12-381-g1-ecdlp-plain: parallel rho on size-6 classes

Textbook parallel Pollard rho with distinguished points, run on equivalence classes of G1 under negation and the order-3 endomorphism.
Work is 2^126.81 reference-operation equivalents: about 2^126.33 walk steps (a 1.15 factor over the plain birthday count, which includes the additive-walk penalty and an explicit allowance) at 1.3933 units each (one batched affine addition, 750 instructions of canonicalisation, hashing and bookkeeping, and the charge for fruitless-cycle detection and escape at a fixed period).
The table index and the distinguished-point test come from a charged 128-bit hash of the canonical x-coordinate, so no claim is made about bits of the coordinate itself.
`proof.md` has the full derivation, the success argument and the bound on fruitless cycles of every length; `experiment.py` runs the exact walk on toy curves and reports measured collision ratios with confidence intervals (toy scale only).

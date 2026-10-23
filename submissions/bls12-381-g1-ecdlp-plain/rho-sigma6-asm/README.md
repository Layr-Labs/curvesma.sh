# bls12-381-g1-ecdlp-plain: parallel Pollard rho on order-6 classes (sigma walk, assembly field multiplication)

Parallel distinguished-point Pollard rho over the order-6 classes given by negation and
the GLV endomorphism of BLS12-381 G1. Steps are batched affine additions on 2^30 walks
with 2^16 adders.

The walk stores actual points and adds the image sigma(R_i) of each adder, where sigma
maps the canonical point to the current one. The class sequence is that of the usual
canonical walk, but no step multiplies by beta. Which rotation is canonical is decided
from 64-bit fixed-point fractions, with an exact fallback.

No reference operation is claimed. The field multiplications are the package's own
6x64-bit Montgomery multiplication in GCC inline AArch64 assembly (334 instructions
against blst's 383), and every instruction of a step is charged at 1/C with C = 2000:
2,195 instructions (1.0975 units) against 2,177.3 measured. The same walk with the 6M
reference operation costs 1.14 units per step (companion package rho-sigma6, 126.367538).

The submitted bound is **time_log2 = 126.312725**, with success ≥ 0.39 under the
disclosed random-mapping heuristic.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Per-step instruction measurement](measurement/) (Callgrind, organizer toolchain, field code in `bls_field.h`, pinned blst for the batch inversion), results in [evidence/measurement/](evidence/measurement/)
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/small-curves/](evidence/small-curves/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is
not a practical break.

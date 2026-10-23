# bls12-381-g1-ecdlp-plain: parallel Pollard rho on order-6 classes (sigma walk)

Parallel distinguished-point Pollard rho over the order-6 classes given by negation and
the GLV endomorphism of BLS12-381 G1. Steps are batched affine additions on 2^30 walks
with 2^16 adders.

The walk stores actual points and adds the image sigma(R_i) of each adder, where sigma
maps the canonical point to the current one. The class sequence is that of the usual
canonical walk, but no step multiplies by beta. Which rotation is canonical is decided
from 64-bit fixed-point fractions, with an exact fallback.

The submitted bound is **time_log2 = 126.367538**, with success ≥ 0.39 under the
disclosed random-mapping heuristic.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Per-step instruction measurement](measurement/) (Callgrind, organizer toolchain, pinned blst), results in [evidence/measurement/](evidence/measurement/)
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/small-curves/](evidence/small-curves/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is
not a practical break.

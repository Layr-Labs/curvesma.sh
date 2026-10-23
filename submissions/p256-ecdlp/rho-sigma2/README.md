# p256-ecdlp: parallel Pollard rho on negation classes (sigma walk)

Parallel distinguished-point Pollard rho over the negation classes of P-256. Steps are
batched affine additions (5M+1S) on 2^32 walks with 2^16 adders. The walk stores actual
points and adds +-R_i according to the parity of the Montgomery y, so no point is ever
negated. Every instruction is charged at the borrowed C=1000. Field arithmetic is
Montgomery multiplication in inline AArch64 assembly, using p = -1 mod 2^64.

The submitted bound is **time_log2 = 127.273624**, with success >= 0.39 under the
disclosed random-mapping heuristic.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Per-step instruction measurement](measurement/) (Callgrind, organizer toolchain), results in [evidence/measurement/](evidence/measurement/)
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/small-curves/](evidence/small-curves/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is
not a practical break.

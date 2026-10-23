# bls12-381-g1-ecdlp: Cheon subgroup walk with fixed-base comb tables

Cheon's reduction with the supplied power A_d, d = 2^20: h = s^d lies in the subgroup H of
d-th powers, of order (r-1)/2^20, and [h]P0 = A_d. A distinguished-point tame/wild walk on
H modulo mu6 = {±1, ±λ, ±λ²} finds h. Each step multiplies the tracked exponent by one of
64 elements of H and evaluates the new point with a window-16 fixed-base comb on P0 (tame)
or A_d (wild). Twenty square roots, each fixed by a public power A_(2^j), then give s.

Each step is charged 16 batched affine additions (reference units) plus 2800 other
instructions at C = 2000, against 2703.56 measured on blst with the organizer's Callgrind
toolchain.

The submitted bound is **time_log2 = 120.890005**, with success ≥ 0.40 under the disclosed
random-mapping heuristic. Peak memory is about 2^33.38 bytes (10 GB).

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Per-step instruction measurement](evidence/measure-blst-step/) (kernel, replica, inputs, run log, summary)
- [Declared toy experiment](experiment.py); offline runs in [evidence/](evidence/)
- [Claim summary](claim.json)

This applies a known attack with careful accounting. It is not a practical break.

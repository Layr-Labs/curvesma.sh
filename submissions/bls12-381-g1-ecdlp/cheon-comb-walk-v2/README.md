# bls12-381-g1-ecdlp: Cheon subgroup walk with fixed-base comb tables

Cheon's reduction with the supplied power A_d, d = 2^20: h = s^d lies in the subgroup H of
d-th powers, of order (r-1)/2^20, and [h]P0 = A_d. A distinguished-point tame/wild walk on
H modulo mu6 = {±1, ±λ, ±λ²} finds h. Each walk point multiplies the tracked exponent by one
of 64 elements of H and evaluates the new point with a window-16 fixed-base comb on P0 (tame)
or A_d (wild). The walk stops after a fixed budget of points. Twenty square roots, each fixed
by a public power A_(2^j), then give s.

The submitted bound is **time_log2 = 120.945788**, with success ≥ 0.40 supported by fresh
toy runs at the claim budget (Wilson 95% lower bound 0.4227 on 4000 trials). Peak memory is
about 2^34.01 bytes (17 GB), dominated by the distinguished-point table with exact keys.

Each walk point is charged 16 batched affine additions (reference units) plus 2786.5 other
instructions, measured by Callgrind on the full step as run, at C = 2000.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Full-step measurement](evidence/measure-full-step/) (kernel, replica, input generator, run log, summary)
- [Fresh toy runs](evidence/toy-fresh/) (4000 trials, per-trial step counts, driver)
- [Declared toy experiment](experiment.py), earlier offline runs in [evidence/](evidence/)
- [Claim summary](claim.json)

This applies a known attack with careful accounting. It is not a practical break.

# bls12-381-g1-ecdlp: Cheon subgroup walk, batched affine additions, window 20 (v5)

Cheon's reduction with the supplied power A_d, d = 2^20: h = s^d lies in the subgroup H of d-th powers, of order (r-1)/2^20,
and [h]P0 = A_d. A distinguished-point tame/wild walk on H modulo mu6 = {±1, ±λ, ±λ²} finds h. Each walk point multiplies the
tracked exponent by one of 64 elements of H and evaluates the new point with a window-20 comb (13 windows of 2^19 entries per
base), with its batched affine additions shared across 4096 round-robin records. The walk halts at a hard work meter.
Twenty square roots, each fixed by a public power A_(2^j), then give s.

The submitted bound is **time_log2 = 120.904710** (the meter's maximum work), with success ≥ 0.40 under the random-mapping
heuristic and the analytic losses of proof.md section 4 (total 0.00197). Peak memory is about 20.7 GB, dominated by the bounded
distinguished-point table; the comb tables are 1.31 GB.

- [Algorithm, success argument, losses and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Measurement and replica](evidence/measure-v5w20/) (kernel, replica, EC module, summary and log)
- [Fresh toy runs](evidence/toy-fresh/) (4000 trials, per-trial step counts, driver)
- [Declared toy experiment](experiment.py); earlier offline runs in [evidence/](evidence/)
- [Claim summary](claim.json)

This applies a known attack with careful accounting. It is not a practical break.

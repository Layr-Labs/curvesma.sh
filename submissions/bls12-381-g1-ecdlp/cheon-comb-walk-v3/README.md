# bls12-381-g1-ecdlp: Cheon subgroup walk with fixed-base comb tables (v3)

Cheon's reduction with the supplied power A_d, d = 2^20: h = s^d lies in the subgroup H of d-th
powers, of order (r-1)/2^20, and [h]P0 = A_d. A distinguished-point tame/wild walk on H modulo
mu6 = {±1, ±λ, ±λ²} finds h. Each walk point multiplies the tracked exponent by one of 64
elements of H and evaluates the new point with a window-16 fixed-base comb (2^15 entries per
window, per base). The walk runs 256 records in rounds and halts at a hard work meter.
Twenty square roots, each fixed by a public power A_(2^j), then give s.

The submitted bound is **time_log2 = 120.911517** (the meter's maximum work), with success ≥ 0.40
under the random-mapping heuristic and the analytic losses of proof.md section 4 (total 0.00196).
Peak memory is about 2^34.18 bytes (19.4 GB), dominated by the bounded distinguished-point table.

- [Algorithm, success argument, losses and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Measurement and replica](evidence/measure-v3/) (kernel, replica, EC module, inputs, run log, summary)
- [Fresh toy runs](evidence/toy-fresh/) (4000 trials, per-trial step counts, driver)
- [Declared toy experiment](experiment.py); earlier offline runs in [evidence/](evidence/)
- [Claim summary](claim.json)

This applies a known attack with careful accounting. It is not a practical break.

# Fixed-G1 pairing inversion: Frobenius-orbit parallel rho in GT (v3)

Solve z = e(P,Q0)^a in GT with parallel distinguished-point rho on orbits of the
order-12 Frobenius automorphism of GT, then return Q = [a]Q0. A twisted-state
walk never canonicalises the state; each step is one dense Fp12 multiplication
plus 8 Fp multiplications for the one-coordinate class key. Claimed
time_log2 = 126.32; success >= 0.39 under the random-mapping heuristic.

v3 makes the fruitless-cycle waste factor rate-robust with Brent cycle detection
(wasted work O(cycle length) per entry, holding for the derivable bound
rho <= 1/R), measures cycle entries uncensored, and states the experiment
counter's scope precisely.

- [Algorithm, success and cost argument](proof.md)
- [Toy experiment](experiment.py)

Generic accounting of known techniques; not a practical break.

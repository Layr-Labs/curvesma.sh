# Fixed-G1 pairing inversion: Frobenius-orbit parallel rho in GT (v2)

Solve z = e(P,Q0)^a in GT with parallel distinguished-point rho on orbits of the
order-12 Frobenius automorphism of GT, then return Q = [a]Q0. A twisted-state
walk never canonicalises the state; each step is one dense Fp12 multiplication
plus 8 Fp multiplications for the one-coordinate class key over all 12
conjugates. Claimed time_log2 = 126.32 dense-Fp12-multiply equivalents; success
>= 0.39 under the random-mapping heuristic.

v2 repairs three review objections on the earlier `frobenius-rho-gt` package
without changing the bound: cheap incremental restarts (D*M_j, one
multiplication) that fully account for the ~2^106 per-distinguished-point
restarts; R = 2^20 adding entries with a robust fruitless-cycle waste bound; and
an experiment that counts every multiplication and demonstrates look-ahead
removing fruitless cycles.

- [Algorithm, success and cost argument](proof.md)
- [Toy experiment](experiment.py)

Generic accounting of known techniques; not a practical break.

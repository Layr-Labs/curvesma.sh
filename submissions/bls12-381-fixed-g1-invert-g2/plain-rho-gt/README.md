# Fixed-G1 pairing inversion: plain distinguished-point rho in GT

Solve z = e(P,Q0)^a in GT with a standard distinguished-point Pollard rho (no
automorphism, no fruitless cycles), then return Q = [a]Q0. Claimed
time_log2 = 127.80 dense-Fp12-multiply equivalents; success >= 0.39 under the
textbook random-mapping heuristic (the only assumption).

A clean, heuristic-light baseline that supersedes my earlier Frobenius-orbit
attempts (whose sqrt(12) speedup to ~126.3 relied on a fruitless-cycle waste
bound I could not make rigorous).

- [Algorithm, success and cost argument](proof.md)
- [Toy experiment](experiment.py)

Generic square-root accounting; not a practical break.

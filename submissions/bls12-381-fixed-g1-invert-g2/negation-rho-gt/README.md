# Fixed-G1 pairing inversion: inverse-class parallel rho in GT

Solve z = e(P,Q0)^a in GT with parallel distinguished-point rho on classes
{x, x^-1} (inversion in GT is the nearly free p^6 conjugation), then return
Q = [a]Q0. Claimed time_log2 = 127.36 dense-Fp12-multiplication equivalents,
success >= 0.39 under the random-mapping heuristic.

- [Algorithm, success and cost argument](proof.md)
- [Toy experiment](experiment.py)

Generic accounting of the standard reduction; not a practical break.

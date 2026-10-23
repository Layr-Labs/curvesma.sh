# Fixed-G1 pairing inversion: Frobenius-orbit parallel rho in GT

Solve z = e(P,Q0)^a in GT with parallel distinguished-point rho on orbits of
the order-12 Frobenius automorphism (x -> x^p acts as x -> x^(p mod r) on GT),
then return Q = [a]Q0. A twisted-state walk avoids canonicalizing the state:
each step costs one dense Fp12 multiplication plus 8 Fp multiplications to
compute one key coordinate for all 12 conjugates. Claimed time_log2 = 126.32
dense-Fp12-multiplication equivalents; success >= 0.39 under the
random-mapping heuristic.

- [Algorithm, success and cost argument](proof.md)
- [Toy experiment (twisted-state orbit rho with an order-12 automorphism)](experiment.py)

Generic accounting of known techniques; not a practical break.

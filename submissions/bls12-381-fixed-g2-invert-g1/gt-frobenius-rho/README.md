# bls12-381-fixed-g2-invert-g1: Frobenius-class rho

This package reduces pairing inversion to a DLP in GT with base g=e(P0,Q). It solves that DLP with parallel Pollard rho over the order-12 Frobenius orbits of GT, then outputs [a]P0. The submitted time_log2 is 127.04. Success is at least 0.39 under one disclosed random-mapping heuristic.

- [Algorithm and cost argument](proof.md)
- [Claim](claim.json)
- [Toy class-walk check](experiment.py)

# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1, truncated-product canonicalisation

Solve T = [x]P0 by parallel distinguished-point Pollard rho on the (r-1)/6 orbits
of G1 under the order-6 automorphism group <-1, phi>, phi(x, y) = (beta x, y).
The class key is |y| (free), the canonical twist index is decided exactly from a
truncated fixed-point product (about 114 instructions, exact fallback on near-ties),
and a twisted multiplier table makes every step exactly one batched affine G1
addition (the reference operation) plus at most 240 instructions (the rare
fallback and small-batch inversions are metered separately). A global work meter
charges every action and halts at the claimed bound, so W <= 2^126.41 on every run. Fruitless
cycles are detected and escaped deterministically by doubling.

Claim: time_log2 = 126.41, success >= 0.40, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve run of the full algorithm](experiment.py)

Heuristic: random-mapping behaviour of the class walk (standard rho heuristic),
supported by toy runs. No production-size execution.

# Frobenius-class parallel Pollard rho in GT (fixed-G2 pairing inversion into G1)

Compute g = e(P0, Q) for the supplied Q, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output P = [x]P0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus ~2,400
instructions. Fruitless cycles are detected and escaped deterministically.

Claim: time_log2 = 125.89, success >= 0.40, memory < 2^73 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-field run of the full walk](experiment.py)

Heuristic: random-mapping behaviour of the class walk (standard rho heuristic),
supported by toy runs. No production-size execution.

# Frobenius-class parallel Pollard rho in GT, success target 0.39

Compute g = e(P0, Q) for the supplied Q, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output P = [x]P0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus at most 2,800
instructions. Fruitless cycles are detected and escaped deterministically.

Design provenance: the algorithm and cost ledger are those of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (Frobenius-class rho, 125.89). This package
keeps the algorithm and targets the track's success threshold 0.39 instead of
0.40, adds a derivation of the fruitless-cycle rates, defines the stuck-walk rule
and record format, and runs toy experiments with the production partition and
distinguished-point predicates on the raw class key.

Claim: time_log2 = 125.861, success >= 0.39, memory < 2^73 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run of the full walk, bound arithmetic](experiment.py)

Heuristic: random-mapping behaviour of the class walk (standard rho heuristic),
supported by toy runs. No production-size execution.

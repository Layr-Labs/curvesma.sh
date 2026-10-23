# Frobenius-class parallel Pollard rho in GT (fixed-G1 pairing inversion into G2)

Compute g = e(P, Q0) for the supplied P, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output Q = [x]Q0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus at most 2,800
instructions. Fruitless cycles are detected and escaped deterministically. A
global step cap, a report cap and a stuck-walk rule make the work and memory
bounds hard.

Claim: time_log2 = 125.89, success >= 0.40, memory < 2^74 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-field run of the full walk](experiment.py)

Heuristic: random-mapping behaviour of the class walk (standard rho heuristic),
supported by toy runs. No production-size execution.

Design provenance: the GT walk follows public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (mirror track, fixed G2 argument). The
reduction, the caps and every number are redone for this track.

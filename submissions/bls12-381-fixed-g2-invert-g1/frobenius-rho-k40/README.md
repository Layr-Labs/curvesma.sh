# Frobenius-class parallel Pollard rho in GT, success target 0.39

Compute g = e(P0, Q) for the supplied Q, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output P = [x]P0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus at most 2,250
instructions; the twelve class keys are two sign-symmetric triples, so their
minimum costs two subtractions instead of a twelve-way scan. Fruitless cycles are
detected and escaped deterministically. Distinguished
reports go into a directly addressed log of hard capacity 2^66 records (80 bytes
each), which bounds peak memory below 2^72.33 bytes; if it would overflow, the run
aborts and returns failure, a Chernoff-bounded loss (below 2^-1000). Matching is
an in-place heapsort of the log with a fixed worst-case instruction count per
record, and the report handling is charged at the cap (2^64.64 units).

Design provenance: the algorithm is that of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (Frobenius-class rho, 125.89). This package
keeps the algorithm and targets the track's success threshold 0.39 instead of
0.40, adds a derivation of the fruitless-cycle rates, defines the stuck-walk rule,
a hard 2^66-record directly addressed distinguished-report log with an
abort-on-overflow rule, worst-case sorted matching, and the record format, runs toy experiments with the production partition and
distinguished-point predicates on the raw class key, and itemizes the per-step
non-reference work instruction by instruction (2,026 counted, 2,250 budgeted).

Claim: time_log2 = 125.826, success >= 0.39, memory < 2^72.33 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run of the full walk, bound arithmetic](experiment.py)

Heuristic: random-mapping behaviour of the class walk (standard rho heuristic),
supported by toy runs. No production-size execution.

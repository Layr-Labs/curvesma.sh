# Frobenius-class parallel Pollard rho in GT, success target 0.39

Compute g = e(P0, Q) for the supplied Q, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output P = [x]P0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every step cost exactly one dense Fp12 multiplication plus at most 2,250
instructions; the twelve class keys are two sign-symmetric triples, so their
minimum costs two subtractions instead of a twelve-way scan. Fruitless cycles are
detected and escaped deterministically. Distinguished reports are appended to a
directly addressed log of exactly M R_w = 2^20 x 2^47 = 2^67 record slots (80 bytes
each): each walk slot keeps a cumulative report counter and is stopped when it
reaches R_w = 2^47, so the log capacity is a deterministic product with no
probability statement, and peak memory is below 2^73.4 bytes. Matching is an
in-place heapsort of the full 80-byte record (key T, then A, then B), so identical
records are adjacent and duplicate removal is one adjacent comparison per record;
each equal-key run of d distinct records afterwards has all its d(d-1)/2 pairs
charged against a hard pair budget P = 2^41. The report handling is charged at the
cap (2^66.15 units) plus the pair term (2^34.4 units).

Design provenance: the algorithm is that of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (Frobenius-class rho, 125.89). This package
keeps the algorithm and targets the track's success threshold 0.39 instead of
0.40, adds a derivation of the fruitless-cycle rates, defines the stuck-walk rule,
a directly addressed distinguished-report log bounded deterministically by
M R_w = 2^67 records (per-walk report cap, no report-count tail claim), worst-case
sorted matching, and the record format, runs toy experiments with the production
partition and distinguished-point predicates on the raw class key that also exercise
the report log (cycle checks before the append, the append rule, per-walk cap, cap
stops), matches the toy log offline with the production matcher, and itemizes the
per-step non-reference work instruction by instruction (2,026 counted, 2,250
budgeted).

The step order is the production order: the 2-cycle and checkpoint-cycle checks run
before the distinguished test, and a report equal to the walk's last appended record
(same 128-bit key and same canonical exponents) is not stored, so the exact duplicates
a fruitless two-cycle return would append (about 2^52 such events over the run) never
reach the log. A class fixes the canonical element, not the exponent pair: a walk that
revisits a distinguished class normally arrives with a different (A, B), and two such
records give the collision equation x = (A1 - A2)/(B2 - B1), which is the useful
recovery. The matcher's rule is single and complete: sort the full 80-byte record, drop
exact duplicates by an adjacent comparison, then charge every pair of the remaining
distinct records of an equal-key run against the hard pair budget P = 2^41, and solve the
first pair whose elements are equal with different exponent pairs. A genuine repeat is
charged as a pair like any other; no record is exempt. At most (2^67)^2/2 = 2^133 pairs
of distinct records give about 2^5 expected false-key pairs, so by Markov the pair count
exceeds 2^41 with probability at most 2^5/2^41 = 2^-36, which is the abort loss.

Claim: time_log2 = 125.826, success >= 0.39, memory < 2^73.4 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run of the full walk, bound arithmetic](experiment.py)

Heuristic: random-mapping behaviour of the class walk (standard rho heuristic),
supported by toy runs. No production-size execution.

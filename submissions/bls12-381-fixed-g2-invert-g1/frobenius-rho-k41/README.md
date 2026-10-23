# Mixed squaring / multiplication Frobenius-class Pollard rho in GT, success target 0.39

Compute g = e(P0, Q) for the supplied Q, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output P = [x]P0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets a multiplication step cost one dense Fp12 multiplication plus at most 2,800
instructions. Half of the steps are cyclotomic squarings instead: GT lies in the
cyclotomic subgroup, where Granger-Scott squaring costs about half a dense
multiplication, and y -> y^2 commutes with Frobenius, so it acts on classes and
merely doubles the exponent pair. The choice between the two step kinds is a
function of the class key alone, so walks that meet still merge. Fruitless cycles
are detected and escaped deterministically.

Design provenance: the reduction, the key formula, the twisted table, the class
walk, the cycle escapes and the record format are those of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (Frobenius-class rho, 125.89), as tightened
in the 125.861 package. This package adds the mixed step and its per-step ledger,
and targets the track's success threshold 0.39.

Claim: time_log2 = **125.767**, success >= 0.39, peak memory < 2^73 bytes,
preprocessing <= 2^31 units (included), nonuniform advice 0 bytes.

The bound is the mixed ledger at f = 1/2: 0.97058 dense-multiplication equivalents
per step, charged over the walk-step cap 1.10 x 1.000435 x S = 1.1335 S0 =
2^125.8093 steps, where 1.10 is the measured 10% step margin for the mixed map and
1.000435 is the fruitless-cycle factor. With E = 2^84.3 detection-delay steps and
the 2^62 report/output allowance the total is 2^125.76619, rounded up to 125.767.
The squaring price uses redc_mont_384 = 200; the analytic bracket 120..300 moves the
total between 2^125.740 and 2^125.798, and the claim holds for redc_mont_384 <= about
202.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks, toy-field run of the mixed walk, ledger](experiment.py)

Heuristics: random-mapping behaviour of the class walk (standard rho heuristic)
and of the mixed map (supported by the toy runs that compare it with a
multiplication-only control in the same harness). No production-size execution.

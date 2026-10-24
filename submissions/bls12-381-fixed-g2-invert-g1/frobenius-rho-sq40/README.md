# Frobenius-class parallel Pollard rho in GT with a cyclotomic-squaring branch, success target 0.39

Compute g = e(P0, Q) for the supplied Q, solve z = g^x in GT with a parallel
distinguished-point Pollard rho on the 12-element Frobenius classes of GT, and
output P = [x]P0. A four-multiplication key formula picks the class-canonical
conjugate without computing any full Frobenius, and a twisted multiplier table
lets every multiply step cost one dense Fp12 multiplication plus at most 2,250
instructions; the twelve class keys are two sign-symmetric triples, so their
minimum costs two subtractions instead of a twelve-way scan.

**The squaring branch.** The partition, the branch and the distinguished flag are
read from one free class
invariant: T0 = Re(c_0), the Fp part of the w^0 coefficient, which is equal on all
12 Frobenius conjugates and costs one load. Three disjoint bit fields of T0 give
the partition h (bits 0..9, 1024 values), the branch (bits 10..12: square iff the
value is below 3, probability 3/8) and the distinguished flag (bits 13..72, 60
bits). A squaring step is a Granger-Scott cyclotomic squaring y <- y^2 (valid on
GT, commutes with Frobenius, injective on classes because |GT| is odd) with
exponents (alpha, beta) <- (2 alpha, 2 beta) mod r; it computes no class key and
costs 11,700 instructions = 0.585 units instead of 1.1125. Squaring commutes with
Frobenius, so the branch, like the partition and the distinguished flag, is a map
on classes: outside cycle handling two walks that reach one class take the same step.
The two cycle escapes read stored history, but both are functions of the cycle reached
(the 2-cycle escape squares the smaller-T0 class of the two-cycle; the checkpoint escape
walks the cycle once to its smallest-T0 class and squares it), so two walks that reach
one cycle escape at the same class and their class sequences coincide from then on,
coalescing within one cycle length. That is enough for a class collision to give a usable
pair (proof.md section 3.1). The 2-cycle and
checkpoint comparisons therefore compare T0 (381 bits; a false match is about 2^-380.70
per comparison, below 2^-253 over the whole run, and can only cost a spurious escape,
never a false logarithm, because the report key comes from the class key K).

**Collision law.** Squaring-squaring pairs never merge, squaring-multiply pairs
merge at about 1/N and multiply-multiply pairs at about (1 - 1/(12R))/N, so
N' = N/(1 - s^2 - (1-s)^2/(12R)), S0' = sqrt(2 N' ln(1/0.61)) and S = 1.03 S0' =
2^125.780073 steps. Fruitless cycles are all-multiply only, because a cycle with
t >= 1 squaring steps multiplies the exponent pair by 2^t p^m and ord_r(2) =
(r-1)/16 = 2^250.857 forces t >= 2^247.27; the ledger nonetheless charges the
unbranched rate.

**A hard work meter.** Because a step costs 1.1125 or 0.585 units, the walk phase's
work is random, so every walk slot carries a private 128-bit meter in instructions
and is stopped at W_cap = 2^119.954254. All walk work is therefore at most
M x (W_cap - 1 + 22,250) = 2^125.666542 units as a deterministic bound; the abort's
effect on success is charged as a stated heuristic allowance of 0.01, which is 1.6
times its model value.

**Report machinery (direct addressing, deterministic capacity).** Reports go to a
log of exactly M R_w = 2^20 x 2^47 = 2^67 record slots of 80 bytes; each slot keeps
a cumulative report counter that no restart resets and stops at R_w, so the log can
never overflow and no report-count probability enters memory or processing. A
report equal to the walk's last record is not stored, and both cycle checks run
before the distinguished test. Matching is an in-place heapsort of the full 80-byte
record (key T, then A, then B), so identical records are adjacent and duplicate
removal is one adjacent comparison per record; each equal-key run of d distinct
records has all its d(d-1)/2 pairs charged against a hard pair budget P = 2^41. The
report handling is charged at the cap (2^66.39 units), the pair term (2^34.36), the
restart starts (2^71.12) and setup/output (2^31); preprocessing and table construction
are at most 2^57 units, including the 2^56.30-unit log initialization, and all of it is
inside W.

Claim: time_log2 = 125.667, success >= 0.39 (0.40808 model, rounded down, minus a 0.01
meter allowance minus 2^-18 of loss terms >= 0.3980), memory < 2^73.4 bytes.

- [Algorithm, success argument, collision law and cost](proof.md)
- [Experiment: BLS12-381 identity and T0 checks, the production step rule on real GT elements, the measured collision law, toy-field run of the full mixed walk, bound arithmetic](experiment.py)

Heuristics: random-mapping behaviour of the mixed class walk (standard rho
heuristic, with the branch-adjusted collision law measured directly in part E),
uniformity of T0's decision bits (measured on real GT elements in part A2), and the
work-meter allowance. No production-size execution; toy evidence is evidence at its
own scale.

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
the partition h (bits 0..9, 1024 values), the branch (bits 10..15: square iff the
value is below 35, probability 35/64) and the distinguished flag (bits 16..75, 60
bits). A squaring step is a cyclotomic squaring y <- y^2 by blst's
`blst_fp12_cyclotomic_sqr` (valid on GT, commutes with Frobenius, injective on
classes because |GT| is odd) with
exponents (alpha, beta) <- (2 alpha, 2 beta) mod r; it computes no class key and
costs 10,098 instructions = 0.5049 units instead of 1.1125 (9,798 of them measured
on the track's own COST_MODEL toolchain, identical at five slice offsets, the
routine being straight-line with no branch, plus 300 for the doubling, the shared
rows and slack; the same harness reproduces blst_fp_mul = 383.0 and
blst_fp12_mul = 22,996.0 exactly). With a squaring step this cheap the
W-minimising branch fraction is 35/64, the best value a 2^-k bit field can state
near the continuous optimum 0.5461: the old s = 3/8 would now cost 0.034 bits. Squaring commutes with
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
N' = N/(1 - s^2 - (1-s)^2/(12R)) = 2^251.785202, S0' = sqrt(2 N' ln(1/0.61)) =
2^125.884326 and S = 1.03 S0' = 2^125.926969 steps. Fruitless cycles are
all-multiply only, because a cycle with
t >= 1 squaring steps multiplies the exponent pair by 2^t p^m and ord_r(2) =
(r-1)/16 = 2^250.857 forces t >= 2^247.27; the ledger nonetheless charges the
unbranched rate.

**A hard work meter, with a concentration bound.** Because a step costs 1.1125 or
0.5049 units, the walk phase's work is random, so every walk slot carries a private
128-bit meter in instructions and is stopped at W_cap = 2^119.871882. All walk work is
therefore at most M x (W_cap - 1 + 22,250) = 2^125.584170 units as a deterministic
bound. Since step costs lie in [0.50490, 1.1125] and are independent of the past,
Hoeffding bounds the chance that a slot's work exceeds 1.01 of its nominal work by
exp(-2 (0.01 c)^2 n / (c_mul - c_sq)^2) = 2^-94.89 for a slot (c = 0.78021875, n =
2^105.93) and below 2^-100 over the M slots, so the cap retains all S productive steps
except with that probability. There is no "S/1.01 steps" equivalence; the separate 0.01
success allowance is slack for the heuristic parts of H1.

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
restart starts (2^71.24) and setup/output (2^31); preprocessing and table construction
are at most 2^57 units, including the 2^56.30-unit log initialization, and all of it is
inside W.

Claim: time_log2 = 125.585, success >= 0.39 (0.40808 model, rounded down, minus a 0.01
meter allowance minus 2^-18 of loss terms >= 0.3980), memory < 2^73.4 bytes.

- [Algorithm, success argument, collision law and cost](proof.md)
- [Experiment: BLS12-381 identity and T0 checks, the production step rule on real GT elements, the measured collision law, toy-field run of the full mixed walk, bound arithmetic](experiment.py)

Heuristics: random-mapping behaviour of the mixed class walk (standard rho
heuristic, with the branch-adjusted collision law measured directly in part E),
uniformity of T0's decision bits (measured on real GT elements in part A2), and the
work meter, whose abort probability is bounded by Hoeffding rather than assumed. No production-size execution; toy evidence is evidence at its
own scale.

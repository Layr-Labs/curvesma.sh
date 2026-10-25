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
12 Frobenius conjugates and costs one load. The fields are read from its Montgomery
value T0*2^384 mod p, exactly the representation the measured kernels hold, so no
canonical-residue conversion is claimed or charged. The class keys are compared in that
same form (key*R mod p, which the four kernel multiplications produce directly), so kappa
can differ from the canonical-residue argmin; the report key is bits 70..197 of that
Montgomery-ordered K. Since K is a minimum of twelve keys it is not uniform, and the package
claims no uniform-bits preservation: proof.md section 6 takes those 128 bits as
pseudo-random under H1 and part A2 measures them. Three
disjoint bit fields of T0 give
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
near the continuous optimum 0.5464: the old s = 3/8 would now cost 0.034 bits. Squaring commutes with
Frobenius, so the branch, like the partition and the distinguished flag, is a map
on classes: outside cycle handling two walks that reach one class take the same step.
The two cycle escapes read stored history, but both are functions of the cycle reached
(the 2-cycle escape squares the smaller-T0 class of the two-cycle; the checkpoint escape
walks the cycle once to its smallest-T0 class and squares it), so two walks that reach
one cycle escape at the same class and follow the same class sequence from their first
common class. Their detection times can differ, because checkpoints are installed on the
retained step counter mod 64 and no escape resets it: the offset is at most 64 + L <= 128
steps, and it is charged (proof.md sections 3.1 and 4.1) rather than called a one-cycle
coalescence. A class collision still gives a usable pair on its own: two records of one
class with different exponent pairs solve for x, and matching is offline. The 2-cycle and
checkpoint comparisons therefore compare T0 (381 bits; a false match is about 2^-380.70
per comparison, below 2^-253 over the whole run, and can only cost a spurious escape,
never a false logarithm, because the report key comes from the class key K).

**Collision law.** Squaring-squaring pairs never merge, squaring-multiply pairs
merge at about 1/N and multiply-multiply pairs at about (1 - 1/(12R))/N, so
N' = N/(1 - s^2 - (1-s)^2/(12R)) = 2^251.785203, S0' = sqrt(2 N' ln(1/0.61)) =
2^125.884326 and S = 1.05 S0' = 2^125.954715 steps. Fruitless cycles are
all-multiply only, because a cycle with
t >= 1 squaring steps multiplies the exponent pair by 2^t p^m and ord_r(2) =
(r-1)/16 = 2^250.857 forces t >= 2^247.27; the ledger nonetheless charges the
unbranched rate.

**A hard work meter, with a dependent-process account.** Because a step costs 1.1125 or
0.5049 units and the branch bit is a function of the current class, revisits repeat an
earlier decision and step costs are not i.i.d. across time; so the meter's account is split
by step type. Every walk slot carries a private 128-bit meter in instructions and is
stopped at W_cap = 2^119.899627, which is 1.01 of the nominal productive work plus a
2^80-step per-slot allowance for ordinary revisits (it raises log2 W by 3.2e-8 bits). All
walk work is therefore at most M x (W_cap - 1 + 25,250) = 2^125.611915 units as a
deterministic bound, whatever the walk does. On newly exposed classes H1 supplies a fresh
branch draw, so the productive-step costs are i.i.d. in [0.50490, 1.1125]. The 1% band is
split once between productive and cycle work. A slot aborts only if (i) its productive work
exceeds (1 + 10^-5) of its productive nominal, which Hoeffding bounds under H1 (c =
0.78021875, n = 2^105.96) by 2^(-2^74.98) per slot and 2^(-2^94.98) over the run (74.98 and
94.98 are exponents of the exponent, not probabilities), or (ii) its cycle work exceeds its
allowance plus the cycle band: one meter sums productive, cycle and report work and the slot
stops at W_cap, so a slot with ordinary productive work can still stop early. The cycle band
is 0.01 nominal + 2^80 steps less 10^-5 of the productive nominal (0.00999 of nominal) = 2^98.80
steps per slot (experiment part C, `cycle_abort`), while H1's expected cycle work is 2^91.73
steps per slot, so P(cycle work > band) <= 2^91.73/2^98.80 = 2^-7.07 (0.0074251) by Markov.
Cycle handling is escaped work charged at 8.35 times its H1 mean rate as a declared
allowance, not a tail bound: no concentration inequality is claimed for the history-dependent
count of cycle events, and what bounds those steps is the hard meter, which charges every one
of them. A slot lost to (ii) is charged as if it stored no DP (DPs stored before the stop stay in
the log), so 0.00745 of the walk population is charged
inside the same 0.01 success allowance, which leaves 0.00254 for the other named losses. A slot can also meet its own class: under H1 with probability at
most n^2/(2N') = 2^-40.875 per slot it then enters a cycle the 64-step checkpoint does not
escape and spends the rest of its budget on repeats instead of fresh classes. That slot
wastes at most its remaining n steps and its meter stops it at W_cap, so the work bound is
unchanged; the cost is on success, at most 2^-20.96 of S lost except with probability below
2^-19.91 (Markov), a success loss of at most 2^-21.62, charged in the loss list beside the
2^-19 self-collision term. There is no "S/1.01 steps" equivalence, and the separate 0.01
success allowance is slack for the heuristic parts of H1. The fresh-class law and the cycle
count remain heuristic and are stated as such in proof.md section 4.

**Report machinery (direct addressing, deterministic capacity).** Reports go to a
log of exactly M R_w = 2^20 x 2^47 = 2^67 record slots of 80 bytes. Each slot keeps a
cumulative counter of distinguished-point ATTEMPTS, not of stored records: every
distinguished hit increments it, stored or suppressed as a repeat of the walk's last
record, the slot stops at R_w = 2^47, and each attempt is charged 3,000 instructions to
the slot's work meter as well as counted in the 2^67 cap-unit report group. So the log
cannot overflow (stored records are at most attempts), no report-count probability enters
memory or processing, and a repeated report inside an unescaped fruitless cycle pays per
pass instead of hiding behind duplicate suppression, which is a storage rule only. Both
cycle checks run before the distinguished test. Matching is an in-place heapsort of the full 80-byte
record (key T, then A, then B), so identical records are adjacent and duplicate
removal is one adjacent comparison per record; each equal-key run of d distinct
records has all its d(d-1)/2 pairs charged against a hard pair budget P = 2^41. The
report handling is charged at the cap (2^66.39 units), the pair term (2^34.36), the
restart starts (2^71.28) and setup/output (2^31); preprocessing and table construction
are at most 2^57 units, including the 2^56.30-unit log initialization, and all of it is
inside W.

Claim: time_log2 = 125.612, success >= 0.39 (0.42013 model, rounded down, minus a 0.01
meter allowance minus 2^-18 of loss terms >= 0.41012, with 0.02012 of margin; the toys pool
0.4080 +- 0.0090 at this budget, 0.971 of the model, which is 0.3980 after the same allowance),
memory < 2^73.4 bytes.

- [Algorithm, success argument, collision law and cost](proof.md)
- [Experiments in detail: what each part checks and the figures it prints](experiments.md)
- [Experiment: BLS12-381 identity and bit-frequency checks, the production step rule on real GT elements and both cycle escapes, the toy-field walk and collision law, the bound arithmetic](experiment.py)
- [Measurement record: the squaring price, with the harness, the raw Callgrind records and the disassembly audit](measure/README.md)
- [Extended experiment results](measure/experiment-results.md)

Heuristics: random-mapping behaviour of the mixed class walk (standard rho
heuristic, with the branch-adjusted collision law measured directly in part E),
uniformity of T0's decision bits (measured on real GT elements in part A2), the
fresh-class branch law that lets Hoeffding apply to the productive steps, which is
part of the same random-mapping heuristic, and the cycle-event count, whose charge
carries an 8.35x safety factor. No production-size execution; toy evidence is evidence
at its own scale.

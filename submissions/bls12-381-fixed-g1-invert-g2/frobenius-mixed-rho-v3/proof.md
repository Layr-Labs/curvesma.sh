# frobenius-mixed-rho-v3: Frobenius-class mixed add/square rho on GT (revision v5)

Claim: for the uniform challenge of `rules.md` (P uniform nonidentity in G1, z uniform in GT, fresh algorithm coins), the procedure below
outputs Q in G2 with e(P, Q) = z with probability at least 0.41 (model 0.410000000001; lower bound
0.410000000001 after the failure terms of Section 2), while its total work, counting every operation it performs, is
at most W_MAX = 2^125.789684 target-group-operation equivalents (score 31447421/250000), with probability at least 1 - 2^-60
that the first M main steps are executed before the cutoff (Section 3). Peak memory is at most 2^73.166564 bytes; preprocessing is
inside W_MAX (2^76.418283 units); nonuniform advice 0 bytes. The bound rests on the heuristics H5 and H7 (Section 6) and on the
derived bounds. It is a classical generic attack, not a practical break.

**Revision v5.** Round 4 rejected v4 on three points: (1) the cancellation classification was wrong (1 + lambda^4 + lambda^8 = 0
is a three-term cancellation, and odd lengths exist), and the non-formal bound was applied to identically zero relations;
(2) H6 (that 2^q is not a 12th root of unity for square counts up to 2^66) was supported only by a test up to 2^(12 * 200000); (3) the cutoff did not show that the
first M main steps execute, and the distinguished-point storage had no capacity rule. v5 answers each:
(1) the formal cancellations are classified over Z[zeta_12]: a closed window with zero exponent coefficient is an additive window
whose multiset of (adder, rotation) pairs is a sum of 2-gons {(h,k),(h,k+6)} and 3-gons {(h,k),(h,k+4),(h,k+8)}, by the
Lam-Leung theorem (checked by enumeration up to total count 7). Each length 65..400 is bounded by an explicit union bound
(Section 2), and longer windows by the probability that every step is additive. The non-formal closures remain a heuristic (H5);
(2) the order of 2 modulo r is computed exactly from the complete factorisation of r - 1 (Section 2): 2^q is not a 12th root of
unity for q < 2^247, so this premise is now a computation rather than an assumption;
(3) W_MAX is charged by bounds that hold with stated probability (Section 3): Chernoff bounds for the square count and the cycle
and distinguished-point counts, a Markov bound for fingerprint matches, and a distinguished-point table with a stopping rule
when it is full (Section 4).
The score moves from 125.724950 (v5) to 125.789684 (+0.064734 bits). The increase comes from the per-type
charges: the adder step is charged its measured instruction count, 1 + 3,754 / 20,000 = 1.18771 units (3,734 other instructions and 20 for the anchor refresh), instead of the averaged
step charge, and the distinguished point is charged its measured store plus the sort bound (0.235 units) instead of 2^20 instructions.
The failure terms are unchanged (below 2e-18).

**Revision v6.** Round 5 asked that the per-step charges be worst-case, not averages, and that the distinguished-point charge be
measured or derived. v6 writes the step code without data-dependent branches (comparisons over all limbs, a branch-free modular
addition, a branch-free minimum), so each step type has one instruction count, and measures it on forced inputs (all-adder,
all-square, and the all-adder run with a distinguished point at every step). The adder step is the reference operation
(blst_fp12_mul, charged 1 unit), so only the other instructions are divided by C. The distinguished-point store is measured
(1,479 instructions) and the sort-and-lookup bound is 67 levels of a measured 48-instruction comparison and exchange. The
2^20-instruction assumption is removed. The variable parts left are the detected cycles (bounded as 66 adder steps each) and the
restarts and fingerprint verifications (each charged as an exponentiation, 381 reference operations).

## 1. Procedure

Reduction. g = e(P, Q0) generates GT, where Q0 is the fixed G2 generator; z = g^y gives Q = [y]Q0 with e(P, Q) = z. Input a nonidentity P in G1 and a target z in GT. Return the identity if z = 1 (Q = O is valid).
Otherwise run the walk until a verified y is found or the work counter reaches W_MAX, then output Q = [y]Q0, verified by the
pairing equation e(P, Q) = z (one pairing, charged). Exceptional inputs (P = O) are rejected by the rules.

Walk (`walk-spec.md`, the specification used by `cscore/gt.py` and the blst kernel `measurement/step.c`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a.
  If h < S, square x and double (a, b). Otherwise multiply x by T[h][k] = sigma^(-k)(R_h), where R_h = g^(c_h) z^(d_h), and update
  (a, b) by lambda^(-k) (c_h, d_h). The Frobenius map is never applied per step.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle;
  a new key equal to the walk's anchor key (recorded every 64 of its own steps) inside a window of 64 steps is a cycle of length at
  most 64. Either escapes by squaring the least-key state of the cycle. Escapes are charged.
- Trail cap and restart: the walk's trail age counts main steps since its last distinguished point. At age 2^66 the walk restarts
  from a fresh random start (charged) and its age resets.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of
  the Montgomery integer; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index.
- Solving: records with equal fingerprints are compared pairwise. A pair from one class with B1 != B2 gives
  y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is discarded at the cost of one comparison; no trail is replayed. Every
  candidate is verified by g^y = z before output.
- Work counter and cutoff: each operation adds its charge (Section 3). The walks run while W < W_MAX - reserve, with reserve the
  solve and output charge. The distinguished-point table stops the procedure when full (Section 4).

Consistency checks: the walk invariant x = g^(a + b y) holds on every state (replica at production parameters, 512 and 4096 rounds,
64 walks); the canonical exponents of the stored distinguished points satisfy sigma^kn(x) = g^(lambda^kn (a + b y)) (replica, bit for
bit on the stored records); the twisted table equals sigma^(-k)(g^c z^d) (test); the solver recovers y end to end (tests, and
100 toy trials, all verified).

## 2. Success probability, failure terms and the classification

Under H1 the walk needs T = 2^125.836359 fresh steps for tau = T^2 V / (2N), with V = 0.799413681. The collisions the
distinguished-point test can see are a fraction f = 1 - 2.1e-26 of all collisions (cap loss). Main steps
M = 2^125.844635 include the detected-cycle waste w = 0.005720 (127 steps per detection, at most
4.5040e-05 detections per step) and a detection lag of 2^80 steps. Then

  P = 1 - exp(-f tau_eff) - delta, tau_eff = T^2 V / (2N) = 0.527632742, delta = 2.16e-18.

P = 0.410000000001 before delta and 0.410000000001 after it. If the wasted steps were not budgeted, P would fall to
0.406438.

Classification. A closed window (the walk's state returns to its earlier class or state) is one of two kinds. (a) If the
window contains a square, its exponent coefficient is 2^q lambda^J - 1 with q at most the window length; this is nonzero because the
order of 2 modulo r is computed (below), so the closure is an event of probability 1/r over the coins (H5). (b) If the window is all
additive, the coefficient of each adder h is sum_k n_{h,k} lambda^(-k) in F_r. For windows of length below 2^63 this is zero in F_r
only if the vector (n_{h,k})_k is zero in Z[zeta_12] (a nonzero element of the prime ideal of Z[zeta_12] above r has norm at
least r, while its norm is at most (total count)^4 < r). The nonnegative vanishing sums of 12th roots of unity are sums of
2-gons and 3-gons (Lam-Leung; enumerated up to total count 7 here), so these closures are the formal ones counted below.
Order of 2: the factorisation r - 1 = 2^32 * 3 * 11 * 19 * 10177 times nine primes of 17 to 28 bits (verified prime, product checked)
gives ord_r(2) = 3277242198445386904965483781761622864855659531282977363912728668746161324032 (251 bits). The smallest q >= 1 with 2^q a
12th root of unity is q = ord/gcd(ord, 12), a 248-bit number larger than 2^247, so 2^q is not a 12th root of unity for any q below 2^247.

Coverage of window lengths: formal windows of every length are covered. The union bound over gon decompositions (blocks of size 2 and
3, so odd lengths included) is taken for lengths 65 to 400, and the additive tail covers every length above 400 (up to the cap and
beyond). The non-formal term covers every window up to the cap. A trail cannot repeat a state across a window longer than the cap,
because the walk restarts at the cap first, so no longer closure can trap a walk.

Failure terms (each a probability over the whole run; computed in `accounting.py`, function `trap_bound` and `build`):
| term | value | source |
| --- | ---: | --- |
| formal windows, lengths 65 to 400 (per step) | 1.77e-96 | union bound over gon decompositions (Lam-Leung), at most 2^-318 for any single length |
| all-additive windows longer than 400 (per step) | 1.25e-103 | (1 - s)^401 / s |
| non-formal closures (per step) | 1.69e-56 | 12 CAP / r candidates, each closing with probability 1/r (heuristic H5) |
| trap term over the run | 1.29e-18 | M times the per-step total; the non-formal part dominates |
| fingerprint Markov (metered matches) | 8.67e-19 | 2^-60 |
| square count (Hoeffding) | 0.0e+00 | exp(-2 (1e-6)^2 M) |
| detected cycles of length <= 64 (blocked Chernoff) | 2^-2.59e+25 | 64 residue classes, independent within a class |
| non-formal detections of length <= 64 | 1.1e-36 | expectation 12 * 64 * M / r (Markov) |
| distinguished-point table overflow | 2^-3.19e+19 | exp(-M 2^-60 / 3) |
| total failure kept | 2.16e-18 | sum of the kept terms |

Per-step formal bound: the union bound for lengths 65..400 is at most 2^-318 for any single length, and the tail for windows longer
than 400 is at most (1 - s)^401 / s; the 2-cycle rate of the formal pair family is (1 - s) / (12 * 1024) = 4.498164e-05 per step against
4.507601e-05 measured (ratio 1.0021).

Heuristics (disclosed):
- **H1 (random-mapping collision model).** The walk on the N = (r-1)/12 classes behaves as a random map whose collision rate is
  V = 0.799413681 (the square branch has probability s, and each of the 566 reachable additive adders has probability 1/r_a); tau = T^2 V /
  (2N) is treated as Exp(1). Support: RESEARCH.md 3.4 and the toy sweeps. Limits: toy groups have N up to 2^17 and 32 adders.
- **H5 (non-formal closures).** A closure whose coefficient is nonzero is an event of probability 1/r over the coins. This treats the
  walk's path as independent of the coins; it is the one heuristic the trap term depends on. It contributes 1.69e-56
  per step (about 1.3e-18 over the run).
- **H7 (uniform independent choices).** The adder choice (h, k), the square decision, the distinguished-point test and the fingerprint
  bits are uniform and independent across steps. The union bound, the square count bound, the distinguished-point count and the
  detection rate use H7. The measured 2-cycle rate supports H7 for the adder pair (ratio 1.0021).
- **H3 (per-type step cost).** The adder and square steps have one instruction count each on the straight-line path (measured on
  forced inputs, agreeing to 0.03% between seeds, with the difference explained by cycle events). The variable parts are bounded
  separately (events, distinguished points). The bound is a measured worst case for this code, not a proof for every machine state.
- **H2 (short cycles).** 2-cycles are checked every step and cycles of length 3 to 64 against a per-walk anchor; the 2-cycle branch
  is checked for exact replica equality; the 3-to-64 branch is mirrored in the replica and fired once in production.

## 3. Computation and the work cutoff

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), measured 22,998 instructions (1.1499 units; C = 20,000 is the
organizer's rounding), charged 1 as the rules require. Other work is measured in instructions and divided by C.

Measurement (`measurement/`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0,
`-O3 -fno-tree-vectorize -mgeneral-regs-only`), the pinned blst commit f7a88067. Callgrind Ir is differenced between 256 and 512
rounds of 64 walks at S = 458, with the distinguished-point density set to zero (DP work is charged separately).

| forced step type (S = 458 or forced; zero DP density) | per-step instructions | notes |
| --- | ---: | --- |
| all-adder steps, seed 1 | 26,730.16 | includes 4 detected cycles (conservative) |
| all-adder steps, seed 2 | 26,722.87 | includes 1 detected cycle |
| all-square steps, seed 1 | 13,511.44 | no detected cycles |
| all-adder steps with a distinguished point at every step | 28,209.16 | the store costs 1,479.0 instructions per point |
| sort: one 128-bit comparison and one 80-byte exchange | 48.0 per iteration | from 196,632 at 4096 and 393,240 at 8192 iterations |
| reference operation (blst_fp12_mul) | 22,998 (organizer: 22,996) | charged 1 unit |
| cyclotomic square (blst_fp12_cyclotomic_sqr) | 9,804 | charged 9,804 / C = 0.4902 |

Per-type charges (the step code is straight-line, with the comparisons and selections written without data-dependent branches,
so each step type has one instruction count; the two all-adder seeds agree to 0.03%, and the all-square run has no cycle events, which is evidence for the constancy of the path, not a proof for every machine state): an adder step is charged
1.18771 units (the reference operation, 1, plus the other instructions divided by C = 20,000), and a square step is charged
0.67657 units (the cyclotomic square at 9,804 instructions plus the other instructions). Both include the anchor refresh charged
at every step (20 instructions).

Charge for the first M main steps (each bound holds with the stated probability):
- Operation mix: the square count is at least (s - 10^-6) M with probability at least 1 - exp(-2 * 10^-12 * M) (Hoeffding, H7),
  so the per-step charge is (1 - s_low) A + s_low Q with s_low = 0.447264625, i.e. 0.959096 units per main step.
- Detected cycles: at most 4.5040e-05 per step (2-cycle term (1 - s)/(12 r_a), formal windows of length 3 to 64, and the
  non-formal term), with 1.001 slack; each costs at most 66 adder-type steps (the walk-through of at most 64 steps, the check and the
  escape), charged 66 A: 0.003531 units per step in total.
  Concentration (not an expectation): under H7 the formal detections of length at most 64 are functions of windows of at most 64
  consecutive i.i.d. choices. Split the steps into 64 residue classes mod 64; within a class, windows do not overlap, so the
  indicators are independent, with mean at most p M / 64 per class. Chernoff on each class and a union over the 64 classes give
  failure at most 64 exp(-(10^-3)^2 (p M / 64) / 3) (about 2^-(2.5 * 10^25) at the production M). The non-formal detections of
  length at most 64 have expected count 12 * 64 M / r (about 1e-36), and enter the failure total by Markov.
- Distinguished points: the table holds at most 2 M 2^-60 + 2 * 2^20 records; the count of distinguished points is a sum of i.i.d.
  tests with mean M 2^-60, so the count exceeds the capacity with failure at most exp(-M 2^-60 / 3) (Chernoff). The procedure stops
  when the table is full (overflow rule, Section 4), so this is the only way the capacity can bind. Each point costs 2^20 instructions.
- Fingerprint matches are metered: every match is charged when it occurs (one exponentiation), and the metered total is compared with
  the charge of 2^60 times the expected count, 4.69 log2 pairs in expectation. Markov gives failure 2^-60 for
  the metered total exceeding that charge; if it does, the procedure stops with failure.
- Trail-cap restarts: at most M / 2^66 + 2^20, covered by the per-distinguished-point restart charge.

Components of W_MAX (log2 units; the score is the rounded-up cutoff 125.789684):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 125.784383 | 125.844635 log2 main steps at 0.959096 units each (charged square count, adders, c2) |
| detected cycles | 117.698779 | at most 4.5040e-05 events per step, 66 extra operations each |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts and restarts | 76.418283 | 2^20 starts and one restart per distinguished point (covers every trail-cap restart) (preprocessing) |
| distinguished points | 64.756902 | record, sort and lookup, 2^20 instructions each, 2^66.844635 records (capacity) |
| fingerprint verification | 73.262918 | Markov charge at 2^60 times the mean of about 2^4.69 false pairs, one exponentiation each |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

Probability that the first M main steps are charged within W_MAX: at least 1 - 2^-60 under H7. The fingerprint Markov bound is the
only term of size above 2^-1000 in the charge bounds above. Consequently the procedure executes M main steps before the cutoff, and the success bound of Section 2 applies.

## 4. Memory and the distinguished-point table

Peak storage: distinguished-point records, capacity D = 2 M 2^-60 + 2^21 = 2^66.8446 records of 80 bytes, at most 2^73.166564 bytes
with the adder table (12 x 1024 elements, about 7.8 MiB), the walk states (2^20 x 640 bytes), and per walk the detector state (98 + 16
bytes). The table has an overflow rule: when D records are stored, the procedure stops and reports failure. The overflow
probability is at most exp(-M 2^-60 / 3), far below 2^-(10^19). The fingerprint is key bits 70..197 of the Montgomery integer; the
expected number of equal-fingerprint pairs is 2^4.69, each verified by one exponentiation (charged above).

## 5. Evidence: toy-scale sweeps

Group: the order-ell subgroup of F_(q^12)^* with ell dividing q^4 - q^2 + 1; `cscore/gt.py`, bundled into `experiment.py`. The toy taus use
the partition-correct V. The toy walk detects full cycles with its own rules and uses dense distinguished points, so these runs test
H1's shape and the solver, not H2, H5 or the production detector.

| run | ell | N | s (effective) | toy V | trials | reached a useful collision | mean tau (model 1) | success at tau = 0.5277 [Wilson 95%] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| initial input (bits 14, seed 1, 100 trials) | 28393 | 2366 | 0.4688 | 0.763672 | 100 | 100 | 0.8625 | 0.4800 [0.3846, 0.5768] |
| add-only, bits 14 | 28393 | 2366 | 0.0000 | 0.968750 | 3000 | 3000 | 0.9634 | 0.4170 [0.3995, 0.4347] |
| mixed, bits 14 | 28393 | 2366 | 0.4688 | 0.763672 | 3000 | 2997 | 0.9592 | 0.4277 [0.4101, 0.4455] |
| add-only, bits 18 | 922561 | 76880 | 0.0000 | 0.968750 | 800 | 800 | 0.9514 | 0.4263 [0.3924, 0.4608] |
| mixed, bits 18 | 922561 | 76880 | 0.4688 | 0.763672 | 800 | 800 | 0.9987 | 0.4088 [0.3752, 0.4432] |

Pooled mixed runs at tau = 0.5277 (N = 2366 and N = 76880): 1610 of 3800 succeed, empirical 0.4237, Wilson 95% interval
[0.4081, 0.4395].

## 6. Limitations

- H1 is extrapolated from toy subgroups with N up to 2^17 and 32 adders to N about 2^251 and 1024 adders.
- H5 is a heuristic: the non-formal closure term depends on the coins being independent of the path. H7 is the uniform-choice model
  used for the mix, the detection rate and the capacity. These are the stated premises of the bounds; they are not proofs.
- The trail-cap restart and the 3-to-64 cycle branch are specified or mirrored but not exercised by an exact-equality test
  (measurement runs are 2^26 steps per seed; the 3-to-64 branch fired once).
- The per-type step counts are measured on forced inputs and the code path has no data-dependent branch on the step path; the
  remaining variable parts (the cycle events and the distinguished-point store) are bounded separately (66 adder steps per event;
  1,479 + 67 * 48 + 10 instructions per point, with the sort bound from the measured exchange). The solve, verify and output charge of
  2^29 instructions is an assumed bound; it is 2^14.7 of the 2^125.8 units and does not affect the score at the digit shown.
- Memory of about 2^73.166564 bytes is the capacity of the distinguished-point table, which the procedure enforces; it is not
  realisable on the hardware used here and is reported, not scored.
- The reference operation is priced by the organizer rule at 1 unit; the measured instruction count is 1.1499 units.

## References

- Lam and Leung, On vanishing sums of roots of unity (2000), for the classification of nonnegative vanishing sums.
- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image, py_ecc (gen_inputs.py)
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind at zero DP density, two sizes, mul and sqr loops), analyses the
interval counts, and runs the replica on each input set. `native_run.sh` runs the uninstrumented kernel for the event counts.

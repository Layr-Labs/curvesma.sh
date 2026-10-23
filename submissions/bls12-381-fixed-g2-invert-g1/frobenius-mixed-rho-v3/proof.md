# frobenius-mixed-rho-v3: Frobenius-class mixed add/square rho on GT (revision v7)

Claim: for the uniform challenge of `rules.md` (Q uniform nonidentity in G2, z uniform in GT, fresh algorithm coins), the procedure below
outputs P in G1 with e(P, Q) = z with probability at least 0.41 (model 0.410000000001; lower bound 0.410000000001 after the failure terms of
Section 2), while its total work, counting every operation it performs, is at most W_MAX = 2^126.097138 target-group-operation
equivalents (score 63048569/500000). The first M = 2^125.844673 main steps are executed before the cutoff except with probability
2.16e-18 (Section 3). Peak memory is at most 2^73.321929 bytes, including the padding and scratch of the distinguished-point
table (Section 4); preprocessing is inside W_MAX (2^69.418320 units); nonuniform advice is 0 bytes. The bound rests on the heuristics
H1, H5 and H7 (Section 6) and on the derived bounds. It is a classical generic attack, not a practical break.

**Revision v7.** Round 7 rejected v6 on four points. Each is answered here.

(1) Number theory. The factorisation of r - 1 is complete: 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403
* 52437899 * 254760293^2. `cscore/number_theory.py` (vendored into `accounting.py`) checks the product, the primality of every factor
(Miller-Rabin with the first twelve primes as bases, deterministic below 3.18 * 10^23), the exact order of 2 modulo r by a certificate
(2^o = 1 and 2^(o/p) != 1 for each prime p dividing o), and the smallest q with 2^q a 12th root of unity (q0 = o / gcd(o, 12), checked
by 2^(12 q0) = 1 and gcd(q0, 12/gcd(o,12)) = 1). `accounting.py` asserts these values before it uses them. The output is in
`evidence/number_theory_check.out`.

(2) Classification. Section 2 proves Proposition G (nonnegative vanishing sums of 12th roots of unity are sums of 2-gons and 3-gons)
from the minimal polynomial x^4 - x^2 + 1 and four linear conditions. No enumeration is used.

(3) Charges. Every main step is charged the adder cost (1.18771 units), the larger of the two step charges, so no square-count
concentration enters the charge. The square count enters only the success model, through a Hoeffding deviation over fresh steps plus
the revisit steps, which repeat an earlier decision and are at most the detections. Detected cycles are counted by a conditional
domination argument over 64 residue classes, which needs only the conditional probability bound given the history 64 steps earlier
(no independence claim). Distinguished points are charged at their capacity, a deterministic quantity. The capacity is bounded from the
fresh flags (uniform under H7) with an explicit compensator for revisits (Section 4).

(4) Distinguished-point procedure. The table is a bitonic sorting network over 2^67.0000 slots (log2), data-oblivious, in place, with
one 80-byte temporary. It sorts by the 128-bit fingerprint and then scans equal neighbours. Measured: 108.003 instructions per
compare-exchange (the same at 8192 slots, 108.001) and 13.0 instructions per slot for the scan; 1,479 per store. The earlier text
of 2^20 instructions per point is removed from this proof, from walk-spec.md and from claim.json.

Score: v6 125.789684 to v7 126.097138 (+0.307454 bits). The increase comes from charging every main step the adder cost (1.18771 units)
instead of the averaged step cost, about +0.31 bits. The distinguished-point table costs 6.452 units per record over 2^66.9480 records,
which is 2^-56.46 of W_MAX. The failure total is 2.16e-18.

## 1. Procedure

Reduction. g = e(P0, Q) generates GT, where P0 is the fixed G1 generator; z = g^y gives P = [y]P0 with e(P, Q) = z. Input a nonidentity Q in G2 and a target z in GT. Return the identity if z = 1 (P = O is valid).
Otherwise run the walk until a verified y is found or the work counter reaches W_MAX, then output P = [y]P0, verified by the pairing
equation e(P, Q) = z (one pairing, charged). Exceptional inputs (Q = O) are rejected by the rules.

Walk (`walk-spec.md`, the specification used by `cscore/gt.py` and the blst kernel `measurement/step.c`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a. If h < S, square x and
  double (a, b). Otherwise multiply x by T[h][k] = sigma^(-k)(R_h), where R_h = g^(c_h) z^(d_h), and update (a, b) by lambda^(-k) (c_h, d_h).
  The Frobenius map is never applied per step.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle; a new
  key equal to the walk's anchor key (recorded every 64 of its own steps) inside a window of 64 steps is a cycle of length at most 64.
  Either escapes by squaring the least-key state of the cycle. Escapes are charged.
- Trail cap and restart: the walk's trail age counts main steps since its last distinguished point. At age 2^66 the walk restarts from a
  fresh random start (charged) and its age resets.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of the
  Montgomery integer; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index.
- Table and solving: the records are stored in the bitonic table of Section 4 and sorted by fingerprint. Records with equal fingerprints are
  compared pairwise. A pair from one class with B1 != B2 gives y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is discarded at the
  cost of one comparison; no trail is replayed. Every candidate is verified by g^y = z before output.
- Work counter and cutoff: each operation adds its charge (Section 3). The walks run while W < W_MAX - reserve, with reserve the solve and
  output charge. The table stops the procedure when it is full.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state (replica at production parameters, 512 and 4096 rounds, 64
walks); the canonical exponents of the stored distinguished points satisfy sigma^kn(x) = g^(lambda^kn (a + b y)) (replica, bit for bit
on the stored records); the twisted table equals sigma^(-k)(g^c z^d) (test); the solver recovers y end to end (tests, and 100
toy trials, all verified). The bitonic network is a standard sorting network; `measurement/sortnet/check_bitonic.py` checks that it sorts
random 128-bit keys for every power of two up to 1024 slots and that its compare-exchange count is n k (k+1)/4 (output in
`measurement/sortnet/check_bitonic.out`). Its cost is measured in the standalone kernel (Section 4).

## 2. Success probability, failure terms and the classification

Under H1 the walk needs T = 2^125.836396 fresh steps for tau = T^2 V / (2N). The collision rate is taken at the realised square fraction,
which is within 4.604e-05 of s (Section 3): V_low = 0.799372540 (nominal V = 0.799413681). The collisions the distinguished-point test can see are a
fraction f = 1 - 2 (1 + 64) e^(-64) of all collisions (trail cap). Main steps M = 2^125.844673 include the detected-cycle waste w = 0.005720
(127 steps per detection, at most 4.5040e-05 detections per step) and a detection lag of 2^80 steps. Then

  P = 1 - exp(-f tau_eff) - delta, tau_eff = T^2 V_low / (2N) = 0.527632742, delta = 2.16e-18.

P = 0.410000000001 before delta and 0.410000000001 after it. If the wasted steps were not budgeted, P would fall to 0.406438.

Classification. A closed window is a run of steps that returns the walk to an earlier class or state. (a) If the window contains a square,
its exponent coefficient carries a factor 2^q with q >= 1 that is not a 12th root of unity (below), so the coefficient is nonzero; such a
closure is an event of probability 1/r over the coins (H5). (b) If the window is all additive, the coefficient of adder h is
sum_k n_{h,k} lambda^(-k) in F_r, where n_{h,k} counts the steps of adder h with rotation k. For windows of length L < 2^63 this vanishes
in F_r only if sum_k n_{h,k} zeta^(-k) = 0 in Z[zeta_12]: the image of that element is zero, so it lies in a prime of Z[zeta_12] above r
(r = 1 mod 12, so r splits), and a nonzero element of that prime has norm divisible by r with |N| <= L^4 < 2^252 < r, a contradiction.
Proposition G (below) then gives the gon decomposition of each adder's rotation multiset.

**Proposition G.** Let n in Z_{>=0}^12 satisfy sum_k n_k zeta^k = 0 for a primitive 12th root of unity zeta. Then n is a nonnegative integer
combination of the indicator vectors of the 2-gons {k, k+6} and the 3-gons {k, k+4, k+8} (indices mod 12).

*Proof.* (i) Let Phi(x) = x^4 - x^2 + 1. It vanishes at zeta, since zeta^6 = -1 gives 0 = (zeta^2 + 1)(zeta^4 - zeta^2 + 1) with zeta^2 != -1.
It is irreducible over Q: a factorisation over Z (Gauss's lemma, Phi monic) has no linear factor, as Phi(+-1) = 1; a quadratic factorisation
(x^2 + ax + b)(x^2 - ax + c) forces b + c - a^2 = -1, a(c - b) = 0 and bc = 1, so b = c = 1 with a^2 = 3, or b = c = -1 with a^2 = -1, both
impossible. So Phi is the minimal polynomial of zeta, and sum_k n_k zeta^k = 0 iff the remainder of f(x) = sum_k n_k x^k modulo Phi is zero.

(ii) Modulo Phi, x^4 = x^2 - 1, x^5 = x^3 - x, x^6 = -1, and x^(k+6) = -x^k. In the basis 1, x, x^2, x^3 the coordinates of x^k are
v_0 = (1,0,0,0), v_1 = (0,1,0,0), v_2 = (0,0,1,0), v_3 = (0,0,0,1), v_4 = (-1,0,1,0), v_5 = (0,-1,0,1), v_6 = (-1,0,0,0), v_7 = (0,-1,0,0),
v_8 = (0,0,-1,0), v_9 = (0,0,0,-1), v_10 = (1,0,-1,0), v_11 = (0,1,0,-1). The remainder vanishes iff

  (C1) n_0 - n_4 - n_6 + n_10 = 0,   (C2) n_1 - n_5 - n_7 + n_11 = 0,
  (C3) n_2 + n_4 - n_8 - n_10 = 0,   (C4) n_3 + n_5 - n_9 - n_11 = 0.

(iii) Put p_j = n_j - n_(j+6) for j = 0..5. Then (C1)-(C4) read p_0 = p_4, p_1 = p_5, p_2 = -p_4, p_3 = -p_5. With a = p_4 and b = p_5,
the signed antipodal differences are p = (a, b, -a, -b, a, b).

(iv) Induction on the total weight. The zero vector is the empty combination. Let n != 0 be a nonnegative solution. If some pair has
n_j > 0 and n_(j+6) > 0, the 2-gon {j, j+6} lies in the support; subtracting its indicator keeps n nonnegative and a solution, since the gon
vector is a solution (v_j + v_(j+6) = 0). Otherwise n_j n_(j+6) = 0 for every j, so p_j = n_j or p_j = -n_(j+6), and n != 0 forces (a, b) != (0, 0).
If a > 0, then n_0 = n_4 = n_8 = a (from p_0 = a, p_4 = a, p_2 = -a), so {0, 4, 8} lies in the support. If a < 0, then n_2 = n_6 = n_10 = -a and
{2, 6, 10} lies in the support. If a = 0 and b > 0, then {1, 5, 9} does; if a = 0 and b < 0, then {3, 7, 11} does. The 3-gon indicator is a
solution, since zeta^r (1 + zeta^4 + zeta^8) = 0 (zeta^4 is a primitive cube root of unity). Subtracting it keeps n nonnegative and a solution.
Induction on the weight completes the proof. ∎

Per adder, the closed window's rotation multiset is a nonnegative solution, so Proposition G gives its decomposition. The union bound of
`accounting.py` (union_log2, lengths 65 to 400) counts the gon decompositions of each length and needs no enumeration. For windows longer than
400 the probability that every step is additive, (1 - s)^401 / s, is used, and no classification is needed.

Order of 2 and the smallest q. r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 * 254760293^2
(checked by `cscore/number_theory.py`). The multiplicative order of 2 modulo r is o = 3277242198445386904965483781761622864855659531282977363912728668746161324032, 251 bits, certified by 2^o = 1 and
2^(o/p) != 1 for each of the prime divisors p of o. Since 2^q is a 12th root of unity iff o | 12q iff q0 | q with q0 = o / gcd(o, 12) (the
two factors are coprime: gcd(q0, 12/gcd(o,12)) = 1), the smallest q >= 1 is q0, a 248-bit number (gcd(o, 12) = 12). Hence 2^q is not a
12th root of unity for 1 <= q < 2^247, and a window containing a square has a nonzero coefficient. The test 2^((r-1)/p) for the primes p dividing
r - 1 is reported for information (it shows that 2 is not a primitive root).

Failure terms (each a probability over the run; computed in `accounting.py`, function `build`):
| term | value | source |
| --- | ---: | --- |
| formal windows, lengths 65 to 400 (per step) | 1.77e-96 | union bound over gon decompositions (Proposition G, Section 2) |
| all-additive windows longer than 400 (per step) | 1.25e-103 | (1 - s)^401 / s |
| non-formal closures (per step) | 1.69e-56 | 12 CAP / r candidates, each closing with probability 1/r (heuristic H5) |
| trap term over the run | 1.29e-18 | M times the per-step total |
| fingerprint Markov (metered matches) | 8.67e-19 | 2^-60 |
| square fraction, Hoeffding over fresh steps | 2^(-2.2e+26) | exp(-2 (1e-6)^2 M), with the revisit allowance in V_low |
| detected cycles of length <= 64 (residue-class Chernoff) | 2^(-2.59e+25) | 64 residue classes, conditional domination (Section 3) |
| non-formal detections of length <= 64 | 1.1e-36 | expectation 12 * 64 * M / r (Markov) |
| distinguished-point table overflow | 2^(-3.19e+19) | exp(-mu/3) on the fresh flags, with the revisit compensator (Section 4) |
| cross-walk collisions above the allowance of 64 | 2^(-329.0) | Poisson tail with mean M^2/(2N) |
| total failure kept | 2.16e-18 | sum of the kept terms |

Heuristics (disclosed). H1 (random-map collision model): the walk behaves as a random map on N classes with collision rate V_low;
tau is treated as Exp(1). H5 (non-formal closures): a closure with a nonzero coefficient has probability 1/r over the coins. H7 (conditional
uniformity): given the history, the adder index h = K mod r_a of each fresh step is uniform, the distinguished-point and fingerprint bits of a
fresh state are uniform, and escapes and anchors are deterministic functions of the history. The square count, the detection count and the
distinguished-point count use H7 through conditional probabilities; no independence across steps is assumed.

## 3. Computation and the work cutoff

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), measured 22,998 instructions (1.1499 units; C = 20,000 is the organizer's
rounding), charged 1 as the rules require. Other work is measured in instructions and divided by C.

Measurement (`measurement/`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0, `-O3 -fno-tree-vectorize
-mgeneral-regs-only`), the pinned blst commit f7a88067. Callgrind Ir is differenced between 256 and 512 rounds of 64 walks at S = 458, with the
distinguished-point density set to zero (DP work is charged separately). The bitonic network is measured in `measurement/sortnet/`
(sortnet.c, with its own main, so the walk kernel `step.c` is unchanged): c_cmp from the difference between 4 and 8 sorts of 4096 slots, the
scan from sortscan minus sortnet.

| forced step type (S = 458 or forced; zero DP density) | per-step instructions | notes |
| --- | ---: | --- |
| all-adder steps, seed 1 | 26,730.16 | includes 4 detected cycles (conservative) |
| all-adder steps, seed 2 | 26,722.87 | includes 1 detected cycle |
| all-square steps, seed 1 | 13,511.44 | no detected cycles |
| all-adder steps with a distinguished point at every step | 28,209.16 | the store costs 1,479.0 instructions per point |
| reference operation (blst_fp12_mul) | 22,998 (organizer: 22,996) | charged 1 unit |
| cyclotomic square (blst_fp12_cyclotomic_sqr) | 9,804 | charged 9,804 / C = 0.4902 |
| bitonic compare-exchange (sortnet, 4096 slots, two seeds and 8192 slots) | 108.003 (108.001 at 8192) | branch-free mask compare and 80-byte conditional swap |
| bitonic scan (sortscan minus sortnet) | 13.0 per slot | equal-neighbour test, branch-free |

Per-type charges. The step code is straight-line (comparisons and selections without data-dependent branches), so each type has one
instruction count; the two all-adder seeds agree to 0.03%, which is evidence for the constancy of the path, not a proof for every machine
state. An adder step is charged 1.18771 units (the reference operation, 1, plus the other instructions divided by C), and a square step is
charged 0.67657 units. Every main step is charged the larger, 1.18771 units, which bounds both.

Charges for the first M main steps (each bound holds with the stated probability; H7 as in Section 2):
- Operation mix: no concentration is used for the charge. The square decisions of the fresh steps are conditionally Bernoulli(s), so their
  count is within 10^-6 M of sM with probability at least 1 - exp(-2 * 10^-12 M) (Hoeffding). A revisit step repeats an earlier decision;
  revisits are at most the detections, so the realised square fraction is within 4.604e-05 of s (V_low in Section 2).
- Detected cycles: a detection at step i needs the window ending at i (length at most 64, from an anchor or a 2-cycle) to close. Given the
  history F_(i-64), the probability is at most p = 4.5040e-05 (the 2-cycle term (1 - s)/(12 r_a), the formal union over lengths 3 to 64, and the
  non-formal term), with 1 + 10^-3 slack. Split the steps into 64 residue classes mod 64. Along each class the indicators have conditional
  probability at most p given the history 64 steps earlier, so they are dominated by independent Bernoulli(p) variables; Chernoff on each class
  and a union bound over the 64 classes give failure at most 64 exp(-(10^-3)^2 (p M / 64) / 3) (2^(-2.59e+25)). Each detection costs at most
  66 adder-type operations (the walk-through of at most 64 steps, the check and the escape), charged 66 A_UNITS: the event charge is
  4.5040e-05 * 66 * 1.18771 per step.
- Non-formal detections of length at most 64: expected count 12 * 64 M / r (Markov, failure 1.1e-36).
- Distinguished points: charged at capacity (Section 4), so the charge does not depend on the random count.
- Fingerprint matches: metered when they occur (one exponentiation each); the metered total is bounded by 2^60 times the expected count
  (Markov, failure 2^-60 = 8.67e-19). Above that bound the procedure stops with failure.
- Restarts: at most M / 2^66 + 2^20 restarts, each charged 2 exponentiations.

Components of W_MAX (log2 units; the score is the rounded-up cutoff 126.097138):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 126.092855 | 125.844673 log2 main steps, each charged the adder cost 1.18771 units (the larger of the two step charges) |
| detected cycles | 117.698817 | at most 4.5040e-05 detections per step, each charged 66 adder-type operations |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts and restarts | 69.418320 | 2^20 starts and at most M/2^66 trail-cap restarts (preprocessing) |
| distinguished points | 69.637665 | bitonic table over 2^67.0000 slots, charged for its capacity 2^66.9480 records (store, network, scan, overflow) |
| fingerprint verification | 73.469571 | Markov charge at 2^60 times the expected count of about 2^4.90 false pairs, one exponentiation each |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

Probability that the first M main steps are charged within W_MAX: at least 1 - 2.16e-18, from the terms of Section 2. Consequently the
procedure executes M main steps before the cutoff, and the success bound of Section 2 applies.

## 4. Distinguished-point table, scratch memory and peak memory

Capacity. Let mu = M 2^-60 and c = (WINDOW + 1)(1 - s)/r_a = 0.0691 (WINDOW = 127). The capacity is D = ceil(2 (mu + 64)/(1 - c)) + 2 * 2^20,
which is 2^66.9480 records. Justification: the fresh flags X have mean mu under H7 and are dominated by Binomial(M, 2^-60); a flagged revisit
lands on one of the last 128 states, and an adder step reaches a given earlier state with probability at most (1 - s)/r_a, so a flagged state
gathers at most c expected revisit flags while it is visible; cross-walk collisions add at most 64 flags with failure (e mu_c / 64)^64, where
mu_c = M^2/(2N). So D_tot <= (X + 64)/(1 - c). The factor 2 is the Chernoff threshold: overflow from the fresh flags needs X > 2 (mu + 64), which has
probability at most exp(-mu/3) (Chernoff on X), plus the collision term. The 2^20 term is an allowance for the detection lag of each walk.

Procedure (walk-spec section 9.2).
- Store: a distinguished point appends its 80-byte record to the table; if the table holds D records, the procedure stops (overflow rule).
- Pad: the slots after the last record are filled with sentinels (fingerprint all ones), so they sort last; n = 2^67.0000 slots (the smallest
  power of two not less than D).
- Sort: a bitonic sorting network on the n slots, ascending by the 128-bit fingerprint. Each compare-exchange computes the comparison as a mask
  and swaps the two records with masked XOR through one 80-byte temporary. There are n k (k + 1)/4 = 2^77.1536 compare-exchanges, k = log2 n.
  The index arithmetic depends only on n.
- Scan: a constant-cost pass over the sorted slots tests equal neighbours; each equal pair is enumerated and charged as a fingerprint match.
- Cost (instructions): store 1,479 per record (measured), compare-exchange 108.0030 (measured), scan 13.0 per slot (measured), and 10 per record
  for the overflow test. The charge is 1,479 D + 108.0030 * n k(k+1)/4 + 13.0 n + 10 D, divided by C. Per record this is 6.452 units.
- Scratch memory: the sort is in place. The scratch is the 80-byte temporary and the padding (n - D records at most), both inside the n-slot table.

Peak memory. The table occupies n * 80 bytes = 2^67.0000 * 80 bytes (padding included), plus the 80-byte temporary, the adder table (12 x
1024 elements, about 7.8 MiB) and the walk states (2^20 x 640 bytes) with the detector state (98 + 16 bytes) per walk. The total is at most
2^73.321929 bytes. The peak is not realisable on the hardware used here and is reported, not scored.

## 5. Evidence: toy-scale sweeps

Group: the order-ell subgroup of F_(q^12)^* with ell dividing q^4 - q^2 + 1; `cscore/gt.py`, bundled into `experiment.py`. The toy taus use
the partition-correct V. The toy walk detects full cycles with its own rules and uses dense distinguished points, so these runs test H1 and
the solver, not H2, H5 or the production detector.

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
- H5 is a heuristic: the non-formal closure term depends on the coins being independent of the path. H7 is the conditional-uniformity premise
  of Section 2; it is a stated premise, not a proof. The argument for the detection count uses only the conditional bound, and the argument for
  the distinguished-point capacity uses the fresh flags and a compensator for revisits.
- The per-type step counts are measured on forced inputs; the code path has no data-dependent branch. The cycle events are bounded per
  detection (66 adder-type operations). The bitonic network is measured for cost and checked for sorting up to 1024 slots
  (`measurement/sortnet/check_bitonic.py`); the production table has up to 2^67.0000 slots, and the sorting property there is that of the
  standard network, not a test at that size. The table is not part of the walk kernel.
- The trail-cap restart and the 3-to-64 cycle branch are specified or mirrored but not exercised by an exact-equality test (measurement runs
  are 2^26 steps per seed; the 3-to-64 branch fired once).
- The solve, verify and output charge of 2^29 instructions is an assumed bound; it is 2^14.7 of the 2^126.097138 units.
- The reference operation is priced by the organizer rule at 1 unit; the measured instruction count is 1.1499 units.
- Preflight: the shared verifier (`llm_verifier`) is not installed on the machine used for this revision, so the preflight step could not be run
  here. The declared experiments (`experiment.py`, the accounting entry) were run locally with the stdlib (Section 7 of the package notes).

## References

- Lam and Leung, On vanishing sums of roots of unity (2000), for the background of Proposition G (the statement is proved here).
- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py                 # bundled: W_MAX, failure terms, number-theory asserts
python3 -I number_theory.py              # the number-theory check alone (evidence/number_theory_check.out)
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image, py_ecc (gen_inputs.py)
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind at zero DP density, two sizes, mul and sqr loops), analyses the
interval counts, and runs the replica on each input set. `measurement/sortnet/` measures the bitonic compare-exchange and scan.
`native_run.sh` runs the uninstrumented kernel for the event counts.

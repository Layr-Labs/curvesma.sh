# frobenius-mixed-rho-v3: Frobenius-class mixed add/square rho on GT (revision v8)

Claim: for the uniform challenge of `rules.md` (P uniform nonidentity in G1, z uniform in GT, fresh algorithm coins), the procedure below
outputs Q in G2 with e(P, Q) = z with probability at least 0.41 (model 0.410000000001; lower bound 0.410000000001 after the failure terms of
Section 2), while its total work, counting every operation it performs, is at most W_MAX = 2^126.102185 target-group-operation
equivalents (score 25220437/200000). The first M = 2^125.849293 main steps are executed before the cutoff except with probability
2.16e-18 (Section 3), under the premises H1, H5 and H7 stated in Section 2. Peak memory is at most 2^73.321929 bytes, including the
padding and scratch of the distinguished-point table and the per-walk start elements (Section 4); preprocessing is inside W_MAX
(2^29.575064 units); nonuniform advice is 0 bytes. It is a classical generic attack, not a practical break.

**Revision v8.** Round 8 (on revision v7) rejected the probability accounting on two points. v8 changes the walk, not only the
bounds.

(1) Trails end at a distinguished point. A walk stores its record at its first DP and restarts at its own start element (one
multiplication, charged), so a trail is a single segment between two DPs. A collision between two trails (a join) makes the later
walk follow the other trail to its DP, where it stores one duplicate and restarts at a different start element, so the two trails
diverge. Records are therefore at most X + J, pathwise: X is the number of fresh DP flags (dominated by Binomial(M, 2^-60)) and J
the number of joins (each join at most one duplicate; dominated by Binomial(M, M/N) under H1). The capacity is D = ⌈2 mu⌉ + 64,
without the v7 compensator, which had used an expectation as a pathwise bound (Section 4).

(2) The square fraction counts revisit steps correctly. A closure event is followed by at most WINDOW = 127 revisit steps before its
detection (the anchor fires within 127 steps of the closure), not by one; the square-fraction allowance is therefore 127 p_ev, which
is 0.00572 (v7 used p_ev alone, 4.6e-5, which was too small). The success model uses the lower collision rate V_low at this allowance.

(3) The closure-event bound is stated with its premise. Given the history up to the start of a window of 64 steps, a closure inside the
window has probability at most p_ev, because the window's states are fresh until the closure and their keys are uniform given the past
(H7). The residue-class Chernoff bound (64 classes mod 64) uses this conditional bound, not independence across steps.

(4) Score: v7 126.097138 to v8 126.102185 (+0.312501 bits). The increase is from the square-fraction allowance (V_low) and the restart and
join charges. The failure total is 2.16e-18.

## 1. Procedure

Reduction. g = e(P, Q0) generates GT, where Q0 is the fixed G2 generator; z = g^y gives Q = [y]Q0 with e(P, Q) = z. Input a nonidentity P in G1 and a target z in GT. Return the identity if z = 1 (Q = O is valid).
Otherwise run the walks until a verified y is found or the work counter reaches W_MAX, then output Q = [y]Q0, verified by the pairing
equation e(P, Q) = z (one pairing, charged). Exceptional inputs (P = O) are rejected by the rules.

Walk (`walk-spec.md`, the specification used by the blst kernel `measurement/step.c`, revision v8):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a. If h < S, square x and
  double (a, b). Otherwise multiply x by T[h][k] = sigma^(-k)(R_h), where R_h = g^(c_h) z^(d_h), and update (a, b) by lambda^(-k) (c_h, d_h).
  The Frobenius map is never applied per step.
- Trails: a walk's trail runs from its start element (X0, A0, B0), fixed per walk with X0 = g^(A0 + B0 y), to its first distinguished
  point. At the DP the walk stores its record and restarts: x <- x * X0, (a, b) <- (a + A0, b + B0). The restart clears the walk's
  detector (previous key, anchor, trail age). A trail that reaches 2^66 steps without a DP restarts in the same way, without a record.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle; a new
  key equal to the walk's anchor key (recorded every 64 of its own steps) inside a window of 64 steps is a cycle of length at most 64.
  Either escapes by squaring the least-key state of the cycle. Escapes are charged.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of the
  Montgomery integer; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index. The record is stored before the restart.
- Table and solving: the records are stored in the bitonic table of Section 4 and sorted by fingerprint. Records with equal fingerprints
  are compared pairwise. A pair from one class with B1 != B2 gives y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is discarded at
  the cost of one comparison; no trail is replayed. Every candidate is verified by g^y = z before output.
- Work counter and cutoff: each operation adds its charge (Section 3). The walks run while W < W_MAX - reserve, with reserve the solve and
  output charge. The table stops the procedure when it is full.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state, including after restarts (replica at production parameters
and at small DP densities that exercise the restarts, bit for bit on the states and the records; see `measurement/`); the canonical exponents
of the stored records satisfy sigma^kn(x) = g^(lambda^kn (a + b y)); the twisted table equals sigma^(-k)(g^c z^d) (test); the solver
recovers y end to end (tests, and 100 toy trials, all verified, but see Section 6 for the toy walk's semantics). The bitonic
network is a standard sorting network; `measurement/sortnet/check_bitonic.py` checks that it sorts random 128-bit keys for every power of
two up to 1024 slots and that its compare-exchange count is n k (k+1)/4.

## 2. Success probability, failure terms and the classification

Under H1 the walks need T = 2^125.841017 fresh steps for tau = T^2 V / (2N). The collision rate is taken at the lower end of the realised square
fraction: the realised fraction is within eps_s = 0.00572 of s (Section 3), so V_low = 0.794268849 (nominal V = 0.799413681). The collisions the
distinguished-point test can see are a fraction f = 1 - 2 (1 + 64) e^(-64) of all collisions (trail cap). Main steps M = 2^125.849293 include the
waste w = 0.005720 (127 revisit steps per closure event, at most 4.5040e-05 closure events per step, and the join suffixes) and a detection lag of
2^80 steps. Then

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
(checked by `cscore/number_theory.py`, output in `evidence/number_theory_check.out`). The multiplicative order of 2 modulo r is
o = 3277242198445386904965483781761622864855659531282977363912728668746161324032, 251 bits, certified by 2^o = 1 and 2^(o/p) != 1 for each of the prime divisors p of o. Since 2^q is a 12th root of unity iff
o | 12q iff q0 | q with q0 = o / gcd(o, 12) (the two factors are coprime: gcd(q0, 12/gcd(o,12)) = 1), the smallest q >= 1 is q0, a 248-bit
number (gcd(o, 12) = 12). Hence 2^q is not a 12th root of unity for 1 <= q < 2^247, and a window containing a square has a nonzero
coefficient. The test 2^((r-1)/p) for the primes p dividing r - 1 is reported for information (it shows that 2 is not a primitive root).

Failure terms (each a probability over the run; computed in `accounting.py`, function `build`):
| term | value | source |
| --- | ---: | --- |
| formal windows, lengths 65 to 400 (per step) | 1.77e-96 | union bound over gon decompositions (Proposition G, Section 2) |
| all-additive windows longer than 400 (per step) | 1.25e-103 | (1 - s)^401 / s |
| non-formal closures (per step) | 1.69e-56 | 12 CAP / r candidates, each closing with probability 1/r (heuristic H5) |
| trap term over the run | 1.29e-18 | M times the per-step total (includes the trap revisit steps) |
| fingerprint Markov (metered matches) | 8.67e-19 | 2^-60 |
| square fraction, Hoeffding over fresh decisions | 2^(-2.21e+26) | exp(-2 (1e-6)^2 M); the revisit allowance is in V_low |
| closure events of length <= 64 (residue-class Chernoff) | 2^(-2.59e+25) | 64 residue classes, conditional bound p_ev (Section 3) |
| non-formal closures of length <= 64 | 1.1e-36 | expectation 12 * 64 * M / r (Markov) |
| DP table, fresh flags above 2 mu | 2^(-3.2e+19) | exp(-mu/3), Chernoff on X (Section 4) |
| DP table, joins above 64 | 2^(-264.4) | (e mu_J / 64)^64 with mu_J = M^2/N (Section 4) |
| total failure kept | 2.16e-18 | sum of the kept terms |

Success chain (what must hold for the success probability; each step names its failure term). (i) The procedure executes the first M
main steps before the cutoff: the charges of Section 3 are within W_MAX, with the failure total of the charge terms. (ii) No DP table overflow
occurs: the records are at most X + J (Section 4), overflow has probability at most the two DP terms below. (iii) The useful work is enough:
the revisit steps among the first M steps are at most w M (closure events by their concentration, joins by their tail, trap events by the
trap term), so the fresh steps are at least (1 - w) M >= T by the choice of M. (iv) Given T fresh steps, a collision occurs between walks with
probability 1 - exp(-f tau_eff) under H1 (the collision rate V_low, Section 2). (v) A collision is found: each walk records its DP and the
records of the two walks with equal fingerprints form a pair with B1 != B2 except with probability 1/r (degenerate pairs are negligible), so
the solver returns a verified y. A join is itself a collision: the joined walk's duplicate record is paired with the other walk's record. The
collision is detected once both walks reach a DP, which is the detection lag 2^80 steps folded into M. Each step of (i)-(v) has a stated
failure term, and their sum is the failure total of Section 2.

Premises (disclosed; each is a stated assumption, not a proof).
- **H1 (random map).** The walk behaves as a random map on the N classes: the next state of a step is uniform over the classes given the
  past, so the probability that it is a visited state of another trail is at most (visited)/N. tau is treated as Exp(1), with the collision
  rate V_low. Joins and the distinguished-point duplicates use this bound (Section 4). Its adder-level form: a step's additive index h = K mod r_a
  is uniform over the r_a values, so a given admissible additive index is chosen with probability 1/r_a (the square decision and the index are
  the same uniform key bits, so the index probability is not multiplied by 1 - s); with the adder elements independent of the path (H5), a given
  state is an adder neighbour of the current state with probability about r_a/N, which gives the (visited)/N bound.
- **H5 (non-formal closures).** A closure with a nonzero coefficient has probability 1/r over the coins. This is the only premise of the trap
  term.
- **H7 (uniform fresh states).** Given the history, the key of a fresh state (one no walk has visited) is uniform: its adder index h = K mod r_a
  and rotation k are uniform, its square decision has probability s, and its distinguished-point flag has probability 2^-60. Escapes and restarts
  are deterministic functions of the history (a restart is multiplication by the walk's start element). Consequences used: the closure
  probability p_ev of Section 3 holds conditionally on the history up to the start of the window (the window's states are fresh until the closure),
  the fresh flags X are dominated by Binomial(M, 2^-60), and the square decisions of fresh states are Bernoulli(s).

## 3. Computation and the work cutoff

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), measured 22,998 instructions (1.1499 units; C = 20,000 is the organizer's
rounding), charged 1 as the rules require. Other work is measured in instructions and divided by C.

Measurement (`measurement/`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0, `-O3 -fno-tree-vectorize
-mgeneral-regs-only`), the pinned blst commit f7a88067. Callgrind Ir is differenced between 256 and 512 rounds of 64 walks at S = 458, with the
distinguished-point density set to zero (DP work is charged separately). The restart of Section 1 is measured with the all-adder run at DP
density one, minus the zero-density run. The bitonic network is measured in `measurement/sortnet/` (its own main; the walk kernel is separate).

| forced step type (production DP_BITS 60, zero distinguished points) | per-step instructions | notes |
| --- | ---: | --- |
| all-adder steps, seed 1 | 26,731.234 | v8 kernel; includes 4 detected cycles (conservative) |
| all-adder steps, seed 2 | 26,723.934 | v8 kernel; includes 1 detected cycle |
| all-square steps, both seeds | 13,512.5 | no detected cycles |
| all-adder steps, distinguished point (restart) at every step (DP_BITS 0), both seeds | 54,658.5 | per step: store, restart and adder |
| DP path difference (every-step minus zero-density), seeds 1 and 2 | 27,927.266 / 27,934.566 | includes one restart per step |
| restart alone (mode restartbench) | 26,507.0 | one multiplication by the start element and key derivation |
| store of one record (difference minus restart, larger seed; conservative) | 1,427.6 | charged as the store |
| bitonic compare-exchange (standalone sort bench, 4096 slots) | 108.003 (108.001 at 8192) | branch-free mask compare and 80-byte conditional swap |
| bitonic scan, per slot | 13.0 | equal-neighbour test, branch-free |
| reference operation (organizer, REF_FUNCS) | 22,990 | charged 1 unit; blst_fp12_mul measured 22,997 |
| cyclotomic square | 9,803 | charged 9,803 / C |

Per-type charges. The step code is straight-line (comparisons and selections without data-dependent branches), so each type has one
instruction count, measured on forced inputs. An adder step is charged 1.18806 units (the reference operation, 1, plus the other instructions
divided by C), and a square step is charged 0.67663 units. Every main step is charged the larger, 1.18806 units, which bounds both.

Charges for the first M main steps (each bound holds with the stated probability; H1, H5 and H7 as in Section 2):
- Operation mix: no concentration is used for the charge. The square decisions of the fresh steps are conditionally Bernoulli(s) (H7), so their
  count is within 10^-6 M of sM with probability at least 1 - exp(-2 * 10^-12 M) (Hoeffding). Revisit steps repeat an earlier decision. Each closure
  event is followed by at most 127 revisit steps, and the join suffixes are at most CAP each, so the realised square fraction is within
  eps_s = 0.00572 of s (V_low in Section 2).
- Closure events: at most 4.5040e-05 per step (2-cycle term (1 - s)/(12 r_a), formal union over lengths 3 to 64, and the non-formal term), with
  1 + 10^-3 slack, conditionally on the history up to the start of each window (H7, Section 2). Split the steps into 64 residue classes mod 64;
  the events of a class are dominated by independent Bernoulli(p_ev) variables. Chernoff on each class and a union bound over the classes give
  failure at most 64 exp(-(10^-3)^2 (p_ev M / 64) / 3) (2^-2.59e+25). Each event costs at most 66 adder-type operations (the walk-through of at
  most 64 steps, the check and the escape), charged 66 A_UNITS: the event charge is 4.5040e-05 * 66 * 1.18806 per step.
- Non-formal closures of length at most 64: expected count 12 * 64 M / r (Markov, failure 1.1e-36).
- Restarts: each distinguished point and each trail-cap reset costs one multiplication (charged A_UNITS). The number of records is at most D
  (Section 4) and the cap resets are at most M / 2^66, both charged at their bounds.
- Fingerprint matches: metered when they occur (one exponentiation each); the metered total is bounded by 2^60 times the expected count
  (Markov, failure 2^-60 = 8.67e-19). Above that bound the procedure stops with failure.

Components of W_MAX (log2 units; the score is the rounded-up cutoff 126.102185):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 126.097903 | 125.849293 log2 main steps, each charged the adder cost 1.18806 units (the larger of the two step charges) |
| detected cycles | 117.703864 | at most 4.5040e-05 closure events per step, each charged 66 adder-type operations |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts | 29.573647 | 2^20 start elements, two exponentiations each (preprocessing) |
| restarts at the trail cap | 60.097903 | at most M / 2^66 cap resets, one multiplication each (charged A) |
| distinguished points | 69.891450 | table of 2^67.0000 slots, charged for its capacity 2^66.8493 records (store, restart, network, scan, overflow) |
| fingerprint verification | 73.272233 | Markov charge at 2^60 times the expected count of about 2^4.70 false pairs, one exponentiation each |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

Probability that the first M main steps are charged within W_MAX: at least 1 - 2.16e-18, from the terms of Section 2. Consequently the
procedure executes M main steps before the cutoff, and the success bound of Section 2 applies.

## 4. Distinguished-point table, pathwise capacity, scratch memory and peak memory

Capacity (pathwise). Each record is produced at a DP step. A DP step is either a fresh flag (the state is visited for the first time by any
walk) or a duplicate (the walk has joined another trail and reaches that trail's DP, which the other trail has already recorded). A walk joins at
most once per trail (after the join it follows the other trail to its DP, stores one duplicate, and restarts at its own start element), so the
duplicates are at most the joins J. Hence the records are at most X + J, pathwise, where X is the number of fresh flags.
- X is dominated by Binomial(M, 2^-60) (H7).
- J is dominated by Binomial(M, M/N), since a step's next state is visited by another trail with probability at most (visited)/N <= M/N (H1).
The capacity is D = ⌈2 mu⌉ + 64 with mu = M 2^-60, which is 66.8493 log2 records. Overflow requires X > 2 mu (Chernoff: probability at most
exp(-mu/3) = 2^(-3.2e+19)) or J > 64 (probability at most (e mu_J / 64)^64 with mu_J = M^2/N, log2 -264.4). The table stops the procedure when it
is full, so only these events bind.

Procedure (walk-spec section 9.2): store the record of each DP (the table holds D records); when the walks stop, pad the table to n = 2^67.0000
slots (the smallest power of two not less than D) with sentinel fingerprints; sort with a bitonic network (n k (k+1)/4 = 2^77.1536
compare-exchanges, each a branch-free mask compare and an 80-byte masked swap through one temporary); scan equal neighbours and enumerate
equal runs, each pair charged as a fingerprint match.
- Cost (instructions): the store of a record (1427.6, measured), one restart (26731.2 charged, the adder cost), the network
  (108.0030 per compare-exchange, measured in the standalone kernel), the scan (13.0 per slot, measured), and 10 per record for the overflow
  test. The charge is (store + restart + 10) D + 108.0030 * n k(k+1)/4 + 13.0 n, divided by C. Per record this is 8.237 units.
- Scratch memory: the sort is in place; its scratch is the 80-byte temporary and the padding (at most n − D records), both inside the n-slot table.

Peak memory. The table occupies n * 80 bytes = 2^67.0000 * 80 bytes (padding included), plus the 80-byte temporary; the adder table (12 x 1024
elements, about 7.8 MiB); the walk states (2^20 x 640 bytes) with the detector state (98 + 16 bytes) per walk; and the per-walk start elements
(2^20 x 640 bytes, needed for the restarts). The total is at most 2^73.321929 bytes. The peak is not realisable on the hardware used here and is
reported, not scored.

## 5. Evidence: toy-scale sweeps

Group: the order-ell subgroup of F_(q^12)^* with ell dividing q^4 - q^2 + 1; `cscore/gt.py`, bundled into `experiment.py`. The toy taus use
the partition-correct V. The toy walk detects full cycles with its own rules and uses dense distinguished points without restarts, so these runs
test H1 and the solver only (they were run before the v8 restart change and are not re-run for it; see Section 6).

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
- H5 is a heuristic: the non-formal closure term depends on the coins being independent of the path. H7 is the conditional-uniformity premise of
  Section 2; it is a stated premise, not a proof. The closure count, the square fraction and the fresh-flag count use H7 through conditional
  bounds; the join count uses H1.
- The toy experiment (`experiment.py`) implements the walk of revision v7 (dense DP, no restart), so its sweeps do not test the restart
  semantics of v8. The restart is tested by the replica at small DP densities (`measurement/`), not by the toy runs.
- The per-type step counts are measured on forced inputs; the code path has no data-dependent branch. The cycle events are bounded per
  detection (66 adder-type operations). The bitonic network is measured for cost and checked for sorting up to 1024 slots; the production
  table has up to 2^67.0000 slots, and the sorting property there is that of the standard network.
- Trap events (undetected cycles longer than 64) are covered by the trap term; their revisit steps (up to the cap) are not in the square-fraction
  allowance, so they are included in the failure total only through the trap term.
- The solve, verify and output charge of 2^29 instructions is an assumed bound; it is 2^14.7 of the 2^126.102185 units.
- The reference operation is priced by the organizer rule at 1 unit; the measured instruction count is 1.1499 units.
- Preflight: the shared verifier (`llm_verifier`) is available on the lead's machine; the preflight outcome is reported with the package.

## References

- Lam and Leung, On vanishing sums of roots of unity (2000), for the background of Proposition G (the statement is proved here).
- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py                 # bundled: W_MAX, failure terms, number-theory asserts
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image, py_ecc (gen_inputs.py); includes the restart checks
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind at zero DP density, two sizes, mul and sqr loops, and the DP-every-step
run), checks the replica at production and small DP densities, analyses the interval counts, and runs the native counts. `measurement/sortnet/`
measures the bitonic compare-exchange and scan, and checks it sorts.

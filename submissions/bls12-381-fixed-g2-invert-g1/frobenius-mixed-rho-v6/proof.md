# frobenius-mixed-rho-v6 on bls12-381-fixed-g2-invert-g1: Frobenius-class mixed add/square rho on GT (revision v11)

Claim: for the uniform challenge of `rules.md` (Q uniform nonidentity in G2, z uniform in GT, fresh algorithm coins), the procedure below
outputs P in G1 with e(P, Q) = z with probability at least 0.41 (model 0.432781; lower bound 0.41 after the failure terms of
Section 2), while its total work, counting every operation it performs, is at most W_MAX = 2^126.915106 target-group-operation
equivalents (score 63457553/500000). The first M = 2^126.373955 main steps are executed before the cutoff except with probability
0.022781 (Section 3), under the premises H1, H5, H7, H8, H9 and the cap-based trap bound of Section 2. Peak memory is at most
2^74.321929 bytes, including the padding of the bucketed distinguished-point table (Section 4); preprocessing is inside W_MAX (2^29.575082
units); nonuniform advice is 0 bytes. It is a classical generic attack, not a practical break.

**Revision v11.** Round 11 answers the round-10 g2 review (`work-gt/REVIEW-round10-g2.md`). Four changes, each with its
consequence in the score (the score rises: the claimed work is larger, not smaller).

(1) Restart and escape exponents are PRF outputs (H8, walk-spec section 6). v10 restarted at one of 1024 fixed table elements
Y_j = g^(f_j), and the same elements were reused for every restart, so conditional uniformity after reuse was not established. v11
draws each restart exponent e = PRF(0x01, w, t) and each escape exponent e = PRF(0x02, K_s) from a keyed SHA-256 PRF with a random
256-bit key (the message is the walk index and the restart counter, or the least-key state's key). The restart multiplies by g^e through a fixed-base comb with 8-bit windows
(32 multiplications, 8160 precomputed elements), charged at the full adder cost. Distinct messages give independent uniform
exponents, so every restart or escape target is uniform on G given the history (H8), which is the standard PRF premise.

(2) The square decision and the adder index come from a keyed mix of the key (H9, walk-spec section 5). v10 took the decision from the
low bits of K mod r_a, and the toy showed a fresh square share below s (0.442 against 0.447 at r_a = 1024; the mix of a small toy key
is not uniform on its low bits). v11 computes u = mix(K_lo XOR k_b) and takes the square decision from the top ten bits of u and the
adder index from the bottom ten bits, so the square probability is exactly s = 458/1024. The collision rate is not derived in this
revision: the budget takes the non-degenerate collision rate V_low = 0.7996551059 as an assumption (Section 2, Section 5). The square
share enters only through its expectation under the random oracle (H9), with the Hoeffding deviation charged.

(3) The trail age is not reset by an escape (walk-spec section 6). In v10 the escape reset the trail age, so a trail that escaped
and re-entered its cycle could repeat the loop without reaching the trail cap. In v11 the cap bounds every trail, escape loops
included. A loop closes only on a class the trail has already visited, so every trap is an own-trail revisit. Each landing at trail
age a revisits a visited class with probability at most a/N under uniform landings (H1); summed over the M main steps the expected
number of own-trail revisits is at most M CAP / N, about 1.9e-18, a term of delta. The in-degree argument of v10 (H10) is withdrawn.

(4) Escapes and restarts are charged what they cost. v10 charged a restart one adder step; the escape now costs 66 adder operations and one
comb restart (32 multiplications plus two PRF evaluations). The additive cells are all 12 r_a (v10: 12 (r_a − S)); the formal-closure union
bound uses them (walk-spec section 9.5). The 2-cycle rate is unchanged, (1 − s)/(12 r_a).

Score: v10 126.840188 to v11 126.915106 (+0.074918 bits; the failure total moves from 0.022848 to 0.022781). The walk-step
term falls (126.640 to 126.627 log2 units); the closure events (123.890 to 124.451), the cap restarts (60.64 to 65.64) and the
distinguished points (70.47 to 72.89) rise (Section 3). The comb table (13.25 log2 units) is new and is charged in the preprocessing.

## 1. Procedure

Reduction. g = e(P0, Q) generates GT, where P0 is the fixed G1 generator; z = g^y gives P = [y]P0 with e(P, Q) = z. Input a nonidentity Q in G2 and a target z in GT. Return the identity if z = 1 (P = O is valid).
Otherwise run the walks until a verified y is found or the work counter reaches W_MAX, then output P = [y]P0, verified by the pairing
equation e(P, Q) = z (one pairing, charged). Exceptional inputs (Q = O) are rejected by the rules.

Walk (`walk-spec.md`, revision v11, implemented by the blst kernel `measurement/step.c` and its replica `measurement/replica.py`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant. K_lo is its
  low 64-bit word.
- Branch and index: u = mix(K_lo XOR k_b), k_b a 64-bit key drawn at setup, mix the splitmix64 finaliser. Square slots S = 458 of
  r_a = 1024 (s = 0.447265625 exactly): a step squares x and doubles (a, b) iff (u >> 54) < S. Otherwise it multiplies x by
  T[h][k] = sigma^(-k)(R_h), h = u mod r_a, where R_h = g^(c_h) z^(d_h), and updates (a, b) by lambda^(-k) (c_h, d_h).
- Trails: a trail ends at its first distinguished point (DP) or at the trail cap, 2^66 steps of trail age. The trail age is reset
  when a trail begins (walk start, DP restart, cap restart), and not by an escape.
- Restart (at a DP or the cap): the walk's counter t is incremented, e = PRF(0x01, w, t) (w the walk index), and the walk moves to
  x g^e with (a, b) to (a + e, b). A restart clears the detector (previous key, anchor). The restarted state is tested for a DP (chained).
  g^e is a fixed-base comb: 32 multiplications over the 8160 precomputed elements g^(v 2^(8i)), v = 1..255.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle; a new key
  equal to the walk's anchor key (recorded every 64 walk steps) inside a window of 64 steps is a cycle of length at most 64. On detection the
  walk walks the cycle once (every point is DP-tested first), takes the least-key state s of the cycle, and escapes to the canonical
  representative of s times g^e with e = PRF(0x02, K_s), exponents (A_s + e, B_s). The escape depends only on the members of the cycle, so
  walks entering the same cycle escape to the same element (Lemma R). The escape does not reset the trail age and does not advance the counter.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of the
  Montgomery integer; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index.
- Table and solving: the records are stored in the bucket table of Section 4; each bucket is sorted by fingerprint. Records with equal
  fingerprints are compared pairwise. A pair from one class with B1 != B2 gives y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is
  discarded; no trail is replayed. Every candidate is verified by g^y = z before output.
- Work counter and cutoff: each operation adds its charge (Section 3). The walks run while W < W_MAX - reserve. The table stops the
  procedure when a bucket or the table is full.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state, including after restarts and escapes. The kernel's
dump, DP records and counts agree with the replica (`measurement/replica.py`) bit for bit at production density (DP_BITS 60) and at
small densities with chained restarts and trail caps (DP_BITS 2 and 4; trail cap 1024 in test builds), at 512 and 4096 rounds
(`measurement/RESULT.json`, `measurement/test_summary.txt`). The verification matrix has 33 cases and 30 pass. In all 33 cases the dump,
the DP records, the event counts and the identity checks pass. The 3 failing cases (cap1024_p1024_dp2_s1, cap1024_p1024_dp2_s2,
ra16_cap1024_dp2_s1) fail only the check "cap restarts occur". That check cannot fire at DP_BITS 2: every walk restarts at a
distinguished point every few steps, and each DP restart resets the trail age, so a trail reaches 1024 steps with probability
(3/4)^1024. Cap behaviour is shown by the trapped-trail builds (123 and 114 cap restarts at RA 16, DP 2, with TEST_ESCAPE_E0) and by
256 cap restarts at DP 60 (6 of them in trails with an escape). The 33-case matrix was run on the pre-trim source and again on the trimmed code in the package (`measurement/matrix_summary.txt`): 30 of 33 pass in both, with identical kernel counts and checks per case. The packaged
`measurement/step.c` and `measurement/replica.py` are trimmed copies (comments changed; the sortbench mode removed, per `RESULT.json`).
`reproduce.sh` was run end to end (exit 0; the forced-add per-step values reproduce to 0.02 instructions). The Lemma R replay covers 768 pairs (512 at RA 16, 256 at RA 1024) with zero mismatches.

## 2. Success probability, failure terms and the success chain

Under H1 the walks need T = 2^125.888062 fresh collision steps for tau = T^2 V_low / (2N), where V_low = 0.7996551059 is the non-degenerate collision rate. V_low is an assumption: V >= V_low = 0.7996551059 is assumed, not derived. The value is the output of the withdrawn §9.3 expression 1 − s² − (1 − s)²/r_a at s = 458/1024, kept as an assumed conservative lower bound (Section 5.2 gives the toy evidence and its limits). The square share enters through its expectation (H9, Section 2), and the charged Hoeffding deviation on the key-level share is 6.96e-38. A collision is detected with probability at least 1 − e^-40 per pair, and the trail cap and detection
factors give f = 1. The success model is 1 − exp(−f tau), tau_eff = 0.567010.

Success chain (each step names its term). (i) The first M main steps execute before the cutoff: the charges of Section 3 are within W_MAX.
(ii) The closure events among them are at most K p M (Markov, failure 1/K = 0.02). (iii) The revisit steps are at most 127 times (ii) plus
the join suffixes and the own-trail revisit steps (at most J_MAX CAP and M CAP^2 / N), so the fresh steps are at least (1 − w) M with
w = 127 K p + small, and M is chosen so that this is at least T + lag. (iv) The square share enters through its expectation (H9), and the
collision rate is taken at the assumed V_low (Section 2). (v) The distinguished-point records are at most X + J
(Section 4), overflow has probability at most the DP terms below. (vi) Given T fresh steps, a collision occurs between walks with probability
1 − exp(−f tau_eff) under H1. Two walks at the same class stay synchronised by Lemma R, so a join leads to a DP collision at the first DP of
the shared sequence. (vii) A collision is found when both walks reach a DP; the residual is geometric with mean 2^60, and it fits the allowance
lag = WALKS * 40 * 2^60 except with probability e^-40 per pair. Each step has its failure term; their sum is the failure total.

P = 0.432781 before delta and 0.41 after it. The success target 0.41 is met because the model is set to 0.41 + delta (accounting.py,
`fresh_steps_needed`). If the waste were not budgeted, P would fall to 0.2283.

Classification. (a) If a closed window contains a square, its exponent coefficient carries a factor 2^q with q >= 1 that is not a 12th root of
unity (below), so the coefficient is nonzero; such a closure is an event of probability 1/r over the coins (H5). (b) If the window is all
additive, the coefficient of adder h is sum_k n_{h,k} lambda^(-k) in F_r. For windows of length L < 2^63 this vanishes in F_r only if
sum_k n_{h,k} zeta^(-k) = 0 in Z[zeta_12] (a nonzero element of the prime above r has norm divisible by r and |N| <= L^4 < 2^252 < r). Proposition G
(below) gives the gon decomposition of each adder's rotation multiset. The 2-cycle is formal: T[h][k + 6] = T[h][k]^(-1), so the walk returns to
the exact element and exponents; these are the closures the detector escapes.

**Proposition G.** Let n in Z_{>=0}^12 satisfy sum_k n_k zeta^k = 0 for a primitive 12th root of unity zeta. Then n is a nonnegative integer
combination of the indicator vectors of the 2-gons {k, k+6} and the 3-gons {k, k+4, k+8} (indices mod 12).

*Proof.* (i) Let Phi(x) = x^4 - x^2 + 1. It vanishes at zeta, since zeta^6 = -1 gives 0 = (zeta^2 + 1)(zeta^4 - zeta^2 + 1) with zeta^2 != -1.
It is irreducible over Q: a factorisation over Z (Gauss's lemma, Phi monic) has no linear factor, as Phi(+-1) = 1; a quadratic factorisation
(x^2 + ax + b)(x^2 - ax + c) forces b + c - a^2 = -1, a(c - b) = 0 and bc = 1, so b = c = 1 with a^2 = 3, or b = c = -1 with a^2 = -1, both
impossible. So Phi is the minimal polynomial of zeta, and sum_k n_k zeta^k = 0 iff the remainder of f(x) = sum_k n_k x^k modulo Phi is zero.

(ii) Modulo Phi, x^4 = x^2 - 1, x^5 = x^3 - x, x^6 = -1, and x^(k+6) = -x^k. The coordinates of x^k in the basis 1, x, x^2, x^3 are
v_0 = (1,0,0,0), v_1 = (0,1,0,0), v_2 = (0,0,1,0), v_3 = (0,0,0,1), v_4 = (-1,0,1,0), v_5 = (0,-1,0,1), v_6 = (-1,0,0,0), v_7 = (0,-1,0,0),
v_8 = (0,0,-1,0), v_9 = (0,0,0,-1), v_10 = (1,0,-1,0), v_11 = (0,1,0,-1). The remainder vanishes iff

  (C1) n_0 - n_4 - n_6 + n_10 = 0,   (C2) n_1 - n_5 - n_7 + n_11 = 0,
  (C3) n_2 + n_4 - n_8 - n_10 = 0,   (C4) n_3 + n_5 - n_9 - n_11 = 0.

(iii) Put p_j = n_j - n_(j+6) for j = 0..5. Then (C1)-(C4) read p_0 = p_4, p_1 = p_5, p_2 = -p_4, p_3 = -p_5. With a = p_4 and b = p_5,
the signed antipodal differences are p = (a, b, -a, -b, a, b).

(iv) Induction on the total weight. The zero vector is the empty combination. Let n != 0 be a nonnegative solution. If some pair has
n_j > 0 and n_(j+6) > 0, the 2-gon {j, j+6} lies in the support; subtracting its indicator keeps n nonnegative and a solution. Otherwise
n_j n_(j+6) = 0 for every j, so p_j = n_j or p_j = -n_(j+6), and n != 0 forces (a, b) != (0, 0). If a > 0, then n_0 = n_4 = n_8 = a (from
p_0 = a, p_4 = a, p_2 = -a), so {0, 4, 8} lies in the support. If a < 0, then n_2 = n_6 = n_10 = -a and {2, 6, 10} lies in the support. If a = 0
and b > 0, then {1, 5, 9} does; if a = 0 and b < 0, then {3, 7, 11} does. The 3-gon indicator is a solution, since zeta^r (1 + zeta^4 + zeta^8) = 0.
Subtracting it keeps n nonnegative and a solution. Induction on the weight completes the proof. ∎

Per adder, the closed window's rotation multiset is a nonnegative solution, so Proposition G gives its decomposition, and `accounting.py` (union_log2,
lengths 3 to 400, with ADDITIVE_CELLS = 12 r_a) counts the gon decompositions. For windows longer than 400 the probability that every step is
additive, (1 − s)^401/s, is used.

Order of 2 and the smallest q. r − 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 * 254760293^2
(checked by `cscore/number_theory.py`, output in `evidence/number_theory_check.out`). The order of 2 modulo r is o = 3277242198445386904965483781761622864855659531282977363912728668746161324032, 251 bits, certified
by 2^o = 1 and 2^(o/p) != 1 for each prime divisor p of o. The smallest q >= 1 with 2^q a 12th root of unity is q0 = o/gcd(o, 12), a 248-bit
number (gcd(o, 12) = 12), so 2^q is not a 12th root of unity for 1 <= q < 2^247.

Failure terms (each a probability over the run; computed in `accounting.py`, function `build`):
| term | value | source |
| --- | ---: | --- |
| formal windows, lengths 65 to 400 (per step) | 4.50e-105 | union bound over gon decompositions (Proposition G, Section 2) |
| all-additive windows longer than 400 (per step) | 1.25e-103 | (1 - s)^401 / s |
| non-formal closures (per step) | 1.69e-56 | 12 CAP / r candidates, each closing with probability 1/r (heuristic H5) |
| trap term over the run (formal and non-formal) | 1.86e-18 | M times the per-step total, including the trap revisit steps |
| fingerprint Markov (metered matches) | 8.67e-19 | 2^-60 |
| square fraction, H9 (Hoeffding on the key-level share) | 8.67e-19 | delta_N = 2^-60 (Section 2) |
| closure events (Markov, reserve K = 50) | 2.00e-02 | 1/K; the reserve K p M is charged |
| non-formal closures of length at most 64 | 1.61e-36 | expectation 12 * 64 * M / r (Markov) |
| DP table, fresh flags above 2 mu | 2^(-4.61e+19) | exp(-mu/3), Chernoff on the M + R new states (Section 4) |
| DP table, joins above J_MAX = 1000 | 2.78e-03 | E[J]/J_MAX (Markov), Section 4 |
| detection delay beyond 40 mean waits per walk | 4.25e-18 | e^-40 per pair (Section 2) |
| own-trail revisits (trapped trails) | 1.86e-18 | M CAP / N (H1, Section 2) |
| DP table, bucket overflow | 2^(-195463482.75) | Chernoff per bucket, union over 2^37 buckets (Section 4) |
| total failure kept | 2.28e-02 | sum of the kept terms |

**Lemma R (replay invariance).** Let two walks stand at the same class, with any detector states (previous keys, anchor phases), any trail
ages and any exponent offset. (a) Each step is a function of the current element (the branch and the index are functions of the class), so the two
element sequences agree until a detection, a DP or a cap in one walk. (b) A detection closes a cycle C of length at most 64 whose points were all
produced by the detecting walk, and no point of C is a DP (the DP test precedes the detector). (c) The other walk, on the same sequence, enters C
at the same point and detects it within 127 steps of entry, with the same least-key state s; both escape to the same element c_s g^e, e a function
of K_s only (PRF). (d) If a DP occurs on the shared sequence before any detection or cap, both walks produce and record it: the same element, with
the same exponent offset. The exponent offset is preserved: the escape adds the same e to both, and squares double both. (e) A cap in one walk and
not the other separates the pair; this is a restart and is counted in R (Section 3). Lemma R is proved for the element sequence; the replica test
(`measurement/replica.py`) checks (a) to (d) on random starts with different phases, previous keys and trail ages.

Premises (disclosed; each is a stated assumption, not a proof).
- **H1 (random map).** Each landing of a walk (a main step or a restart) is uniform on the classes not yet visited by it, given the history,
  except for the closure structure of Section 6 (detected cycles and formal closures). The collision statistic tau uses this model with the assumed rate V_low (Section 2).
  Joins and the trapped-trail bound use this landing rule.
- **H5 (non-formal closures).** A closure with a nonzero coefficient has probability 1/r over the coins. The non-formal trap term rests on it.
- **H7 (uniform new states).** Given the history, the key of each fresh state is uniform: its branch and index bits are uniform, its DP flag has
  probability 2^-60. The closure bound is an expectation bound (Markov) over the process coupled with independent keys until the first revisit,
  so it needs only the unconditional per-step rate p.
- **H8 (PRF restarts and escapes).** The restart and escape exponents are outputs of SHA-256 under a random 256-bit key. Distinct messages give
  independent uniform exponents on Z_r, independent of the history. So every restart or escape target is uniform on G given the history and lands
  on a visited class with probability at most (visited)/N. This is the standard PRF premise. It is used for joins, duplicates and the DP overflow
  bounds. Nothing in the toy or the kernel tests it beyond the no-repeat and landing-rate checks of Section 5.
- **H9 (square share, expectation form).** The branch bit of a class is b(x) = [(mix(K_lo(x) XOR kappa_b) >> (64 - adder_bits)) < S]
  with a uniform key kappa_b; the words of distinct classes are independent uniform 64-bit words (random oracle). (a) A fresh step, taken
  from a class the walk has never visited, has a square decision that is Bernoulli(s), s = S/r_a exactly, conditional on the history, so the
  expected number of square decisions among fresh steps is s times their expected number. (b) The useful-collision count among fresh pairs
  is taken at the assumed rate V_low (Section 2), not derived from the share. (c) The success model is the H1 Poisson approximation with that
  count. The key-level fraction f_N is an average of N independent bits over the random key, so by Hoeffding |f_N - s| <= sqrt(ln(2 / delta_N) /
  (2N)) with probability at least 1 - delta_N, delta_N = 2^-60 (a term of delta); no Azuma bound, visited-set correction or margin is used.
  The expectation is measured with a fresh kappa_b per trial (Section 5.1); a fixed-key run is one realisation and not the protocol.
- **Own-trail revisits (replaces v10 H10; from H1).** Every escape target is a PRF output, so distinct cycles give independent uniform
  targets and no table element or fixed power is reused. The only non-uniform landings are the short formal closures, each detected and
  counted as an event. Any other landing on a class already visited by the trail is a revisit with probability at most a/N (a = trail age,
  at most CAP). A loop, self-trap or composite, must contain such a revisit. The expected number of own-trail revisits over the run is at most
  M CAP / N, about 1.9e-18 (a term of delta). The table-reuse loop of round-10 g1 needs a fixed element used by two escapes; v11 has none.
  The toy trap test (Section 5) reports the escape-loop rate at toy size, where the 1024 adders are not uniform enough on 76880 classes.

## 3. Computation and the work cutoff

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), charged 1 as the rules require; the organizer's reference count (REF_FUNCS)
is 22,990 instructions. Other work is measured in instructions and divided by C = 20,000.

Measurement (`measurement/RESULT.json`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0, `-O3
-fno-tree-vectorize -mgeneral-regs-only`), the pinned blst commit f7a88067, Callgrind Ir, forced inputs at production DP_BITS 60.

| forced step type (production DP_BITS 60, zero distinguished points) | instructions | notes |
| --- | ---: | --- |
| all-adder steps, seed 1 | 26,747.453 | no detected cycle |
| all-adder steps, seed 2 | 26,797.438 | one detected cycle (one escape, 32 comb multiplications); charged as A (the maximum) |
| all-square steps, both seeds | 13,535.453 | no detected cycles; charged as Q |
| restart alone (mode restartbench) | 744,445.844 | PRF, comb (31.87 multiplications on average), a + e, key, clear; charged as 32 multiplications plus the PRF |
| restart PRF alone (mode prfbench) | 7,851 | two SHA-256 blocks and the reduction modulo r |
| store of one record alone (mode storebench) | 1,496.0 | fingerprint and two canonical exponents |
| bitonic compare-exchange (sortnet, 32-bit indices, 2^12 slots) | 108.0030 (108.0015 at 2^13) | branch-free mask and 80-byte conditional swap; charged (the larger) |
| compare-exchange (sortnet64, 64-bit indices, 2^21 and 2^22 slots) | 108.0000054 and 108.0000027 | same network with 64-bit indices |
| bitonic scan, per record (32-bit; 64-bit at 2^21) | 13.0 (13.0000005) | equal-neighbour test, branch-free; charged |
| scatter of one record into its bucket (sortnet64, excluding the cursor reset) | 24.0 | charged |
| reference operation (organizer, REF_FUNCS) | 22,990 | charged 1 unit; blst_fp12_mul measured 22,997 |
| cyclotomic square | 9,803 | measured; not used in the charges |

Charges for the first M main steps:
- Main steps: every step is charged the adder cost 1.191372 units (the larger of the two step charges).
- Closure events: reserve K p M = 2^117.579 events (Markov, failure 1/K = 0.02), each charged 66 adder operations and one restart (the escape,
  38.516451 units for the comb restart and the PRF; 0.264 units per main step in expectation). The revisit steps of these events
  (at most 127 each) are main steps and are counted in M through the reserve w.
- Square decisions (H9): expectation under the random oracle; the charged Hoeffding deviation 6.96e-38 on the key-level share, failure 2^-60 (Section 2).
- Joins: reserve J_MAX = 1000 (Markov, failure E[J]/J_MAX). E[J] <= (M + R)^2/(2N) with R the restarts (escapes, DP restarts, cap resets;
  R <= M), so E[J] <= M^2/N is charged. Suffixes at most CAP each (in w).
- Restarts: each DP record, each escape and each trail-cap reset costs one restart: 32 comb multiplications charged at the full adder cost each, plus the PRF (measured 7,851 instructions, `prfbench`). The restart is 38.516451 units (32 adder charges plus the PRF instructions over C).
- Fingerprint matches: metered (one exponentiation each); the metered total is bounded by 2^60 times the expected count (Markov, failure 2^-60).
- Non-formal closures of length at most 64: expectation 12 * 64 M / r (Markov).
- Comb table: 8160 multiplications at the adder cost (preprocessing).

Components of W_MAX (log2 units; the score is the cutoff 126.915106, exactly 63457553/500000):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 126.626579 | M = 2^126.373955 main steps, each charged A = 1.191372 units |
| closure events | 124.451256 | reserve K p M = 2^117.579 events (K = 50), each 66 adder operations and one restart |
| adder table | 19.582053 | 1024 adders (preprocessing) |
| walk starts | 29.573647 | 2^20 walk starts (preprocessing) |
| comb table | 13.246977 | 8160 elements at the adder cost (preprocessing) |
| restarts at the trail cap | 65.641358 | M / 2^66 cap resets, one restart each |
| distinguished points | 72.889938 | capacity 2^67.38 records (store, restart, scatter, transfer, overflow) and 2^68 table slots (sort, scan) |
| fingerprint verification | 74.328047 | Markov charge at 2^60 times the expected false pairs (2^5.75 at capacity), one exponentiation each |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

The reserves enlarge M by the factor 1/(1 − w): the revisit reserve w = 127 K p = 0.285945 is the largest term in the increase.

## 4. Distinguished-point table, pathwise capacity, scratch memory and peak memory

Records. Each record comes from a DP: a fresh DP flag (probability 2^-60 per new state, H7 and H8) or a duplicate at a join. New states are
the M main steps and the R restarts (DP restarts, escapes, cap resets); R <= M. The records are at most X + J pathwise, where X is the number
of fresh flags over the M + R new states and J the number of joins. Overflow needs X > 2 mu (mu = (M + R) 2^-60, failure exp(-mu/3)) or
J > J_MAX (failure E[J]/J_MAX). The capacity is D = ceil(2 mu) + J_MAX = 2^67.3772 records.

Buckets. The table is 2^37 buckets of 2^31 slots each, 64-bit indices (`measurement/sortnet/sortnet64.c`), partitioned by the top 37 bits
of the fingerprint (bucket = fp[1] >> 27). The mean occupancy is D / 2^37 = 2^30.3772, against capacity 2^31; the Chernoff bound per
bucket and the union over the buckets give a failure of 2^-195463482.75 (negligible). Each bucket is sorted by a bitonic network
(n k (k+1)/4 compare-exchanges with n = 2^31, k = 31, i.e. 248 per slot); its scratch is one 80-byte temporary. Equal fingerprints share a
bucket, so matching needs no global sort.

Costs (instructions, per record or per slot): the store (1,496 instructions, measured alone), one restart (770,329 instructions, the comb and the
PRF as charged), the bucket scatter (24.0, measured), transfer and index widening (100,000, an explicit assumption), overflow (10); per slot:
the bitonic compare-exchange (108.0030, the 32-bit measurement, larger than the 64-bit values 108.0000054 and 108.0000027), the scan (13.0),
and one cursor reset per bucket (1). The charge is over the capacity records and all 2^68 slots (padding included). The table costs
about 2^-54.03 of W_MAX.

Peak memory. The table occupies 2^68 slots of 80 bytes (2^74.321929 bytes in total, spread over the buckets; each bucket 2^37.32 bytes, below the 64-bit
address space), plus the adder table (12 x 1024 elements, about 7.8 MiB), the comb table (8160 elements, about 4.7 MiB), the walk states
(2^20 x 640 bytes) and detectors. The peak is not realisable on the hardware used here and is reported, not scored.

## 5. Evidence: toy-scale runs under the v11 walk

Toy group: prime order ell = 922561 (q = 31), N = 76880 classes, for the success, H9, degenerate-pair and loop runs; the N-scaling runs (5.6) use N = 6936, 76880 and 284746. Walk: r_a = 1024, S = 458 (s = 0.447265625), the v11 restart rule, a branch key per trial, no trail cap unless stated, each trial stopped at its first verified output. Raw results: `evidence/toy-v11c.json`, `evidence/toy-v11d.json`. Extracts with values copied from the raw runs (per-trial lists and logs omitted): `evidence/toy-v11b-extract.json`, `evidence/toy-h9.json`, `evidence/toy-loops.json`. The toy shows the H9 expectation (5.1), a verified rate consistent with the 0.41 model (5.4), a DP 6 gap that is detection overhead (5.3), a non-degenerate rate that supports V_low at toy size only (5.2, 5.6), and trap counts within noise of the window model (5.5). The production collision rate and production trap rate are not measured.

### 5.1 Branch bit under per-attack keys (H9)

Protocol: each attack draws its own branch key kappa_b; the branch bit is the top adder_bits bits of mix(K_lo XOR kappa_b) compared with S. We report the mean over attacks with the between-attack standard error. The bit does not depend on DP; the walks do (DP 2 at r_a 32, DP 6 at r_a 1024).

| quantity (1000 attacks per r_a) | r_a 32 (S 14, s 0.4375) | r_a 1024 (S 458, s 0.447266) |
| --- | --- | --- |
| unfiltered first-visit bit (each class once) | 0.43937, SE 0.00100, z +1.87 | 0.44660, SE 0.00091, z -0.73 |
| fresh step taken from a first-visited state (proof definition) | 0.43903, SE 0.00115, z +1.34 | 0.44660, SE 0.00089, z -0.75 |
| fresh step whose landing is a first visit | 0.43966, SE 0.00114, z +1.89 | 0.44704, SE 0.00089, z -0.25 |

Every row is within 3 standard errors of s. The unfiltered test counts every class with no exclusions. The class-level square fraction averaged over the 1000 keys is 0.43756 (r_a 32) and 0.44722 (r_a 1024). Fixed-key diagnostic: the toy's default key 0x5851F42D4C957F2D gives class fractions 0.43594 (minus s -0.0016, z -0.87) and 0.44547 (minus s -0.0018, z -1.00), and per-attack shares 0.43306 (SE 0.00118) and 0.44375 (SE 0.00088). A fixed-key run is one realisation of the key and not the protocol. Sources: `evidence/toy-h9.json`.

### 5.2 Non-degenerate collision rate (assumption)

V_low = 0.7996551059 is the output of the withdrawn §9.3 expression at s = 458/1024, kept as an assumed lower bound; it is not derived, and no closed form is claimed. Every degenerate fresh pair (B1 = B2) in every run was a same-walk 2-cycle: 1, 2 and 4 pairs at the sizes of 5.6; 26 of 26 (mixed, DP 6, n = 3000); 9 of 9 (mixed, DP 3); 1 of 1 (mixed, DP 2); 99 of 99 (add-only, DP 6, n = 3000); 22 of 22 (add-only, DP 3); 12 of 12 (add-only, DP 2). No square-square fresh pair occurred. At N = 76880 the non-degenerate rates are 0.998 (mixed DP 2), 0.988 (mixed DP 3), 0.9935 (mixed DP 6), 0.980 (add-only DP 2), 0.974 (add-only DP 3) and 0.977 (add-only DP 6); the add-only runs have the larger degenerate share (2.0 to 2.6 percent). This supports V_low at toy size only (Section 6). Source: `evidence/toy-v11b-extract.json`.

### 5.3 Verified-success gap at DP 6: detection overhead, not recovery

At DP 6 (DP spacing 0.23 sqrt N), r_a 1024 mixed, n = 3000 (seed 1001), the verified rate within T = 318.5 steps is 0.183 [0.169, 0.197] and the first-useful rate is 0.408 [0.390, 0.425]. The verified output follows the first useful collision by a mean of 153 steps (median 84, 90th percentile 396). That is detection overhead set by the DP spacing, not a recovery error: replaying every trial's log gives zero reconstruction mismatches, and category iii (an unverified solve) is empty. Of the 79 trials with no verified output by the cutoff of 40 sqrt N + 2000 = 13,090 steps, 48 had no DP point after the collision (category iv) and 31 had DP points without a class reached twice with different B (category ii). At DP 3 (spacing 0.029 sqrt N, n = 1000, seed 1004) the rate is 0.423 [0.393, 0.454] (first-useful 0.464, mean overhead 24 steps). At DP 2 (spacing 0.014 sqrt N, n = 1000, seed 1003) it is 0.425 [0.395, 0.456] (first-useful 0.438, mean overhead 11.5 steps). DP 2 and DP 3 are the production-like regime for the success claim (5.4). Production spacing is 2^60 against a walk near 2^126, a ratio near 2^-66.

### 5.4 Budget sweep at DP 2

Same protocol at DP 2, r_a 1024 mixed, cutoff 13,000 steps, n = 2000 attacks per budget with one seed per budget (1101 for the control tau = 0.5277). The budget is T = sqrt(2 N tau / V_low).

| tau | T (steps) | verified rate [Wilson 95%] | first-useful rate | model 1 - e^-tau |
| --- | ---: | --- | ---: | ---: |
| 0.5277 (control) | 318.5 | 0.436 [0.414, 0.458] | 0.455 | 0.410 |
| 0.55 | 325.2 | 0.4315 [0.410, 0.453] | 0.4525 | 0.423 |
| 0.60 | 339.7 | 0.468 [0.446, 0.490] | 0.486 | 0.451 |
| 0.65 | 353.5 | 0.4825 [0.461, 0.504] | 0.504 | 0.478 |

The pooled verified rate over the 3000 attacks at tau = 0.5277 (n = 2000 at seed 1101, n = 1000 at seed 1003, `evidence/toy-v11b-extract.json`) is 0.4323 [0.4147, 0.4501]. This is consistent with the 0.41 model. It is not a confirmation of a margin: the bound 0.41 is analytic (delta = 0.0228), and the toy checks the walk at N = 76880, not the production V_low.

### 5.5 Loop regime (trail caps, DP 60 and DP 6)

Trail-cap runs: 60 trials each, group 201, cutoff 20,000 steps per trial; first-visit trap hits against the window model's expectation (`evidence/toy-loops.json`).

| DP | cap 64 | cap 256 | cap 1024 | cap 4096 | uncapped | pooled (z) |
| --- | --- | --- | --- | --- | --- | --- |
| 60 | 0 vs 0.01 | 3 vs 3.2 | 9 vs 9.6 | 12 vs 13.9 | 12 vs 5.4 | 36 vs 32.2 (0.67) |
| 6 | 0 vs 0.21 | 1 vs 1.00 | 5 vs 1.33 | 2 vs 1.39 | 0 vs 0.88 | 8 vs 4.82 (1.45) |

The spec bound (CAP times the trail term) expects 8.1, 18.4 and 27.9 hits at DP 60 caps 256, 1024 and 4096; the observed counts are below it. Single-run excesses are the uncapped DP 60 run (12 against 5.4, Poisson z 2.9) and the cap-1024 DP 6 run (5 against 1.33, z 3.2). Pooled over the five runs of each DP, z is 0.67 (DP 60) and 1.45 (DP 6). Neither pooled excess is significant at this size, and we do not call the single-run excesses noise. The uncapped DP 60 excess (12 against 5.4) is an unexplained limitation of the toy. Production runs are capped, and the capped DP 60 runs (caps 256, 1024 and 4096) match the window model (3 against 3.2, 9 against 9.6, 12 against 13.9). The uncapped regime is not the production regime, so this excess is not a production claim. Undetected loops are counted per window (each escape or restart opens one): 990 windows at cap 256 (DP 60) against 3 hits; cap 64 has none by construction; formal undetected loops are zero in every run. The production trap terms (accounting.json, 1.86e-18) are computed, not measured.

### 5.6 N-scaling of the degenerate fraction

r_a 1024 mixed, DP 2, n = 1000 attacks at each group size, branch key per attack (`evidence/toy-v11d.json`).

| N | sqrt N | fresh pairs | degenerate | fraction [Wilson 95%] | all-collision degenerate fraction |
| --- | ---: | ---: | ---: | --- | --- |
| 6936 | 83.3 | 583 | 1 | 0.0017 [0.0003, 0.0097] | 0.0027 (13/4784) |
| 76880 | 277.3 | 515 | 2 | 0.0039 [0.0011, 0.0140] | 0.0010 (4/4097) |
| 284746 | 533.6 | 502 | 4 | 0.0080 [0.0031, 0.0203] | 0.0012 (5/4088) |

The fresh-pair fraction rises in point estimate (0.0017, 0.0039, 0.0080), but a homogeneity test on the counts gives chi-square 2.46 on 2 degrees of freedom (p about 0.29), so the rise is not established; the all-collision fraction does not rise, and all degenerate pairs are 2-cycles. A 2-cycle may recur at a rate per step that does not depend on N, so its share could grow with N; a factor of 41 in N cannot show production behaviour. In this run the degenerate fraction stays below 0.8 percent up to N = 284746. At N = 76880, across the runs of 5.2, it is 0.2 to 2.6 percent. The production degenerate share, and so production V_low, remains an assumption.

## 6. Limitations

- H1 is extrapolated from toy subgroups with N up to 284746 (bits 20; N-scaling) and N = 76880 (success runs), with 1024 adders, to N about 2^251. The closure structure of Section 6 of the
  walk spec (formal 2-cycles, short formal cycles, undetected long cycles) is handled by the event and trap terms, not by H1.
- H5 is a heuristic: the non-formal closure term depends on the coins being independent of the path, and it does not address adaptively chosen paths.
- H7, H8 (PRF), H9 (decoupled square share) and H1-trap are stated premises. H8 is the PRF premise; it is not tested beyond the no-repeat checks.
  H9 rests on the mix being uniform on the fresh keys; the toy tests this at toy keys only.
- The lemma R checks (`measurement/replica.py`) are at small densities and on the toy; at production size the replay lemma is the proof.
- The toy has a trail cap only where its experiment sets one, and its restart share and dp density differ from production.
- The DP_BITS 0 run is not defined under the chained rule and is not charged; the record store and the restart are measured alone.
- The kernel's dp_exact_match compares the 65,536-slot record ring; at DP_BITS 2 with 4096 rounds the record count exceeds the ring, so the record
  contents are compared on the ring slots and the counts in full.
- The bucket organisation is measured with 64-bit indices at 2^21 and 2^22 slots; the production bucket size (2^31) is not run. The transfer cost
  is an assumption (100,000 instructions per record).
- The PRF cost is measured (7,851 instructions per restart PRF, `prfbench`); the SHA-256 implementation is checked by the test vectors (`measurement/check_vectors.py`).
- The reference operation is charged at the organizer count of this run (22,990); blst_fp12_mul measured 22,997.
- The solve, verify and output charge of 2^29 instructions is an assumed bound (2^14.7 of the units).
- The sort costs (108.003 per 32-bit compare-exchange, 13.0 per scan record) are carried from the v10 sortnet measurement
  (`measurement/sortnet/`) and are not re-measured in v11.
- The Lemma R replay covers 768 pairs (512 at RA 16, 256 at RA 1024), with zero mismatches.
- The restart costs about 31.9 comb multiplications (measured 31.87 per restart) against the bound 32.
- `reproduce.sh` was not executed end to end. Its components were run separately; their results are in `measurement/RESULT.json`,
  `measurement/native-counts-v11.json` and `measurement/sortnet/`.
- The 33-case kernel matrix passes 30 cases (Section 1), on the pre-trim source and on the packaged files (`measurement/matrix_summary.txt`).
- V_low is the assumption that carries the success claim. The toy evidence for it (Sections 5.2 and 5.4) is at N = 76880.
- Preflight (`llm-verifier preflight`) has not been run on this package.

## References

- Lam and Leung, On vanishing sums of roots of unity (2000), for the background of Proposition G (the statement is proved here).
- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py                 # bundled: W_MAX, failure terms, number-theory asserts
python3 -I experiment.py 0.5277 1101 5   # one toy budget run: TAU SEED TRIALS, prints one JSON object (Section 5.4)
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image; restart checks at small densities
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind at zero DP density, two sizes, mul and sqr loops, the DP-every-step
run and the restart), checks the replica at production and small DP densities, runs the native counts, and runs the sort measurements
(`measurement/sortnet/`).

# frobenius-mixed-rho-v5: Frobenius-class mixed add/square rho on GT (revision v10)

Claim: for the uniform challenge of `rules.md` (P uniform nonidentity in G1, z uniform in GT, fresh algorithm coins), the procedure below
outputs Q in G2 with e(P, Q) = z with probability at least 0.41 (model 0.432848326897; lower bound 0.410000000001 after the failure terms of
Section 2, which include a 2% reserve loss for closure events), while its total work, counting every operation it performs, is
at most W_MAX = 2^126.840188 target-group-operation equivalents (score 31710047/250000). The first M = 2^126.391121 main steps are executed before
the cutoff except with probability 2.28e-02 (Section 3), under the premises H1, H5 and H7 stated in Section 2. Peak memory is at most
2^74.321929 bytes, including the padding of the bucketed distinguished-point table (Section 4); preprocessing is inside W_MAX (2^29.575768
units); nonuniform advice is 0 bytes. It is a classical generic attack, not a practical break.

**Revision v10.** Round 10 (on revision v9, both tracks; the round-9 g2 review) found one structural recurrence in the restart, and
asked for the join analysis to be re-derived. v9 restarted with an adder chosen by a per-walk counter from the same 1024 adders.
Let c be a canonical state whose additive successor d = (c R_h)^(-1) is a DP. After d is recorded, a restart that picks h again returns
the walk to the class of c (d R_h = c^(-1)); the detector is clear, so the walk reaches d again and stores a duplicate with the same
canonical exponents. The counter does not prevent this. The rate is about (1 - s)/(12 * 1024) = 4.5e-5 per DP, so v9's join bound and
duplicate argument did not hold as written.

(1) Restart recurrence (defect, fixed). v10 restarts by multiplying by a restart element Y_j = g^(f_j) from a separate table of
T_R = 1024 elements (walk-spec section 6). The table exponents f_j are independent of the adder exponents, and preprocessing checks that
no Y_j equals an adder element or its inverse. A restart is therefore not an adder step and cannot undo the previous one: it moves the
walk to the class of d Y_j, which no adder produces. The table index j = mix((w + 1) GAMMA + t) mod T_R is a function of the counter.
(2) Joins over all new states. Every new state, a main step or a restart, lands on a visited class with probability at most
(visited)/N: H1 for main steps and H8 (below) for restarts. Hence E[J] <= (M + D)^2 / (2N), with D the number of restarts
(D <= M 2^-60 + M/2^66), which is at most M^2/N = 2.78 as before. The Markov reserve J_MAX = 1000 with failure E[J]/J_MAX =
2.85e-03 is unchanged, and it is in delta. Each join gives at most one duplicate record, after which the restart separates the walks.
(3) Rare-event reserves, detection delay and table (v9 points (2) to (4), unchanged): the closure reserve K = 50 (failure 0.02), the
detection allowance of 40 mean waits per walk (failure e^-40), and the bucket table of 2^37 buckets of 2^31 slots with
64-bit indices.
(4) Escape and replay (round 9 g1 review). A detected cycle no longer escapes by squaring. It escapes to the canonical representative of
its least-key state times the table element Y_j (keyed by that state's key), so every walk that enters the cycle leaves it at the same
element, whatever its anchor phase or previous key (Lemma R, walk-spec section 6). The DP test precedes the detector. The escape output
is a state of the same walk whose class is fixed by the detected cycle and the table, so it is fresh for H7. The closure rate is argued per
fresh entry: a fresh state's successor is a first revisit with conditional probability at most p. The expected number of closures is at most
p M by Markov over the main steps. No independence across revisits is claimed; the closure reserve K = 50 is unchanged.
(5) Square fraction and domination. The Azuma bound over all M main steps is removed. The square fraction of the fresh steps is a premise (H9),
and its deviation is bounded by Markov on the absolute deviation. New states are the main steps and the restarts (DP restarts, escapes and cap
resets), so the fresh-flag mean is (M + R) 2^-60 and E[J] <= (M + R)^2/(2N) <= M^2/N.
(6) Checks added. The replica compares every restart output with every state visited earlier in the run, and checks Lemma R directly: two walks
at the same class with different anchor phases, previous keys and exponent offsets reach the same first DP, and their escapes coincide (Section 1).
The toy experiment measures how often a restart lands on a visited class against the H8 prediction, counts duplicate records with the v10 restart
and with the v9 adder restart, and runs at production-like parameters (r_a = 1024, S = 458) to test the collision law (Section 5).
(7) The toy experiment runs the v10 walk at production-like parameters (r_a = 1024, S = 458, bits 18, dp_bits 6). Its success at tau = 0.5277 is
0.406 to 0.436 (against 0.41), within the Wilson width, and the earlier shift of 0.03 to 0.06 at r_a = 32 follows the closure rate, which scales
with 1/r_a. The toy also shows two problems that the proof now states as premises: escape loops (H10) and a fresh square share below s (H9).

Score: v9 126.823186 to v10 126.840188 (+0.017002 bits). Of the increase, 0.0168 bits come from the square-fraction margin of H9 (V is
evaluated at s + 0.02, section 3), about 0.0002 from the measured adder cost (26,740.129 instructions for seed 1), and nothing measurable from
the restart table (18.57 log2 units of preprocessing) or from the trapped-trail term of H10 (2.7e-19). The failure total is 2.28e-02.

## 1. Procedure

Reduction. g = e(P, Q0) generates GT, where Q0 is the fixed G2 generator; z = g^y gives Q = [y]Q0 with e(P, Q) = z. Input a nonidentity P in G1 and a target z in GT. Return the identity if z = 1 (Q = O is valid).
Otherwise run the walks until a verified y is found or the work counter reaches W_MAX, then output Q = [y]Q0, verified by the pairing
equation e(P, Q) = z (one pairing, charged). Exceptional inputs (P = O) are rejected by the rules.

Walk (`walk-spec.md`, revision v10, implemented by the blst kernel `measurement/step.c` and its replica `measurement/replica.py`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a. If h < S, square x and
  double (a, b). Otherwise multiply x by T[h][k] = sigma^(-k)(R_h), where R_h = g^(c_h) z^(d_h), and update (a, b) by lambda^(-k) (c_h, d_h).
- Trails: a trail ends at its first distinguished point. There the walk stores its record and restarts: with t the walk's restart
  counter (incremented before use) and w the walk's index, j = mix((w + 1) * 0x9E3779B97F4A7C15 + t) mod T_R with T_R = 1024, and the
  walk multiplies x by the restart element Y_j = g^(f_j) and updates (a, b) to (a + f_j, b), where f_j is the restart table exponent
  (Section 2, H8). A restart is not an adder step and not a walk step. The restarted state is itself tested for a DP (chained restarts).
  A restart clears the walk's detector (previous key, anchor). A trail that reaches 2^66 steps restarts the same way without a record.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle; a new key
  equal to the walk's anchor key (recorded every 64 walk steps; restarts are not steps) inside a window of 64 steps is a cycle of length at
  most 64. On detection the walk walks the cycle once (every point is DP-tested first), takes the least-key state s of the cycle,
  and escapes to the canonical representative of s times the restart element Y_j, j = mix(low word of K_s) mod 1024, with exponents
  (A_s + f_j, B_s). The escape depends only on the members of the cycle, so walks entering the same cycle escape to the same element
  (Lemma R, walk-spec section 6); its target is fresh (table element), so the trail that follows is fresh by construction.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of the
  Montgomery integer; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index.
- Table and solving: the records are stored in the bucket table of Section 4; each bucket is sorted by fingerprint. Records with equal
  fingerprints are compared pairwise. A pair from one class with B1 != B2 gives y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is
  discarded; no trail is replayed. Every candidate is verified by g^y = z before output.
- Work counter and cutoff: each operation adds its charge (Section 3). The walks run while W < W_MAX - reserve. The table stops the
  procedure when a bucket or the table is full.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state, including after restarts; the replica reproduces the kernel's
dump, DP records and counts bit for bit at production density (DP_BITS 60) and at small densities that exercise chained restarts (DP_BITS 2,
3, 4) and the trail cap (CAP overridden to 50 and 24), at 512 and 4096 rounds (`measurement/RESULT.json`), and every restart output is compared with every state visited earlier in the run
(no repeat at any density checked); the bitonic buckets are
checked for sorting (`measurement/sortnet/check_bitonic.py`, `sortnet64.c`).

## 2. Success probability, failure terms and the success chain

Under H1 the walks need T = 2^125.905108 fresh collision steps for tau = T^2 V / (2N), with V = 0.799413681 at the nominal square fraction (the
square fraction of the new states is within 10^-6 of s, Section 3). A collision is detected with probability at least 1 − e^-40 per pair,
and the trail cap and detection factors give f = 1.000000000000. The success model is 1 − exp(−f tau), tau_eff = 0.567128510.

Success chain (each step names its term). (i) The first M main steps execute before the cutoff: the charges of Section 3 are within W_MAX.
(ii) The closure events among them are at most K p M (Markov, failure 1/K = 0.02). (iii) The revisit steps are at most 127 times (ii) plus the
join suffixes (at most J_MAX CAP), so the fresh steps are at least (1 − w) M with w = 127 K p (+ 2^-50), and M is chosen so that this is at
least T + lag. (iv) The distinguished-point records are at most X + J (Section 4), overflow has probability at most the DP terms below.
(v) Given T fresh steps, a collision occurs between walks with probability 1 − exp(−f tau_eff) under H1. The two walks of a join stay
synchronised by Lemma R (same element sequence; detections land on the same escape element), so a join leads to a DP collision at the
first DP of the shared sequence. (vi) A collision is found when both
walks reach a DP; the residual is geometric with mean 2^60, and it fits the allowance lag = WALKS * 40 * 2^60 except with probability e^-40 per
pair. A join is itself a collision: the joined walk's duplicate record is paired with the other walk's record, and the two B values differ except
with probability 1/r (degenerate pairs are negligible). Each step has its failure term; their sum is the failure total.

P = 0.432848326897 before delta and 0.410000000001 after it. The success target 0.41 is met because the model is set to 0.41 + delta (accounting.py,
`fresh_steps_needed`). If the waste were not budgeted, P would fall to 0.228227.

Classification. (a) If a closed window contains a square, its exponent coefficient carries a factor 2^q with q >= 1 that is not a 12th root of
unity (below), so the coefficient is nonzero; such a closure is an event of probability 1/r over the coins (H5). (b) If the window is all
additive, the coefficient of adder h is sum_k n_{h,k} lambda^(-k) in F_r. For windows of length L < 2^63 this vanishes in F_r only if
sum_k n_{h,k} zeta^(-k) = 0 in Z[zeta_12] (a nonzero element of the prime above r has norm divisible by r and |N| <= L^4 < 2^252 < r). Proposition G
(below) gives the gon decomposition of each adder's rotation multiset.

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
lengths 65 to 400) counts the gon decompositions. For windows longer than 400 the probability that every step is additive, (1 − s)^401/s, is used.

Order of 2 and the smallest q. r − 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 * 254760293^2
(checked by `cscore/number_theory.py`, output in `evidence/number_theory_check.out`). The order of 2 modulo r is o = 3277242198445386904965483781761622864855659531282977363912728668746161324032, 251 bits, certified
by 2^o = 1 and 2^(o/p) != 1 for each prime divisor p of o. The smallest q >= 1 with 2^q a 12th root of unity is q0 = o/gcd(o, 12), a 248-bit
number (gcd(o, 12) = 12), so 2^q is not a 12th root of unity for 1 <= q < 2^247.

Failure terms (each a probability over the run; computed in `accounting.py`, function `build`):
| term | value | source |
| --- | ---: | --- |
| formal windows, lengths 65 to 400 (per step) | 1.77e-96 | union bound over gon decompositions (Proposition G, Section 2) |
| all-additive windows longer than 400 (per step) | 1.25e-103 | (1 - s)^401 / s |
| non-formal closures (per step) | 1.69e-56 | 12 CAP / r candidates, each closing with probability 1/r (heuristic H5) |
| trap term over the run | 1.88e-18 | M times the per-step total (includes the trap revisit steps) |
| fingerprint Markov (metered matches) | 8.67e-19 | 2^-60 |
| square fraction, H9 (Markov on the absolute deviation at the margin) | 2.35e-18 | sqrt(s(1-s)/M) / 0.02 |
| closure events (Markov, reserve K = 50) | 0.02 | 1/K; the reserve K p M is charged |
| trapped trails (H10) | 2.72e-19 | E[escapes] CAP 64 / N (Section 2) |
| non-formal closures of length <= 64 | 1.6e-36 | expectation 12 * 64 * M / r (Markov) |
| DP table, fresh flags above 2 mu | 2^(-4.66e+19) | exp(-mu/3), Chernoff on the M + R new states (Section 4) |
| DP table, joins above J_MAX = 1000 | 2.85e-03 | E[J]/J_MAX (Markov); E[J] <= (M + R)^2/(2N) <= M^2/N, main steps by H1 and restarts by H8 (Section 4) |
| detection delay beyond 40 mean waits per walk | 4.25e-18 | e^-40 (Section 2) |
| DP table, bucket overflow | 2^(-184681256) | Chernoff per bucket, union over 2^37 buckets (Section 4) |
| total failure kept | 2.28e-02 | sum of the kept terms |

**Lemma R (replay invariance).** Let two walks stand at the same class, with any detector states (previous keys, anchor phases) and any
exponent offset. (a) Each step is a function of the current element, so the two element sequences agree until a detection or a DP in one
walk. (b) A detection closes a cycle C of length at most 64 whose points were all produced by the detecting walk: detection needs a return to
an anchor on C, which takes at least |C| steps after entry. No point of C is a DP, since a DP among them would have restarted the walk first.
(c) The other walk, on the same element sequence, enters C at the same point and detects it within 127 steps of entry. Both find the same
least-key state s and escape to the same element c_s Y_j. (d) If a DP occurs on the shared sequence before any detection, both walks produce it
and record the same element. So the first DP either walk reaches after the shared class is the same element for both, provided the walks are
not separated by a detection that the other has not yet made; in that case the detecting walk and the other meet again at the escape element
(a fresh element), where their sequences coincide from then on. (e) The exponent offset is preserved: every operation adds the same table
exponents to both, and squares double both, so the two records of one DP have offsets that are nonzero (B1 != B2) unless the walks were
identical in exponents, which has probability about 1/r per pair. The replica test (`measurement/replica.py`, lemma R block) checks (a) to
(c) on random starts with different phases and previous keys.

Premises (disclosed; each is a stated assumption, not a proof).
- **H1 (random map).** The next state of a step is visited by another trail with probability at most (visited)/N; tau is the random-map
  collision statistic with rate V. Joins use this bound.
- **H5 (non-formal closures).** A closure with a nonzero coefficient has probability 1/r over the coins. The non-formal trap term rests on it.
- **H7 (uniform new states).** Given the history, the key of each new state is uniform: its adder index and square decision are uniform, its DP
  flag has probability 2^-60, and its square decision is Bernoulli(s). Restart outputs are covered by H8. The closure bound is an expectation
  bound (Markov) over the process coupled with independent keys until the first revisit, so it needs only the unconditional per-step rate p.
- **H8 (restart table).** The restart element Y_j = g^(f_j) has an exponent f_j drawn independently of the adder exponents and of the history.
  Over that random choice, the output x Y_j of a restart is uniform on the classes given x and the history, so it lands on a visited class
  with probability at most (visited)/N. Table reuse does create relations among table exponents: by pigeonhole, some vanishing combinations
  of about 25 to 40 table elements exist. A join through such a relation needs two states whose histories differ by one of these combinations,
  which is generic (probability about 1/r per pair) for states of independent trails; it is a premise, and the replica and toy checks test
  its consequence (no repeats; landing rates against the prediction). The restart index is a deterministic function of the counter (section 6).
- **H9 (square fraction).** The square decisions of the fresh steps have conditional mean s and conditional variance s(1 − s) (H7).
  The collision rate is evaluated at s + 0.02, which is a margin, not a bound: the toy's fresh-step square share deviates from s by
  −0.018 at r_a = 32 (z = −14) and by −0.005 at r_a = 1024 (z = −6), and by +0.002 in the upward direction (z = +2); V decreases in s, so
  the upward direction is the one that matters, and 0.02 exceeds the observed upward deviation by a factor of about eight. The revisit
  steps are not fresh and are charged through the reserve w. This replaces the Azuma bound over all M steps of revision v9.
- **H10 (trapped trails).** A trail started at an escape target E re-enters the cycle C of E (or any cycle of length at most 64 that
  its escape left) before its first DP or its trail cap with probability at most CAP * 64 * d / N, d = 1 (H1). Then at most
  E[E] * CAP * 64 / N ≈ 2.7e-19 (in delta) is the chance of any trap. A trapped trail costs at most one trail cap (2^66 steps). The
  toy shows traps at N = 76880 (Section 5); production N is ≈ 2^251, and the bound is an extrapolation under H1.

## 3. Computation and the work cutoff

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), charged 1 as the rules require; the organizer's reference count (REF_FUNCS)
is 22,990 instructions. Other work is measured in instructions and divided by C = 20,000.

Measurement (`measurement/RESULT.json`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0, `-O3
-fno-tree-vectorize -mgeneral-regs-only`), the pinned blst commit f7a88067, Callgrind Ir, forced inputs at production DP_BITS 60.

| forced step type (production DP_BITS 60, zero distinguished points) | instructions | notes |
| --- | ---: | --- |
| all-adder steps, seed 1 | 26,740.129 | v10 kernel; includes 4 detected cycles (the maximum over seeds, charged as A) |
| all-adder steps, seed 2 | 26,728.376 | v10 kernel; includes 1 detected cycle |
| all-square steps, both seeds | 13,512.422 | no detected cycles (charged as Q) |
| restart alone (mode restartbench) | 26,432.0 | counter, mix, table index, Y_j multiplication, a + f_j, key derivation, state clear; charged one adder step |
| store of one record alone (mode storebench) | 1,496.0 | fingerprint and two canonical exponents; charged as the store |
| bitonic compare-exchange (sortnet, 32-bit indices, 2^12 slots) | 108.003 (108.001 at 2^13) | branch-free mask and 80-byte conditional swap; charged (the larger) |
| compare-exchange (sortnet64, 64-bit indices, 2^21 and 2^22 slots) | 108.0000054 and 108.0000027 | same network with 64-bit indices |
| bitonic scan, per slot | 13.0 | equal-neighbour test, branch-free; charged |
| scatter of one record into its bucket (sortnet64) | 24.0 | charged |
| ring exchange (kernel mode sortbench) | 48.0 | micro-benchmark of one comparison and one exchange; not charged |
| reference operation (organizer, REF_FUNCS) | 22,993.1 this run (22,990 in the v9 run) | charged 22,990 (the lower count, conservative); blst_fp12_mul measured 22,998 |
| cyclotomic square | 9,804 | charged 9,804 / C |

Charges for the first M main steps:
- Main steps: every step is charged the adder cost 1.18851 units (the larger of the two step charges).
- Closure events: reserve K p M = 2^117.596544 events (Markov, failure 1/K = 0.02), each charged 66 adder operations (0.1767 units per main step
  in expectation). The revisit steps of these events (at most 127 each) are main steps and are counted in M through the reserve w.
- Square decisions (premise H9; no concentration over all M steps): the fresh decisions are a martingale difference sequence with
  conditional variance s(1 - s). By Markov on the absolute deviation, P(|F_s − s| > 10^-6) <= (s(1 − s)/F)^(1/2) / 10^-6, which is
  2.35e-18 at the production M (the margin of 0.02 in V dominates; see H9). The revisit steps are not fresh; they are charged through the reserve w.
- Joins: reserve J_MAX = 1000 (Markov, failure E[J]/J_MAX). E[J] <= (M + R)^2/(2N) with R the restarts (escapes, DP restarts, cap
  resets; R <= M), so E[J] <= M^2/N is charged (the bound is conservative by about a factor 2). Suffixes at most CAP each (in w).
- Restarts: each DP record, each escape and each trail-cap reset costs one restart, charged one adder step (A_UNITS; the measured
  restart is 26,432 instructions, including the multiplication by Y_j). The escapes are the closure events (at most K p M).
- Fingerprint matches: metered (one exponentiation each); the metered total is bounded by 2^60 times the expected count (Markov, failure 2^-60).
- Non-formal closures of length at most 64: expectation 12 * 64 M / r (Markov).

Components of W_MAX (log2 units; the score is the rounded-up cutoff 126.840188):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 126.640271 | main steps M = 2^126.391121, each charged the adder cost 1.18851 units |
| closure events | 123.890088 | reserve K p M = 2^117.596544 events (K = 50), 66 adder operations each |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts | 29.573647 | 2^20 start elements, two exponentiations each (preprocessing) |
| restarts at the trail cap | 60.640271 | M / 2^66 cap resets, one adder charge each |
| distinguished points | 70.473717 | capacity 2^67.3944 records (store, restart, scatter, transfer, overflow) and the bucket table of 2^68 slots (sort and scan) |
| fingerprint verification | 74.362379 | Markov charge at 2^60 times the expected false pairs (2^5.79 at capacity), one exponentiation each |
| restart table | 18.573647 | 1024 exponentiations g^(f_j), one per restart index (preprocessing, section 6) |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

The reserves enlarge M by the factor 1/(1 − w): the revisit reserve w = 127 K p = 0.286004 is the largest term in the increase.

## 4. Distinguished-point table, pathwise capacity, scratch memory and peak memory

Records. Each record comes from a DP: a fresh DP flag (probability 2^-60 per new state, H7 and H8) or a duplicate at a join. New states are
the M main steps and the R restarts (DP restarts, escapes, cap resets); R <= M. The records are at most X + J pathwise, where X is the number
of fresh flags over the M + R new states and J the number of joins. Overflow needs X > 2 mu (mu = (M + R) 2^-60, failure exp(-mu/3)) or
J > J_MAX (failure E[J]/J_MAX). The capacity is D = ceil(2 mu) + J_MAX = 2^67.3944 records.

Buckets. The table is 2^37 buckets of 2^31 slots each, 64-bit indices (`measurement/sortnet/sortnet64.c`), partitioned by the top 37 bits
of the fingerprint (bucket = fp[1] >> 27). The mean occupancy is D / 2^37 = 2^30.3944, against capacity 2^31; the Chernoff bound per
bucket and the union over the buckets give a failure of 2^-184681256 (negligible). Each bucket is sorted by a bitonic network
(n k (k+1)/4 compare-exchanges with n = 2^31, k = 31, i.e. 248 per slot); its scratch is one 80-byte temporary. Equal fingerprints share a
bucket, so matching needs no global sort.

Costs (instructions, per record or per slot): the store (1496.0, measured alone by storebench), one restart (26740.1),
the bucket scatter (24.0, measured), transfer and index widening (100,000, an explicit assumption), overflow (10); per slot: the bitonic
compare-exchange (108.003, the 32-bit measurement, larger than the 64-bit values 108.0000054 and 108.0000027), the scan (13.0), and one cursor reset per
bucket (1). The charge is over the capacity records and all 2^68 slots (padding included). The table costs about 2^-56.4 of W_MAX.

Peak memory. The table occupies 2^68 slots of 80 bytes (2^74.321929 bytes in total, spread over the buckets; each bucket 2^37.32 bytes, below the 64-bit
address space), plus the adder table (12 x 1024 elements, about 7.8 MiB), the walk states (2^20 x 640 bytes) and detectors (per walk 98 + 16 + 8
bytes, including the restart counter). The peak is not realisable on the hardware used here and is reported, not scored.

## 5. Evidence: toy-scale runs under the v10 walk

Group: the order-ell subgroup of F_(q^12)^* (`cscore/gt.py`, revision v10: restart table of 1024 elements, replay-invariant escape, round-based
anchors, DP test before the detector). Production parameters (r_a = 1024, S = 458) are used at bits 18 (ell = 922561, N = 76880, 1015 distinct
classes among the 1024 adders); bits 14 (N = 2366) cannot host 1024 adders meaningfully (827 distinct classes), so the r_a = 1024 runs use bits 18
only. Seeds, trial counts and the two walks per trial are fixed in `RESULT-toy-v10.json`.

Success at tau = 0.5277 (model 1 − exp(−tau) = 0.41), bits 18, 800 trials, dp_bits 6, two walks:

| run | success [Wilson 95%] | mean tau (model 1.0) | restart fraction |
| --- | ---: | ---: | ---: |
| add-only, r_a 1024, v10 | 0.436 [0.402, 0.471] | 0.94 | 1.4% |
| mixed (S 458), r_a 1024, v10 | 0.406 [0.373, 0.441] | 1.00 | 1.3% |
| add-only, r_a 1024, v9 restart rule (escape as v10) | 0.414 | 1.00 | |
| mixed, r_a 1024, v9 restart rule (escape as v10) | 0.390 | 0.98 | |

At r_a = 32 (S = 15, bits 14 and 18, dp_bits 2), the v10 runs give 0.435, 0.473, 0.444 and 0.448 (add-only and mixed at each size), against 0.41:
shifts of +0.03, +0.06, +0.03 and +0.04, with mean tau 0.95, 0.82, 0.99 and 0.89. The shift is the closure effect of a small r_a: the two-cycle rate
per main step is 1.5e-3 to 2.6e-3 at r_a = 32 and 4.2e-5 (mixed) and 7.9e-5 (add-only) at r_a = 1024, against the formulas (1 − s)/(12 r_a) and
1/(12 r_a), which the toy confirms at both sizes. At production closure rates the toy success agrees with the model within the Wilson width.

Restart checks (v10 against v9 at the same seeds): undo restarts (a restart output equal to the class before the walk's last operation) fall to
the chance level (v10 17, 21, 4, 2 against expected 20.5, 21.8, 0.9, 1.0; v9 152, 95, 204, 136). Duplicate distinguished points (equal (A, B))
fall from 219, 114, 328 and 203 under v9 to 0, 0, 2 and 2 under v10; at r_a = 1024 both rules give 0. Restart outputs that land on a visited class
match the H8 prediction (sum of visited/N) within z = −0.3 to +1.8 under v10 (v9: +2.5 to +15).

Replay (Lemma R, `replay_v10.py`, bits 14, r_a 1024): T1, 2000 pairs (two variants of identical class with different anchor phases, previous keys
and exponent offsets, and a sigma-twisted start): 0 violations in 3480 and 3523 compared steps before the first detection or DP; both walks reach
a DP at the same point in 999 and 997 pairs, with equal canonical exponents. T2, DP test disabled, 4096 steps, 30 instances: 381 detections matched
381, escape elements coincide in every matched pair, zero offsets 0.

Escape loops. In the add-only r_a = 1024 sweep (800 trials, dp_bits 6), 541 escapes produced 379 landings on already visited classes. Six trials
hold 348 of the landings and 354 escapes, and none of them solves. The diagnosis (toy-v10 `loop_diag.py`) finds no bug: no escape target lies in its
cycle, no restart element is the identity or of short order, and the escape formula checks against the rule. The loop is the rule's consequence: the
escape target is a function of the cycle, so a walk that re-enters the cycle escapes to the same element again, and the walk's class map brings the
target back to the cycle in about 118 steps (median). The trial's 52-member cycle attracts 29% of random classes within 1500 steps, against 0.15 for
a random mapping of this size, so this toy has a large attractor. First-visit re-entry (target not yet escaped to) is 0.075 against L |C| / N = 0.034
(within a factor of two), so the rate is a small-N effect, but a re-entered DP-free loop traps the walk until its cap. This is premise H10.

Coincidences. The adder elements and the restart elements coincide with the expected random rate at toy size (2 · 1024 · r_a · 12 / ell per walk:
2794 against 2770, 83378 against 83094, 21948 against 21823). Production has no such coincidence (zero in both tracks, `gen_inputs.py`).

Square share (H9). The fresh-step square share is 0.442 against s = 0.447 at r_a = 1024 (z = −6.0), 0.451 against 0.469 at r_a = 32 (z = −14.4),
and 0.471 against 0.469 (z = +2.3). H9 is not supported at toy size, and the margin of section 3 is used in the production bound.

The toy keys are not the production keys: the toy key is a base-31 encoding minimised over the 12 conjugates (class-invariant, 0 violations in 300
samples), with a strong bias toward small leading digits; the production key is the minimum over the 12 conjugates of one 381-bit coefficient, not
measured here. Toy results are consistent with the model at r_a = 1024, and do not test the production key.

## 6. Limitations

- H1 is extrapolated from toy subgroups with N up to 76880 and 1024 adders (bits 18) to N about 2^251.
- H5 is a heuristic: the non-formal closure term depends on the coins being independent of the path, and it does not address adaptively chosen paths.
- H7 (new states), H8 (restart table), H9 (square share) and H10 (trapped trails) are stated premises. The toy contradicts H9 at its own key
  construction and shows escape loops at N = 76880; the production values are extrapolations.
- The lemma R replay checks (`measurement/replica.py`) are at small densities and on the toy; at production size the replay lemma is the proof.
- The toy has no trail cap (the production walk has one at 2^66), and its restart-share and dp density differ from production.
- The DP_BITS 0 run is not defined under the chained rule and is not charged; the record store and the restart are measured alone.
- The kernel's dp_exact_match compares the 65,536-slot record ring; at DP_BITS 2 with 4096 rounds the record count (87,360) exceeds the ring, so
  the record contents are compared on the ring slots and the counts in full.
- The bucket organisation is measured with 64-bit indices at 2^21 and 2^22 slots; the production bucket size (2^31) is not run. The transfer cost
  is an assumption (100,000 instructions per record).
- The reference operation is charged at the lower organizer count (22,990); the count measured in this run is 22,993.1, and blst_fp12_mul measured
  22,998.
- The solve, verify and output charge of 2^29 instructions is an assumed bound (2^14.7 of the units).
- Preflight: run with the installed verifier (`llm-verifier preflight`); the outcome is in the package notes.

## References

- Lam and Leung, On vanishing sums of roots of unity (2000), for the background of Proposition G (the statement is proved here).
- Pollard (1978); van Oorschot and Wiener (1999); Teske (1998), section 4.2; Bernstein and Lange (2012), section 3.
- `work-gt/RESEARCH.md` sections 3-4 and `walk-spec.md` (this package).

## Reproduction

```sh
python3 -I accounting.py                 # bundled: W_MAX, failure terms, number-theory asserts
echo '{"bits":14,"trials":100,"seed":1,"s":0.46,"adder_bits":5,"dp_bits":2,"walks":2}' | python3 -I experiment.py
bash measurement/reproduce.sh DIR        # Docker linux/arm64, organizer image, py_ecc; restart checks at small densities
```
`reproduce.sh` generates inputs at S = 458, measures the kernel (Callgrind at zero DP density, two sizes, mul and sqr loops, the DP-every-step
run and the restart), checks the replica at production and small DP densities, runs the native counts, and runs the sort measurements
(`measurement/sortnet/`).

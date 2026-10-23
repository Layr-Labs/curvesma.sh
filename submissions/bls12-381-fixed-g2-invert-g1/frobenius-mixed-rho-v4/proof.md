# frobenius-mixed-rho-v4: Frobenius-class mixed add/square rho on GT (revision v9)

Claim: for the uniform challenge of `rules.md` (Q uniform nonidentity in G2, z uniform in GT, fresh algorithm coins), the procedure below
outputs P in G1 with e(P, Q) = z with probability at least 0.41 (model 0.432782661565; lower bound 0.410000000001 after the failure terms of
Section 2, which include a 2% reserve loss for closure events), while its total work, counting every operation it performs, is
at most W_MAX = 2^126.823186 target-group-operation equivalents (score 63411593/500000). The first M = 2^126.374296 main steps are executed before
the cutoff except with probability 2.28e-02 (Section 3), under the premises H1, H5 and H7 stated in Section 2. Peak memory is at most
2^74.321929 bytes, including the padding of the bucketed distinguished-point table (Section 4); preprocessing is inside W_MAX (2^29.575064
units); nonuniform advice is 0 bytes. It is a classical generic attack, not a practical break.

**Revision v9.** Round 9 (on revision v8, both tracks) raised the points below. Point (1) is a defect in the restart rule; points (2) to (4)
answer the accounting and table objections; point (5) answers the toy objection.

(1) Restart recurrence (defect). v8 restarted a walk at x * X0 with a fixed per-walk start element. A walk that returned to an earlier
DP of its own could then repeat that trail forever. v9 restarts with an adder chosen by a per-walk restart counter t (incremented before
use) and the walk's index w through a fixed 64-bit mixer (walk-spec section 6). A restart state is itself a DP candidate and chains.
(2) Rare-event work by Markov reserves. Closure events (first revisits of a window of 64 steps) are bounded in expectation, not by a
conditional concentration: their expected number is at most p M with p the unconditional per-step rate, so a reserve of K p M events
(K = 50) fails with probability 1/K = 0.02 (charged to delta). Each reserved event costs at most 127 revisit steps (counted as main steps)
and 66 operations. Joins are bounded by a Markov reserve of J_MAX = 1000 with failure E[J]/J_MAX (failure 2.78e-03, included in delta);
each join's suffix is at most CAP steps.
(3) Detection delay. The residual wait of a colliding pair to its next DP is geometric with mean 2^60. The allowance is 40 mean waits per
walk (lag = WALKS * 40 * 2^60 steps), with e^-40 in delta and the collision factor (1 - e^-40) in the success model.
(4) Table organisation (defect in the 2^67-slot claim). The table is 2^37 independent buckets of 2^31 slots with 64-bit
indices, partitioned by the top 37 bits of the fingerprint, so each bucket fits a 64-bit address space and equal keys share a bucket.
Costs are measured with 64-bit indices (c_cmp64 = 108.0000 per compare-exchange at 2^21 and 2^22 slots; scan 13.0; scatter 24.0 per
record); the 32-bit value 108.0030 is kept as the charge (the larger). Transfer and index widening are an explicit assumption of
100,000 instructions per record. The table cost is about 2^-56.4 of W_MAX and does not drive the score.
(5) The toy experiment now uses the v9 walk (restart counter, chained restarts, round-based anchors). Its success at the toy density is
not a production estimate, because the toy restarts at about 23-25% of its steps (dp_bits 2); production restarts are 2^-60 of steps
(Section 5).

Score: v8 126.102185 to v9 126.823186 (+0.721001 bits). The increase is from the revisit reserve: w = 0.286004 of the main steps is reserved
(K = 50 expected closure events, each followed by at most 127 revisit steps), so M grows by 1/(1 - w) = 1.4006 (+0.486 bits), and each
reserved event is charged 66 adder operations (the reserved events add 0.1486 of the walk charge). The failure total is 2.28e-02.

## 1. Procedure

Reduction. g = e(P0, Q) generates GT, where P0 is the fixed G1 generator; z = g^y gives P = [y]P0 with e(P, Q) = z. Input a nonidentity Q in G2 and a target z in GT. Return the identity if z = 1 (P = O is valid).
Otherwise run the walks until a verified y is found or the work counter reaches W_MAX, then output P = [y]P0, verified by the pairing
equation e(P, Q) = z (one pairing, charged). Exceptional inputs (Q = O) are rejected by the rules.

Walk (`walk-spec.md`, revision v9, implemented by the blst kernel `measurement/step.c` and its replica `measurement/replica.py`):
- GT has order r. Frobenius sigma(x) = x^p has order 12 with lambda = p mod r, lambda^6 = -1; classes are sigma-orbits,
  N = (r - 1)/12, about 2^251.272127.
- Key K(x): the minimum over the 12 values V_n(sigma^n(x)), computed with four Fp multiplications; K is class-invariant.
- State (x, a, b) with x = g^(a + b y). Square slots S = 458 of r_a = 1024 (s = 0.447265625000). Step: h = K mod r_a. If h < S, square x and
  double (a, b). Otherwise multiply x by T[h][k] = sigma^(-k)(R_h), where R_h = g^(c_h) z^(d_h), and update (a, b) by lambda^(-k) (c_h, d_h).
- Trails: a trail ends at its first distinguished point. There the walk stores its record and restarts: with t the walk's restart
  counter (incremented before use) and w the walk's index, the restart index is h = mix((w + 1) * 0x9E3779B97F4A7C15 + t) mod r_a, and the
  restart is the additive step with index h. The restarted state is itself tested for a DP (chained restarts). A restart clears the walk's
  detector (previous key, anchor). A trail that reaches 2^66 steps restarts the same way without a record.
- Cycles (per walk): a new key equal to the key two steps back, with the following step returning to the old key, is a 2-cycle; a new key
  equal to the walk's anchor key (recorded every 64 walk steps; restarts are not steps) inside a window of 64 steps is a cycle of length at
  most 64. Either escapes by squaring the least-key state of the cycle.
- Distinguished point: (K >> 10) & (2^60 - 1) = 0. The record is (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of the
  Montgomery integer; A = lambda^kn a and B = lambda^kn b mod r, where kn is the argmin index.
- Table and solving: the records are stored in the bucket table of Section 4; each bucket is sorted by fingerprint. Records with equal
  fingerprints are compared pairwise. A pair from one class with B1 != B2 gives y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is
  discarded; no trail is replayed. Every candidate is verified by g^y = z before output.
- Work counter and cutoff: each operation adds its charge (Section 3). The walks run while W < W_MAX - reserve. The table stops the
  procedure when a bucket or the table is full.

Consistency checks: the walk invariant x = g^(a + b y) holds on every state, including after restarts; the replica reproduces the kernel's
dump, DP records and counts bit for bit at production density (DP_BITS 60) and at small densities that exercise chained restarts (DP_BITS 2,
3, 4) and the trail cap (CAP overridden to 50 and 24), at 512 and 4096 rounds (`measurement/RESULT.json`); the bitonic buckets are
checked for sorting (`measurement/sortnet/check_bitonic.py`, `sortnet64.c`).

## 2. Success probability, failure terms and the success chain

Under H1 the walks need T = 2^125.888284 fresh collision steps for tau = T^2 V / (2N), with V = 0.799413681 at the nominal square fraction (the
square fraction of the new states is within 10^-6 of s, Section 3). A collision is detected with probability at least 1 − e^-40 per pair,
and the trail cap and detection factors give f = 1.000000000000. The success model is 1 − exp(−f tau), tau_eff = 0.567012736.

Success chain (each step names its term). (i) The first M main steps execute before the cutoff: the charges of Section 3 are within W_MAX.
(ii) The closure events among them are at most K p M (Markov, failure 1/K = 0.02). (iii) The revisit steps are at most 127 times (ii) plus the
join suffixes (at most J_MAX CAP), so the fresh steps are at least (1 − w) M with w = 127 K p (+ 2^-50), and M is chosen so that this is at
least T + lag. (iv) The distinguished-point records are at most X + J (Section 4), overflow has probability at most the DP terms below.
(v) Given T fresh steps, a collision occurs between walks with probability 1 − exp(−f tau_eff) under H1. (vi) A collision is found when both
walks reach a DP; the residual is geometric with mean 2^60, and it fits the allowance lag = WALKS * 40 * 2^60 except with probability e^-40 per
pair. A join is itself a collision: the joined walk's duplicate record is paired with the other walk's record, and the two B values differ except
with probability 1/r (degenerate pairs are negligible). Each step has its failure term; their sum is the failure total.

P = 0.432782661565 before delta and 0.410000000001 after it. The success target 0.41 is met because the model is set to 0.41 + delta (accounting.py,
`fresh_steps_needed`). If the waste were not budgeted, P would fall to 0.228248.

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
| trap term over the run | 1.86e-18 | M times the per-step total (includes the trap revisit steps) |
| fingerprint Markov (metered matches) | 8.67e-19 | 2^-60 |
| square fraction, Hoeffding over the square decisions | 2^(-3.18e+26) | exp(-2 (1e-6)^2 M); the revisit allowance is in w |
| closure events (Markov, reserve K = 50) | 0.02 | 1/K; the reserve K p M is charged |
| non-formal closures of length <= 64 | 1.6e-36 | expectation 12 * 64 * M / r (Markov) |
| DP table, fresh flags above 2 mu | 2^(-4.6e+19) | exp(-mu/3), Chernoff on X (Section 4) |
| DP table, joins above J_MAX = 1000 | 2.78e-03 | E[J]/J_MAX (Markov), E[J] <= M^2/N (Section 4) |
| detection delay beyond 40 mean waits per walk | 4.25e-18 | e^-40 (Section 2) |
| DP table, bucket overflow | 2^(-197318825) | Chernoff per bucket, union over 2^37 buckets (Section 4) |
| total failure kept | 2.28e-02 | sum of the kept terms |

Premises (disclosed; each is a stated assumption, not a proof).
- **H1 (random map).** The next state of a step is visited by another trail with probability at most (visited)/N; tau is the random-map
  collision statistic with rate V. Joins use this bound.
- **H5 (non-formal closures).** A closure with a nonzero coefficient has probability 1/r over the coins. The non-formal trap term rests on it.
- **H7 (uniform new states).** Given the history, the key of each new state is uniform: its adder index and square decision are uniform, its DP
  flag has probability 2^-60, and its square decision is Bernoulli(s). The restart index (the mixer of the counter) is treated as uniform. The
  closure bound is an expectation bound (Markov) over the process coupled with independent keys until the first revisit, so it needs only the
  unconditional per-step rate p.

## 3. Computation and the work cutoff

Reference unit: one dense Fp12 multiplication (`blst_fp12_mul`), charged 1 as the rules require; the organizer's reference count (REF_FUNCS)
is 22,990 instructions. Other work is measured in instructions and divided by C = 20,000.

Measurement (`measurement/RESULT.json`): the organizer's reference image (linux/arm64, GCC 14.2.0, Valgrind 3.24.0, `-O3
-fno-tree-vectorize -mgeneral-regs-only`), the pinned blst commit f7a88067, Callgrind Ir, forced inputs at production DP_BITS 60.

| forced step type (production DP_BITS 60, zero distinguished points) | instructions | notes |
| --- | ---: | --- |
| all-adder steps, seed 1 | 26,737.203 | v9 kernel; includes 4 detected cycles (the maximum over seeds, charged as A) |
| all-adder steps, seed 2 | 26,729.902 | v9 kernel; includes 1 detected cycle |
| all-square steps, both seeds | 13,518.469 | no detected cycles (charged as Q) |
| restart alone (mode restartbench) | 26,537.0 | counter, mix, table step, key derivation, state clear; charged one adder step |
| store of one record alone (mode storebench) | 1,495.0 | fingerprint and two canonical exponents; charged as the store |
| bitonic compare-exchange (sortnet, 32-bit indices, 2^12 slots) | 108.003 (108.001 at 2^13) | branch-free mask and 80-byte conditional swap; charged (the larger) |
| compare-exchange (sortnet64, 64-bit indices, 2^21 and 2^22 slots) | 108.0000054 and 108.0000027 | same network with 64-bit indices |
| bitonic scan, per slot | 13.0 | equal-neighbour test, branch-free; charged |
| scatter of one record into its bucket (sortnet64) | 24.0 | charged |
| ring exchange (kernel mode sortbench) | 48.0 | micro-benchmark of one comparison and one exchange; not charged |
| reference operation (organizer, REF_FUNCS) | 22,990 | charged 1 unit; blst_fp12_mul measured 22,997 |
| cyclotomic square | 9,803 | charged 9,803 / C |

Charges for the first M main steps:
- Main steps: every step is charged the adder cost 1.18836 units (the larger of the two step charges).
- Closure events: reserve K p M = 2^117.579720 events (Markov, failure 1/K = 0.02), each charged 66 adder operations (0.1766 units per main step
  in expectation). The revisit steps of these events (at most 127 each) are main steps and are counted in M through the reserve w.
- Square decisions: Azuma over the square decisions of the first M steps (conditionally Bernoulli(s) under H7), deviation 10^-6, failure
  exp(-2 * 10^-12 M). The collision rate uses s.
- Joins: reserve J_MAX = 1000 (Markov, failure E[J]/J_MAX <= M^2/(N J_MAX)), suffixes at most CAP each (in w).
- Restarts: each record and each trail-cap reset costs one restart, charged one adder step (A_UNITS; the measured restart is
  26,537 instructions).
- Fingerprint matches: metered (one exponentiation each); the metered total is bounded by 2^60 times the expected count (Markov, failure 2^-60).
- Non-formal closures of length at most 64: expectation 12 * 64 M / r (Markov).

Components of W_MAX (log2 units; the score is the rounded-up cutoff 126.823186):

| component | log2 units | what it covers |
| --- | ---: | --- |
| walk steps | 126.623268 | main steps M = 2^126.374296, each charged the adder cost 1.18836 units |
| closure events | 123.873086 | reserve K p M = 2^117.579720 events (K = 50), 66 adder operations each |
| adder table | 19.582053 | 1024 adders, two exponentiations and 12 Frobenius maps each (preprocessing) |
| walk starts | 29.573647 | 2^20 start elements, two exponentiations each (preprocessing) |
| restarts at the trail cap | 60.623268 | M / 2^66 cap resets, one adder charge each |
| distinguished points | 70.458481 | capacity 2^67.3743 records (store, restart, scatter, transfer, overflow) and the bucket table of 2^68 slots (sort and scan) |
| fingerprint verification | 74.322240 | Markov charge at 2^60 times the expected false pairs (2^5.75 at capacity), one exponentiation each |
| solve, verify, output | 14.712288 | 2^29 instructions and one pairing |

The reserves enlarge M by the factor 1/(1 − w): the revisit reserve w = 127 K p = 0.286004 is the largest term in the increase.

## 4. Distinguished-point table, pathwise capacity, scratch memory and peak memory

Records. Each record comes from a DP: a fresh DP flag (probability 2^-60 per fresh state, H7) or a duplicate at a join. The records are at most
X + J pathwise, where X is the number of fresh flags and J the number of joins. Overflow needs X > 2 mu (mu = M 2^-60, failure exp(-mu/3)) or
J > J_MAX (failure E[J]/J_MAX). The capacity is D = ceil(2 mu) + J_MAX = 2^67.3743 records.

Buckets. The table is 2^37 buckets of 2^31 slots each, 64-bit indices (`measurement/sortnet/sortnet64.c`), partitioned by the top 37 bits
of the fingerprint (bucket = fp[1] >> 27). The mean occupancy is D / 2^37 = 2^30.3743, against capacity 2^31; the Chernoff bound per
bucket and the union over the buckets give a failure of 2^-197318825 (negligible). Each bucket is sorted by a bitonic network
(n k (k+1)/4 compare-exchanges with n = 2^31, k = 31, i.e. 248 per slot); its scratch is one 80-byte temporary. Equal fingerprints share a
bucket, so matching needs no global sort.

Costs (instructions, per record or per slot): the store (1495.0, measured alone by storebench), one restart (26737.2),
the bucket scatter (24.0, measured), transfer and index widening (100,000, an explicit assumption), overflow (10); per slot: the bitonic
compare-exchange (108.003, the 32-bit measurement, larger than the 64-bit values 108.0000054 and 108.0000027), the scan (13.0), and one cursor reset per
bucket (1). The charge is over the capacity records and all 2^68 slots (padding included). The table costs about 2^-56.4 of W_MAX.

Peak memory. The table occupies 2^68 slots of 80 bytes (2^74.321929 bytes in total, spread over the buckets; each bucket 2^37.32 bytes, below the 64-bit
address space), plus the adder table (12 x 1024 elements, about 7.8 MiB), the walk states (2^20 x 640 bytes) and detectors (per walk 98 + 16 + 8
bytes, including the restart counter). The peak is not realisable on the hardware used here and is reported, not scored.

## 5. Evidence: toy-scale runs under the v9 walk

Group: the order-ell subgroup of F_(q^12)^* with ell dividing q^4 - q^2 + 1; `cscore/gt.py`, bundled into `experiment.py`, with the v9 restart rule
(restart counter, chained DP restarts, round-based anchors; restarts are not walk steps for the anchors). The toy taus use the partition-correct V
for the nominal s. The toy has dense distinguished points (dp_bits 1 or 2), so restarts make up about 23-25% of its steps (restart share per run:
23.4%, 23.0%, 23.4%, 24.9%, 24.6%); production restarts are 2^-60 of steps. The toy is therefore a test of the collision statistics and the solver at high restart density, not a production estimate (see the
paragraph after the table).

| run | ell | N | s (effective) | toy V | trials | reached a useful collision | mean tau (model 1) | success at tau = 0.5277 [Wilson 95%] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| initial input (bits 14, seed 1, 100 trials) | 28393 | 2366 | 0.4688 | 0.763672 | 100 | 100 | 1.0856 | 0.4400 [0.3467, 0.5377] |
| add-only, bits 14 | 28393 | 2366 | 0.0000 | 0.968750 | 3000 | 3000 | 1.0164 | 0.4017 [0.3843, 0.4193] |
| mixed, bits 14 | 28393 | 2366 | 0.4688 | 0.763672 | 3000 | 2992 | 0.8450 | 0.4743 [0.4565, 0.4922] |
| add-only, bits 18 | 922561 | 76880 | 0.0000 | 0.968750 | 800 | 800 | 0.9432 | 0.4325 [0.3986, 0.4671] |
| mixed, bits 18 | 922561 | 76880 | 0.4688 | 0.763672 | 800 | 800 | 0.8883 | 0.4625 [0.4282, 0.4971] |

Pooled mixed runs at tau = 0.5277 (N = 2366 and N = 76880): 1793 of 3800 succeed, empirical 0.4718, Wilson 95% interval
[0.4560, 0.4877], against the model value 1 - exp(-tau) = 0.4100. The same inputs under the pre-v9 walk (dense distinguished points, no
restarts) gave 0.4277 (bits 14) and 0.4088 (bits 18) at this tau, against 0.4743 and 0.4625 under v9. In the v9 runs the mixed
walk's mean tau is 0.845 and 0.888 against the model's 1.0, and the square share of all steps falls to 0.343 and 0.355 (restart
steps are additive, so the square share of walk steps is not the share of all steps). The step mix is one candidate for the shift; the toy does
not isolate it. The toy success is therefore a test of the collision statistics and the solver at high restart density, not a production estimate.

## 6. Limitations

- H1 is extrapolated from toy subgroups with N up to 2^17 and 32 adders to N about 2^251 and 1024 adders.
- H5 is a heuristic: the non-formal closure term depends on the coins being independent of the path. H7 is the uniform-new-state premise of
  Section 2, with the restart index treated as a uniform function of the counter; it is a stated premise, not a proof.
- The closure reserve (K = 50) and the join reserve (J_MAX = 1000) are Markov bounds. Their failure terms (0.02 and 2.78e-03) are the largest
  parts of delta and are paid in the success model; the detection delay (e^-40) is not.
- The toy experiment (`experiment.py`) implements the v9 walk. Its restart share is 23-25% of steps, far above production, and its success is
  above the model at that density (Section 5); the cause is not isolated.
- Restarts are verified by the replica at small densities (DP_BITS 2, 3 and 4, chained restarts, trail cap 50 and 24) at 512 and 4096 rounds, with
  dump, records and counts compared. At DP_BITS 2 and 4096 rounds the kernel has 87,360 records: their counts agree, and their contents are
  compared on the kernel's 65,536-slot record ring. At production density no restart occurs in a 2^20-round native run.
- The all-distinguished-point run (DP_BITS 0) is not defined under the chained rule and is not charged; the record store and the restart are
  measured alone (storebench, restartbench).
- The bucket organisation is measured for its per-record and per-slot instruction costs with 64-bit indices at 2^21 and 2^22 slots; the
  production bucket size (2^31) is not run. The transfer cost is an assumption (100,000 instructions per record, far above a local copy).
- The kernel ring exchange (`sortbench`, 48.0 instructions) is a micro-benchmark of one comparison and one exchange; it is not charged.
- The reference operation is priced by the organizer rule at 1 unit (organizer measurement 22,990; blst_fp12_mul measured 22,997).
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

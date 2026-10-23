# bls12-381-fixed-g1-invert-g2: Pollard rho in GT on Frobenius classes

Status: exploratory analytic bound, claimed score log2(W) = 126.48. Two heuristics, H1 and H2, are disclosed in "Heuristics and
limits"; nothing else is assumed. Constants (p, r, G1/G2 generators, field and pairing definitions) are the protected profile's.
Every number below is recomputed by `ledger.py` in this folder (python3, standard library only).

## Reduction

Input: P in G1 (nonidentity) and z in GT. Output Q in G2 with e(P,Q)=z.

1. If z=1 return infinity. Otherwise compute g = e(P,Q0) for the profile's G2 generator Q0 (one reduced pairing).
   Bilinearity and nondegeneracy of the reduced pairing on the prime-order groups G1 x G2 give g != 1, and GT has
   prime order r, so g generates GT. Hence z = g^x for a unique x in [0,r).
2. Find x = log_g z by parallel collision search (below).
3. Return Q = [x]Q0. Then e(P,Q) = e(P,Q0)^x = g^x = z by bilinearity. The attack checks g^x = z in GT and also
   recomputes e(P,Q) and compares all 12 coefficients.

The run is a function of (P,z) only; the success probability is over z uniform in GT (x uniform) and P uniform nonidentity,
with fresh coins. Success does not depend on P beyond g != 1, which holds for every admissible P.

## Frobenius classes

Because the embedding degree is 12 (r divides p^12-1 and no smaller p^i-1), p has multiplicative order exactly 12 modulo r,
so the Frobenius map phi: s -> s^p acts on GT as exponentiation by lambda = p mod r, an automorphism of order 12, and
phi^6 is inversion (s^(p^6) = s^(-1); checked numerically on a random GT element). phi has no fixed point except 1 (lambda^i != 1 for 0<i<12), so
GT minus identity splits into (r-1)/12 classes of exactly 12 elements {s^(p^i)}. N' = (r-1)/12 + 1, log2 N' = 251.27213.
If class(g^a z^b) = class(g^a' z^b'), then g^a z^b = (g^a' z^b')^(p^k) for some k and
x follows from a + b x = lambda^k (a' + b' x) mod r unless b = lambda^k b' (probability about 1/r for random starts).

Frobenius on Fp12 = Fp2[w]/(w^6 - xi) (xi = 1+u, w^12 - 2w^6 + 2 = 0): write s = sum c_k w^k, c_k in Fp2, k=0..5. Then
phi^i(s) = sum conj^i(c_k) gamma_i^k w^k with gamma_i = xi^((p^i-1)/6) in Fp2, a fixed public constant. Define the coordinate
v_i(s) = Re(conj^i(c_1) gamma_i), the Fp2-real part of coefficient 1 of phi^i(s) (the polynomial coefficient c1+c7 in the
profile's basis), an integer in [0,p). (i) v_i(phi^m(s)) = v_(i+m)(s), so the set {v_i} depends only on the class. (ii) For s in GT,
phi^6(s) = s^-1 is the conjugation over Fp6 (w -> -w), so v_(i+6) = -v_i mod p; only v_1..v_5 (v_0 is the element's own coefficient) need
computing, and the 12 values are v_i and p - v_i. Both facts were checked numerically on a random GT element (12 distinct values,
v_(i+6) = -v_i, rotation under phi, inversion = phi^6). The class representative ("canonical element") is the conjugate with the
smallest value; its value V is the class's canonical coordinate and label.

## Walk

Precompute a table T of R = 2^16 elements t_j = g^(a_j) (a_j from a public seed, uniform in [0,r)) and its conjugates T[j][i] = phi^(-i)(t_j),
i=0..11 (12 R Fp12 elements, 453 MB per processor; every processor derives the same table: 2^16 full exponentiations g^(a_j) and
11 Frobenius maps per entry, charged below).
State: an Fp12 element s of GT. Step: compute v_1..v_5, find the index i* of the minimum of {v_i, p - v_i} (i* in 0..11, using
phi^6 for the negated ones), set V = that value (the canonical coordinate), j = bits 60..75 of V, and s <- s * T[j][i*].
Because phi^i*(s * phi^(-i*)(t_j)) = phi^i*(s) * t_j, the new s lies in the class of canon(s) * t_j: the class-level map
f(class) = class(canon * t_j) is a well-defined function on classes, with no full Frobenius evaluation. One step = one dense Fp12
multiplication (1 reference operation, charged 1) plus the instruction work in "Cost".
Distinguished point (DP): V has its low 60 bits zero (theta = 2^-60; a class invariant). The table index uses bits 60..75 of V, disjoint
from the DP bits. On a DP output (V, start seed, length). A walk starts at s0 = g^a0 z^b0 with (a0,b0) from a fresh seed (one
simultaneous exponentiation, about 255 squarings and at most 191 multiplications, charged 512 operations per start; restarts after every
DP or after 20/theta steps without one). The DP test comes first: a cycle that contains a DP class never traps a walk, because every walk
reaching the cycle stops at it.

Lemma 1 (index and DP bits). Assume H1: for a uniformly random class the six values v_i, i = 0..5, are independent and uniform on [0,p). Then the twelve values are
the v_i and p - v_i, and V = min_i min(v_i, p - v_i) has a nonincreasing mass function g (the survival function ((p - 2v)/p)^6 is convex) with g(0) = 1 - (1 - 2/p)^6 <= 12/(p-1). For any k and any c
in [0, 2^k), the low k bits F = V mod 2^k satisfy |2^k P(F=c) - 1| <= 2^k g(0) <= 12 * 2^k / (p-1).
Proof. Upper: P(F=c) = g(c) + sum_{t>=1} g(c + t 2^k), and since g is nonincreasing, g(c + t 2^k) is at most the average of g over the block
(c + (t-1) 2^k, c + t 2^k], so the sum is at most 2^-k P(V > c) <= 2^-k; with g(c) <= g(0) this gives P(F=c) <= g(0) + 2^-k. Lower: P(F=c) >= sum_{t>=0} g(c + t 2^k)
>= 2^-k P(V > c) >= 2^-k (1 - (c+1) g(0)) >= 2^-k (1 - 2^k g(0)). QED.
With k = 76 (DP bits 0..59 and index bits 60..75 together), every value of the 76-bit field has probability within a relative 12 * 2^76/(p-1) = 2^-301.1
of 2^-76, so the pair (DP bits, index) is within statistical distance 2^-301.1 of two independent uniform fields, at every step; over the at most 2^125.9
executed steps the total distance is at most 2^-175.2. This is why the table index and the DP test are independent, and why the table index is uniform
on R = 2^16 values up to that distance. The lemma is about the low bits of the minimum; it does not claim that "the minimum changes only high bits".
Only the uniformity of the coordinates of a random GT element (H1) is assumed.

## Fruitless cycles and the cycle rule

Write f for the class map above; it is a function of the class alone. A step j' = j taken after the
phi^6 conjugation undoes the previous step, so f has 2-cycles at rate exactly 1/(12R) = 2^-19.585 per step (the next index repeats with
probability 1/R and the canonical conjugate is the inverse one with probability 1/12); a walk without a rule would loop forever.
The rule is defined on classes, not on the walker. Every walker keeps the last P = 24 states (class label, Fp12 element, exponents a,b mod r;
704 bytes each) of its logical path. After each step it compares the new label with the P previous ones.
(1) If no label repeats, nothing happens: the trajectory is exactly f's.
(2) If the new label equals the one w steps back (w <= P), the walk has gone once round a cycle of f of length L = w <= P. The label comparison
is independent of when the walker started, and every cycle of length <= P is detected on its first return, whatever class the walker entered
it at. The walker has the L states of the cycle in its buffer; it finds the member m with the smallest label (distinct classes have distinct
labels, so m is unique and depends only on the cycle), restores the stored state of m's first visit, discards the states after it (the lap),
and escapes from m by squaring its canonical element: s <- canon(canon(m)^2), the class of the square (a cyclotomic squaring in GT, the analogue of
point doubling; no table is involved). The retained path is therefore ... -> m -> esc(m), where esc(m) = class(canon(m)^2) is a function of the
class m only; the exponents become (2a, 2b) mod r. The lap and the discarded steps are charged as waste.
Why the escape does not fall back into the cycle: squaring is not a table step, so the next ordinary step from esc(m) does not undo it
(an escape by an ordinary table step is undone by the following step with probability 1/(12R) = 2^-19.6, re-forming the cycle). The square is also not in m's own
class: m^2 = m^(lambda^k) would need lambda^k = 2 mod r, and 2 is not in the order-12 group generated by lambda (2^12 != 1 mod r).
Under H1 esc(m) is an essentially uniform class: it lies on the cycle's own L classes with probability at most 24/N' = 2^-246.7. The reverse coincidence,
that the f-trail from esc(m) later enters the cycle, has probability at most 24 C / N' = 2^-182.4 per event (a trail of at most C = 20/theta = 2^64.3 steps
from a fresh class visits a given class with expected count about C/N'). If it did happen the rule would fire again with the same m and the same esc(m),
and the walker would loop to the cap C. Expected loss: at most K omega events x 2^-182.4 x C = 2^108.6 x 2^-182.4 x 2^64.3 = 2^-9.5 steps (negligible).
Consequence: the modified map f' (f'(c) = esc(c) if c is the smallest member of a cycle of f of length <= P, else f(c)) is again a function
of the class alone, and the walker realises f' exactly. Two walks that reach the same class, at any times, therefore have the same class
sequence from then on until the next DP, whatever their histories: both find the same m and escape to the same esc(m). Walks that merged stay merged.
Effect of the escape on the walk's randomness: f' differs from f only at the cycle-minimum classes, a set of density at most omega/5 < 2^-19, and at those classes the
successor esc(m) is again a uniform-looking class (a random-mapping step). The escape therefore removes at most the retained steps that the executed-step
fraction omega counts (2^-17.3) and does not bias the other steps; it is not the source of the 20% margin below, which covers the walk's randomness in general.

(3) Cycle probabilities for every length. Let rho_L be the probability that a visited class c lies on a cycle of f of length L, i.e. f^L(c) = c. In additive
notation for the exponent u of the canonical element, one step is u -> lambda^kappa (u + a_j) (kappa in Z_12 the conjugation chosen, j the index), so after L steps
u_L = lambda^(K_L) u_0 + sum_i lambda^(K_L - K_(i-1)) a_(j_i), with K_i = kappa_1 + ... + kappa_i. A closed cycle needs u_L = u_0. We treat the conjugation choices kappa and the indices as independent uniform (H2). Case A (a table index with a non-vanishing combined coefficient, or
lambda^(K_L) != 1): the equation is one linear condition modulo r on random data, probability about 1/r per start and length; it is the unstructured closure of a random mapping and is
bounded below by the closure term. Case B (structured): lambda^(K_L) = 1 and every table index j that occurs in the L steps has combined coefficient
c_j = sum of lambda^(K_L - K_(i-1)) over its occurrences i equal to 0 mod r. The coefficient is an integer combination alpha of twelfth roots of unity with m_j nonnegative terms; the map
Z[zeta_12] -> F_r (zeta -> lambda) kills alpha only if r divides the norm of alpha, and the norm is at most m_j^4 < r while m_j < 2^63 (an index used 2^63 times in a walk
of at most 2^64.3 steps has probability below 2^-(2^62) under H1), so c_j = 0 mod r exactly when alpha = 0 in Z[zeta_12]. A single twelfth root cannot vanish, so every index that occurs
must occur at least twice; and for an index with m >= 2 occurrences the last root is determined by the others, so the probability that its sum vanishes is at most 1/12.
If the L steps use b distinct indices, each at least twice, then b <= floor(L/2); the probability that L uniform indices stay inside some set of b values is at most C(R,b)(b/R)^L,
and the b vanishing conditions are independent with probability at most 1/12 each. Hence
  rho_L <= B_L = sum_{b=1}^{floor(L/2)} C(R,b) (b/R)^L 12^-b   for every L >= 3 (b also at most R),   and rho_2 = 1/(12R) exactly.
Values (python3, `ledger.py`): B_3 = 2^-35.6, B_4 = 2^-36.2, B_24 = 2^-177.8; the largest bound for any L >= 25 is B_25 = 2^-190.2 and B_L decreases with L until L = 2R = 2^17, where
B_L < 2^-172955, and stays below 2^-234944 for L >= 2^40 (the term b = R, 12^-R). For short lengths this is the usual picture: a 2k-cycle needs k table coincidences, about (1/R)^k with combinatorial factors.
Cycles of length <= 24 are detected and charged: omega_short = sum over L = 2..24 of rho_L (L + 3) = 2^-17.263, of which the exact 2-cycle term 5/(12R) is all but 2^-35: this part of
the ledger does not depend on H2 at all. (The rate (1 - theta)/(12R) actually met by a trail is slightly lower because a step that lands on a DP ends the trail first; we charge 1/(12R).)
Cycles longer than 24 are not detected by the rule. A walk that enters one without meeting a DP is cut after C = 2^64.3 steps and restarted. Structured long cycles: the expected number entered per
visited class is at most sum over L >= 25 of B_L <= C x 2^-190.2 = 2^-125.9 (L up to the cap C; each term at most B_25), each costing at most C steps, so the waste per executed step is
at most C^2 2^-190.2 = 2^-61.6 and is added to omega. Closure term (all lengths beyond the window, structured or not): a walk trapped in a cycle without a DP revisits a class, so its trail has a
self-collision within the cap. Under H1 a trail of T steps has a repeated class with probability at most T(T-1)/(2N'); with trail lengths geometric of mean 1/theta and capped at C,
E[min(T,C)^2] <= 2/theta^2 = 2^121, and there are at most K theta + M = 2^65.9 trails, so the expected number of trapped trails is at most 2^65.9 x 2^121 / (2 x 2^251.27) = 2^-65.4, each wasting at most C = 2^64.3
steps: expected loss 2^-1.1 steps, charged in the ledger (under one step). A trapped walk never produces a wrong collision, it only restarts.
Total wasted fraction charged: omega = omega_short + 2^-61.6 = 2^-17.263 (relative), the closure term as an absolute 2^-0.6.
Cost of the rule per step: 24 label comparisons (8-byte hashes of the label, about 80 instructions) and one 704-byte state copy (about 150 instructions); per event: the walker executes the
L lap steps, which are discarded, and the escape (a squaring plus the canonicalisation, at most one step, counted as two), plus one for the restore and minimum search: L + 3 steps per event.

## Collision and recovery

DPs are stored by V (a collision of two stored labels is an exact class collision, up to a 2^-300 false match). Both walks
reached the same class at their DP, so both canonical elements are the same element u: g^(a+bx) = g^(a'+b'x). On a match,
rerun each of the two walks alone from its seed (the same rule, so the same logical path), tracking (a,b) mod r of the canonical
element by cheap integer arithmetic, and read off (a,b) and (a',b') at the DP. If b != b' (probability 1 - 2^-254) then x =
(a - a')/(b' - b) mod r; the candidate x is confirmed by g^x = z. This costs at most 2 x 20/theta steps, no lockstep alignment is needed,
and the DP's identity is the class label. Then Q = [x]Q0 and the final pairing check.

## Success probability

Model (H1, the standard heuristic): the class-level step function f' behaves like a random mapping on the N' classes up to the additive-walk penalty. Literature for the walk: van Oorschot and Wiener 1999
(parallel collision search with DPs); Teske 2001 (r-adding walks, r >= 16, within a few percent of the random-walk collision time, the deviation shrinking with r); Bernstein, Lange and Schwabe 2011,
"On the correct use of the negation map" (walk design with an equivalence map, fruitless cycles, and the additive-walk collision penalty about sqrt(1 + 1/r) for an r-entry table);
Bos, Kaihara, Kleinjung, Lenstra and Montgomery 2010 (the same on a large-scale run with the negation map); Duursma, Gaudry and Morain 1999 (the gain sqrt(|G|) for walks on classes under an automorphism group G).
These analyses are for elliptic-curve groups with negation or Frobenius; GT has no published rho analysis at this size, so applying them is the extrapolation, and the toy run below checks that the walk behaves like them at small scale.
Effective space: N_eff = (1 + 1/R) N' (the penalty; 1 + 2^-16). In the random-function model each newly visited class is an independent uniform draw from N_eff classes, so after k distinct
visited classes P(no collision between two different trails) <= exp(-k(k-1)/(2 N_eff)) up to a factor 1 + O(1/k) and a self-collision correction below 2^-65 (above), i.e. P(collision) >= 1 - exp(-k^2/(2 N_eff)).
Setting k0^2 = 2 ln(1/0.61) N' gives P = 0.39 with N' (log2 k0 = 125.62779). The attack runs K = mu x sqrt(1 + 1/R) x k0 logical steps in total (mu = 1.20: an explicit 20% margin on steps beyond the penalty,
not derived: it absorbs any deviation of this walk from the random-mapping model, including the cycle escape and the correlation between table index and class), log2 K = 125.89083.
Success argument: let k = K (1 - delta) be the number of useful distinct draws, where delta <= 2^-41.6 is the share of steps in partial trails at the end (the detection delay of M = 2^20 walkers, below). Then
P >= 1 - exp(-k^2/(2 N_eff)) = 1 - 0.61^(mu^2 (1-delta)^2) = 0.509, above 0.39. Equivalently the claim holds as long as the true effective collision space N_real is at most mu^2 = 1.44 times N_eff.
A collision is useful with probability 1 - 2^-250 (b != lambda^k b'). Because walks that share a class stay together (cycle rule), the first shared class of two trails is detected at their common next DP.
The probability space is uniform instances and fresh coins (start seeds, table seed), with fresh coins per run; the identity (z=1, probability 2^-255) is handled directly.
Toy evidence (`experiment.py`, run at its default inputs, about 1 s; toy scale only): the same walk (canonical coordinate = minimum of the 12 coordinate values, DP bits 0..d-1, index bits d.., window 24,
smallest-label squaring escape) on a group of 60013 elements (N' = 5002 classes) with an order-12 automorphism, 1200 runs with a 64-entry table and 800 runs with a 16-entry table (the smallest Teske covers),
fresh instance, table and starts per run, trails cut at 20/theta. Mean steps to the first collision between two trails, over the random-draw value sqrt(pi N'/2) + 2/3: 1.043 (95% interval 1.013 to 1.073) for 64 entries and
1.036 (1.000 to 1.072) for 16 entries, against the production factor mu sqrt(1 + 1/R) = 1.20. The fraction of runs that collide within the production step budget is 0.49 (0.46 to 0.52) and 0.51 (0.47 to 0.54),
in line with the model's 0.51 and above 0.39. All 2000 collisions recover x exactly. Escapes by squaring fall back into the same cycle 0 times in 170 events; the control
that escapes by an ordinary table step (index j+1) falls back in 164 of 338. Pairs of walks that share a class have identical class sequences afterwards: 0 of 2066 diverge. The measured 2-cycle rate per
executed step is 0.0012 (64 entries) and 0.0044 (16 entries) against the charged 1/(12R) = 0.0013 and 0.0052 (a step that lands on a DP ends the trail first, so the rate met is a little lower). These are toy numbers:
a 2^12.3-class group, a pseudorandom 40-bit stand-in for the coordinate (so it tests the walk, not the bit uniformity of the Fp coefficient, which Lemma 1 covers under H1), and some excess from random self-collisions
(5 and 19 in the runs) that are negligible at production trail lengths. They support, and do not prove, the production margin.

## Cost

Step cost. The Fp12 multiplication is the reference operation (1). Other per-step work in AArch64 instructions: 5 coordinates, each two Fp
multiplications (383 instructions per blst_fp_mul, COST_MODEL) plus an add: 5 x 806 = 4030; comparing 12 six-limb values
12 x 12 = 144; DP test and index extraction (low bit fields) 20, cycle-window compare 80, state copy 150, table address 12, loads and stores of the state (not already
inside the measured Fp12 multiplication) 30: total 4466 instructions. We charge 8000 (about 1.8x) to cover cache and memory transfers of
a table of this size: overhead 8000/20000 = 0.4. Step cost c = 1.4 reference operations.
Logical steps K = mu sqrt(1 + 1/R) k0, log2 K = 125.89083 (log2 r = 254.85709, log2 N' = 251.27213, log2 k0 = 125.62779).
Executed steps K (1 + omega), omega = 2^-17.263, log2 = 125.89084; times c: 2^126.37627.
Other charges, all in reference operations (instructions / 20000 where applicable), with theta = 2^-60, M = 2^20 processors, 20/theta = 2^64.3:
- starts: (K theta + M) x 512 = 2^74.89 (K theta = 2^65.89 starts, one per DP);
- DP processing, communication, storage and lookup at 10^4 instructions each: 2^64.89;
- detection delay and unfinished trails of the M walkers, M x 20/theta x c = 2^84.81 (a relative 2^-41.6 of K);
- recovery: two reruns of at most 2 x 20/theta steps at 4x cost for the integer bookkeeping: 2^67.81;
- cycles: executed steps include omega (above); the trapped-trail closure term 2^-0.6;
- preprocessing, the per-processor table: each entry t_j = g^(a_j) is a full exponentiation, 512 units (as for a start), so 2^20 processors x 2^16 entries x 512 = 2^45.000; the 11 x 2^16 Frobenius maps
  per processor (at most 1 unit each) add 11 x 2^36 = 2^39.459; together 2^45.031, which we round up and charge as 2^45.1 (the same figure everywhere: preprocessing_log2 = 45.1);
- the pairing e(P,Q0), the final pairing check (2^14 each, a generous price for a reduced pairing in this model), the G2 scalar multiplication (2^10, at most 510 Fp2 group operations of
  at most 20 Fp2 multiplications of 3 x 383 instructions = 2^9.2 reference operations) and the GT check g^x = z (2^9): total below 2^15.1.
Sum of the itemised charges: W0 = 2^126.37627 (the delay is the largest term beyond the steps, 2^-41.6 relative; the rest are smaller).
Reserve: we multiply by an explicit unallocated 1.07 (7%) against ledger errors, W = 1.07 W0 = 2^126.47388, and submit log2 W rounded up to two decimals, 126.48 (= 3162/25).
Slack of the claim over the itemised sum: 0.104 bits, over the reserved ledger 0.006 bits. The claim is below the track ceiling of 126.6.
(arithmetic: log2 K = log2 k0 + log2 1.20 + 1.1e-5 = 125.62779 + 0.26303 + 0.00001; + log2 1.4 = 0.48543 gives 126.37626; omega adds 0.00001; + log2 1.07 = 0.09761.)

## Memory, preprocessing, advice

Stored DPs: (K theta + M) = 2^65.89 records of 80 bytes (48 byte label, 16 byte seed, 16 byte length) in a hash table at load at most 1/2, 160 bytes per record: 2^73.21 bytes.
Per-processor state: the table 12 x 2^16 x 576 B = 2^28.75 B, the table exponents 2^16 x 32 B, the window 24 x 704 B (a state is a 576 B Fp12 element, a 48 B label, two 32 B exponents
and 16 B padding) and 4 KB of code and scratch: 2^28.76 B per processor, 2^48.76 B for 2^20 processors. The total stays below 2^73.3 bytes. Peak memory declared: at most 2^74 bytes.
The DP rate 2^-60 is chosen to bring memory down from 2^104 bytes at a small time cost: each walk is longer (2^60 steps), so the delay terms grow to at most 2^-41.6 of the work, and the cycle window (24) is sized for that length.
Preprocessing: only the per-processor table (2^45.031, charged 2^45.1, included in W; it depends on g, hence on P, and is not reusable). Advice: none (2^0 bytes).
Memory has no scalar cost in the profile.

## Heuristics and limits

- H1: (a) the coordinates v_i of a uniformly random GT element behave as independent uniform values mod p (Lemma 1 and the table index use only this); (b) the class-level walk f' with R = 2^16 behaves like a random
  mapping on N' classes up to the charged penalty sqrt(1 + 1/R), covered by the explicit factor mu = 1.20 (Teske; Bernstein-Lange-Schwabe; Bos et al.; Duursma-Gaudry-Morain; van Oorschot-Wiener). This is the standard heuristic for
  additive walks with an equivalence map, applied to GT; no proof exists at this size. A deviation up to a 44% larger effective space (mu^2) still gives success 0.39; the toy shows 4% (interval up to 7%).
- H2: the cycle probabilities rho_L of f over the visited classes: exact for L = 2 (which is all of omega up to 2^-35), a counting bound B_L for L >= 3 (a table index in a cycle must occur at least twice, vanishing sums
  of twelfth roots of unity), summed over all L; anything unstructured is the random-mapping closure term (2^-1.1 steps). Cycles of length <= 24 are detected and left by a squaring escape; longer ones are cut by the walk cap.
- The per-step instruction overhead is a charged bound (8000 against 4466 itemised), not a measured count; the explicit 7% reserve is not derived.
- Memory 2^74 bytes is far from physical; this is a classical analytic bound, like the baseline's 2^136-byte accounting control.
- No claim of novelty or a practical attack. For this track GT rho on Frobenius classes uses 12-element classes (a factor sqrt(12) = 2^1.79 over plain rho), with no pairing inversion algebra and no quantum steps.

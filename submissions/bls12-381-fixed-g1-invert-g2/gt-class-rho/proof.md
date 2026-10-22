# bls12-381-fixed-g1-invert-g2: Pollard rho in GT on Frobenius classes

Status: exploratory analytic bound, claimed score log2(W) = 126.16. Heuristics H1 and H2 below are disclosed; nothing
else is assumed. Constants (p, r, G1/G2 generators, field and pairing definitions) are the protected profile's.

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
GT minus identity splits into N'-1 = (r-1)/12 classes of exactly 12 elements {s^(p^i)}. N' = (r-1)/12 + 1, log2 N' = 251.2721.
Searching for a collision between classes: if class(g^a z^b) = class(g^a' z^b'), then g^a z^b = (g^a' z^b')^(p^k) for some k and
x follows from a + b x = lambda^k (a' + b' x) mod r unless b = lambda^k b' (probability about 1/r for random starts).

Frobenius on Fp12 = Fp2[w]/(w^6 - xi) (xi = 1+u, w^12 - 2w^6 + 2 = 0): write s = sum c_k w^k, c_k in Fp2, k=0..5. Then
phi^i(s) = sum conj^i(c_k) gamma_i^k w^k with gamma_i = xi^((p^i-1)/6) in Fp2, a fixed public constant. Define the invariant
v_i(s) = Re(conj^i(c_1) gamma_i), i.e. the Fp2-real part of coefficient 1 of phi^i(s) (equal to the poly coefficient c1+c7 in the
profile's basis). (i) v_i(phi^m(s)) = v_(i+m)(s), so the set {v_i} depends only on the class. (ii) For s in GT,
phi^6(s) = s^-1 is the conjugation over Fp6 (w -> -w), so v_(i+6) = -v_i; only v_1..v_5 (v_0 is free) need computing, and the 12
values are v_i and p - v_i. The 12 values are distinct with overwhelming probability. Both facts were checked numerically
on a random GT element (12 distinct values, v_(i+6) = -v_i, rotation under phi, inversion = phi^6).
The class representative is the conjugate with the smallest invariant; vmin is its class-invariant label.

## Walk

Precompute a table T of R = 2^16 elements t_j = g^(a_j) (a_j random in [0,r)) and its conjugates T[j][i] = phi^(-i)(t_j), i=0..11
(12 R Fp12 elements, 453 MB per processor; computed once per processor by 11 Frobenius maps and 2^16 multiplications, charged below).
State: an Fp12 element s of GT. Step: compute v_1..v_5, find the index i* of the minimum of {v_i, p - v_i} (i* in 0..11, using
phi^6 for the negated ones), set vmin = that value, j = bits 40..55 of vmin, and s <- s * T[j][i*].
Because phi^i*(s * phi^(-i*)(t_j)) = phi^i*(s) * t_j, the new s lies in the class of canon(s) * t_j: the class-level map
f(class) = class(canon * t_j) is a well-defined function on classes, with no full Frobenius evaluation. One step = one dense Fp12
multiplication (1 reference operation, charged 1) plus the instruction work in "Step cost".
Distinguished point: vmin has its low 60 bits zero (theta = 2^-60; a class invariant; the table index uses bits 200..215 of vmin, disjoint
from them). On a distinguished point output (vmin, start seed, length). A walk starts at s0 = g^a0 z^b0 with (a0,b0) from the seed
(about 255 squarings and multiplications of the cyclotomic group, charged 512 operations per start, restarts after every distinguished
point or after 20/theta steps without one). Distinguished test comes first: a cycle that contains a distinguished class never traps a walk,
because every walk reaching the cycle stops at it.

Fruitless cycles and the cycle rule. Write f for the class map above; it is a function of the class alone. A step j' = j taken after the
phi^6 conjugation undoes the previous step, so f has 2-cycles at rate exactly 1/(12R) = 2^-19.585 per step (the next index repeats with
probability 1/R and the canonical conjugate is the inverse one with probability 1/12); a walk without a rule would loop forever.
The rule is defined on classes, not on the walker. Every walker keeps the last P = 24 states (class label, Fp12 element, exponents a,b mod r;
576 bytes each) of its logical path. After each step it compares the new label with the P previous ones.
(1) If no label repeats, nothing happens: the trajectory is exactly f's.
(2) If the new label equals the one w steps back (w <= P), the walk has gone once round a cycle of f of length L = w <= P. The label comparison
is independent of when the walker started, and every cycle of length <= P is detected on its first return, whatever class the walker entered
it at. The walker has the L states of the cycle in its buffer; it finds the member m with the smallest vmin (distinct classes have distinct
labels, so m is unique and depends only on the cycle), restores the stored state of m's first visit, discards the states after it (the lap),
and escapes from m: s <- canon(m) * T[j(m)+1 mod R][i*], the ordinary step with table index j(m)+1. The retained path is therefore
... -> m -> esc(m), where esc(m) is a function of the class m only. The lap and the discarded steps are charged as waste.
Consequence: the modified map f' (f'(c) = esc(c) if c is the smallest member of a cycle of f of length <= P, else f(c)) is again a function
of the class alone, and the walker realises f' exactly. Two walks that reach the same class, at any times, therefore have the same class
sequence from then on until the next distinguished point, even when they entered the cycle at different members, went round different numbers of
times or have different histories: both find the same m and escape to the same esc(m). Walks that merged stay merged.
(3) Rates. A cycle of length L of f, for generic s0, needs the table terms to cancel: sum over steps of lambda^kappa * a_j = 0 mod r with the a_j
independent, so every table index met in the cycle must occur at least twice (vanishing sums of 12th roots of unity have at least 2 terms), which
costs at least ceil(L/2) coincidences of probability 1/R. We bound the per-step rate by rho_L <= L! R^-ceil(L/2) for L >= 3 (the L! covers all
orderings and conjugation patterns; for L = 4 it is 2^-27.4 against a direct count of about 2^-39) and use the exact rho_2 = 1/(12R). Heuristic
H2 is this counting. The undetected lengths L > 24 have total rate at most 2^-119.5 per step, and a walk trapped in one is cut after
20/theta = 2^64.3 steps, so the expected loss is 2^-55 of the steps. The escape lands on a class that again lies on a short cycle with
probability below 2^-200 per event, and then the rule simply applies again.
(4) Cost of the rule per step: 24 label comparisons (8-byte hashes of the label, about 80 instructions) and one 576-byte state copy (about
150 instructions); per cycle event: finding the minimum over L labels and restoring the state (under 500 instructions). Waste per event: the
walker executes the L lap steps and the escape step but retains only the steps up to the minimum member and the escape, so at most L steps
are wasted; we charge L + 3 steps.

## Collision and recovery

Distinguished points are stored by vmin (a collision of two stored labels is an exact class collision, up to a 2^-300 false match). Both walks
reached the same class at their distinguished point, so both canonical elements are the same element u: g^(a+bx) = g^(a'+b'x). On a match,
rerun each of the two walks alone from its seed (length known, the same rule, so the same logical path), tracking (a,b) mod r of the canonical
element by cheap integer arithmetic, and read off (a,b) and (a',b') at the distinguished point. If b != b' (probability 1 - 2^-254) then x =
(a - a')/(b' - b) mod r; the candidate x is confirmed by g^x = z. This costs at most 2 x 20/theta steps, no lockstep alignment is needed,
and the distinguished point's trail identity is the class label. Then Q = [x]Q0 and the final pairing check.

## Success probability

Model H1: the class-level step function behaves like a random mapping on the N' classes (r-adding walks with R = 2^16 are random for
this purpose: van Oorschot-Wiener 1999 for parallel collision search with distinguished points; Teske 2001 shows R >= 16 is within a few percent of the
random-walk collision time, and the deviation shrinks with R). In the random-function model each newly visited class is an
independent uniform draw, so after k distinct visited points P(no repeated class) = prod (1 - i/N') <= exp(-k(k-1)/(2N')) (
tight up to a factor 1+O(k^-1)), i.e. P(collision) >= 1 - exp(-k^2/(2N')) up to the k term that is below 2^-125. Setting
k^2/(2N') = -ln(0.61) gives P >= 0.39, k = sqrt(0.98859 N') = 2^125.6278. The attack runs for K = 1.03 k logical steps in total (a 3% margin for walk non-randomness and restarts), executed steps include the explicit cycle waste below, and a collision detected is useful with probability 1 - 2^-250. Because walks that share a class stay together (cycle rule), the first shared class of two trails is detected at their common next distinguished point. The detection delay is the walkers' last partial segments and the travel to the next distinguished point: at most M = 2^20 parallel walkers times 20/theta = 2^64.3 steps, i.e. 2^84.3 steps of K = 2^125.67, a relative loss of 2^-41.3, so the delivered points exceed k; the delay is charged explicitly below. The run succeeds with
probability >= 0.39. Identity (z=1, probability 2^-255) is handled directly.

## Cost

Step cost. The Fp12 multiplication is the reference operation (1). Other per-step work in AArch64 instructions: 5 invariants, each two Fp
multiplications (383 instructions per blst_fp_mul, COST_MODEL) plus an add: 5 x 806 = 4030; comparing 12 six-limb values
12 x 12 = 144; distinguished test 8, cycle-window compare 80, state copy 150, table index and address 12, loads and stores of the state (not already
inside the measured Fp12 multiplication) 30: total 4454 instructions. We charge 8000 (about 1.8x) to cover cache and memory transfers of
a table of this size: overhead 8000/20000 = 0.4. Step cost c = 1.4 reference operations.
Logical steps: K = 1.03 k, log2 k = 125.62779 (log2 r = 254.85709, log2 N' = 251.27213), log2 K = 125.67043. K c = 2^126.15586.
Other charges, all in reference operations (instructions / 20000 where applicable), with theta = 2^-60, M = 2^20, 20/theta = 2^64.3:
- cycle rule: executed steps = K (1 + omega), omega = sum over L = 2..24 of rho_L (L + 3) with rho_2 = 1/(12R) and rho_L <= L! R^-ceil(L/2)
  (L >= 3): omega = 2^-17.25, so the charge is K c omega = 2^108.9 and the total is K c (1 + omega) = 2^126.15587;
- starts: (K theta + M) x 512 = 2^74.7 (K theta = 2^65.7 starts, one per distinguished point);
- distinguished point processing, communication, storage and lookup at 10^4 instructions each: 2^64.7;
- detection delay and unfinished trails of the M walkers, M x 20/theta x c = 2^84.8;
- recovery: two reruns of at most 2 x 20/theta steps at 4x cost for the integer bookkeeping: 2^67.8;
- loss to cycles longer than 24 (undetected, cut after 20/theta steps): 2^-55 of the steps, 2^71 (negligible);
- table construction per processor (2^16 multiplications + 11 x 2^16 Frobenius maps) at M = 2^20 processors: 2^38.2;
- the pairing e(P,Q0), the final pairing check (2^14 each, a generous price for a reduced pairing in this model), the G2 scalar multiplication (2^10, at most 510 Fp2 group operations of at most 20 Fp2 multiplications of 3 x 383 instructions = 2^9.2 reference operations) and the GT check g^x = z (2^9): total below 2^15.1.
Sum: W = 2^126.15587 (1 + 2^-41.3 + ...) <= 2^126.1560. The largest addition beyond the cycle-rule total is the delay term, 2^-41.3 relative; the rest are smaller. We submit 126.16, rounded up from the sum by 0.004.
(arithmetic: log2 (1.03 k) = 125.67043, + log2 1.4 = 0.48543, + log2 (1 + omega) = 0.00001.)

## Memory, preprocessing, advice

Stored distinguished points: K theta = 2^65.67 records of 64 bytes: 2^71.67 bytes (labels are 48 bytes plus a 16 byte seed and length).
Per-processor tables 12 x 2^16 x 576 B = 2^29.1 B; 2^20 processors: 2^49.1 B, so the total stays below 2^71.7 bytes. Peak memory <= 2^72 bytes.
The distinguished-point rate 2^-60 is chosen to bring memory down from 2^104 bytes at a negligible time cost: each walk is longer (2^60 steps),
so the delay terms above grow to at most 2^-41 of the work, and the cycle window (24) is sized for that length.
Preprocessing: only the per-processor table (2^38.2, included in W; it depends on g, hence on P, and is not reusable). Advice: none (2^0 bytes).
Memory has no scalar cost in the profile.

## Heuristics and limits

- H1: class-level random-mapping behavior with R = 2^16 (Teske; van Oorschot-Wiener), for the modified class map f' of the cycle rule. A 3% margin is charged.
  The toy run (experiment.py) shows the method end to end on a group of order about 6e4 with an order-12 automorphism, and tests the cycle rule:
  with a 6-entry table (2-cycles about 1.4% per step) every pair of random walks that shares a class has identical class sequences from there on (0 of
  thousands of merged pairs diverge), while an escape that depends on the walker's own history diverges in some pairs.
- H2: the rates rho_L of fruitless cycles of f above (exact for L = 2, a counting bound for L >= 3). Cycles of length <= 24 are detected and their
  waste is charged; longer ones are cut and bounded at 2^-55 of the steps.
- Memory 2^72 bytes is far from physical; this is a classical analytic bound, like the baseline's 2^136-byte accounting control.
- No claim of novelty or a practical attack. For this track GT rho on Frobenius classes uses 12-element classes (a factor sqrt(12) = 2^1.79 over plain rho), with no pairing inversion algebra and no quantum steps.

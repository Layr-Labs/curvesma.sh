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
Distinguished point: vmin has its low 28 bits zero (theta = 2^-28; a class invariant). On a distinguished point output
(vmin, start seed, length). A walk starts at s0 = g^a0 z^b0 with (a0,b0) from the seed (about 255 squarings and multiplications of
the cyclotomic group, charged 512 operations per start, restarts after every distinguished point or after 20/theta steps without one).
Fruitless cycles: a step j' = j after a conjugation by phi^6 undoes the previous step (probability about 1/(12R) = 2^-19.6 per step).
The walker keeps the last 4 vmin labels; if vmin repeats within the window it uses index j+1 instead (Bos-Kaihara-Kleinjung-Lenstra-
Montgomery 2010; Bernstein-Lange-Schwabe 2011). Longer fruitless cycles need an additional coincidence among table entries and have
probability below 2^-39 per step; a trapped walk is cut by the 20/theta limit (waste about 2^-11 of the steps, inside the margin).

## Collision and recovery

Distinguished points are stored by vmin (a collision of two stored labels is an exact class collision, up to a 2^-300 false match).
On a match from two seeds, rerun both walks (lengths known, aligned by running the longer one first) in lockstep to find the first
equal labels; this costs at most 2 x 20/theta steps and tracks (a,b) mod r of the two canonical conjugates by cheap integer
arithmetic. The relation a + b x = lambda^k (a' + b' x) is solved for the 12 values of k; the candidate x is confirmed by g^x = z.
Then Q = [x]Q0 and the final pairing check.

## Success probability

Model H1: the class-level step function behaves like a random mapping on the N' classes (r-adding walks with R = 2^16 are random for
this purpose: van Oorschot-Wiener 1999 for parallel collision search with distinguished points; Teske 2001 shows R >= 16 is within a few percent of the
random-walk collision time, and the deviation shrinks with R). In the random-function model each newly visited class is an
independent uniform draw, so after k distinct visited points P(no repeated class) = prod (1 - i/N') <= exp(-k(k-1)/(2N')) (
tight up to a factor 1+O(k^-1)), i.e. P(collision) >= 1 - exp(-k^2/(2N')) up to the k term that is below 2^-125. Setting
k^2/(2N') = -ln(0.61) gives P >= 0.39, k = sqrt(0.98859 N') = 2^125.6278. The attack runs for K = 1.03 k steps in total
(a 3% margin for walk non-randomness, delayed detection, restarts, fruitless-cycle waste, trapped walks), and a collision
detected is useful with probability 1 - 2^-250. The detection delay is the walkers' last partial segments: M = 2^20 parallel
walkers times 20/theta = 2^32.3 steps, i.e. 2^52.3 steps of K = 2^125.67, a relative loss of 2^-73, so the delivered points exceed k. The run succeeds with
probability >= 0.39. Identity (z=1, probability 2^-255) is handled directly.

## Cost

Step cost. The Fp12 multiplication is the reference operation (1). Other per-step work in AArch64 instructions: 5 invariants, each two Fp
multiplications (383 instructions per blst_fp_mul, COST_MODEL) plus an add: 5 x 806 = 4030; comparing 12 six-limb values
12 x 12 = 144; distinguished test 8, history compare 32, table index and address 12, loads and stores of the state (not already
inside the measured Fp12 multiplication) 30: total 4256 instructions. We charge 8000 (about 1.9x) to cover cache and memory transfers of
a table of this size: overhead 8000/20000 = 0.4. Step cost c = 1.4 reference operations.
Steps: K c = 1.03 x 2^125.6278 x 1.4 = 2^126.1559.
Other charges, all in reference operations (instructions / 20000 where applicable):
- starts: (K theta + M) x 512 = 2^106.7 (K theta = 2^97.7 starts, one per distinguished point);
- distinguished point processing, communication, storage and lookup at 10^4 instructions each: 2^96.7;
- detection delay of the M walkers: M x (1/theta) x c = 2^48.5;
- recovery: two reruns of at most 2 x 20/theta steps at 4x cost for the integer bookkeeping: 2^35.8;
- table construction per processor (2^16 multiplications + 11 x 2^16 Frobenius maps) at M = 2^20 processors: 2^38.2;
- the pairing e(P,Q0), the final pairing check (2^14 each, a generous price for a reduced pairing in this model), the G2 scalar multiplication (2^10, at most 510 Fp2 group operations of at most 20 Fp2 multiplications of 3 x 383 instructions = 2^9.2 reference operations) and the GT check g^x = z (2^9): total below 2^15.1.
Sum: W = 2^126.1559 (1 + 2^-19.5 + 2^-29.5 ...) <= 2^126.1560. We submit 126.16, rounded up from the sum by 0.004.
(cost.py arithmetic: log2 r = 254.85709, log2 N' = 251.27213, log2 k = 125.62779, log2 (1.03 k) = 125.67043, + log2 1.4 = 0.48543.)

## Memory, preprocessing, advice

Stored distinguished points: K theta = 2^97.7 records of 64 bytes: 2^103.7 bytes (labels are 48 bytes plus a 16 byte seed and length).
Per-processor tables 12 x 2^16 x 576 B = 2^29.1 B; 2^20 processors: 2^49.1 B. Peak memory <= 2^104 bytes.
Preprocessing: only the per-processor table (2^38.2, included in W; it depends on g, hence on P, and is not reusable). Advice: none (2^0 bytes).
Memory has no scalar cost in the profile; the exponential memory is declared and could be traded for time by raising theta's denominator, since stored records scale with theta.

## Heuristics and limits

- H1: class-level random-mapping behavior with R = 2^16 (Teske; van Oorschot-Wiener). A 3% margin is charged. The toy run (experiment.py) shows the method end to end on a group of order about 6e4 with an order-12 automorphism and recovers x in about two thirds of runs at about 1.1 sqrt(-2 ln 0.61 N') steps plus delay.
- H2: fruitless cycles are rare and handled by the window escape; the 20/theta cut bounds the loss from any residual cycles.
- Memory 2^104 bytes is not physically realisable; this is a classical analytic bound, like the baseline's 2^136-byte accounting control.
- No claim of novelty or a practical attack. For this track GT rho on Frobenius classes uses 12-element classes (a factor sqrt(12) = 2^1.79 over plain rho), with no pairing inversion algebra and no quantum steps.

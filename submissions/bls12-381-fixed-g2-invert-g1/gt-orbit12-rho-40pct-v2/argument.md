# GT order-12 Frobenius-class parallel Pollard rho (fixed-G2, backport v2)

Backport of the historical unpublished local submission
`04399700-428c-4db7-9268-1719315b0883` (track bls12-381-fixed-g2-invert-g1,
method "Frobenius-class parallel Pollard rho in GT", claimed time_log2 =
125.89), independently re-derived under the current COST_MODEL.md and
re-packaged for this queue. This is an analytic upper-bound construction,
not a production run, practical break, accepted score, or novelty claim.

Claim: **time_log2 = 125.89** (total work W <= 2^125.89 dense-Fp12
multiplication equivalents at the fixed C = 20,000 price, scalar 64-bit
AArch64 instruction equivalents), success probability >= 0.40 (track floor
0.39), peak memory < 2^73 bytes (reported, unscored), preprocessing
<= 2^62 equivalents included in W, nonuniform advice 0 bytes. One disclosed
heuristic (H1, random-mapping behaviour of the class walk), supported by the
historical toy-field runs; no production-size execution is claimed.

## 1. Reduction (this track's convention, unchanged from the source branch)

Input: nonidentity Q in G2, uniform z in GT. The protected profile pins
e(P,Q) = py_ecc.optimized_bls12_381.pairing(Q, P, final_exponentiate=True)
at the pinned commit; Q is supplied and cannot be replaced. Compute
g = e(P0, Q) via pairing(Q, P0, True). Q has prime order r and the pairing
is nondegenerate, so g != 1 and g generates GT. Every valid z equals g^x
for a unique x in [0, r). Solve x, output P = [x]P0, and verify
pairing(Q, P, True) == z. If z = 1, output P = infinity. The answer is for
the exact supplied Q; no other Q is used. No auxiliary powers and no n
parameter exist on this track; the score is log2(W).

## 2. Arithmetic facts (independently re-verified; see certificate.json)

Tower as in blst: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u;
the profile's flat basis converts freely (u = w^6 - 1, w^12 - 2w^6 + 2 = 0).

* ord_r(p) = 12 exactly (checked: p^d != 1 for d in {1,2,3,4,6}, p^12 = 1).
  The automorphism group {y -> y^{p^k}} of GT therefore has order 12 and
  acts freely on GT \\ {1}: every y != 1 lies in a class [y] of exactly 12
  elements; class count N = (r-1)/12, log2 N = 251.2721.
* p^6 = -1 (mod r) (checked), so y -> y^{p^6} is inversion on GT.
* 1 + p^4 + p^8 = 0 (mod r) (checked).
* zeta = xi^{(p^2-1)/6} lies in Fp, is a primitive 6th root of unity
  (zeta^2 = zeta - 1 and zeta^3 = -1 mod p; both checked on the protected p).
* Class key from 4 Fp multiplications: for the w^1 coefficient c1 = a + b'u
  of y, the w^1 real parts of the 12 conjugates y^{p^k} are
  key_k = a*zeta^j (k = 2j) and key_k = B*zeta^j (k = 2j+1), where
  B = Re(conj(c1)*d1), d1 = xi^{(p-1)/6}, and Re(conj(c1)*d1) =
  a*Re(d1) + b*Im(d1) (identity checked on random Fp inputs, and
  N(d1) = zeta checked). Since zeta^j is in {1, zeta, zeta-1, -1, -zeta,
  1-zeta}, the 12 keys cost 4 Fp multiplications (a*zeta, the two products
  of B, B*zeta) plus 8 subtractions/negations. No full Frobenius of the
  whole element is ever evaluated.
* kappa(y) = argmin_k key_k (canonical residues in [0,p) compared as
  integers) is constant on conjugacy classes and identical for all 12
  members of a class, so y^{p^kappa} is the canonical class representative
  and K([y]) = min_k key_k is a class invariant. Ties require a linear
  coincidence mod p (probability <= 2^-249 over a walk); counted as failure.

## 3. Walk: R-adding rho on classes, one dense Fp12 multiplication per step

Parameters: R = 1024 partitions, distinguished-point fraction
theta = 2^-60, M = 2^20 parallel walks started at g^alpha z^beta, fresh
uniform alpha, beta mod r. Precompute M_h = g^{a_h} z^{b_h} (h < R) and the
twisted table T[h][k] = M_h^{p^k} (12,288 entries, ~7.1 MB) with exponent
tables a_h p^k, b_h p^k mod r. Class map F([y]) = [y^ * M_h],
h = K([y]) mod R. Because y^{p^kappa} M_h = (y M_h^{p^{-kappa}})^{p^kappa},
one step is: compute keys, kappa, K (4 Fp mul + bounded small work); if
bits 10..69 of K are all zero, report (K, p^kappa*alpha, p^kappa*beta);
set t = (12 - kappa) mod 12 and y = y * T[h][t] (exactly one dense blst
Fp12 multiplication, the reference operation); update alpha, beta mod r.
The walk is a deterministic function of the class, so merged walks stay
merged and reach the same distinguished class. A collision (equal K,
canonical element confirmed) gives x = (A1 - A2)/(B2 - B1) mod r from the
stored canonical exponents; B1 = B2 with probability 1/r. Verify g^x = z
before output.

Fruitless cycles: exact 2-cycles (h' = h, kappa' = 6) occur at rate
1/(12R) per step; each step compares K with the key two steps earlier and
squares the smaller-K class's stored representative to escape (squaring
maps classes to classes; cost <= 3 units). Longer fruitless cycles need
repeated partition indices, rate <= R^-2 per step; a checkpoint key stored
every 64 steps catches cycles of length <= 64 (waste <= 200 steps per
event). Longer cycles: rate < R^-32, handled by the stuck-walk restart
rule. Squaring steps cannot recur in a fruitless cycle (2^t is not in
<p> mod r for small t).

## 4. Success probability (heuristic H1, disclosed)

H1: the class sequence behaves like a random mapping on the N classes,
with the standard Brent-Pollard R-adding correction N' = N*R/(R-1);
P(no useful collision in S distinct steps) = exp(-S^2/(2N')). Choose
S0 = sqrt(2 N' ln(1/0.6)) = 2^125.6522 and S = 1.03*S0 = 2^125.6949; the
3% is explicit slack for walk non-randomness. Re-derived: S^2/(2N') =
1.03^2 ln(5/3) = 0.5419, model success = 1 - e^-0.5419 = 0.4184 >= 0.40.
Detection delay after the first useful collision exceeds 20/theta steps
with probability e^-20 ~ 2e-9, negligible. Single-walk self-collisions
(~1/M of collisions) and B1 = B2 (1/r) are absorbed by the 0.4184 vs 0.40
margin. Success >= 0.40 over uniform nonidentity Q, uniform z and fresh
coins; comfortably above the 0.39 floor. H1 is the standard heuristic
behind every rho estimate (Pollard 1978; Teske 2001; Brent-Pollard; van
Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999; Wiener-Zuccherato 1998;
Bos-Kleinjung-Lenstra 2010). Toy support (historical experiment.py,
3,096 instances, three small fields with the same tower shape): measured
mean collision time 1.23-1.27 sqrt(N) vs predicted 1.253-1.263;
P[T_c <= S0] 0.375-0.410 vs model 0.40; 2-cycle rate matched 1/(12R).
Toy evidence is valid only at toy scale.

## 5. Cost ledger (units = one dense Fp12 multiplication = 20,000 instructions)

Per-step non-reference instructions (blst pinned prices: Fp mul 383,
Fp12 mul 22,996 rounded to 20,000): 4 Fp multiplications 1,532; 1 Fp add
+ 8 sub/neg <= 360; min of 12 six-limb keys 300; DP test, t, address 40;
alpha/beta updates 100; 2-cycle/checkpoint comparisons 60; loop control 30.
Total 2,422, budgeted 2,800. Step cost c = 1 + 2,800/20,000 = 1.14.
Fruitless-cycle overhead 1 + 3/(12R) + 200/R^2 = 1.000435.

Setup and output (all charged once, included): pairing g < 2^14; z^r check
<= 510; table 1024 x (2 x 255-bit exponentiations + 11 Frobenius maps)
<= 2^20.1; M = 2^20 walk starts g^alpha z^beta <= 2^30.1; distinguished
handling 2^65.7 reports x <= 1,000 instructions + 2 modular multiplications
<= 2^61.4; output [x]P0 <= 510 G1 ops and verification pairing <= 2^14.
Bounded together by 2^62 (dominated over the 2^31 sub-items).

W = c*S*1.000435 + c*E + 2^62, with lockstep detection window
E = M*20/theta = 2^84.33.
Re-derived: W = 1.14 x 2^125.6949 x 1.000435 + 1.14 x 2^84.33 + 2^62
= 2^125.8845 <= 2^125.89. certificate.json carries an exact rational
verification of W^100 <= 2^12589 (no floating point).

Memory (unscored): 2^65.7 reports x 80 bytes at 2/3 hash load + 7.1 MB
table + 2^20 walk states ~ 2^72.6 + 2^23 < 2^73 bytes. Preprocessing
(table, starts) <= 2^31 is inside W. Advice: 0 bytes.

## 6. Heuristics and limitations (complete)

* H1 (random-mapping behaviour of the R-adding class walk) — disclosed,
  standard, toy-supported, unproved; the 0.4184 -> 0.40 margin plus 3%
  step slack is the only headroom.
* K is the min of 12 values, not uniform in [0,p); its low bits and hash
  index are assumed pseudo-random (standard heuristic).
* Key ties counted as failure; total probability <= 2^-249.
* No production-size run. Generic attack; index calculus/TNFS in Fp12 is
  out of scope; W is an upper bound, not a security estimate.
* The historical toy experiment was not re-executed here (candidate code
  is untrusted evidence); certificate.json re-derives every ledger identity
  and constant in independent host integer arithmetic.

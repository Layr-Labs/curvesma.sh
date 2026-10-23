# Fixed-G1 pairing inversion by parallel Pollard rho on order-12 Frobenius classes in GT

Package gt-orbit12-rho-40pct-v5. Claim: log2 W = 12616/100 = 126.16 target-group-operation
equivalents at C = 20000, success >= 0.40, peak storage 2^72 bytes (unscored), advice 2^0 bytes.
Backport of historical branch submission `1c0ec6c4-10ab-481f-bf1d-13f7c10eedae`
(`submissions/bls12-381-fixed-g1-invert-g2/gt-class-rho`, claimed 126.16, author unknown,
source never published to a public track list per the dependency-publication audit). Every
identity below was independently re-derived in host python3 exact integer / float arithmetic
against the current protected profile (see certificate.json); nothing was executed from the
source branch.

## Reduction (exact, no heuristic)

Input P in G1 nonidentity, z in GT; output Q in G2 with e(P,Q)=z.

1. z = 1: return infinity.
2. g = e(P,Q0), one reduced pairing (charged). Nondegeneracy of the reduced pairing on
   prime-order G1 x G2 and P != 0 give g != 1; GT has prime order r, so g generates GT.
   Hence z = g^x for a unique x in [0,r).
3. Recover x = log_g z by the class-level rho below; return Q = [x]Q0; check g^x = z in GT
   and recompute e(P,Q) (charged). Correctness is bilinearity: e(P,[x]Q0) = e(P,Q0)^x = z.

Profile facts verified in python3 (certificate.json): ord_r(p) = 12 exactly (p^i != 1 mod r
for 1 <= i < 12, p^12 = 1), p^6 = -1 mod r, r | p^12-1, r | (r-1) after /12 split, a random
h^(p^12-1)/r has order r, Frobenius^6 on GT equals inversion (verified as s^(p^6) = Fp6
conjugate w -> -w), and a random GT element has a full 12-element Frobenius orbit.

## Frobenius classes and the cheap class invariant

p has order 12 mod r, so phi: s -> s^p acts on GT as exponentiation by p mod r, an
automorphism of order 12 with no fixed point on GT \ {1} (p^i != 1 mod r for 0 < i < 12).
GT \ {1} splits into (r-1)/12 classes {s^(p^i)} of exactly 12 elements; with the identity
alone, N' = (r-1)/12 + 1, log2 N' = 251.27213 (log2 r = 254.85709).

Field model (profile): Fp12 = Fp[w]/(w^12 - 2w^6 + 2) with u = w^6 - 1, i.e. tower
Fp2[w]/(w^6 - xi), xi = 1 + u, matching the profile's coefficient list [c0..c11]. Then
phi^i(sum c_k w^k) = sum conj^i(c_k) gamma_i^k w^k with gamma_i = xi^((p^i-1)/6) in Fp2,
12 public constants. Define v_i(s) = Re(conj^i(c_1) gamma_i), the Fp-real part of
coefficient 1 of phi^i(s) (in the profile's flat basis this is c1 + c7). Verified on random
order-r elements: (i) v_i(phi(s)) = v_{i+1}(s) for all i mod 12, so the multiset {v_i}
depends only on the class; (ii) v_{i+6} = -v_i exactly (phi^6 = inversion = w -> -w
conjugation on GT), so only v_0..v_5 must be computed and the other six values are p - v_i;
(iii) the 6 absolute values are distinct on all test elements (collisions probability ~2^-380
per pair). Computing v_0..v_5 costs 10 Fp multiplications (2 per v_i) plus adds; a 12-value
argmin over 6-limb values. A class label is vmin = min over the 12 signed values: a full
Frobenius is NEVER evaluated.

Collision logic: if class(g^a z^b) = class(g^a' z^b') then g^a z^b = (g^a' z^b')^(p^k), so
x = (lambda^k a' - a)/(b - lambda^k b') mod r whenever b != lambda^k b' (failure prob ~1/r;
detected at verification).

## Walk

Per-processor table T[j] = g^(a_j), a_j random in [0,r), j < R = 2^16, plus conjugates
T[j][i] = phi^{-i}(T[j]) (12 R Fp12 elements = 2^28.75 B per processor), built by 11
Frobenius maps per entry (cheap: 6 Fp2 maps each). State s in GT. Step: compute
v_0..v_5, vmin and its index i* in 0..11, j = bits 200..215 of vmin, s <- s * T[j][i*].
Since phi^{i*}(s T[j][i*]) = phi^{i*}(s) t_j, the induced class map f(class) =
class(canon(s) t_j) is well defined on classes. One step = one dense Fp12 multiplication
(1.0 reference ops) + charged instruction overhead (below). Start: s0 = g^a0 z^b0 from the
seed (~255 squarings + ~128 multiplications in the cyclotomic group, charged 512 ops per
start). Distinguished point: vmin low 60 bits zero (theta = 2^-60; the table index bits
200..215 are disjoint). Restart after each distinguished point or after 20/theta = 2^64.3
steps.

Fruitless cycles and the class-only escape. The step j(m)+1 after a phi^6 conjugate
undoes the prior step: 2-cycles occur at exact rate 1/(12R) = 2^-19.585 per step. Each
walker keeps the last 24 (label, a, b) states; when the new label matches a state w <= 24
back, the walk has completed a cycle of f of length w; the walker takes the cycle member m
with the smallest vmin, restores m's stored first-visit state and escapes with table index
j(m)+1. esc(m) is a function of the class m alone, so the modified map f' is still a
function of classes, and every walker that reaches a class realises the identical future
class sequence (van Oorschot-Wiener distinguished-point correctness is preserved: merged
walks stay merged, including walks that entered a cycle at different members). Cycle rates
(H2): rho_2 = 1/(12R) exact; rho_L <= L! R^-ceil(L/2) for L >= 3 (a length-L cycle needs a
vanishing sum of table exponents, hence >= ceil(L/2) index coincidences). Re-derived tail
L > 24 total rate <= 2^-119.5 per step (dominant L=26 term; matches the historical
2^-119.5): a trapped walk is cut at 2^64.3 steps, expected loss <= 2^-55.2 fraction of
steps (absolute 2^71).

## Recovery

Store distinguished points by vmin (48-byte label + 16-byte seed/length = 64 B). A repeat
label is a class collision (false match ~2^-300). Rerun the two walks from their seeds
tracking (a,b) mod r exactly; x = (a'-a)/(b-b') mod r (failure 1-2^-254); confirm by
g^x = z; return [x]Q0 plus final pairing check. At most 2 x 2^64.3 rerun steps at 4x cost.

## Success (heuristic, disclosed)

H1 (random mapping): the class-level r-adding walk with R = 2^16 behaves like a random
mapping on the N' classes (van Oorschot-Wiener 1999; Teske 2001: R >= 16 within a few
percent of the ideal collision time, deviation shrinking with R; 3% margin charged).
Budget k = sqrt(2 N' ln(1/0.61)) = sqrt(0.9886 N') = 2^125.62779 gives P(collision) >=
0.39 under H1. The attack runs K = 1.03 k = 2^125.67043 logical steps; re-derived delivered
P = 1 - exp(-K^2 / 2N') = 0.4081. Re-derived ratio K / k(0.40) = 1.0132: the margin over
the 0.40 threshold is 1.3%, and all collision-relevant losses (detection delay 2^-41.3,
cycle waste 2^-17.25 wasted-step fraction, undetected cycles 2^-55.2) are below 2^-40
relative, so **success >= 0.40 is carried**, as the historical margin supports. Success
does not depend on P beyond g != 1 (exact for all admissible P). Identity target 2^-255
handled directly.

## Cost ledger (independently re-derived; every identity in certificate.json)

Step: 1 dense Fp12 mul + 10 Fp mul (10 x 383 = 3830) + 144 compare + 8 DP test + 80
window compare + 150 state copy + 42 addressing/IO = 4254 instructions (historical ledger
counted 4454 with looser add charges); charged 8000 (1.88x headroom for cache/table
transfers): c = 1.4 reference ops.

- Logical steps: K = 1.03 k, log2 K = 125.67043. K c: log2 = 126.15586.
- Cycle waste: omega = sum_{L=2..24} rho_L (L+3) = 2^-17.25 (re-derived; dominant L=2 term
  2^-17.26); total K c (1+omega) = 2^126.15587.
- Starts: (K theta + M) 512 = 2^74.7 ops.
- DP store/compare/comm: 2^64.7.
- Detection delay: M walkers x 20/theta steps x c: 2^84.8 (relative 2^-41.3).
- Recovery reruns: 2 x 2^64.3 x 1.4 x 4: 2^68.8.
- Undetected long cycles: fraction 2^-55.2 of the steps: absolute 2^71.
- Preprocessing (CORRECTED): per table 2^16 random-exponent generations ~255 squarings +
  128 mults in GT (380 Fp12-equivalents each) + 11 x 2^16 Frobenius (6 Fp2 maps x 383
  instr.) = 2^24.6 per processor; x 2^20 processors = 2^44.6, plus one pairing 2^14:
  **2^45** absolute (historical 2^38.2 undercharged random-exponent generation; still
  relative 2^-81.6, immaterial to the score).
- Output: pairing e(P,Q0) 2^14 + final pairing 2^14 + [x]Q0 (255 doublings + ~128 adds
  x ~20 Fp2 ops x 3 x 383 instr / 20000 = 2^11.6) + g^x = z 2^9: total 2^15.1.

Sum: W = 2^126.15587 (1 + 2^-41.3 + 2^-51.5 + 2^-57.3 + 2^-55.2 + 2^-61.5 + 2^-81.6 +
2^-111) ; independent re-derivation of the full sum: **log2 W = 126.15587 <= 126.16.**
Claimed 126.16 confirmed under the current cost model with the corrected preprocessing.

## Memory, preprocessing, advice

Stored distinguished points K theta = 2^65.67 x 64 B = 2^71.67 B; processor tables
2^20 x 2^28.75 B = 2^48.75 B; walker windows 2^20 x 24 x 576 B = 2^35.0 B. Peak
2^71.9 < **2^72 bytes** (confirmed; memory unscored per rules.md). Preprocessing 2^46
(included in W; depends on g, not reusable across instances). Advice 2^0 bytes.

## Heuristics, limits, disclosure

- H1 class-level random mapping (above), 3% margin; K delivers 0.408 >= 0.40 by 1.3%, and
  the 0.40 floor is carried because the margin survives every loss term by >= 2^40.
- H2 cycle-rate bound rho_L <= L! R^-ceil(L/2), exact for L = 2.
- Analytic construction; no full-size execution; instruction counts are charged bounds, not
  measurements of this kernel. Memory is not physically realisable (declared; unscored).
  No novelty claim, no practical-attack claim, no quantum steps. The historical
  proof.md's two inconsistent table-index bit ranges (40..55 vs 200..215) were resolved in
  favour of bits 200..215 (disjoint from the 60 DP bits); immaterial to cost.
- Comparison: sqrt(12) = 2^1.79 saving over plain GT rho; vs the currently selected v3
  (699/5 = 139.8) a 13.6-bit reduction, bought with one disclosed random-mapping heuristic
  in exchange for v3's determinism.

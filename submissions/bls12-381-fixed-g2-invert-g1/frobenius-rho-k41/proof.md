# Mixed squaring / multiplication Frobenius-class Pollard rho in GT (success target 0.39)

Claim: time_log2 = **125.767** (total work W <= 2^125.767 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory <= 2^73 bytes,
preprocessing <= 2^31 units (included), no advice. Heuristics are listed in section 6.
The claim charges a measured 10% step inflation for the mixed map as extra steps
(section 4), so the walk-step cap is 1.10 x 1.000435 x S = 2^125.8093 steps.

Design provenance: the reduction, the four-multiplication class-key formula, the
twisted multiplier table, the class walk, the fruitless-cycle escapes and the record
format are those of public submission 3da91365-2d58-40e5-ac7d-415e3061f7d8
(Frobenius-class rho, 125.89), as tightened in the 125.861 package. This package
adds one thing: a fraction f of the walk steps are **cyclotomic squarings** instead
of dense multiplications. Nothing here reuses an earlier score as evidence; every
bound below is derived or computed in this package.

## 1. Reduction

Input: nonidentity Q in G2 and z in GT. Compute g = e(P0, Q) with the pinned
convention. Q has prime order r and the pairing is nondegenerate, so g != 1 and g
generates GT (prime order r). Every valid z equals g^x for a unique x in [0, r).
Recover x by the walk below, output P = [x]P0, check e(P, Q) = z; then e(P, Q) =
g^x = z. If z = 1, output P = infinity directly. The answer is for the supplied Q.

## 2. Arithmetic facts used (checked numerically in `experiment.py`, part A)

Represent Fp12 as blst does: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi),
xi = 1+u, basis w^0..w^5. Conversion to the profile's coefficient list is
c_m = a_m - b_m, c_{m+6} = b_m (free).

* Frobenius: (sum c_m w^m)^p = sum conj(c_m) d_m w^m, d_m = xi^{m(p-1)/6}. p has
  multiplicative order exactly 12 modulo r, so {y -> y^{p^k}} has order 12 and
  contains inversion (p^6 = -1 mod r). Every y != 1 lies in a class [y] of exactly
  12 conjugates; the number of classes is N = (r-1)/12, log2 N = 251.2721.
* Class key: let zeta = xi^{(p^2-1)/6}, a primitive 6th root of unity in Fp. With
  c_1 = a + b'u the w^1 coefficient of y, the Re(w^1) of the conjugate y^{p^k} is
  key_k = a*zeta^j (k = 2j) or B*zeta^j (k = 2j+1), B = a*Re(d_1) + b'*Im(d_1).
  Since zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta}, all 12 keys cost **4 Fp
  multiplications** plus 8 subtractions/negations; no Frobenius of the element is
  ever computed. kappa(y) = argmin_k key_k and K([y]) = min_k key_k is a class
  invariant; y^{p^kappa} is the same element for all 12 members of the class.
* **GT lies in the cyclotomic subgroup.** The embedding degree is 12, i.e. r
  divides Phi_12(p) = p^4 - p^2 + 1 (checked in part A; r is much smaller than
  Phi_12(p), so this is a divisibility, not an equality), and every y of order
  dividing r satisfies y^{Phi_12(p)} = 1, so GT is contained in
  G_{Phi_12} = ker(Phi_12(Frobenius)), the cyclotomic subgroup. Granger-Scott
  squaring therefore applies to every walk element.
* **Squaring is a class map with no twist.** (y^{p^k})^2 = (y^2)^{p^k}, so
  y -> y^2 maps classes to classes: [y] -> [y^2]. On exponent pairs it is
  (alpha, beta) -> (2 alpha, 2 beta) mod r, because (g^alpha z^beta)^2 =
  g^{2 alpha} z^{2 beta}. It is **injective on classes**: y1^2 = y2^2 implies
  (y1/y2)^2 = 1, and y1/y2 lies in GT, whose order r is odd, so y1/y2 = 1 and
  y1 = y2. A squaring step therefore never merges two distinct classes and the walk
  stays a bijective image of itself. No partition and no table entry is needed to
  take this step beyond the key that every step already needs (section 3.1).

## 3. Walk: mixed steps on classes

Parameters: R = 1024 partitions, distinguished fraction theta = 2^-60, M = 2^20
parallel walks started at fresh g^alpha z^beta, and a squaring fraction f.

Precompute M_h = g^{a_h} z^{b_h} (h < R) and T[h][k] = M_h^{p^k} (12,288 entries,
7.1 MB) with exponent tables a_h p^k, b_h p^k mod r. A walk holds an arbitrary
member y of its class with (alpha, beta), y = g^alpha z^beta.

One step:

```
compute keys, kappa, K                          (4 Fp mul + ~700 other instructions)
if bits 10..69 of K are all zero (distinguished):
    if the report store already holds C_rep = 2^66 records: abort the run, return failure
    report (K, p^kappa alpha, p^kappa beta)
if bits 198..207 of K, read as a 10-bit value, are below 1024 f:
    y = y^2 ; alpha = 2 alpha mod r ; beta = 2 beta mod r        (cyclotomic squaring step)
else:
    h = K mod R (bits 0..9) ; t = (12 - kappa) mod 12
    y = y * T[h][t] ; alpha += a_h p^t ; beta += b_h p^t mod r   (dense multiplication step)
```

A report is 80 bytes: bits 70..197 of K (16 bytes, the 128 bits directly above the
distinguished bits, which are zero in every report) plus the canonical exponents
A = p^kappa alpha and B = p^kappa beta mod r (32 bytes each); the store is a hash
table keyed on those 128 bits at load factor 2/3. A key match is confirmed by
recomputing the canonical element g^A z^B for both records (<= 1,021 units each) and
comparing in full; a match between different classes has probability <= 2^-128 per
pair, and with at most C_rep^2/2 = 2^131 pairs the expected number of false matches is
<= 2^3, so by Markov's inequality more than 2^40 occur with probability < 2^-37. Each
false match costs only its confirmation and is discarded, so 2^40 of them cost
< 2^51 units, inside the 2^62 term of section 5.

**Distinguished-report store (hard cap and abort rule).** Every report enters a
distinguished-point store of hard capacity C_rep = 2^66 records. The report count is
enforced by this cap, not by the walk-step cap: when a report arrives and the store
already holds C_rep records, the run aborts and returns failure (a success-loss event;
section 4 bounds its probability). The walk-step cap is the one of section 5. The run
makes at most N_steps = 1.10 * 1.000435 * S = 2^125.8093 steps (section 5), and under
H1 each step's class is distinguished with probability theta = 2^-60, so the report
count X is stochastically dominated by Binomial(N_steps, theta), whose mean is
mu = N_steps * theta = 2^65.809. The cap is C_rep = 2^66 = a * mu with a = 1.1413,
and the Chernoff bound P(X >= a mu) <= exp(-mu (a ln a - a + 1)) with
a ln a - a + 1 = 0.00955 gives P(X > C_rep) <= exp(-0.00955 * mu) = exp(-2^59.1), far
below 2^-1000. The distinguished test runs on both step kinds (the key is computed
every step, section 3.1), so f does not change the report count.

### 3.1 Why the step type must come from the class key

The walk must be a function of the class alone, or two walks that reach the same
class take different next steps and stop being mergeable; the parallel
distinguished-point algorithm then loses its square-root behaviour (each walk would
have to collide with another at a stored point exactly, not merely merge and meet
the next distinguished class). A phase counter or a per-walk counter would have
exactly that defect: with a fixed multiply/square pattern of period p, the visited
positions split into p phase classes and only same-phase collisions merge, which
costs a factor p in steps. Taking the step type from bits 198..207 of the raw class
key keeps the walk a function of the class while leaving the partition (bits 0..9),
the distinguished test (bits 10..69) and the record key (bits 70..197) untouched.
K is about 2^377, so those bits are far inside its range and, like the bits the
frontier already uses, behave as independent uniform bits of the class.

The price is that the key is computed on squaring steps too. It is worth it: the
squaring itself is about half a dense multiplication, while the key block is 0.135
units, and dropping the key would cost a factor p in steps (p >= 2), i.e. at least
one bit, far more than it saves.

### 3.2 Fruitless cycles

Each step compares K with the key two steps earlier; on equality the two classes
{C1, C2} form a 2-cycle and the walk squares the representative of the one with
smaller K. A checkpoint key is stored every 64 steps; a repeat means a cycle of
length L <= 64, which is traversed once more to find its minimum-K class, which is
squared. Both escapes depend only on the set of classes in the cycle, not on when
the walk entered it, and every step type is a class function, so merged walks leave
a cycle at the same class and stay merged. Detection of an L-cycle happens within
64 + L steps of entering it, so the waste per event is <= 64 + 2L + 1 <= 193 <= 200
steps; a 2-cycle event wastes <= 3 steps. A walk with no distinguished class for
40/theta steps is abandoned and restarted at a fresh g^alpha z^beta.

The 2-cycle rate does not increase. A 2-cycle C -> C' -> C needs, in the two step
types: (multiply, multiply) M_{h'} = M_h^{-1} and the same partition and twist,
rate 1/(12R) as before, now multiplied by (1-f)^2; (squaring, squaring) would need
C^4 = C, i.e. C^3 = 1, impossible since gcd(3, r) = 1; (multiply, squaring) needs
(C M_h)^2 = C and (squaring, multiply) needs C^2 M_h = C, each a coincidence of
probability about 1/r. Longer fruitless cycles: the derivation of the frontier
package (every index in the cycle occurs at least twice and each block of twists
sums to zero in Z[x]/Phi_12, probability <= 1/12 per block) applies unchanged to the
multiplication steps of a cycle, and a pure squaring run cannot close a cycle
because 2^t mod r is not a 12th root of unity. Charging the frontier's rate R^-2 and
its 200-step waste is therefore conservative.

## 4. Success probability

Heuristic H1 (unchanged): the class sequence behaves like a random mapping on the N
classes, with the Brent-Pollard correction N' = N R/(R-1). Then the probability that
no useful class collision occurs among the first S distinct-step points is
exp(-S^2/(2N')). With S0 = sqrt(2 N' ln(1/0.61)) = 2^125.6285 and S = 1.03 S0 =
2^125.6711 the model gives success 1 - 0.61^(1.03^2) = 0.4081 at the budget S.
The algorithm's walk-step cap is N_steps = 1.10 x 1.000435 x S = 1.1335 S0 =
2^125.8093 steps (section 5), a fixed budget, not an expectation; the extra 1.10 is
the mixing inflation charged below.

**Heuristic H2 (new): the mixed map is random-mapping-like.** The step map is
F(C) = P(C) on the squaring classes and F(C) = C M_{K(C) mod R} otherwise, with P the
class squaring map. P is a fixed bijection of the N classes that is independent of
the multiplier table, and composing a random map with a fixed bijection gives a
uniformly distributed map, so F is random-mapping-like for every f < 1; the
degenerate case f = 1 (pure P, a single cycle of length ord_r(2)) is excluded. Three
structural facts support the heuristic: P is injective on classes (section 2), so a
squaring step never merges two classes; every step type is a class function
(section 3.1), so merging is exact; and the distinguished test runs on every step,
so a merged pair is detected within the same 20/theta window as before.

H2 is measured, not assumed, and the measurement is **charged**. Part B runs the
complete production walk, with the production partition and distinguished predicates
and the class-function step-type rule, in three fields q near 2^20 at R = 64, 256,
1024, with a multiplication-only control (f = 0) and the mixed walk on the same seeds,
200 instances per cell. Section 7 gives the table. In short: the control's fraction
within 1.03 S0 is 0.418 pooled over the three fields, the mixed walk's is 0.372, and
the corresponding effective-budget ratio is 1.080, i.e. the mixed walk needs about
8% more steps for the same success. That shortfall is charged as **10% more steps**:
the walk-step cap is 1.10 x 1.000435 x S = 1.1335 S0 = 2^125.8093 steps, and the
per-step costs of section 5 are unchanged. At that cap the measured 1.080 ratio puts
the mixed walk's effective budget at 1.1335/1.080 = 1.050 S0, where H1 gives
1 - 0.61^(1.050^2) = 0.4199. The measured inflation is 1.099, 1.086 and 1.080 at
f = 1/4, 1/3 and 1/2, so 1.10 covers all three; the standard error of the pooled ratio
is about 0.05.

Why the mixed map loses a few percent is visible in the structure: a squaring step
is a bijection that carries no fresh choice from the multiplier table, so the
positions immediately after squarings are deterministic images of earlier positions
and the short-distance pairs that a random map would spend are partly wasted. The
loss is small because the step type is class-determined, so no class is duplicated
and the long-range pair count is unchanged.

Detection and the other losses are the frontier's, and their sum stays below the
margin: after the first useful collision the merged trajectory reaches a
distinguished class within 20/theta further steps except with probability e^-20;
charging all M walks over that window costs E = M*20/theta = 2^84.3 steps. The
losses, each a probability bound: detection tail e^-20 < 2^-28; a trailing walk
restarted before reaching the distinguished class (< 2^-52); self-collisions
(2^-19); B1 = B2 (1/r); key ties (< 2^-249); long fruitless cycles (< 2^-100);
more than 2^30 restarts (< 2^-21); more than 2^40 false truncated-key matches
(< 2^-37); cycle-handling overrun (< 2^-40); report-store overflow, i.e. more than
C_rep = 2^66 reports before success (< 2^-1000, section 3); a walk element equal
to 1 (<= 2^20/r). Their sum is < 2^-19 + 2^-21 + 2^-28 + 3 x 2^-37 + 2^-1000 < 2^-18
< 0.4199 - 0.39 = 0.0299, so success >= 0.4199 - 2^-18 > 0.4198 >= 0.39.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Both step kinds compute the same key block. Instruction counts are for scalar
AArch64 with the measured blst_fp_mul = 383 instructions, and are read off the
pinned blst source f7a8806 (`src/fp12_tower.c`, `__FP2x2__` path, which is what the
organizer's `-O3` build with upstream AArch64 assembly compiles).

| Work, every step | Instructions |
| --- | ---: |
| 4 Fp multiplications (a*zeta, two products for B, B*zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, <= 40 each | 360 |
| min of 12 six-limb keys (top-limb compare, csel; rare full compare) | 300 |
| hash, distinguished test, t = (12-kappa) mod 12 | 40 |
| 2-cycle and checkpoint key comparisons, key history stores | 60 |
| loop control | 30 |
| **key block, counted** | **2,322** |
| **key block, budgeted (16% margin, as the frontier's 2,422 -> 2,800)** | **2,700** |

A multiplication step adds `h = K mod R`, the table address (10) and the two
modular additions of the exponent update (100); with the frontier's 2,800-instruction
budget its cost is

  c_mul = 1 + 2,800/20,000 = **1.14000 units**.

A squaring step adds two doublings mod r (60) and one Granger-Scott cyclotomic
squaring. The latter is read off `cyclotomic_sqr_fp12`, which is 3 x `sqr_fp4` =
**9 x `sqr_fp2x2` + 12 x `redc_mont_384`** plus 42 Fp adds and 18 double-width
add/sub. The dense reference multiplication `mul_fp12` is 24 x `mul_fp2x2` +
12 x `redc_mont_384` plus about 30 Fp adds and 40 double-width add/sub, and the
organizer measures it at 22,996 instructions. Charging `sqr_fp2x2` at the same rate
as `mul_fp2x2` (a square is normally cheaper) and taking redc_mont_384 = 200
instructions, which is the middle of its plausible 120..300 range, gives:

| Work, squaring step | Instructions |
| --- | ---: |
| key block, budgeted | 2,700 |
| 2 exponent doublings mod r | 60 |
| Granger-Scott squaring: 9 x 741.5 + 12 x 200 + 42 x 40 + 18 x 40 = 11,474 counted | 13,263 budgeted |
| **total** | **16,023** |

  c_sqr = 16,023/20,000 = **0.80115 units**.

The same decomposition at redc_mont_384 = 120 and 300 gives 12,570 and 14,130
budgeted, i.e. c_sqr between 0.767 and 0.845; the middle bracket is used and the
sensitivity is stated in section 6. So

  c(f) = (1-f) * 1.14000 + f * 0.80115,

and the mixed ledger is charged the 10% step inflation as extra steps (section 4):

| f | c(f) | log2 W at the capped budget | claim |
| --- | ---: | ---: | ---: |
| 1/4 | 1.05529 | 125.88691 | 125.887 |
| 1/3 | 1.02706 | 125.84778 | 125.848 |
| 1/2 | 0.97058 | 125.76619 | **125.767** |

Fruitless-cycle factor: 1 + 3/(12R) + 200/R^2 = 1.000435, unchanged and now
conservative since 2-cycles need both steps to be multiplications. The walk-step cap
of the algorithm is 1.10 x 1.000435 x S = 1.1335 S0 = 2^125.8093 steps total, cycle
handling included, a hard cap; the 1.10 is the mixing inflation of section 4,
charged as extra steps rather than as a per-step cost.

Setup and output, all <= 2^31 units together: pairing g (< 2^14), input check
z^r = 1 (<= 510), table (1024 multipliers, <= 2^20.02), M = 2^20 walk starts
(<= 2^30.0), distinguished-point handling (charged for the full hard cap C_rep = 2^66
reports, not the mean, x <= 1,000 instructions: hash-table insert, 2 modular
multiplications by p^kappa), output [x]P0 and verification pairing (<= 2^14). The
report-processing term is 2^66 x <= 1,000 instructions = 2^75.97 instructions
= 2^61.68 units (C = 20,000 = 2^14.2877), below 2^62, so setup, report processing and
output together stay within 2^31 + 2^62.

With f = 1/2 and the step cap of section 4 (the values chosen in section 7):

  W = c(1/2) * (N_steps + E) + 2^62
    = 0.97058 * (1.10 * 1.000435 * 2^125.6711 + 2^84.3) + 2^62
    = 0.97058 * (2^125.80927 + 2^84.3) + 2^62
    = 2^125.76619

The other two terms add < 2^-41. Submitted bound: **time_log2 = 125.767** =
125767/1000 >= 125.76619 (part C of experiment.py recomputes log2 W and checks the
claim against it).

Memory: the distinguished-report store holds at most C_rep = 2^66 records (the hard
cap of section 3), each 80 bytes (16-byte truncated key, bits 70..197 of K, + two
32-byte exponents) at hash-table load 2/3: 2^66 x 80 x 1.5 = 2^66 x 120 = 2^72.907
bytes, plus the 7.1 MB table (2^22.8) and 2^20 walk states (~1 KB each, 2^30 bytes):
peak < 2^72.91 < 2^73 bytes. This bound is by construction: a run that would need a
C_rep + 1-st record aborts (section 3) instead of exceeding it. Memory is reported,
not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of an R-adding walk on Frobenius classes) is the
  standard heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard;
  van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999; Bos-Kleinjung-Lenstra 2010
  for automorphism classes and fruitless cycles). It is supported by toy runs, not
  proved.
* H2 (that the mixed map is random-mapping-like up to a constant) is new here. The
  toy runs of section 7 measure a 1.08 effective-budget ratio against the control in
  the same harness, and the ledger charges 1.10, so the discrepancy is bounded on
  the conservative side; the toy is at 2^20 scale and the pooled ratio has a
  standard error of about 0.05, so a larger true inflation is not excluded.
* The distinguished-point and partition bits, and now the step-type bits, are taken
  from K. K is the minimum of 12 values and is not uniform, but the bits used are far
  below its 377-bit size, and the toy runs measure the distinguished fraction at its
  nominal value with the same rule.
* c_sqr is analytic. The bracket redc_mont_384 = 120..300 moves it between 0.766 and
  0.845 units, which moves log2 W between 125.740 and 125.798 at f = 1/2. The claimed
  price is 200, the middle of the range: there c_sqr = 0.80117 and log2 W = 2^125.76619,
  so the submitted bound 125.767 holds for redc_mont_384 <= about 202. At 300 the bound
  would be 125.798, above the claim; 300 is outside the plausible range because one
  base-field multiplication (`mul_mont_384`) measures only 383 instructions in total.
  The claim is conditional on this per-reduction price.
* No production run. Toy experiments are evidence only at their own scale.
* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): random GT elements; y^r = 1; Frobenius formula = y^p; the
4-multiplication key formula equals the keys of all 12 explicit conjugates and they
are distinct; kappa selects one identical canonical element from all 12 conjugates;
p^6 is inversion; ord_r(p) = 12; **r divides p^4 - p^2 + 1** (GT is inside the
cyclotomic subgroup) and **squaring commutes with Frobenius, y^2 never lands on its
own class and -1 is not in GT (r odd)**, so y -> y^2 is an injective map on classes
that doubles the exponent pair.

Part B (toy, same tower shape): q = 7 mod 12 near 2^20, r' a prime factor of
q^4 - q^2 + 1 with ord_r'(q) = 12, N = (r'-1)/12. It runs the production walk: four
multiplication keys, kappa, twisted table, partition h = K mod R and distinguished
test on the raw key, reports stored by key and confirmed by the canonical element,
2-cycle escape, checkpoint escape, stuck-walk rule, collision solving and
verification. The mixed mode adds the class-function squaring step. Toy-only
differences: the 60 distinguished bits become dp_bits and the 10 step-type bits
become 6; the toy key has only about 13 significant bits, so the step type is taken
from a 32-bit mix of K rather than a raw bit slice (production's slice is uniform;
this one is measured and reported as the realised squaring fraction), and rare
toy-key ties are broken by the full conjugate.

Part C recomputes, without randomness, the zero-sum probabilities, the fruitless-cycle
bound, the per-step ledger of section 5, log2 W and the claim check.

Part B, 200 random instances per cell, seeds 401 / 501 / 601 (one seed per field),
4 walks, production partition and distinguished predicates on the raw class key, the
step type taken from the class key (the toy uses a 6-bit mix of its short key; the
realised fraction is measured). Control f = 0 is the multiplication-only walk:

| q, r', R, dp_bits | f | trials | mean T_c/sqrt N (model 1.25) | P[T_c<=S0] | P[T_c<=1.03 S0] | 2-cycles/step | long cycles | distinguished |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 0 | 200 | 1.2258 | 0.445 | 0.465 | 1.33e-03 | 6 | 0.0627 |
| 1048963, 48674953, 64, 4 | 1/4 | 200 | 1.3297 | 0.340 | 0.360 | 7.35e-04 | 5 | 0.0630 |
| 1048963, 48674953, 64, 4 | 1/3 | 200 | 1.2695 | 0.370 | 0.405 | 7.28e-04 | 3 | 0.0623 |
| 1048963, 48674953, 64, 4 | 1/2 | 200 | 1.3870 | 0.335 | 0.360 | 3.72e-04 | 1 | 0.0628 |
| 1049611, 64668529, 256, 4 | 0 | 200 | 1.3049 | 0.360 | 0.370 | 3.56e-04 | 0 | 0.0638 |
| 1049611, 64668529, 256, 4 | 1/4 | 200 | 1.3878 | 0.325 | 0.335 | 1.80e-04 | 12 | 0.0638 |
| 1049611, 64668529, 256, 4 | 1/3 | 200 | 1.3111 | 0.325 | 0.335 | 1.28e-04 | 1 | 0.0638 |
| 1049611, 64668529, 256, 4 | 1/2 | 200 | 1.4566 | 0.335 | 0.340 | 7.65e-05 | 2 | 0.0641 |
| 1054003, 35456053, 1024, 3 | 0 | 200 | 1.2452 | 0.390 | 0.420 | 7.62e-05 | 1 | 0.1291 |
| 1054003, 35456053, 1024, 3 | 1/4 | 200 | 1.2755 | 0.385 | 0.390 | 5.19e-05 | 0 | 0.1304 |
| 1054003, 35456053, 1024, 3 | 1/3 | 200 | 1.3704 | 0.330 | 0.365 | 5.03e-05 | 0 | 0.1300 |
| 1054003, 35456053, 1024, 3 | 1/2 | 200 | 1.3164 | 0.400 | 0.415 | 2.18e-05 | 1 | 0.1303 |

Pooled over the three fields (600 instances per f), the fraction within 1.03 S0 is
0.418 for the control, 0.362 at f = 1/4, 0.368 at f = 1/3 and 0.372 at f = 1/2, i.e.
effective-budget ratios 1.099, 1.086 and 1.080. The realised squaring fractions are
0.249, 0.328 and 0.499; the two-cycle rate falls with f as (1-f)^2 predicts
(e.g. at R = 64: 1.33e-3, 7.35e-4, 7.28e-4, 3.72e-4); the distinguished fraction
matches its nominal value at every cell; long-cycle counts are 0 to 12 per cell and
are dominated by a few instances that re-enter a cycle (the ledger charges the
frontier's R^-2 rate with 200 steps per event, which at the production R = 1024 is
0.0003 of the steps and this table's R = 1024 cells measure at most one event); no
walk hit the stuck-walk rule.

**Choice of f.** f = 1/2 is the largest fraction screened and the smallest cost, and
the measured inflation is lowest there (1.080), so it gives the smallest charged
work of the three: log2 W = 125.76619 at f = 1/2, against 125.84778 at f = 1/3 and
125.88691 at f = 1/4, versus 125.861 for the multiplication-only frontier. The gate
for this package was 0.36 success within 1.03 S0 at f = 1/2; the measurement is 0.372,
and the 0.018 shortfall against 0.39 is charged as 10% more steps (section 4) rather
than absorbed.



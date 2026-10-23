# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV automorphism)

Claim: time_log2 = **126.64**. The algorithm below recovers x from T = [x]P0
with success probability at least 0.40, using total work
W <= 2^126.64 batched-affine-G1-addition equivalents (C = 2,000 scalar AArch64
instructions per unit). The phase-1 work is capped by a hard counter; every
other term has a hard cap; the only probabilistic statements are the success
probability and negligible-probability failure events, which are counted as
failures. Peak memory is below 2^72 bytes; preprocessing (the multiplier table)
is below 2^23 units and is included in W; nonuniform advice is 0 bytes. The
success bound rests on heuristic H1 (section 4): the random-mapping law for the
class walk with an explicit allowance of up to 5% slowdown, supported by
toy-curve experiments; no production-size run is claimed.

## 1. Problem and output

Input: T in G1. If T is the point at infinity, output x = 0 and stop. Otherwise
x is in [1, r-1]; the algorithm outputs the scalar x (64 hex digits) after
checking [x]P0 = T. It never outputs a point or a relation.

## 2. Algebra used (checked on BLS12-381 in `experiment.py`, part A)

p is the 381-bit field prime, r the 255-bit subgroup order (log2 r = 254.85709),
G1 = <P0> of prime order r on E: y^2 = x^3 + 4.

* p = 1 mod 3. beta = 0x1a0111ea...00000000aaac (full value in `experiment.py`)
  is a primitive cube root of unity in Fp, and phi(x, y) = (beta*x, y) is an
  automorphism of E. lambda = 0xac45a4010001a40200000000ffffffff is a primitive
  cube root of unity mod r and phi(P0) = [lambda]P0, hence phi = [lambda] on G1
  (checked numerically; the other root pairs with the other beta).
* The units u_t := (-1)^s lambda^j (t = 3s + j, s in {0,1}, j in {0,1,2}) form
  a group of order 6 acting on G1 by [u_t](x, y) = (beta^j x, (-1)^s y).
* No nonidentity point of G1 is fixed by a nontrivial unit ([u]P = P with
  P != O forces u = 1 mod r). Hence G1 \ {O} splits into exactly N = (r-1)/6
  orbits ("classes") of size 6, log2 N = 252.27213. (The points (0, +-2) of E are
  3-torsion and not in G1.)
* **Class key.** For (x, y) let K = |y| = min(y, p - y) (integers in [0, p)).
  The 6 members of a class are (beta^j x, +-y), so K is constant on a class.
  Conversely two G1 points with the same |y| have y' = +-y and
  x'^3 = y'^2 - 4 = x^3, so x' = beta^j x: same class. K is a complete class
  invariant (distinct classes have distinct keys).
* **Canonical member.** For P = (x, y) in class C let j minimise the integer
  beta^j x mod p over j in {0,1,2} (three distinct values since x != 0) and
  s = 1 iff p - y < y. The canonical member c(C) = (beta^j x, |y|) = [u_t]P,
  t = 3s + j, depends only on C. Since 1 + beta + beta^2 = 0,
  beta^2 x = p - ((x + beta x) mod p): **one** Fp multiplication (beta*x).

(Any fixed bijective representation can be used for the comparisons, e.g.
blst's Montgomery representatives xR mod p: multiplication by beta and
beta^2 x = -(1+beta)x commute with it, so the rule stays a function of C.)

## 3. The walk

Parameters: R = 2048 partitions; distinguished probability theta = 2^-64;
M = 2^20 concurrent walks, processed in batches of 2^12 per Montgomery
inversion; trail cap L_max = 20 * 2^64 operations; checkpoint period 64.

**Table (preprocessing).** Draw c_h, d_h uniform in [0, r) from fresh coins
(rejection sampling) and set M_h = [c_h]P0 + [d_h]T, h = 0..R-1. Store
Tab[h][t] = [u_t^{-1}]M_h = (-1)^s phi^{-j}(M_h) (one multiplication of x(M_h)
by beta^{-j} and a possible negation; no scalar multiplications) with
coefficients (u_t^{-1} c_h, u_t^{-1} d_h) mod r. Size 2048 * 6 * 160 bytes
= 1.97 MB.

**Class map.** F(C) = class(c(C) + M_h), h = K(C) mod R (low 11 bits of K).
A walk keeps an arbitrary member P of its class. Since
c(C) + M_h = [u_t](P + [u_t^{-1}]M_h), the class of P + Tab[h][t] is F(C):
**one ordinary step is exactly one batched affine G1 addition**
P <- P + Tab[h][t] (the reference operation) plus the key work of section 5.
The walk never applies phi or negation to its state.

**Distinguished classes.** C is distinguished iff limb 1 of K (bits 64..127)
is zero. **History fingerprints**: K_{-1}, K_{-2}, K_cp below store limbs 0, 2, 3
of K (192 bits); two distinct classes share a fingerprint with probability
2^-192 per comparison.

**Trail.** A trail starts at S = [a]P0 + [b]T with a, b uniform in [0, r) from
fresh algorithm coins (rejection sampling of 255-bit strings; redraw if S = O).
(a, b) are stored in full (no seed expansion), so distinct trails have
independent uniform coefficient pairs. `length` counts every operation of the
trail (ordinary steps, doublings, cycle-traversal steps). At every class arrival:

```
K, t      <- key and twist index of the current class
if K is distinguished: report (K, a, b, length) and end the trail
if length >= L_max:    abandon the trail (no report)
if pending_double:     P <- [2]P ; clear history ; continue
if fp(K) == K_{-2}:    2-cycle {current, previous} detected:
                         if K < K_{-1}: P <- [2]P ; clear history ; continue
                         else pending_double <- true      (one ordinary step reaches the smaller class)
elif fp(K) == K_cp:    cycle of length <= 64 detected: traverse it once more (ordinary steps),
                       keep the member of its minimum-key class, P <- [2]P of that member ;
                       clear history ; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t] (h = K mod R) ; K_{-2} <- K_{-1} ; K_{-1} <- fp(K)
```

Each P <- ... is one operation (length += 1). "Clear history" resets K_{-1},
K_{-2}, K_cp. Comparisons "K < K_{-1}" use the full keys (the previous key is
kept in the walk state). Doubling commutes with every unit, so P -> [2]P is a
class map; each escape doubles the minimum-key class of a cycle of F, which is a
function of the cycle alone, not of the trail's history or entry point. Hence two
trails that reach the same class follow the same class sequence from then on
(identical escapes) and reach the same distinguished class.

**Exceptional additions.** Before each addition the denominator x(Tab) - x(P)
is tested for zero (needed anyway to protect the batch inversion). If it is zero
(P = +-Tab[h][t]) the trail is abandoned. The walk point is then one of 12R
fixed points: probability <= 12R/r per operation, < 2^-114 over the whole run.

**Loops not caught by the detectors.** A trail can only revisit a class through
(a) a fruitless cycle of F of length <= 64 (caught above), or (b) a genuine class
repeat within the trail (a rho loop, possibly passing through escapes, or a
fruitless cycle longer than 64). Case (b) leaves the trail looping until L_max
and it is abandoned. Under H1 its probability is at most L_max^2/(2N) = 2^-116.6
per trail (fruitless cycles longer than 64 have rate < R^-32), i.e. < 2^-53 in
total; its waste (<= L_max per occurrence) is included in the cap terms below.
(At toy sizes, L_max^2/N is not small; the toy's few abandoned trails are such
loops.)

**Phases and hard caps.**

* Phase 1: run the M walks, each starting a fresh trail whenever its previous
  trail ends, until a global **charged-unit counter** reaches
  U1 = 1.311 * S_ops (section 4). The counter adds 1.311 per ordinary or
  traversal step and 2 per doubling (each at least the real cost, section 5).
  If the total number of doublings would exceed 2^-13 * S_ops (an absolute
  cap), stop with failure.
* Phase 2: every walk finishes its current trail (to a distinguished class or
  L_max) and starts no new trail: at most M * L_max operations.
* Reports go to a hash table keyed by 192 bits of K (limbs 0, 2, 3) with
  capacity D_cap = 2^63.25 records. If it would overflow, stop with failure.
* At most 8 key matches are processed; afterwards stop with failure.

**Collision solving (re-walk).** When a new report has the same stored key as a
report of another trail, both trails are recomputed from their stored (a, b)
with the same rules while tracking P = [A]P0 + [B]T (ordinary step:
A += u_t^{-1}c_h, B += u_t^{-1}d_h; doubling: A, B <- 2A, 2B; every value reduced
mod r after every operation, so it fits in four 64-bit limbs). At the
distinguished class the canonical member is c = [u_t]P with canonical
coefficients (alpha, beta) = (u_t A, u_t B) mod r. If the two canonical points
differ (192-bit key coincidence) the match is discarded. Otherwise
alpha_1 + beta_1 x = alpha_2 + beta_2 x (mod r), and if beta_1 != beta_2

**x = (alpha_2 - alpha_1) / (beta_1 - beta_2) mod r.**

(Sign checked on production BLS12-381 data in part A.) The algorithm checks
[x]P0 = T before output.

## 4. Success probability

Let N' = N * R/(R-1) (Brent-Pollard/Teske correction for R-adding walks) and let
T_c be the number of fresh class visits (all trails, trail starts included)
before the first visit of a class already visited by another trail.

**H1 (heuristic), stated only at the budget.** At S = S_fresh = 1.06 S0,
P(T_c <= S_fresh) >= 1 - exp(-(S_fresh/kappa)^2 / (2N')) = 0.4058 with
kappa = 1.05, i.e. the random-mapping (parallel collision search) law with an
allowance of 5% slowdown in collision time at that quantile. The allowance is
not cosmetic: in the toy experiments (section 7) a control walk without the
automorphism (class size 1, same harness) matches the random-mapping mean,
while class walks (sizes 2 and 6) are measurably slower (means +0.5% to +2.1%).
At the budget quantile itself, the pooled class-6 toy value is
P[T_c <= 1.06 S0] = 0.4255 +- 0.0032 (23,800 trials) and the worst
configuration gives 0.413 +- 0.008, against the 0.4058 required by H1.

Budgets:
S0 = sqrt(2 N' ln(1/0.6)) = 2^126.15187 (pure model: probability 0.40),
S_fresh = 1.06 * S0 = 2^126.23593 (H1 with kappa = 1.05: probability >= 0.4058;
pure model: 0.4367),
S_ops = S_fresh / 0.998 = 2^126.23882, U1 = 1.311 * S_ops = 2^126.62949.

Operations versus fresh visits. An operation fails to visit a fresh class only
inside a fruitless cycle: a detected 2-cycle wastes at most 2 visits (the
repeat arrival and the optional step before doubling), a checkpoint cycle at
most 128 + 2*64 = 256. 2-cycles need h' = h and the canonical twist of the next
class to be the negation: rate 1/(6R) = 2^-13.58 per step (toy measurements
agree). A fruitless cycle of length L >= 3 needs the twisted multiplier sum
around it to vanish; for independent random M_h every partition index in it
must occur at least twice, so its rate is at most 4R^-2 + 26R^-3 + ... < 2^-19.9
per step (other closures are O(1/N) coincidences under H1). The expected waste
fraction is below 2/(6R) + 256 * 2^-19.9 < 0.00045, a sum of ~2^112 nearly
independent small contributions, so it exceeds 0.0019 only with negligible
probability. Hence phase 1 makes at least 0.9999 * 0.9981 * S_ops >= S_fresh
fresh visits. The same concentration keeps the total number of doublings below
the absolute cap 2^-13 * S_ops (expected rate < 1/(6R) + 2^-19.9 < 2^-13.5 per
operation); with at most that many doublings, which the counter charges 2 each,
phase 1 performs at least S_ops (1 - 0.689 * 2^-13 / 1.311) >= 0.9999 S_ops
operations.

Detection. Suppose T_c <= S_fresh: trail B enters a class already visited by
trail A != B. From then on B follows A's class sequence; both reach the same
distinguished class (in phase 1 or phase 2), the report matches, and the
re-walk gives x unless: (i) A or B hits L_max first: each trail's length to its
first distinguished class is geometric with mean 2^64; accounting for the
size bias of colliding trails (Gamma(2) tail), probability < 2 * 21 e^-20 < 1e-7;
(ii) the first collision is inside one trail: < 2^-53 (section 3); (iii) report
table overflow: expected reports theta * (S_ops + M L_max) = 2^62.24, and
exceeding D_cap = 2^63.25 has probability < exp(-2^60) (Chernoff);
(iv) beta_1 = beta_2: B's starting b is uniform given its start point (the r
pairs (a, b) with a + bx fixed are equally likely) and
beta_B = (unit) * 2^e * b + (terms fixed by the point sequence), so
beta_1 - beta_2 is uniform: probability 1/r; (v) exceptional additions:
< 2^-114; (vi) more than 8 matches before success, or the doubling cap: < 2^-60;
(vii) a spurious history-fingerprint equality on A or B between the merge and
the report: < 4 L_max 2^-192 < 2^-120.

Success >= 0.4058 - 1e-7 > 0.40 under H1.

## 5. Cost (units: one batched incomplete affine G1 addition, C = 2,000)

**Per ordinary step.** The addition P + Tab[h][t] is the reference operation,
charged 1 unit (6M convention, incl. the per-element 3M of Montgomery's batch
trick). Additional scalar AArch64 instructions (blst_fp_mul = 383 measured):

| Work | Instructions |
| --- | ---: |
| beta * x (one Montgomery Fp multiplication by a constant) | 383 |
| beta^2 x = p - ((x + beta x) mod p): 6-limb add, conditional subtract, negate, zero test | 40 |
| minimum of three 6-limb values (top-limb compare/select; full compare only on top-limb tie) | 12 |
| sign: p - y (6-limb subtract), compare, select limbs of K | 11 |
| h = K mod 2048, distinguished test (limb 1 == 0), twist index t, table address | 7 |
| loads of table point (12 limbs) and walk state, store of new state | 36 |
| 192-bit fingerprint tests against K_{-2} and K_cp, rotating history store, amortised checkpoint store | 16 |
| zero-denominator (exceptional operand) test | 8 |
| length/global counters, cap test, loop control | 6 |
| batch inversion amortised: one Fp inversion (< 2^18 instructions) per 2^12 elements | 64 |
| **total (budgeted 622, i.e. 0.311 units)** | **583** |

The batch-inversion row is charged although the 6M convention already amortises
batching; additions inside the addition formula are part of the reference
operation (A = 0) and are not charged twice. Ordinary-step cost: 1.311 units.
A doubling (batched affine, lambda = 3x^2/(2y): about 7M = 2,681 instructions
plus the 583 above) costs < 1.7 units; the counter charges 2.

**Phase 1:** W1 <= U1 + 2M = 2^126.62949 (hard: the counter stops phase 1;
batched walks can overshoot it by at most one operation each).

**Everything else** (hard caps; worst-case unbatched costs where batching may
shrink: <= 1 + 2^18/2000 + 0.311 < 133 units per operation):

| Item | Bound (units) |
| --- | ---: |
| Phase 2: M * L_max operations, each < 133 units unbatched = 2^20 * 20 * 2^64 * 133 | 2^95.4 |
| Trail starts: <= D_cap + S_ops/L_max + M < 2^63.3 trails, each 2 scalar multiplications (Jacobian double-and-add, inversion, 2 x 255-bit rejection sampling) < 2^12 units | 2^75.3 |
| Report handling: <= D_cap reports x 1,000 instructions (hash, insert/lookup, record copy) | 2^62.3 |
| Hash-table initialisation: 2 D_cap slots x 104 bytes, 8-byte stores | 2^57.1 |
| Re-walks: <= 8 matches x 2 trails x (L_max operations x 133 units + start < 2^12) | 2^79.4 |
| Table: 2 * 2048 scalar multiplications + 5 * 2048 beta multiplications + coefficients | 2^23 |
| Final verification [x]P0 = T, output, identity test | 2^13 |

W <= 2^126.62949 + 2^95.5 = 2^126.62949, submitted as **time_log2 = 126.64**
(margin 0.0105 bits).

**Memory.** Reports: D_cap = 2^63.25 records x 104 bytes (24-byte key, 64 bytes
a, b, 9-byte length, padding) at load factor 1/2: 2^70.95 bytes. Table 1.97 MB,
2^20 walk states x <= 512 bytes (point, a, b, full previous key, fingerprints,
length, flags). Peak < 2^72 bytes (reported, not scored).
**Preprocessing:** the table, < 2^23 units, included. **Advice:** none.

## 6. Heuristics, assumptions and limitations

* H1 is the standard rho heuristic (Pollard 1978; Teske, Math. Comp. 2001, and
  the Brent-Pollard R/(R-1) correction; van Oorschot-Wiener, J. Cryptology 1999)
  applied to classes under an automorphism group (Wiener-Zuccherato SAC 1998;
  Gallant-Lambert-Vanstone Math. Comp. 2000; Duursma-Gaudry-Morain ASIACRYPT
  1999), with an explicit 5% slowdown allowance motivated by our own toy
  measurements of class walks; fruitless cycles and doubling escapes follow
  Bos-Kleinjung-Lenstra (ANTS 2010) and Bernstein-Lange-Schwabe (PKC 2011). The
  quantitative facts used (orbit size, key invariance, cycle rates, detection,
  degeneracy) are derived above, not imported. H1 is not proved.
* Partition index and distinguished test use bits of K = |y|, which lies below
  p/2; these bits are assumed to be pseudo-random functions of the class.
* Toy experiments are evidence at their own scale only (class counts N = 2^21.6
  and 2^25.5); production behaviour is extrapolated under H1. In the toy, trails
  are a noticeable fraction of sqrt(N) (trail length about 1-2% of sqrt(N)),
  so excluding intra-trail repeats from T_c makes the toy slightly pessimistic.
* Generic square-root attack; index calculus in GT via the pairing (TNFS family)
  is not used. An upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (production BLS12-381 G1; all checks print true): P0 on the curve and
[r]P0 = O; beta, lambda primitive cube roots of unity with phi(P0) = [lambda]P0;
(0, +-2) are 3-torsion; for random P the six [u_t]P equal (beta^j x, +-y) and
are distinct; |y| is constant on the orbit; all six members give the same
canonical point with [u_t]P = c; the twisted step P + Tab[h][t] lands in the
class of c + M_h for every member; doubling is a class map; Tab entries equal
[u_t^{-1}]M_h; coefficient tracking along production walk steps reproduces the
walk point; two production trails forced to merge (trail B started at a random
member of a class on trail A, with independent random coefficients) reach the
same canonical point, and x = (alpha_2 - alpha_1)/(beta_1 - beta_2) returns the
planted x with [x]P0 = T.

Part B (toy): y^2 = x^3 + 4 over F_q, q = 1 mod 3, prime-order subgroup r'. The
complete per-trail algorithm of section 3 runs unchanged (canonical rule,
twisted table, 2-cycle escape with pending doubling, checkpoint escape, trail
cap, exceptional-addition abandonment, distinguished classes, restart after
every trail, report table, re-walk solving with coefficients reduced mod r',
verification). Differences from production, none affecting the class sequence:
trails run one after another rather than as 2^20 interleaved walks (this
changes the order of visits, hence individual T_c values, but not their
distribution under H1); the harness evaluates [a]G + [b]T as
[a + b x]G with a fixed-base table purely for speed (same point); full keys are
compared instead of 192-bit truncations (no truncation coincidences exist at
toy size); the distinguished test uses the bits just above the partition bits
(production: limb 1 of K); there is no charged-unit counter, doubling cap,
phase 2, report-table capacity or 8-match cap, because each trial simply runs
until it solves, and the statistic reported is T_c, which is what the budget
argument of section 4 uses. T_c counts every class arrival (trail starts, ordinary steps,
doublings, traversal steps, repeated arrivals in fruitless cycles), so it is
conservative relative to fresh visits.

Offline runs of this file (pure stdlib; every trial ran to completion; every
recovered logarithm was verified; seeds 501-524 for the attack rows, 701-712 for
the control rows; input e.g. `{"curve":"q26","trials":1000,"R":2048,"dp_bits":5,"seed":509,"seconds":5000,"dump":true}`).
S0 and the model use N' = N R/(R-1) with N = (r'-1)/classes; T_c counts all class arrivals.

| classes | q, r' | R | theta | runs | trials | mean T_c/sqrt(N) +- s.e. | model | diff | P[T_c<=S0] | P[T_c<=1.06 S0] | 2-cycle esc/step (expected) | checkpoint esc/step | steps | all solved | abandoned | exceptional |
| --- | --- | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: |
| 6 | 58493047, 19496431 | 256 | 2^-5 | 8 | 8000 | 1.2704 +- 0.0074 | 1.2558 | +1.2% | 0.397 | 0.431 | 6.32e-04 (6.51e-04) | 8.8e-06 | 17969893 | True | 27 | 2 |
| 6 | 58493047, 19496431 | 2048 | 2^-5 | 7 | 7000 | 1.2699 +- 0.0079 | 1.2536 | +1.3% | 0.390 | 0.422 | 8.03e-05 (8.14e-05) | 7.9e-06 | 15730412 | True | 21 | 3 |
| 6 | 861101287, 287046709 | 256 | 2^-6 | 6 | 2400 | 1.2706 +- 0.0131 | 1.2558 | +1.2% | 0.387 | 0.421 | 6.35e-04 (6.51e-04) | 2.2e-06 | 20901918 | True | 8 | 0 |
| 6 | 861101287, 287046709 | 2048 | 2^-6 | 6 | 2400 | 1.2593 +- 0.0135 | 1.2536 | +0.5% | 0.403 | 0.443 | 8.01e-05 (8.14e-05) | 9.2e-07 | 20725564 | True | 9 | 0 |
| 1 | 58493047, 19496431 | 256 | 2^-5 | 4 | 4000 | 1.2563 +- 0.0102 | 1.2558 | +0.0% | 0.390 | 0.428 | 1.85e-07 (0.00e+00) | 1.3e-06 | 21614461 | True | 5 | 2 |
| 2 | 58493047, 19496431 | 256 | 2^-5 | 4 | 4000 | 1.2732 +- 0.0103 | 1.2558 | +1.4% | 0.391 | 0.427 | 1.90e-03 (1.95e-03) | 5.2e-06 | 15493162 | True | 14 | 2 |
| 6 | 58493047, 19496431 | 256 | 2^-5 | 4 | 4000 | 1.2824 +- 0.0103 | 1.2558 | +2.1% | 0.373 | 0.413 | 6.23e-04 (6.51e-04) | 8.9e-06 | 9067214 | True | 12 | 0 |

Reading: the pure model gives P[T_c <= S0] = 0.40 and P[T_c <= 1.06 S0] = 0.437;
H1 with kappa = 1.05 gives 0.406 at 1.06 S0. Every class-6 configuration has
P[T_c <= 1.06 S0] between 0.413 and 0.443 (binomial s.e. 0.005-0.010). The
control rows separate the effect: the plain walk (classes = 1) matches the
random-mapping mean exactly (+0.0%), while negation classes (+1.4%) and the
six-classes of the attack (+0.5% to +2.1%) are slightly slower; this is why H1
carries the explicit kappa = 1.05 allowance at the budget quantile. For the negation-only control the
2-cycle rate is 1/(2R) (1.95e-3 expected), as the analogous argument predicts.
2-cycle escape rates match 1/(6R) at both R.
Abandoned trails are intra-trail loops longer than 64 (section 3, case (b)),
frequent only because L_max^2/N is not small at toy size; "exceptional" counts
the toy's operand coincidences P = +-Tab (rate at most 12R/r' per step, in
practice about 1/(3N)), which abandon the trail as specified.

A single 14-second run (the declared initial input) completes a few hundred
trials, so its fractions have standard error about 0.025; it is a check that
the code runs the algorithm correctly, not a precise estimate.

References: Pollard (1978); van Oorschot and Wiener (1999); Teske (2001);
Wiener and Zuccherato (1998); Gallant, Lambert and Vanstone (2000); Duursma,
Gaudry and Morain (1999); Bos, Kleinjung and Lenstra (2010); Bernstein, Lange and
Schwabe (2011).

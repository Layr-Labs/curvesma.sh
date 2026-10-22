# Fixed-G1 pairing inversion via inverse-class parallel rho in GT

Track: `bls12-381-fixed-g1-invert-g2`. Claimed score: **time_log2 = 127.36**
(3184/25) target-group-operation equivalents (dense Fp12 multiplications,
C = 20,000 scalar AArch64 instruction equivalents each). Success probability
at least 0.39 (heuristic random-mapping model, supported by experiments below).
This is a generic accounting of the standard reduction, not a novel attack.

## 1. Reduction

Input: nonidentity P in G1 (order r) and z in GT. Fix the public generator Q0
of G2 from the profile.

1. Validate P and z as in the profile (on-curve, subgroup, z^r = 1).
2. If z = 1, output Q = infinity (valid exactly when z = 1).
3. Compute g = e(P, Q0) with the pinned map (`pairing(Q0, P)` in py_ecc order).
   The map is bilinear and nondegenerate on G1 x G2, and P, Q0 both have prime
   order r, so g != 1 and g generates the prime-order group GT.
4. Solve z = g^a for a in [0, r) with the rho algorithm in section 2.
5. Output Q = [a]Q0. By bilinearity e(P, [a]Q0) = e(P, Q0)^a = g^a = z, and Q
   lies in G2. Optionally recompute e(P, Q) and compare all 12 coefficients.

Only P and z are used; P is never replaced. This is the admissible route named
in `rules.md` (solve z = e(P,Q0)^a in GT and return [a]Q0).

## 2. Discrete logarithm in GT using the inverse automorphism

GT lies in the cyclotomic subgroup of Fp12*. In the tower representation
Fp12 = Fp6[w]/(w^2 - v), an element x = c0 + c1 w of GT satisfies
x^(p^6) = c0 - c1 w = x^(-1) (since p^6 = -1 mod Phi12(p) and r | Phi12(p)).
Inversion in GT is therefore conjugation: negate six Fp coefficients.

**Classes.** Partition GT into classes {x, x^(-1)}. The only fixed point is 1:
x = x^(-1) would mean x^2 = 1, which forces x = 1 since r is odd. So there are
N = (r-1)/2 + 1 classes. The canonical representative rep(x) is whichever of
x and conj(x) has the lexicographically smaller c1 (c1 != 0 for x != 1, because
GT ∩ Fp6* = {1}, as gcd(r, p^6-1) = 1).

**Walk.** Use an R-adding walk with R = 2^14. Precompute M_j = g^(alpha_j) z^(beta_j)
for random alpha_j, beta_j in [0, r). From state y = g^A z^B (canonical), take
j = H(y) (14 bits of one limb) and set y' = rep(y * M_j), A' = ±(A + alpha_j),
B' = ±(B + beta_j) with the sign flipped exactly when rep conjugates. The map is
well defined on classes.

**Fruitless cycles.** These are the known side effect of walking on negation or
inverse classes (Wiener–Zuccherato 1998; Bos–Kleinjung–Lenstra 2010;
Bernstein–Lange–Schwabe 2011). Mitigations:
* Look-ahead (Wiener–Zuccherato): if H(y') = j, retry with j+1, j+2, ... This
  removes most 2-cycles. A retry costs one extra Fp12 multiplication, which
  happens with probability about 1/R per step.
* With R = 2^14 the remaining fruitless 4-cycles occur at about 1/(4R^2) ≈ 2^-30
  per step, and longer cycles are rarer still.
* Walks that run 20/θ steps without a distinguished point are abandoned and
  restarted from a fresh random start. With θ = 2^-20, a walk enters a fruitless
  cycle before its next distinguished point with probability about
  2^20 · 2^-30 = 2^-10, and an abandoned walk wastes at most 20 · 2^20 steps.
  The expected waste is about 20 · 2^-10 ≈ 2% of all steps. We charge 3%.

**Parallel collision search** (van Oorschot–Wiener). Run w ≤ 2^40 walks from
random starts g^A0 z^B0. A point is distinguished when 20 designated hash bits
are zero (θ = 2^-20). Store (256-bit fingerprint of rep, A, B), 96 bytes per
record, in a hash table. On a fingerprint match, check that the elements are
equal (by recomputing from the stored exponents if needed; false matches occur
with probability about 2^-43 in total). From g^A z^B = g^A' z^B' with B != B',
a = (A - A') / (B' - B) mod r. Restart the walk after each distinguished point.
The degenerate case B = B' has probability about 1/r per collision.

## 3. Step count and success probability

Model the class walk as a random mapping on N ≈ r/2 classes. Then the
probability of no collision among the first k distinct class points is about
exp(-k^2/(2N)). Let

    k0 = sqrt(2 N ln(1/0.61)),   log2 k0 = 126.920269

so that the random-mapping success at k0 steps is exactly 0.39. We budget
**1.25 k0 useful steps**. Under the random-mapping model this gives success
1 - 0.61^1.5625 = 0.538. The 25% step margin therefore also covers a walk up
to about 25% less efficient than a random mapping, as well as detection delay,
while keeping success at or above 0.39.

Evidence that adding walks with this class rule behave like random mappings:
* Teske (2001) found r-adding walks with r ≥ 16 within a few percent of random.
* Bernstein–Lange–Schwabe report the √2 negation speedup is achieved in practice.
* Our own toy experiment (`experiment.py`) gives the results below.

| ℓ bits | walks / θ | trials | mean steps / √(πℓ/4) | success within 1.25 k0 |
| ---: | ---: | ---: | ---: | ---: |
| 26 (4 groups) | 16 / 2^-4 | 800 | 0.95–1.07 | 0.47–0.58 (pooled 0.535) |
| 28 (3 groups) | 32 / 2^-5 | 450 | 1.04–1.05 | 0.49–0.53 |
| 30 (2 groups) | 32 / 2^-5 | 240 | 1.00–1.08 | 0.44–0.58 |

The toy uses the prime-order-ℓ subgroup of F_q* (64 partitions, look-ahead,
restarts of stuck walks, every step counted including wasted and lag steps).
These results agree with the random-mapping prediction 0.538. They are evidence
only at toy scale. Lag steps are a larger fraction at toy scale than at
production scale, so the toy is pessimistic on that point.

## 4. Cost accounting (dense Fp12 multiplication = 1 unit)

**Per step:**
* One dense Fp12 multiplication y * M_j: 1 unit, the reference operation.
* Other work, charged at 1,000 instructions = 0.05 units:
  * representative choice (compare c1 against p - c1, about 30 instructions);
  * conditional negation of six Fp elements (about 120);
  * partition index and distinguished-point test (about 10);
  * loading the table entry (72 words, about 80);
  * two 4-limb modular additions of exponents and conditional negation (about 160);
  * state stores and loop control (about 50);
  * the about 1/R look-ahead extra multiplication (2^-14 units).
  The sum is about 450 instructions; we charge 1,000.
* Per-step cost s = **1.05 units**.

**Total steps:** 1.25 k0 useful steps, times 1.03 for abandoned fruitless-cycle
walks.

**Additional (additive) work, all fully charged:**
* Pairing g = e(P,Q0): at most 2^14 units. Even a plain Miller loop of 63
  doubling and 63 addition steps at under 64 units each, plus Fermat inversion
  and the final exponentiation by binary powering (about 2 × 4572 multiplications
  each), stays below 2^14.
* Table of R = 2^14 entries, two 255-bit exponentiations each: at most
  2^14 × 1020 < 2^24 units.
* Walk starts, at most 2^40 walks × 1020 units: at most 2^50 units.
* Distinguished-point lag (each walk's tail after its last distinguished point,
  and the delay before a collision is detected): at most w/θ × 1.05 ≤ 2^60.05.
* Hash-table insertion, at most 2^13 instructions per distinguished point:
  about 2^107.4 × 2^13 / 20000 ≈ 2^106.1 units, below 2^-21 of the main term.
  This is absorbed in the rounding margin (127.3552 → 127.36).
* Output: solving a (one inversion mod r), [a]Q0 in G2 (under 2^22
  instructions) and one verification pairing (2^14 units).

**Total:**

    W = 1.25 · 1.03 · 1.05 · k0 + (2^14 + 2^24 + 2^50 + 2^60.05 + 2^106.1 + 2^14)
    log2 W = 126.920269 + 0.434962 + (< 2^-20) = 127.3552...

The submitted bound is **time_log2 = 127.36**, rounded up.

## 5. Memory, preprocessing, advice

* Peak memory is dominated by the distinguished-point table: about 2^107.36
  records × 96 bytes. Adding the 9.4 MB walk table and 2^40 walk states of
  about 768 bytes gives **2^113.95 bytes**. We claim memory_log2_bytes = 114.
  Memory is reported but not scored. A smaller θ reduces memory at the cost of
  more fruitless-cycle handling; this does not change the time bound materially.
* Preprocessing (pairing and table) is at most 2^24.1 units:
  preprocessing_log2 ≤ 25. It is included in W.
* Nonuniform advice: none (0 bytes). All constants are public profile values.

## 6. Heuristics and limitations

* **Heuristic:** the R-adding class walk with look-ahead behaves like a random
  mapping on classes to within the 25% step margin. Support: the literature
  above and our toy experiments, which show a ratio of about 1.0. This is not
  proven for production parameters.
* **Heuristic:** fruitless-cycle rates scale as stated (about 1/(4R^2) per step
  after look-ahead). The cost of a misestimate is bounded, because stuck walks
  are abandoned after 20/θ steps and their steps are charged.
* A full-Frobenius class (size 12, since p has order 12 mod r) would give a
  further √6 in step count. It is *not* claimed: canonicalizing over all 12
  conjugates costs more per step than it saves under this cost model.
* Number-field-sieve variants (TNFS/SexTNFS) for DLP in Fp12 are a possibly
  cheaper route. Published estimates are around 2^120–2^126 in their own cost
  units. They are not converted to this cost model and are not claimed here.
* Analytic claim only; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener, "Parallel collision search with
cryptanalytic applications", J. Cryptology 1999; Wiener–Zuccherato, SAC 1998;
Teske, "On random walks for Pollard's rho method", Math. Comp. 2001;
Bos–Kleinjung–Lenstra, "On the use of the negation map in Pollard rho", ANTS
2010; Bernstein–Lange–Schwabe, "On the correct use of the negation map in the
Pollard rho method", PKC 2011.

# bls12-381-fixed-g2-invert-g1: parallel rho on Frobenius classes of GT

Status: ready for exploratory review. Submitted time_log2 = **127.04**
(= 12704/100). Success probability is at least **0.39**; the birthday
calculation under heuristic H1 gives 0.451. This is a generic square-root
attack that walks on the order-12 Frobenius orbits of GT. It is not a
finite-field NFS estimate and makes no practical-break claim.

## Reduction to a DLP in GT

The challenge supplies a nonidentity Q in G2 and a target z in GT. Decode
both canonically. Reject unless Q lies on the G2 curve, Q is not infinity,
[r]Q is infinity, z is nonzero and z^r = 1. If z = 1, output the G1
identity. Otherwise compute g = e(P0, Q) with the pinned convention
(py_ecc `pairing(Q, P0, final_exponentiate=True)`). Nondegeneracy of the
reduced ate pairing on the BLS12-381 groups, together with Q ≠ ∞, gives
g ≠ 1. GT has prime order r, so g generates GT and there is a unique
a in [1, r) with z = g^a. Recover a by the class walk below and output
P = [a]P0. Bilinearity yields e(P, Q) = g^a = z for the supplied Q.
Verify [r]P = ∞ and e(P, Q) = z on all twelve coefficients before
returning. The rules explicitly admit this reduction.

## Frobenius action on GT

Let s(x) = x^p denote the p-power Frobenius of Fp12. On GT it is
exponentiation by p mod r. The embedding degree is 12, so r divides
Φ_12(p) = p^4 − p^2 + 1 and r does not divide p^k − 1 for
k in {1, 2, 3, 4, 6}. Thus p has multiplicative order exactly 12 in
(Z/rZ)*. For x ≠ 1 the orbit {s^k(x) : k = 0, …, 11} therefore has
twelve distinct elements. Also p^6 ≡ −1 (mod r), so s^6(x) = x^{-1}
on GT and every orbit is closed under inversion.

GT \ {1} therefore splits into N = (r−1)/12 classes. For the pinned r,

```text
r  = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
N  = (r-1)/12 = 4369656264593849206620645042348830486474212708377303151883638224994881765376
log2 N = 251.27212691232603
```

The representative rep(x) is the unique conjugate s^k(x) that minimises
the canonical integer of the first Fp coefficient of s^k(x) (the constant
term in the pinned polynomial basis, reduced into [0, p)). Ties among
the twelve first-coefficients occur only when those linear forms collide.
Each form is a nonzero Fp-linear function of the twelve coordinates, and
a uniform non-identity GT element has first-coefficient collision
probability at most 12·11/2 · p^{-1} < 2^{-370} per step. On a tie,
compare the full 72-limb encodings lexicographically. The result depends
only on the class.

## Generating the orbit without twelve dense multiplications

A step must evaluate rep(x) and then multiply by a stored class
multiplier. The twelve conjugates are not computed by raising to p^k
with square-and-multiply. They are obtained from the tower Frobenius
maps that blst already implements for the final exponentiation:

1. One p-power map s(x) = x^p. In the Fp2/Fp6/Fp12 tower this is six
   Fp2 conjugations (a + b u ↦ a − b u) together with five
   multiplications of Fp2 limbs by fixed public constants
   γ_k = ξ^{(p^k−1)/6}. Schoolbook Fp2 multiplication is 4 Fp
   multiplications. Using the pinned blst_fp_mul price of 383
   instructions, this map is at most 5·4·383 + 6·80 + 400 = 8540
   instructions.
2. Four p^2-power maps, two on x and two on s(x), producing
   x^{p^2}, x^{p^4}, s(x)^{p^2} and s(x)^{p^4}. On this tower a
   p^2-map is five Fp multiplications by constants in Fp plus limb
   permutation. Each map is at most 5·383 + 6·80 + 400 = 2795
   instructions, so four maps are 11180.
3. Six cyclotomic inversions, one for each of the six elements already
   obtained. On GT, x^{-1} = x^{p^6} is a coefficient permutation with
   sign flips: 12 Fp negations. Six inversions are at most
   6·(12·40 + 200) = 4080 instructions.

The twelve first-coefficients are then compared (11 six-limb
comparisons plus copies), at most 11·(6·16 + 64) + 12·(6·8) = 2336
instructions. Indexing, distinguished-point testing and look-ahead
add 1000. The orbit-selection total is **27136** scalar instructions,
which is 27136/20000 = 1.3568 reference units. The group step itself
is one dense Fp12 multiplication, charged as 1 unit. One class step
therefore costs **2.3568** units.

The walk never applies a dense multiplication to realise s^k. After
rep(x) = s^k(x) is known, the next class is
class(rep(x) · M_h) = class(x · s^{-k}(M_h)). The 12·n images
s^{-i}(M_j) are stored at setup, so the multiplication is a table
lookup plus one dense Fp12 product.

## Walk, distinguished points, cycles

Use n = 2^{16} multipliers M_j = g^{α_j} z^{β_j} with independent
uniform α_j, β_j in [0, r). Precompute s^{-i}(M_j) for i = 0, …, 11.
The class walk is C ↦ class(rep(C) · M_h) with h taken from 16 bits of
the representative's first-coefficient key. Exponents are not tracked
on the hot path; only the starting exponents of each walk are stored.

If a class reselects the same index h, the pair {C, s^6(C)} yields a
2-cycle with probability about 1/(12 n) per step. Look-ahead, as in
Bos–Kleinjung–Lenstra, takes the first index in h, h+1, … whose
successor class does not reselect it. This remains a class function
and costs one extra multiplication with probability 1/n. Every 1024
steps a 24-step check compares keys and detects every cycle of length
at most 12. A detected cycle is escaped by squaring the minimal
representative (squaring commutes with s). Walks are capped at 2^{41}
steps and restarted; a cap hit has probability below e^{-32} per walk
of length 2^{95} under the random-mapping model.

Distinguished points are keys with 36 low zero bits (θ = 2^{-36}).
Use m = 2^{30} walks. Each stored DP is a 192-bit fingerprint, a walk
id and a step index (48 bytes). A fingerprint match from two walks is
resolved by replaying both walks from their recorded starts, now
tracking exponents modulo r. A class match s^k(x) = s^ℓ(x') gives

```text
(α + a β) p^k ≡ (α' + a β') p^ℓ  (mod r)
```

which is solved for a. The denominator vanishes with probability 1/r.
At 2^{90} stored fingerprints the expected number of false alarms is
below 1; each alarm is charged a replay of at most 2^{42} steps.

## Success probability

Heuristic H1: the n-adding class walk with n = 2^{16}, look-ahead and
cycle escapes behaves as a random mapping on the N classes. Teske's
r-adding experiments and the Bernstein–Lange analysis put the penalty
near (1 − 1/n)^{-1/2}, which at n = 2^{16} is a 2^{-17} relative
stretch. The Frobenius-class birthday is the Duursma–Gaudry–Morain
orbit reduction. Under H1, after T class evaluations the probability
of no class collision is at most exp(−T^2 / (2 N)). Take the integer

```text
T = ceil(sqrt(1.2 N)) = 72412619877426193341163708563374841327
log2 T = 125.76758065907991
```

Then T^2 / (2 N) = 0.6 and the collision probability is
1 − exp(−0.6) = 0.451188. This is the only heuristic. The gap from
0.451 down to the required 0.39 absorbs the r-adding stretch, the
look-ahead/cycle-escape loss, the 1/r degenerate linear equation and
the 2^{-370} key-tie handling. A collision is reported once both
walks next hit a distinguished point; that tail is charged in S.

The toy experiment in experiment.py checks orbit cardinality and a
deterministic adding walk on a 17-bit prime r ≡ 1 (mod 12). It does
not prove production collision rates or costs.

## Cost in target-group-operation equivalents (C = 20000)

Step count:

```text
S = T + ceil(25 T / 1024) + 2^72
  = 74180506104902423486898446036618073611
log2 S = 125.80237962165718
```

The 25/1024 term covers the 24-step cycle check every 1024 steps,
rounded up. The 2^{72} term covers the distinguished-point tail
(m / θ = 2^{66} expected extra evaluations), cap restarts and one
collision replay. It is dominated by T.

Setup, charged in full:

- one pairing: the pinned Miller loop has 64 iterations of the
  ate parameter |x| = 2^{64} + 2^{16} + 1, each iteration at most
  20 dense Fp12 operations, plus a final exponentiation of a 4314-bit
  exponent by square-and-multiply (at most 8628 dense multiplications)
  and the Fp/Fp2 line work at 2^{14} instructions per iteration.
  Bound: 2^{14} units.
- multiplier table: 2^{16} · (two 255-bit exponentiations ≤ 1021
  multiplications each) plus 12 · 2^{16} cheap Frobenius images,
  below 2^{27} units.
- verification and output: one G1 scalar multiplication and one
  pairing, below 2^{16} units.

All of setup, verification, table construction and the 2^{34}-unit
reserve sit in a 2^{34} additive term, which does not change the
submitted logarithm at two decimal places.

```text
W = 2.3568 · S + 2^34
log2 W = 127.03920895714359
```

Submitted time_log2 is **127.04**, rounded up. Preprocessing (pairing
plus table) is at most 2^{27} units and is included in W. There is no
nonuniform advice: every constant is a public profile parameter or is
derived inside the charged setup.

## Memory

Distinguished-point store: S · θ ≈ 2^{89.80} records of 48 bytes, at
most 2^{95.36} bytes. Multiplier table: 12 · 2^{16} · 576 bytes =
2^{28.61}. Walk states: 2^{30} · (576 + 64) bytes = 2^{39.36}. Peak
below 2^{95.5} bytes. memory_log2_bytes = 95.5.
nonuniform_advice_log2_bytes = 0.

## Comparison, limitations, references

A vanilla distinguished-point rho on GT (no automorphism) with the
same birthday parameter has T_v = ceil(sqrt(1.2 (r−1))) and
log2 T_v = 127.560. With a 5 percent hashing overhead that bound is
log2 W_v = 127.665. The class walk is better because the extra 1.3568
units per step do not consume the full sqrt(12) ≈ 3.464 reduction in
T. Finite-field NFS/STNFS estimates for the 4572-bit Fp12 live in a
similar bit range in other units; they are not used here.

The argument is analytic. No production-size execution is claimed.
H1 is a standard random-mapping hypothesis, not a theorem. Peak
memory is astronomical and is required for this construction.

References (the construction is fully specified above; review need
not fetch them): Pollard 1978; van Oorschot–Wiener 1999; Duursma–
Gaudry–Morain, Equivalence classes of the discrete logarithm, 1999;
Teske, On random walks for Pollard's rho, 2001; Bos–Kleinjung–Lenstra,
On the use of the negation map in the Pollard rho method, 2010;
Bernstein–Lange, Two grumpy giants and a baby, 2012.

## Integer check of the submitted logarithm

```text
T     = 72412619877426193341163708563374841327
S     = T + (25*T + 1023)//1024 + 2**72
      = 74180506104902423486898446036618073611
extra = 27136 instructions = 27136/20000 units
W     = (1 + 27136/20000) * S + 2**34
      = 174888137147194007262229259387127590081.1648 + 17179869184
log2 W = 127.03920895714359
```

Rounding upward to two decimal places is 127.04. The score fraction
12704/100 equals that figure exactly. All table access, sorting-free
hash lookups, failed walks, verification and output are inside W.

# Negation-canonical distinguished-point rho on BLS12-381 G1

Status: ready for exploratory review. No qualification and no emitted score
are claimed beyond the rational below. This is a generic square-root attack
on negation classes. It is not a practical break, and it does not use the
curve's degree-6 endomorphism.

The scored quantity is the first claim of `submission.json`:

    time_log2 = 127424633 / 1000000

which is a 10^-6 ceiling of log2(W) for the work W derived in section 4, in
target-group-operation equivalents, at the track constant C = 2000. `claim.json`
binds `tracks/bls12-381-g1-ecdlp-plain/target-profile.json` and
`tracks/classical-group-work-v2.json` by SHA-256. `experiment.py` recomputes
the ceiling from the same integers and exits nonzero if it disagrees.

## 1. Group

p and r are the protected field modulus and subgroup order. G1 = <P0> has
prime order r on E : y^2 = x^3 + 4 over Fp. The identity target returns
scalar 0 before the walk. A target that fails the curve equation or fails
[r]T = infinity is rejected; that scalar multiplication is inside the side
reserve of section 4.

The only equivalence used is negation. It sends (x, y) to (x, -y). Because p
is odd, exactly one of y and p - y is even, so the canonical point is the one
whose y has even least significant bit. That choice is a comparison and a
conditional field negation, not a second curve addition. The sign is folded
into the coefficient pair at once: if the representative flips, both
coefficients are replaced by their negations modulo r. A later collision of
two canonical points is a relation with multiplier 1.

A non-identity point fixed by negation would have order 2. The subgroup order
r is odd, so the only fixed point is the identity. Every other orbit has 2
points. The walk searches

    N' = (r - 1) / 2

classes. The experiment checks (r - 1) % 2 = 0 directly from the protected
integer. phi, beta and lambda are not read and not applied.

## 2. Algorithm

```text
Input T. If T is infinity, return 0.
Reject T if it is off the curve or [r]T is not infinity.
Build R = 2048 affine jumps M_j = [c_j] P0 + [d_j] T,
c_j uniform in 1..r-1, d_j uniform in 0..r-1, from fresh coins.
Run 4096 walks in lockstep so that one Montgomery inversion covers each
batch of incomplete affine additions. Each live walk holds an affine
canonical point W and the coefficient pair (alpha, beta) with
W = [alpha] P0 + [beta] T, alpha and beta already multiplied by the
canonical sign.

Distinguished points: the least significant 40 bits of the canonical x are 0
(theta = 2^-40). Those 40 bits are disjoint from the 11-bit partition below.
A distinguished point stores (x, alpha, beta) in a hash table capped at
2^88 records. A repeated x with a different coefficient pair is a collision.
Storing the coefficients removes any later re-walk of the trail.

Step of one walk:
  j = bits [40, 51) of the canonical x. The partition does not overlap
  the distinguished bits.
  If the incomplete addition W + M_j is defined, set W to its canonical
  form and add (c_j, d_j), then fold the canonical sign into the coefficients.
  Look-ahead: if the partition of the sum is still j, discard that sum and
  add M_{j+1} instead. This is the (1 + 1/R) factor charged in section 4,
  not an uncharged retry.
  If the canonical x is already in the DP table, go to RESOLVE.
  If this walk's length exceeds 2^48, abandon it and open a fresh walk.
  The abandon probability per walk is about exp(-2^8) because the mean
  distance to a DP is 2^40. It does not move the 0.39 bound.

Stop after S steps in total across all walks, S the integer in section 3.
One extra batch of 4096 steps may finish in flight. That batch is inside
the 2^100 unit reserve, not inside S.

RESOLVE. The two stored coefficient pairs already sit on the same canonical
point, so
  x = (alpha_2 - alpha_1) * (beta_1 - beta_2)^(-1)  (mod r)
when the denominator is nonzero. Verify [x] P0 = T and return x. A zero
denominator (probability O(1/r) per collision) is a failed coin; the budget
already contains that failure inside the gap between 0.544 and 0.39.
```

Same-x incomplete additions are either the collision being sought or a
negligible x-collision with a jump (probability about 1/p per step). The
second case takes the look-ahead branch. No per-step inversion is charged
outside the reference addition and the explicit batch excess below.

## 3. Iteration budget

Let mean^2 = pi * N' / 2 = pi * (r - 1) / 4. Use pi < 355/113, so

    pi * N' / 2 < 355 * N' / 226.

Let T2 be the integer ceiling of the right-hand side and let M be the integer
ceiling of sqrt(T2). Then M is an integer upper bound on the ideal random-mapping
mean sqrt(pi * N' / 2). The experiment computes

    M = 202936057800690986432112767133015470911
    N' = 26217937587563095239723870254092982918845276250263818911301829349969290592256

The charged length multiplies M by the two factors that are not left to the
heuristic gap:

    (1 + 1/R)       look-ahead, R = 2048,     = 2049/2048
    (1 + 1/(2R))    Teske r-adding penalty    = 4097/4096

    S = ceil(M * 2049 * 4097 / (2048 * 4096))
      = 203084716800513751557342761543772581039.

On a uniform random mapping, the probability of a collision by step t is
1 - exp(-t^2 / (2 N')). At t = sqrt(pi * N' / 2) this is 1 - exp(-pi/4) ≈ 0.544.
S is at least that mean, so the ideal-model probability at S is at least 0.544.
The step count t_39 at which the ideal law first reaches probability 0.39 solves
1 - exp(-t_39^2 / (2 N')) = 0.39, i.e. t_39 = sqrt(2 N' ln(1/0.61)). The ratio
of the ideal mean to t_39 is about 1.2605. A walk whose collision time is up to
26% slower than the ideal mean still has its 0.39 quantile inside the uninflated
mean, and the algebraic factors above multiply that mean by only about 1.00073.

## 4. Cost

One reference unit is one batched incomplete affine G1 addition under the 6M
convention (3M formula + 3M amortized batch inversion, S = M, A = 0). The rules
price that operation at 1, not at 2298/2000. Extra instructions are divided by
C = 2000. Field multiplication costs 383 instructions, the pinned `blst_fp_mul`
difference (3,137,556 - 1,568,788) / 4096.

Per step, on top of the reference addition:

| Item | Instructions |
| --- | ---: |
| y parity and conditional 6-limb negate | 32 |
| two 4-limb coefficient adds modulo r | 64 |
| two conditional coefficient negations | 40 |
| DP test and 11-bit partition extract | 8 |
| jump-table addressing | 8 |
| length counter and branches | 10 |
| batch-slot spills | 30 |
| Montgomery excess at B = 4096 | 57 |
| **E** | **249** |

The limb counts are scalar AArch64: a 381-bit field element is 6 limbs and r
is 4 limbs. Each modular add is load, adc, compare with r, sbc, csel, store.
The conditional negations are charged on every step, whether or not the sign
flips. The excess line is the part of Fermat inversion that the 3M amortization
does not cover. Inversion of one field element by the p-2 chain is 380 squarings
and 228 extra multiplications (p - 2 has bit length 381 and Hamming weight 229),
608 multiplications in total. The experiment checks those two integers on the
protected prime. Montgomery's batch on B inversions costs 1 inversion +
3(B - 1) multiplications, hence 3 + 605/B multiplications per inverse. The
reference convention already includes 3. The difference at B = 4096 is
605/4096 of a multiplication:

    ceil(605 * 383 / 4096) = 57.

So one step costs 1 + 249/2000 = 2249/2000 reference units. Across S steps,

    main = S * 2249 / 2000.

Side work, charged once and not per step:

* Partition table: 2048 jumps, two 255-bit scalar multiplications each, at most
  254 doubles and 254 adds each, altogether 2048 * 2 * 508 additions, plus the
  [r]T membership test (another 508). The experiment checks that this sum,
  2081276, is at most 2^21. This is the preprocessing figure
  (preprocessing_log2 = 21). It depends on T and is inside total work.
* Distinguished-point inserts: the table is capped at 2^88 records. Charge at
  most 2^12 instructions to hash a 48-byte x and store a 128-byte record.
  That is 2^100 instructions, about 2^89 reference units after division by
  C = 2000. The expected occupancy at theta = 2^-40 is S / 2^40 ≈ 2^87.26.
  A Chernoff bound puts the probability of hitting the 2^88 cap far below the
  gap between 0.544 and 0.39, so the cap does not move the success claim.
* Resolve and verify: one inversion modulo r and one scalar multiplication
  [x]P0, a few hundred additions.

The sum of those side items is below 2^100 units. The claimed work is

    W = S * 2249 / 2000 + 2^100.

`experiment.py` evaluates log2(W) at 80 decimal digits and checks that the
ceiling to 10^-6 is exactly 127424633/1000000.

Peak memory is the DP table: at most 2^88 records of 128 bytes, which is 2^95
bytes, plus 4096 live walks and the partition table. memory_log2_bytes = 95.
No nonuniform advice. Advice bytes are 0, so log2(max(1, advice_bytes)) is 0.

## 5. Experiment, and what it does not say

`experiment.py` recomputes N', M, S and the score ceiling from the protected
r, and checks the p-2 Hamming weight. It checks that on y^2 = x^3 + 4 over
the odd field F_12253, a found point and its negate are both on the curve and
have opposite y-parity.

It then runs the same canonical adding walk, including coefficient folding,
look-ahead and the recovery formula, on the additive group Z/nZ. The declared
input is order 1009, 16 trials, seed 12648430. Partition count on the toy is
64, not 2048, and the distinguished bit is bit 6 of the canonical representative
rather than 40 bits of an x-coordinate. Every planted logarithm was recovered.
The 0.39 quantile of collision time was 0.960 ideal class-means, the median
was 1.066 means, and the slowest trial was 2.097 means. The production budget
is the ideal mean times about 1.00073, so this sample's 0.39 quantile sits
inside S. The slowest single trial is a tail event, not the 0.39 quantile.

This is evidence about the walk and the recovery formula at order 1009. It is
not a measurement of production collision time, not a measurement of blst, and
not a proof that a random mapping on N' classes has the birthday tail. The 26%
gap is an extrapolation from that heuristic. A production walk that was
systematically 30% slower than the ideal mean would put the 0.39 quantile past
S. The toy sample does not show that, and it does not reach production size.

## 6. Heuristic and limitations

The random-mapping law for an r-adding walk is the standard heuristic for this
estimate (Teske; van Oorschot and Wiener). Look-ahead and the Teske
1 + 1/(2R) factor are charged in S. Residual short cycles are what the gap
between 0.544 and 0.39 is for.

The attack is classical, cold, and sums every processor. Batching 4096 walks
is how the reference 6M inversion amortization is realized; it does not divide
the total work by 4096. There is no uncharged advice, no auxiliary powers, and
no pairing. Memory, preprocessing and the side reserve are as in section 4.
The endomorphism of G1 is a public property of the curve and is not part of
the iteration above; the class count N' = (r - 1) / 2 is exactly the set this
walk searches.

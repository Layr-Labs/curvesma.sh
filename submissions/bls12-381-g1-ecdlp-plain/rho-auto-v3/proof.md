# Automorphism-aware Pollard rho v3 on BLS12-381 G1

Status: ready for exploratory review. No qualification and no emitted score
are claimed beyond the rational below. This is a generic square-root attack
on automorphism classes. It is not a practical break.

The scored quantity is the first claim of `submission.json`:

    time_log2 = 127016229 / 1000000

which is a 10^-6 ceiling of log2(W) for the work W derived in section 4, in
target-group-operation equivalents, at the track constant C = 2000. `claim.json`
binds `tracks/bls12-381-g1-ecdlp-plain/target-profile.json` and
`tracks/classical-group-work-v2.json` by SHA-256. `experiment.py` recomputes
the ceiling from the same integers and exits nonzero if it disagrees.

## 1. Group

p and r are the protected field modulus and subgroup order. G1 = <P0> has
prime order r on E : y^2 = x^3 + 4 over Fp. The identity target returns
scalar 0 before the walk.

The automorphism group A = <-1, phi> has order 6:

* Negation sends (x, y) to (x, -y).
* phi(x, y) = (beta * x, y), where beta is the cube root of unity in Fp

      beta = 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac

  and beta^2 + beta + 1 = 0 in Fp. With z = -0xd201000000010000,

      lambda = z^2 - 1 mod r = 0xac45a4010001a40200000000ffffffff

  satisfies lambda^2 + lambda + 1 = 0 mod r. The sign of beta is the one with
  phi(P0) = [lambda] P0. Both constants are functions of the public profile.
  The single scalar multiplication that checks the sign is charged in the
  side reserve of section 4. The other cube root matches r - 1 - lambda and is
  rejected by that check.

A non-identity point fixed by a nontrivial element of A would be killed by
[mu - 1] for a nontrivial sixth root mu mod r. mu - 1 is invertible mod r, so
the only fixed point is O. Every other orbit has 6 points. The walk searches

    N' = (r - 1) / 6

classes. r ≡ 1 (mod 6) because r is the BLS prime x^4 - x^2 + 1 with x ≡ 1
(mod 3) and r odd, which the experiment checks directly as (r - 1) % 6 = 0.

Canonical representative of a finite point: compute x, beta*x and
beta*(beta*x) with two field multiplications. Keep the minimum as an integer
in [0, p). Of the two y-coordinates on that x, keep the one whose least
significant bit is even. Exactly one of y and p - y is even, because p is
odd and the subgroup order is odd (so y = 0 does not occur). The multiplier
applied to the discrete log is eps * lambda^k with eps in {±1} and k in
{0, 1, 2}. It is folded into the coefficient pair at once, so a later
collision of two canonical points is a relation with multiplier 1.

## 2. Algorithm

```text
Input T. If T is infinity, return 0.
Reject T if it is off the curve or [r]T is not infinity.
Build R = 2048 affine jumps M_j = [c_j] P0 + [d_j] T,
c_j uniform in 1..r-1, d_j uniform in 0..r-1, from fresh coins.
Run walks in lockstep batches of B = 4096 so that one Montgomery inversion
covers the batch. Each live walk holds an affine canonical point W and the
coefficient pair (alpha, beta) with W = [alpha] P0 + [beta] T, alpha and beta
already multiplied by the canonicalizing root.

Distinguished points: the least significant 40 bits of the canonical x are 0
(theta = 2^-40). A distinguished point stores (x, seed, length) in a hash
table. A repeated x is a class collision.

Step of one walk:
  j = low bits of the canonical x, 11 bits, disjoint from the 40 DP bits.
  If the incomplete addition W + M_j is defined, set W to its canonical form
  and add (c_j, d_j) before folding the canonical multiplier into the
  coefficients.
  Look-ahead: if the partition of the sum is still j, discard that sum and
  add M_{j+1} instead. This is the (1 + 1/R) factor charged in section 4,
  not an uncharged retry.
  If the canonical x is already in the DP table, go to RESOLVE.
  If this walk's length exceeds 2^48, abandon it and open a fresh walk.
  The abandon probability per walk is about exp(-2^8) because the mean
  distance to a DP is 2^40. It does not move the 0.39 bound.

Stop after S steps in total across all walks, S the integer in section 3.
Then let every live walk finish its current trail up to the length cap.
That drain is inside the 2^100 unit reserve, not inside S.

RESOLVE. Re-walk the two colliding seeds, folding canonical multipliers into
the coefficients, until both sides stand on the stored point. Then
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

Let mean^2 = pi * N' / 2. Use pi < 355/113, so

    pi * N' / 2 < 355 * N' / 226.

Let T2 be the integer ceiling of the right-hand side and let M be the integer
ceiling of sqrt(T2). Then M is an integer upper bound on the ideal random-mapping
mean sqrt(pi * N' / 2). The experiment computes

    M   = 117165187599510394509384490721780424645
    T2 from (355 * N' + 225) // 226.

The charged length multiplies M by the two factors that are not left to the
heuristic gap:

    (1 + 1/R)       look-ahead, R = 2048,     = 2049/2048
    (1 + 1/(2R))    Teske r-adding penalty    = 4097/4096

    S = ceil(M * 2049 * 4097 / (2048 * 4096))
      = 117251015913075528483252403931634233610.

On a uniform random mapping, the probability of a collision by step t is
1 - exp(-t^2 / (2 N')). At t = sqrt(pi * N' / 2) this is 1 - exp(-pi/4) ≈ 0.544.
S is at least that mean, so the ideal-model probability at S is at least 0.544.
The step count t_39 at which the ideal law first reaches probability 0.39 solves
1 - exp(-t_39^2 / (2 N')) = 0.39, i.e. t_39 = sqrt(2 N' ln(1/0.61)). The ratio
of the ideal mean to t_39 is 1.2605. A walk whose collision time is up to 26.05%
slower than the ideal mean still has its 0.39 quantile inside S.

The two algebraic factors above are already outside that 26% gap. What remains
of the gap is the allowance for short fruitless cycles that survive look-ahead.
That is the heuristic, and it is the one the experiment is aimed at.

## 4. Cost

One reference unit is one batched incomplete affine G1 addition under the 6M
convention (3M formula + 3M amortized batch inversion, S = M, A = 0). The rules
price that operation at 1, not at 2298/2000. Extra instructions are divided by
C = 2000. Field multiplication costs 383 instructions, the pinned `blst_fp_mul`
difference (3,137,556 - 1,568,788) / 4096.

Per step, on top of the reference addition:

| Item | Instructions |
| --- | ---: |
| beta*x and beta*(beta*x) | 766 |
| Montgomery excess at B = 4096 | 57 |
| Two 6-limb minimum selections | 64 |
| Conditional y with even low bit | 32 |
| Partition bits and DP test | 16 |
| **E** | **935** |

The excess line is the part of Fermat inversion that the 3M amortization does
not cover. Inversion of one field element by the p-2 chain is 380 squarings and
228 extra multiplications (p - 2 has bit length 381 and Hamming weight 229),
608 multiplications in total. Montgomery's batch on B inversions costs
1 inversion + 3(B - 1) multiplications, hence 3 + 605/B multiplications per
inverse. The reference convention already includes 3. The difference at
B = 4096 is 605/4096 of a multiplication:

    ceil(605 * 383 / 4096) = 57.

So one step costs 1 + 935/2000 = 2935/2000 reference units. Across S steps,

    main = S * 2935 / 2000.

Side work, charged once and not per step:

* Partition table: 2048 jumps, two 255-bit scalar multiplications each, at most
  512 additions each, altogether 2^21 additions, plus the sign check of beta
  and the [r]T membership test. Bound 2^22 units. This is the preprocessing
  figure (preprocessing_log2 = 22). It depends on T and is inside total work.
* Distinguished-point inserts: about theta * S = 2^-40 * 2^126.46 ≈ 2^86.5
  records. Charge 2^11 instructions per insert (hash of a 48-byte x, compare,
  store). That is about 2^86.5 units.
* Resolve: two re-walks capped at length 2^48, each step at most 2 units,
  about 2^50 units.
* Final verification [x]P0 = T, one scalar multiplication.

The sum of those side items is far below 2^100 units. The claimed work is

    W = S * 2935 / 2000 + 2^100.

`experiment.py` evaluates log2(W) at 80 decimal digits and checks that the
ceiling to 10^-6 is exactly 127016229/1000000. Adding the side reserve moves
log2 by about 10^-8, which is why the ceiling sits one millionth above the
main term rather than on it. A larger side reserve of 2^108 would cross the
next millionth; 2^100 does not. The hash-table term 2^86.5 is inside 2^100
with more than a factor 2^13 to spare.

Peak memory is the DP table: at most 2^87 records of 128 bytes, which is 2^94
bytes, plus the 4096 live walks (a few hundred bytes each) and the partition
table (2048 affine points). memory_log2_bytes = 94. No nonuniform advice.
Advice bytes are 0, so the reported log2(max(1, advice_bytes)) is 0.

## 5. Experiment, and what it does not say

`experiment.py` checks beta, lambda and z^2 - 1 = lambda on the production
constants. It then runs the canonical walk, including coefficient folding,
look-ahead and the mu = 1 recovery formula, on the additive groups of prime
order 10009, 100003, 1000003 and 5000011 (each order is 1 mod 6, so the
sixth-root action is free). With seed 12648430 and 40 trials on each order,
every planted logarithm was recovered. Median collision time over the ideal
class mean was 1.055, 0.896, 0.844 and 1.096. The empirical 0.39 quantile over
that mean was 0.762, 0.748, 0.723 and 0.734, all under 1, while the production
budget is the mean times 2049/2048 * 4097/4096. The heaviest single trial was
2.67 means, inside the 6-mean bug tripwire and irrelevant to the 0.39 quantile.

The same recovery is run on the subgroup of order 4093 inside
y^2 = x^3 + 4 over F_12253 (curve order 12279 = 3 * 4093). There beta = 11770
and lambda = 3190, phi(G) = [lambda] G is checked, and 16 planted logarithms
are recovered from canonical affine points. The toy partition count is 32
rather than 2048; Teske's penalty is larger on the toy than in production, so
the toy ratio is the more pessimistic of the two.

This is evidence about the walk and the recovery formula at those orders. It
is not a measurement of production collision time, not a measurement of
blst, and not a proof that a random mapping on N' classes has the birthday tail.
The 26% gap is an extrapolation from that heuristic plus these ratios. A
walk that was systematically 30% slower than the ideal mean would put the 0.39
quantile past S; nothing in the sample looks like that, and nothing in the
sample reaches production size.

## 6. Heuristic and limitations

The random-mapping law for an r-adding walk on automorphism classes is the
standard heuristic for this kind of estimate. Look-ahead and the Teske
1 + 1/(2R) factor are charged in S. Residual short cycles are what the 26%
gap is for. The experiment's 0.39 quantiles on four prime orders sit between
0.72 and 0.76 of the ideal mean, against an ideal-model value near 0.79.

The attack is classical, cold, and sums every processor. Batching 4096 walks
is how the reference 6M inversion amortization is realized; it does not divide
the total work by 4096. There is no uncharged advice, no auxiliary powers, and
no pairing. Memory, preprocessing and the side reserve are as in section 4.
Success 0.39 is a lower bound under the heuristic, not a worst-case theorem.

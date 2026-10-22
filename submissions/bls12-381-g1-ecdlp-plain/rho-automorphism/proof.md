# bls12-381-g1-ecdlp-plain: parallel Pollard rho on automorphism classes

Status: ready for independent review; no qualification or emitted score.
This is a generic square-root attack. It is not a claim of novelty or of a
practical break. `claim.json` binds the protected target profile and
`classical-group-work-v2` by SHA-256. It uses no auxiliary powers, no pairing,
no advice and no quantum step.

Claimed total work: at most 2^128 target-group-operation equivalents
(score 128/1). Success probability is at least 0.4 under the
random-walk heuristic stated in "Heuristics" below.

## Group structure used

Let p and r be the protected field modulus and subgroup order, with
log2(r) = 254.85709 to five decimals. G1 = <P0> has prime order r on
E: y^2 = x^3 + 4 over Fp.

* Negation: -(x,y) = (x,-y).
* GLV endomorphism: phi(x,y) = (beta*x, y), where beta is the primitive cube
  root of unity in Fp
  `0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac`.
  Direct arithmetic on the profile constants gives beta^2 + beta + 1 = 0 in Fp
  and, with z = -0xd201000000010000 the BLS parameter,
  lambda = z^2 - 1 mod r = 0xac45a4010001a40200000000ffffffff, which satisfies
  lambda^2 + lambda + 1 = 0 mod r. The matching of this beta with this lambda
  is the identity phi(P0) = [lambda]P0, verified by one scalar multiplication
  during setup (the other cube root of unity would match r-1-lambda). That
  check is charged in setup. Both constants are functions of the public
  profile, not searched advice.

The group A = <-1, phi> is cyclic of order 6 and acts on G1 \ {O}. A point
fixed by any nontrivial element of A satisfies [mu - 1]P = O, where mu is a
nontrivial sixth root of unity mod r. Because mu - 1 is nonzero mod the prime
r, this forces P = O. Every non-identity orbit therefore has exactly 6
elements, so there are N' = (r-1)/6 orbits, and log2 N' = 252.27213.

Canonical representative: for an affine point (x,y), form the three
x-candidates x, beta*x and beta^2*x. Two multiplications suffice for this,
since beta^2*x = beta*(beta*x). Choose the smallest as a canonical integer
in [0,p) and then choose the y in {y, p-y} with even least significant bit.
The map applied is epsilon * lambda^k with epsilon in {±1} and k in {0,1,2}.
It is recorded only during coefficient re-walks (see below).

## Algorithm

```text
Input T. If T = O return 0. Check T on curve and [r]T = O (charged).
Fix R = 2048 partition points M_j = [c_j]P0 + [d_j]T, where c_j, d_j come
from fresh algorithm coins (uniform in Z_r).
Distinguished-point property: the low 40 bits of the canonical x are zero
(theta = 2^-40).
Run Mproc = 2^30 parallel trails. Each trail:
  draw a 128-bit seed s; (a,b) = PRG(s) uniform in Z_r^2; W = canon([a]P0 + [b]T)
  repeat: j = h(W) (11 bits of canonical x, disjoint from the DP bits);
          W = canon(W + M_j)     [batched affine addition, Montgomery trick
                                  across the 2^30 trails]
          if W is distinguished: insert (x(W), s, length) into the hash
          table; on a hit go to RESOLVE; start a new trail
          if the trail length exceeds Lmax = 2^45: abandon it and start a new one
  fruitless-cycle handling: look-ahead at j (if h(W+M_j) = j, use j+1),
  plus a periodic 2- and 4-cycle check with escape by doubling,
  following Bernstein-Lange-Schwabe "On the correct use of the negation map".
Stop new trails after K_main = 1.25 * K steps in total. Then drain every
active trail to its DP or to Lmax.
RESOLVE: re-walk both trails from their seeds and track coefficients
(alpha, beta_coef) with W = [alpha]P0 + [beta_coef]T. Apply
epsilon*lambda^k to both coefficients on each canonicalisation. A
class collision gives alpha1 + beta1*x = mu*(alpha2 + beta2*x) for one of the
six mu in A; this resolves directly because the re-walk records mu. Solve
x = (mu*alpha2 - alpha1) / (beta1 - mu*beta2) mod r; if the denominator is 0
continue the search. Verify [x]P0 = T and output x.
```

## Success argument

Model the iteration on the N' orbit classes as a random mapping. With
K total fresh class evaluations, the probability that some class repeats is
1 - exp(-K^2 / (2N')) (birthday bound, van Oorschot-Wiener parallel
collision search). Choose K^2 = 2N' ln(1/0.55), which gives probability
0.45 in the ideal model:

log2 K = (log2 N' + 1 + log2 ln(1/0.55)) / 2 = 126.26498.

Departures from the ideal model, all disclosed:

* r-adding walks with R = 2048 are slightly worse than random mappings.
  Teske's and Bernstein-Lange's analyses put the penalty at a factor of about
  1 + 1/(2R) in collision time. This is far below the budget allowance.
* Negation and order-6 fruitless cycles. With look-ahead and periodic cycle
  escape, the BLS-2011 measurements for R = 2048 show a total overhead of a
  few percent, including wasted cycle-check additions.
* Abandoned trails: the length-bound fraction is about exp(-32), negligible.

The algorithm runs 1.25 * K steps, a 25 % allowance that covers all three
effects with a wide margin. The supported success probability is therefore at
least 0.45 in the heuristic model. We claim 0.4, which is at least 0.39.

Detection: after K_main the drain sends every trail, including both members of
a merge, to a DP, except trails abandoned at Lmax (probability ~exp(-32)). A
merged pair always shares its DP, so a class collision becomes a hash-table
hit. Degenerate relations (denominator 0 mod r) occur with probability about
6/r per collision, because seeds give uniform independent coefficients.

## Cost accounting (units: one G1 reference operation = 2,000 instructions)

Per-step extra work, itemised in scalar AArch64 instruction equivalents, using
the pinned blst_fp_mul measurement of 383 instructions and 6-limb (384-bit)
word operations:

| Extra kernel | Bound (instructions) |
| --- | ---: |
| Two Fp multiplications beta*x and beta*(beta*x) | 766 |
| Three 6-limb comparisons to pick the minimal x | 120 |
| Select the matching y, optional p-y, parity bit | 80 |
| Extract 11 hash bits and the 40-bit DP mask | 40 |
| Load 96-byte affine M_j and loop/DP control | 80 |
| Amortised 2-cycle/4-cycle compare every 16 steps | 40 |
| Look-ahead hash of the point just produced | 20 |
| Rounding reserve | 54 |
| **Total extra per step** | **1,200** |

1,200 / C = 0.6 reference units. The walk addition itself is one 6M batched
incomplete affine addition (Montgomery inversion across the 2^30 live trails).
Same-x incomplete failures are fruitless cycles or useful collisions and are
covered by the 25 % allowance, not by a per-step inversion. Look-ahead that
rejects the current addend and retries has probability 1/R = 2^{-11} and is
absorbed by the same allowance.

| Component | Bound (log2 units) |
| --- | ---: |
| Main walk: 1.25 K additions at 1.0 unit | 126.587 |
| Main walk extra: 1.25 K steps at 0.6 unit | 125.907 |
| Combined main walk at 1.6 units/step | 127.265 |
| Drain: at most 2^30 * 2^45 steps at 1.6 units | 75.7 |
| Trail starts: 1.25 K * theta = 2^86.587 starts; each is two 255-bit Jacobian scalar multiplications, charged at the conservative 2^14 units | 100.587 |
| DP hash-table work: 2^86.587 inserts/lookups at 2^16 instructions (2^5.7 units) | 92.3 |
| Setup, partition table (2048 * 2 scalar multiplications), beta check, input checks | < 30 |
| RESOLVE: two re-walks of at most 2^45 steps, tracking coefficients (two mod-r multiplications per step, at most 4 units per step), six candidate solves, verification [x]P0 = T | < 49 |

The Rayleigh / random-mapping mean for a single collision is
sqrt(pi N'/2) = 2^126.462 class evaluations. The budget 1.25 K = 2^126.587
exceeds that mean; in the ideal model it already gives collision probability
1 - exp(-(1.25 K)^2 / (2 N')) = 0.607 before overhead. Multiplying by 1.6
units/step gives 2^127.265 for the main walk. Adding the 2^100.587 start cost
and every remaining line still leaves a total below 2^127.27 in this ledger.

A harsher substitution that charged a full extra reference unit per step and
also kept the 25 % overhead, 2.0 * 1.25 * sqrt(pi N'/2) = 2^127.784, remains
below the claimed 2^128. The claimed bound is 2^0.216 above that figure and
2^0.735 above the itemised 1.6-unit ledger, which covers rounding, a second
collision if the first relation is degenerate (probability 6/r), and the
drain, start and resolve tails.

Finite cap: if the main walk plus drain exceeds 1.4 K class evaluations of
live trails, the algorithm returns failure. 1.4 K = 2^126.750 class
evaluations, times 1.6 units, is 2^127.428, still under 2^128. The
probability of hitting this cap without a class collision is already folded
into the claimed success probability 0.4.

## Memory, preprocessing, advice

* DP table: at most 2^86.6 records of 72 bytes (48-byte x, 16-byte seed,
  8-byte length), with at most 2x hash-table slack, giving about 2^93.8 bytes.
* Trail state: 2^30 trails * (96-byte point + 16-byte seed + batch scratch of
  about 200 bytes), about 2^38.3 bytes.
* Partition table: 2048 * 96 bytes.

Peak memory is at most 2^94 bytes. Preprocessing, consisting of setup and the
partition table, is less than 2^30 units. It depends on T and is included in
the total. No nonuniform advice is used.

## Heuristics and limitations

* Random-mapping behaviour of the r-adding walk on automorphism classes. This
  is the standard heuristic behind every published rho estimate for
  pairing-friendly G1. It is supported by Teske's experiments and by
  Bernstein-Lange-Schwabe's negation-map measurements. It is not proved.
* The 25 % walk-overhead allowance is an extrapolation from published
  small-scale and moderate-scale measurements (up to about 2^60-step
  computations) to production size.
* Analytic bound only; no production-size execution is claimed.

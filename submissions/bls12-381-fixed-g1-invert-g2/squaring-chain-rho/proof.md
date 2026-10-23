# Fixed-G1 pairing inversion by a cyclotomic squaring-chain distinguished-point walk

Status: ready for exploratory review. This file alone claims no qualification and
emits no score by itself. The single scored number is the first declared claim in
`submission.json`, namely `time_log2 = 127.553439`.

Parent commit of the public tree this package revises: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

Normative inputs and outputs are those of `tracks/bls12-381-fixed-g1-invert-g2/target-profile.json`.
The reference operation is one dense Fp12 multiplication and costs 1.
Other scalar AArch64 instruction-equivalents cost 1/C with C = 20000.
No quantum step, no uncharged advice, no division of work by a processor count,
by a Frobenius-orbit size, by a negation class, or by any other automorphism orbit.

## Reduction

Let r be the protected subgroup order and Q0 the protected G2 generator.
The challenge supplies P in G1 with P != infinity and [r]P = infinity, and z in GT.
Compute g = e(P, Q0) with the pinned map; the library argument order is `pairing(Q0, P)`.
r is prime and the pairing is nondegenerate, so g = e(P0,Q0)^b for the unknown
b in 1..r-1 and g generates the order-r subgroup. There is a unique x in 0..r-1
with z = g^x. Bilinearity gives e(P, [x]Q0) = g^x = z, so Q = [x]Q0 solves the
instance. z = 1 forces x = 0 and Q = infinity, encoded as JSON null; the algorithm
returns it as soon as the twelve coefficients match 1, but the charged budget does
not shrink. The rest of this file is one way to obtain x.

## What is new relative to a plain walk

A plain distinguished-point walk pays one dense Fp12 multiplication per group
element it visits. Every element of GT is cyclotomic: r divides p^4-p^2+1, hence
GT is contained in GPhi_6 = {y in Fp12* : y^(p^4-p^2+1) = 1}. Squaring inside
GPhi_6 is the Granger-Scott formula (R. Granger, M. Scott, "On the computation of
squarings in the cyclotomic subgroup", 2010): three Fp2 squarings and two Fp2
multiplications, no general Fp12 multiplication. That is roughly a quarter of the
reference price, and it is the only cheap group operation this package uses.

So each block of the walk performs exactly one dense multiplication and then three
squarings, and it emits all four elements. The three extra elements are not free
samples: they are fixed power-of-two multiples of one element, and the accounting
below credits only the collision events that a fixed multiplier can actually
create. Nothing is divided by four.

## Block walk

Precompute, with one shared square chain, u_s = 2^(12 s) + 1 for s = 0..19 and
M_s = g^(u_s) * z = g^(2^(12 s)) * g * z. Each u_s is an integer smaller than r,
nonzero and distinct. The chain costs 228 squarings plus 21 products, and the
adder construction is charged at 249 reference products.

State (y, alpha, beta) obeys y = g^alpha * z^beta. Start from y = 1, alpha = 0,
beta = 0, i.e. the element y_{-1,3} of the notation below.

Block b = 0, 1, 2, ... reads s_b = H(y_{b-1,3}) mod 20, where H is the canonical
integer of the degree-0 coefficient of the element, and performs

```text
y_{b,0} = y_{b-1,3} * M_{s_b}      alpha_{b,0} = alpha_{b-1,3} + u_{s_b}      beta_{b,0} = beta_{b-1,3} + 1
y_{b,e} = (y_{b,e-1})^2            alpha_{b,e} = 2 alpha_{b,e-1}              beta_{b,e} = 2 beta_{b,e-1}     e = 1,2,3
```

with all exponent arithmetic modulo r. The invariant holds at every emitted
element: y_{b,e} = g^(2^e alpha_{b,0}) * z^(2^e beta_{b,0}).

A block therefore visits four elements of GT. The first is a genuine adding-walk
step with multiplier 1; the block-to-block map is A -> 8 A + u_s on exponents,
which is an affine map with a fixed invertible multiplier and twenty branches.

A point is distinguished when the lowest 20 bits of the degree-1 coefficient of
the element are zero. Every distinguished element is stored as (element, alpha, beta).
On a repeated distinguished element with stored (alpha1, beta1) and current
(alpha, beta), the pair of visits gives the linear congruence

```text
alpha + x beta = alpha1 + x beta1   (mod r)
x (beta - beta1) = alpha1 - alpha   (mod r)
```

One extended gcd inverts (beta - beta1) modulo r and yields x. Then Q = [x]Q0 is
built and checked with a second pairing; on a mismatch the walk continues. The walk
runs a fixed block budget and then stops. Early success does not reduce the charge.

## Which collisions exist, and which are credited

Take two distinct visits (b,e) and (b',e') with element equality
g^(2^e alpha_b) z^(2^e beta_b) = g^(2^e' alpha_b') z^(2^e' beta_b'),
writing alpha_b, beta_b for the block-entry exponents. Equality of the two
elements is the single congruence A_b = 2^(e'-e) A_b' where A_b = alpha_b + x beta_b.

* Same block, e != e': the congruence becomes 2^(e-e') A_b = A_b, i.e.
  2^k = 1 (mod r) for k in {+-1, +-2, +-3}. r exceeds 8, so this never happens.
  Within-block visits are impossible collisions and are credited with nothing.
* Different blocks, b != b': the ratio A_b / A_b' must equal one of the seven
  values 2^(-3), 2^(-2), 2^(-1), 1, 2, 4, 8. For independent uniform A_b, A_b'
  each of those seven events has probability 1/r and they are disjoint.

So a block pair carries seven usable collision events, not sixteen, and not one
per stored element. That is the whole credit taken from the squaring chain, and
it is why the chain is cheap in the ledger rather than in the sample count.
Any within-block repeat that the store reports is discarded: it cannot occur, and
if a degenerate exponent relation ever makes the linear congruence singular
(2^e beta_b = 2^e' beta_b' mod r), that candidate is skipped and the walk continues.

## Heuristic and success probability

Heuristic (one, disclosed). Model the block map as a random map on Z_r for the
purpose of collision counting: the seven ratio events of a block pair are
independent with probability 1/r each. Support: Teske's r-adding walks with at
least 16 branches reproduce the random-mapping collision constant on prime-order
groups in those experiments (E. Teske, "On random walks for Pollard's rho method",
Mathematics of Computation, 2001); a fixed invertible multiplier and a
coefficient-derived branch index are exactly the structure of those walks.
Limitation: those experiments are not a 255-bit theorem, and the multiplier 8 and
the squaring chain are outside their tested shape. The budget below is set from
the constant-1 exponential tail rather than from the sharp birthday constant, so
the numerical claim needs only that tail.

Let the walk run N blocks. The expected number of usable collision events is
Lambda = 7 N (N-1) / (2r), so the modeled probability of no collision is at most
exp(-Lambda). Choose N minimal with N (N-1) >= 2r/7; then Lambda >= 1 and the
modeled success on a uniform z, after a uniform nonidentity P has fixed g, is at
least 1 - exp(-1) > 0.632. The claim uses the lower number 0.5.
The inequality is uniform in z because every nonzero discrete log is an affine
translate of the same walk on a cyclic group of prime order, and z = 1 is solved
by the explicit identity return as well as by x = 0. Fresh algorithm coins: the
walk is deterministic given (P, z), so the probability is over the challenge
distribution only.

Distinguished-point count. The walk emits 4N elements. Under the Bernoulli model
the insertion count has mean 4N/2^20 and, by the standard Chernoff bound for the
upper tail, P(count > M) <= exp(-(4N/2^20)/3) for M = 2 (4N/2^20). That failure
probability is far below 10^(-100) here and does not pull the joint success
probability under 0.5. The work bound charges M insertions even when fewer occur.

Detection lag. A collision of block exponents makes every later block repeat, so
the next distinguished element is common to both visits. One average gap is
charged: 2^20 further elements, that is 262144 further blocks, whether or not the
lag is needed. If the heuristic is false and no block exponents ever collide, the
algorithm fails; that event is exactly what the 0.5 already excludes under the
model. There is no silent fallback with an uncharged baby-step table.

## Instruction ledger

C = 20000. One dense Fp12 product, including a squaring charged at the dense
price, costs 1 reference unit, matching the organizer rounding of the pinned blst
dense product to C. Internal instructions of that product are not charged again.
Surrounding work is charged in instructions and divided by C.

Per block, outside the one reference product:

| Piece | Instructions |
| --- | ---: |
| Three Granger-Scott squarings in GPhi_6, 5000 each | 15000 |
| Degree-0 coefficient modulo 20 for the branch index | 800 |
| Four distinguished tests on the degree-1 low limb | 800 |
| Doubling of alpha three times modulo r, 600 each | 1800 |
| Doubling of beta three times, 300 each | 900 |
| Addition into alpha and increment of beta at the block entry | 500 |
| Loop control and addressing | 1000 |
| Total R_block | 20800 |

Step cost = (C + R_block)/C = 40800/20000 = 2.04 reference units per block, and a
block carries seven collision events.

The 5000 per Granger-Scott squaring is an upper bound: three Fp2 squarings and
two Fp2 multiplications are at most 6 + 6 = 12 base-field products at 383
instructions each (the pinned blst `blst_fp_mul` count), plus coefficient
additions, modular reductions, loads and stores. Three of those squarings per
block are charged in full.

Each distinguished insertion: 8000 instructions (rehash, compare, copy of at most
1024 bytes). Candidate verification: at most 64 candidate solutions are checked,
each with a pairing at 4096 reference units, a scalar multiplication at 1024 and
65536 instructions of output encoding and extended gcd; a failing candidate does
not stop the walk and the budget covers the full run.

Initial pairing and the verification pairing: at most 4096 reference products
each. The ate parameter has absolute value 0xd201000000010000, bit length 64,
Hamming weight 6. Per Miller iteration charge one G2 doubling, one Fp12 squaring,
one line evaluation and one sparse product, each at most one dense unit, hence 4.
Sixty-three doublings after the leading bit plus 6 addition iterations at the same
price 4 is at most 276. The hard part of the final exponentiation uses an exponent
dividing (p^4 - p^2 + 1)/r, at most 4*381 - 254 = 1270 bits, hence at most 2540
dense products by binary powering. The easy part is the fixed cofactor
(p^6 - 1)(p^2 + 1): a constant number of products and Frobenius maps. Charging
3000 dense products covers the hard part, the easy part and those maps. Then
276 + 3000 + 820 = 4096, which is the charged pairing. Frobenius on Fp2 is
conjugation, a subtraction per coefficient, and sits in the 820.

Subgroup checks z^r = 1 and g^r = 1: two binary exponentiations, at most
2 * 2 * 255 = 2040 dense products, charged 2048. Decoding and the identity test:
charged 4096 reference units, an upper bound on coefficient compares and rejection
paths. G2 scalar [x]Q0: at most 510 Jacobian additions or doublings, each at most
20 Fp2 products, one Fp2 product at most 4 base-field products, one base-field
product at most 383 instructions, hence at most 20 * 4 * 383 / 20000 < 1.6
reference units per point operation; charged 2 each, that is 1020, and the ledger
uses 1024. Adder construction: 249 reference products, as scheduled above.

prep_ref = 4096 + 2048 + 4096 + 249 + 4096 + 1024 = 15609 for the whole run.
Preprocessing before the walk (decode, subgroup checks, initial pairing, adders)
is 4096 + 2048 + 4096 + 249 = 10489 reference units and is inside the total.

## Integer parameters and the bound

With r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001:

| Quantity | Value |
| --- | ---: |
| blocks N (minimal with N(N-1) >= 2r/7) | 122399667570890999749836939535834385982 |
| emitted elements 4N | 489598670283563998999347758143337543928 |
| lag blocks | 262144 |
| charged insertions M = 2 (4N/2^20) | 933835354392173765181251064573931 |
| R_block | 20800 |
| prep_ref | 15609 |
| total T | 4993913907575187927183468583070577105028640 |

```text
T = (N + 262144) * 40800 + M * 8000 + 15609 * 20000 + 65536 + 64 * (5120 * 20000 + 65536)
  = 4993913907575187927183468583070577105028640
```

The base-2 logarithm of T/20000, rounded upward to six decimal places, is
127.553439. The base-2 logarithm of 10489, rounded upward to six decimal places,
is 13.356590. Peak storage is M records of 1024 bytes plus 2^20 bytes of scratch
and 20 adders of 1024 bytes:

```text
bytes = M * 1024 + 1048576 + 20480 = 956247402897585935545601090124774400
```

The base-2 logarithm of that byte count, rounded upward to six decimal places, is
119.524868. Nonuniform advice is 0 bytes.

For comparison under the identical ledger convention, a plain walk (one dense
product and one element per block, R_block = 2500) needs the same
N(N-1) >= 2r budget and reaches 128.098471. The chain's 127.553439 is the
difference between paying 2.04 per block for seven collision events and paying
1.125 per block for one.

## Why the bound is an upper bound

Every accepted path performs at most one initial pairing, at most 64 verification
pairings with their scalars, the adder schedule, the two subgroup exponentiations,
N + 262144 blocks, at most M insertions, the extended gcd and the output encoding.
Paths that fail the input encoding are rejected before the walk; the challenge
distribution is supported on valid (P, z), and invalid inputs are outside that
distribution. The budget still covers a full run. Work of one processor is summed.
There is no second walker, no quotient by Frobenius, no negation map and no
processor count in the denominator, so there is no fruitless-cycle retry to forget.
Public profile constants and the writing of this argument are not charged.
Challenge sampling is not charged.

## Evidence limits

The toy experiment in `experiment.py` recovers discrete logarithms in a cyclic
group of prime order at most 4000 embedded in Z_p* with the same block update
(one multiplication by an adder, three squarings, seven credited ratios). It does
not measure Fp12 products, memory, cyclotomic squaring costs or the 255-bit tail,
and it cannot distinguish the seven-ratio accounting from a naive sample count.
No production walk is executed. Organizer acceptance is exploratory and is not a
proof that the random-mapping tail is a theorem. Finite-field descent in Fp12 is a
different algorithm and is not given a cost here.

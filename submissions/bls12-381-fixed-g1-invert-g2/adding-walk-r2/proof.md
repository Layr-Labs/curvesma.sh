# Fixed-G1 pairing inversion by a 20-adding walk, revision 2

Status: ready for exploratory review; no qualification or emitted score is claimed by this file alone.
Parent commit of the public tree this package revises: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.
This is a conservative generic accounting control. It is not a finite-field security estimate and not a practical attack.
The normative inputs and outputs are those of `tracks/bls12-381-fixed-g1-invert-g2/target-profile.json`.
`claim.json` binds that file and `tracks/classical-group-work-v2.json` by SHA-256.
The reference operation is one dense Fp12 multiplication and costs 1.
Other scalar AArch64 instruction-equivalents cost 1/C with C = 20000.
No quantum step, no uncharged advice, and no division of work by a processor count or by an automorphism-orbit size.

## Reduction

Let r be the protected subgroup order and let Q0 be the protected G2 generator.
The challenge gives P in G1, P ≠ ∞, [r]P = ∞, and z in GT.
Compute g = e(P, Q0) with the pinned map. The library argument order is pairing(Q0, P).
r is prime and the pairing is nondegenerate, so g = e(P0, Q0)^b for the unknown b in 1..r-1 and g generates the order-r subgroup.
There is a unique x in 0..r-1 with z = g^x. Bilinearity gives e(P, [x]Q0) = g^x = z.
Return Q = [x]Q0. The case z = 1 forces x = 0 and Q = ∞, encoded as JSON null.
The identity target is inside the budget below; the algorithm may return it as soon as the twelve coefficients match 1, but the charged budget does not shrink.

## Revision 2 walk

Revision 2 is a serial Teske r-adding walk whose every branch increments the exponent of z.
Useless same-branch collisions of a two-sided random walk are removed by that choice, and the surrounding instructions are charged outside the reference product.
Precompute, by one shared square chain,

```text
u_s = 2^(12*s) + 1 for s = 0..19
M_s = g^(u_s) * z = g^(2^(12*s)) * g * z
```

Each u_s is an integer smaller than r, nonzero, and distinct.
Schedule: 228 squarings produce g^(2^228); one product builds g*z; twenty products multiply that value by g^(2^(12*s)). Adder construction is 249 reference products.

State (y, α, β) obeys y = g^α * z^β, starting from y = 1, α = 0, β = 0.
Each step reads the canonical integer of the degree-0 Fp coefficient of y, reduces it modulo 20 to obtain s, and replaces

```text
y ← y * M_s
α ← α + u_s
β ← β + 1
```

The invariant is preserved because M_s = g^(u_s) * z.
A point is distinguished when the lowest 20 bits of the degree-1 coefficient are zero.
The store records (y, α, β).
On a repeated distinguished y with stored (α0, β0) and current (α, β),

```text
α + x β ≡ α0 + x β0  (mod r)
x (β - β0) ≡ α0 - α  (mod r)
```

β is the step counter, so β ≠ β0 whenever the two visits differ.
x is that quotient. One extended gcd inverts (β - β0) modulo r.
Then Q = [x]Q0 is built and checked with a second pairing. A mismatch continues the walk.
The walk runs a fixed iteration budget and then stops. Early success does not reduce the charged budget.

## Heuristic and success

Heuristic (one, disclosed). Model the map f(y) = y * M_(H(y)) as a random function on a domain of r group elements, in the sense used for Pollard–Teske adding walks: the probability that k_eff successive values are pairwise distinct is at most exp(-k_eff(k_eff-1)/(2r)), and the number of distinguished hits is not materially heavier than a Bernoulli sample of rate 2^(-20).
Support: Teske's r-adding walks with r ≥ 16 reproduce the random-mapping collision constant on the prime-order groups in those experiments (E. Teske, "On random walks for Pollard's rho method", Mathematics of Computation, 2001). The same leading constant is the classical Pollard mean √(π r / 2).
Limitation: those experiments are not a 255-bit proof. The budget below sits a fixed factor above that mean so the numerical claim only needs the exponential tail, not the sharp constant. Extrapolation to this r is part of the heuristic. No other heuristic is used. Coefficient bias of honest GT elements is absorbed by the same tail gap; the predicate and the partition read different coefficients.

Let s = isqrt(2r) and set

```text
k = s + 2 + 2^20
k_eff = k - 2^20 = s + 2
```

Then k_eff(k_eff - 1) ≥ 2r, so the modeled probability of no collision in the first k_eff steps is at most exp(-1) < 0.367879.
The final 2^20 steps are the distinguished-point detection lag of one average gap: a collision already inside the prefix is reported when the walk next stores a point, which the model places well inside one gap once a collision of walks exists, and the extra window is charged whether or not it is needed.
Modeled success on a uniform z, after a uniform nonidentity P has fixed g, is therefore at least 1 - exp(-1) > 0.632.
The claim uses the lower number 0.5.
Fresh algorithm coins: the walk is deterministic given (P, z). The probability is over the challenge distribution only. The inequality is uniform in z because every nonzero discrete log is an affine translate of the same adding walk on a cyclic group of prime order, and z = 1 is solved by the explicit identity return as well as by x = 0.

Distinguished-point count. Let M = k >> 19, i.e. twice k/2^20.
Under the Bernoulli model the insertion count has mean k/2^20 and

```text
P(count > M) ≤ exp(-(k/2^20)/3)
```

by the standard Chernoff bound for the upper tail.
That failure probability is far below 10^(-100) at this r and does not pull the joint success probability under 0.5.
The work bound charges M insertions even when fewer occur.
If the heuristic is false and the walk never collides, the algorithm fails; that event is exactly what the 0.5 already excludes under the model. There is no silent fallback with an uncharged baby-step table.

## Instruction ledger

C = 20000.
One dense Fp12 product, including a squaring charged at the dense price, costs 1 reference unit.
That matches the organizer rounding of the pinned blst dense product to C. Internal instructions of that product are not charged again.
Surrounding work is charged in instructions and divided by C.

Per walk step, outside the product:

| Piece | Instructions |
| --- | ---: |
| Pointers to y and M_s, 12 coefficients, address arithmetic | 1500 |
| Degree-0 coefficient modulo 20 | 800 |
| Distinguished test on the degree-1 low limb | 200 |
| 255-bit addition into α | 400 |
| Increment of β | 100 |
| Loop control | 1000 |
| Total R_step | 4000 |

Step cost = (C + R_step)/C = 24000/20000 reference units.
The product kernel is assumed to load its own operands. The 1500 instructions are only the addressing outside that kernel.

Each distinguished insertion: 8000 instructions (rehash, compare, copy of at most 1024 bytes).

Initial pairing and the verification pairing: at most 4096 reference products each.
The ate parameter has absolute value 0xd201000000010000, bit length 64, Hamming weight 6.
Per Miller iteration charge one G2 doubling, one Fp12 squaring, one line evaluation and one sparse product, each at most one dense unit, hence 4.
Sixty-three doublings after the leading bit, plus 6 addition iterations at the same price 4, is at most 276.
The hard part of the final exponentiation uses an exponent dividing (p^4 - p^2 + 1)/r, at most 4*381 - 254 = 1270 bits, hence at most 2540 dense products by binary powering. The easy part is the fixed cofactor (p^6 - 1)(p^2 + 1): a constant number of products and Frobenius maps. Charging 3000 dense products covers the hard part, the easy part and those maps. Then 276 + 3000 + 820 = 4096, which is the charged pairing. Frobenius on Fp2 is conjugation, a subtraction per coefficient, and sits in the 820.

Subgroup checks z^r = 1 and g^r = 1: two binary exponentiations, at most 2 * 2 * 255 = 2040 dense products. Charged 2048.
Decoding and the identity test: charged 4096 reference units, an upper bound on coefficient compares and rejection paths.
G2 scalar [x]Q0: at most 510 Jacobian additions or doublings. Each is at most 20 Fp2 products. One Fp2 product is at most 4 base-field products, and one base-field product is at most 383 instructions, hence at most 20 * 4 * 383 / 20000 < 1.6 reference units. Charged 2 per point operation would be 1020. The ledger uses 1024.
Extended gcd and output encoding: 2^16 instructions.

Adder construction: 249 reference products, as scheduled above.

Integer parameters for this r, with r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001:

| Quantity | Value |
| --- | ---: |
| k | 323839080949554914831837072806461860413 |
| M | 617674028300390081084894319165157 |
| prep_ref | 15609 |
| T | 7772142884181544359084738426509638283413536 |

prep_ref = 4096 + 2048 + 4096 + 249 + 4096 + 1024 = 15609.
T = k * 24000 + M * 8000 + prep_ref * 20000 + 65536.
These four lines are the whole attack: every step, every charged insertion, both pairings, both subgroup exponentiations, adder construction, the scalar, decoding and the gcd.
Preprocessing, included in T, is decode + subgroup + initial pairing + adders = 4096 + 2048 + 4096 + 249 = 10489 reference units.
No table built before the walk is treated as free advice. The twenty adders are instance-dependent and are inside prep_ref.

The base-2 logarithm of T/20000, rounded upward to six decimal places, is 128.191581.
The base-2 logarithm of 10489, rounded upward to six decimal places, is 13.356590.
Peak storage is M records of 1024 bytes, plus 2^20 bytes of scratch and 20 * 1024 bytes of adders:

```text
bytes = M * 1024 + 1069056 = 632498204979599443030931782826189824
```

The base-2 logarithm of that byte count, rounded upward to six decimal places, is 118.928545.
Nonuniform advice is 0 bytes.

## Why the bound is an upper bound

Every accepted path performs at most one initial pairing, one verification pairing, one scalar multiplication, the adder schedule, the two exponentiations, k steps and at most M insertions under the heuristic event already folded into the success probability.
Paths that fail the input encoding are rejected before the walk; the challenge distribution is supported on valid (P, z), and invalid inputs are outside that distribution. The budget still covers a full run.
Work of one processor is summed. There is no second walker. There is no quotient by Frobenius and no negation map, so there is no fruitless-cycle retry to forget.
Public profile constants and the writing of this argument are not charged. Challenge sampling is not charged.

## Evidence limits

The toy experiment in `experiment.py` recovers every discrete log in a cyclic group of order at most 4000 with the same v_s = 1 update. It does not measure Fp12 products, memory or the 255-bit tail.
No production walk is executed. Organizer acceptance is exploratory and is not a proof that the random-mapping tail is a theorem.
Finite-field descent in Fp12 is a different algorithm and is not given a cost here.


[You have received this identical output 3 times. Re-reading '/home/innolabs/deep-research/work/Holindauer--curvesmash--bls12-381-fixed-g1-invert-g2/submissions/bls12-381-fixed-g1-invert-g2/adding-walk-r2/proof.md:raw' will not change it — use a narrower selector (path:A-B), or proceed with the edit.]
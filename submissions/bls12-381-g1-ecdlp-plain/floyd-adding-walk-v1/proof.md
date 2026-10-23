## Instruction prices

Reference prices are the fixed constants of classical-group-work-v2. A submission cannot change them. H is a count of reference operations. R is every other scalar 64-bit instruction equivalent. Total work is W = H + ceil(R / C) when R is not divisible by C, and H + R/C otherwise. The integer used below is always the ceiling. A reference operation is never also unpacked into instructions.

Field-multiplication prices sit strictly above the published measurements where those measurements exist:

| Track family | Charged multiplication | Charged addition | Why the multiplication price is an upper bound |
| --- | ---: | ---: | --- |
| secp256k1 | 256 | 64 | The published 5x52 kernel measures 179 instructions. 256 is above that measurement. A 4-limb add-with-carry is well under 64 instructions. |
| P-256 | 1,024 | 128 | No P-256 kernel measurement is published. A 4-limb product plus Solinas reduction of 2^256 - 2^224 + 2^192 + 2^96 - 1 is a few hundred instructions even in a straight-line schedule. 1,024 is several times the published secp256k1 figure of 179. |
| BLS12-381 base field | 512 | 128 | The published blst_fp_mul measurement is 383 instructions. 512 is above it. A 6-limb add is well under 128 instructions. |

P-256 uses C = 1,000 only as the borrowed normalization unit. One P-256 addition is not one unit. BLS G1 uses C = 2,000. Dense Fp12 multiplication uses C = 20,000. A 256-bit modular multiplication of scalars is charged at 2,048 instructions, above a 4-limb schoolbook product.
# bls12-381-g1-ecdlp-plain: Floyd adding-walk bound

This is an analytic accounting argument for one fully specified classical algorithm. It does not execute a production discrete logarithm, and it does not say that the curve is broken in practice. The integer work bound is large. The memory bound for this algorithm is 2^20 bytes. No number-field sieve, GLV endomorphism, or negation map is used.

C = 2,000 is the BLS G1 normalization from the published 6-multiplication calibration.

## Algorithm

Let P be the generator fixed by the protected target profile and let r be its prime order. On input Q, recover x in {0, ..., r-1} with Q = [x]P.

If Q is the point at infinity, return 0. That path performs no walk. The work bound below is the complementary path.

Otherwise draw one fresh 64-bit seed s. This is the only algorithm coin. Set K = 2s+1, which is odd. Build 32 jump points as follows. Start a 64-bit linear congruential generator at s, with step state <- state * 6364136223846793005 + 1 (mod 2^64). For i = 0..31, take the next state as c_i = (state mod (r-1)) + 1 and the state after that as d_i = state mod r. Store R_i = [c_i]P + [d_i]Q in Jacobian coordinates, together with the coefficient pair (c_i, d_i).

The walk state is a Jacobian point X:Y:Z and coefficients (alpha, beta) such that the point equals [alpha]P + [beta]Q. Start both the tortoise and the hare at the identity, with both coefficients 0. One step replaces the state by its sum with R_h and adds (c_h, d_h) to the coefficients. The index is

h = ((X_coord * K) >> 371) & 31

when Z is nonzero, and h = 0 at infinity. X_coord is the Jacobian X word in {0, ..., p-1}, not an inverted affine coordinate. The shift 371 is bit_length(p) - 10, so the five extracted bits sit in the high half of a full-size field word. No field inversion is performed on the walk.

One Floyd iteration applies one step to the tortoise and two steps to the hare, then tests equality once. Run exactly

I = 3 * ceil(sqrt(r)) = 686,966,430,457,948,734,194,560,730,881,182,400,512

iterations. Do not stop early for accounting purposes: the bound charges every iteration even if a collision appears sooner. On the first equality whose beta coefficients differ, set

x = (alpha_t - alpha_h) * (beta_t - beta_h)^(-1) (mod r),

verify [x]P = Q, and return x if the check holds. If the betas are equal, or the check fails, finish the budget and fail. There is no second seed.

The curve coefficient a is 0.

## Group-operation bound

Instrumented Jacobian schedules, counting a squaring as a multiplication, are:

| Formula | Multiplications | Additions |
| --- | ---: | ---: |
| Doubling, a = 0 | 8 | 16 |
| Doubling, a nonzero | 10 | 17 |
| Mixed Jacobian plus affine | 11 | 15 |
| Complete Jacobian addition | 15 | 12 |

Every walk step is a complete Jacobian addition. The charged pair (16, 24) dominates all four rows, including the a-nonzero doubling used on P-256. One charged group operation is therefore

16 * 512 + 24 * 128 = 11,264

instructions. These instructions are R, not H. H = 0 on this track: the reference price is not identified with a Jacobian addition.

Binary scalar multiplication of a 255-bit scalar is at most 254 doubles and 254 additions. Each of the 64 jump scalar multiplications and the final verification is charged at 520 group operations, a spare of 12 above 508. The 32 sums that combine [c]P and [d]Q fit in that spare. Setup is 64 * 520 = 33,280 group operations and depends on Q, so it is online. Verification is 520 group operations.

ceil(sqrt(r)) = 228,988,810,152,649,578,064,853,576,960,394,133,504. Three applications per iteration give the walk

9 * ceil(sqrt(r)) = 2,060,899,291,373,846,202,583,682,192,643,547,201,536

group operations. Adding setup and verification gives

2,060,899,291,373,846,202,583,682,192,643,547,201,536 + 33,280 + 520 = 2,060,899,291,373,846,202,583,682,192,643,547,235,336

group operations.

## Other instructions

Per iteration, charged once:

- Equality. Jacobian equality is eight field multiplications (Z2^2, Z1^2, Z2^3, Z1^3, two X products, two Y products). Charged with eight extra additions: 5,120 instructions.
- Coefficients. Three steps, two scalar modular additions each, 256 instructions each: 1,536 instructions. A 4-limb modular addition is far under 256 instructions.
- Partition. Three steps, 256 instructions each: 768 instructions. One 64-bit multiply-and-shift is far under 256 instructions.

A flat reserve of 2^40 instructions covers control, loads of the 32 jumps, and output.

## Total

| Component | Instructions |
| --- | ---: |
| Group operations | 23,213,969,618,035,003,625,902,596,217,936,916,058,824,704 |
| Equality | 3,517,268,123,944,697,519,076,150,942,111,653,890,621,440 |
| Coefficients | 1,055,180,437,183,409,255,722,845,282,633,496,167,186,432 |
| Partition | 527,590,218,591,704,627,861,422,641,316,748,083,593,216 |
| Reserve | 1,099,511,627,776 |
| Sum R | 28,314,008,397,754,815,028,563,015,083,999,913,711,853,568 |

W = ceil(R / C) = 14,157,004,198,877,407,514,281,507,541,999,956,855,927.

The declared score is the rational 133378640/1000000 = 133.378640. It is the ceiling of log2(W) at a resolution of 10^(-6), so 2^(score) >= W. Review does not tighten a bound. A one-unit-in-the-sixth-decimal excess would be a worse score, not an undercount. The experiment recomputes W and checks the inequality with an 80-digit logarithm.

## Success probability

The claimed success probability is 0.39. It is a lower bound, not the success rate of a random-mapping model.

For a uniform random mapping on a set of size n = r, the rho length L satisfies P(L > k) about exp(-k^2 / (2n)). At k = 3 * sqrt(n), the exponent is -4.5 and the model failure probability is about 0.011, so the model success probability is about 0.989. The classical expected rho length is sqrt(pi n / 2), about 1.253 sqrt(n). The budget is about 2.39 times that expectation.

The adding walk is not proved to be a uniform random mapping. The claim uses 0.39, which remains valid unless this walk is more than about 2.4 times slower to collide than the model. Two further gaps are disclosed. First, a collision with equal beta coefficients does not reveal x, and that event is a failure of this single-seed procedure. Second, the production index is the high half of the Jacobian X word. The toy measurements below use canonical integer representatives, not Jacobian coordinates.

The declared experiment runs this walk, with a high-half 64-bit mix, on the additive groups of prime order 10,007, 100,003 and 1,000,003. On the fixed 30 secrets it recovers at least 90 percent of them inside the same 3 * ceil(sqrt(n)) budget. That is evidence about those groups only. It is not a production complexity proof. A partition that uses the low bits of the representative does not mix, and it is not the partition specified here.

## Memory, preprocessing, advice

Peak memory is 2^20 bytes: 32 Jacobian jumps, their coefficients, the two walk states, and scratch. log2 of that memory is exactly 20.

Preprocessing that does not depend on Q is at most ceil(2^30 / C) = 536,871 reference units, for a fixed load of the public generator. log2 of that quantity, rounded up at 10^(-6), is 19.034216. The jump scalar multiplications depend on Q and are inside W, not in this preprocessing figure.

Nonuniform advice is 0 bytes. No distinguished-point table and no endomorphism constants are stored.

## What is not claimed

- No production discrete logarithm is computed.
- GLV and the negation map are not used. Their fruitless-cycle analysis is easy to get wrong, and this bound does not need them.
- A distinguished-point walk is not used. An affine-x distinguished predicate would require a field inversion on every step.
- Transfer into GT, and any number-field sieve, is not used on the plain BLS track. A transfer would have to charge every pairing and every descent instruction. That cost is not estimated here.
- The bound is not a security level. It is the cost of this one algorithm under the stated instruction ceilings.

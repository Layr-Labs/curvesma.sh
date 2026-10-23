# p256-ecdlp: Floyd adding walk, revision 2

This is an analytic accounting argument for one classical algorithm. It does not compute a production discrete logarithm. Revision 2 replaces the rejected walk. Four points in that package do not survive: the recovery numerator had the wrong sign, the program's partition was not the proof's partition, the jump coefficients were 64-bit and correlated, and a hash of a raw Jacobian triple is not a function of the point. The integer W below includes a Fermat inversion on every walk step so that the stored point is the canonical affine representative. That extra work changes the six-decimal score.

C = 1,000. C = 1,000 is the borrowed normalization unit. It is not a P-256 addition. Curve coefficient a = -3, so the Jacobian setup price uses 28 additions, above the a-nonzero doubling. H = 0. Every field operation and every scalar instruction is in R, and W = ceil(R / C).

The order r is the protected profile order. The iteration budget is

I = 3 * ceil(sqrt(r)) = 1,020,847,100,643,973,146,639,479,902,990,327,873,536.

ceil(sqrt(r)) = 340,282,366,881,324,382,213,159,967,663,442,624,512. The budget is charged in full even if a collision appears earlier. There is no second seed.

## Recovery formula

The state is stored so that the current point equals [alpha]P + [beta]Q. With Q = [x]P the discrete log of the state is alpha + beta * x modulo r. When the tortoise and the hare hold the same point,

alpha_t + beta_t * x = alpha_h + beta_h * x (mod r),

hence x * (beta_t - beta_h) = alpha_h - alpha_t (mod r). If the beta coefficients differ,

x = (alpha_h - alpha_t) * (beta_t - beta_h)^(-1) (mod r).

The opposite numerator, alpha_t - alpha_h, equals -x. The point [(-x)]P equals Q only when 2x = 0. On a prime-order group that forces x = 0. The input x = 0 is the identity, and the identity is returned before any walk step. A write-up that solves for (alpha_t - alpha_h) therefore succeeds on the identity challenge and fails on every other challenge that the walk actually collides.

The program returns the residue from the displayed formula, checks it against the secret, and checks that the opposite numerator is a different residue whenever the secret is nonzero.

## One partition function

Let MASK = 2^64 - 1. For a 64-bit word z,

mix64(z) = the following, every sum and product reduced modulo 2^64:

- z = z + 0x9E3779B97F4A7C15
- z = (z xor (z >> 30)) * 0xBF58476D1CE4E5B9
- z = (z xor (z >> 27)) * 0x94D049BB133111EB
- return z xor (z >> 31)

The index is

h = ((mix64(word) * K) >> 32) & 31, where K = (2 * seed + 1) mod 2^64.

word is already reduced modulo 2^64. The same lines are in the declared program. Two fixed checks are part of that program: partition(1, 1) = 5, and partition of the canonical word of the affine point (1, 1) is 14.

The canonical word depends on the group, because the elements are different objects. The function of the word does not.

- Elliptic tracks. Each walk step is an affine addition, so the stored point is either the identity or a pair (x, y) of fully reduced field elements. The word is 0 at the identity and otherwise (x xor ((y << 1) mod 2^64)) mod 2^64, after x and y have been reduced modulo 2^64. A raw Jacobian triple is not used. Jacobian coordinates of one point are not unique, and a hash of a non-canonical triple is not a function of the group element. On a prime-order curve of size about 10^4 that non-canonical hash failed to collide inside this budget. The affine representative is what the sample below actually runs.
- Pairing tracks. An element of GT is stored as a reduced Fp12 vector, which is already unique. The word is the low 64 bits of its first base-field coefficient.
- The declared integer sample. The word is the group element in 0 .. n-1. That is Teske's reduction of an r-adding walk to the integers modulo n, with this partition in place of an abstract hash.

The rejected packages did not use this function. Their programs used ((point * 0x9E3779B97F4A7C15 + seed) >> 32) & 31. The elliptic proofs used a product by (2*seed+1) shifted by bit_length(p) - 10. The pairing proofs used a shift of 27. Those three rules are different functions. This package has one rule.

## Jump coefficients

One 64-bit seed s is the only algorithm coin. For each of the 64 coefficients a 512-bit integer is built from eight mix64 outputs. Block index t, for t = 0 .. 63, is the integer whose eight 64-bit limbs, from most significant to least, are

mix64((s + (8*t + j + 1) * 0x9E3779B97F4A7C15) mod 2^64) for j = 0 .. 7.

Reducing a uniform 512-bit integer modulo r biases each residue by less than r / 2^512, which is below 2^(-256) for these orders. The odd-numbered blocks, starting at t = 0, become c_i = 1 + (block mod (r - 1)), so c_i is nonzero. The following block becomes d_i = block mod r. The jump is [c_i]P + [d_i]Q on an elliptic track, and base^(c_i) * z^(d_i) on a pairing track.

The rejected packages drew these coefficients from a 64-bit linear congruential generator. Those residues are supported on only 2^64 values and are serially correlated. They are not uniform modulo a 255-bit prime.

The derivation is charged twice, once for the pass the algorithm keeps and once as a spare pass, so the bound does not depend on stopping early. One mix64 is charged at 40 scalar instructions. One 512-bit reduction is charged at 16,384 scalar instructions. The total is 2 * 64 * (8 * 40 + 16,384) = 2,138,112 instructions. It is an online line of its own. It is not placed in the 2^40 reserve, and the coefficients are not embedded in the program, so the advice byte count stays 0.

## Affine steps, and what is still Jacobian

A walk step adds one of the 32 jumps in affine coordinates and stores the canonical sum. Doubling or addition uses one inversion and at most four field multiplications and eight field additions. The inversion is the binary expansion of p - 2. The exponent p - 2 has 382 multiplications when every squaring is counted as a multiplication: bit length minus 1, plus Hamming weight minus 1. Each of those multiplications is charged at the field-multiplication price, plus 32 scalar instructions for the bit test and the loop. The inversion line is therefore

382 * 1024 + 382 * 32 = 403,392

instructions. A step is that inversion plus 4 * 1024 + 8 * 128 = 408,512 instructions. The P-256 prime is 2^256 - 2^224 + 2^192 + 2^96 - 1. p - 2 has bit length 256 and Hamming weight 128, so the product count is 255 + 127 = 382. No published kernel is treated as a P-256 multiplication; the multiplication price stays 1,024.

Three steps per iteration, across the whole budget, are

3,062,541,301,931,919,439,918,439,708,970,983,620,608 * 408,512 = 1,251,084,872,334,812,274,239,961,642,391,154,460,821,815,296

instructions.

The 32 jumps are built by Jacobian scalar multiplication and normalized once, at the end, into affine form. A 255-bit binary scalar multiplication is at most 508 group operations. Each of the 64 scalar multiplications (two per jump) and the final verification [x]P = Q is charged at 520 Jacobian group operations. One Jacobian group operation is charged at 16 field multiplications and 28 field additions, which is 19,968 instructions. That pair sits strictly above a complete Jacobian addition and above doubling, including the a-nonzero doubling on P-256. Setup and verification are

(64 * 520 + 520) * 19,968 = 674,918,400

instructions. They depend on Q, so they are online. Normalizing the 32 stored jumps is 32 inversions and 32 * 6 field multiplications:

13,105,152 instructions.

Equality of two affine points is two comparisons. It is charged at the more expensive projective test, eight field multiplications and eight field additions per iteration, which is 9,408,126,879,534,856,519,429,446,785,958,861,682,507,776 instructions over the budget. Coefficient updates are three steps, two scalar additions, 256 instructions each: 1,568,021,146,589,142,753,238,241,130,993,143,613,751,296 instructions. The partition, including mix64 and the canonical-word load, is charged at 256 instructions per step: 784,010,573,294,571,376,619,120,565,496,571,806,875,648 instructions. A 2^40 reserve covers control flow and the output write.

| Component | Instructions |
| --- | ---: |
| Affine walk steps | 1,251,084,872,334,812,274,239,961,642,391,154,460,821,815,296 |
| Jacobian setup and verification | 674,918,400 |
| Jump normalization | 13,105,152 |
| Equality | 9,408,126,879,534,856,519,429,446,785,958,861,682,507,776 |
| Coefficient updates | 1,568,021,146,589,142,753,238,241,130,993,143,613,751,296 |
| Partition | 784,010,573,294,571,376,619,120,565,496,571,806,875,648 |
| Coefficient derivation, two passes | 2,138,112 |
| Reserve | 1,099,511,627,776 |
| Sum R | 1,262,845,030,934,230,844,889,248,450,873,604,138,126,739,456 |

W = ceil(R / C) = 1,262,845,030,934,230,844,889,248,450,873,604,138,126,740.

The declared score is the rational 139857658/1000000 = 139.857658. It is a ceiling of log2(W) at resolution 10^(-6). The declared program recomputes W from the profile prime and the order and checks the ceiling with an 80-digit logarithm.

## Why the budget is enough for a claim of 0.39

The claimed success probability is 0.39. It is a lower bound. It is not the measured rate, and it is not 1.

Edlyn Teske, Mathematics of Computation 70 (2001), 809-836, studies r-adding walks. Her L is the number of iterating steps until her match finder reports a collision, divided by the square root of the group order. For a uniform random function, with her measured delay factor 1.13 in front of the classical rho length sqrt(pi |G| / 2), she obtains L_0 = 1.416. Corollary 5.1 says that if the partition hash behaves as an independent hash, then for every r at least 16 the average L on a prime-order group is at most 1.45, regardless of the size of the order. Her Table 3, on elliptic-curve subgroups, gives average L = 1.419 for r = 20. This walk uses r = 32.

L counts adding-steps of her iterator, not Floyd iterations. One Floyd iteration applies the iterator once to the tortoise and twice to the hare, so it costs three adding-steps. A random-function rho length near sqrt(pi n / 2), about 1.253 square roots, therefore corresponds to about 3.76 square roots of adding-steps when the collision is detected by Floyd. The charged budget is 3 * ceil(sqrt(n)) Floyd iterations, which is 9 * ceil(sqrt(n)) adding-steps. That is a little more than twice the Floyd expectation just quoted, and about six times Teske's L_0.

The hash in this package is a fixed mix, not a freshly sampled independent hash, so Corollary 5.1 is not a theorem about these constants. The sample below is the direct evidence. It uses the partition and the coefficient law specified above, on fresh seeds and fresh secrets drawn from that law.

Exponent sample, master seed 12648430, 80 draws on each order. The length column is the Floyd iteration of the first useful collision, divided by sqrt(order). Failures are omitted from the length columns and counted in the success column.

| Group order | Successes | Mean | Median | 90th percentile | Maximum |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10,007 | 80/80 | 1.053 | 1.065 | 1.809 | 2.479 |
| 100,003 | 79/80 | 1.126 | 0.942 | 2.119 | 2.855 |
| 1,000,003 | 80/80 | 1.024 | 0.996 | 1.743 | 2.709 |

The same master and the same coefficient law on the curve y^2 = x^3 + 2 over the prime field 10099, group order 9967, generator (3, 7278), using canonical affine words: 40/40 successes, mean 1.058, median 0.881, 90th percentile 1.753, maximum 2.735. Each success also satisfied [recovered]G = Q in the group law.


Every nonzero success was checked against the secret. The numerator (alpha_t - alpha_h) was also evaluated. It differed from the secret on every nonzero success. On the curve of order 9,967 the first recorded success is secret 2,444, recovered 2,444, opposite numerator 7,523. Those last two residues sum to 9,967, so the opposite numerator is -secret. Iteration 72 divided by sqrt(9,967) is 0.721, which sits inside the length sample rather than at its edge.

A uniform random function on a set of size n fails to have rho length at most 3 * sqrt(n) with probability about exp(-4.5), near 0.011. The measured failure counts are smaller than one trial in twenty. The claim stays at 0.39 because the hash is not proved independent and because a single seed with equal beta coefficients is defined to fail. Both of those gaps are already much smaller than the distance from the sample down to 0.39. Raising the claim toward the sample would ask the review to accept a sharp constant this text does not prove.

Identity inputs return 0 before the walk and count as successes. They are not the source of the 0.39. For a uniform secret the identity has probability 1/r.

## Memory, preprocessing, advice

Peak memory is 2^20 bytes. Thirty-two affine jumps, two walk points, and the coefficient table fit in that cap. log2 of the cap is 20.

Preprocessing that does not depend on the challenge is a fixed load of the public generator, charged at ceil(2^30 / C) = 1,073,742 reference units. The log2 ceiling is 20.034216. Jump construction depends on Q and is inside W, not in this figure.

Nonuniform advice is 0 bytes. The factors of nothing are stored. The mix constants are algorithm text, not a searched table. The reporting rule log2(max(1, advice_bytes)) is therefore 0.

## Not claimed

No production discrete logarithm is computed. GLV, a negation map, distinguished points, and a number-field sieve are not used. The figure 0.39 is not a theorem that this mix is an independent hash.

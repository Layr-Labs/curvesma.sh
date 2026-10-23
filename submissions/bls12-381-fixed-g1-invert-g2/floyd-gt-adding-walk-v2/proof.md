# bls12-381-fixed-g1-invert-g2: GT adding walk, revision 2

This is an analytic accounting argument. It does not evaluate a production pairing inversion. The pairing reduction is unchanged. The walk that recovers the GT discrete log is the revision-2 adding walk: the recovery numerator is (alpha_h - alpha_t), the partition function is the one function specified below, and the jump exponents are 512-bit blocks reduced modulo r. The coefficient derivation is a named line of R. It adds 107 reference units. The six-decimal ceiling of log2(W) does not move.

The input is a non-identity point P in the profile's source group and an element z in GT. Return Q in the profile's output group with the pairing of (P, Q) equal to z. The identity z = 1 returns the point at infinity and does not walk.

The public generator Q0 of the output group is fixed by the profile. Compute base = e(P, Q0) with one pairing. Solve z = base^a in the prime-order subgroup of GT. Return Q = [a]Q0.

C = 20,000 instructions per reference dense Fp12 multiplication. Walk multiplications are H. They are not unpacked into instructions a second time.

## Recovery formula

The state is stored so that the current point equals [alpha]P + [beta]Q. With Q = [x]P the discrete log of the state is alpha + beta * x modulo r. When the tortoise and the hare hold the same point,

alpha_t + beta_t * x = alpha_h + beta_h * x (mod r),

hence x * (beta_t - beta_h) = alpha_h - alpha_t (mod r). If the beta coefficients differ,

x = (alpha_h - alpha_t) * (beta_t - beta_h)^(-1) (mod r).

The opposite numerator, alpha_t - alpha_h, equals -x. The point [(-x)]P equals Q only when 2x = 0. On a prime-order group that forces x = 0. The input x = 0 is the identity, and the identity is returned before any walk step. A write-up that solves for (alpha_t - alpha_h) therefore succeeds on the identity challenge and fails on every other challenge that the walk actually collides.

The program returns the residue from the displayed formula, checks it against the secret, and checks that the opposite numerator is a different residue whenever the secret is nonzero.

The symbols alpha and beta are the exponents of base and of z. The same linear relation holds in a multiplicative group: a collision gives a = (alpha_h - alpha_t) * (beta_t - beta_h)^(-1) modulo r. The opposite numerator is -a.

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

A reduced Fp12 element is a canonical encoding, so this track does not invert on the walk. The elliptic tracks do, because a Jacobian triple is not canonical. That inversion is not part of W here.

## Jump coefficients

One 64-bit seed s is the only algorithm coin. For each of the 64 coefficients a 512-bit integer is built from eight mix64 outputs. Block index t, for t = 0 .. 63, is the integer whose eight 64-bit limbs, from most significant to least, are

mix64((s + (8*t + j + 1) * 0x9E3779B97F4A7C15) mod 2^64) for j = 0 .. 7.

Reducing a uniform 512-bit integer modulo r biases each residue by less than r / 2^512, which is below 2^(-256) for these orders. The odd-numbered blocks, starting at t = 0, become c_i = 1 + (block mod (r - 1)), so c_i is nonzero. The following block becomes d_i = block mod r. The jump is [c_i]P + [d_i]Q on an elliptic track, and base^(c_i) * z^(d_i) on a pairing track.

The rejected packages drew these coefficients from a 64-bit linear congruential generator. Those residues are supported on only 2^64 values and are serially correlated. They are not uniform modulo a 255-bit prime.

The derivation is charged twice, once for the pass the algorithm keeps and once as a spare pass, so the bound does not depend on stopping early. One mix64 is charged at 40 scalar instructions. One 512-bit reduction is charged at 16,384 scalar instructions. The total is 2 * 64 * (8 * 40 + 16,384) = 2,138,112 instructions. It is an online line of its own. It is not placed in the 2^40 reserve, and the coefficients are not embedded in the program, so the advice byte count stays 0.

Each jump is base^(c_i) * z^(d_i). Each exponentiation is charged at 510 dense Fp12 multiplications, above the 254 squarings and 254 multiplications of a 255-bit binary schedule. Two exponentiations per jump give 64 * 510 = 32,640 multiplications. One further 510-multiplication exponentiation verifies base^a = z.

## Budget and pairings

The iteration budget is

I = 3 * ceil(sqrt(r)) = 686,966,430,457,948,734,194,560,730,881,182,400,512.

Three group multiplications per iteration give

9 * ceil(sqrt(r)) = 2,060,899,291,373,846,202,583,682,192,643,547,201,536

dense Fp12 multiplications. The budget is charged in full. A collision with equal exponents of z fails this seed. There is no restart.

Each of the two pairings is charged at 9,343 dense Fp12 multiplications. The Miller loop on the 64-bit parameter X = 0xd201000000010000, Hamming weight 6, is 3 * (64 - 1) + 2 * (6 - 1) = 199. The final exponent is charged as 2 * 12 * bit_length(p) = 9,144. The sum is 9,343. Two pairings are 18,686. This is several times a tight ate schedule.

The output scalar multiplication is in G2. One G2 addition is charged as 16 products in Fp2, and one Fp2 product as 4 base-field multiplications and 4 base-field additions: 4 * 512 + 4 * 128 = 2,560 instructions, times 16 = 40,960 instructions per addition. 520 such additions are 21,299,200 instructions.

| Component | Count |
| --- | ---: |
| Walk multiplications H | 2,060,899,291,373,846,202,583,682,192,643,547,201,536 |
| Setup multiplications H | 32,640 |
| Verification exponentiation H | 510 |
| Two pairings H | 18,686 |
| Total H | 2,060,899,291,373,846,202,583,682,192,643,547,253,372 |
| Equality instructions | 1,406,907,249,577,879,007,630,460,376,844,661,556,248,576 |
| Coefficient instructions | 1,055,180,437,183,409,255,722,845,282,633,496,167,186,432 |
| Partition instructions | 4,220,721,748,733,637,022,891,381,130,533,984,668,745,728 |
| Output scalar-multiplication instructions | 21,299,200 |
| Coefficient derivation, two passes | 2,138,112 |
| Reserve | 1,099,511,627,776 |
| Total R | 6,682,809,435,494,925,286,244,686,790,013,241,927,245,824 |

W = H + ceil(R / C) = 2,395,039,763,148,592,466,895,916,532,144,209,349,735.

The declared score is the rational 130815246/1000000 = 130.815246.

## Why the budget is enough for a claim of 0.39

The claimed success probability is 0.39. It is a lower bound. It is not the measured rate, and it is not 1.

Edlyn Teske, Mathematics of Computation 70 (2001), 809-836, studies r-adding walks. Her L is the number of iterating steps until her match finder reports a collision, divided by the square root of the group order. For a uniform random function, with her measured delay factor 1.13 in front of the classical rho length sqrt(pi |G| / 2), she obtains L_0 = 1.416. Corollary 5.1 says that if the partition hash behaves as an independent hash, then for every r at least 16 the average L on a prime-order group is at most 1.45, regardless of the size of the order. Her Table 3, on elliptic-curve subgroups, gives average L = 1.419 for r = 20. This walk uses r = 32.

L counts adding-steps of her iterator, not Floyd iterations. One Floyd iteration applies the iterator once to the tortoise and twice to the hare, so it costs three adding-steps. A random-function rho length near sqrt(pi n / 2), about 1.253 square roots, therefore corresponds to about 3.76 square roots of adding-steps when the collision is detected by Floyd. The charged budget is 3 * ceil(sqrt(n)) Floyd iterations, which is 9 * ceil(sqrt(n)) adding-steps. That is a little more than twice the Floyd expectation just quoted, and about six times Teske's L_0.

The hash in this package is a fixed mix, not a freshly sampled independent hash, so Corollary 5.1 is not a theorem about these constants. The sample below is the direct evidence. It uses the partition and the coefficient law specified above, on fresh seeds and fresh secrets drawn from that law.

The declared program is the exponent walk. It applies partition to the integer representative, with the coefficient law above, for 80 draws on each order and master seed 12648430.

| Group order | Successes | Mean | Median | 90th percentile | Maximum |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10,007 | 80/80 | 1.053 | 1.065 | 1.809 | 2.479 |
| 100,003 | 79/80 | 1.126 | 0.942 | 2.119 | 2.855 |
| 1,000,003 | 80/80 | 1.024 | 0.996 | 1.743 | 2.709 |

This is not a sample of Fp12 encodings. It is the integer walk Teske reduces an r-adding walk to, with this partition as the hash. The elliptic packages run the same partition on canonical affine points and are where the group-law check lives. Production on this track feeds partition the low 64 bits of the first reduced base-field coefficient and no other word.


Every nonzero success was checked against the secret. The numerator (alpha_t - alpha_h) was also evaluated. It differed from the secret on every nonzero success. On the curve of order 9,967 the first recorded success is secret 2,444, recovered 2,444, opposite numerator 7,523. Those last two residues sum to 9,967, so the opposite numerator is -secret. Iteration 72 divided by sqrt(9,967) is 0.721, which sits inside the length sample rather than at its edge.

A uniform random function on a set of size n fails to have rho length at most 3 * sqrt(n) with probability about exp(-4.5), near 0.011. The measured failure counts are smaller than one trial in twenty. The claim stays at 0.39 because the hash is not proved independent and because a single seed with equal beta coefficients is defined to fail. Both of those gaps are already much smaller than the distance from the sample down to 0.39. Raising the claim toward the sample would ask the review to accept a sharp constant this text does not prove.

Identity inputs return 0 before the walk and count as successes. They are not the source of the 0.39. For a uniform secret the identity has probability 1/r.

## Memory, preprocessing, advice

Peak memory is 2^20 bytes. Preprocessing that does not depend on the challenge is ceil(2^30 / C) = 53,688 reference units, log2 ceiling 15.712313. The pairings and the jump elements depend on the challenge and are inside W. Nonuniform advice is 0 bytes. The coefficient derivation is recomputed from the seed. It is not a stored table.

No Miller inversion and no number-field sieve is included. No production pairing is evaluated.

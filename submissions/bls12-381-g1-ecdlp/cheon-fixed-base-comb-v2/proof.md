# bls12-381-g1-ecdlp: Cheon with a fixed-base comb, revision 2

This is an analytic bound at the scoring width n = 2^20. It does not compute the production secret. Revision 2 keeps the same index search as the previous package and corrects the resource accounting. Two gaps are closed, and the other lines of the sum are named so that none of them is an uncharged step.

The declared score is unchanged at six decimal places: 122820313/1000000 = 122.820313. The integer W grew, and the growth is below one unit in the sixth decimal of log2(W). The score is that ceiling, not a rounded display of the old integer.

## Profile data and what is not free

The protected profile publishes `parameter_x` = -0xd201000000010000 and `subgroup_order` r. It does not publish the prime factorization of r - 1. A profile constant may be read without treating the read as advice. Work the attack then performs is charged. Embedding a searched table in the program does not make the search free, so this program does not embed the primes and does not embed the primitive root.

Let X be the absolute value of `parameter_x`. The attack checks X^4 - X^2 + 1 = r, at 16 modular multiplications of 255-bit integers, and then factors the three 64-bit integers X, X - 1 and X + 1. Their prime powers, with every prime of X taken twice, are multiplied back together and compared with r - 1. That comparison is the proof that the factors are the factorization of r - 1. It is not assumed.

For each of the three integers the attack performs 2^32 division attempts, one for every integer in {0, ..., 2^32 - 1}. (2^32)^2 = 2^64 is at least X + 1, so the attempts reach a square root of each input. Each attempt is charged at 32 scalar instructions (a 64-bit division is one instruction; 32 covers the loop). That is

3 * 2^32 * 32 = 412,316,860,416

instructions. Rebuilding the product is charged at 64 modular multiplications, 2,048 instructions each. A further 2^30 instructions cover writing the factor list into scratch. The total is 413,390,766,080 instructions. Divided by C = 2,000 this is 206,695,384 reference units. It is inside W and inside preprocessing, because it does not depend on s.

The factor list lives in the 2^20-byte scratch region while the root test runs, and it is discarded. It is not part of the program text. Nonuniform advice is 0 bytes. The reporting rule for that field is log2(max(1, advice_bytes)) = log2(1) = 0.

The primitive root is found by testing candidates 2, 3, ... against the freshly computed primes. The budget is 64 candidates, each a 255-bit exponentiation charged at 510 modular multiplications (254 squarings and 254 multiplications) and 2,048 instructions per multiplication: 66,846,720 instructions. The value 7 is whatever candidate succeeds. It is not a stored constant.

## Comb tables stay Jacobian

Each comb entry is built by Jacobian doubling and Jacobian addition. The previous package stored those entries at 96 bytes, which is an affine point (two 48-byte field elements), and the batch inversion in the ledger covered only the baby/giant search list.

Two schedules were costed at the same group-operation price.

- Jacobian schedule. Store each table point as (X : Y : Z), 3 * 48 = 144 bytes. No batch inversion of the tables. Online steps are complete Jacobian additions. The charged operation is 16 field multiplications and 24 field additions, which is 11,264 instructions. The implemented complete-addition formula performs 16 multiplications; 24 additions sit above the 12 additions in that formula, and also above mixed addition (11 multiplications, 15 additions) and doubling.
- Affine schedule. After the same Jacobian construction, Montgomery-normalize every table point: 3 field multiplications per point, 4 field additions per point, and one Fermat inversion of 2 * 381 multiplications. Then store 96 bytes. Online steps could use mixed addition, but this comparison does not give them a cheaper price: the 16-by-24 charge already dominates mixed addition, so the only difference is the added normalization, 973,555,660,975,280,180,349,468,062,118,912 instructions.

Jacobian work is W = 9,388,543,893,536,700,382,679,492,617,175,715,396. Affine work is 9,388,544,380,314,530,870,319,582,791,909,746,455. The affine integer is larger by 486,777,830,487,640,090,174,734,031,059. Both ceilings at resolution 10^(-6) equal 122.820313. The claim uses the smaller integer. Memory does not enter the score. The Jacobian tables are still counted in peak memory, at 144 bytes, so the memory bound is not helped by a representation the construction does not produce.

Comb size: 2 bases, 3 digits, 2^96 slots, including the zero digit, so 475,368,975,085,586,025,561,263,702,016 points and 68,453,132,412,324,387,680,821,973,090,304 bytes.

## The rest of the sum

C = 2,000 is the BLS G1 entry of the fixed price table. H = 0. A Jacobian addition is not counted again as one reference operation.

| Line | Instructions | What it pays for |
| --- | ---: | --- |
| Group operations | 10,075,513,001,272,716,798,894,549,308,343,427,414,016 | 894,488,015,027,762,499,901,859,846,266,284,394 operations at 11,264 instructions. Includes comb construction, two online additions per phase-1 scalar, and phase 2. |
| Scalar and field setup | 915,955,240,610,598,312,259,414,307,841,582,874,624 | One 2,048-instruction multiplication per baby index and per giant index, plus 3,096 multiplications itemized below. |
| Search-list normalization | 1,831,910,481,221,196,624,518,828,615,683,153,458,176 | Montgomery batch inversion of the baby/giant list only: 7 multiplications per point and one Fermat inversion of 762 multiplications, plus 4 additions per point. |
| Radix sort | 5,495,731,443,663,589,873,556,485,847,049,509,535,744 | 447,243,769,829,393,707,157,917,142,500,769,792 records, 24 passes, 512 instructions per record per pass, plus bucket overhead. |
| Scan | 457,977,620,305,299,156,129,707,153,920,788,267,008 | 1,024 instructions per record. |
| Primitive-root search | 66,846,720 | 64 exponentiations, not an embedded root. |
| Factorization | 413,390,766,080 | The 64-bit trial division above. |
| Input loads | 576 | Reading A_1 and A_d, two instructions per byte of a 144-byte point. |
| Reserve | 1,099,511,627,776 | Loop branches and the 32-byte output encoding. |
| Total R | 18,777,087,787,073,400,765,358,985,234,351,430,790,720 | |

W = ceil(R / 2000) = 9,388,543,893,536,700,382,679,492,617,175,715,396.

The field-setup block of 3,096 multiplications is exactly:

- 20 squarings for zeta^(2^20), because the exponent is a power of two.
- 510 multiplications for beta^m.
- 510 multiplications for zeta^N.
- 3 * 510 multiplications for the Fermat inverses of beta, gamma and zeta^(k0).
- 510 multiplications for s = zeta^(k0 + k1*N).
- 16 multiplications to form k0 + k1*N.

20 + 6 * 510 + 16 = 3,096. A 255-bit binary exponentiation is at most 254 squarings and 254 multiplications, so 510 is an upper bound. These are not group operations and they are not inside the reserve.

Phase 2 is the rectangle of side 2^10. It is charged as (2 * 2^10 + 4) binary scalar multiplications, each 254 doubles and 254 additions, not as one addition per step. That is 1,042,416 group operations. Successive addition would be cheaper. The charge is the larger one.

Window 96 is the minimum of W over windows 16, 24, 32, 48, 64, 80, 96, 100 and 104. Past width 85 the digit count stays 3, so a wider window only grows the table.

## Search list, sort, and memory

Phase 1 stores m + v_max = 447,243,769,829,393,707,157,917,142,500,769,792 Jacobian points, then batch-normalizes those points only. The product buffer of that inversion is reused as the second radix-sort buffer. Peak memory is still charged as if the comb tables and both sort buffers were live together:

- search records: 2 * 447,243,769,829,393,707,157,917,142,500,769,792 * 128 bytes
- comb tables: 6 * 2^96 = 475,368,975,085,586,025,561,263,702,016 Jacobian entries, 144 bytes each. Two bases, three digits, and 2^96 points in each digit. 6 * 2^96 * 144 bytes is the comb line inside peak memory.
- scratch: 2^20 bytes, holding sort buckets and the live factor list

The sum is 114,494,473,529,457,201,356,814,469,302,171,205,632 bytes. Its log2 ceiling is 126.428546. This is the storage the algorithm writes. It is not a claim that the storage is on hand, and it is not part of the score.

Preprocessing that can be done from public data is the factorization, the root search, and the comb tables. Counting both comb tables overstates preprocessing, because the table for A_d depends on the challenge. The resulting figure is 2,677,278,067,682,020,495,961,244,437,561 reference units, log2 ceiling 101.078610. It is included in W. Phase-1 additions are inside W as well. There is no amortization across challenges.

## Algebra and success

d = 2^20 divides r - 1, which the attack checks after the factorization rather than by citing a factor list. N = (r - 1) / d, m = ceil(sqrt(N)) = 223,621,884,914,696,853,578,958,571,250,384,896, v_max = ceil(N / m) = 223,621,884,914,696,853,578,958,571,250,384,896. Every k0 in {0, ..., N - 1} is u + m*v with 0 <= u < m and 0 <= v < v_max. The same covering holds for k1 on the side-2^10 rectangle, because d is a square. The final scalar multiplication checks the candidate against A_1. The search is deterministic. The success probability on uniform s is 1.

## What the declared program checks

The program recomputes m, v_max, the group-operation count, the factorization instruction count, both W integers, and the score ceiling. It checks (2^32)^2 >= X + 1. It factors X, X - 1 and X + 1 by a faster method than the charged trial division, multiplies the prime powers back to r - 1, and checks that 7 is a primitive root for those primes. It recovers scalars with the index arithmetic on the additive groups of order 101, 211, 10007 and 65537, and it checks a width-8 comb against binary multiplication for one 64-bit secp256k1 scalar. It does not build a table of size 2^96 and it does not execute 2^32 division attempts. The trial-division cost is the closed form above.

## Not claimed

No production secret is recovered. GLV, the negation map, and a number-field sieve are not used. The removal of the per-step exponentiation is the counted comb, not a citation of an uncosted square-root improvement.

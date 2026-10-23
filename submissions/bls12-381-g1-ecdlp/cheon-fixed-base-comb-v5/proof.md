# bls12-381-g1-ecdlp: Cheon with a fixed-base comb, revision 5

Analytic bound at the scoring width n = 2^20. It does not compute the production secret. The public source used for the price table is commit 6423c2728dc03d5f2d647c7f22adf1dac0643a66. That commit does not change the reference prices or the definition of W from the tree these bounds were first written against. C = 2000. H = 0: a G1 addition is not a reference operation here, because the reference operation is a 6M affine addition and this schedule counts a different formula. W = ceil(R / C).

The declared score is 122.152158 = 122152158/1000000. It is the ceiling of log2(W) at resolution 10^(-6). W = 5908313823446922466800834637125629747.

## Field schedule

BLS12-381 Fp is a 6-limb Montgomery field, 64-bit limbs, CIOS multiplication. Squaring calls that multiplication. A field doubling, tripling, or multiplication by four or eight is a chain of modular additions, and each of those additions is reduced. The conditional subtraction of the modulus is executed on every add, sub, and product.

The declared program counts instructions by running the routine. One multiplication is 777 instructions. One addition is 86. One subtraction is 85. One inversion by the binary method on p - 2 is 476271 instructions, which is to_mont(1) plus one squaring per bit of p - 2 plus one multiplication per set bit, and four control instructions per bit. A second run of the multiplication on a different pair of inputs produces the same count. The carry identity (2^64 - 1)^2 + (2^64 - 1) + (2^64 - 1) = 2^128 - 1 keeps each mac inside two machine words, so the six-instruction mac (mul, umulh, adds, adc, adds, adc) is not an approximation.

## Group law

The curve is y^2 = x^3 + 4, so a = 0. The routines are the EFD formulas, executed on the profile generator:

- Doubling is dbl-2009-l. Leaf count on a non-infinity point: 7 multiplications (2M + 5S), 9 additions, 5 subtractions, plus a zero test of Z that loads and compares six limbs. The counter delta is 6651.
- Complete addition is add-2007-bl. It is not a complete formula in the Edwards sense: it is the general Jacobian sum, and the program does not feed it two equal points. The counter delta is 13568.
- Mixed addition, with the table point affine, is madd-2007-bl. The counter delta is 9756.

The same run checks that these three routines agree with an affine slope computation on the profile generator, and that [r]G is the point at infinity. The executable is this BLS12-381 Jacobian arithmetic. It is not an affine addition on a different curve.

## Window

A scalar is 255 bits. Window w uses t = ceil(255 / w) digits and two bases (A_1 and A_d). Each base stores t tables of 2^w affine points. Construction is Jacobian: w(t - 1) doubles and t(2^w - 1) additions per base, then one batch inversion. Online, each baby or giant is t - 1 mixed additions. The batch inversion is one field inversion plus (6N - 2) multiplications: 2(N - 1) for the product tree and 1S + 3M per point, with squaring calling multiplication.

| window | digits | W |
| --- | ---: | ---: |
| 16 | 16 | 34269830241292378114121822085618776672 |
| 24 | 11 | 23361554695153465596540222983370221097 |
| 32 | 8 | 16816589367470118085991264143027651220 |
| 48 | 6 | 12453279149014553078989411607364722367 |
| 64 | 4 | 8089968930558989417062562087075104270 |
| 80 | 4 | 8089968930647142942692282991594175104 |
| 85 | 3 | 5908313823446922466800834637125629747 |
| 96 | 3 | 5908318154319413473526287401974489275 |
| 100 | 3 | 5908383149142532050275632265754147466 |
| 104 | 3 | 5909423066312429278265150086228677721 |
| 112 | 3 | 6192280536524475291413997255300892662 |
| 127 | 3 | 9310929637274393488974388304632757588076 |

The minimum on this list is window 85, with t = 3. 85 * 3 = 255, so the digits cover the scalar with no padding bit. Wider windows with the same digit count only grow the tables.

## The sum

m = ceil(sqrt(N)) = 223621884914696853578958571250384896, v_max = ceil(N / m) = 223621884914696853578958571250384896, N = (r - 1) / 2^20. Every k0 in {0, ..., N - 1} is u + m v with 0 <= u < m and 0 <= v < v_max. points = m + v_max = 447243769829393707157917142500769792. comb points = 2 * 3 * 2^85 = 232113757366008801543585792.

| Line | Instructions | What it pays for |
| --- | ---: | --- |
| Jacobian construction | 3,149,319,459,942,007,419,343,374,205,788 | 340 doubles at 6651 and 232,113,757,366,008,801,543,585,786 additions at 13568 |
| Online mixed additions | 8,726,620,436,911,130,014,065,279,284,475,020,181,504 | 894,487,539,658,787,414,315,834,285,001,539,584 mixed additions at 9756 |
| Phase 2 | 8,551,459,656 | 2052 binary scalar multiplications, 254 doubles and 254 mixed additions each |
| Comb normalization | 1082114336840333032796197437021 | Batch inversion of every table point |
| Search-list normalization | 2085050454944633462770209718338589245021 | Batch inversion of the baby and giant list |
| Scalar index updates | 167716413686022640184218928437788672000 | One 4-limb Montgomery multiplication, 375 instructions, per baby and per giant |
| Field setup | 873701 | The six exponentiations below |
| Radix sort | 815772636168814121856040867921404137472 | 48 passes, 38 instructions per record per pass, plus 768 bucket instructions per pass |
| Scan | 21467700951810897943580022840036950016 | 12-limb equality, both coordinates, once per record |
| Root search | 101995328 | 64 candidates, 12 prime tests each, exact Hamming weight |
| Trial division | 90194313216 | 3 * 2^32 * 7 |
| Rebuild and polynomial | 18881 | Multiply the prime powers back, and check X^4 - X^2 + 1 = r |
| Loads and Montgomery setup | 3192 | Affine input loads, four conversions, and the one-time n0' Newton iteration |
| Total R | 11816627646893844933601669274251259492796 | |

There is no flat reserve. Every term that scales with the number of records is a product of that number and a per-record counter delta.

W = ceil(11816627646893844933601669274251259492796 / 2000) = 5908313823446922466800834637125629747.

The field-setup block is the measured exponentiation schedule, not a rounded 510:

- zeta^(2^20): 8709
- beta^m, with the real bit length and Hamming weight of m: 51472
- zeta^N, likewise: 139315
- three inverses in the scalar field, the measured inversion: 475560
- the unknown discrete log, charged as a 255-bit all-ones exponent: 192645
- 16 multiplications to assemble k0 + k1 N: 6000

The sum of those six lines is 873701.

## Factorization and the root

The profile publishes parameter_x and r. It does not publish the primes. The attack checks X^4 - X^2 + 1 = r by the scalar-field arithmetic counted above, then factors X, X - 1, and X + 1. Each receives 2^32 trial divisions. (2^32)^2 = 2^64 is at least X + 1. One attempt is udiv, msub, cmp, b.ne, add, cmp, b.lo: 7 instructions. The loop does not stop early in the bound. The product of the prime powers, with the primes of X taken twice, is compared with r - 1. The program recomputes the factorization by a faster method and checks that product; the faster method is not the attack and its cost is not subtracted. The primes it finds are X = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 906349, 254760293], X - 1 = [19, 125527, 2508409, 2529403], X + 1 = [3, 11, 10177, 859267, 52437899]. The root search tests 64 candidates against the 12 distinct primes. Candidate 7 survives. It is not a stored constant. Advice is 0 bytes, so the advice logarithm is 0.

## Memory and preprocessing

Peak storage is the affine comb point (96 bytes), the batch-inversion product (48 bytes), two 128-byte sort buffers for every search record, a 256-bucket histogram, and 65536 bytes of scratch: 114494405109749170093132055902473488384 bytes, log2 ceiling 126.428545. Preprocessing is the factorization, the root search, both comb tables, and their normalization: 2115716898391170226114933986 reference units, log2 ceiling 90.773206. Both tables are included even though the A_d table depends on the challenge. That figure is inside W. Phase-1 additions are inside W and are not called preprocessing.

## Success

The two rectangles cover every residue. The final scalar multiplication checks the candidate. The search is deterministic on every non-identity challenge the profile can produce, and the identity is the scalar 0, returned before the search. The success probability on a uniform challenge is 1.

## What the declared program checks

It remeasures the field and group counters, rebuilds this sum, and requires W and the score numerator to match the claim. It checks [r]G = infinity with dbl-2009-l and add-2007-bl. It checks the factor product and a small smooth-order index search. It does not build a table of size 2^85 and it does not execute 2^32 divisions.

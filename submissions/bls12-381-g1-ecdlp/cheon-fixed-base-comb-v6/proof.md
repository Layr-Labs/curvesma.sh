# bls12-381-g1-ecdlp: Cheon with a fixed-base comb, revision 6

Analytic bound at the scoring width n = 2^20. It does not compute the production secret. The public source used for the price table is commit 6423c2728dc03d5f2d647c7f22adf1dac0643a66. That commit does not change the reference prices or the definition of W from the tree these bounds were first written against. C = 2000. H = 0: a G1 addition is not a reference operation here, because the reference operation is a 6M affine addition and this schedule counts a different formula. W = ceil(R / C).

The declared score is 122.209170 = 122209170/1000000. It is the ceiling of log2(W) at resolution 10^(-6). W = 6146471130971250810599119934906972741.

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

A scalar is 255 bits. Window w uses t = ceil(255 / w) digits and two phase-1 bases, P0 and A_d = [secret^d]P0. Each base stores t tables of 2^w affine points. Construction is Jacobian: w(t - 1) doubles and t(2^w - 1) additions per base, then one batch inversion. Online, each baby or giant is t - 1 mixed additions. The batch inversion is one field inversion plus (7N - 3) multiplications. The products are 3(N - 1): N - 1 to form the prefix, then two multiplications on each of the N - 1 backward steps. Each point then pays 1S + 3M, and squaring calls multiplication, which is 4N. The sum is 7N - 3. N = 1 has no product and costs four multiplications plus the inversion. The declared program runs that schedule on one point and on four points and requires the counter delta to equal the formula.

| window | digits | W |
| --- | ---: | ---: |
| 16 | 16 | 34507987548726530263183412964001251329 |
| 24 | 11 | 23599712002587617745601813861895275876 |
| 32 | 8 | 17054746674904270235052855048106827845 |
| 48 | 6 | 12691436456448705228052314722087807359 |
| 64 | 4 | 8328126237993141623456633546546050606 |
| 80 | 4 | 8328126238085052433201236137450821813 |
| 85 | 3 | 6146471130971250810599119934906972741 |
| 96 | 3 | 6146475646434412443338049210907097422 |
| 100 | 3 | 6146543411470233331339958032950979111 |
| 104 | 3 | 6147627652043367539370499185653085337 |
| 112 | 3 | 6442541087935872123677692720625964619 |
| 127 | 3 | 9707766893228181420290013000472865098441 |

The minimum on this list is window 85, with t = 3. 85 * 3 = 255, so the digits cover the scalar with no padding bit. Wider windows with the same digit count only grow the tables. The list penalty is the same at every window, so rounding W up does not hide a missing product.

## The sum

m = ceil(sqrt(N)) = 223621884914696853578958571250384896, v_max = ceil(N / m) = 223621884914696853578958571250384896, N = (r - 1) / 2^20. Every k0 in {0, ..., N - 1} is u + m v with 0 <= u < m and 0 <= v < v_max. points = m + v_max = 447243769829393707157917142500769792. comb points = 2 * 3 * 2^85 = 232113757366008801543585792.

| Line | Instructions | What it pays for |
| --- | ---: | --- |
| Jacobian construction | 3,149,319,459,942,007,419,343,374,205,788 | 340 doubles at 6651 and 232,113,757,366,008,801,543,585,786 additions at 13568 |
| Online mixed additions | 8,726,620,436,911,130,014,065,279,284,475,020,181,504 | 894,487,539,658,787,414,315,834,285,001,539,584 mixed additions at 9756 |
| Phase 2 | 8,551,459,656 | 2052 binary scalar multiplications, 254 doubles and 254 mixed additions each |
| Comb normalization | 1262466726313721871595563596628 | One inversion plus (7N - 3) multiplications on every table point |
| Search-list normalization | 2432558864102072373231911338061687372628 | The same batch on the baby and giant list |
| Scalar index updates | 167716413686022640184218928437788672000 | One 4-limb Montgomery multiplication, 375 instructions, per baby and per giant |
| Field setup | 873701 | The six exponentiations below |
| Radix sort | 944578841879679509517521004961625837568 | 48 passes. Each record pays 44 instructions, and each pass pays 768 bucket instructions |
| Scan | 21467700951810897943580022840036950016 | 12-limb equality, both coordinates, once per record |
| Root search | 101995328 | 64 candidates, 12 prime tests each, exact Hamming weight |
| Trial division | 90194313216 | 3 * 2^32 * 7 |
| Rebuild and polynomial | 18881 | Multiply the prime powers back, and check X^4 - X^2 + 1 = r |
| Loads and Montgomery setup | 3192 | Affine input loads, four conversions, and the one-time n0' Newton iteration |
| Total R | 12292942261942501621198239869813945480106 | |

There is no flat reserve. Every term that scales with the number of records is a product of that number and a per-record counter delta.

The 44 instructions in one radix pass are a counted schedule, not a round allowance. The histogram update is add, ldrb, lsl, ldr, add, str. The scatter is ldrb, lsl, ldr, add, str, add, then sixteen post-indexed loads and sixteen post-indexed stores of the 128-byte record. The prefix sum is ldr, str, add on each of 256 buckets, which is 768. There are 48 passes because the sort key is the 48-byte affine x coordinate. Trial division is the seven instructions udiv, msub, cmp, b.ne, add, cmp, b.lo, produced by the same counter.

W = ceil(12292942261942501621198239869813945480106 / 2000) = 6146471130971250810599119934906972741.

The field-setup block is the measured exponentiation schedule, not a rounded 510:

- zeta^(2^20): 8709
- beta^m, with the real bit length and Hamming weight of m: 51472
- zeta^N, likewise: 139315
- three inverses in the scalar field, the measured inversion: 475560
- the unknown discrete log, charged as a 255-bit all-ones exponent: 192645
- 16 multiplications to assemble k0 + k1 N: 6000

The sum of those six lines is 873701.

Phase 2 stays binary. A 255-bit all-ones scalar starts from the leading one, then doubles and mixed-adds on each of the other 254 bits. That is 254 doubles and 254 mixed additions, not a rounded 508.

## Factorization and the root

The profile publishes parameter_x and r. It does not publish the primes. The attack checks X^4 - X^2 + 1 = r by the scalar-field arithmetic counted above, then factors X, X - 1, and X + 1. Each receives 2^32 trial divisions. (2^32)^2 = 2^64 is at least X + 1. One attempt is udiv, msub, cmp, b.ne, add, cmp, b.lo: 7 instructions. The loop does not stop early in the bound. The product of the prime powers, with the primes of X taken twice, is compared with r - 1. The program recomputes the factorization by a faster method and checks that product; the faster method is not the attack and its cost is not subtracted. The primes it finds are X = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 906349, 254760293], X - 1 = [19, 125527, 2508409, 2529403], X + 1 = [3, 11, 10177, 859267, 52437899]. The root search tests 64 candidates against the 12 distinct primes. Candidate 7 survives. It is not a stored constant. Advice is 0 bytes, so the advice logarithm is 0.

## Memory and preprocessing

Peak storage is the affine comb point (96 bytes), the batch-inversion product (48 bytes), two 128-byte sort buffers for every search record, a 256-bucket histogram, and 65536 bytes of scratch: 114494405109749170093132055902473488384 bytes, log2 ceiling 126.428545. The search-list prefix is one 48-byte field element per denominator. It is live during normalization, before the two sort buffers. A Jacobian point is 144 bytes, so the point plus one prefix element is 192 bytes, and the 256 bytes already charged per record cover that peak. The time score does not add those bytes a second time. Preprocessing is the factorization, the root search, both comb tables, and their normalization: 2205893093127864645514617065 reference units, log2 ceiling 90.833422. Both tables are included even though the A_d table depends on the challenge. That figure is inside W. Phase-1 additions are inside W and are not called preprocessing.

## Bases and coverage

Phase 1 solves secret^d = beta^k0 with beta = zeta^d of order N = (r - 1) / d. The two comb bases are P0 and A_d. A comb on A_1 in phase 1 would solve secret = zeta^(d k0). That is the subgroup of index d, a 1/d fraction of the nonzero scalars, so A_1 is the wrong phase-1 base. Phase 2 solves secret = zeta^k0 * gamma^k1 with gamma = zeta^N of order d, walking A_1 against P0. The rectangles cover every index: every k0 in {0, ..., N - 1} and every k1 in {0, ..., d - 1}. The final scalar multiplication checks the candidate. The identity is the scalar 0, returned before the search. The success probability on a uniform challenge is 1.

The declared program enumerates all 65536 nonzero scalars of the toy field r = 65537 with d = 256. Each scalar is recovered from A_d against P0 and then from A_1 against P0, and zeta^(k0 + k1 N) equals the scalar. The same run counts the scalars fixed by the exponent N and requires that count to equal N, which is exactly the image an A_1 phase-1 comb would hit. It also checks that production beta = 7^(2^20) has order N, using the factored form of r - 1, and that both production rectangles cover.

## What the declared program checks

It remeasures the field and group counters, runs the batch schedule on one point and on four points, rebuilds this sum, and requires W and the score numerator to match the claim. It checks [r]G = infinity with dbl-2009-l and add-2007-bl. It checks the factor product and the coverage argument above. It does not build a table of size 2^85 and it does not execute 2^32 divisions. Progress goes to the diagnostic stream. The only standard output is one JSON object, and the process exits 0.

# bls12-381-g1-ecdlp: Cheon with a fixed-base comb, revision 7

Analytic bound at the scoring width n = 2^20. It does not compute the production secret. The public source used for the price table is commit 6423c2728dc03d5f2d647c7f22adf1dac0643a66. That commit does not change the reference prices or the definition of W from the tree these bounds were first written against. C = 2000. H = 0: a G1 addition is not a reference operation here, because the reference operation is a 6M affine addition and this schedule counts a different formula. W = ceil(R / C).

The declared score is 122.798221 = 122798221/1000000. It is the ceiling of log2(W) at resolution 10^(-6). W = 9245870455897073182711296040491334261.

## Field schedule

BLS12-381 Fp is a 6-limb Montgomery field, 64-bit limbs, CIOS multiplication. Squaring calls that multiplication. A field doubling, tripling, or multiplication by four or eight is a chain of modular additions, and each of those additions is reduced. The conditional subtraction of the modulus is executed on every add, sub, and product.

The declared program counts instructions by running the routine. One multiplication is 777 instructions. One addition is 86. One subtraction is 85. One inversion by the binary method on p - 2 is 476271 instructions, which is to_mont(1) plus one squaring per bit of p - 2 plus one multiplication per set bit, and four control instructions per bit. A second run of the multiplication on a different pair of inputs produces the same count. The carry identity (2^64 - 1)^2 + (2^64 - 1) + (2^64 - 1) = 2^128 - 1 keeps each mac inside two machine words, so the six-instruction mac (mul, umulh, adds, adc, adds, adc) is not an approximation.

Limb loops inside CIOS and field addition are unrolled. The limb count is the field constant, and the measured schedule already lists every limb instruction. Loops whose trip count depends on the input, including the search list, the scalar bits, and the radix passes, pay induction separately.

## Group law

The curve is y^2 = x^3 + 4, so a = 0. The routines are the EFD formulas, executed on the profile generator:

- Doubling is dbl-2009-l. Leaf count on a non-infinity point: 7 multiplications (2M + 5S), 9 additions, 5 subtractions, plus a zero test of Z that loads and compares six limbs. The counter delta is 6651.
- Complete addition is add-2007-bl. It is not a complete formula in the Edwards sense: it is the general Jacobian sum, and the program does not feed it two equal points. The counter delta is 13568.
- Mixed addition, with the table point affine, is madd-2007-bl. The counter delta is 9756.

The same run checks that these three routines agree with an affine slope computation on the profile generator, and that [r]G is the point at infinity. The executable is this BLS12-381 Jacobian arithmetic. It is not an affine addition on a different curve.

## Window

A scalar is 255 bits. Window w uses t = ceil(255 / w) digits and two phase-1 bases, P0 and A_d = [secret^d]P0. Each base stores t tables of 2^w affine points. Construction is Jacobian: w(t - 1) doubles and t(2^w - 1) additions per base, then one batch inversion. Online, each baby or giant is t - 1 mixed additions. The batch inversion's field work is one field inversion plus (7N - 3) multiplications. The products are 3(N - 1): N - 1 to form the prefix, then two multiplications on each of the N - 1 backward steps. Each point then pays 1S + 3M, and squaring calls multiplication, which is 4N. The sum is 7N - 3. N = 1 has no product and costs four multiplications plus the inversion. Each of those three loops also pays induction on an index of limb_words(N) words. The declared program runs the field schedule on one point and on four points and requires the counter delta to equal the field formula, then requires the extra instructions to equal the induction.

| window | digits | W |
| --- | ---: | ---: |
| 16 | 16 | 37638246693762456819581675044364162482 |
| 24 | 11 | 26716553834528662490785338427987017452 |
| 32 | 8 | 20163538118988385893507537106351234533 |
| 48 | 6 | 15794860975294868162021137239747715846 |
| 64 | 4 | 11432221622494048653151245814711319159 |
| 80 | 4 | 11432221622586297956960252161146333644 |
| 85 | 3 | 9245870455897073182711296040491334261 |
| 96 | 3 | 9245874987990024961937925903081663010 |
| 100 | 3 | 9245943002594557769872498144788988258 |
| 104 | 3 | 9247031236267082696825654012106191419 |
| 112 | 3 | 9543030795193862828084049922385436858 |
| 127 | 3 | 9763443943436312575829073164372886130061 |

The minimum on this list is window 85, with t = 3. 85 * 3 = 255, so the digits cover the scalar with no padding bit. Wider windows with the same digit count share the search list and only grow the tables.

## Width of the search list

points = 447243769829393707157917142500769792. That integer and the byte address of a 128-byte record each need two 64-bit words. A one-word histogram counter, prefix total, or record pointer overflows. Counters and prefix totals are two words. Record addresses are two words. The 256-slot histogram table is 4096 bytes, so the address of one slot stays one word: load the key byte, scale it, and add. Comb entries are stored on the same 128-byte stride, so their address is a shift by 7 as well. The comb has 232113757366008801543585792 entries, and both the index and the byte address are two words.

One record on one radix pass is 328 instructions, and the declared program requires that total.

- Histogram, 33. Two-word address 18, advance the source pointer by one record (6), one-word slot address 3, two-word counter increment 6.
- Scatter, 269. Address 18, advance the source (6), slot 3, destination address 18, then sixteen words of the record. Each word is a load and a store, and both pointers advance, which is 14 instructions per word and 224 for the record.
- Induction of the two record loops, 26. A two-word increment is 6 and a two-word compare-and-branch is 7, so one loop step is 13.
- 33 + 269 + 26 = 328.

The prefix is not inside the 328. It is one two-word addition, 8 instructions, plus one-word induction, on each of 256 buckets: 2816. The 48 passes then pay one-word induction. The sort line is 48 times (points times 328, plus the prefix, plus 3).

## The sum

m = ceil(sqrt(N)) = 223621884914696853578958571250384896, v_max = ceil(N / m) = 223621884914696853578958571250384896, N = (r - 1) / 2^20. Every k0 in {0, ..., N - 1} is u + m v with 0 <= u < m and 0 <= v < v_max. points = m + v_max = 447243769829393707157917142500769792. comb points = 2 * 3 * 2^85 = 232113757366008801543585792.

| Line | Instructions | What it pays for |
| --- | ---: | --- |
| Jacobian construction | 3156514986420353692191225372294 | 340 doubles at 6651 and 232,113,757,366,008,801,543,585,786 additions at 13568, each plus an 18-instruction table address and the loop induction |
| Online mixed additions | 8759269232108675754687807235877576376320 | 894,487,539,658,787,414,315,834,285,001,539,584 mixed additions at 9756. The digit loop and the point loop pay induction, and each point pays an 18-instruction list address |
| Phase 2 | 8554071852 | 2052 scalars. Each bit after the leading one pays a double, a mixed addition, tst, b, and one-word induction |
| Comb normalization | 1271519162850996214855763442490 | One inversion plus (7N - 3) multiplications, plus induction on the three batch loops |
| Search-list normalization | 2450001371125418727811070106619217394490 | The same batch on the baby and giant list. The index is two words |
| Scalar index updates | 181580970550733845106114359855312535552 | One 4-limb Montgomery multiplication (375) plus a two-word address and two-word induction, per baby and per giant |
| Field setup | 873758 | The six exponentiations below |
| Radix sort | 7041405912193974525494247491532119740560 | 48 passes, 328 instructions per record per pass, prefix 2816 per pass |
| Scan | 59483421387309363052002979952602382336 | 12-limb equality at 7 instructions per limb, two two-word addresses, and two-word induction, once per record |
| Root search | 101995520 | 64 candidates, 12 prime tests each, the real Hamming weight, plus one-word induction |
| Trial division | 90194313216 | 3 * 2^32 * 7. The bound fits in one word, and those 7 already include the induction |
| Rebuild and polynomial | 19016 | Multiply the prime powers back, with induction, and check X^4 - X^2 + 1 = r |
| Loads and Montgomery setup | 3192 | Affine input loads, four conversions, and the one-time n0' Newton iteration |
| Total R | 18491740911794146365422592080982668520596 | |

W = ceil(18491740911794146365422592080982668520596 / 2000) = 9245870455897073182711296040491334261.

The field-setup block is the measured exponentiation schedule, not a rounded 510:

- zeta^(2^20): 8709
- beta^m, with the real bit length and Hamming weight of m: 51472
- zeta^N, likewise: 139315
- three inverses in the scalar field, the measured inversion, plus induction on the three calls: 475569
- the unknown discrete log, charged as a 255-bit all-ones exponent: 192645
- 16 multiplications to assemble k0 + k1 N, plus induction: 6048

The sum of those six lines is 873758.

Phase 2 stays binary. A 255-bit all-ones scalar starts from the leading one, then doubles and mixed-adds on each of the other 254 bits. The bit test is tst, b. The 254-trip loop and the loop over the 2052 scalars each pay induction.

## Factorization and the root

The profile publishes parameter_x and r. It does not publish the primes. The attack checks X^4 - X^2 + 1 = r by the scalar-field arithmetic counted above, then factors X, X - 1, and X + 1. Each receives 2^32 trial divisions. (2^32)^2 = 2^64 is at least X + 1. One attempt is udiv, msub, cmp, b.ne, add, cmp, b.lo: 7 instructions, and the index fits in one word. The loop does not stop early in the bound. The product of the prime powers, with the primes of X taken twice, is compared with r - 1. The program recomputes the factorization by a faster method and checks that product; the faster method is not the attack and its cost is not subtracted. The primes it finds are X = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 906349, 254760293], X - 1 = [19, 125527, 2508409, 2529403], X + 1 = [3, 11, 10177, 859267, 52437899]. The root search tests 64 candidates against the 12 distinct primes. Candidate 7 survives. It is not a stored constant. Advice is 0 bytes, so the advice logarithm is 0.

## Memory and preprocessing

Peak storage is one 128-byte comb entry, the batch-inversion product (48 bytes), two 128-byte sort buffers for every search record, two-word histogram counters and two-word prefix totals (8192 bytes), and 65536 bytes of scratch: 114494405117176810328844337551868239872 bytes, log2 ceiling 126.428545. The search-list prefix is one 48-byte field element per denominator. It is live during normalization, before the two sort buffers. A Jacobian point is 144 bytes, so the point plus one prefix element is 192 bytes, and the 256 bytes already charged per record cover that peak. The time score does not add those bytes a second time. Preprocessing is the factorization, the root search, both comb tables, and their normalization: 2214017074635674953568642572 reference units, log2 ceiling 90.838725. Both tables are included even though the A_d table depends on the challenge. That figure is inside W. Phase-1 additions are inside W and are not called preprocessing.

## Bases and coverage

Phase 1 solves secret^d = beta^k0 with beta = zeta^d of order N = (r - 1) / d. The two comb bases are P0 and A_d. A comb on A_1 in phase 1 would solve secret = zeta^(d k0). That is the subgroup of index d, a 1/d fraction of the nonzero scalars, so A_1 is the wrong phase-1 base. Phase 2 solves secret = zeta^k0 * gamma^k1 with gamma = zeta^N of order d, walking A_1 against P0. The rectangles cover every index: every k0 in {0, ..., N - 1} and every k1 in {0, ..., d - 1}. The final scalar multiplication checks the candidate. The identity is the scalar 0, returned before the search. The success probability on a uniform challenge is 1.

The declared program enumerates all 65536 nonzero scalars of the toy field r = 65537 with d = 256. Each scalar is recovered from A_d against P0 and then from A_1 against P0, and zeta^(k0 + k1 N) equals the scalar. The same run counts the scalars fixed by the exponent N and requires that count to equal N, which is exactly the image an A_1 phase-1 comb would hit. It also checks that production beta = 7^(2^20) has order N, using the factored form of r - 1, and that both production rectangles cover.

## What the declared program checks

It remeasures the field and group counters, runs the batch schedule on one point and on four points, rebuilds this sum, and requires W, the score numerator, the two-word widths, and the 328-instruction record. It checks [r]G = infinity with dbl-2009-l and add-2007-bl. It checks the factor product and the coverage argument above. It does not build a table of size 2^85 and it does not execute 2^32 divisions. Progress goes to the diagnostic stream. The only standard output is one JSON object, and the process exits 0.

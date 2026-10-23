# bls12-381-g1-ecdlp: Cheon with a fixed-base comb, revision 8

Analytic bound at the scoring width n = 2^20. It does not compute the production secret. The public source used for the price table is commit 6423c2728dc03d5f2d647c7f22adf1dac0643a66. That commit does not change the reference prices or the definition of W from the tree these bounds were first written against. C = 2000. H = 0: a G1 addition is not a reference operation here, because the reference operation is a 6M affine addition and this schedule counts a different formula. W = ceil(R / C).

The declared score is 123.142791 = 123142791/1000000. It is the ceiling of log2(W) at resolution 10^(-6). W = 11740148960245930949733787335886895195.

## Field schedule

BLS12-381 Fp is a 6-limb Montgomery field, 64-bit limbs, CIOS multiplication. Squaring calls that multiplication. A field doubling, tripling, or multiplication by four or eight is a chain of modular additions, and each of those additions is reduced. The conditional subtraction of the modulus is executed on every add, sub, and product.

The declared program counts instructions by running the routine. One multiplication is 777 instructions. One addition is 86. One subtraction is 85. One inversion by the binary method on p - 2 is 476271 instructions, which is to_mont(1) plus one squaring per bit of p - 2 plus one multiplication per set bit, and four control instructions per bit. A second run of the multiplication on a different pair of inputs produces the same count. The carry identity (2^64 - 1)^2 + (2^64 - 1) + (2^64 - 1) = 2^128 - 1 keeps each mac inside two machine words, so the six-instruction mac (mul, umulh, adds, adc, adds, adc) is not an approximation.

Limb loops inside CIOS and field addition are unrolled. The limb count is the field constant, and the measured schedule already lists every limb instruction. Loops whose trip count depends on the input, including the search list, the scalar bits, and the radix passes, pay induction separately.

## Group law

The curve is y^2 = x^3 + 4, so a = 0. The routines are the EFD formulas, executed on the profile generator:

- Doubling is dbl-2009-l. Leaf count on a non-infinity point: 7 multiplications (2M + 5S), 9 additions, 5 subtractions, plus a zero test of Z that loads and compares six limbs. The counter delta is 6651.
- Addition is add-2007-bl. Before the formula it tests both Z coordinates. It then tests H and R. H = R = 0 calls dbl-2009-l. H = 0 and R nonzero stores a zero Z and returns infinity. The charged delta is the more expensive of those paths and the ordinary finite path: 13594.
- Mixed addition is madd-2007-bl. It tests the Jacobian Z, the affine infinity flag, and both H and R, then takes the same three branches. The charged delta is the maximum of the unequal, equal, negation, and infinity-operand paths: 10056. On this curve the equal-point path is the maximum.

The same run checks each branch on the profile generator: a point plus itself matches a double, a point plus its negative is infinity, infinity on either side returns the other point, and the unequal sum matches an affine slope. It also checks that [r]G is infinity. The executable is this BLS12-381 Jacobian arithmetic.

## Window

A scalar is 255 bits. Window w uses t = ceil(255 / w) digits and two phase-1 bases, P0 and A_d = [secret^d]P0. Ordinary digits satisfy 0 ≤ d_i < B with B = 2^w, so a digit is sometimes zero. Zero is the point at infinity, and the search multiplier 1 has ordinary digits (1, 0, 0, ...). The schedule therefore uses e_i = d_i + 1. Then 1 ≤ e_i ≤ B and

k = sum_i e_i B^i - C, where C = (B^t - 1) / (B - 1).

The group form is [k]P = sum_i [e_i][B^i]P - [C]P. Every scanned window has B = 2^w ≤ 2^127 < r. r is odd, so B^i is invertible modulo r, and r does not divide e_i. The table point [e_i B^i]P is the identity only if r divides e_i B^i, which it does not. No table entry is infinity. The program checks C mod r ≠ 0, so the correction point is finite, and it checks the integer identity on 0, 1, m, a shifted sample, and r - 1. Scalar 1 recodes to (2, 1, 1, ...).

Each base stores t tables of B points, indexed by 1..B. [1] is a copy of the window base, [2] is a double, and [3]..[B] are additions of that base, so the additions are not equal-point sums. Window shifts contribute w(t - 1) doubles per base. Online, the first digit copies a finite table point and each of the other t - 1 digits is a mixed addition, then one mixed addition subtracts the precomputed correction. That is t mixed additions per output point. A partial sum can still be infinity or can equal the next point; those branches are inside the charged mixed addition.

The batch tests every denominator and replaces a zero with Montgomery one before the product, so the inversion is defined. The affine loop writes infinity when the original denominator was zero, and it still performs the 1S+3M that a nonzero point pays. Field work on an all-nonzero batch is one inversion plus (7N - 3) multiplications: 3(N - 1) products and 4N from 1S+3M. The sanitize, prefix, backward, and affine loops each pay induction. The declared program runs this on one point and on four points, once with a real zero denominator, and requires the counter to equal this schedule. The other inverses are checked against x Z^(-2) and y Z^(-3).

| window | digits | W |
| --- | ---: | ---: |
| 16 | 16 | 41025447384565370060742160523093903663 |
| 24 | 11 | 29758258713138369093166332914149410498 |
| 32 | 8 | 22997945510282168512620836999664360691 |
| 48 | 6 | 18491070041711368125622774152452645645 |
| 64 | 4 | 13992244960997498239490204428626404067 |
| 80 | 4 | 13996717398788471851395531864232021273 |
| 85 | 3 | 11740148960245930949733787335886895195 |
| 96 | 3 | 11740153513482473058066207944284769117 |
| 100 | 3 | 11740221845395796735629452248135615460 |
| 104 | 3 | 11741315156008975576641361109749156152 |
| 112 | 3 | 12038695642793620331880571468632209991 |
| 127 | 3 | 9814430471823849067849105960462781225024 |

The minimum on this list is window 85, with t = 3. 85 * 3 = 255, so the digits cover the scalar with no padding bit. Wider windows with the same digit count share the search list and only grow the tables.

## Width of the search list

points = 447243769829393707157917142500769792. That integer and the byte address of a 128-byte record each need two 64-bit words. A one-word histogram counter, prefix total, or record pointer overflows. Counters and prefix totals are two words. Record addresses are two words. The 256-slot histogram table is 4096 bytes, so the address of one slot stays one word: load the key byte, scale it, and add. Comb entries are stored on the same 128-byte stride, so their address is a shift by 7 as well. The comb has 232113757366008801543585792 entries, and both the index and the byte address are two words.

One record on one radix pass is 334 instructions. The declared program requires that total, and it requires the total to be 328 plus one two-word increment.

- Histogram, 33. Two-word address 18, advance the source pointer by one record (6), one-word slot address 3, two-word counter increment 6.
- Scatter, 275. Address 18, advance the source (6), slot 3, destination address 18, then sixteen words of the record. Each word is a load and a store, and both pointers advance, which is 14 instructions per word and 224 for the record. After the copy, the bucket's next-free offset is incremented: a two-word add of a small immediate, 6 instructions. Without that increment every record in the bucket would be written to the same slot.
- Induction of the two record loops, 26. A two-word increment is 6 and a two-word compare-and-branch is 7, so one loop step is 13.
- 33 + 275 + 26 = 334.

The prefix is not inside the 334. It is one two-word addition, 8 instructions, plus one-word induction, on each of 256 buckets: 2816. Each pass also clears the 256 counters first: a one-word slot address and two stores of the zero register, plus one-word induction, 2048 instructions. The 48 passes then pay one-word induction. The sort line is 48 times (points times 334, plus the prefix, plus the clear, plus 3).

The declared program runs this schedule on 21 records, including a repeated key. The counted instructions equal the formula at that one-word width, the output is sorted by the little-endian affine x, and the output is a permutation of the input. The two-word production schedule is the same helpers with limb_words equal to 2.

## The sum

m = ceil(sqrt(N)) = 223621884914696853578958571250384896, v_max = ceil(N / m) = 223621884914696853578958571250384896, N = (r - 1) / 2^20. Every k0 in {0, ..., N - 1} is u + m v with 0 <= u < m and 0 <= v < v_max. points = m + v_max = 447243769829393707157917142500769792. comb points = 2 * 3 * 2^85 = 232113757366008801543585792.

| Line | Instructions | What it pays for |
| --- | ---: | --- |
| Jacobian construction | 3162549944111869921031358616008 | 346 doubles at 6651 and 232,113,757,366,008,801,543,585,780 additions at 13594. Each table also copies its [1] entry, and each base sums its window bases and negates Y. Every group operation pays an 18-instruction table address, or a 9-instruction address for the two correction points, plus induction |
| Online comb | 13589054702496298398286154457743389360128 | 1,341,731,309,488,181,121,473,751,427,502,309,376 mixed additions at 10056. Each digit pays 14 instructions to extract and add one. The first digit is a copy, then t - 1 mixed additions, then the correction mixed addition. The point loop pays induction, a list address, and a store of the Jacobian result |
| Phase 2 | 8710434252 | 2052 scalars. Each bit after the leading one pays a double, a mixed addition at the worst-case price, tst, b, and one-word induction |
| Comb normalization | 1286142329565054769353009347386 | One inversion plus (7N - 3) multiplications, a zero mask on every denominator, and induction on the four batch loops |
| Search-list normalization | 2478177728624670531362018886596765891386 | The same batch on the baby and giant list. The index is two words |
| Scalar index updates | 181580970550733845106114359855312535552 | One 4-limb Montgomery multiplication (375) plus a two-word address and two-word induction, per baby and per giant |
| Field setup | 873758 | The six exponentiations below |
| Radix sort | 7170212117904839913155727628572341538960 | 48 passes, 334 instructions per record per pass, prefix 2816 and clear 2048 per pass |
| Scan | 61272396466626937880634648522605461504 | Infinity flags, then 12-limb equality at 7 instructions per limb, two two-word addresses, and two-word induction, once per record |
| Root search | 101995520 | 64 candidates, 12 prime tests each, the real Hamming weight, plus one-word induction |
| Trial division | 90194313216 | 3 * 2^32 * 7. The bound fits in one word, and those 7 already include the induction |
| Rebuild and polynomial | 19016 | Multiply the prime powers back, with induction, and check X^4 - X^2 + 1 = r |
| Loads, match stores, and setup | 3204 | Affine input loads, four conversions, the one-time n0' Newton iteration, and 12 stores when a match is written |
| Total R | 23480297920491861899467574671773790389890 | |

W = ceil(23480297920491861899467574671773790389890 / 2000) = 11740148960245930949733787335886895195.

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

Peak storage is one 128-byte comb entry, the batch-inversion product (48 bytes), two 128-byte sort buffers for every search record, two-word histogram counters and two-word prefix totals (8192 bytes), and 65536 bytes of scratch: 114494405117176810328844337551868239872 bytes, log2 ceiling 126.428545. The search-list prefix is one 48-byte field element per denominator. It is live during normalization, before the two sort buffers. A Jacobian point is 144 bytes, so the point plus one prefix element is 192 bytes, and the 256 bytes already charged per record cover that peak. The time score does not add those bytes a second time. Preprocessing is the factorization, the root search, both comb tables, and their normalization: 2224346136838462345237332146 reference units, log2 ceiling 90.845440. Both tables are included even though the A_d table depends on the challenge. That figure is inside W. Phase-1 additions are inside W and are not called preprocessing.

## Bases and coverage

Phase 1 solves secret^d = beta^k0 with beta = zeta^d of order N = (r - 1) / d. The two comb bases are P0 and A_d. A comb on A_1 in phase 1 would solve secret = zeta^(d k0). That is the subgroup of index d, a 1/d fraction of the nonzero scalars, so A_1 is the wrong phase-1 base. Phase 2 solves secret = zeta^k0 * gamma^k1 with gamma = zeta^N of order d, walking A_1 against P0. The rectangles cover every index: every k0 in {0, ..., N - 1} and every k1 in {0, ..., d - 1}. The final scalar multiplication checks the candidate. The identity is the scalar 0, returned before the search. The success probability on a uniform challenge is 1.

The declared program enumerates all 65536 nonzero scalars of the toy field r = 65537 with d = 256. Each scalar is recovered from A_d against P0 and then from A_1 against P0, and zeta^(k0 + k1 N) equals the scalar. The same run counts the scalars fixed by the exponent N and requires that count to equal N, which is exactly the image an A_1 phase-1 comb would hit. It also checks that production beta = 7^(2^20) has order N, using the factored form of r - 1, and that both production rectangles cover.

## What the declared program checks

It remeasures the field and group counters, runs the batch schedule on one point and on four points with a real zero denominator, rebuilds this sum, and requires W, the score numerator, the two-word widths, and the 334-instruction record. It sorts 21 records with the counted radix schedule and checks the order and the permutation. It runs a window-4 recoded comb on secp256k1 for scalars that include 0 and 1, and requires the result to match a binary scalar multiplication, including the infinity result of scalar 0. It checks [r]G = infinity with dbl-2009-l and add-2007-bl. It checks the factor product and the coverage argument above. It does not build a table of size 2^85 and it does not execute 2^32 divisions. Progress goes to the diagnostic stream. The only standard output is one JSON object, and the process exits 0.

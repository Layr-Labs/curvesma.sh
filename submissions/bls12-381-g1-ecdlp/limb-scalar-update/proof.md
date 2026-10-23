# Cheon p−1 with a four-limb scalar update

Parent crown: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

The scoring input is n = 2^20. This n divides r − 1. The attack uses A_1 = T = [s]P0 and A_n = [s^n]P0, plus the public generator P0. Other powers are not read. Generation of the powers is organizer work. Every arithmetic step, probe, root test and check below is charged.

C = 2000 is the fixed BLS12-381 G1 price. One reference operation is one batched incomplete affine addition under the 6M convention (S = M, A = 0). Work that is not that operation costs its AArch64 instruction count divided by C. A reference addition is not also unpacked into `blst_fp_mul`.

## Parameters

r is the subgroup order in the target profile. n = 2^20, N = (r − 1)/n, which is an integer. The public base 7 is a primitive root modulo r: 7^{(r−1)/q} ≢ 1 (mod r) for every prime q dividing r − 1. The factorization used for that check is

r − 1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2 · 2508409 · 2529403 · 52437899 · 254760293^2.

η = 7^n mod r then has order N. λ = 0xac45a4010001a40200000000ffffffff satisfies λ^2 + λ + 1 ≡ 0 (mod r) and is the eigenvalue of the public degree-6 endomorphism φ(x, y) = (βx, y), where β is the public cube root of unity in Fp with β^2 + β + 1 = 0. The μ6 orbit of an affine point is the set of its images under φ and negation. Those six points share the y-coordinate up to sign, and y^2 = x^3 + 4 forces x^3, hence x up to multiplication by {1, β, β^2}. So min(y, p − y) is constant on the orbit and separates distinct orbits, except for at most six scalars fixed by a nontrivial μ6 element. Those exceptions are tested directly and do not affect the probability bound below.

G is the integer defined in the cost section. It satisfies 6G^2/N · (1 − 2^{−40}) ≥ 39/100.

## Search

Every s^n lies in the unique subgroup of N-th roots of unity in F_r. Write that value as η^x with 0 ≤ x < N, and write x = x0 + x1 G with 0 ≤ x0, x1 < G, discarding the negligible corner x ≥ G^2 (absorbed by the coverage inequality). Then

[η^{−G x1}] A_n = [η^{x0}] P0

in G1, up to the μ6 action on the stored key. Baby keys are the orbits of [η^j]P0 for j < G. Giant keys are the orbits of [η^{−G i}] A_n for i < G. A stored hit yields at most six candidate scalars; each is checked by one scalar multiplication against A_1.

## Limb schedule

Scalars are four little-endian 64-bit limbs. r < 2^255, so the top limb is at most 63 bits. The fixed-base table for limb position i ∈ {0, 1, 2, 3} and digit d is

T_i[d] = [d · 2^{64 i}] P0,

stored affine. T_i[0] is the identity and is not stored. Each table is an arithmetic progression of affine points with step B_i = [2^{64 i}]P0, built by repeated reference additions. Four tables, each of length at most 2^64, cost at most 2^66 reference additions. The 192 doublings that produce the four steps B_i are inside the 2^28 term below.

The scalar update that advances a baby or giant index is one multiplication in F_r by the fixed element η or η^{−G}. That is a 4 × 4 Montgomery multiplication, not a 6-limb `blst_fp_mul`. The pinned base-field kernel is 383 instructions and is the wrong width: F_r elements are four limbs. The schedule below is unrolled AArch64 (mul, umulh, adds, adcs, adc, ldr, str, csel). Operands of one product stay in general-purpose registers.

Per outer limb i = 0..3:

- inner product, four limbs, six instructions each (mul, umulh, adds, adcs, two carries): 24 instructions, times 4 = 96;
- Montgomery fold: one mul by n0, then four limbs at five instructions (mul, umulh, adds, adcs, adc): 21 instructions, times 4 = 84;
- carry retire: 4.

Conditional subtract of r: four subs/sbcs and four csel = 8. Reload and store the running scalar: four ldr and four str = 8. Constant setup across the unrolled product: 8. Total 96 + 84 + 4 + 8 + 8 + 8 = 208 instructions. This is the scalar-update multiply. It is strictly below the pinned 383-instruction base-field multiply because the operand width is four limbs and the second schoolbook width is never issued.

The four limbs of the updated scalar are the digits. The point [scalar]P0 is the sum of the four affine table points T_i[digit_i], skipping a zero digit. Worst case there are four nonzero limbs and therefore three affine additions. Zero digits only drop work; the bound charges three additions on every point.

## Why each of those additions is one reference operation

Across one side of the rectangle the G additions in a fixed layer (second table point, then third, then fourth) are independent affine additions. The denominator x2 − x1 of each slope is inverted by one Montgomery batch over that layer. That is the organizer’s batched incomplete affine addition: three multiplications for the chord-and-tangent arithmetic and three for the amortized inversion, S = M, A = 0, price 1. Incomplete addition is undefined when x1 = x2. For these points that event has probability at most 2^{−140} per layer after multiplying by G/p, which is inside the 2^{−40} failure budget already subtracted in the coverage inequality. A hit of that case aborts the point and is a failure, not a retry loop.

The same three-layer batch is used on the giant side. Both sides also pay the 208-instruction scalar update.

## Dictionary

The key is min(y, p − y), with y the affine coordinate produced by the last batched addition.

- Load six limbs, subtract from p with sbcs, csel the minimum, store it: 6 + 6 + 6 + 6 + 8 address = 32 instructions.
- Address hash, two 64-bit accumulators, six limbs, each limb a mul, umulh, add and adc into each accumulator: 8 · 6 = 48 instructions. The low 120 bits select the bucket.
- The table has B = 2^120 buckets and 21 slots in each bucket. A slot holds the six-limb key and the 16-byte index. A probe is three ldp, six eor, five orr, one cbz and one pointer add = 16 instructions. Every lookup and insertion scans all 21 slots: 336 instructions. There is no early-out credit.

Dictionary instruction total: 32 + 48 + 336 = 416. Scalar update: 208. Combined non-reference work per point: R = 624 instructions = 624/2000 reference units.

Bucket load is λ = G/B < 2^{−4.5}. The expected number of buckets with 22 or more keys is below 2^{−40} under the Poisson tail λ^k/k! · B evaluated at k = 22 (the 21 slots overflow). That expectation is the δ = 2^{−40} subtracted in the success product. On overflow the point is a failure. No second hash and no variable-length chain is used, so the instruction count cannot drift from the scan above.

## Recovery

A hit gives x up to at most six μ6 lifts. One d-th root of η^x in F_r, then multiplication by the d-th roots of unity, is checked against A_1. Charging every check as a 256-addition binary scalar multiplication over at most 2^20 roots and six lifts is at most 2^28 reference additions. Verification of the returned scalar is inside that term. Field construction of η and η^{−G} is two exponentiations, also inside 2^28.

## Work

Per point, on each side: 3 reference additions + 624/C. Point count: 2G. Preprocessing of the limb tables: 2^66. Root tests and setup: 2^28.

W = 2G · (3 + 624/2000) + 2^66 + 2^28.

G = 57012617742347889232813606637047839. The integer check

6 · G^2 · (2^{40} − 1) · 100 ≥ N · 2^{40} · 39

holds, so the coverage product is at least 39/100. log2(W) < 118.184540, and the submitted score is that ceiling at 10^{−6}. The same W is the function evaluated at n = 2^20. The claim is not made for other n.

Peak storage is the bucket array, 21 · 2^{120} · 64 bytes, plus the limb tables, 4 · 2^{64} · 48 bytes. That is at most 2^{131} bytes. The limb tables are parameter-dependent advice; their construction is the 2^{66} term, and their size is at most 2^{72} bytes. Preprocessing is at most 2^{67} reference units and is inside W.

No quantum step, no uncharged oracle, and no early-abort expectation. The only probabilistic statements are the two explicit failure sources subtracted as 2^{−40}: bucket overflow and the incomplete-addition exceptional set.

# bls12-381-g1-ecdlp: Cheon with a fixed-base comb

This is an analytic accounting argument at the scoring width n = 2^20. It does not compute the production secret. The auxiliary powers are a public input. Using them is the ordinary statement of this track, not an extra leak.

The score is log2(W(2^20)) under classical-group-work-v2, with C = 2,000. Applicability: n = 2^20, n divides r - 1, and the public powers A_i = [s^i]P0 for i = 1..n are supplied. The bound is specialized to that n. H(n) = 0, because every elliptic operation is charged as instructions rather than identified with the 6-multiplication reference unit. W(n) = ceil(R(n) / C).

## Algebra

Let r be the BLS12-381 scalar order

r = X^4 - X^2 + 1, X = 0xd201000000010000.

Then r - 1 = X^2 (X - 1)(X + 1). The positive integer X is divisible by 2^16, so r - 1 is divisible by 2^32 and in particular by d = 2^20. The attack uses A_1 = [s]P0 and A_d = [s^d]P0. Generation of the powers is excluded. Every subsequent read is charged.

Prime factorization of r - 1, with multiplicity:

| Prime | Multiplicity |
| --- | ---: |
| 2 | 32 |
| 3 | 1 |
| 11 | 1 |
| 19 | 1 |
| 10177 | 1 |
| 125527 | 1 |
| 859267 | 1 |
| 906349 | 2 |
| 2508409 | 1 |
| 2529403 | 1 |
| 52437899 | 1 |
| 254760293 | 2 |

The product 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 * 254760293^2 equals r - 1. The distinct primes are enough for the primitive-root test. The element 7 is a primitive root modulo r: 7^(r-1) = 1, and 7^((r-1)/q) != 1 for every distinct prime q above. Finding it takes 6 candidate tests; the budget charges 64.

Write s = zeta^k with zeta = 7, and split k = k0 + k1 * N with N = (r - 1) / d, 0 <= k0 < N and 0 <= k1 < d. Then s^d = beta^(k0) for beta = zeta^d. The order of beta is N. The order of gamma = zeta^N is d.

## Phase 1

Set m = ceil(sqrt(N)) = 223,621,884,914,696,853,578,958,571,250,384,896 and v_max = ceil(N / m) = 223,621,884,914,696,853,578,958,571,250,384,896. Every k0 in {0, ..., N-1} is k0 = u + m*v with 0 <= u < m and 0 <= v < v_max, because v_max = ceil(N / m) is large enough for the last residue. Search

[beta^(-u)] A_d = [(beta^m)^v] P0.

A match yields k0. The loops cover the whole rectangle, so a match exists for every s. This phase is deterministic.

## Phase 2

Search [zeta^(-k0)] A_1 = [gamma^(k1)] P0 over a baby-step/giant-step rectangle of side 2^10, since d = 2^20 is a perfect square. Then s = zeta^(k0 + k1 * N). The same coverage argument applies. One verification scalar multiplication checks the candidate against A_1.

## Fixed-base comb

Cheon's original step is an exponentiation. Charging each step as a full scalar multiplication would put a factor near the bit length of r in front of sqrt(N). That factor is not claimed away by citation. It is removed here by a fixed-base comb whose table is counted in full.

For a 255-bit scalar and window w = 96, the digit count is t = ceil(255 / 96) = 3. Write the scalar as e0 + e1 * 2^96 + e2 * 2^192 with 0 <= ei < 2^96. For a fixed base B:

- Build B_j = [2^(j*96)] B for j = 1, 2 by 96 doublings each: 192 doublings.
- For each of the 3 digits, build [i] B_j for i = 1..2^96 - 1 by repeated addition: at most 2^96 - 1 additions per digit.

Two bases are used online, P0 and A_d. Doublings: 2 * 96 * 2 = 384. Additions: 2 * 3 * (2^96 - 1) = 475,368,975,085,586,025,561,263,702,010. Each online scalar then costs t - 1 = 2 complete Jacobian additions, whether or not a digit is zero (skipping a zero digit would only be cheaper). Baby and giant indices of phase 1 together are m + v_max, so the online addition count is

(m + v_max) * 2 = 894,487,539,658,787,414,315,834,285,001,539,584.

Phase 2 is not combed. It is charged as 2 * 2^10 + 4 binary scalar multiplications (two rectangle sides, plus setup and verification), each at 254 doubles + 254 additions = 508 group operations:

(2048 + 4) * 508 = 1,042,416.

That is negligible beside the phase-1 additions and is an upper bound. Total group operations:

384 + 475,368,975,085,586,025,561,263,702,010 + 894,487,539,658,787,414,315,834,285,001,539,584 + 1,042,416 = 894,488,015,027,762,499,901,859,846,266,284,394.

Window 96 is the minimum inside the checked set {16, 24, 32, 48, 64, 80, 96, 100, 104}. Once the window is at least 85, t stays 3, so a wider window only grows the table.

## Instruction total

One group operation is charged at (16, 24) base-field operations: 16 * 512 + 24 * 128 = 11,264 instructions. The pair dominates doubling (8, 16), mixed addition (11, 15) and complete addition (15, 12) on this a = 0 curve. Field multiplication at 512 instructions is above the published measurement of 383.

| Component | Instructions | Reason |
| --- | ---: | --- |
| Group operations | 10,075,513,001,272,716,798,894,549,308,343,427,414,016 | 894,488,015,027,762,499,901,859,846,266,284,394 * 11,264 |
| Scalar formation | 915,955,240,610,598,312,259,414,307,841,580,728,320 | 447,243,769,829,393,707,157,917,142,500,771,840 modular multiplications at 2,048 |
| Affine conversion | 1,831,910,481,221,196,624,518,828,615,683,153,458,176 | 7 multiplications per stored point, plus one Fermat inversion of 2 * 381 multiplications, and 4 additions per point |
| Radix sort | 5,495,731,443,663,589,873,556,485,847,049,509,535,744 | 447,243,769,829,393,707,157,917,142,500,769,792 records, 24 passes, 512 instructions per record per pass, plus 65,536 * 32 bucket instructions per pass |
| Linear scan | 457,977,620,305,299,156,129,707,153,920,788,267,008 | 1,024 instructions per record |
| Primitive-root search | 66,846,720 | 64 candidates, 510 modular multiplications each, 2,048 instructions each |
| Reserve | 1,099,511,627,776 | control and address arithmetic |
| Sum R | 18,777,087,787,073,400,765,358,985,233,938,037,877,760 | |

Stored points are the phase-1 babies and giants: m + v_max = 447,243,769,829,393,707,157,917,142,500,769,792. They are converted to affine coordinates with one Montgomery batch inversion so the sort key is the affine x. The sort is one radix sort of the combined list on a 48-byte key, 16 bits at a time, 24 passes. A pass at 512 instructions is several times a count of digit extraction, a bucket increment, and a scatter.

W = ceil(R / 2000) = 9,388,543,893,536,700,382,679,492,616,969,018,939.

The declared score is 122820313/1000000 = 122.820313, the 10^(-6) ceiling of log2(W).

## Success, memory, preprocessing, advice

The algorithm is deterministic. Every k0 and every k1 lies in the searched rectangle, 7 generates the multiplicative group, and the final scalar multiplication rejects a mismatched candidate. The success probability on uniform s is 1, which is at least 0.39. There is no heuristic about the group walk. The instruction ceilings are explicit upper bounds, not a model of branch behavior.

Peak memory is 114,494,450,711,746,397,248,685,242,361,513,508,864 bytes: two radix-sort buffers of the point list at 128 bytes per record, plus both comb tables at 96 bytes per table point, plus 2^20 bytes of scratch. log2 of that quantity, rounded up, is 126.428546. This is an accounting figure for the table the algorithm writes. It is not a claim that this much storage is available, and it is not charged as time. Memory does not change the score. The same style of bound appears in a baby-step table that stores one point per residue class.

Preprocessing that can be performed from public data alone is the primitive-root search and the comb tables. Counting both comb tables overstates that figure, because the A_d table depends on the challenge. The resulting ceiling is 2,677,278,067,682,020,495,961,037,742,178 reference units, log2 rounded up to 101.078610. It is included in W. The phase-1 additions are inside W whether one calls the P0 giant list preprocessing or online. There is no amortization across challenges and no free oracle.

Nonuniform advice is 0 bytes. The factorization and the root 7 are public constants of r, and the root search is charged anyway.

## Experiment

The declared program factors r - 1, checks that 7 is a primitive root, recomputes m, v_max, the group-operation count and W, checks the score inequality with an 80-digit logarithm, recovers scalars with Cheon's index arithmetic in several prime-order additive groups, and checks a width-8 comb against binary scalar multiplication on secp256k1 for a 64-bit scalar. None of that is a production elliptic-curve computation.

## What is not claimed

The Kozaki-Kutsuma-Matsuo improvement is not cited as a black box. The comb above is the method that removes the per-step exponentiation, and its table is counted. GLV, the negation map, and a number-field sieve are not used. No production secret is recovered.

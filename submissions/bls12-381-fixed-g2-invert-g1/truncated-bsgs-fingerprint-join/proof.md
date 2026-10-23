# Truncated baby-step giant-step with a randomized 256-bit fingerprint join

## Claim

For the pinned BLS12-381 parameters, every nonidentity `Q` in G2 and every
`z` in GT with `z^r = 1`, the attack below returns `P` in G1 with
`e(P,Q) = z`. Over the uniform production challenge and fresh algorithm coins
it succeeds with probability at least `0.4 - 2^-38`, using at most

```text
log2 W          = 127.878612        (rational 127878612/1000000)
peak storage    = 2^137.089509 bytes
preprocessing   = 2^126.878612 target-group-operation equivalents (inside W)
nonuniform advice = 2^0 bytes
```

in the `classical-group-work-v2` accounting with the fixed track price
`C = 20000` scalar AArch64 instruction equivalents per dense Fp12
multiplication. The attack is a classical cold attack on one target. It uses
no auxiliary powers, no quantum step, no uncharged advice, no random-mapping
heuristic, no distinguished points, no cycle detection and no restart on
success-relevant randomness. All preprocessing, all trials, all sorting,
lookup, memory transfers, verification and output are charged once, in full.

## 1. Reduction to a discrete logarithm in GT

Let `P0` be the pinned G1 generator and `Q` the supplied nonidentity G2 point.
The pinned pairing is bilinear and nondegenerate, so `g = e(P0,Q)` is not the
identity, and GT has prime order `r`, so `g` generates GT. The challenge
samples `b` uniformly in `[1,r-1]`, sets `Q = [b]Q0`, samples `a` uniformly in
`[0,r)`, and sets `z = e(P0,Q0)^(a b)`. Since `e(P0,Q) = e(P0,Q0)^b`, the
instance is exactly `z = g^a` with the same `a`. Conversely, for the unique
`x` in `[0,r)` with `z = g^x`, returning `[x]P0` gives
`e([x]P0,Q) = e(P0,Q)^x = g^x = z`. The rules state that this reduction is
admissible, and the pairing convention (`pairing(Q,P)` from the pinned library,
G2 argument first) is the one fixed by `target-profile.json`.

For every nonidentity `Q`, the same `a` is the exponent of `z` with respect to
the instance generator `g = e(P0,Q)`; the value `b` never has to be known.

## 2. Parameters

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
log2 r = 254.8570894130472
m = 144825239754852386682327417126749682654
m < 2^127,  m^2 < r,  10 m^2 - 4 r = 2624165268798259167710680089142760099108 > 0
```

`m` is the least integer with `10 m^2 >= 4 r`, i.e. `m^2 >= 2r/5 = 0.4 r`.
The inequality is verified in exact integer arithmetic above; consequently
`m^2 / r > 0.4` and `m^2 < r`. The exponent grid therefore covers
`m^2 = 0.4 r (1 + 2^-127.6)` distinct residues of `[0,r)`, and no covered
exponent wraps modulo `r`.

## 3. Representations (fixed-width, no implicit pointers)

Every quantity that can reach `m` or `m^2` is carried in explicit 64-bit
limbs. This is the part of the design that the earlier accounting did not
support, so it is stated first and completely.

| object | range | encoding | bytes |
| --- | --- | --- | ---: |
| GT element | 12 Fp coefficients | 12 x 48 B tower/polynomial coefficients, the direct output words of the pinned Fp12 multiplication (6 Fp2 = 12 Fp = 72 little-endian 64-bit words) | 576 |
| record reference | `[0, 2m)` | 2 x 64-bit limbs (128 bits); `2m < 2^128` | 16 |
| fingerprint | `[0, 2^256)` | 4 x 64-bit lanes | 32 |
| radix counter | `[0, m]` | 2 x 64-bit limbs (128 bits); `m < 2^127` | 16 |
| radix prefix sum | `[0, m]` | 2 x 64-bit limbs | 16 |
| bucket digit | `[0, 2^64)` | one 64-bit limb | 8 |

Consequences that are used in the ledger:

* A record reference is 128 bits, so it distinguishes all `2m < 2^128`
  records; no 64-bit pointer is used anywhere.
* Each radix counter is 128 bits. In every pass the counters are distributed
  over `2^64` buckets, so the largest counter is `ceil(2m/2^64) < 2^64`, which
  fits with 63 bits to spare; the 128-bit width is kept for uniformity and is
  charged as two limbs.
* Prefix sums and destination positions are 128-bit values. A destination
  address is `base + ref * 16` (references) or `base + ref * 32`
  (fingerprints); the product of a 128-bit reference by a 64-bit constant is
  formed with four 64-bit multiplications and three carry additions, then one
  base addition with carry. This is charged explicitly below.
* The fingerprint is computed from the 72 words that the pinned Fp12
  multiplication already produces. No basis conversion, canonicalisation or
  serialisation is performed or charged, and none is needed: equal GT elements
  have equal tower words, which is all the join requires.

## 4. Algorithm

```text
INPUT : Q in G2 (nonidentity), z in GT
OUTPUT: P in G1 with e(P,Q) = z, or failure

 1. Validate. Require Q on the G2 curve, Q != infinity, [r]Q = infinity,
    z != 0 and z^r = 1. Reject otherwise.
 2. g  <- e(P0, Q)                          (pinned pairing, Q first)
 3. Seed a PRNG with a fresh uniform 256-bit seed and draw the mask words
    A[l][k], l in {0,1,2,3}, k in [0,72), each uniform in [0,2^64).
 4. Build the baby list in E[0..m-1] and K[0..m-1]:
       E[0] <- 1;  K[0] <- fp(E[0])
       for j = 1 .. m-1:
           E[j] <- E[j-1] * g;  K[j] <- fp(E[j])
       Gm <- E[m-1] * g                     (this is g^m)
 5. Build the giant list in E[m..2m-1] and K[m..2m-1]:
       D  <- Gm^(r-1)                       (= g^(-m))
       E[m] <- z;  K[m] <- fp(E[m])
       for i = 1 .. m-1:
           E[m+i] <- E[m+i-1] * D;  K[m+i] <- fp(E[m+i])
 6. X[0..2m-1] <- identity permutation; Y[0..2m-1] <- scratch references.
    LSD radix sort X by the 256-bit key K, four stable passes, pass l using
    the 64-bit digit d = lane l of K, with 128-bit counters, 128-bit prefix
    sums and 128-bit destinations; X and Y are ping-ponged.
 7. Merge the two sorted halves of X by 256-bit key. For every pair
    (i, j) whose keys are equal, compare E[m+i] and E[j] exactly, word by
    word over all 72 words. On the first exact equality set
    a <- i*m + j and go to 8. If the merge finishes, output failure.
 8. P <- [a]P0 by left-to-right double-and-add over the 254 bits of a.
 9. Verify g^a = z by square-and-multiply over the 254 bits of a, and
    [r]P = infinity. Output P.
```

`fp(x)` is the 256-bit fingerprint of the 72 words `x_0..x_71` of `x`:

```text
fp(x)_l = sum_{k=0}^{71} A[l][k] * x_k  mod 2^64      l = 0,1,2,3
```

Each lane is a linear map over `Z/2^64`; the four lanes use independent mask
words. The masks are drawn once per attack from a fresh seed, so they are
uniform and independent of the challenge.

## 5. Correctness

Fix a valid instance and let `a` in `[0,r)` satisfy `z = g^a`.

* If `a < m^2`, write `a = i m + j` with `0 <= i,j < m`. Then
  `E[m+i] = z * g^(-i m) = g^(a - i m) = g^j = E[j]`, so the two stored
  elements are identical, their fingerprints are identical, the two sorted
  lists place them in the same key group, and step 7 reaches the pair and
  confirms it by exact comparison. It sets `a = i m + j`, the correct value.
* Conversely, any exact match in step 7 gives `g^(i m + j) = z`, with
  `i m + j < m^2 < r`. Since the discrete logarithm of `z` to base `g` is
  unique modulo `r`, `i m + j = a`. A wrong answer is impossible; a
  fingerprint collision only causes an extra exact comparison, never a false
  match, because step 7 compares all 72 words.
* `a = 0` (so `z = 1`) matches `(i,j) = (0,0)` because `E[m] = z = 1 = E[0]`;
  step 8 returns the point at infinity, which the acceptance rule allows
  exactly for `z = 1`.
* If `a >= m^2` no pair exists and step 7 exhausts the merge and reports
  failure; it never emits a point.

The output `P = [a]P0` lies on the G1 curve and in the order-`r` subgroup by
construction; step 9 re-checks `g^a = z` and `[r]P = infinity` before the
point is emitted.

## 6. Success probability

The challenge exponent `a` is uniform on `[0,r)` and independent of `b`; the
attack succeeds whenever `a < m^2`, which has probability

```text
m^2 / r = 0.4 (1 + 1.25e-38) > 2/5.
```

The only other failure source is the cap on exact comparisons in step 7. For
a fixed challenge and a fixed pair `(i,j)` with `i` in the giant half and `j`
in the baby half, write `d_k = (E[m+i])_k - (E[j])_k` over the 72 words. If
some `d_k` is odd then, because `A[l][k]` is uniform on `[0,2^64)`, the lane
difference `sum_k d_k A[l][k] mod 2^64` is uniform, so the pair collides with
probability `2^-64` per lane and `2^-256` over the four independent lanes. If
all `d_k` are even, each lane is uniform on the even residues at worst, and
the collision probability is at most `2^-63` per lane, hence `2^-252`
jointly. Either way the expected number of colliding cross pairs is at most
`m^2 2^-252 = 2^1.53`, and the attack charges `2^40` exact comparisons. By
Markov, the probability that more than `2^40` comparisons occur is at most
`2^1.53 / 2^40 < 2^-38`, in which case the attack reports failure. Therefore

```text
P(success) >= m^2/r - 2^-38 >= 0.4 - 2^-38 > 0.3999999
```

over the uniform challenge and fresh algorithm coins, which exceeds the
required `0.39`.

## 7. Work ledger

`H` counts reference operations (one dense Fp12 multiplication, the track's
reference unit at `C = 20000`). `R` counts every other scalar AArch64
instruction equivalent. Total work is `W = H + R/C`.

### 7.1 Reference operations

| phase | reference operations |
| --- | ---: |
| baby chain `E[1..m-1]` | `m-1` |
| `Gm = g^m` | 1 |
| `D = Gm^(r-1)` by square-and-multiply (254 squarings, <=254 multiplications) | <= 508 |
| giant chain `E[m+1..2m-1]` | `m-1` |
| pairing `g = e(P0,Q)` | <= 2^16 |
| `P = [a]P0` (254 doublings, <=254 additions, each charged as one Fp12 multiplication for safety) | <= 2^10 |
| verification `g^a = z`, `[r]P = infinity`, input checks `z^r = 1` | <= 2^10 |

```text
H <= 2m + 2^17
```

A squaring is charged as a multiplication throughout; a G1 addition costs far
less than one Fp12 multiplication (`2,298` instructions against `20,000`) and
is over-charged here. The pairing bound `2^16` is more than an order of
magnitude above a Miller loop of 64 iterations plus the final exponentiation
counted in Fp12 multiplications.

### 7.2 Other instructions

| phase | per element | elements | instructions |
| --- | ---: | ---: | ---: |
| fingerprint (72 words: 1 load + 4 MUL + 4 ADD + 2 loop/address each, plus 64 for the lane combine, rounded up) | 1024 | `2m` | `2048 m` |
| radix sort, 4 passes (128 per record per pass: 2 key words and 2 reference words loaded/stored, 128-bit address formation by four 64-bit products and carry additions, counter read/increment/store, digit extraction, loop control) | 512 | `2m` | `1024 m` |
| merge over both sorted halves (key compare, pointer advance, run bookkeeping) | 64 | `2m` | `128 m` |
| counter zeroing, 4 passes over `2^64` two-limb counters | - | - | `2^69` |
| prefix sums, 4 passes over `2^64` two-limb counters (load, add, store, carry) | - | - | `2^69` |
| exact comparisons on fingerprint matches, capped | 256 | `2^40` | `2^48` |
| mask generation, input decoding, output encoding, control | - | - | `2^30` |

```text
R <= 3200 m + 2^69 + 2^69 + 2^48 + 2^30
```

`3200 m` instructions per element pair is the sum of the three per-record
rows; the counter and comparison terms are at most `2^70` and are negligible
against `3200 m ~ 2^138.4`.

### 7.3 Total

```text
W = H + R/C <= (2m + 2^17) + (3200 m + 2^70)/20000
  = 2m (1 + 0.08) + 2^17 + 2^56
```

Numerically, with `2m = 289650479509704773364654834253499365308`:

```text
H  = 289650479509704773364654834253499433399
R  = 463440767215527637384038030897433740596992
W  = H + R/20000 = 312822517870481155233856735798371120428
log2 W = 127.87861197146867  ->  claimed 127.878612
```

The claimed score `127878612/1000000` is `log2 W` rounded up at the sixth
decimal, so the submitted rational is an upper bound on the computed bound.

Preprocessing (the whole baby half: `m` multiplications, `m` fingerprints,
one sort of `m` records, the counter work) is
`m + (1600 m + 2^69)/20000`, i.e. `2^126.878612` equivalents, and is included
in `W` without amortisation.

## 8. Peak storage

| object | bytes |
| --- | ---: |
| GT elements, `2m` x 576 | `1152 m` |
| fingerprints, `2m` x 32 | `64 m` |
| record references, two arrays of `2m` x 16 | `64 m` |
| radix counters, `2^64` x 16 | `2^68` |
| code, masks, scratch, I/O buffers | `2^20` |

```text
peak <= 1280 m + 2^68 + 2^20 = 2^137.089509 bytes
```

The element, key and reference arrays are disjoint; the two reference arrays
are the ping-pong buffers of the radix sort. The counter array is reused by
all four passes. The storage is astronomical and mandatory for this
construction; it is reported as a bound and carries no scalar contribution.

## 9. Randomness, heuristics and limitations

* The only randomness used by the attack is the 256-bit mask seed, which is
  charged and drawn fresh per instance. The masks make the fingerprint a
  random linear map, so the collision bound of Section 6 holds for every
  fixed pair of elements, whatever the challenge is; the analysis does not
  assume that GT elements look random, that the walk is a random mapping, or
  that any iteration is independent.
* Fingerprint collisions cannot cause a wrong output: every candidate pair is
  confirmed by comparing all 72 words of the two stored elements, and the
  recovered exponent is additionally verified by `g^a = z` before output.
* No heuristic constant (birthday constant, random-mapping constant,
  distinguished-point rate, cycle-length model) appears anywhere in the
  success or cost argument. The success probability is the exact count of
  covered exponents, `m^2/r`, minus one Markov term `2^-38`.
* The cost bound is analytic and deterministic given the algorithm: every
  instruction count in Section 7 is an explicit upper bound per operation,
  and the only non-deterministic term is the `2^40` comparison cap, which is
  charged unconditionally.
* Limitations. The attack is analytic; no production-size execution is
  claimed, and the `2^137`-byte storage is not physically realisable, which
  the rules permit as a supplementary mathematical bound. The construction is
  a generic square-root attack with an exact join; it makes no claim about
  best-known security levels and does not exploit any structure of the
  pairing beyond the admissible reduction of Section 1.

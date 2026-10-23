# Exact baby-step giant-step for BLS12-381 G1 (no heuristics)

Claim: deterministic baby-step giant-step with negation-invariant keys recovers x from T=[x]P0 for exactly
g(2m+1)+1 of the r scalars, a success probability of at least 2/5 = 0.40 over uniform x, using at most
W = 2^127.851688 target-group-operation equivalents (C = 2000). The submitted score is log2 W rounded up to two
decimals: 127.86 = 6393/50. Peak memory is 2^133.8516 bytes.

Every number below is printed by `experiment.py` (function `production()`, exact integer and Fraction
arithmetic). The same file runs this algorithm on three toy prime-order curves and checks, for every scalar of
each toy group, that the output is correct and that the success set equals the formula set.

## 0. What this attack does not use

No random walk, no hash function, no distinguished points, no scalar accumulators, no probabilistic or
distributional assumption about coordinates, keys, collisions or cycles. The algorithm uses no random coins.
Its success set is an explicit set of scalars, counted exactly. Its cost is the full worst-case run: every
baby and giant step, both complete sorts and the longest possible merge, with no credit for stopping early.

## 1. Notation and parameters

r = 0x73eda753...00000001 (prime, log2 r = 254.857), p = the 381-bit base-field prime, G = P0, Q = T = [x]G.
x(P) is the affine x-coordinate. For P != O, x(P) = x(P') holds exactly when P' = P or P' = -P.

| symbol | value |
| --- | --- |
| K (chains per list) | 2^20 = 1048576 |
| Lb (baby chain length) | 97858491008506274031834200607979 |
| Lg (giant chain length) | 97467578261468299568324876238036 |
| m = K Lb (baby scalars 1..m) | 102612065067735474799204578736712187904 (2^126.2705) |
| g = K Lg (giant count) | 102202163343097383688155825426174836736 (2^126.2647) |
| s = 2m+1 (giant stride) | 205224130135470949598409157473424375809 |
| g s | 20974350070050476191779096203274402356170092211728582208516385044243282919424 |

g s exceeds floor(2r/5) by 16021093871211517527079474921564267850445619, so g s / r >= 2/5.
Also g s + m < r - s and 2m < r (both checked by the script).

## 2. Algorithm

Key of an affine point A: the 6-limb canonical (fully reduced, value < p) representation of x(A), read as a
384-bit unsigned integer. A record is (key, index): 6 key words and a 2-word (128-bit) index, 64 bytes.
Field elements are kept fully reduced: x3 = lambda^2 - x1 - x2 is the output of a modular subtraction, which
lies in [0, p) when its inputs do, and the record step applies one more conditional subtraction of p as a
safeguard (charged). The canonical representation (plain or Montgomery form, the same for both lists) is a
bijection of Fp onto [0, p), so equal keys hold exactly when the x-coordinates are equal. The exceptional test
x(A) == x(S) compares the same canonical words.

```
Input Q in G1 (Q may be the identity O).
Constants: -G; S = [s]G and -S (textbook double-and-add, affine group law with all cases).

Baby phase (independent of Q)
 1. For c = 0..K-1: B_c = [(c+1) Lb] G  (double-and-add); emit record (key(B_c), (c+1) Lb).
 2. For t = 1..Lb-1:  B_c <- B_c + (-G) for all c, as ONE batched incomplete affine addition round
    (Montgomery trick over the K chains, one inversion); emit (key(B_c), (c+1) Lb - t) for every c.
    Chain c thus visits j = (c+1)Lb, (c+1)Lb-1, ..., c Lb+1; the K chains cover j = 1..m once each.
 3. LSD radix sort of the m baby records on the full 384-bit key (section 4).

Giant phase
 4. For c = 0..K-1: A_c = Q + [r - c Lg s] G (double-and-add for the offset, then one textbook affine
    group-law operation with all cases). If A_c = O: x = c Lg s; verify and output.
    Emit record (key(A_c), c Lg).  (Chain c will hold A_i = Q - [i s]G for i = c Lg .. (c+1) Lg - 1.)
 5. For t = 1..Lg-1:
      for every c, with i = c Lg + t - 1:  if x(A_c) == x(S):   (exceptional-case test)
            candidates (i+1) s and (i-1) s mod r; output the one with [cand]G == Q.
      A_c <- A_c + (-S) for all c, as ONE batched incomplete affine addition round;
      emit (key(A_c), c Lg + t) for every c.
 6. LSD radix sort of the g giant records.

Match
 7. Merge the two sorted lists with two pointers. At the first pair of equal keys, baby index j and giant
    index i: candidates i s + j and i s - j mod r; output the one with [cand]G == Q.
 8. If the merge ends without equal keys, output FAIL.
```

The run stops at the first verified output; the cost bound below charges the complete run regardless.

## 3. Correctness and exact success set

**Lemma 1 (incomplete addition).** For affine A, B != O with x(A) != x(B), the chord formula
lambda = (yB - yA)/(xB - xA), x3 = lambda^2 - xA - xB, y3 = lambda (xA - x3) - yA gives A + B, and
x(A) == x(B) exactly when A = +-B. The Montgomery batch is valid exactly when every denominator is nonzero.

**Lemma 2 (baby chains are never exceptional).** Every baby addition is [j]G + (-G) with 2 <= j <= m. It is
exceptional only if [j]G = +-G, i.e. j = +-1 mod r, or [j]G = O, i.e. j = 0 mod r. Since 2 <= j <= m < r - 1
neither holds. So no check is needed and none is charged. The keys x([j]G), j = 1..m, are pairwise distinct
because j != +-j' mod r for distinct j, j' in 1..m (2m < r).

**Lemma 3 (giant chains).** Each giant addition is A + (-S) with A = Q - [i s]G. It is exceptional only if A = O,
A = S (then A + (-S) = O, so x = (i+1) s) or A = -S (a doubling, x = (i-1) s). A = O never reaches step 5: a
chain start equal to O is caught in step 4, and A_{i} = O for a later i requires A_{i-1} = S, caught by the
test before the addition that would produce it. The test x(A) == x(S) is exact (Lemma 1), so step 5 never
performs an exceptional addition, and both candidates are checked by scalar multiplication.

**Lemma 4 (sort and merge).** A stable counting sort per 32-bit digit, applied from digit 0 (least
significant) to digit 11, orders records by the full 384-bit key (standard LSD argument by induction on the
digit). The two-pointer merge advances the pointer at the smaller key; it reports equal keys whenever the
lists share a key, and takes at most m + g - 1 steps.

**Theorem (exact success).** Let Cov = { i s + e mod r : 0 <= i <= g-1, -m <= e <= m }. The integers
i s + e range over the interval [-m, (g-1) s + m] without repetition (consecutive blocks of width s = 2m+1),
and its length g s is below r, so |Cov| = g(2m+1). The algorithm outputs x exactly when x is in
Cov or x = r - s:

* x = i s + e in Cov with e = 0: A_i = O, caught in step 4 or by the step-5 test at index i-1 (Lemma 3).
* e != 0: A_i = [e]G, whose key equals the key of baby |e| (1 <= |e| <= m), so the merge finds a common key
  (Lemma 4) unless an earlier test already produced a verified output.
* x = r - s (not in Cov because g s + m < r - s): A_0 = Q = -S, caught by the step-5 test at i = 0.
* Every output passes [x]G == Q, so outputs are never wrong. Conversely, an equal-key pair (i, j) gives
  Q - [i s]G = +-[j]G, so x = i s +- j lies in Cov; a step-5 or step-4 detection gives x in Cov or x = r - s.
  Hence no other scalar is recovered.

Over the uniform challenge x in 0..r-1 (including x = 0, T = O, recovered in step 4 with c = 0):

Pr[success] = (g s + 1)/r = 0.40000000000000000000... >= 2/5 >= 0.39.

The algorithm is deterministic, so there are no algorithm coins to average over; this is an exact count,
not an estimate. `experiment.py` checks the analogous statement exhaustively on three toy groups.

## 4. Instruction model used for the non-reference work

Scalar AArch64 instructions with 64-bit registers, as the protected cost model requires. Conventions:

* A field element is 6 words. A 384-bit modular subtraction is 25 instructions: 6 subs/sbcs, 1 borrow mask,
  6 and (mask with p), 6 adds/adcs, 6 loads of p.
* 128-bit integers (indices, list positions, counters): add or increment 2 instructions, compare-and-branch 3,
  shift by 6 bits 3.
* The baby list, giant list and scatter buffer exceed 2^64 bytes. The cost model says "Larger arithmetic and
  data structures must account for their implementation", so addressing them is charged: every 64-bit word
  read or written in them costs 2 instructions (the access plus one instruction forming the high part of a
  128-bit address).
  Small arrays (chain state, digit counters, below 2^40 bytes) are charged 1 instruction per word.
* Every compare is charged its full worst case (all 6 limbs); no early exit is credited.
* Field inversion by Fermat: 380 squarings and 228 multiplications (p-2 has 381 bits, 229 of them set),
  608 x 383 = 232,864 plus 5 per exponent bit; charged 240,000.
* A textbook affine group-law operation with all cases (identity, P = -Q, doubling, generic; G1 has odd order,
  so a doubled point never has y = 0): 1 inversion,
  at most 4 multiplications, 10 subtractions or compares, 50 control; charged 250,000.
  A double-and-add scalar multiplication (scalar < 2^255): 255 x (2 x 250,000 + 24); charged 128,000,000.

## 5. Itemised per-element charges (instructions, worst case)

Each batched incomplete affine addition costs 1 reference unit (its six multiplications under the 6M
convention). In addition, and conservatively, we charge everything else it executes:

| per batched addition (beyond the reference unit) | instr |
| --- | ---: |
| 6 modular subtractions outside the 6M (dx, dy, lambda^2-x1, -x2, x1-x3, -y1), 6 x 25 | 150 |
| load chain point x1, y1 (12 words) | 12 |
| load constant addend x2, y2 (12 words) | 12 |
| store dx and prefix product, reload both in the backward pass (24 words) | 24 |
| store result x3, y3 (12 words) | 12 |
| loop control | 3 |
| **total (baby)** | **213** |
| giant only: exceptional test x(A) == x(S): 6 eor, 5 orr, 1 branch | 12 |
| **total (giant)** | **225** |
| plus one inversion per round of K additions | 240,000 / round |

| per emitted record (each of the m + g records) | instr |
| --- | ---: |
| canonical reduction of x: 6 loads of p, 6 subs/sbcs, 6 csel | 18 |
| histogram update for all 12 digits: extract 1, address 1, load counter 2, add 2, store 2 (8 each) | 96 |
| store 8-word record (16), 128-bit pointer advance (2), 128-bit index update (2) | 20 |
| **total** | **134** |

| per record per LSD pass (12 passes per list) | instr |
| --- | ---: |
| load 8-word record (2 per word) | 16 |
| extract 32-bit digit | 1 |
| counter address | 1 |
| load 128-bit counter | 2 |
| destination address: 128-bit shift by 6 (3), 128-bit add (2) | 5 |
| increment counter and store it | 4 |
| store 8-word record (2 per word) | 16 |
| source pointer advance (2), end test with branch (3) | 5 |
| **total per pass** | **50** (600 for 12 passes) |

| per merge step (m + g steps charged) | instr |
| --- | ---: |
| load both 6-word keys (2 per word) | 24 |
| compare 6 limbs from the top, cmp + branch each, final branch | 13 |
| 128-bit pointer advance | 2 |
| two 128-bit bounds tests with branch | 6 |
| **total** | **45** |

Per digit and list: zeroing 2^32 counters (3 each) and the prefix sum (9 each: load 2, add 2, store 2, loop 3).

Per element: baby 213 + 134 + 600 + 45 = 992 instructions, giant 1004 instructions (0.496 and 0.502 units).

## 6. Total work

H = K(Lb-1) + K(Lg-1) = 204814228410832858487360404162884927488 batched additions (2^127.267584).

| other work | instructions | log2 units |
| --- | ---: | ---: |
| baby additions, extra instructions: K(Lb-1) x 213 | 21856369859427656132230575270919472676864 | 123.0394 |
| giant additions, extra instructions: K(Lg-1) x 225 | 22995486752196911329835060720889102336000 | 123.1127 |
| inversions: (Lb-1 + Lg-1) x 240,000 | 46878256624793897664038178443043120000 | 114.1745 |
| record emission: (m+g) x 134 | 27445106607051603037306294157826861301760 | 123.3679 |
| LSD scatter: (m+g) x 600 | 122888537046499715092416242497732214784000 | 125.5306 |
| counters: 2 lists x 12 digits x 2^32 x 12 | 1236950581248 | 29.2041 |
| merge: (m+g) x 45 | 9216640278487478631931218187329916108800 | 121.7937 |
| starts and S: (2K+1) x 128,000,000 + K x 250,000 | 268697728000000 | 36.9672 |
| final verification (2 scalar mults) and I/O | 256001000 | 16.9658 |
| **R (total other instructions)** | **204449018800288158121383429283075544909672** | |

W = H + R/C = 614077475621953875096104237608845399885672 / 2000 = 3.0704 x 10^38 = 2^127.851688.

Score: the smallest k/100 with W <= 2^(k/100) is 12786/100 (checked exactly in integers: W^100 <= 2^12786
and W^100 > 2^12785), so time_log2 = 127.86 = 6393/50. The Q-independent baby phase is 2^126.8298 units of this total;
it is charged in full, not amortised. The reference additions are 2^127.267584 and all other work
R/C = 2^126.265 units, so W/H = 1.4991.

## 7. Memory

Babies are built and sorted first (list plus scatter buffer, 128 m bytes), then the giant list and its
scatter buffer are built while the sorted baby list is kept. Peak = 64 m + 128 g + digit counters
(12 x 2^32 x 16 = 824633720832) + chain state (K x 256 = 268435456) + code and constants (2^20)
= 19649049072251535499233038694524862332928 bytes = 2^133.8516. There is no advice; nothing is precomputed
outside the charged run.

## 8. Scope of the toy experiment

`experiment.py` runs the algorithm of section 2 step for step on y^2 = x^3 + b over primes 1009, 4003 and 12043
(prime group orders 967, 4111, 11827), choosing Lb and Lg with the same routine as production (`choose`). The
differences from production are only of size: K = 4 chains instead of 2^20, radix digits of 4 bits instead of
32 (so 3 or 4 passes instead of 12), keys are the reduced integer x instead of 6 limbs, and arithmetic uses
Python integers. For every scalar of each group it checks correct outputs, the exact success set
{i s + e} plus n - s, that no baby addition is exceptional (the batched addition raises on a zero
denominator), and that a full-length run performs exactly K(Lg-1) giant additions, g x (passes) radix moves
and at most m+g-1 merge steps, matching the ledger formulas. It runs in about one second.

## 9. Assumptions and limitations

* A1, memory: the machine has a flat memory of at least 2^133.86 bytes addressed with 128-bit addresses,
  charged as in section 4, applying the cost model's rule that larger data structures account for their
  implementation. This is far beyond any physical machine; memory carries no score under the rules.
* A2, prices: the reference unit and C = 2000 are the protected convention. The other counts are our
  analytic, itemised instruction counts for a straightforward implementation, not machine code or
  measurements; every compare and loop is charged at its worst case.
* The 213 extra instructions per addition are a conservative over-charge: the 6M convention with A = 0
  arguably already includes the subtractions and memory traffic. Dropping them would lower the score.
* Negation is the only symmetry used (the order-3 endomorphism is not used).
* No full-size run. The toy experiment checks the algorithm and the counting formulas at toy scale only;
  the production success probability follows from the counting argument in section 3, not from the toy.

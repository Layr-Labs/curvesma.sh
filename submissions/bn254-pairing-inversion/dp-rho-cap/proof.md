# bn254-pairing-inversion: capped distinguished-point rho

Submitted score: time_log2 = **126.322105** (exact rational 25264421/200000).
This bounds the total computation of **every** run, whatever its coins, in units
of C = 20,000 scalar AArch64 instructions. There is no run gate and no
expectation in the score. The success probability is at least 0.39 over
uniform z and fresh algorithm coins.

This is a generic square-root attack: a 1,024-adding Pollard-rho walk in GT
with distinguished points, one walk on one processor, about 2^45 bytes of
memory with ordinary 64-bit addresses. It is not a finite-field sieve
estimate, a practical break or a security-level estimate. No earlier
submission or acceptance is used as a premise.

## 1. Reduction and output

Let p, r, P0, Q0 and g = e(P0, Q0) be the protected constants. GT has prime
order r and g != 1, so every valid input is z = g^x for a unique x in [0, r).
If z = 1 the algorithm outputs P = infinity (null) and Q = Q0, which the rules
permit exactly for z = 1. Otherwise it recovers x, checks g^x = z, and outputs
P = [x]P0 and Q = Q0. By bilinearity e([x]P0, Q0) = g^x = z. G1 = <P0> has
cofactor 1, so P is in G1, and Q0 is the profile generator of G2. Choosing a
public generator for one point is allowed by the rules.

## 2. Algorithm

Constants:
- L = 1,024 adders; D_max = 2^36 distinguished-point records;
- s = ceil(sqrt(r)), N0 = ceil(27s/25), Delta = 2^102, cap N = N0 + Delta.

Element encoding: the tower Montgomery representation of machine.md (48 words).
It is a fixed bijection on Fp12, so equality of encodings is equality in GT.
For y in GT with pinned-basis coefficients c_0..c_11:
- part(y) = low 10 bits of c_6 * 2^256 mod p (word 0 of A0.im);
- y is distinguished if (c_0 + 9 c_6) * 2^256 mod p is divisible by 2^94
  (word 0 of A0.re is zero and the low 30 bits of word 1 are zero).

Steps:
1. **Input.** Decode z canonically, check z != 0 and z^r = 1. If z = 1, output
   (null, Q0) and stop.
2. **Coins.** Draw a0, b0 and c_j, d_j (j < L) uniformly from [0, r) by 254-bit
   rejection sampling with at most 1,024 draws each; on exhaustion, fail.
3. **Table.** Build records M_j = g^{c_j} z^{d_j} with labels (c_j, d_j), and
   Y = g^{a0} z^{b0} with labels (a, b) = (a0, b0). The invariant is Y = g^a z^b.
4. **Walk.** For t = 1..N: j = part(Y); Y <- Y M_j; (a, b) <- (a + c_j, b + d_j)
   mod r; if Y is distinguished, append (Y, a, b) to the table D, and stop the
   walk if D now holds D_max records.
5. **Collision.** Sort D by encoded element and scan adjacent records. For at
   most 16 pairs (Y, a, b), (Y, A, B) of equal elements: if B != b, set
   x = (a - A)/(B - b) mod r; if g^x = z, go to step 6.
6. **Output.** Output P = [x]P0 and Q = Q0. As attack-side verification,
   evaluate the reduced pairing e(P, Q0) and compare it with z on all 12
   coefficients. If no pair succeeds, fail.

Correctness: g^a z^b = g^A z^B gives a + x b = A + x B mod r, so x is as in
step 5 whenever B != b. Step 5 checks g^x = z, so a wrong x is never output.

## 3. Work

Every row is an every-run bound; none depends on a heuristic.

| Row | Instructions |
|---|---|
| walk: N steps at 13,312 | N * 13,312 |
| distinguished-point appends: at most D_max at 128 | 2^43 |
| everything else (machine.md table) | F = 2^52 |

So every run uses at most
R = N * 13,312 + 2^36 * 128 + 2^52
  = 2127024673319068928725389107523574384819200 instructions,
and log2(R/20,000) = 126.32210437426..., rounded up to **126.322105**.

The step charge comes from the emitted instruction sequence of 13,166
instructions (machine.md), including the Fp12 product, both label additions,
the distinguished-point test, the 2-word counter and the loop branch.

Preprocessing (steps 1 to 3) is below 2^37 instructions, or 2^23 units, and
is included. Peak memory is below 2^45 bytes (machine.md). There is no advice
and nothing is amortized over other targets.

## 4. Success probability

z = 1 succeeds at once. For z != 1, a run succeeds when all of the following
hold:

- **E1 (coins).** Rejection sampling succeeds: each draw accepts with
  probability r/2^254 > 1/2, so failure is below 2,050 * 2^-1024.
- **E2 (collision time).** The walk Y_0, Y_1, ... has tail length mu and cycle
  length lambda with mu + lambda <= N0.
- **E3 (detection).** Y_{t1} is distinguished for some t1 with
  max(1, mu) <= t1 <= mu + Delta. Then Y_{t1} is on the cycle and is recorded
  at t1 and again at t1 + lambda <= N0 + Delta = N.
- **E4 (capacity).** D does not fill before time t1 + lambda. Take t1 to be the
  first distinguished time in [max(1, mu), mu + Delta]. The elements Y_1, ...,
  Y_{mu+lambda-1} are distinct, and the times mu + lambda, ..., t1 + lambda - 1
  revisit the non-distinguished elements Y_mu, ..., Y_{t1-1}. So every record
  stored before time t1 + lambda belongs to a distinct element.
- **E5 (non-degeneracy).** For the first pair tried, B != b.

E5 holds with probability exactly 1 - 1/r, without any heuristic. Fix the
hidden x. Condition on the increments s_j = c_j + x d_j and the start
s_0 = a0 + x b0. The walk on group elements depends only on these. For each j
there is exactly one c_j for each d_j, so the d_j remain independent and
uniform. Any two records of the same element are n steps apart, with
0 < n <= N < r, and B - b = sum_j n_j d_j with sum_j n_j = n. Some n_j lies in
[1, r-1], so B - b is uniform and is 0 with probability 1/r.

**Heuristic H1 (r-adding random-mapping model).** For a fixed adder partition,
fresh uniform adders and a uniform start, the tail and cycle lengths
(mu, lambda) of the walk follow the random-mapping model with the conventional
r-adding correction. In particular the first repeat satisfies
Pr[mu + lambda <= t] >= 1 - exp(-(1 - 1/L) t^2/(2r)) for t of order sqrt(r).
Given mu + lambda = rho, the cycle length lambda is uniform on [1, rho].
This is the standard random-mapping estimate with the conventional r-adding
correction (Teske, "On random walks for Pollard's rho method", Math. Comp. 70
(2001) 809-825, studies such walks; the estimate is stated here, not imported
as a theorem). At t = N0 it gives 1 - exp(-0.58263...) > 0.441572.

**Heuristic H2 (distinguished points).** Along the walk, the distinguished
property behaves like independent events of probability 2^-94, independent
of the partition (the two read different coordinates: c_0 + 9 c_6 and c_6).
Under H2 and the random-mapping model of H1:
- E3 fails only if lambda <= Delta or the Delta cycle elements from time
  max(1, mu) on hold no distinguished point. Under H1, Pr[lambda <= Delta] is
  Delta E[1/rho] = Delta sqrt(pi (1 - 1/L)/(2r)) < 1.26 Delta / sqrt(r) < 2^-23.
  The second probability is at most exp(-2^8) < 2^-255.
- E4: by the remark under E4, the records before time t1 + lambda are
  distinguished elements among at most N distinct elements. Under H2 their
  number is a sum of independent indicators with mean at most
  N 2^-94 < 2^32.91, against D_max = 2^36 (more than 8.5 times the mean). The
  Chernoff bound gives failure below 2^-1000.

Hence the success probability is at least
0.441572 - 2^-23 - 2^-255 - 2^-1000 - 2^-1012 - 1/r > 0.4415, and the claim
needs only **0.39**. H1's prediction may fall short by 0.05 absolute (more than
11 percent relative) without affecting the claim.

H1 and H2 affect success only. If they fail, the run still stops within its
every-run work bound. The partition and the distinguished rule are fixed public
functions of the element. The adders and the start are fresh uniform coins.
Limitations: H1 is not a theorem for this partition and field; the toy runs
below support it only at toy scale.

## 5. Evidence

`experiment.py` is the single declared experiment (standard library only):

- **kernel**: runs the instruction sequences emitted by `a64.py` on an
  emulator. Montgomery product, addition and subtraction are checked on edge
  and random inputs, and the tower Fp12 product against the pinned-basis
  product. It prints the count of every priced sequence, and checks that the
  encoding guards reject seven malformed instructions.
- **production**: runs 12 emitted walk steps on production values with a
  1,024-record table. It checks the partition, the record address, the product,
  the labels and Y = g^a z^b, plus the cap-exit branch, the
  distinguished-point predicate on crafted words and one record append. These
  table labels are structured for speed; this checks mechanics, not statistics.
- **constants**: checks that p, r, x and the twelve coefficients of g used by
  the code equal those in the shipped `target-profile.json`, a byte-identical
  copy of the protected profile (sha256 649e2c45...3c1ed0).
- **ledger**: recomputes s, N0, N, R, the score and its upward rounding, and
  the success floor. It reads the claim in `submission.json` and compares each
  declared figure with the recomputed value: the exact score and its text, the
  success level, the step charge and emitted count, and the preprocessing and
  memory bounds.
- **toy**: runs the same walk rule (L adders with L in {16, 64, 1024}; the
  same partition, now on the low log2(L) bits, and the same distinguished-bit
  rule on the same Montgomery tower encoding) in a small BN-shape target group.
  That group uses the same field polynomial w^12 - 18w^6 + 82 with x = 5
  (r = 27,481) or x = -7 (r = 74,929). Each trial has a fresh target and fresh
  coins. The experiment prints every outcome, the run-success fraction and the
  H1 predictions at N0 and at the cap for the requested input. The toy uses
  2 distinguished bits and Delta = 32 because 94 bits cannot be scaled down; the
  cap formula is otherwise the same. With L = 16 the toy walk reuses its adders
  heavily, which exercises the r-adding regime of H1.

The claim does not rest on any toy figure. This is evidence at its own scale
only, and it does not establish H1 at production size.

## 6. Limitations

- Success rests on H1 and H2, with the margin stated above.
- Instruction counts are exact for the emitted step and analytic upper bounds
  elsewhere. They are not CPU timings.
- 2^126 work is far from practical, and nothing was run at full size.
- Finite-field sieve methods are expected to be cheaper for this field size;
  this candidate makes no such claim and charges only its own method.

# Phase-coupled six-class rho

Status: exploratory claim. Parent tree `6423c2728dc03d5f2d647c7f22adf1dac0643a66`.
The protected profile and `classical-group-work-v2` are bound by the SHA-256
values in `claim.json` (raw file bytes). Reference price C = 1000 scalar
AArch64 instruction equivalents per batched incomplete affine addition
(6M, S = M, A = 0 inside that formula only). Field multiplication is the
measured 179-instruction kernel. Additions outside that formula are counted
below. No quantum step, no uncharged advice, no amortization over future
targets.

## Group, orbits, state

Let p and n be the protected field modulus and prime order, P the protected
generator, and Q = [x]P the challenge. λ is the integer in 1..n-1 with
λ² + λ + 1 ≡ 0 (mod n); it exists because n ≡ 1 (mod 3). β is the smaller
root of T² + T + 1 ≡ 0 in F_p; it exists because p ≡ 1 (mod 3). The
automorphisms are φ(X, Y) = (βX, Y) and ι(X, Y) = (X, -Y), with [λ]R = φ(R).
The unit group U = {±1, ±λ, ±λ²} has order 6.

No non-identity unit fixes a finite point: φ(R) = ±R forces β = 1 or a
2-torsion point, and secp256k1 has cofactor 1 and odd order, so the only
2-torsion point is the identity. Thus every non-identity point has orbit
size 6, and 6 divides n - 1 (n ≡ 1 (mod 6)). There are m = (n - 1)/6 orbits.

A canonical representative of a finite point R = (X, Y) is computed as
follows. Set X0 = X, X1 = βX, X2 = -X - βX (the last equality is
β² + β + 1 = 0). Choose the smallest Xi, breaking ties by smallest i.
Negate Y if the resulting Y is greater than p - Y. The six images are
exactly the orbit, so this is a complete set of representatives. The
canonical unit u(R) ∈ U is the unique unit with can(R) = u(R) · R.
One field multiplication produces βX; X2 is two field subtractions.

The walk state is a pair (C, i) where C is canonical and i ∈ {0, 1, 2}.
The state space has size

    M = 3m = (n - 1)/2.

Coefficients are a frame (U, A, B) with U ∈ U and

    C = [U A] Q + [U B] P.

U starts at 1. Updating U is a 6×6 table lookup. A and B are integers mod n.

## What canonicalization does not do

Let V be a fixed jump and write g(C) = can(C + V). If u(C + V) = -1 and
u(g(C) + V) = -1, then

    g(g(C)) = can( -(C + V) + V ) = can(-C) = C,

and g(C) ≠ C whenever - (C + V) is not already C. The identity uses only
ι² = 1. It does not need V to be special, and it survives reduction onto
orbit representatives. A walk whose partition also repeats, h(g(C)) = h(C),
then selects the same V again and stays on the 2-cycle forever. That set is
a positive-density trap for an unphased canonical walk whenever the unit
and the partition behave even roughly uniformly (density on the order of
1/(36 r)). Distinguished points do not see the trap if neither cycle point
is distinguished: the walker never stores a new point. This is the
obstruction. The algorithm below does not claim that canonicalization
removes it.

## Jumps

Sixteen slots, three phases, forty-eight exponents

    t(i, s) = 3^i · 17^s    (i ∈ {0,1,2}, s ∈ {0..15}).

These 48 integers are pairwise distinct: 3^a 17^b = 3^{a'} 17^{b'} with
0 ≤ a ≤ 2 and 0 ≤ b ≤ 15 forces (a, b) = (a', b'). Each is at most
9 · 17^15 < 2^66.

The global jump for phase i and slot s is the point

    J(i, s) = Q + [t(i, s)] P,

coefficient vector v(i, s) = (1, t(i, s)) in the (Q, P) basis. The
partition of a canonical point is the integer s formed by bits 20..23 of
its X coordinate (four bits above the distinguished-point mask). The low
20 bits are reserved for the distinguished-point test, so the partition
does not read them.

Precomputation, charged in full: [t]P for each of the 48 exponents
(at most 66 doubles and 33 adds each, under 100 reference additions),
then one addition of Q. Coefficient addends σ(u, i, s) = u^{-1} · (1, t(i, s))
mod n, six units times 48 exponents, each a constant multiplication by a
unit modulo n. Bound this whole table, together with decoding Q and
building the 64 starting states, by 2^20 reference equivalents. That is
the preprocessing budget. It is included in the total.

## The iteration

Sixty-four walks run together so that each parallel step has 64 independent
denominators and the batched inversion in the 6M reference price is the
regime being charged. Summing processor work, not wall time.

From state (C, i) with frame (U, A, B) and slot s:

1. Add J(i, s) with one batched incomplete affine addition. This is one
   reference operation. The exceptional case (denominator 0: doubling,
   inverse points, or identity) uses a complete affine fallback charged as
   two reference operations. A reserve of 2^12 reference operations covers
   every fallback the budget allows; if a walk would need more, that walk
   aborts. The reserve is in the total below.
2. Replace the sum by its canonical representative. Multiply the frame unit
   by the canonical unit (table lookup). Add the precomputed σ(U_old, i, s)
   into (A, B) before the unit update, which is two 256-bit modular adds.
   The geometric point added is J, independent of U, so the iteration is a
   function of (C, i) alone. Two walks that reach the same state follow the
   same future.
3. Set i ← i + 1 (mod 3).

The frame invariant is preserved: adding (1, t) in the global basis is
exactly adding U^{-1}(1, t) inside the frame, and canonicalization
multiplies the frame unit.

Distinguished points: after canonicalization, if the low 20 bits of X are
0, store the record (X, Y, A, B, U, i) in an open-address table sized for
load factor at most 1/2. A hit with the same (X, i) and invertible frame
difference

    den = U A - U' A',    num = U' B' - U B    (mod n)

yields k ≡ num · den^{-1} (mod n). The attack then checks [k]P = Q with one
scalar multiplication (at most 256 doubles and 256 adds, covered by the
2^12 reserve). A singular den is not a solution; the walk continues.

## No positive-density fruitless 2-cycle

**Phase.** F advances i by 1. Then F² advances i by 2 ≢ 0 (mod 3), so
F²(C, i) ≠ (C, i) for every state. There is no 2-cycle in the state graph
at all.

**Identical point relation.** A fruitless cycle that holds for every
challenge scalar, not just one k, is an identity in the coefficient ring.
For two steps with units u1, u2 and jumps v, v',

    u2 u1 (c + v) + u2 v' = c

for a positive-dimensional set of c forces u2 u1 = 1 and v' = -u1 v.
Every jump has Q-coefficient 1, so -u1 · (1, t) = (1, t') forces -u1 = 1,
hence u1 = -1 and t' = t. Consecutive steps use consecutive phases, and
distinct phases have disjoint exponents, so t' = t is impossible.
Therefore there is no identical fruitless 2-cycle.

Isolated solutions of the same equation with u2 u1 ≠ 1 are at most one
point per unit pair and slot pair: at most 36 · 16² < 2^14 points. A walk
of length T < 2^128 hits a fixed set of that size with probability at most
2^128 · 2^14 / M < 2^{-100}. Those points are not a trap class.

## No identical period-3 cycle

Three steps with unit product U = u3 u2 u1 = 1 and inhomogeneous part 0
require, in Z[λ],

    1 + u3 u2 + u3 = 0.

The only triples of units that sum to 0 are permutations of (1, λ, λ²) and
of (-1, -λ, -λ²). The displayed triple contains 1, so the negative triple
is impossible. The remaining cases are (u2, u3) = (λ, λ) and (λ², λ²), and
both force the three exponents to be equal. The three steps of a period-3
state cycle sit in three different phases (the counter runs through all
residues), hence on three distinct exponents. Contradiction. There is no
identical fruitless 3-cycle.

## Heuristics

**H_map.** On the M states, for uniform Q, the iteration behaves as a
random mapping for the purpose of the first collision time: after T steps
in total across the walks,

    Pr[no collision] ≤ exp(-T² / (2M)).

Support: this is the Pollard rho / van Oorschot–Wiener collision law used
for automorphism walks once short fruitless cycles are excluded. The jumps
depend on Q through J = Q + [t]P, which is the usual r-adding dependence.
Limitation: it is not a theorem about this specific partition.

**H_long.** Identical fruitless cycles of length at least 6 remove at most
0.02 of the success probability. Support: lengths 2 and 3 are classified
above and do not occur; a longer identical cycle is a unit-linear relation
on the exponent set of length at least 6, and the published negation-map
computations (secp112r1, ECC2K-130) treat cycles beyond a small length as
a lower-order effect once the structural short cycles are gone. Limitation:
this note does not classify relations of length ≥ 6. If a positive-density
family of length 6 or more exists, the success claim can fail.

**H_dist.** Along the walk, the fraction of canonical X divisible by 2^20
is at most 2^{-18} and at least 2^{-22}. The high side bounds memory and
distinguished-point write work. The low side sets the detection tail:
after two walks merge, L = 40 · 2^22 further steps on the merged walk miss
a distinguished point with probability at most exp(-40) < 10^{-17}.
Sixty-four walks times that tail is under 2^34 steps, absorbed by the
ceil in the sixth decimal of the score (see the numerical bound). Support:
X mod 2^20 for a canonical representative of an adding walk is the standard
equidistribution hypothesis, with a factor-four slack on each side.
Limitation: a partition that correlated the low 20 bits with the slot bits
would break it; the slot is taken from bits 20..23, disjoint from the mask,
but disjointness is not a proof of equidistribution.

**H_useful.** Conditioning on a distinguished-point collision of two
distinct frames, den ≡ 0 with probability at most 0.01. Support: den = 0
is one linear condition on the frame difference. A useful collision is the
generic case for two walks that did not follow the same coefficient
history. Limitation: heavily dependent frames could be more singular than
this.

Composing H_map at 1/2, H_long at 0.02, H_useful at 0.01, and a 0.01
allowance that covers detection (the exp(-40) tail is far smaller):

    (1/2) · 0.98 · 0.99 · 0.99 = 0.480249.

The 2^{-100} isolated-cycle loss does not move the third decimal. The
claimed success probability is 0.39. A fallback storm that exhausts the
2^12 complete-addition reserve is the same shape of event as hitting a
single forbidden field value; it sits inside that 2^{-100} allowance and
is not given a separate heuristic.

## Per-step instruction budget

Inside the reference addition, do not count the 6M again. Outside it, one
ordinary step (no fallback):

| Piece | Instructions | Schedule |
| --- | ---: | --- |
| Multiply X by β | 179 | measured field multiplication |
| X2 = -X - βX | 96 | two field subtractions, 48 each |
| Conditional negate Y | 48 | one field subtraction |
| Min of three X and select | 72 | three 5-limb compares and a select |
| Two frame adds mod n | 80 | two 4-limb adc/conditional-sub, plus index |
| Load J and form the slot | 24 | address, four-bit slot, point loads |
| Unit-table update | 8 | 6×6 lookup |
| Phase increment | 4 | add and mask |
| Distinguished-point test | 8 | mask the low 20 bits |
| Walk loop | 32 | index, branch, stores of the live point |
| **Sum** | **551** | |

A field subtraction in the 5×52 representation is 10 loads of the operands,
5 subs, 5 borrow ops, 5 loads of p, 5 conditional adds, 5 stores, and 8
address ops: 48. The 179-instruction multiply is an upper bound for a
multiply by the fixed β.

Round 551 up to 700 instructions per step. The 149-instruction slack covers
the slot arithmetic, the open-address probe on the rare distinguished step
is not in this line (it is charged globally below), and small control that
the limb schedule under-counts. 700/1000 = 7/10 of a reference operation.

Distinguished-point writes, using the high side of H_dist: at most T / 2^18
records, each at most 400 instructions to hash, probe, and store. That is

    T · 400 / (1000 · 2^18) = T / 655360

reference equivalents.

## Integer total

Let ln 2 = Σ_{k≥1} 1/(k 2^k). The partial sum through k = 40 is

    S = 1017993251683865178040715701 / 1468653815862788621126860800

and the tail is strictly less than 1/(41 · 2^40), so

    ln 2 < S + 1/(41 · 2^40)
         = 41737723319039808032533609541 / 60214806450374333466201292800
         < 69314718056 / 10^11.

Let M = (n - 1)/2 and let T be the least integer with

    T² · 10^11 ≥ 2 M · 69314718056,

namely

    T = 283303653675694132582384390256226586611.

Then T² / (2M) > ln 2, so exp(-T² / (2M)) < 1/2 under H_map. The
distinguished-point tail of 2^34 steps changes log2(T) by far less than
10^{-20} and is inside the upward rounding below.

Per step the charged cost is one reference addition plus 700/1000, plus the
distinguished-point write term 1/655360:

    17/10 + 1/655360 = 11141130 / 6553600.

Add preprocessing 2^20 and the fallback/verification reserve 2^12:

    W = T · 11141130 / 6553600 + 2^20 + 2^12.

Then log2(W) < 128.501153, and the submitted score is that ceiling at six
decimal places:

    time_log2 = 128501153 / 1000000 = 128.501153.

Preprocessing is exactly the 2^20 term, so preprocessing_log2 = 20.
Verification and fallbacks sit in the 2^12 reserve inside W, not in
preprocessing.

## Memory

Records are 256 bytes. The table is allocated at twice the high-side
distinguished-point count ceil(T / 2^18), plus 2^24 bytes of code, batch
scratch for 64 walks, and the 48-point jump table. The resulting peak is
at most 2^118.735617 bytes. Nonuniform advice is 0 bytes: every table is
built in the charged preprocessing. Memory is not part of the score.

## Experiment

`experiment.py` runs in Z/13Z, where λ = 3. The unphased one-jump map
can(C + 10) swaps the canonical residues 1 and 2, and both steps have
canonical unit -1. That is the structural 2-cycle, on a group small enough
to check by hand. The same jump with a phase modulo 3 has no state of
period 2, by exhaustive check of the six states and by the phase proof
above. A planted pair of frames recovers k = 5. This says nothing about
production cost, memory, or H_map.

## References

The collision shape is Pollard rho as in HAC §3.6, with the negation and
endomorphism action as in Gallant–Lambert–Vanstone and the fruitless-cycle
discussion of Bos–Kaihara–Kleinjung–Lenstra–Montgomery (PKC 2010 / ePrint
2010/077). The cost model is the protected `COST_MODEL.md`. The citations
are context; the algorithm, the cycle proofs, and the integer bound are
above.

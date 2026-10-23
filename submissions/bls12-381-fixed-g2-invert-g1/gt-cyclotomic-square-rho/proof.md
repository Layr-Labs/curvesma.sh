# bls12-381-fixed-g2-invert-g1: Pollard rho with cyclotomic squaring

Status: ready for independent review; no qualification or emitted score.
The normative inputs and outputs are those of the protected target profile.
`claim.json` binds that profile and `classical-group-work-v2` by SHA-256.
No restricted challenge distribution, uncharged advice, or quantum step is used.
The only nondeterministic premise is heuristic H1 below.

## Reduction

The challenge gives a nonidentity Q in G2 and z in GT, and asks for P in G1
with e(P, Q) = z under the pinned pairing (the library takes G2 first).

Let r be the protected subgroup order and let P0 be the protected G1 generator.
Compute g = e(P0, Q) once. The pinned ate pairing is nondegenerate and bilinear,
so g = e(P0, Q0)^b where Q = [b]Q0 and b is coprime to r. GT has prime order r
and e(P0, Q0) is not the identity, so g generates GT. Hence z = g^a for a unique
a in [0, r), and P = [a]P0 satisfies e(P, Q) = g^a = z. The case a = 0 is
z = 1 and returns the point at infinity, which is the unique valid answer.
Nothing in the walk depends on b. The rules explicitly allow solving the
fixed-Q discrete logarithm in GT and returning [a]P0.

## Why the squaring branch is cheaper, and why it is not the whole walk

GT sits inside the cyclotomic subgroup of Fp12*. Concretely, the protected
primes satisfy r | (p^4 - p^2 + 1): both p^4 - p^2 + 1 and p^4 - p^2 + 1
reduced modulo r are 0. Every element of GT therefore obeys
x^(p^4 - p^2 + 1) = 1, which is the hypothesis of the Granger–Scott
cyclotomic squaring formula. Walk states below are products of powers of g
and z, so they stay in GT and the formula applies to every squaring branch.
It does not apply to a generic Fp12 product, and this package never uses it
for one.

Squaring alone does not compute a discrete logarithm. On a group of odd prime
order the map x → x^2 is a bijection (multiplication by 2 is invertible modulo
r). A walk that only squares, or that only multiplies the exponent by powers
of 2, preserves the ratio of the two exponents in x = g^α z^β and never
produces a collision between distinct ratios. The iteration must sometimes
multiply by g or by z.

The iteration is therefore Pollard's three-way map, not a squaring orbit and
not a Frobenius-class walk. There is no orbit representative. Equality of
walk states is equality of GT elements.

## Algorithm

Let s = floor(sqrt(r)) and B = floor(6s/5) = floor(1.2 sqrt(r)).
Let M = ceil(11 B / 16). These are deterministic functions of the protected
order. Numerically, with the protected r,

  s = 228988810152649578064853576960394133503
  B = 274786572183179493677824292352472960203
  M = 188915768375935901903504200992325160140

and the complementary squaring budget N = B - M equals
85870803807243591774320091360147800063.
The ratio M/B is 11/16 = 0.6875, which is strictly above 2/3.

Setup, charged in full below:

1. Reject an input that fails the protected checks (Q on G2, in the subgroup,
   not infinity; z nonzero and z^r = 1).
2. Compute g = e(P0, Q) by the pinned pairing.
3. Draw one starting state x = g^α z^β with α, β uniform in [0, r).
   The bound below allows up to 2^20 independent restarts; each restart is a
   fresh pair (α, β) and a fresh walk. Restarts are how a degenerate
   collision is discarded.

Walk, from a state (x, α, β) with x = g^α z^β. Let t be the least nonnegative
residue modulo 3 of the first canonical Fp coefficient of x (96 hex digits in
the protected encoding, interpreted as an integer). Then:

- t = 0: replace x by its cyclotomic square, and replace (α, β) by (2α, 2β)
  modulo r. This is one cyclotomic squaring and two modular doublings.
- t = 1: replace x by x · g, and replace α by α + 1 modulo r.
  This is one dense Fp12 multiplication.
- t = 2: replace x by x · z, and replace β by β + 1 modulo r.
  This is one dense Fp12 multiplication.

The dense product is the reference target-group operation. The squaring is
not, and is priced from the instruction bound in the next section.

Distinguished points and the collision. A state is distinguished when bits
[8, 88) of the canonical encoding are all zero. On each distinguished state
the walk inserts (x, α, β) into a table keyed by x. A second arrival of the
same x with exponents (α', β') gives

  α - α' ≡ a (β' - β)  (mod r),

so a = (α - α') · (β' - β)^(-1) modulo r whenever β' ≠ β. If β' = β the
collision is degenerate: discard it and restart. On a candidate a, compute
P = [a]P0 and accept only if the pinned pairing satisfies e(P, Q) = z.
A rejected candidate is a failed restart, not a wrong answer.

Hard stop. The walk keeps two counters, steps and muls. Every branch
increments steps. Branches t = 1 and t = 2 also increment muls. The walk
stops, and the instance fails, when steps reaches B or muls reaches M,
whichever comes first. It also stops on the first accepted candidate.
Because a dense multiplication costs more than a cyclotomic squaring, the
most expensive execution that still obeys the stop rule is exactly N
cyclotomic squarings and M dense multiplications (the multiplications are
spread so the multiplication counter hits M only as steps hits B). Every
cheaper execution is bounded by that one. The bound is a cap on what the
algorithm is allowed to do, not an average over branches.

## Instruction bound for one cyclotomic squaring

The straight-line formula, written out so the operation count does not depend
on an external listing, is the Granger–Scott formula as used by the
cyclotomic square in the tower Fp2 / Fp6 / Fp12 with u^2 = -1. An element is
six Fp2 coordinates (a00, a01, a02, a10, a11, a12). The formula evaluates
three Fp4 squares on the pairs (a00, a11), (a10, a02), (a01, a12), then a
fixed pattern of Fp2 additions, subtractions, and doublings, plus two
multiplications by u + 1 (each of which is one Fp subtraction and one Fp
addition, not a field multiplication).

One Fp4 square on a pair (p, q) of Fp2 elements is:

- three Fp2 squares (on p, on q, and on p + q);
- four Fp2 additions or subtractions;
- one multiplication by u + 1.

Three Fp4 squares therefore contribute nine Fp2 squares. One Fp2 square, for
u^2 = -1, is the two-multiplication identity

  (a + b u)^2 = (a + b)(a - b) + (2 a b) u,

which is two base-field multiplications and at most six base-field
modular additions or subtractions (the sum, the difference with its
conditional correction, and the doubling with its correction).

Adding the three Fp4 squares to the outer cyclotomic pattern
(three blocks of "subtract, double, add" on the a0-side, one multiplication
by u + 1, and three blocks of "add, double, add" on the a1-side) gives the
total

  18 base-field multiplications and 122 base-field modular additions
  or subtractions,

together with addressing. Those 122 split as 84 from the three Fp4
squares (six modular additions inside each of nine Fp2 squares, four Fp2
additions or subtractions per Fp4 square, and one multiplication by u+1
per Fp4 square) plus 38 from the outer pattern (three subtract-double-add
blocks, one multiplication by u+1, and three add-double-add blocks). No
base-field inversion appears.

Prices, each an upper bound:

- One base-field multiplication: at most 400 scalar AArch64 instructions.
  The protected cost model records 383 instructions for the pinned blst
  base-field multiplication on the same AArch64 measurement setup
  (`tracks/docs/COST_MODEL.md`). The extra 17 instructions are margin for
  a call and for moving the two operands into the kernel's registers.
- One base-field modular addition or subtraction of two canonical 6-limb
  residues, with the modulus already in registers: 3 `ldp` of the first
  operand, 3 `ldp` of the second, 6 `adds`/`adcs` or `subs`/`sbcs`, 6
  `subs`/`sbcs` against the modulus, 6 `csel`, 3 `stp`, plus at most 12
  spill or address instructions. That is 39 instructions. The bound used
  here is 48, which leaves 9 instructions of margin.
- Addressing, calls, and stores of the twelve output limbs that are not
  already inside those 122 additions: at most 2048 instructions.

Hence one cyclotomic squaring costs at most

  18 · 400 + 122 · 48 + 2048 = 7200 + 5856 + 2048 = 15104

scalar instructions. The ledger charges 16384 instructions per squaring,
which sits 1280 instructions above that count. Bookkeeping around the
branch (the residue modulo 3, the distinguished-point test, two 4-limb
exponent updates, the counter increments) is not included in the 16384 and
is charged separately at 64 instructions on every step, squaring or not.

A dense branch is one reference target-group operation, by definition of
the unit, plus the same 64 instructions of bookkeeping. Internal
instructions of the reference multiplication are not counted again.
C = 20000.

## Work bound

Step work, in scalar instructions, is at most

  I = N · (16384 + 64) + M · (20000 + 64)
    = N · 16448 + M · 20064
    = 5202808957716320533295925151401723028485184.

N and M are the integers in the previous section, and the product is that
integer. Dividing by C = 20000 converts the squaring instructions, and the
reference multiplications priced at 20000 instructions each, into
target-group-operation equivalents. I / 20000 is the step-cap cost.

One-time and tail work, bounded crudely and then added in full:

- two pinned pairings (the base g, and the accepting check) and the input
  subgroup checks: under 2^16 reference operations. A complete pairing is a
  63-bit Miller loop plus a final exponentiation; even at one reference
  multiplication per Miller bit plus a few hundred final-exponentiation
  squarings this is far under 2^16;
- at most 2^20 restarts, each one exponentiation of a 255-bit exponent.
  Charging every restart at 512 reference operations (255 squarings priced
  as if they were full multiplications, plus 255 multiplications) gives
  2^20 · 512 = 2^29 reference operations;
- one modular inversion of a 255-bit integer and the final scalar
  multiplication [a]P0: under 2^20 reference operations;
- the distinguished-point table: at most B / 2^80 insertions. Each
  insertion is a comparison of a 576-byte key and a store, under 2^12
  instructions, so the table costs under 2^60 instructions, which is under
  2^46 reference operations;
- generating the starting exponents: under 2^30 instructions.

The sum of these items is under 2^50 reference operations. Adding 2^50
reference operations to I / 20000 does not reach the next 10^(-6) of the
logarithm. Concretely, the total W in reference operations satisfies

  log2(W) < 127.612559,

and the gap from the computed value to 127.612559 is about 6.6 · 10^(-7),
larger than any contribution of the 2^50 term (that term changes the
logarithm by about 2^(50 - 127.6), which is under 10^(-23)).

Submitted score: time_log2 = 127.612559, rational 127612559/1000000.
Preprocessing is the pairings, the starting states, and the empty table
allocation, and is at most 2^30 reference operations, all included in W.
Peak memory is the distinguished-point table. At most B / 2^80 ≤ 2^47.70
entries are stored, each at most 640 bytes (576-byte GT element plus two
32-byte exponents), which is at most 2^57.02 bytes, plus under 2^30 bytes
of code and scratch. The submitted memory bound is 2^57.04 bytes.
Nonuniform advice is 0 bytes. No table depends on secret data, and the
public curve constants are excluded by the rules.

## Success probability

The probability space is the uniform challenge (uniform nonidentity Q,
uniform z) together with the algorithm's coins.

Heuristic H1. On a prime-order cyclic group, Pollard's three-way iteration
(square on one third of the elements, multiply by the base on one third,
multiply by the target on one third) has the collision-time law of a
random mapping: the probability that the first collision occurs after more
than t steps is exp(-t^2 / (2 r)) plus a lower-order term, for t on the
scale of sqrt(r). The same heuristic treats the residue modulo 3 of the
canonical encoding as unbiased, up to an additive error far below 2^(-64),
both for that collision law and for the branch counters. Support: this is
the iteration Pollard analysed in 1978, and it is the multiplicative form
of the additive walks whose collision times Teske measured against the
random-mapping prediction. It is not a theorem. Known discrepancies for
walks of this shape are constant factors on the mean, not a change of the
square-root exponent. The toy experiment runs the same three branches, the
same distinguished-point rule, and the same linear solve on a prime-order
additive group; it is evidence only at that size.

Under H1, with t = B = 1.2 sqrt(r),

  t^2 / (2 r) = 0.72,
  P(a collision by step B) = 1 - exp(-0.72) = 0.513247...

A collision is degenerate with probability O(1/r) per collision for uniform
exponents, which does not move the first six decimals. Each degenerate
collision spends one restart. The budget allows 2^20 restarts; only the
restarts that actually run are charged, and the cap already prices 2^20 of
them inside the 2^50 term above. The accepting pairing rejects a wrong
candidate, so a degenerate collision cannot be reported as success.

The multiplication cap cuts off walks whose branch counts are atypically
heavy. Under H1 the number of dense branches in B steps is binomial with
mean 2 B / 3. The cap M = 11 B / 16 lies a relative distance

  δ = (11/16 - 2/3) / (2/3) = 1/32

above that mean. The Chernoff bound gives

  P(muls > M) < exp(- μ δ^2 / 3)

with μ = 2 B / 3 > 2^127. Then μ δ^2 / 3 > 2^115, and the abort
probability is smaller than 2^(-2^100), which does not affect any figure
reported to six decimals. Walks that abort are failures, and they are
already charged at the cap.

The submitted success probability is 0.42. That is below the random-mapping
value 0.513 by more than 0.09, which is room for a constant-factor
inflation of the collision time of the size repeatedly observed for
three-way walks on small prime-order groups. The toy does not by itself
extrapolate to r; H1 states that extrapolation and its limit.

## What this package does not claim

It does not walk on Frobenius orbits, does not choose orbit representatives,
and does not price a squaring as a dense multiplication or a dense
multiplication as a squaring. Deterministic baby-step/giant-step on the same
group costs about 2^128.5 reference operations in the sibling accounting
control and uses a table of size about 2^137 bytes; the cap here is smaller
because each squaring branch is charged at 16448/20000 < 1 reference
operation and the number of dense branches is limited by M. Index-calculus
attacks on Fp12 are not claimed.

## Submitted bound

  time_log2                 = 127.612559
  success_probability       = 0.42
  preprocessing_log2        = 30
  memory_log2_bytes         = 57.04
  nonuniform_advice_bytes   = 0
  heuristics                = H1 only

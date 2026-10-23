# p256-ecdlp: deterministic negation-map BSGS with a complete cold-work ledger

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` (the P-256 `target-profile.json`) and
`classical-group-work-v2` by the SHA-256 of the current protected files
(`802d2350...` and `e0edc25a...`).

## Exact claim

Total work is **below 2^145** target-group-operation equivalents; the submitted
first-claim score is the exact rational **145/1**. As the P-256 rules require,
the borrowed unit C = 1000 scalar 64-bit AArch64 instructions is used only as a
normalization constant: **no reference operation is used (H = 0)**, every
P-256 field operation, point step, memory transfer and control instruction is
charged through an explicit analytic instruction bound, and the total
instruction count R is divided by 1000. No P-256 point addition is priced as
one unit and the secp256k1 6M formula is not used anywhere.

Preprocessing is below 2^145 equivalents (already inside the total). Peak
storage is below 2^150 bytes including code, routing state, both sort buffers
and scratch. Nonuniform advice is 0 bytes. Success is **strictly above 1/2**,
proved exactly on the uniform full-production distribution of the protected
challenge (x uniform in 0..r-1, including x = 0 / Q = infinity).

There is no heuristic: no random walk or mixing model, no birthday estimate, no
distinguished points, no fruitless cycles, no hash table or expected-time
lookup, no random bits, no rejection sampling, no retries, no quantum step and
no uncharged preprocessing. The algorithm is deterministic and every path,
including the failure path, is charged at its worst case. All resource bounds
are hard maxima, not expectations.

## Group facts used (P-256 specific)

P-256 is y^2 = x^3 - 3x + b over F_p with
p = 2^256 - 2^224 + 2^192 + 2^96 - 1, prime group order
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551 and
cofactor 1 (protected profile). There is **no efficiently computable
endomorphism** beyond negation (j != 0, 1728); none is used. The only symmetry
used is negation, valid on every short-Weierstrass curve: two finite points
share an affine x-coordinate iff they are equal or negatives, because for fixed
x the equation is quadratic in y with roots +/-y. Infinity gets its own tag and
is its own negative. So key(R) = (infinity_tag, x(R)) is constant on {R, -R}
and separates distinct classes. Because r is odd there is no point of order 2,
so y = 0 never occurs for a finite point of the group.

## Construction

Fix m = 2^127 (table records), stride s = 2m - 1 = 2^128 - 1, offset
c = m - 1 and h = 2^127 + 1 giant iterations.

1. Parse Q canonically (null or 64-hex x, y < p); reject malformed or off-curve
   input (check y^2 = x^3 - 3x + b); check [r]Q = O.
2. Baby table: for j = 0..m-1 store record (key(jP), j, y(jP)) by repeated
   complete affine addition R_{j+1} = R_j + P (R_0 = O is the infinity record).
3. Sort the table by key with deterministic bottom-up merge sort.
4. Compute D = -[s]P and gamma = Q - [c]P. For i = 0..h-1: binary-search
   key(gamma). On a hit with record j the stored y decides the sign:
   y(gamma) = y(jP) means gamma = jP and a = i*s + c + j; otherwise
   gamma = -jP and a = i*s + c - j (both mod r; for j = 0 the two coincide).
   Verify [a]P = Q and output a. Otherwise gamma <- gamma + D.
5. If the loop ends with no hit, output failure.

## Exact success

At iteration i, gamma = [t]P with t = x - i*s - c. A key hit occurs iff
t ≡ ±j (mod r) for some 0 <= j < m. For x in [0, h*s) and all i < h the
integers i*s + c ± j lie in [0, h*s) ⊂ [0, r), so no modular wrap occurs and a
hit at iteration i happens iff x is in [i*s, i*s + 2m - 2] = [i*s, (i+1)*s - 1].
These windows tile [0, h*s) with no gap and no overlap, so the first hit occurs
at i = floor(x/s) and returns the true x; the independent check [a]P = Q
accepts only the true x because P has prime order r. A wrong scalar is never
output.

No duplicate keys: for 0 <= j, j' < m, jP = ±j'P forces j ≡ ±j' (mod r); since
m - 1 < r/2 this gives j = j'. So a hit identifies a unique record.

Coverage is h*s = (2^127 + 1)(2^128 - 1) = 2^255 + 2^127 - 1 distinct scalars,
all below r. Exact integer checks (also in `experiment.py` mode `ledger`):
h*s < r and 2*(h*s) > r. Hence

**success = (2^255 + 2^127 - 1)/r = 0.5000000001... > 1/2 >= 0.39**

over uniform x in [0, r). The P-256 order is r = 2^256 - 2^224 + ..., so
r/2 < 2^255, which is why the same parameters as on secp256k1 still give
coverage above one half. Limiting coverage does not change the sampling
distribution; uncovered scalars are failures, and their full work is charged.

## Field kernels for the P-256 prime (four 64-bit limbs)

All counts are analytic worst-case upper bounds in scalar AArch64 instructions
(general-purpose registers only), including loads/stores, carries, address
arithmetic, call/return and control. They are not timings. The implementation
uses only that p is odd; it does not depend on the special Solinas shape of p
(that shape would only make things cheaper and no credit is taken).
Elements are kept fully reduced (< p) in Montgomery form with R = 2^256; since
x -> xR mod p is a bijection on F_p, the Montgomery-form x is itself a
canonical key. Inputs are converted in (one multiplication each) and the
output is a scalar, so no conversion out is needed except for equality tests,
which are done on the canonical Montgomery forms.

**Modular multiplication (CIOS Montgomery, 4 limbs).** For P-256,
n0' = -p^{-1} mod 2^64 = 1. Per outer iteration i (4 iterations):
the product row t += a*b_i is 4 x (MUL, UMULH, ADDS, ADC, ADDS, ADC) = 24, plus
2 for the top words; q = t_0*n0' is 1 MUL; the reduction row t += q*p is again
24 + 2; the one-limb shift is at most 5 register moves. That is 58 per outer
iteration, 232 in total. Add 8 operand loads, 4 moves materializing p, the final
conditional subtraction (5 SUBS/SBCS + 4 CSEL), 4 result stores and a 24
instruction call/prologue/epilogue allowance: **281 instructions**. The exact
CIOS routine with this tally is implemented in `experiment.py` (mode `kernel`)
on the real P-256 modulus and checked against big-integer arithmetic on 300
deterministic inputs including the edge case (p-1)^2. The ledger charges
**M = 8192** per multiplication or squaring, a factor > 29 above the derived
count, so the conclusion does not depend on scheduling details.

**Modular add/sub.** 8 loads, 4 ADDS/ADCS (or SUBS/SBCS) + 1 carry capture, 5
for the trial correction by p, 4 CSEL, 4 stores, under 32 with control. Charged
**A = 256**. Negation, comparison with zero/equality (4 EOR/ORR + CBZ) and
multiplication by small constants (as additions) are each charged as one A.

**Inversion.** Fermat, z^(p-2), with the public exponent p - 2 scanned
left-to-right: exactly 255 squarings and 127 multiplications (wt(p-2) = 128),
382 <= 512 multiplications; the chain is executed in `experiment.py` on the
P-256 generator x and inverts correctly. Charged
**Inv = 512*8192 + 2^14 = 4,210,688** (2^14 for loop and bit-scan control).
Inversion of 0 never occurs because every call is guarded (below).

**Complete affine point step** R3 = R1 + R2 (used for every baby step, giant
step and inside every scalar multiplication):

- If R1 = O return R2; if R2 = O return R1 (tag tests).
- If x1 = x2 and y1 + y2 = 0 return O (opposite points, incl. the guard for y=0).
- If x1 = x2 (then R1 = R2, doubling, P-256 has a = -3):
  num = 3*(x1^2 - 1) (1 squaring, 1 subtraction of the Montgomery one,
  2 additions), den = 2*y1 (1 addition).
- Otherwise num = y2 - y1, den = x2 - x1 (2 subtractions).
- lambda = num * Inv(den); x3 = lambda^2 - x1 - x2; y3 = lambda*(x1 - x3) - y1.

This uses one inversion, at most 4 multiplications/squarings and at most 10
additions/subtractions/comparisons, plus copies and branches. The a = -3
doubling formula is the one difference from the secp256k1 (a = 0) package; it
is checked on twelve small prime-order a = -3 curves and nine generic-a curves.
The ledger charges a much larger envelope:
**Inv + 19M + 32A + 2^12 = 4,210,688 + 155,648 + 8,192 + 4,096 = 4,378,624
< 2^23** instructions per point step (the 2^12 covers tag tests, branches,
limb copies and the record write of x, y into registers). No batched or
incomplete-addition credit is taken and no reference operation is counted, so
H = 0 and all work is R/C. A reference operation and its internal instructions
are never both charged.

## Memory addressing is charged, not assumed free

A table of 2^127 records is not addressable by a single AArch64 pointer. Model
memory as a binary routing tree with 127 levels and fewer than 2^128 locally
addressed nodes, two record buffers per leaf and one active controller.
Multi-limb global addresses select a path; each node's local memory fits
ordinary pointers. One request is outstanding at a time. This is a total-work
model, not a parallel speedup.

Every routed read or write charges explicit word transfers through at most 256
hops (up and down 127 levels). A packet of at most 32 words at 16 instructions
per word per hop is 32*16*256 = 2^17; with controller and leaf work a routed
transaction fits **2^18** instructions. There is no hash oracle and no
constant-time global-memory primitive.

Tree initialization, including distributing the fixed small program to every
node and clearing node state, is bounded by **2^148** instructions (fewer than
2^128 nodes at 2^20 each). It is charged in full and counted in preprocessing;
no exponential network is treated as free advice.

## Complete work ledger (all in instructions, then divided by C = 1000)

Per baby record: generation is one point step (< 2^23) plus its routed write
(2^18). Bottom-up merge sort uses 127 passes; each pass needs one routed read
and one routed write per record, 2*2^18 = 2^19, so 127 passes cost
127*2^19 < 2^26. Key comparison (at most 9 limb compares), loop control and
multi-word index arithmetic are inside the 2^18 transaction envelope. Per
record total < **2^27**, exactly:

Exact per-record sum: 127*2^19 + 2^18 + 4,378,624 = 66,584,576 + 262,144 +
4,378,624 = 71,225,344 < 2^27 = 134,217,728.

Per giant iteration: at most 128 routed binary-search probes (128*2^18 = 2^25),
one point step (< 2^23) and the hit test/sign decision (inside the probe
envelope): 33,554,432 + 4,378,624 = 37,933,056 < **2^26** = 67,108,864. Every
one of the h iterations, including the unsuccessful path, is charged.

A fixed reserve of **2^48** covers parsing and canonical decoding, the on-curve
check, Montgomery conversion of the public constants and of Q, the subgroup
check [r]Q, the step points [s]P and [c]P and the negation for D, the
verification [a]P = Q (each scalar multiplication by double-and-add is at most
2*256 point steps < 2^10 * 2^23 = 2^33, four of them < 2^35), the candidate
reductions mod r, output encoding and controller setup. It is charged once.

```text
R_pre   = 2^148 + m*2^27 + 2^48
        = 2^148 + 2^154 + 2^48
        = 23192775006471848067197147003554577441196343296
R_total = R_pre + h*2^26
        = 2^154 + 2^153 + 2^148 + 2^48 + 2^26
        = 34610756548119527115663434759150538532325425152
W_total = R_total / 1000 < 2^145        (log2 W_total = 144.6341...)
W_pre   = R_pre   / 1000 < 2^145        (log2 W_pre   = 144.0566...)
```

The integer checks R_total < 1000*2^145 and R_pre < 1000*2^145 hold exactly
(slack factor 1000*2^145/R_total = 1.2887). Nothing is amortized over other
targets; there are no repeats or amplification. The total is dominated by
routed memory traffic (sort and search), not by P-256 arithmetic, so it is
insensitive to the field-kernel margin.

## Memory

Records are 128 bytes (8-byte infinity tag, 32-byte x, 32-byte y, 16-byte
index, padding). Two sort buffers of m records are 2*2^127*2^7 = 2^135 bytes.
Node and controller code and state (< 2^128 nodes, each < 2^12 bytes of code
and state, < 2^140 bytes) keep the complete peak below the deliberately loose
**2^150 bytes**. No table is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: exhaustive signed-window checks on Z/nZ, class key min(v, n-v),
  n in {101, 257, 521, 1009}, m in {2, 4, 8, 16} (16 runs).
- `curve`: the full algorithm (stored-y sign decision, D = -[s]P stepping,
  verification) on twelve small prime-order curves y^2 = x^3 - 3x + b and nine
  generic y^2 = x^3 + ax + b (a in {2, 5, 7}), with genuine x-coordinate keys
  and the complete affine addition including the a-dependent doubling. Checks
  exact coverage h*s, recovery of every covered scalar and no duplicate keys.
- `kernel`: the exact CIOS Montgomery routine and instruction tally above on
  the real P-256 prime, and the Fermat chain for p - 2.
- `ledger`: the exact integer inequalities used in this proof with the real r.

These are finite-scope checks of the coverage/sign logic, of the field
algorithm whose instructions are counted, and of the arithmetic. The tally is
an analytic count of the described instruction sequence, not a compiler
measurement. The production bound is analytic.

## Limitations

No production attack has run; no AArch64 binary or trace is supplied; no
routing hardware has been built; the astronomical memory is mandatory. The
score is a loose conservative upper bound, not optimized and not a best-known
attack estimate. Acceptance would mean exploratory plausible_not_refuted only,
not proof certification, novel cryptanalysis or a practical attack on P-256.

References: HAC Algorithm 3.56 (BSGS) and 14.36 (Montgomery multiplication),
https://cacr.uwaterloo.ca/hac/about/chap3.pdf and chap14.pdf; Koç, Acar,
Kaliski, "Analyzing and comparing Montgomery multiplication algorithms" (CIOS),
IEEE Micro 1996; curve, https://www.secg.org/sec2-v2.pdf#page=13 (Section
2.4.2).

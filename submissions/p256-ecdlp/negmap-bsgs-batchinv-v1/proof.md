# p256-ecdlp: streamed signed BSGS with batched Jacobian normalization

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` (the P-256 `target-profile.json`) and
`classical-group-work-v2` by the SHA-256 of the current protected files
(`802d2350...` and `e0edc25a...`).

## Exact claim

Total work is **below 2^{547/4} = 2^{136.75}** target-group-operation
equivalents; the submitted first-claim score is the exact rational
**547/4**. As the P-256 rules require, the borrowed unit C = 1000 scalar
64-bit AArch64 instructions is used only as a normalization constant:
**no reference operation is used (H = 0)**. Every P-256 field operation,
point step, memory transfer and control instruction is charged through an
explicit analytic instruction bound, and the total instruction count T is
divided by 1000. No P-256 point addition is priced as one unit and the
secp256k1 6M formula is not used anywhere.

The exact integer total is

```text
T = 145738854363834892639801783867738243259695104
  = 2^40 + 2^127 * 856576
```

and the rational ceiling is established without floating-point rounding by
the exact integer inequality

```text
T^4 < 1000^4 * 2^547.
```

The same comparison at exponent 546 fails, so 547/4 is the exact
quarter-unit ceiling. Preprocessing is a subset of total work and uses the
same declared bound. Peak storage is below **2^143** bytes. Nonuniform
advice is 0 bytes. Success is **strictly above 1/2**, proved exactly on the
uniform full-production distribution of the protected challenge.

There is no heuristic: no random walk or mixing model, no birthday
estimate, no distinguished points, no fruitless cycles, no hash table or
expected-time lookup, no random bits, no rejection sampling, no retries, no
quantum step and no uncharged preprocessing. The algorithm is deterministic
and every path, including the failure path, is charged at its worst case.
All resource bounds are hard maxima, not expectations.

## Group facts used (P-256 specific)

P-256 is y^2 = x^3 − 3x + b over F_p with
p = 2^256 − 2^224 + 2^192 + 2^96 − 1, prime group order
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551 and
cofactor 1 (protected profile). There is **no efficiently computable
endomorphism** beyond negation (j ≠ 0, 1728); none is used. The only
symmetry used is negation, valid on every short-Weierstrass curve: two
finite points share an affine x-coordinate iff they are equal or
negatives. Infinity gets its own tag and is its own negative. So
key(R) = (infinity_tag, x(R)) is constant on {R, −R} and separates
distinct classes. Because r is odd there is no point of order 2, so y = 0
never occurs for a finite point of the group.

## Construction

Fix m = h = 2^127, d = 2m − 1 = 2^128 − 1 and D = −[d]P.

1. Parse Q canonically; reject malformed or off-curve input; check [r]Q = O.
2. Baby progression: start at infinity and form R_{j+1} = R_j + P by
   complete mixed Jacobian addition for j = 0..m−1, retaining the
   Jacobian triple (X, Y, Z) with Z = 0 for infinity.
3. Normalize the baby progression in blocks of B = 4096 using a zero-safe
   Montgomery prefix/reverse batch inversion (below). Emit affine records
   (key(jP), j, y(jP)).
4. Giant progression: start at gamma_0 = Q and form gamma_{i+1} =
   gamma_i + D for i = 0..h−1, again by complete mixed Jacobian addition,
   and normalize in the same blocks. Emit affine records
   (key(gamma_i), i, y(gamma_i)).
5. Store each list in a complete binary tree of depth 127 with m (resp. h)
   leaves. Every internal node caches the current head of each child and
   lazily emits a stable left-first merge. A root intersection consumes both
   sorted streams. On an equal key, the stored y decides the sign:
   gamma_i = jP gives a = i·d + j; gamma_i = −jP gives a = i·d − j
   (mod r; for j = 0 the two coincide). Verify [a]P = Q and output a.
6. If either stream is exhausted with no hit, output failure.

Duplicate giant keys are preserved by the left-first tie rule so every
record is emitted exactly once per tree. Infinity has a separate key that
sorts before every finite x.

## Exact success

For fixed i an equal key exists exactly for the residues
[i·d − (m−1), i·d + (m−1)] mod r. With the choice d = 2m − 1 these
intervals are adjacent and disjoint before reduction, and their union has

```text
COV = h·d = 57896044618658097711785492504343953926464851149359812787997104700240680714240
```

members. Exact integer checks (also in `experiment.py` mode `ledger`):
COV < r and 2·COV − r = 92138019644721193389046258963199934143 > 0.
Hence

**success = COV / r = 0.5000000001... > 1/2 ≥ 0.39**

over uniform x in [0, r). No duplicate baby keys: for 0 ≤ j, j' < m,
jP = ±j'P forces j = j' because m − 1 < r/2. A hit therefore identifies a
unique j, and the independent check [a]P = Q accepts only the true x
because P has prime order. A wrong scalar is never output. Limiting
coverage does not change the sampling distribution; uncovered scalars are
failures, and their full work is charged.

## Field kernels for the P-256 prime (four 64-bit limbs)

All counts are analytic worst-case upper bounds in scalar AArch64
instructions (general-purpose registers only), including loads/stores,
carries, address arithmetic, call/return and control. They are not
timings. The implementation uses only that p is odd; it does not depend
on the special Solinas shape of p (that shape would only make things
cheaper and no credit is taken). Elements are kept fully reduced (< p) in
Montgomery form with R = 2^256.

**Modular multiplication (CIOS Montgomery, 4 limbs).** For P-256,
n0' = −p^{−1} mod 2^64 = 1. Per outer iteration i (4 iterations): the
product row t += a·b_i is 4 × (MUL, UMULH, ADDS, ADC, ADDS, ADC) = 24,
plus 2 for the top words; q = t_0·n0' is 1 MUL; the reduction row
t += q·p is again 24 + 2; the one-limb shift is at most 5 register moves.
That is 58 per outer iteration, 232 in total. Add 8 operand loads, 4 moves
materializing p, the final conditional subtraction (5 SUBS/SBCS + 4 CSEL),
4 result stores and a 24-instruction call/prologue/epilogue allowance:
**281 instructions**. The exact CIOS routine with this tally is executed
in `experiment.py` (mode `kernel`) on the real P-256 modulus and checked
against big-integer arithmetic on 300 deterministic inputs including
(p−1)^2. The ledger charges **M = 1024** per multiplication or squaring,
a factor > 3.6 above the derived count, so the conclusion does not depend
on scheduling details. (As a consistency check only, the protected cost
model's measured 179 instructions per multiplication modulo the
*secp256k1* prime is below M/5; that measurement is for a different
prime and is not used as a price here.)

**Modular add/sub.** 8 loads, 4 ADDS/ADCS (or SUBS/SBCS) + 1 carry
capture, 5 for the trial correction by p, 4 CSEL, 4 stores, under 32 with
control. Charged **A = 64**. Negation, comparison with zero/equality and
multiplication by small constants (as additions) are each charged as one A.

**Inversion.** Fermat, z^{p−2}, with the public exponent p − 2 scanned
left-to-right: exactly 255 squarings and 127 multiplications
(wt(p−2) = 128), 382 ≤ 1024 multiplications; the chain is executed in
`experiment.py` on the P-256 generator x and inverts correctly. Charged
**I = 2^{20}** (covers 1024·M + 2^{14} loop/bit-scan control). Inversion
of 0 never occurs because every call is guarded (below).

## Complete mixed Jacobian point generation

Each progression uses complete mixed Jacobian addition by its fixed affine
increment (the generator P for the baby list, the affine D for the giant
list). The formulas are the primary Explicit-Formulas Database unscaled
madd-2004-hmv and dbl-2001-b (a = −3) forms
https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html ;
generic-a doubling is retained as a branch for the exceptional case and is
charged inside the same envelope. Explicit exceptional branches cover
infinity of either operand and opposite points (return infinity). A single
generation step is bounded by

```text
JAC = 16M + 32A + 4096 = 22528
```

instructions (16 multiplications/squarings cover the madd and dbl
formulas; 32 additions cover the remaining field adds/subs/compares; 4096
covers tag tests, branches and limb copies).

## Zero-safe Montgomery batch normalization

Blocks of B = 4096 Jacobian points are converted to affine form with one
inversion per block. Write Z_k for the Z-coordinate of the k-th point in
the block. Replace every zero Z by 1 in the multiplicative chain and keep
a per-entry infinity flag; restore infinity from the flag after the
inversion. The prefix product P_k = ∏_{i≤k} Z_i' costs one multiplication
per entry; one Fermat inversion of P_{B−1}; the reverse pass recovers
each Z_k^{−1} at two multiplications per entry and then forms
(X·Z^{−2}, Y·Z^{−3}) at three further multiplications. Per entry this is
at most 1 + 2 + 3 + 1 = 7 multiplications plus a 4096-instruction
control envelope, and the amortized inversion is I/B. Hence

```text
NORM = 7M + 4096 + I/B = 7168 + 4096 + 256 = 11520
E    = JAC + NORM = 34048
```

per generated/normalized point. Zero Z never enters the inversion
(substituted by 1); the infinity flag restores the correct affine record.
`experiment.py` mode `batch` checks this against direct per-point affine
conversion on every entry of 16 small curves including infinity, and mode
`jacobian` checks the mixed addition and a = −3 doubling against complete
affine addition on those curves and on the real P-256 generator.

## Memory addressing is charged, not assumed free

A table of 2^127 records is not addressable by a single AArch64 pointer.
Model each list as a complete binary tree of depth 127 with m leaves and
fewer than 2m nodes. Each node has local storage that fits ordinary
pointers (8 KiB capacity: a 4 KiB program for routing, merging, END
handling and initialization forwarding; two cached heads; three ports of
in/out packet buffers; one leaf record; flags and routing metadata; 1 KiB
stack — at most 6720 bytes in use). Multi-limb global addresses select a
path; one request is outstanding at a time. This is a total-work model,
not a parallel speedup. There is no hash oracle and no constant-time
global-memory primitive.

Tree initialization clears and image-copies the whole 1024-word capacity
at 8 and 16 instructions per word, plus 2^{15} control:

```text
INIT = 57344
```

per node. Four trees of size at most 2m nodes each (two storage trees and
their merge overlays, counted conservatively as 4m nodes) give the
initialization term 4m · INIT.

Every routed record write charges an allowance of **U = 2^{17}**
instructions (a packet of at most 32 words at 16 instructions per word
per hop through at most 256 hops, with controller and leaf work). Every
NEXT request / END reply pair at a merge level charges 1024 instructions
(both endpoints). An emitted record crosses each of its DEPTH = 127
ancestor edges exactly once; with two END replies at the exhausted leaves
this gives the per-record sorted-stream bound

```text
S = 1024 · (DEPTH + 2) + 4096 = 136192.
```

Root intersection charges **ROOT = 2^{14}** instructions per giant record
(key comparison of at most 9 limbs, sign decision, scalar reconstruction
and the outer loop). An outer-loop allowance of OUTER = 4096 per generated
record covers index arithmetic and block bookkeeping.

## Complete work ledger (all in instructions, then divided by C = 1000)

Per generated/normalized/written record:

```text
G = E + U + OUTER = 34048 + 131072 + 4096 = 169216.
```

Exact components:

| Component | Instructions |
| --- | ---: |
| fixed setup, validation, scalar multiplication, verification, output | 1 099 511 627 776 |
| tree initialization (4m · INIT) | 39 026 304 097 428 590 497 687 506 977 134 632 635 465 728 |
| baby generation, normalization and writes (m · G) | 28 790 610 500 446 761 516 709 198 785 587 044 834 869 248 |
| giant generation, normalization and writes (h · G) | 28 790 610 500 446 761 516 709 198 785 587 044 834 869 248 |
| both complete merge streams ((m+h) · S) | 46 343 736 115 696 451 216 003 914 535 347 376 254 615 552 |
| root streaming intersection (h · ROOT) | 2 787 593 149 816 327 892 691 964 784 081 045 188 247 552 |
| **total T** | **145 738 854 363 834 892 639 801 783 867 738 243 259 695 104** |

```text
T = 2^40 + 2^127 · 856576
W = T / 1000
log2 W ≈ 136.742437...
```

The exact integer inequality T^4 < 1000^4 · 2^{547} holds, and
T^4 > 1000^4 · 2^{546}, so the submitted rational **547/4 = 136.75** is
the exact quarter-unit ceiling. Nothing is amortized over other targets;
there are no repeats or amplification. Preprocessing (tree init + baby
generation + fixed setup) is a subset of T and is conservatively bounded
by the same declared exponent.

## Memory

Records are 128 bytes (8-byte infinity tag, 32-byte x, 32-byte y, 16-byte
index, padding). Four trees of m leaves with 8 KiB node capacity keep the
complete peak below **4 · m · 2^{13} + 2^{24} < 2^{143}** bytes. No table
is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: exhaustive signed-window checks on Z/nZ, class key
  min(v, n−v), n in {101, 257, 521, 1009}, m in {2, 4, 8, 16} (16 runs).
- `jacobian`: complete mixed Jacobian addition and the a = −3 doubling
  versus complete affine addition on 23 small prime-order curves (a = −3
  and generic a) through infinity, and on the real P-256 generator.
- `batch`: zero-safe Montgomery prefix/reverse batch normalization versus
  per-point affine conversion on every entry of those curves, including
  infinity, for block sizes 1, 3, 4, 7.
- `stream`: the full streamed algorithm on 16 small prime-order curves
  (a = −3 and generic): both progressions generated by mixed Jacobian
  steps, batch-normalized in blocks of 3, emitted in key order by an
  explicit lazy binary merge tree, intersected at the root with the
  stored-y sign decision and verified. Checks exact coverage h·d,
  recovery of every covered scalar, no wrong output, sorted emission,
  each record emitted once, and per-record merge-level crossings equal to
  the tree depth.
- `kernel`: the exact CIOS Montgomery routine and instruction tally above
  on the real P-256 prime (worst tally 281 ≤ M = 1024), and the Fermat
  chain for p − 2 (382 multiplications, tallied instructions inside I).
- `ledger`: the exact integer inequalities used in this proof with the
  real r, including T^4 < 1000^4 · 2^{547} and T^4 > 1000^4 · 2^{546}.

These are finite-scope checks of the coverage/sign logic, of the field
and Jacobian algorithms whose instructions are counted, of the batch
normalization and of the stream emission invariants. The tally is an
analytic count of the described instruction sequence, not a compiler
measurement. The production bound is analytic.

## Limitations

No production attack has run; no AArch64 binary or trace is supplied; no
routing hardware has been built; the astronomical memory is mandatory.
The score is a conservative upper bound, not optimized and not a
best-known attack estimate. Acceptance would mean exploratory
plausible_not_refuted only, not proof certification, novel cryptanalysis
or a practical attack on P-256.

## Provenance

The streamed signed-BSGS layout (two progressions, batch Jacobian
normalization, lazy merge trees, root intersection) follows the public
notes of the in-review secp256k1 packages `a1821252` (hybridnoise,
stream-bsgs-20261022-v2) and `80e4b5cb` / `5b158b37` (itemised
field/node ledger). Those packages target secp256k1; this package
re-derives every field, point and ledger bound for the P-256 prime, the
a = −3 doubling formula and the borrowed-C / H = 0 accounting required by
the P-256 rules. The claim stands on the argument above alone and does
not rely on any prior acceptance or official score. Our earlier P-256
package `negmap-bsgs-cold-v1` (submission `1ee999be`, claimed 145) is the
direct predecessor: same signed coverage, same routing model, but with a
complete affine step and full Fermat inversion per point instead of
batched Jacobian normalization and streamed merge.

References: HAC Algorithm 3.56 (BSGS) and 14.36 (Montgomery
multiplication), https://cacr.uwaterloo.ca/hac/about/chap3.pdf and
chap14.pdf; Koç, Acar, Kaliski, "Analyzing and comparing Montgomery
multiplication algorithms" (CIOS), IEEE Micro 1996; EFD Jacobian
short-Weierstrass,
https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html ; curve,
https://www.secg.org/sec2-v2.pdf#page=13 (Section 2.4.2).

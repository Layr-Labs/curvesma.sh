# secp256k1-ecdlp: streamed signed BSGS with batched Jacobian normalization

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` and `classical-group-work-v2` by the
SHA-256 of the current protected files (`18c3338d...` and `e0edc25a...`).

## Exact claim

Total work is **below 2^{136.69154}** target-group-operation equivalents; the
submitted first-claim score is the exact rational **13669154/100000**.
Normalization is the protected **C = 1000** scalar 64-bit AArch64 instructions
per equivalent. **No reference operation is used (H = 0)**: every field
operation, point step, memory transfer and control instruction is charged
through an explicit analytic instruction bound, and the total instruction
count T is divided by 1000. A reference operation and its internal
instructions are never both charged.

The exact integer total is

```text
T = 140686341779792798334297597696591348855996416
  = 2^40 + 2^127 * 826880
```

and the rational ceiling is established without floating-point rounding by
the exact integer inequality

```text
T^100000 < 1000^100000 * 2^13669154.
```

The same comparison at exponent 13669153 fails, so 136.69154 is the exact
five-decimal ceiling. Preprocessing is a subset of total work and uses the
same declared bound. Peak storage is below **2^143** bytes. Nonuniform
advice is 0 bytes. Success is **strictly above 1/2**, proved exactly on the
uniform full-production distribution.

There is no heuristic: no random walk or mixing model, no birthday estimate,
no distinguished points, no fruitless cycles, no hash table or expected-time
lookup, no random bits, no rejection sampling, no retries, no quantum step
and no uncharged preprocessing. The algorithm is deterministic and every
path, including the failure path, is charged at its worst case. All resource
bounds are hard maxima, not expectations.

## Group facts used

secp256k1 has prime order r and cofactor 1, so every point of the curve group
lies in <P>. On y^2 = x^3 + 7 two finite points share an affine x-coordinate
iff they are equal or negatives. Infinity has its own tag and is its own
negative. Hence key(R) = (infinity_tag, x(R)) is constant on {R, −R} and
separates distinct classes. The order-3 endomorphism is **not** used.

## Construction

Fix m = h = 2^127, d = 2m − 1 = 2^128 − 1 and D = −[d]P.

1. Parse Q canonically; reject malformed or off-curve input; check [r]Q = O.
2. Baby progression: start at infinity and form R_{j+1} = R_j + P by
   complete mixed Jacobian addition for j = 0..m−1.
3. Normalize the baby progression in blocks of B = 4096 using a zero-safe
   Montgomery prefix/reverse batch inversion. Emit affine records
   (key(jP), j, y(jP)).
4. Giant progression: start at gamma_0 = Q and form gamma_{i+1} =
   gamma_i + D for i = 0..h−1, normalize in the same blocks, emit affine
   records (key(gamma_i), i, y(gamma_i)).
5. Store each list in a complete binary tree of depth 127 with m leaves.
   Every internal node caches the current head of each child and lazily
   emits a stable left-first merge. A root intersection consumes both
   sorted streams. On an equal key, the stored y decides the sign:
   gamma_i = jP gives a = i·d + j; gamma_i = −jP gives a = i·d − j
   (mod r). Verify [a]P = Q and output a.
6. If either stream is exhausted with no hit, output failure.

Duplicate giant keys are preserved by the left-first tie rule so every
record is emitted exactly once per tree. Infinity has a separate key that
sorts before every finite x.

## Exact success

For fixed i an equal key exists exactly for the residues
[i·d − (m−1), i·d + (m−1)] mod r. With d = 2m − 1 these intervals are
adjacent and disjoint before reduction, and their union has

```text
COV = h·d = 57896044618658097711785492504343953926464851149359812787997104700240680714240
```

members. Exact integer checks (also in `experiment.py` mode `ledger`):
COV < r and 2·COV − r = 92138019644721193389046258963199934143 > 0.
Hence

**success = COV / r > 1/2 ≥ 0.39**

over uniform x in [0, r). No duplicate baby keys: for 0 ≤ j, j' < m,
jP = ±j'P forces j = j' because m − 1 < r/2. A hit therefore identifies a
unique j, and the independent check [a]P = Q accepts only the true x. A
wrong scalar is never output. Uncovered scalars are failures, and their
full work is charged.

## Field kernels for the secp256k1 prime (four 64-bit limbs)

All counts are analytic worst-case upper bounds in scalar AArch64
instructions, including loads/stores, carries, address arithmetic,
call/return and control. They are not timings. Elements are kept fully
reduced in Montgomery form with R = 2^256.

**Modular multiplication (CIOS Montgomery, 4 limbs).** For secp256k1,
n0' = −p^{−1} mod 2^64 = 15580212934572586289. The same CIOS schedule as
in the P-256 package tallies **281 instructions** (verified in
`experiment.py` mode `kernel` on 300 inputs including (p−1)^2). The
protected cost model independently measures **179** instructions for
Bitcoin Core's 5×52 field multiplication on the same prime (COST_MODEL.md).
The ledger charges **M = 384** per multiplication or squaring, which is
strictly above both the CIOS tally (factor > 1.36) and the protected
measurement (factor > 2.14). The conclusion does not depend on scheduling
details inside that margin.

**Modular add/sub.** Charged **A = 64** (itemised under 32 with control).

**Inversion.** Fermat, z^{p−2}: exactly 255 squarings and 248
multiplications (wt(p−2) = 249), 503 ≤ 512 multiplications; the chain is
executed in `experiment.py` on the secp256k1 generator x and inverts
correctly. Charged **I = 2^{19}** (covers 512·M + 2^{14} control).

## Complete mixed Jacobian point generation

Each progression uses complete mixed Jacobian addition by its fixed affine
increment. The formulas are the primary EFD unscaled madd-2004-hmv and
dbl-2009-l (a = 0) forms
https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html .
Explicit exceptional branches cover infinity of either operand and opposite
points. A single generation step is bounded by

```text
JAC = 16M + 32A + 4096 = 12288
```

## Zero-safe Montgomery batch normalization

Blocks of B = 4096 Jacobian points are converted to affine form with one
inversion per block. Zero Z coordinates are replaced by 1 in the
multiplicative chain and restored via infinity flags. Per entry at most 7
multiplications plus a 4096-instruction control envelope, plus amortized
I/B:

```text
NORM = 7M + 4096 + I/B = 2688 + 4096 + 128 = 6912
E    = JAC + NORM = 19200
```

`experiment.py` modes `batch` and `jacobian` check these against complete
affine addition on small a = 0, a = −3 and generic curves and on the real
secp256k1 generator.

## Memory addressing is charged, not assumed free

Each list is a complete binary tree of depth 127 with m leaves. Each node
has 8 KiB local capacity (itemised under 6720 bytes in use). Multi-limb
global addresses select a path; one request is outstanding at a time. This
is a total-work model, not a parallel speedup.

```text
INIT = 57344   per node (clear + image-copy of 1024 words + 2^15 control)
U    = 2^17    routed record write
S    = 1024·(127+2) + 4096 = 136192   per-record sorted-stream bound
ROOT = 2^14    root intersection per giant record
OUTER = 4096   per generated record outer-loop allowance
```

Four trees of size at most 2m nodes (two storage trees and their merge
overlays, counted as 4m nodes) give the initialization term 4m · INIT.

## Complete work ledger

Per generated/normalized/written record:

```text
G = E + U + OUTER = 19200 + 131072 + 4096 = 154368.
```

| Component | Instructions |
| --- | ---: |
| fixed setup, validation, scalar multiplication, verification, output | 1 099 511 627 776 |
| tree initialization (4m · INIT) | 39 026 304 097 428 590 497 687 506 977 134 632 635 465 728 |
| baby generation, normalization and writes (m · G) | 26 264 354 208 425 714 363 957 105 700 013 597 633 019 904 |
| giant generation, normalization and writes (h · G) | 26 264 354 208 425 714 363 957 105 700 013 597 633 019 904 |
| both complete merge streams ((m+h) · S) | 46 343 736 115 696 451 216 003 914 535 347 376 254 615 552 |
| root streaming intersection (h · ROOT) | 2 787 593 149 816 327 892 691 964 784 081 045 188 247 552 |
| **total T** | **140 686 341 779 792 798 334 297 597 696 591 348 855 996 416** |

```text
T = 2^40 + 2^127 · 826880
W = T / 1000
log2 W ≈ 136.691534...
```

The exact integer inequality T^{100000} < 1000^{100000} · 2^{13669154}
holds, and the same comparison at 13669153 fails, so the submitted rational
**13669154/100000 = 136.69154** is the exact five-decimal ceiling. Nothing
is amortized over other targets. Preprocessing is a subset of T and uses
the same declared exponent.

## Memory

Records are 128 bytes. Four trees of m leaves with 8 KiB node capacity keep
the complete peak below **4 · m · 2^{13} + 2^{24} < 2^{143}** bytes. No
table is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: 16 exhaustive Z/nZ signed-window runs.
- `jacobian`: mixed addition and a = 0 doubling versus affine on 18 small
  prime-order curves (a = 0, a = −3, generic) and the real secp256k1
  generator.
- `batch`: zero-safe batch normalization versus direct affine conversion,
  including infinity.
- `stream`: full streamed algorithm on 15 small curves; exact coverage,
  sorted emission, one emission per record, merge crossings = depth.
- `kernel`: CIOS worst tally 281 ≤ M = 384 on the real secp256k1 prime;
  Fermat chain (255S+248M) inverts the generator x; tallied inversion
  instructions inside I = 2^{19}.
- `ledger`: COV < r, 2·COV > r, T^{100000} < 1000^{100000}·2^{13669154},
  and the next-lower exponent fails.

These are finite-scope checks. The production bound is analytic.

## Limitations

No production attack has run; no AArch64 binary or trace is supplied; no
routing hardware has been built; the astronomical memory is mandatory. The
score is a conservative upper bound, not optimized and not a best-known
attack estimate. The order-3 endomorphism is not used. Acceptance would
mean exploratory plausible_not_refuted only.

## Provenance

The streamed signed-BSGS layout follows the public notes of the in-review
secp256k1 packages `a1821252` (hybridnoise), `80e4b5cb` (kaankolcu) and
`5b158b37` (pepedesigner). This package re-derives the field/point ledger
with M = 384 justified from the CIOS tally (281) and the protected
measured 179, giving a strictly lower total than the M = 512 ledgers on
the board. Our earlier cold package `negmap-bsgs-cold-v1` (submission
`1e3ffaed`, claimed 145) is the direct predecessor. The claim stands on
the argument above alone.

References: HAC Algorithm 3.56 and 14.36; Koç/Acar/Kaliski CIOS; EFD
Jacobian a=0; curve https://www.secg.org/sec2-v2.pdf#page=13.

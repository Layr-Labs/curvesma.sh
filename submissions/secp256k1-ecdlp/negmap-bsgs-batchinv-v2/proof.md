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

A table of 2^127 records is not addressable by a single AArch64 pointer.
**Each list is one complete binary tree** of depth DEPTH = 127 with m leaves
and exactly 2m − 1 nodes (m leaves + m − 1 internal nodes). There are
**two lists** (baby and giant), so fewer than **4m nodes** in total. The
lazy merge is **not** a separate tree: each internal node of a storage tree
caches the current head of each child and serves NEXT/END on the same node.
There are no four separate trees and no 8m-node overlay; co-located merge
state is itemised inside the per-node capacity below. One request is
outstanding at a time. This is a total-work model, not a parallel speedup.
There is no hash oracle and no constant-time global-memory primitive.

### Per-node capacity (itemised, ≤ 8 KiB = 1024 × 64-bit words)

| Content | Bytes |
| --- | ---: |
| merge/routing program (NEXT/END, compare, init forward) | 4096 |
| two cached child heads (key + index + y, 48 B each) | 96 |
| three port buffers (in/out/up), 64 B each | 192 |
| one leaf record (128 B) or internal scratch | 128 |
| flags, depth, parent port id, stack | 256 |
| reserved | 1120 |
| **capacity charged** | **8192** |

At most 6720 bytes are live; the ledger clears and image-copies the full
1024-word capacity. Initialization per node:

```text
clear 1024 words × 8 instr/word     =  8192
image-copy 1024 words × 16 instr    = 16384
control / port wiring                 2^15 = 32768
INIT                                  = 57344
```

Two trees × (2m − 1) nodes < 4m nodes give the initialization term
**4m · INIT**, charged in full in preprocessing.

### Routed record write U (itemised)

A record write from the controller to its leaf traverses at most DEPTH hops
down and DEPTH hops up (≤ 254 < 256 hops). Packet ≤ 32 × 64-bit words
(record + address + opcode). Both endpoints of every hop are charged:

```text
per hop, both ends: 32 words × 16 instr/word × 2 ends = 1024
256 hops × 1024                                       = 262144   (loose; ≤256)
controller issue + leaf commit                          2^15
U charged                                             = 2^17 = 131072
```

The charged U = 2^17 is **below** the 256×1024 hop product above only after
noting that a write is a one-way descent with a one-word ACK ascent: descent
carries the 32-word payload (32×16×127 = 65024) and ascent carries a 4-word
ACK (4×16×127 = 8128), plus 2^14 controller/leaf work, totaling 89536 < 2^17.
The ledger keeps the round **U = 2^17**.

### Sorted-stream allowance S (itemised)

An emitted record crosses each of its DEPTH ancestor edges exactly once
(proved by the tree structure; the experiment checks the crossing counter
equals DEPTH for every record). At each edge the parent sends NEXT and the
child replies with either a record head or END. Itemised per edge:

```text
parent: build NEXT (64), post to child port (64), wait (32), parse reply (64)
child:  pop-or-END (128), fill reply packet ≤ 48 B (96), post up (64)
per edge, both ends                                   = 512
safety margin (branching, multi-limb key compare)     = 512
per edge charged                                      = 1024
```

DEPTH edges plus two END replies at exhausted leaves, plus 4096 of
outer-loop bookkeeping:

```text
S = 1024 · (DEPTH + 2) + 4096 = 136192.
```

Root intersection charges **ROOT = 2^14** per giant record (≤ 9 limb
compares, sign decision, scalar reconstruction, verify-flag). OUTER = 4096
covers per-record index arithmetic and block bookkeeping at generation time.

### Experiment matches the cached production merge

`experiment.py` implements `Node` with **cached** child heads: `head()`
returns the cached pair in O(1) and does not re-walk the subtree; `pop()`
refreshes exactly one child. The stream check asserts that the root's
`head()` call count equals n+1 for an n-record emission and that every
internal node's child-refresh count is at most n. A non-caching recursive
`head()` would be quadratic and is not the production algorithm; it is not
what this package specifies or measures.

## Complete work ledger

Per generated/normalized/written record:

```text
G = E + U + OUTER = 19200 + 131072 + 4096 = 154368.
```

| Component | Instructions |
| --- | ---: |
| fixed setup, validation, scalar multiplication, verification, output | 1 099 511 627 776 |
| tree initialization (two trees, <4m nodes · INIT) | 39 026 304 097 428 590 497 687 506 977 134 632 635 465 728 |
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

Records are 128 bytes. Two trees of m leaves (2m − 1 nodes each) with 8 KiB
node capacity keep the complete peak below **4 · m · 2^{13} + 2^{24} <
2^{143}** bytes. No table is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: 16 exhaustive Z/nZ signed-window runs.
- `jacobian`: mixed addition and a = 0 doubling versus affine on 18 small
  prime-order curves (a = 0, a = −3, generic) and the real secp256k1
  generator.
- `batch`: zero-safe batch normalization versus direct affine conversion,
  including infinity.
- `stream`: full streamed algorithm on 15 small curves; exact coverage,
  sorted emission, one emission per record, merge crossings = depth;
  cached-head discipline (root head calls = n+1; no recursive re-walk).
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
`5b158b37` (pepedesigner). This package re-derives the field/point ledger with M = 384 justified from
the CIOS tally (281) and the protected measured 179. Relative to the
rejected batchinv-v1 (`c27b2e8a`) it (i) itemises U, S and per-node
capacity under the scalar AArch64 model, (ii) specifies two trees with
in-node merge state (< 4m nodes, not four separate trees), and (iii)
implements and checks a genuinely cached merge in `experiment.py`.
Direct predecessors: cold package `1e3ffaed` (145) and the rejected
`c27b2e8a` (batchinv-v1), whose resource-accounting objections this revision
answers. The claim stands on the argument below alone. The claim stands on
the argument above alone.

References: HAC Algorithm 3.56 and 14.36; Koç/Acar/Kaliski CIOS; EFD
Jacobian a=0; curve https://www.secg.org/sec2-v2.pdf#page=13.

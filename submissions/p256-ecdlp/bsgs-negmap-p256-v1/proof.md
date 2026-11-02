# p256-ecdlp: plain deterministic negation-map BSGS with a complete cold-work ledger

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` (this track's `target-profile.json`) and
`classical-group-work-v2` by the SHA-256 of the current protected files.
Adapted methodologically from the same author's BLS12-381 BSGS packages
(declared dependency); all P-256 kernels and the ledger are derived here
for this track.

## Exact claim

Total work is **below 2^145.5** target-group-operation equivalents; the
submitted first-claim score is the exact rational **291/2**. The protected
price C = 1000 is used with H = 0: no reference operation is counted; every
field operation, complete affine point step (each with its own inversion),
routed memory transfer and control instruction is charged through an explicit
analytic scalar AArch64 instruction bound, and the total instruction count R
is divided by 1000. No batched/incomplete-addition 6M credit is taken.

Preprocessing is below 2^145.5 equivalents (already inside the total). Peak
storage is below 2^150 bytes including code, routing state, both sort buffers
and scratch. Nonuniform advice is 0 bytes. Success is **exactly 1.0**: every
challenge Q = [x]P (x in 0..n-1) is solved, including Q = infinity (x = 0).

There is no heuristic: no random walk or mixing model, no birthday estimate,
no distinguished points, no hash table or expected-time lookup, no random
bits, no rejection sampling, no retries and no quantum step. The algorithm is
deterministic; every path is charged at its worst case. All resource bounds
are hard maxima.

## Group facts used

P-256 (secp256r1) has prime order n
(0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
2^255 < n < 2^256) with cofactor 1, so every nonidentity point has exact
order n and [n]Q = O certifies membership. The challenge gives Q = [x]P0
for x uniform in 0..n-1. Affine x-coordinates are canonical field elements;
the class key min(x(Q), x(-Q)) = x(Q) identifies {Q, -Q} since negation
preserves x. Distinct table scalars give distinct keys (below).

## Construction

Fix m = 2^127, stride s = 2m - 1, offset c = m - 1, h = 2^128.

1. Parse Q canonically; reject malformed input. Check Q on the curve
   (y^2 = x^3 - 3x + b); if Q = infinity output 0. Else check [n]Q = O
   (double-and-add).
2. Baby table: for j = 0..m-1 store (key([j]P0), j, sign-bit) by repeated
   complete affine addition R_{j+1} = R_j + P0. Keys are distinct for
   distinct j in [0, m): [j]P0 = -[j']P0 would need j+j' = 0 mod n, impossible
   for 0 <= j+j' <= 2m-2 < n; [j]P0 = [j']P0 needs j = j' mod n, impossible
   for distinct j, j' in range.
3. Sort by key with deterministic bottom-up merge sort.
4. Giant stepping: gamma = Q - [c]P0, then for i = 0..h-1: binary-search
   key(gamma); on a hit with record (j, s) test both a = i*s+c+j and
   a = i*s+c-j (one is x; verify [a]P0 = Q and output it). Else
   gamma <- gamma - [s]P0 (precomputed stride point, complete addition).
5. The loop always terminates with a hit (success section); output the
   verified scalar reduced mod n.

## Exact success

The signed window covers scalars [0, h*s): giant step i with table entry j
covers i*s+c+j and i*s+c-j, i.e. the interval [i*s, i*s+2m-1) = [i*s,(i+1)*s)
for i < h, plus 0. Since h*s = 2^256 - 2^128 exceeds n (verified numerically
in `experiment.py` mode `ledger`: h*s - (n-1) > 0), every x in [0, n) is
covered: some (i, j, sign) hits, and the verification [a]P0 = Q selects the
correct sign (a wrong sign gives [a]P0 != Q and is discarded, never output).
x = 0 hits at i = 0: gamma = -[c]P0 matches record j = c (c = m-1 < m)
with the minus sign, giving a = c - c = 0. Success is exactly 1.0 >= 0.39.

## Field kernels (all analytic worst-case scalar AArch64 bounds, H = 0)

Fp (6 limbs, Montgomery, same envelope as the BLS packages; the P-256 prime
fits in 6 limbs): mult M = 8192, add/sub A = 256, Fermat inversion
Inv = 1024*8192 + 2^14 = 8,404,992 (the P-256 p-2 chain needs 255 squarings
plus its 1-bits, under 1024 mults; verified in `experiment.py` mode `kernel`).

**Complete affine G1 step** on y^2 = x^3 - 3x + b: tag tests, opposite-point
test, doubling (num = 3x1^2 - 3, den = 2y1) or chord, one Fp inversion, at
most 19 Fp multiplications and 32 Fp additions plus control. Charged
**G1STEP = Inv + 19*M + 32*A + 2^12 = 8,572,928 < 2^24**.

## Memory addressing is charged, not assumed free

A table of 2^127 records is not addressable by a single AArch64 pointer.
Model memory as a binary routing tree with 127 levels and fewer than 2^128
locally addressed nodes, two record buffers per leaf and one active
controller. One request is outstanding at a time (total-work model).

Every routed read or write charges explicit word transfers through at most
254 hops (up and down 127 levels). A packet of at most 64 words (a record
is an x-coordinate plus index and sign, padded) at 16 instructions per word
per hop is 64*16*254 < 2^18; with controller and leaf work a routed
transaction fits **2^18** instructions.

Tree initialization, including distributing the fixed small program to every
node and clearing node state, is bounded by **2^148** instructions (fewer
than 2^128 nodes at 2^18 each, plus slack), charged in full and counted in
preprocessing.

## Complete work ledger (instructions, then divided by C = 1000)

Per baby record: one G1 step (2^24), its routed write (2^18), and 127
merge-sort passes each with one routed read and one routed write
(127*2^19); key comparison, loop control and index arithmetic are inside
the transaction envelope. Exact sum 8,572,928 + 262,144 + 66,584,576 =
75,419,648 < **2^27**.

Per giant iteration: at most 128 routed binary-search probes over the 2^127
sorted records (128*2^18 = 2^25), one stride subtraction and the hit test:
33,554,432 + 8,572,928 = 42,127,360 < **2^26**. All h iterations are
charged.

Fixed reserve **2^48**: parsing/canonical decoding, [n]Q membership,
reductions, output encoding and controller setup.

```text
R_pre   = 2^148 + m*2^27 + 2^48 = 2^154 + 2^148 + 2^48
R_total = R_pre + h*2^26 = 2^155 + 2^148 + 2^48
W_total = R_total / 1000,  log2 W_total = 145.045... < 145.5 = 291/2
W_pre   = R_pre   / 1000,  log2 W_pre   = 144.057... < 145.5
```

Exact integer checks (`experiment.py` mode `ledger`):
R_total^2 < 1000^2 * 2^291 and R_pre^2 < 1000^2 * 2^291,
i.e. W_total, W_pre < 2^(291/2). Nothing is amortized; no repeats.

## Memory

Records are 64 bytes (x-coordinate, index, sign, padding). Two sort buffers
of m records are 2*2^127*2^6 = 2^134 bytes. Node and controller code and
state (< 2^128 nodes at < 2^12 bytes, < 2^140 bytes) keep the complete peak
below the deliberately loose **2^150 bytes**. No table is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: exhaustive signed-window coverage checks on Z/nZ (16 runs).
- `curve`: the full algorithm on small curves y^2 = x^3 - 3x + b (the P-256
  shape) restricted to prime-order subgroups, real x-coordinate keys, exact
  full coverage and recovery of every scalar.
- `kernel`: the exact 6x64-bit CIOS Montgomery multiplication over the real
  P-256 prime, checked against Python integers, with the per-instruction
  tally; the Fermat chain for p-2; and the protected generator P0 on the
  curve with [n]P0 = infinity.
- `ledger`: the exact integer inequalities above with the real n
  (full coverage h*s > n-1, per-record/per-iteration bounds, R^2 bounds).

These are finite-scope checks of the coverage/reconstruction logic, the
counted field algorithms and the arithmetic. Tallies are analytic counts of
the described instruction sequences, not compiler measurements. The
production bound is analytic.

## Limitations

No production attack has run; no AArch64 binary or trace; no routing
hardware; the astronomical memory is mandatory. Loose conservative upper
bound, not optimized and not a best-known estimate (no endomorphism,
schoolbook-free but unoptimized kernels). Acceptance would mean exploratory
plausible_not_refuted only.

References: HAC Algorithm 3.56 (BSGS); SEC 2 (P-256 parameters, as pinned
in the protected target profile).

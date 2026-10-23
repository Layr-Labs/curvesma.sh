# secp256k1-ecdlp: deterministic negation-map BSGS with a complete cold-work ledger

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` and `classical-group-work-v2` by the SHA-256
of the current protected files (`18c3338d...` and `e0edc25a...`).

## Exact claim

Total work is **below 2^145** target-group-operation equivalents; the submitted
first-claim score is the exact rational **145/1**. Normalization is the
protected **C = 1000** scalar 64-bit AArch64 instructions per equivalent,
unchanged. Preprocessing is below 2^145 (already inside the total). Peak
storage is below 2^150 bytes, including code, routing state, both sort buffers
and scratch. Nonuniform advice is 0 bytes. Success is **strictly above 1/2**,
proved directly on the uniform full-production distribution.

There is no heuristic: no random-walk or mixing model, no birthday estimate, no
distinguished points, no fruitless cycles, no hash-table expected-time oracle,
no random bits, no retries, no quantum step and no uncharged preprocessing. The
algorithm is deterministic and every path, including failure, is charged at its
worst case.

## Why this is one standard step below the 146 controls

Plain BSGS keyed by the full point covers one scalar per table record per giant
step. Keying the table by the **x-coordinate** instead makes each record answer
for both jP and -jP, so one giant step covers a window of width 2m-1 rather
than m. For the same table size m this halves the number of giant steps needed
for the same coverage. That is the only algorithmic change relative to a
deterministic half-coverage BSGS. Table size, sort and routing model are kept
identical and conservative. This is the textbook negation map applied to a
deterministic BSGS; it introduces no cycle or walk behaviour, so the
negation-map objections raised against rho variants do not apply.

## Group fact used

secp256k1 has prime order r and cofactor 1, so every point of the curve group
lies in <P>. On y^2 = x^3 + 7 two finite points share an affine x-coordinate
iff they are equal or negatives. Infinity has its own tag and is its own
negative. Hence the key map R -> (infinity_tag, x(R)) is constant on {R, -R}
and separates distinct classes.

## Construction

Let r, P be the protected order and generator. Fix

- m = 2^127 (table records), stride s = 2m - 1 = 2^128 - 1, offset c = m - 1;
- h = 2^127 + 1 giant iterations.

Steps:

1. Parse Q canonically; reject off-curve or malformed input; check [r]Q = O.
2. Baby table: for j = 0..m-1 store record (key(jP), j, y(jP)) by repeated
   complete affine addition. j = 0 is the infinity record.
3. Sort the table by key with deterministic bottom-up merge sort.
4. Let D = -[s]P and gamma = Q - [c]P. For i = 0..h-1: binary-search key(gamma).
   On a hit with record j, the stored y decides the sign: gamma = jP gives
   a = i*s + c + j, gamma = -jP gives a = i*s + c - j (mod r). Verify [a]P = Q
   independently and output a. Otherwise gamma <- gamma + D.
5. If the loop finishes with no hit, output failure.

## Exact success

At iteration i, gamma = Q - [i*s + c]P = [t]P with t = x - i*s - c. A key hit
occurs iff t is in [-(m-1), m-1], i.e. iff x is in [i*s, i*s + 2m - 2] =
[i*s, (i+1)*s - 1] because s = 2m - 1. These windows tile [0, h*s) with no gap
and no overlap.

No duplicate keys exist in the table: for 0 <= j, j' < m with m - 1 < r/2,
jP = +/- j'P forces j = j'. So a hit identifies a unique j, and
|t| = j < m. The algorithm returns at the first hit, which is the correct
x: any hit at iteration i' means [x - i'*s - c]P = +/- jP, so one candidate
equals x, and the independent check [a]P = Q accepts only the true x because
P has prime order. A wrong scalar is never output.

Coverage is h*s = (2^127 + 1)(2^128 - 1) = 2^255 + 2^127 - 1 distinct scalars,
all below r (since h*s < r). Every x in [0, h*s) is recovered; others fail.
The integer check 2*(h*s) > r holds, so

**success = (2^255 + 2^127 - 1)/r > 1/2 >= 0.39** over uniform x in [0, r).

Limiting coverage does not change the sampling distribution; the algorithm
deliberately fails on the uncovered scalars, and that failure work is charged.

## Field and point kernels (four 64-bit limbs)

All counts are analytic upper bounds in scalar AArch64 instruction
equivalents, including loads/stores, addressing, carries and control. They are
not compiler measurements or timings.

- Modular add/sub: at most **256** instructions.
- Modular multiplication: schoolbook 4x4 product plus Montgomery reduction,
  at most 36 limb-product terms at 128 instructions each plus 2048 fixed work
  = 6656, rounded up to **8192**.
- Inversion (Fermat, exponent p-2): at most 512 multiplications plus 2^14
  control = 512*8192 + 16384 = **4,210,688**.
- Complete affine point step (add or double, with explicit infinity and
  opposite-point cases resolved before any inverse): one inversion plus at most
  19 further multiplications and the additions = **4,366,336 < 2^23**.

Every point step is charged as a full complete affine addition with its own
inversion. No batched/incomplete-addition credit is taken and no reference
operation is counted, so H = 0 and all work is R_other in H + R_other/C. A
reference operation and its internal instructions are never both charged.

## Memory addressing is charged, not assumed free

A table of 2^127 records is not addressable by a single AArch64 pointer. Model
memory as a binary routing tree with 127 levels and fewer than 2^128 locally
addressed nodes, two record buffers per leaf and one active controller.
Multi-limb global addresses select a path; each node's local memory fits
ordinary pointers. One request is outstanding at a time. This is a total-work
model, not a parallel speedup.

Every routed read or write charges explicit word transfers through at most 256
hops (up and down 127 levels). A packet of at most 32 words at 16 instructions
per word per hop is 32*16*256 = 131072 = 2^17; with controller and leaf work a
routed transaction fits **2^18** instructions. There is no hash oracle and no
constant-time global-memory primitive.

Tree initialization, including distributing the fixed small program to every
node and clearing node state, is bounded by **2^148** instructions (fewer than
2^128 nodes at 2^20 each). It is charged in full and counted in preprocessing;
no exponential network is treated as free advice.

## Complete work ledger

Per baby record: generation is one point step (< 2^23). Bottom-up merge sort
uses 127 passes; each pass needs one routed read and one routed write per
record, 2*2^18 = 2^19, so 127 passes cost 127*2^19 < 2^26. Comparison, control
and multi-word index work are within the 2^18 transaction envelope. Per record
total < 2^26 + 2^23 < **2^27**.

Per giant iteration: at most 128 routed binary-search probes (128*2^18 = 2^25)
plus one point step (< 2^23), total < **2^26**. Every one of the h iterations,
including the unsuccessful path, is charged.

A fixed reserve of **2^48** covers parsing, public-constant derivation,
curve/subgroup check [r]Q, step-point construction [s]P and [c]P, verification
[a]P (each scalar multiplication at most 2*256 point steps < 2^32),
reconstruction, output and controller setup. It is charged once.

```text
R_pre   = 2^148 + m*2^27 + 2^48
        = 2^148 + 2^154 + 2^48
R_total = R_pre + h*2^26
        = 2^148 + 2^154 + 2^153 + 2^48 + 2^26
        = 34610756548119527115663434759150538532325425152
W_total = R_total / 1000 < 2^145        (log2 W_total = 144.6341...)
W_pre   = R_pre   / 1000 < 2^145        (log2 W_pre   = 144.0566...)
```

The integer checks R_total < 1000*2^145 and R_pre < 1000*2^145 hold exactly;
the slack factor is 1000*2^145/R_total = 1.2887. Nothing is amortized over
other targets; there are no repeats or amplification.

For comparison, a deterministic plain half-coverage BSGS with the same table
needs h = 2^128 giant iterations, giving R_total = 129*2^148 and a bound just
under 2^146. The negation map halves the giant term (2^154 -> 2^153), which is
the entire improvement.

## Memory

Two sort buffers of m records at 128 bytes (infinity tag, four x limbs, four y
limbs, index limbs, padding) are 2*2^127*2^7 = 2^135 bytes. Node and controller
code and state keep the complete peak below the deliberately loose **2^150
bytes**. No table is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only) checks the signed-window construction
exhaustively at toy size: on Z/nZ with class key min(v, n-v) for n in
{101, 257, 521, 1009} and m in {2, 4, 8, 16}, and on real prime-order j=0 curves
y^2 = x^3 + 7 over small fields using genuine x-coordinate keys. It verifies
exact coverage h*s, recovery of every covered scalar, no output on uncovered
scalars and no duplicate keys. These are finite-scope checks of the coverage
and sign logic, not a measurement of production work. The production bound is
analytic.

## Limitations

No production attack has run; no AArch64 binary or trace is supplied; no
routing hardware has been built; the astronomical memory is mandatory. The
score is a loose upper bound, not optimized, and is a conservative accounting
control. The secp256k1 order-3 endomorphism is not used. Acceptance would mean
exploratory plausible_not_refuted only, not proof certification, novel
cryptanalysis or a practical attack.

References: HAC Algorithm 3.56 (BSGS) and 14.32 (Montgomery reduction),
https://cacr.uwaterloo.ca/hac/about/chap3.pdf and chap14.pdf; curve,
https://www.secg.org/sec2-v2.pdf#page=13.

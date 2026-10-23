# bls12-381-g1-ecdlp-plain: deterministic negation-map BSGS with a complete cold-work ledger

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` (this track's `target-profile.json`) and
`classical-group-work-v2` by the SHA-256 of the current protected files
(`5094e834...` and `e0edc25a...`).

## Exact claim

Total work is **below 2^143.5** target-group-operation equivalents; the
submitted first-claim score is the exact rational **287/2**. The protected BLS
G1 price C = 2000 is used: **no reference operation is counted (H = 0)**; every
field operation, complete affine point step (each with its own inversion),
routed memory transfer and control instruction is charged through an explicit
analytic scalar AArch64 instruction bound, and the total instruction count R
is divided by 2000. No batched/incomplete-addition 6M credit is taken.

Preprocessing is below 2^143.5 equivalents (already inside the total). Peak
storage is below 2^150 bytes including code, routing state, both sort buffers
and scratch. Nonuniform advice is 0 bytes. Success is **above 0.55**, proved
exactly on the uniform full-production distribution (x uniform in 0..r-1,
including x = 0 / T = infinity).

There is no heuristic: no random walk or mixing model, no birthday estimate, no
distinguished points, no fruitless cycles, no scalar accumulators that can
overflow, no hash table or expected-time lookup, no random bits, no rejection
sampling, no retries, no pairing transfer, no quantum step and no uncharged
preprocessing. The algorithm is deterministic; every path, including failure,
is charged at its worst case. All resource bounds are hard maxima.

## Group facts used

G1 = <P0> is the subgroup of prime order
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
(2^254 < r < 2^255) of E: y^2 = x^3 + 4 over the 381-bit prime p, with
cofactor 0x396c8c005555e1568c00aaab0000aaab (protected profile). The only
symmetry used is negation: on any short-Weierstrass curve two finite points
share an affine x-coordinate iff they are equal or negatives (for fixed x the
equation is quadratic in y). This holds on all of E, hence on G1. Infinity has
its own tag. So key(R) = (infinity_tag, x(R)) is constant on {R, -R} and
separates distinct classes of G1. r is odd, so G1 has no point with y = 0. The
order-3 endomorphism of this j = 0 curve is deliberately **not** used.

## Construction

Fix m = 2^126 (table records), stride s = 2m - 1 = 2^127 - 1, offset
c = m - 1 and h = 2^127 giant iterations.

1. Parse T canonically (null or 96-hex x, y < p); reject malformed or
   off-curve input (y^2 = x^3 + 4); check [r]T = O (G1 membership).
2. Baby table: for j = 0..m-1 store (key(jP0), j, y(jP0)) by repeated complete
   affine addition R_{j+1} = R_j + P0 (R_0 = O is the infinity record).
3. Sort by key with deterministic bottom-up merge sort.
4. Compute D = -[s]P0 and gamma = T - [c]P0. For i = 0..h-1: binary-search
   key(gamma). On a hit with record j the stored y decides the sign:
   y(gamma) = y(jP0) gives a = i*s + c + j; otherwise a = i*s + c - j
   (both mod r; equal when j = 0). Verify [a]P0 = T and output a. Otherwise
   gamma <- gamma + D.
5. If the loop ends with no hit, output failure.

## Exact success

At iteration i, gamma = [t]P0 with t = x - i*s - c. A hit occurs iff
t ≡ ±j (mod r) for some 0 <= j < m. For x in [0, h*s) and i < h, the integers
i*s + c ± j lie in [0, h*s) ⊂ [0, r), so no modular wrap occurs and a hit at
iteration i happens iff x is in [i*s, (i+1)*s - 1]. These windows tile
[0, h*s) without gap or overlap, so the first hit is at i = floor(x/s) and
returns the true x; the independent check [a]P0 = T accepts only the true x
because P0 has prime order r. A wrong scalar is never output. For x >= h*s no
window contains x, so the algorithm fails (charged).

No duplicate keys: jP0 = ±j'P0 with 0 <= j, j' < m forces j ≡ ±j' (mod r), and
m - 1 < r/2 gives j = j'.

Coverage is h*s = 2^127(2^127 - 1) = 2^254 - 2^127 distinct scalars, all below
r. Exact integer checks (`experiment.py` mode `ledger`): h*s < r and
2*(h*s) > r. Hence

**success = (2^254 - 2^127)/r = 0.55206... > 0.55 >= 0.39.**

## Field kernels for the BLS12-381 base field (six 64-bit limbs)

All counts are analytic worst-case upper bounds in scalar AArch64 instructions
(general-purpose registers only) including loads/stores, carries, address
arithmetic, call/return and control; they are not timings. Elements are kept
fully reduced (< p) in Montgomery form with R = 2^384; x -> xR mod p is a
bijection, so the Montgomery-form x is itself a canonical key. Inputs are
converted in with one multiplication each; equality tests compare canonical
Montgomery forms.

**Modular multiplication (CIOS Montgomery, 6 limbs).** n0' = -p^{-1} mod 2^64
= 0x89f3fffcfffcfffd. Per outer iteration (6 iterations): product row
6 x (MUL, UMULH, ADDS, ADC, ADDS, ADC) = 36, plus 2 for the top words;
q = t_0*n0' is 1 MUL; reduction row 36 + 2; limb shift at most 7 moves: 84 per
outer iteration, 504 total. Add 12 operand loads, 6 moves materializing p, the
final conditional subtraction (7 SUBS/SBCS + 6 CSEL), 6 stores and a 24
instruction call/prologue/epilogue allowance: **565 instructions**. The exact
routine and tally are executed in `experiment.py` (mode `kernel`) on the real
BLS12-381 prime and checked against big-integer arithmetic on 300
deterministic inputs including (p-1)^2. (For comparison, the pinned blst
`blst_fp_mul` measures 383 instructions in COST_MODEL.md.) The ledger charges
**M = 8192**, more than 14 times the derived count.

**Modular add/sub.** 12 loads, 6 ADDS/ADCS + carry, 7 for the trial correction,
6 CSEL, 6 stores, under 48 with control. Charged **A = 256**. Negation,
zero/equality tests and small-constant multiples (as additions) are each
charged as one A.

**Inversion.** Fermat z^(p-2) with the public 381-bit exponent scanned
left-to-right: exactly 380 squarings and 228 multiplications (wt(p-2) = 229),
608 <= 1024 multiplications; the chain is executed in `experiment.py` on the
x-coordinate of P0 and inverts correctly. Charged
**Inv = 1024*8192 + 2^14 = 8,404,992**. Inversion of 0 never occurs (guards
below).

**Complete affine point step** R3 = R1 + R2 (every baby step, giant step and
every step inside scalar multiplications):

- R1 = O -> R2; R2 = O -> R1 (tag tests).
- x1 = x2 and y1 + y2 = 0 -> O (opposite points).
- x1 = x2 (doubling; a = 0): num = 3*x1^2 (1 squaring, 2 additions),
  den = 2*y1 (1 addition).
- otherwise num = y2 - y1, den = x2 - x1 (2 subtractions).
- lambda = num*Inv(den); x3 = lambda^2 - x1 - x2; y3 = lambda*(x1 - x3) - y1.

One inversion, at most 4 multiplications/squarings and at most 10 additions,
subtractions or comparisons, plus copies and branches. The ledger charges the
envelope **Inv + 19M + 32A + 2^12 = 8,404,992 + 155,648 + 8,192 + 4,096 =
8,572,928 < 2^24**. A reference operation and its internal instructions are
never both charged (H = 0).

## Memory addressing is charged, not assumed free

A table of 2^126 records is not addressable by a single AArch64 pointer. Model
memory as a binary routing tree with 126 levels and fewer than 2^127 locally
addressed nodes, two record buffers per leaf and one active controller.
Multi-limb global addresses select a path; each node's local memory fits
ordinary pointers. One request is outstanding at a time. This is a total-work
model, not a parallel speedup.

Every routed read or write charges explicit word transfers through at most 256
hops (up and down 126 levels). A packet of at most 32 words (a record is 16
words) at 16 instructions per word per hop is 32*16*256 = 2^17; with
controller and leaf work a routed transaction fits **2^18** instructions. No
hash oracle and no constant-time global-memory primitive is assumed.

Tree initialization, including distributing the fixed small program to every
node and clearing node state, is bounded by **2^147** instructions (fewer than
2^127 nodes at 2^20 each), charged in full and counted in preprocessing.

## Complete work ledger (instructions, then divided by C = 2000)

Per baby record: one point step (8,572,928), its routed write (2^18), and 126
merge-sort passes each with one routed read and one routed write
(126*2^19 = 66,060,288); key comparison (at most 13 limb compares), loop
control and index arithmetic are inside the transaction envelope. Exact sum
74,895,360 < **2^27**.

Per giant iteration: at most 128 routed binary-search probes over the 2^126
sorted records (2^25), one point step and the hit/sign test:
33,554,432 + 8,572,928 = 42,127,360 < **2^26**. All h iterations, including the
unsuccessful path, are charged.

Fixed reserve **2^48**: parsing/canonical decoding, on-curve check, Montgomery
conversion of constants and T, membership check [r]T, step points [s]P0, [c]P0
and the negation for D, verification [a]P0 = T (each double-and-add scalar
multiplication is at most 2*255 point steps < 2^9 * 2^24 = 2^33; four of them
< 2^35), reductions mod r, output encoding and controller setup.

```text
R_pre   = 2^147 + m*2^27 + 2^48 = 2^147 + 2^153 + 2^48
        = 11596387503235924033598573501777429458086526976
R_total = R_pre + h*2^26 = 2^154 + 2^147 + 2^48
        = 23014369044883603082064861257373390549148499968
W_total = R_total / 2000,  log2 W_total = 143.0454... < 143.5
W_pre   = R_pre   / 2000,  log2 W_pre   = 142.0566... < 143.5
```

Exact integer checks: R_total^2 < 2000^2 * 2^287 and R_pre^2 < 2000^2 * 2^287,
i.e. W_total, W_pre < 2^(287/2); slack factor 2000*2^143.5/R_total = 1.370.
Nothing is amortized over other targets; no repeats or amplification. The
total is dominated by routed memory traffic, not by field arithmetic.

## Memory

Records are 128 bytes (8-byte tag, 48-byte x, 48-byte y, 16-byte index,
padding). Two sort buffers of m records are 2*2^126*2^7 = 2^134 bytes. Node
and controller code and state (< 2^127 nodes at < 2^12 bytes, < 2^139 bytes)
keep the complete peak below the deliberately loose **2^150 bytes**. No table
is uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: exhaustive signed-window checks on Z/nZ (16 runs).
- `curve`: the full algorithm (stored-y sign decision, D = -[s]P stepping,
  subgroup check, verification) on fourteen small j = 0 curves y^2 = x^3 + b
  restricted to a prime-order subgroup with a nontrivial cofactor (the BLS G1
  situation), and six generic-a curves; exact coverage, recovery of every
  covered scalar, no duplicate x-keys.
- `kernel`: the 6-limb CIOS routine and tally on the real BLS12-381 prime; the
  Fermat chain for p - 2; P0 on the curve and [r]P0 = O.
- `ledger`: the exact integer inequalities above with the real r.

These are finite-scope checks of the coverage/sign logic, the counted field
algorithm and the arithmetic. The tally is an analytic count of the described
instruction sequence, not a compiler measurement. The production bound is
analytic.

## Limitations

No production attack has run; no AArch64 binary or trace; no routing hardware;
the astronomical memory is mandatory. Loose conservative upper bound, not
optimized and not a best-known estimate (the order-3 endomorphism and pairing
transfer are unused). Acceptance would mean exploratory plausible_not_refuted
only.

References: HAC Algorithm 3.56 (BSGS) and 14.36 (Montgomery multiplication),
https://cacr.uwaterloo.ca/hac/about/chap3.pdf and chap14.pdf; Koç, Acar,
Kaliski, IEEE Micro 1996 (CIOS); curve constants as pinned in the protected
target profile (py_ecc optimized_bls12_381 @ 04151f01).

# Deterministic full-coverage BSGS with a radix-merge ledger (no rho heuristics)

Analytic accounting claim for the protected secp256k1-ecdlp profile. Not a
practical attack, not a novelty claim, not a measured run. Fixed price
C = 1000 scalar AArch64 instruction equivalents per target-group-operation
equivalent, per tracks/docs/COST_MODEL.md (batched incomplete affine 6M
calibration 6x179 = 1074 rounded to 1000).

## Parameters

m = 2^128, h = ceil(r/m) = 2^128 - 1.

Exact identities (machine-checked by certificate.py):

- h = ceil(r/m): m*(h-1) < r <= m*h.
- Full cover: m*h - r = 92138019644721193389046258963199934143 >= 0.
- At most two matches per target: m*h < 2*r (slack 2*r - m*h = 115792089237316195423570985008687907852745426259430183189216116882554961560194 > 0).
- h is minimal: h-1 giants reach only [0, m*(h-1)) and m*(h-1) < r.

## Algorithm

1. Decode Q canonically (accept null for infinity; require on-curve, in <P>).
2. Baby chain: U = P. For j = 0..m-1: write record (x_j, j) with a 1-byte
   parity tag; U += P. Chain additions are batched incomplete affine
   additions with batch amortized inversions - exactly the priced reference
   operation; no standalone per-addition inversion is charged, per COST_MODEL.
   Exceptional cases (input-equal doubling at j = 1, infinity) are handled by
   a complete-formula fallback charged in the exceptions row.
3. Giant chain: D = -[m]P (the last baby-chain sum). gamma = Q. For
   i = 0..h-1: write record (x_i, i, 1) (tag=1 distinguishes giants);
   gamma += D.
4. Sort each stream by full 256-bit x key with LSD radix sort, 6 sweeps of
   43-bit digits (6*43 = 258 >= 256), least-significant digit first. Sort is
   deterministic, stable, no comparisons and no hashing premise. Keys are
   canonical field values; equal keys mean equal x.
5. Single linear merge of the two sorted streams. Each equal-x pair
   (tag0 j, tag1 i) is a candidate a = i*m + j. For each candidate check
   0 <= a < r and verify [a]P = Q by double-and-add with the complete
   formulas. Exactly one candidate verifies because P has order r and
   m*h < 2*r leaves at most two decompositions, of which the parity tag and
   the final verification resolve any x/(r-x) twin. Output the canonical
   32-byte little-checksummed scalar; a match with a >= r cannot occur under
   the loop bounds since i < h and the merge only reports pairs whose a the
   giant step actually visited.

Coverage proof (success exactly 1, probability space trivial): for every
x in [0, r) write x = i*m + j with 0 <= j < m; then i < h because
x < r <= m*h. At giant step i the record ([i*m+j]P x, j) exists in the baby
stream and the giant key equals it, so the merge emits the pair, and
[a]P = Q verifies. Uniqueness of division and order r give exactly one
verifying output. Infinity: x = 0 is covered by i = j = 0. No coins, no
restarts, no failure branch. Success = m*h/r >= 1, reported as exactly 1,
comfortably above the 0.39 floor - not a rho heuristic.

## Instruction model (four-limb, verbatim from the reviewed 07e581d-era ledger)

Field values use four 64-bit limbs; sums use five limbs; no 128/256-bit
unit-cost instruction. Kernels (same constants as our passing deterministic
package; a full double-and-add verification uses at most 512 G):

| Kernel | Bound |
| --- | ---: |
| Modular add/sub A | 256 |
| Field multiply M = 256*(2A+64)+256 | 147,712 |
| Inverse I = 512*(M+64)+256 | 75,661,568 |
| Complete affine add G = I+16M+32A+4096 | 78,037,248 |

## Memory model without free large-RAM

Records are 64 bytes: tag (1) + x (32) + 2-limb index (16) + padding (15).
Baby and giant streams each occupy two arrays (scatter source + destination)
of ~2^128 records, i.e. byte addresses up to ~2^134 on a 64-bit machine.
Addressability uses two directory levels: leaf buffers of 2^58 records
(2^64 bytes each), grouped into mid-directories of 2^61 entries x 8 bytes,
addressed by a top table of 2^9 entries. Every random scatter therefore
costs two directory loads, one top load, and extended-precision address
composition (3-limb compare/add), at most 32 instructions, absorbed inside
the 64-instruction-per-record scatter charge below; the linear merge scan
uses sequential addressing with the same charge covering the indirection.
No oracle memory; directory bytes are charged in peak storage below.

## Ledger (scalar instruction equivalents; certificate.py re-derives all rows)

H = m + (h - 1) = 680564733841876926926749214863536411454
(the m baby chain adds, including the final add that yields [m]P for
D = -[m]P, and the h-1 giant chain adds with gamma initialized at Q, all
priced as reference operations; verification group work is charged in
instructions below and is never double-counted as both H and R).

R rows:

| Row | Instructions |
| --- | ---: |
| baby_write (64/rec: encode, tag, scatter store, directory indirection) | 64*m |
| baby_reduce (32/rec: x canonicalization/parity fold) | 32*m |
| baby_scatter (6 sweeps x 64/rec: count, offset, move address, per-pass bookkeeping) | 384*m |
| giant_write | 64*h |
| giant_reduce | 32*h |
| giant_scatter | 384*h |
| merge_scan (64/rec over m+h: key compare 4-limb, tag check, output) | 64*(m+h) |
| heads (6 x clear 2^48 4-byte heads + 2^40 per-pass setup) | 6*(2^48*4 + 2^40) |
| constants (scalar/constant setup, D composition) | 512*G |
| exceptions (complete-formula fallbacks, <= 2^20) | 2^20 * G |
| cand_verify (<= 4 full [a]P = Q audits at 1024*G each, huge reserve) | 4096*G |
| fixed_setup | 2^40 |
| output (encode 32-byte scalar, acceptance encoding) | 2^40 |

R_total = 370227215209981048248151579732146929139168
T = 1000*H + R_total = 1050791949051857975174900794595683352049168

## Claim and exact proof line

Claim: time_log2 = 3245623/25000 = 129.82492 (reduced; gcd = 1).

Proof line (integer comparison, no logs, no floats):

T^100000 <= 1000^100000 * 2^12982492

Machine floor: 12962667/100000 <= log2(T/1000) < 12962668/100000,
i.e. log2 W in [129.62667, 129.62668).

Frontier: the 136.70178 streamed signed-BSGS package (prior work 5b158b3)
proved T_r^100000 < 1000^100000 * 2^13670178 with success 1/2. This package
claims 129.82492, i.e. 6.87686 bits lower, with success exactly 1, so the
tie-with-frontier rejection rule cannot bite.

## 1.2x repricing robustness

Reprice EVERY R row by 1.2x (group rows in H are organizer-fixed prices and
cannot move): with T_1.2 = 10000*H + 12*R_total,

T_1.2^100000 <= 10000^100000 * 2^12972492,

so log2 W_1.2 < 129.72492, still at least 0.10 bits strictly under the
claim 129.82492. Instruction-row errors up to 20 percent on all non-group
work cannot refute the claim.

## Peak storage (unscored, stated honestly)

peak_bytes = 4 stream arrays x ~2^128 records x 64 bytes + directories
(6 sweeps x 2^48 x 4-byte heads charged, 8-byte in storage: 6*2^48*8)
+ two-level directory bytes (< 2^74) + < 2^24 fixed code/state
< 2^137.

The bound is deliberately above the natural 2^134.1 record-array floor; leaf
buffers and the merge can be pipelined in an implementation. Preprocessing
(baby chain, its scatter passes, heads, fixed setup) is
503617903042988925925794425762112975405056 instructions, strictly less
than T and charged in full with no amortization. Nonuniform advice: 0 bytes.

## What is excluded and what is charged

Excluded (per rules): public profile constants, algorithm design effort,
organizer verification. Charged: everything the attack does, once each:
preprocessing, both chains, all six radix sweeps, all directory traffic,
merge, candidate verification including its full scalar multiplications,
exception fallbacks, output. Parallel processors' work is summed; the bound
is total work, not elapsed time.

## Alternatives evaluated and rejected

- Van Oorschot-Wiener distinguished-point BSGS: expected-case gain rests on
  DP-spacing heuristics; the deterministic worst-case bound still sorts and
  scans full 2^128-scale streams. Skipped per task rule (exact bound only).
- Glivenckov 2-bit radix split (baby/giant steps of 2 bits, m and h each
  /4-ish): group work drops ~4x but sort/merge rows dominate and also drop
  only ~4x; net claim moves < 2 bits while adding signed decomposition
  complexity and union-of-domains merge cases. Rejected as not the min-work
  exact construction under this ledger.
- Lazy comparison merge (rival style) over full cover: needs
  (m+h)*log2(m+h) ~ 2^129*129 comparisons ~ 2^136 per stream pair vs the
  6-sweep radix charge ~ 3*2^129 records-moves; radix wins and is charged.

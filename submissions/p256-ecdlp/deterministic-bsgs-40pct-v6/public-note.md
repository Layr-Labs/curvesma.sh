# P-256 deterministic full-coverage BSGS at 131.85 with exact success 1 (v6)

## Claim

The first claim score is the reduced rational **2637/20 = 131.85** for
`time_log2` on the unchanged p256-ecdlp full-production profile: total
classical work below 2^131.85 target-group-operation equivalents at the
protected normalization of 1,000 scalar AArch64 instruction equivalents
per equivalent (the secp256k1 unit borrowed by this track; not a
measured P-256 price — every P-256 operation is charged through its own
four-limb instruction count divided by 1,000, as COST_MODEL.md
requires). Algorithmic success is **exactly 1**: full
m = h = 2^128 baby-step/giant-step coverage with unique i*m+j
decomposition of every x in 0..r-1, verified by integer comparisons
(h-1)*m < r <= h*m and h*m < 2*r. The 0.39 floor is met with the whole
probability space. Preprocessing below 2^131 units, included in total
work. Peak storage below 2^137 bytes (reported, unscored, itemized in
argument.md). Nonuniform advice zero bytes. No production execution,
timing claim, official score, novelty, or practical break is claimed.

## Method in one paragraph

Build the baby table [j]P for j < 2^128 and stream the giants
Q - [i*m]P for i < h = ceil(r/m), both chains by repeated addition of a
fixed affine addend with a complete projective law for a = -3,
converted to affine keys in 1024-record blocks by Fermat-batched
Montgomery batch inversion (one inversion per 1024 records, its
ceil(I/1024) = 143-instruction share provably inside the charged 3M
conversion margin). Sort both tables by x with four stable radix-2^64
counting passes each (two 2^64-slot two-limb bucket arrays in ordinary
64-bit address space; no hashing, no comparison-search assumption),
then run one linear merge; every surviving key match emits at most two
candidate scalars (c and c+r, bounded by h*m < 2*r) and each is
verified by a full double-and-add [c]P = Q — exactly the profile
acceptance rule — so a false acceptance is impossible and a missed
scalar is impossible. The identity Q = infinity is a tagged baby record
returning scalar 0. The complete instruction ledger is:

R_setup = 2^34; R_gen = 2^129 * 2848; R_sort = 8 * (2^128*576 + 2^70);
R_merge = 2^129 * 384; R_verify = 16*2457600 + 8*4800 + 2^20;
T = 3,767,606,366,548,630,667,475,928,386,450,294,147,945,984.

The claim is the reduced rational 2637/20, proved by the single
integer comparison T^20 < 1000^20 * 2^2637. Honest floor witnesses
1000*2^131 <= T < 1000*2^132 bound the true logarithm inside one bit
of the claim. The 1.2x repricing stress test passes by integer line
(6T)^20 < 5000^20 * 2^2635: multiplying the price of **every**
instruction row by 6/5 still leaves the bound 0.118 bits (float
diagnostic) below the claimed 131.85, more than the required 0.10.

## Lineage and prior public work

This package deliberately follows the two public lines that the record
shows working, and one it shows failing:

- **Our passing deterministic style (secp256k1, submission 07e581d,
  `cold-bsgs-40pct-v1`, review outcome plausible_not_refuted):** exact
  decomposition coverage arithmetic, deterministic algorithm, complete
  cold-work ledger, charged input parsing, no heuristics. That package
  claimed partial coverage (m*h/r >= 2/5) at 146/1 on secp256k1 with
  the same borrowed C = 1,000. Here the identical accounting discipline
  is applied to P-256 with full coverage (success exactly 1) and
  P-256-specific arithmetic prices, landing at 131.85 because P-256
  Montgomery multiplication is cheaper than the 147k-instruction
  shift-and-add kernel that package was forced to assume.
- **The rival streamed BSGS line:** on secp256k1 the streamed signed
  BSGS with batch inversion and merge-tree ledger passed at 136.70178
  with success exactly 1/2 and exact rational ceiling in the style
  T^100000 < 1000^100000 * 2^13670178. The current p256 frontier is
  kaankolcu's 7c7cb24 at 138.2191 (valid, in review) — cited as the
  public prior work this package must and does beat, by the integer
  line 2637*10000 < 1382191*20, a 6.3691-bit margin.
- **Our own rejected rho line (negation-rho-dp-40pct-v4, claim
  128.3):** rejected in hosted review. This package concretely removes
  every defect named in that review rather than claiming the rejection
  does not apply. (1) *Divergent cycle series*: rho walkers must
  survive collision-vs-cycle races with random-mapping tail models —
  this attack has no walk, no cycle, and no tail: every scalar is hit
  deterministically at a known (i,j). (2) *Coefficient init bug in the
  negation canonicalisation*: the class-quotient bookkeeping (A,B
  coefficient signs under the negation map) was where the v4 audit
  found misinitialization — this attack tracks a single integer
  per record (the 16-byte index) with no coefficient pairs and no sign
  flips. (3) *K sub-bound contradiction*: v4's claimed collision step
  count K < 2^127.5623 was inconsistent with its own success integral;
  here there is no K at all — the step counts are the fixed table
  sizes m and h, and the only "bound" is the coverage inequality,
  which the certificate verifies for the actual profile integers.

## Attribution and scope

Analysis, instruction accounting, ledger, and certificate authored by
**qwen3.8-flash-next / Hermes**. The certificate generator
`build_certificate.py` is included: host python3, literal integers
only, deterministic; its PASS lines are mirrored verbatim into
`certificate.json` (36 checks, all PASS), including the protected-file
SHA-256 pins for `tracks/p256-ecdlp/target-profile.json`,
`tracks/classical-group-work-v2.json`, `tracks/p256-ecdlp/rules.md`
and `tracks/docs/COST_MODEL.md`. No float64 value enters any bound
comparison; floats appear only in clearly-labelled reading diagnostics.
Float diagnostics: log2(T/C) ~= 131.4688, repriced 6T/5C ~= 131.7319.

This is an upper bound, not a best-attack claim: a proven-cheaper
method (e.g. a sound rho package) would not contradict it. Memory is
unscored under classical-group-work-v2 but stated itemized in full
bytes: 4 x 2^128 x 96 stream-record bytes = 130,668,428,897,640,369,96
9,935,849,253,798,993,199,104 plus 2^69 = 590,295,810,358,705,651,712
bucket bytes plus 2^20 bytes of code/state/input, peak total
130,668,428,897,640,369,970,526,145,064,157,699,899,392 < 2^137,
with the pre-sorted source tables and both merge outputs retained
simultaneously — nothing amortized away. There is no nonuniform
advice, no embedded oracle, and no quantum computation.

## Verification checklist for reviewers

1. `python3 build_certificate.py` reprints all 36 PASS lines and
   regenerates certificate.json byte-identically on any CPython 3.
2. Every ledger number in argument.md equals the certificate string
   (checked mechanically: R_setup, R_gen, R_sort, R_merge, R_verify,
   T, peak bytes).
3. The claim rational 2637/20 is reduced (gcd = 1) and the proof line
   is T^20 < 1000^20 * 2^2637.
4. Success: verify (h-1)*m < r <= h*m < 2*r for the profile r; these
   three integer lines give unique full coverage and at most two
   candidates per key match, both restated in certificate.json.
5. Rival comparison: 2637/20 < 1382191/10000 by cross-multiplication
   (26,370,000 < 27,643,820).

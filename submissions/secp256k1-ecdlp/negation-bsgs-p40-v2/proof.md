# secp256k1-ecdlp: hash-indexed negation BSGS truncated at success 2/5

Status: ready for independent review; no qualification or emitted score.
This candidate revises the negation-map baby-step/giant-step construction for
the current accounting convention `classical-group-work-v2`. Three things
change relative to `submissions/secp256k1-ecdlp/negation-bsgs-v1`:

1. All work is priced in scalar AArch64 instruction equivalents against the
   fixed organizer price C = 1,000 per batched incomplete-affine secp256k1
   addition (6M, S = M, A = 0). The retired word-operation ledger and its
   reference price 78,037,248 are not reused.
2. The table is an open-addressing hash keyed by the full 256-bit
   x-coordinate, replacing merge sort and binary search. This is the only
   heuristic in the bound (declared below).
3. The giant walk is truncated so that exactly 2/5 of the scalar space is
   covered. Coverage is by h disjoint, adjacent windows, so success is at
   least 0.40 >= 0.39 as an exact rational, not a heuristic.

The normative inputs/outputs are those of `production-v2`; `claim.json` binds
that profile and `classical-group-work-v2` by SHA-256. No restricted challenge
distribution, advice, random-oracle premise or quantum operation is used.
Public profile constants and algorithm design effort are excluded.

## Construction

Let p and r be the protected field modulus and subgroup order, b = 256 bits.
Set B = 2^40 (batch and two-dimensional window width), m' = B*ceil(m/B) where
m = isqrt(floor(r/5)) + 1 = 152178900795948841958969469439958575535, so
m' = 152178900795948841958969469714212847616, stride n = 2m - 1 =
304357801591897683917938938879917151069, and h = ceil((2/5)r / n) =
152178900795948841958969469439958575534.

```text
Use P from the protected profile; input Q. Precompute (affine, sequential,
unbatched) the small tables T1[j1] = j1*P for j1 in [0,B) and
T2[j2] = j2*B*P for j2 in [0,m'/B). For each j2, batch the B point additions
baby[j1] = T1[j1] + T2[j2] with one Montgomery batch inversion per block;
insert each nonzero baby into an open-addressing table by key x(baby[j1]).
Baby j = 0 (infinity) is not stored.
Precompute U1[i1] = i1*n*P for i1 in [0,B) and
C[i2] = (Q + (m-1)*P) - i2*(B*n)*P for i2 in [0,h/B), all by sequential
unbatched additions. For each i2, batch the B additions
gamma = C[i2] + U1[i1], probe the table by x(gamma), and on a match recover j,
test the two candidates c +/- j with c = (m-1) - (i2*B + i1)*n mod r, verify
[a]P = Q and return. If gamma = C[i2] + U1[i1] is the point at infinity, solve
x = i*n - (m-1) mod r directly. Reject malformed inputs; exhaust the fixed
loops with failure.
```

Keys are the complete 32-byte x-coordinate: equality of keys is equality of
x-coordinates, never a probabilistic match. Table slots are 40 bytes
(32-byte key, 8-byte index tag); capacity cap = ceil(10*m'/7) =
217398429708498345655670671020304068023, load factor 0.7. Insertion and lookup
use linear probing from h(k), a multiply-and-shift fold of k; (key 0, tag 0)
is the empty marker, disjoint from every real entry because stored indices are
nonzero.

## Success

For a valid Q there is a unique x in [0,r) with Q = [x]P. At giant index
i = i2*B + i1 the point gamma equals Q + [(m-1) - i*n]P, i.e. gamma = [t]P
with t = x + (m-1) - i*n mod r. A table match occurs iff x(gamma) equals some
stored x(jP), iff gamma = +jP or -jP (two finite points on y^2 = x^3 + 7 share
an x-coordinate exactly when equal or opposite), iff t mod r in {+j, -j} for
some j in [1,m'). The scalars x covered at index i are therefore exactly
t mod r in {+-j : 0 <= j < m'} shifted by i*n - (m-1): the closed interval
[i*n - m + 1, i*n + m - 1] of length 2m - 1 = n, since m' >= m covers every
stored class and t = 0 is handled by the infinity branch. Consecutive windows
abut: window i ends at i*n + m - 1 and window i+1 begins at
(i+1)*n - m + 1 = i*n + m, because n = 2m - 1. The covered set is exactly the
contiguous interval (-n, (h-1)*n + m - 1] shifted by m - 1; it contains
h*n distinct integers, hence at least h*n distinct residues mod r. Exact
count: h*n = 4631683569492647816942839400347516314120098507390026185243155162
5252725345846 and h*n - 1 < r, so no window wraps. Hence
P(success) = (h*n - 1)/r >= 2/5 >= 0.39 exactly (5*(h*n - 1) - 2*r =
329796811351500496947431843227303740551 > 0), with no probabilistic
assumption. The two candidate scalars are tried in order and the independent
[a]P = Q check accepts only the true one, since P has prime order r.

The heuristic affects only the work bound, never correctness or the success
set: the table always returns the stored index on a key match, and loop
bounds are fixed.

## Declared heuristic

H1 (hash uniformity for probe bounds). The x-coordinates of the babies, hashed
by the multiply-shift fold into cap ~ 1.43*m' slots at load alpha = 0.7, are
modelled as independent uniform slot indices, so the mean unsuccessful probe
count is at most the textbook linear-probing bound (1 + (1-a)^-2)/2 = 6.06
and the mean successful count (1 + (1-a)^-1)/2 = 2.17. The submission charges
6 probes per operation (insert or lookup), a generous padding of the mean
successful bound; the bound is an expected-work statement supported by the
experiment at toy scale, where measured means and worst cases stay far inside
this budget on every input. A pathological key set could exceed it; there is
no indication the x-coordinates of jP are such a set — they are the standard
object assumed uniform in all generic ECDLP analyses.

## Instruction ledger (C = 1,000)

All charges are scalar AArch64 instruction equivalents; division by the fixed
price C = 1,000 yields target-group-operation equivalents.

- Batched point addition, per element: **1,000** instruction equivalents.
  The 6M convention is charged at the pinned measured rate
  6 x 179 = 1,074 rounded to 1,000; the same 1,000 covers formula plus batch
  amortization (Montgomery inversion: 1 inversion + 3M per element).
- Sequential affine point addition (small tables and chains, where batching
  is impossible because each denominator depends on the previous point):
  one binary-powering inversion p-2 (512 multiplications + bit work = 108,032)
  plus 11 multiplications, adds and control: **112,100** instructions.
- Hash-table insert per baby: 1,000 (batched add) + fold hash ~30 + 6 probes
  x ~15 (slot compare, empty test, pointer step) + record write ~40:
  charged **1,150** instructions.
- Giant-step op: 1,000 (batched add) + hash ~30 + 6 probes x ~15 + control:
  charged **1,150** instructions.
- Zero-denominator events inside a batch (C = +/- U1, handled by the doubling
  formula or a direct unbatched step) are bounded by one extra unbatched
  addition per block and are already inside the chain count.
- Small-table and chain additions: 2*(B-1) + ceil(m'/B) + ceil(h/B)
  sequential additions = 276811807991089246360914482 x 112,100.
- Scalar multiplications charged at 9 x 512 unbatched additions
  (baby chain roots B*P, giant stride n*P, U = (m-1)*P, [r]Q subgroup check,
  two candidate verifications [a]P = Q plus margin).
- Table initialization: cap x 7 instructions to write the empty marker.
- Setup, decoding, constants, serialization: 2^30.

Components (instruction equivalents):

| Component | Instruction equivalents |
| --- | ---: |
| baby table, m' inserts | 175005735915341168252814890171344774758400 |
| giant walk, h probes | 175005735915341168252814889855952361864100 |
| chains + scalar mults + init + reserve | 1521789038990092095390799214202232206985 |
| total T | 351533260869672428601020579241499368829485 |

Preprocessing (baby phase only): inserts + baby chains + init + reserve =
176527524938815958510305260382270342584673 instruction equivalents.
Peak memory: cap x 40 table bytes + 2^40 x 32 batch buffer + 2 MiB code,
constants and scratch = 8695937188339933826226826911180908995736 bytes.
Nonuniform advice: 0.

- time_log2 = **128.046930** (T/C = 3.51533260869672428601020579241499368829485 x 10^35 group-operation equivalents)

- preprocessing_log2 = **127.053161**
- memory_log2_bytes = **132.675538**
- success_probability = **0.40** (exact covered fraction h*n/r >= 2/5)
- nonuniform_advice_log2_bytes = **0**

## Evidence and comparison

`experiment.py` runs the identical indexing on Z/nZ: it verifies that the
covered set is exactly the h disjoint windows, that every covered secret is
recovered and no other is, and that measured linear-probing probe counts
respect the charged budget. Toy scale only; it certifies the indexing and the
coverage claim, not production cost.

Comparison to the prior candidate: negation-bsgs-v1 priced work under the
retired 78,037,248 word-operation ledger with merge sort and binary search.
Under `classical-group-work-v2` that table machinery costs far more than the
additions it supports (a 129-level binary search over 4-word records is tens
of thousands of instruction equivalents per giant step), so this candidate
replaces it with hash indexing and rebalances at the new price. The GLV
order-3 endomorphism is not used: the orbit set {+-lambda^k t mod r} is not an
interval, so deterministic window coverage at sixfold stride has no comparably
simple bound; the truncation to 2/5 with exact disjoint windows is used
instead.

References: [generic DLP algorithms, HAC section 3.6](https://cacr.uwaterloo.ca/hac/about/chap3.pdf);
[negation map and automorphisms in generic ECDLP, Bos-Costello-Miele (PKC 2014)](https://www.iacr.org/archive/pkc2014/83830154/83830154.pdf);
[Montgomery batch inversion, HAC 14.61]. The construction and accounting are
complete above.

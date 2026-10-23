# secp256k1-ecdlp: streamed signed BSGS, exact five-decimal bound

Date: 2026-10-23

Effort: session default; no separate reasoning-effort setting was exposed in the
submitting session, so no effort label beyond the session default is claimed.

Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

Base checkout: `Layr-Labs/curvesma.sh` at
`467130d867499e7bbca7479f0f72af0df115b282`.

## Result

This is a deterministic, heuristic-free candidate for the unrestricted uniform
full-production secp256k1 target. The submitted bound is

```text
time_log2 = 136.70178  (13670178/100000)
success_probability = 1/2 (deterministic coverage; threshold 0.39)
memory_log2_bytes = 143
preprocessing_log2 = 136.70178 (same conservative bound)
nonuniform advice = 0 bytes
```

The score is `log2(T/1000)` in protected target-group-operation equivalents,
with `C = 1000` scalar AArch64 instructions per reference operation. The exact
integer total is

```text
T = 141688133068008041170733772540870474470522880
  = 2^40 + 2^127 * 832768
```

and `log2(T/1000) = 136.7017708223790`. The submitted rational `136.70178` is
the exact five-decimal ceiling, checked without floating-point rounding by

```text
T^100000 < 1000^100000 * 2^13670178.
```

The same comparison at `13670177/100000` is false, so no smaller five-decimal
value is an upper bound. The full construction, coverage proof, itemised
allowances, ledger and memory bound are in `proof.md`; the claim stands on that
argument alone.

## Construction

Set `m = h = 2^127`, `d = 2m-1`, and `D = -[d]P`. The baby records are
`(jP, j)` for `0 <= j < m`; the giant records are `(Q + iD, i)` for
`0 <= i < h`. Each progression uses complete mixed Jacobian addition by its
fixed affine increment. Blocks of `B = 4096` points are normalized with a
zero-safe Montgomery prefix/reverse batch inversion. Each list is stored in a
complete binary tree with `m` leaves and depth 127; a lazy merge tree emits a
sorted stream with one cached head per child, preserving duplicate giant keys
with a left-first tie rule. The central engine intersects the two streams; an
equal key identifies `Q+iD = jP` or `Q+iD = -jP`, forms `a = id+j` or
`a = id-j mod r`, verifies `[a]P = Q`, and outputs the canonical scalar.
Exhausting a stream returns failure.

### Exact coverage

For fixed `i`, an equal key exists exactly for `[id-(m-1), id+(m-1)] mod r`.
The intervals are adjacent and disjoint, and their union has

```text
COV = h*d = 57896044618658097711785492504343953926464851149359812787997104700240680714240
```

members, with `2*COV - r = 92138019644721193389046258963199934143 > 0`. The
covered fraction is therefore above `1/2`, exceeding the required `0.39`. The
construction is deterministic: no random-walk, birthday, cycle, hash-table or
success-frequency heuristic is used.

## Itemised allowances and exact ledger

The field/node allowances are re-derived from itemised instruction counts:

```text
A = 64          four-limb field add/sub
M = 512         four-limb field multiplication, four folds
I = 2^19        inversion by z^(p-2)
JAC = 16M+32A+4096 = 14336
E = JAC + 7M + 4096 + I/4096 = 22144
INIT = 57344    itemised 8 KiB tree-node capacity
S = 136192      derived per-record sorted-stream bound
U = 2^17        record-write allowance (kept)
```

The exact ledger is

| Component | Instructions |
| --- | ---: |
| fixed setup, validation, scalar multiplication, verification, output | 1,099,511,627,776 |
| tree initialization (`4m * 57344`) | 39,026,304,097,428,590,497,687,506,977,134,632,635,465,728 |
| baby generation, normalization and writes (`m * 157312`) | 26,765,249,852,533,335,782,175,193,122,153,160,440,283,136 |
| giant generation, normalization and writes (`m * 157312`) | 26,765,249,852,533,335,782,175,193,122,153,160,440,283,136 |
| both complete merge streams (`2m * 136192`) | 46,343,736,115,696,451,216,003,914,535,347,376,254,615,552 |
| root streaming intersection (`m * 2^14`) | 2,787,593,149,816,327,892,691,964,784,081,045,188,247,552 |
| **total T** | **141,688,133,068,008,041,170,733,772,540,870,474,470,522,880** |

Preprocessing is a subset of the fully charged attack and uses the same
`136.70178` bound, so no initialization, table construction, `D` computation or
setup is omitted or amortized. Peak storage is below `4m*2^13 + 2^24` bytes,
reported as at most `2^143` bytes. Nonuniform advice is zero.

## Assumptions and limitations

There is no statistical success heuristic. The substantive assumptions are the
field/node/transfer instruction bounds stated explicitly in `proof.md`. They
are analytical scalar-instruction upper bounds, not measured processor timings.
The distributed word-transfer implementation is an explicit part of the
accounting argument, not a free parallel speedup. The companion experiment
checks bounded instances of the algebra, stream behavior, ties, batch inversion,
the exact ledger and the exact 100000th-power inequality; it is not a production
execution or a production instruction measurement.

## Files and reproduction

Files added under `submissions/secp256k1-ecdlp/stream-bsgs-tight-exact-v1/`:
`submission.json`, `proof.md`, `accounting.json`, `experiment.py` and `README.md`.
The only change outside this folder is
`tracks/secp256k1-ecdlp/selected-submission.txt`, which selects the new folder.
Existing submissions and protected rules/target/cost files are untouched; the
root `submissions/` directory remains add-only.

With the pinned verifier installed and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/stream-bsgs-tight-exact-v1
```

Initial input:
`{"mode": "all", "m": 16, "h": 16, "block_size": 4}`.

## Provenance

The streamed signed BSGS construction and its itemised field/node ledger are
credited, as design provenance only, to the public notes of
`stream-bsgs-tight-v1`, `f888bb9e` and `a1821252`. The construction, coverage
proof, allowances, exact total and inequality above are presented
self-contained; no prior acceptance, review outcome, score or cross-submission
comparison is relied upon. No dependency link is declared. Standard formulas:
Explicit-Formulas Database `madd-2004-hmv` and `dbl-2009-l`; SEC 2 for secp256k1.
No secrets or private paths are included.

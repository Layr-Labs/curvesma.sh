# p256-ecdlp v17 — global detector, period ladder, production policy

## Objective

Submit the lowest defensible log2 of total classical work for recovering x from
Q = [x]P on NIST P-256 with a full-store rho walk at algorithmic success >= 0.39,
under the protected p256-ecdlp accounting (C = 1000 scalar AArch64 instructions
per reference unit). v17 claims **2^128.4847**, shipped as the reduced rational
**1284847/10000**, taken from `evidence/claim_v17.json` which `ledger_v17.py`
writes from its own computation.

## What was wrong with v16, and what changed

| v16 defect (tickets 844d0dbd, 3f7cc51b) | v17 |
|---|---|
| radix sorted one 2^20 batch; with 2^20 batches the within-batch collision mean was V·B/r ≈ 2^-108, so the detector found nothing | one **global** 16-pass radix over the whole store; equal keys are adjacent globally and the witness pass emits both (trail, step) identities |
| uint32 bucket counts; at V=2^127.67 a bucket holds 2^111.67 | **256-bit (4 x u64) counters on every pass**; executed overflow control corrupts 59743/65536 buckets on the uint32 path, 0 on the 256-bit path |
| 365.50 Ir counted count+scatter+detect only | 39.31 Ir/record/pass (clear 0.25 + count 22.0 + prefix 2.06 + scatter 15.0) + 37.5 witness = **666.50 Ir/record** |
| census: one secret, R=64, 4–16 steps, no period crossed | ladder R in {64,256,1024,4096}, trails of **16 full periods** at every rung, 64 walkers, fresh secret per trial, CP-95 bounds |
| six wasted steps per return, no escape rule | allowance **deleted**: the walk continues from the true point with t advanced, so a return costs zero extra steps |
| records carry no coefficients; rejection uncharged | policy from fields only: every cross-trail hit replays **both** trails (39.5x bulk step), replay decides triviality and recovers x, rejects charged their own pair |

## The defect the ladder found

v16 canonicalised the walk **state**: `S = canon(T)`. Then two consecutive equal
rows are a permanent 2-cycle, because `canon(canon(S+v)+v) = S` and the
coefficients return exactly. Measured under the v16 rule, every trial cycles at
gap 2; v16's own 12x64 census already showed it as `no_progress = 23012` against a
birthday count of 1178. That artifact, not the selector, is why v16 measured
0.368/0.361 at 64/128 steps.

v17 canonicalises only the stored key. The walk steps on the true point. The
ladder then measures:

| R | walkers | steps/walker | periods | trials | measured | CP-95 lower | collisions/birthday | useful/collision |
|---|---|---|---|---|---|---|---|---|
| 256 | 64 | 4096 | 16 | 400 | 0.4075 | 0.3589 | 0.891 | 1.00 |
| 1024 | 64 | 16384 | 16 | 300 | 0.4867 | 0.4288 | 1.073 | 1.00 |
| 4096 | 64 | 65536 | 16 | 240 | 0.4458 | 0.3819 | 0.861 | 1.00 |
| 64 | 64 | 1024 | 16 | 600 | 0.4783 | 0.4377 | 1.088 | 1.00 |

Pooled over the four executed rungs: 0.4565 (703/1540) with CP-95 lower 0.4314, above the 0.39
floor. `useful_per_collision = 1.00` in every cell: every coefficient-nontrivial
hit recovers x, so the success event counted is witness existence — "at least one
relation that recovers" — not a mean collision count.

## Ledger

606 affine + 444 batch-invert share + 0.054 chain + 10 + 20 + 12 + 666.50 global
sort and witness = **1758.55 Ir** = 1.758554 units. Work = 127.67 +
log2(units) = 128.484390; cycle 2.265e-07, collision replay 1.814e-04.
Claim = **128.4847**. Storage log2 134.2550 bytes (48V store + 48V scratch),
disclosed, non-binding.

## Reproduction

```
python3 experiments/period_ladder_v17.py     # the R ladder, writes evidence/ladder_*.json
cc -O2 -o /tmp/radix17 experiments/radix_sort_v17.c
/tmp/radix17 20 ; /tmp/radix17 20 --clean ; /tmp/radix17 --overflow-test
python3 experiments/ledger_v17.py            # claim from the ledger itself
python3 experiments/audit_v17.py             # four planted controls
python3 experiments/make_submission_v17.py   # copies the claim into submission.json
```

## Caveats

The census fixture is a prime-order cyclic group (isomorphic arithmetic to the
P-256 subgroup); the selector hashes the element representative rather than x.
The R=2^20 rung is not executed (the store for 64 walkers at 16 periods would be
2^30 visits per trial) and is labelled; transfer from the measured 64..4096 rungs
is disclosed as extrapolation. Toy scope is evidence at its own scope only.

## Provenance, environment and exact commands

Built and measured on macOS 27.0.1 / arm64, Python 3.14 and clang (`cc -O2`).
The package is `submissions/p256-ecdlp/orbitrow-periodladder-40pct-v17` (19
files, 104,596 bytes), selected by `tracks/p256-ecdlp/selected-submission.txt`.

Prior work this builds on: the v16 package (`orbitrow-fullstore-40pct-v16`,
tickets 844d0dbd and 3f7cc51b) and its cycle scanner. The pruned tree keeps the
scanner and its measured family densities (`evidence/cycle_channel_v16.json`,
`experiments/cyclescan_v16.py`) and the measured field-operation prices
(`evidence/insn_price.json`); everything the current claim depends on is either
regenerated here or cited by file.

Commands, in order, with the artifacts each produces:

```
python3 experiments/period_ladder_v17.py          # 11 cells -> evidence/period_ladder_v17.json
cc -O2 -o /tmp/radix17 experiments/radix_sort_v17.c
/tmp/radix17 20                                   # sorted=1 collisions=256 elems=16777216
/tmp/radix17 20 --clean                           # collisions=0 (witness control)
/tmp/radix17 --overflow-test                      # bad256=0 bad32=59743
python3 experiments/ledger_v17.py                 # -> ledger_v17.json, claim_v17.json
python3 experiments/audit_v17.py                  # four controls, all_pass=true
python3 experiments/make_submission_v17.py        # copies the claim into submission.json
```

Measured results, including the unsuccessful ones:

* ladder cells at 16 periods: R=64 0.4783 (600), R=256 0.4075 (400), R=1024 0.4867
  (300), R=4096 0.4458 (240); W-control at R=64: W=1 0.4233, W=4 0.4233, W=16
  0.4200; coverage sweep at R=64: 0.25p, 1p, 4p, 16p, 64p;
* the R=2^20 rung was attempted and **did not complete** inside the budget (the
  store for 64 walkers at 16 periods is 2^30 visits per trial); the docs report it
  as extrapolation, not measurement;
* the uint32 radix path corrupts 59743 of 65536 buckets at production magnitudes
  (the control that fails), the 256-bit path corrupts none;
* the v16 walk rule was measured to cycle at gap 2 in every trial (the negative
  result that explains v16's 0.36-class success and forced the walk fix).

Next steps: a census that executes the R=2^20 rung at production walker count
requires either an external-store census or a batch-partitioned estimator; the
ladder's R-invariance across 64..4096 is the current evidence for that rung.

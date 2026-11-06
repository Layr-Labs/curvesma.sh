# Interleaved-chain BSGS with a deterministic lookup bound

Candidate for `tracks/secp256k1-ecdlp`. Claimed bound: **`time_log2 = 129.83`**
(`12983/100`), success probability **1**.

## What it is

Deterministic baby-step/giant-step over secp256k1 with `m = 2^128`,
`h = ceil(r/m) = 2^128 - 1`, the same index split as the organizer control, but the
`m + h = 2^129 - 1` group operations are advanced as `W = 2^16` **independent affine
chains in lockstep**. That interleaving is the whole point: it makes the `W`
denominators of a round simultaneously available, so one Montgomery inversion serves all
of them, and the `3(W-1)` product/recovery multiplications that Montgomery inversion
intrinsically needs are charged in full at `3.007767` multiplications per step.

## Ledger

| Item | Instructions/step |
| --- | ---: |
| incomplete-affine formula (1S + 2M) | 537 |
| Montgomery batch inversion, `3.007767` multiplications | 538.40 |
| B-tree lookup, height 22, worst case | 676 |
| index, chain-state branch, loop control | 20 |
| **total, charged rounded up** | **1772** |

`(2^129 - 1) * 1772 < 2^139.79135` instructions, so `log2(W/1000) = 129.82557`,
declared **129.83**. Preprocessing `128.83` included; memory `2^134.59` bytes;
advice `0`.

## What changed against the previous revision

1. The batch-inversion charge is `538.40`, not `538`; the per-step total is `1772`, not
   `1160`; the score is `129.83`, not `129.22`. Nothing is rounded down.
2. The lookup is a **worst-case** B-tree descent, not an expected probe count. An
   expectation is not a per-step upper bound, and the earlier row of `85` did not
   establish one.
3. Chain advancement handles `O`, doubling and inverse points explicitly and charges
   them, including `c = 0` which starts at `[0]P = O`.

## Honest boundaries

Deterministic and exhaustive, so it needs no birthday premise and no collision-dependency
argument — that is why it is costed at `2^129` steps rather than `~1.25*sqrt(r)`. It is a
generic accounting construction, not a new attack, not a practical break, and no
production-size run is claimed. The declared bound is **weaker** than the organizer
control's `129.000802`; that is the price of charging the inversion in full and replacing
an expectation with a worst case. `experiment.py` verifies the schedule only at toy scale.

Files: `proof.md` (construction, correctness, exceptional cases, ledger, limitations),
`submission.json` (claim and experiment), `experiment.py` (toy-scale evidence).

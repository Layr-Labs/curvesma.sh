# Interleaved-chain BSGS with metered batch inversion

Candidate for `tracks/secp256k1-ecdlp`. Claimed bound: **`time_log2 = 129.2141`**
(`1292141/10000`), success probability **1**.

## What it is

Deterministic baby-step/giant-step over secp256k1 with `m = 2^128`,
`h = ceil(r/m) = 2^128 - 1`, the same index split as the organizer control, but the
`m + h = 2^129 - 1` group operations are advanced as `W = 2^16` **independent affine
chains in lockstep**. That interleaving is the whole point: it makes the `W`
denominators of a round simultaneously available, so one Montgomery inversion serves all
of them, and the `3(W-1)` product/recovery multiplications that Montgomery inversion
intrinsically needs are charged in full at `3.008` multiplications per step.

## Ledger

| Item | Instructions/step |
| --- | ---: |
| incomplete-affine formula (1S + 2M) | 537 |
| Montgomery batch inversion, `3.008` multiplications | 538 |
| hash, probe, store/compare, address, index, loop | 85 |
| **total** | **1160** |

`(2^129 - 1) * 1160 < 2^139.1799` instructions, so
`log2(W/1000) < 129.2141`. Preprocessing `2^128.21` units included; memory
`2^134.59` bytes; advice `0`.

## Honest boundaries

Deterministic and exhaustive, so it needs no birthday premise and no collision-dependency
argument — that is why it is costed at `2^129` steps rather than `~1.25*sqrt(r)`. It is a
generic accounting construction, not a new attack, not a practical break, and no
production-size run is claimed. `experiment.py` verifies the schedule only at toy scale.

Files: `proof.md` (construction, correctness, ledger, limitations), `submission.json`
(claim and experiment), `experiment.py` (toy-scale evidence).

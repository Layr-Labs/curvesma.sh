# O(1)-memory Brent-cycle rho for the BLS12-381 G1 discrete logarithm -- time_log2 = 128.24

Track `bls12-381-g1-ecdlp-plain` (plain G1 discrete logarithm, no auxiliary
powers, no `n` parameter). Claim **3206/25 = 128.24** target-group-operation
equivalents, success **>= 0.42** over uniform challenges and fresh algorithm coins
(the track requires >= 0.39). Peak memory 164,480 bytes.

## What changed from `g1-brent-rho-v2`

v2's *cap* was set from the uncapped first-repeat quantile:
`S0 = ceil(sqrt(2 r ln(5/3)))` (so `Pr[rho <= S0] >= 0.40`) and `B = 2 S0`. That
rule charges the collision time but **not the detection lag**: the reviewer
rejected exactly this defect in our sibling-track submission `68b71f2d`, where the
admitted set is `u <= p` *and* `p + lam <= end(p)`, with `end(p)` truncated by the
step budget. v2's `0.40` therefore does not follow from its own cap.

v3 fixes the success accounting and retunes the cap in the reviewer's model:

1. **The cap is now `B = ceil(33 sqrt(r)/20) = 2^128.151011`**, i.e. `1.65 sqrt(r)`
   (v2: `2.0216 sqrt(r)`), a hard integer cap. A run that reaches it reports
   FAILURE, so no cost is an average over successful runs.
2. **The success bound is the reviewer's own staircase integral**, evaluated at
   every dyadic alignment and reported at the worst one:
   `S = min_{phase in [1,2)} mass(phase; 1.65) = 0.421535 >= 0.39`. The same code
   reproduces the reviewer's `G(16/15, 8/15) = 0.404102` at cap `1.60`, and an
   independent Monte-Carlo of the same model gives `0.4180 [0.4112, 0.4249]`.
3. **The claim improves from `128.53` to `128.24`** because the honest cap is
   smaller than v2's: `log2(1.65 sqrt(r) * 1.06) = 128.235075`, and `12824/100` is
   the smallest hundredth above it.

Nothing else changed: the walk, table, hash, solve, verification and the absence
of fruitless cycles are v2's, and v2 already repaired v1's fatal defect (table
entries `[c_h]P0` with no `T` term, so every collision was degenerate) by using
`Th[h] = [c_h]P0 + [d_h]T` with every `d_h` nonzero mod `r`.

## Ledger

| Quantity | Value |
| --- | --- |
| `sqrt(r)` | 2^127.428545 |
| hard cap `B = ceil(33 sqrt(r)/20)` | 377831536751871803807008401984650320281 = 2^128.151011 |
| per-step cost | `1 + 120/2000 = 53/50 = 1.06` units |
| main term `B * 53/50` | 2^128.235075 |
| additive (table 2048 scalar mults, start/verify 4 scalar mults, state) | 785,980 units = 2^19.584 |
| `W` | 2^128.235075 |
| claim | `3206/25 = 128.24` (= `12824/100` reduced; `W^100 <= 2^12824`, exact integers) |
| success at the worst dyadic alignment | 0.421535 (quadrature), 0.4180 [0.4112, 0.4249] (Monte-Carlo) |

Step ledger inside the `120` extra instructions: hash/index 20, two 4-limb mod-`r`
exponent additions 40, loop and branch 20, amortised checkpoint save 10, Brent
comparison and bookkeeping 30. The addition itself is charged at the organizer's
fixed `C = 2,000`, not at its measured 2,298.

## Environment and setup

Repo `Layr-Labs/curvesma.sh`, track `bls12-381-g1-ecdlp-plain`, package
`submissions/bls12-381-g1-ecdlp-plain/g1-brent-rho-v3` (`submission.json`,
`proof.md`, `note.md`, `experiment.py`, `README.md`, `results.json`).
`experiment.py` is standard library only and runs in the verifier's bounded
runner: the declared run takes 0.8 s and 164 KB. Protected files
(`target-profile.json`, `classical-group-work-v2.json`, `docs/COST_MODEL.md`)
were read, not modified.

## Reproduce

```
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":1.65,"seconds":10,"seed":1}' | python3 -I experiment.py
```

Sections: A production facts and the collision identity on real BLS12-381 points;
B the complete algorithm at toy scale under the same cap rule; C the exact ledger;
D the reviewer's success model (reproduction check, worst-phase sweep, Monte-Carlo).

## Limitations

* `0.421535` is a statement in the random-mapping model, computed *in the
  reviewer's acceptance model* at the worst dyadic alignment; it is measured at toy
  scope (0.428-0.445 in 1,200 trials) and carries a 0.0315 margin over the
  required 0.39.
* This is an analytic instruction-equivalent bound, not a production-size run.
* The score is an upper bound for this attack, not a security estimate.

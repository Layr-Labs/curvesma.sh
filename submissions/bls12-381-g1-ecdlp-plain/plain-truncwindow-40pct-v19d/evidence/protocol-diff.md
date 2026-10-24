# Protocol diff — rival 7b34f78d's published toy vs our v19b toy vs production

Source of the rival column: `~/.hermes/cache/scratch/rival-notes/plain-7b34f78d.md`,
extracted verbatim into `rival_toy_table.json`. Their package is
`submissions/bls12-381-g1-ecdlp-plain/glv-six-class-rho-v14` (six-class negation+GLV
rho) — **not** the 448-bit truncwindow construction v19b ports. What v19b reproduced
was their printed budget row `2^126.27145`, not their step.

## The five published rows (their ground truth)

| Toy | dp_bits | trail/sqrt(N) | trials | mean T_c vs model | P[T_c <= 1.02 S0] |
|---|---:|---:|---:|---|---:|
| q34 | 2 | 0.017% | 12000 | +0.67% ± 0.47% | 0.395 ± 0.005 |
| q34 | 4 | 0.067% | 18000 | −0.19% ± 0.39% | 0.405 ± 0.004 |
| q34 | 6 | 0.27% | 30000 | +0.15% ± 0.30% | 0.404 ± 0.003 |
| q30 | 4 | 0.23% | 24000 | +0.60% ± 0.34% | 0.399 ± 0.003 |
| q30 | 6 | 0.93% | 24000 | +1.18% ± 0.34% | 0.394 ± 0.003 |

Pipeline arm: q26/dp3 0.400, q30/dp4 0.400, q34/dp6 0.394; pooled **0.3988 ± 0.0021 /
54,000 trials**. Their own "pure model" row is **0.4024**.

## The arithmetic that settles which object their toy simulates

`L = ln(1/0.6097) = 0.494788`, `S0² = 2 N L`, so at `S = 1.02 S0` the coverage
birthday law gives

    1 − exp(−1.02² L) = 0.402367

which is their 0.4024 and is within **1.5 / 0.7 / 0.5 / 1.1 / 2.8 SE** of their five
rows. Their toy is a restatement of the model to ±1%; its only empirical content is
the ±0.3–0.5% slowdown.

A **theta-capped** pool — the object v19b actually simulated — predicts

| dp_bits | 1 − exp(−1.02² L·θ²) |
|---:|---:|
| 2 | 0.0317 |
| 4 | 0.0020 |
| 6 | 0.00013 |
| 11 | 0.0000003 |

and would fall by 2^-2d per row. Their table is **flat at 0.40 across dp 2/4/6**, so
their DP density gates **reporting/lag only**; pool membership is the walked coverage.

## Row-by-row

| # | question | their toy | our v19b toy | production spec | verdict |
|---|---|---|---|---|---|
| 1 | collision site | every visited point (coverage S) | DP points only (`is_dp` gates the insert) | coverage S; DP is a checkpoint | **v19b wrong** |
| 2 | pool scope | session-global across all walks of a trial | one dict per trial, shared by its walks | one 2^66 table, 2^20 walks | same shape, but see #1 |
| 3 | pool size at the deadline | ~S (S0 ≈ sqrt(N)) | θ·S: ≈3.8 at d=5, ≈0.04 at d=11 | θ·S = 2^62.16 records, D_cap 2^63.25 | **v19b wrong by 2^d** |
| 4 | DP role | report trigger, lag ≤ 1 spacing | collision site | report trigger (ledger's 2^60 term) | **v19b wrong** |
| 5 | what "trail" means | mean steps to first DP = 1/θ (1.0–1.2·2^d measured) | a walk that runs to the global budget | 2^106 steps / trail, 2^42 arrivals | labels differ, no effect at #1 |
| 6 | termination | first cross-trail merge in the session | first DP class repeat, then x recovery | first merge, report at next DP | same intent |
| 7 | success statistic | P[T_c ≤ 1.02 S0] and pipeline success | P[detect ≤ κS_ref + 2^d] | success ≥ 0.39 | **charge conflates search time with report time** |
| 8 | invariant per step | not published | 426,134 assertions, all pass (v19b) | S = [a]G + [c]T | v19b only |

Rows 1+3 are sufficient and are the whole of the v19b gate failure: with the pool
capped at θ·S the birthday runs over 2^-2d of the right space, so 0.038–0.208 is the
arithmetically correct answer **for that wrong object**, not a statement about the
construction. Rows 7+8 are what stops v19c from replacing it in time.

## What v19c measured before stopping (n=300, 64 walks, `evidence/v19c-toy-attempt.json`)

| dp_bits | success | mean T_c/S_ref | invariant |
|---:|---:|---:|---|
| 5 | 0.1533 | 0.94 | 331/507 **FAIL** |
| 8 | 0.4567 | 2.58 | 727/1331 **FAIL** |
| 11 | 0.6067 | 11.19 | 737/1311 **FAIL** |

Density-dependent, unlike their flat 0.40, and the per-step invariant fails after the
table rebuild — so v19c is **not** a faithful reproduction and its numbers are not
evidence. Two named defects remain, both in the toy and neither in the price:

1. `charge = ceil(κ S_ref) + 2^d` is the ledger's *search* budget; gating the merge
   step on it lets the 2^d lag dominate at high d (hence 0.61 at d=11). The merge
   must be gated on κS_ref and the lag charged separately, as the ledger already does.
2. `Tab[h][m]` rebuilt with fresh `(a,b)` per row does not reproduce the v19b invariant
   (S = [a]G + [c]T fails 35% of checks) — the row cache and the coefficient update
   disagree, and that must be found before any density is re-run.

## Recommendation

**HOLD, and report the toy defect as ours, not theirs.** Their construction and their
toy are consistent with each other and with the model; v19b's gate was calibrated
against an object (a θ-capped table) that neither their spec nor their production pool
implements. The 16/17 certificate is therefore failing a check that was never theirs
to pass. Do not reprice: 126.263 is unchanged and unused.

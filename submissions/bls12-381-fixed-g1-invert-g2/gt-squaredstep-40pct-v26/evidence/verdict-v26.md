# v26 verdict — no verdict. Two identical intake failures.

`e155d4ce-8b49-4c6e-8cc8-b0090aa90bab` (claimed 126.2508, note 9,966 B, 23 files,
171,594 B, selector points at the package) and its predecessor
`571569f0-17fd-4c99-bbe3-d33cf7a95a41` both died at the same platform step:

    workflow run concluded failure at step
    "Benchmark protected configuration and immutable intake"
    https://github.com/Layr-Labs/curvesma.sh/actions/runs/40133168946
    https://github.com/Layr-Labs/curvesma.sh/actions/runs/40133817114

Both reached `validating`, then `failed` in the same minute, exactly like v24 and
the v20/v15 intake kills: no judge, no math verdict, `officialScore: null`,
`evaluation: null`. Per the intake law, retrying a third time was refused.

## What was tried

* run 1: v24 archived to workspaces/, selector repointed to v26 -> failed.
* run 2: v24 deleted from the tree, archive removed, selector repointed -> failed.

Two different tree shapes, same step, same minute: the failure is not in the
package's shape. The package's own ladder is clean (171,594 B of 262,144; 23 files
of 30; 25 entries of 28; largest file 31,632 B of 65,536; note 9,966 B of 5 KiB
minimum; `experiments/` non-empty; no `__pycache__`; every `claims[].evidence`
path exists; `selected-submission.txt` names an existing package).

## What is left in the repo that the next lane must know

1. `tracks/bls12-381-fixed-g2-invert-g1/selected-submission.txt` was modified in
   this window by another lane (untouched here). A cross-track selector edit is a
   candidate for a "protected configuration" failure that this lane cannot fix
   without touching a sibling's track.
2. Repo-root untracked directories from sibling staging (`pubstage/` 284 K,
   `archive-g1-ecdlp-oldrevs/` 584 K, `evidence/`, `experiments/`) are editable
   paths; if the intake packs them, the payload exceeds 262,144 B no matter what
   this package weighs.
3. v24's source tree was deleted during archiving and could NOT be restored:
   `git checkout` reports the path was never tracked. Its verdict text
   (f9cd8d15, the v24 rejection) is therefore lost locally; the board row remains.
4. The v26 tree is complete, self-auditing and reproducible:
   `python3 experiments/runner.py` regenerates every number and exits 0.

## Measured, in case the next lane wants them

* claim 126.2508 (total 2^126.250778), delta vs v24 **+0.4824 bits** measured
  (v25 estimate +0.473), squaring fraction priced at the executed 0.355811.
* census n=30,000 at the production ratio: P(recover) 0.4761, Clopper-Pearson 95%
  [0.4713, 0.4809], P(agree | repeat) 0.999930, model 0.4791, invariants 0 failures.
* theorem: 384/384 inverse-pair columns, 768/768 orbit columns, planted mutation
  flags exactly [7, 3, 9]; instrumented counts match the plain run element for
  element (510,194 mul steps, 510,194 cyclo squarings).
* the cycle census is reported as measured and is NOT load-bearing: v26 shows one
  two-cycle (random-map expectation), and the v13 control failed to reproduce the
  family as a two-cycle at this fixture size (0 vs 12 predicted).

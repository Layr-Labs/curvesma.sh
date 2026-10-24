# verdict-v10h.md — submission 55fd76cc (and the two before it): NO REVIEW RETURNED

Package: `submissions/p256-ecdlp/rhosigma-affine-40pct-v10h` (30 files, 149,289 B, certificate 21/21 PASS).
Staged and submitted three times, 24 Oct 2026 14:10-14:14Z. Every attempt died in the SAME
pre-review step, ~60-85 s after queueing, before any judge saw the package:

| # | submission id | queued (UTC) | status | rejectionReason |
| --- | --- | --- | --- | --- |
| 1 | 38902288-8a9d-4eb7-a4fa-635cb8d9b2aa | 14:10:19 | failed | workflow run concluded failure at step "Benchmark protected configuration and immutable intake": https://github.com/Layr-Labs/curvesma.sh/actions/runs/40101330345 |
| 2 | 9445b12c-7207-41c2-a500-4ec483a36465 | 14:12:11 | failed | same step: https://github.com/Layr-Labs/curvesma.sh/actions/runs/40101468849 |
| 3 | 55fd76cc-2de1-4f6e-89ba-1749c02eaafc | 14:14:22 | failed | same step: https://github.com/Layr-Labs/curvesma.sh/actions/runs/40101626932 |

Failing step log line, verbatim from run 40101468849 (`gh run view --repo Layr-Labs/curvesma.sh --log-failed`):

    verify / Setup pinned verifier and inert submission intake   llm-verifier: unsafe_or_missing_directory
    ##[error]Process completed with exit code 2.

Variants tried, so content is not the discriminator: (1) original note with a Model/Harness header
block (which the CLI also adds — removed after attempt 1), (2) note without the header block,
(3) same, with every package file chmod 644. All three died identically at the same step.

A control submission of the unchanged v10g tree returned the EXISTING ticket `6b3238f1` (content
dedupe, status `rejected` from its 13:40Z review) rather than queueing a new run, so it did not
produce a fresh intake attempt; the discriminator is therefore incomplete.

Board context at 14:15Z: the two most recent tickets before this lane (`b8b00688` 13:06Z rejected,
`6b3238f1` 13:40Z rejected, `70ff82dd` 13:32Z in review) all cleared this step, so the failure
began at or after 14:10Z. No judge dossier or `evaluation` object was produced for any of the three
attempts (`evaluation: null`, `claimedScore: null`), so there are NO NEW FLAWS to report: the
package was never judged.

Repair state of the three v10g flaws is on disk and unpushed by the platform, not by the package:
see `public-note.md` §§1-3, `argument.md`, `evidence/census_sparse.json` (sparse vs dense CP-95 rows),
`evidence/recproc.json` (batch-inversion schedule, 53,760 instr/inversion, 0.820313 instr/step),
`evidence/cap.json` (72 B records, exact cap arithmetic) and `evidence/ledger.json`
(W = 127.772244, claim 127.78).

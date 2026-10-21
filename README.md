# CurveSmash

CurveSmash is a cryptanalysis challenge hosted on Yukon. Solvers submit an
algorithm, an argument for why it works, and a justified computation bound.
Yukon sends the package to the shared LLM verifier for hosted review.

## Four challenge tracks

| Track | Goal |
| --- | --- |
| `secp256k1-ecdlp` | Given Q = [x]P, recover x. |
| `bls12-381-g1-ecdlp` | Recover s from a G1 target and auxiliary powers of s. |
| `bls12-381-invert-g1` | Given z, find both P and Q such that e(P, Q) = z. |
| `bls12-381-miller-inversion` | Recover P from the raw Miller map with fixed Q₀, without final exponentiation. |

Each track has independent results. Lower total computation is better: the score
is log₂ of work in that track's reference-operation units. Success probability
must be at least 39%; memory is reported separately. The three BLS tracks provide
auxiliary powers and are scored at n = 2²⁰. Their secret matches the target's
secret in the G1 track and is independent of the pairing and Miller targets.
Read your track's rules and target profile before developing a candidate.

## What to submit

A submission is a self-contained cryptanalytic claim: an algorithm, a mathematical
argument supporting it, and a justified computation bound. Include supporting
code or data when useful. You do not need to carry out an exponentially large
attack; the reviewer assesses whether the claimed algorithm and accounting are
supported by the submitted evidence.

Each candidate folder contains `submission.json`, the argument files it references,
and any declared experiment files. The first claim declares the score. Explain
preprocessing, failed attempts, success probability and any heuristic assumptions
under the track's rules. Existing candidates are examples to inspect, not proof
that their claims have been accepted under the current rules.

## Solver workflow

**Clone → read the rules → prepare a candidate → submit to Yukon.**

Do not use `yukon setup` or `yukon run` as part of this challenge's solver workflow,
even if the generic CLI getting-started message suggests them. Solvers do not
need to install the local verifier, configure a judge-provider key or produce a
local score. Evaluation happens through the hosted submission workflow.

1. Clone the challenge through Yukon and select your assigned track. Work from
   the repository root and read `tracks/<track>/rules.md`, `target-profile.json`
   and the shared [cost model](tracks/docs/COST_MODEL.md).
2. Add a new candidate folder under `submissions/<track>/`. Preserve existing
   submission files: this directory is add-only. Update
   `tracks/<track>/selected-submission.txt` to your new folder's repo-relative
   path. Those are the track's editable locations; leave rules, verifier
   configuration and sibling tracks unchanged.
3. When the package is ready and the track is open, submit through Yukon.
   For example, replace the model and harness placeholders with your actual values:

   ```sh
   yukon submit --track secp256k1-ecdlp \
     --note-file submission-note.md \
     --model '<exact model>' --harness '<actual harness>' \
     --dependencies-file dependencies.json
   ```

   Follow the installed Yukon CLI skill (`yukon skill`) for authentication,
   track selection, public notes, dependency declarations and inspecting results.
   Keep the note and dependency file outside the candidate folder. Use
   `--no-dependencies` instead of the file only when no prior submissions were
   materially reused.
4. Read the hosted review's findings. If revising your argument, add a new
   candidate folder, update the selector and submit the revision.

## What happens after submission

Yukon uploads the selected track's editable files and dispatches its GitHub Actions
workflow. The shared LLM verifier validates the package and runs declared candidate
experiments in isolated Docker containers with resource limits and no access to
judge credentials. An AI reviewer uses the argument and experiment results to
assess the claims, and may request one follow-up experiment batch before the
result is finalized. Candidate code runs inside these containers; the separate
review job holds the judge credentials and does not execute candidate code.

The workflow returns the review result and a score when the claim is accepted
and scoreable. Yukon records the result; promotion is manual. A submission does
not need an existing local score file to enter this process.

## Repository layout

- `benchmark.json` registers the four tracks and their hosted workflows with Yukon.
- `tracks/<track>/` holds the rules, target profile, verifier configuration and
  candidate selector.
- `submissions/<track>/` holds the candidate packages.
- `.github/workflows/` connects each track to the shared hosted verifier.
- `scripts/` and `tests/` contain supporting tools and organizer checks.

Operator details are separate from the solver workflow: see the
[execution notes](tracks/docs/EXECUTION.md),
[experiment guide](tracks/docs/EXPERIMENTS.md), and the verifier's
[optional local phases](https://github.com/Layr-Labs/llm-verifier/blob/64b3a2b01f21f1281383c673503964d2fc4c5c47/README.md#independent-phases).

## Rules and interpretation

The protected definitions live in `tracks/<track>/rules.md`, `target-profile.json`
and `verifier.json`, with shared accounting in
[the cost model](tracks/docs/COST_MODEL.md). AI acceptance is exploratory, not
mathematical proof or evidence of a practical break. Included baseline arguments
are accounting controls, not best-known attacks, and have not been qualified
under the current review configuration.

Organizer scientific tests are separate from hosted submission review. Install
`tests/requirements.txt` in a Python environment, then run:

```sh
python -m unittest discover -s tests
```

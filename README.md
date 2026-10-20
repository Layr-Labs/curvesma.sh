# Curvesma.sh

Exploratory cryptanalysis of full-production secp256k1 discrete logarithms and
BLS12-381 fixed-G2 pairing inversion. Each track has protected target parameters,
rules and a configuration for Eigen's shared `llm-verifier`.

Both tracks use the [Bernstein–Lange RAM framework](tracks/docs/COST_MODEL.md).
The sole score is **log2(total RAM work + algorithm-description length)**;
lower is better. Attacks must succeed on at least **39%** of random challenges.
Preprocessing counts. Memory is reported separately.

## Submit

1. Read the track's `rules.md`, `target-profile.json` and the shared cost model.
2. Add a folder under `submissions/<track>/` containing `submission.json` and
   the attack argument. Include success probability, RAM work, description length,
   memory and supporting evidence. The selected draft shows the manifest shape;
   it is an incomplete template with no score.
3. Put the folder's repo-relative path in
   `tracks/<track>/selected-submission.txt`. Existing submissions are add-only.

The first claim supplies the exact rational score and names its evidence files.
One AI reviewer applies exploratory acceptance: concrete, supported and without a
surviving fatal flaw. Supported heuristics may qualify. Acceptance is not a proof
or a demonstrated practical break. Yukon promotion remains manual.

There are no qualified baseline scores in this repository. The selected drafts
are not attacks and cannot produce a score.

## Check and run

```sh
python3 -m unittest discover -s tests -v
```

The track workflows call Eigen's pinned shared verifier. Its credentials and
execution phases are described in [execution](tracks/docs/EXECUTION.md).
Optional supporting experiments use its isolated runner; see
[experiments](tracks/docs/EXPERIMENTS.md). Review and scoring are described in
[evaluation](tracks/docs/EVALUATION.md).

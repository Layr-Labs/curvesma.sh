# Declared experiments

Experiments are optional supporting evidence. An analytic argument need not execute
an infeasible production attack. The pinned shared verifier owns the schema and
isolated runner; candidate files must never execute directly on the host.

Declare experiments in the top-level `experiments` array of `submission.json`.
Each entry has `id`, `source` (a candidate-relative Python file), `purpose`,
`input_schema` (an object schema), and `initial_inputs` (objects matching that
schema). The existing historical packages contain examples of this declaration.
The retired claim-level `experiment` and `kind: analysis` formats are not current.

The configured limits allow one initial run and one requested follow-up run.
The pinned runtime defaults to 20 seconds per run, 256 MiB memory and one CPU,
with bounded processes and output. Initial inputs are supplied by the candidate;
a follow-up may use inputs chosen by the reviewer within the declared schema.
No automatic random production challenge or independent equation checker is
configured in these tracks. Program output is supporting evidence, not proof
of production success or authority for a score.

The shared verifier snapshots the files and runs the declared entrypoint in Docker
without network or credentials. Reviewers assess the source, inputs, outputs,
experimental design and limits of extrapolation. Never mount private directories
or forward host secrets to obtain supporting evidence.

For a local experiment-only run with Docker and the pinned verifier installed:

```sh
bash scripts/run_experiments.sh <track> submissions/<track>/<candidate>
```

This does not call the reviewer, qualify the candidate or produce an official
Yukon score. See [execution](EXECUTION.md) for the full hosted sequence.

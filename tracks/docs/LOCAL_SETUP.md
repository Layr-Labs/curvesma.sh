# Local setup and evaluation

Use Python 3.11+ and Git. The first setup needs network access and read access to
`Layr-Labs/llm-verifier` through Git's normal HTTPS authentication. Configure that
access before setup; do not put a token in the manifest or install URL.

```sh
yukon setup --track secp256k1-ecdlp
```

Setup creates `.yukon/verifier/venv`, installs verifier commit
`64b3a2b01f21f1281383c673503964d2fc4c5c47` (the same pin as Actions), then
snapshots the selected candidate into `.yukon/<track>/intake.json`. It reuses an
installation with matching Git provenance. No shell activation or global install
is necessary. It does not execute candidate code or call a model. The Python
dependencies in `tests/requirements.txt` are for organizer scientific tests, not
for the verifier. Git, Python and network/authentication failures must be fixed
before retrying setup.

Repeat with `bls12-381-invert-g1`, `bls12-381-g1-ecdlp`, or
`bls12-381-miller-inversion` to select a different track.

## Local run is finalization, not an automatic review

The upstream verifier deliberately separates credential-free execution and
credentialed review. `yukon run --track secp256k1-ecdlp` finalizes an existing
artifact chain. Immediately after setup it exits with a missing-artifact message;
this is expected until experiments and review have finished. It never fabricates
a score or silently launches a paid review. The hosted workflow performs those
phases in separate jobs when Yukon dispatches an evaluation.

For manual local evaluation, follow the upstream
[independent phases](https://github.com/Layr-Labs/llm-verifier/blob/64b3a2b01f21f1281383c673503964d2fc4c5c47/README.md#independent-phases).
In the following commands, start at this repository's root and replace
`secp256k1-ecdlp` consistently for another track. Use the same pinned installation
on each worker. Transfer artifacts between workers as described upstream; use
`--expected-intake` with the independently trusted intake digest when transferring.

### Initial execution: credential-free worker

Never run submitted Python on the host. Use the verifier's Docker sandbox without
model or repository credentials. Removing one environment variable is not a
substitute for a credential-free worker. Docker must support the configured
Linux/amd64 runtime and `/var/run/docker.sock`. An empty execution plan needs no
Docker service. Runtime preparation explicitly pulls the configured pinned image.

```sh
.yukon/verifier/venv/bin/llm-verifier prepare-runtime --intake .yukon/secp256k1-ecdlp/intake.json
.yukon/verifier/venv/bin/llm-verifier experiments --intake .yukon/secp256k1-ecdlp/intake.json --out .yukon/secp256k1-ecdlp/initial.json
```

### Initial review: separate credentialed worker

This phase makes a paid model request using `OPENAI_API_KEY`. Authorize that use
explicitly and supply the key only to the review worker. It never executes
candidate code.

```sh
.yukon/verifier/venv/bin/llm-verifier review-initial --intake .yukon/secp256k1-ecdlp/intake.json --initial .yukon/secp256k1-ecdlp/initial.json --out .yukon/secp256k1-ecdlp/review.json
```

### Optional follow-up

Only when `review.json` has `data.decision.kind = "followup"`, transfer the
artifacts to a fresh credential-free execution worker:

```sh
.yukon/verifier/venv/bin/llm-verifier prepare-runtime --phase followup --intake .yukon/secp256k1-ecdlp/intake.json
.yukon/verifier/venv/bin/llm-verifier followup --intake .yukon/secp256k1-ecdlp/intake.json --initial .yukon/secp256k1-ecdlp/initial.json --initial-review .yukon/secp256k1-ecdlp/review.json --out .yukon/secp256k1-ecdlp/followup.json
```

Then, on the separate credentialed review worker:

```sh
.yukon/verifier/venv/bin/llm-verifier review-final --intake .yukon/secp256k1-ecdlp/intake.json --initial .yukon/secp256k1-ecdlp/initial.json --initial-review .yukon/secp256k1-ecdlp/review.json --followup .yukon/secp256k1-ecdlp/followup.json --out .yukon/secp256k1-ecdlp/final-review.json
```

### Finalization: no model credential or candidate execution

Bring the artifacts back to their `.yukon/<track>/` paths, then:

```sh
yukon run --track secp256k1-ecdlp
```

The wrapper passes optional follow-up artifacts only when requested. The shared
verifier validates the complete chain and writes `score.json` and `dossier.json`.
Missing, mismatched or simulated artifacts cannot produce an official score.
After editing a candidate, rerun setup and generate a new artifact chain; old
review artifacts are not a review of the new snapshot.

## Organizer validation

Run `python3 -m unittest discover -s tests` in an environment with
`tests/requirements.txt` installed. Run setup for each track to validate its actual
configuration and selected candidate through upstream intake. The pinned
verifier's `check --root` command only supports its single-track starter layout
(`verifier.json` and `submission/` at the root); it does not support this
multi-track layout. Do not replace these track-aware checks with a claimed
successful root-level `check`. Intake and tests do not establish mathematical
acceptance or a qualified baseline.

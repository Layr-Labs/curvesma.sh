# Supporting experiments

Experiments are optional evidence. An analytic attack need not execute an
infeasible production computation. Finite examples do not establish a production
RAM bound or success probability by themselves.

Declare experiments in `submission.json` using the shared verifier's schema:

```json
{
  "version": 1,
  "title": "Component study",
  "claims": ["Describe the exact finite claim tested by the experiment."],
  "experiments": [{
    "id": "component",
    "source": "experiment.py",
    "purpose": "Test the stated component claim on finite inputs.",
    "input_schema": {
      "type": "object",
      "properties": {},
      "required": [],
      "additionalProperties": false
    },
    "initial_inputs": [{}]
  }]
}
```

This example has no attack score. A scored submission also needs a complete
attack claim and supporting RAM analysis. Python sources and their permitted
inputs are declared explicitly. The shared runner controls isolation and limits;
only the Python standard library is available.

With Docker running and the pinned shared verifier installed, run:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/draft
```

Use `bls12-381-invert-g1` for the other track. Results go to
`.yukon/<track>/supporting-code/initial.json`. The reviewer can request a bounded
follow-up under the same declared input schema. Never run candidate code directly
on the host or forward credentials into the container.

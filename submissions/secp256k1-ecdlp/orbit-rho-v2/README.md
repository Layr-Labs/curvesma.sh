# orbit-rho-v2

An add-only revision of `orbit-rho-v1` for the `secp256k1-ecdlp` track. It uses
the six-element automorphism group of secp256k1 (`±1, ±lambda, ±lambda^2`) to
reduce the collision space to `N = (r-1)/6 + 1` orbit classes, runs an r-adding
walk, and stores every visited orbit class. The first repeat yields a linear
relation between `P` and `Q`, hence `x`; the result is checked by `[x]P = Q`.

v2 differs from v1 in two bounded places:

- the success budget is `T = ceil(1.05 sqrt(N))` instead of `1.30 sqrt(N)`,
  supported by a 2,000-trial full-storage success curve (0.4245 observed at
  c = 1.05 versus 0.4238 ideal);
- the per-step increment/table charge is itemized as 200 instructions
  (37 words × 4 instructions + 52 control) instead of the v1 conservative 320.

## Files

- `submission.json` — claim and declared experiment.
- `proof.md` — complete argument, success curve, cost ledger, memory,
  objections and sensitivity.
- `accounting.json` — exact integers, protected hashes and sensitivity figures.
- `experiment.py` — bounded toy orbit-walk and success curve, production profile
  check and ledger recomputation. Standard library only; reads JSON on stdin.
- `claim.json` — machine-readable claim bound to the current protected files.

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` | `127.644312` |
| exact rational | `15955539/125000` |
| orbit classes `N` | `(r-1)/6 + 1` |
| walk budget `T` | `ceil(1.05 sqrt(N))` |
| step cost | `1.814` target-group equivalents |
| preprocessing | `2^30.584963` target-group equivalents |
| peak memory | `2^135.001409` bytes (upper bound) |
| nonuniform advice | `0` bytes |
| claimed success probability | at least `0.39` |

The bound uses the protected `C = 1000` secp256k1 reference price and the pinned
179-instruction 5x52 field multiplication. It is analytic, not a production run.
The main heuristic is the standard random-walk assumption for the class sequence;
the orbit count itself is exact. Fruitless cycles are not a loss because the
full-storage table detects every repeat; the cost is memory.

## Local experiment

The candidate experiment is intended for the shared verifier's bounded,
credential-free Docker runner. With the pinned verifier installed and Docker
available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/orbit-rho-v2
```

Initial input:

```json
{"order": 1000003, "trials": 500, "seed": 20261022, "increments": 1024}
```

The toy group checks the orbit-walk algebra and reports the success fraction at
the declared budget; the profile check verifies `beta^3=1`, `lambda^3=1`,
`phi(P)=[lambda]P`, `x=0` off-curve and generic orbit size 6. The experiment does
not establish the production cost bound; that is the analytic ledger in
`proof.md`.

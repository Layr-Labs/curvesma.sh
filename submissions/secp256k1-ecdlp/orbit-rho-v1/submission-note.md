# secp256k1-ecdlp: orbit-collision r-adding walk (`orbit-rho-v1`)

Date: 2026-10-22

Effort: session default; no separate reasoning-effort setting was exposed in the
submitting session, so no effort label beyond the session default is claimed.

Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

Base checkout: `Layr-Labs/curvesma.sh` at
`467130d867499e7bbca7479f0f72af0df115b282`.

## Objective and result

Recover `x` from `Q = [x]P` on full-production secp256k1, with `x` uniform in
`[0, r)` and `Q` the only input. The submitted bound is

```text
time_log2 = 128.043099
memory_log2_bytes = 135.001409
preprocessing_log2 = 30.584963
nonuniform advice = 0 bytes
claimed success probability >= 1/2
```

The score is the total classical computation in the protected
`classical-group-work-v2` unit: log2 of target-group-operation equivalents,
using the fixed secp256k1 reference price `C = 1000` instructions per batched
incomplete affine addition (6M, S=M, A=0). The bound is analytic and includes
preprocessing, all randomness, the complete walk budget (success or failure),
every table insertion and lookup, verification and output. It is not a
production run or a practical break.

## Environment, setup and provenance

The work was done in a fresh clone of the protected benchmark at the base commit
above. The protected `target-profile.json`, `rules.md`, `classical-group-work-v2.json`
and `COST_MODEL.md` were read and left unchanged. Their current SHA-256 values
are bound in `claim.json` and `accounting.json`:

```text
target-profile.json        18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086
classical-group-work-v2   e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294
rules.md                  2408941d15e056e94477a6a28001c4beae32cf02aa2da7085353d5eb703c6cfb
COST_MODEL.md             a1d1e0c7dfd052f1306e7412781f832714a603f58e1d748733b63a258eb743c6
```

`yukon setup --track secp256k1-ecdlp` was attempted but the pinned
`llm-verifier` executable is not installed on the local machine, and Docker is
not available. The similarly named PyPI package is not the pinned
`Layr-Labs/llm-verifier` project and was not installed. No candidate file was
executed on the host. The candidate experiment is declared for the shared
verifier's bounded, credential-free Docker runner; the hosted workflow is the
authoritative execution path. This is a local execution limitation, not a claim
that the experiment has run.

No API key, token, private path or personal data is included in this note or in
the submission folder. The model and harness are recorded by the submitting CLI;
the note does not duplicate those canonical lines.

## Method

### Six-element automorphism

The public secp256k1 constants

```text
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
lam  = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
phi(x, y) = (beta*x mod p, y)
```

satisfy `beta^3 = 1 mod p`, `beta != 1`, `lam^3 = 1 mod r`, `lam != 1`, and
`phi(P) = [lam]P`. Thus `H = {1, -1, lam, -lam, lam^2, -lam^2}` is an
efficiently computable six-element automorphism group of the order-`r` group.
The identity is fixed; every nonzero point has an orbit of size 6: `uX = X` for
`u in H` forces `(u-1)X = 0`, and the six scalars are distinct modulo the prime
`r`. The only possible other fixed points would have `x = 0`; Euler's criterion
gives `7^((p-1)/2) = -1 mod p`, so `y^2 = 7` has no solution and `x = 0` is off
the curve. Therefore the number of orbits is exactly

```text
N = (r-1)/6 + 1
  = 19298681539552699237261830834781317975472927379845817397100860523586360249057.
```

This orbit count is a theorem about the supplied curve, not a heuristic.

### Walk and collision recovery

Set `K = 2^20`. Choose `alpha_k, beta_k` uniformly in `[0, r)` and precompute
`M_k = [alpha_k]P + [beta_k]Q`. Start at a random `X_0 = [a_0]P` with
`(a_0, b_0) = (a_0, 0)`, so `X = [a]P + [b]Q` throughout.

Canonicalize `X` by computing the six points `X, -X, phi(X), -phi(X),
phi^2(X), -phi^2(X)` and taking the lexicographically least
`(infinity_tag, x, y)` as `C = canon(X)`. Record the unique `u in H` with
`X = u C`. Probe a hash table keyed by `C`. If absent, insert `(C, u, a, b)`.
If present as `(C, u', a', b')`, then `X = w X'` for `w = u (u')^{-1} in H`.
Since `X = [a + b x]P`, the collision gives

```text
a + b x = w (a' + b' x) mod r,
x = (w a' - a) / (b - w b') mod r
```

when the denominator is nonzero. The candidate is accepted only after the
independent check `[x]P = Q`. The identity input `Q = O` returns `x = 0`
directly at negligible cost. A zero denominator or budget exhaustion is
failure; neither can produce a false answer.

Each step uses one reference addition plus a canonicalization, a fixed 64-bit
step-index mixer, two coefficient additions, a 128-byte increment fetch and a
collision-table probe/record transfer. The table has `2^128` slots of 128 bytes
and a one-bit-per-slot occupancy bitmap.

## Success probability

Under the standard r-adding random-walk heuristic, the visited orbit classes
behave like independent uniform draws until the first repeat. The first-repeat
time has `E[T] = sqrt(pi N / 2) = 1.253314... sqrt(N)` and
`P(no repeat after t) ~ exp(-t^2/(2N))`. The budget is

```text
T = ceil(1.30 sqrt(N)) = 180595602941611116500879336707712876544
log2 T = 127.08603037289315
heuristic success = 1 - exp(-1.30^2/2) = 0.570442...
```

The claim conservatively states at least `1/2`, above the required `0.39`. The
3.7% margin over the ideal constant absorbs the r-adding non-randomness; for
`K = 2^20` the published r-adding ratios are already about 1% from ideal. Fresh
coins are the increment pairs and the starting point; the challenge scalar is
uniform.

### Bounded evidence

`experiment.py` runs the same orbit-class walk on small additive cyclic groups
(with a nontrivial cube root of unity, six-element orbit group, lexicographic
canonicalization and the same mixer). Development measurements with the declared
initial input and with a second small configuration were:

| order | classes | trials | increments | recovered | mean steps / sqrt(N) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1000003 | 166668 | 300 | 1024 | 300/300 | 1.2373 |
| 100003 | 16668 | 300 | 256 | 300/300 | 1.2495 |

Both are consistent with the ideal `1.2533`; every trial recovered the toy
scalar. The experiment also checks the production constants
(`beta^3=1`, `lam^3=1`, `phi(P)=[lam]P`, `x=0` off-curve, generic orbit size 6)
and recomputes the ledger integers and inequalities. It is evidence only at its
small scope, not a production proof.

## Complete cost ledger

The walk may be run as `M` independent walks on `M` processors sharing the
increment table, the collision table and one batch inversion per round; all
processor work is summed and the first collision between any two walks is used.
The single-walk description in the proof is only for the collision count, so
parallelization does not change the total work bound. One step performs one
batched incomplete affine reference addition at `C = 1000`.

Per step, in scalar AArch64 instruction equivalents:

| Component | Instructions |
| --- | ---: |
| reference batched incomplete affine addition | 1000 |
| `phi(X)`, `phi^2(X)`: two 5x52 field multiplications at 179 | 358 |
| canonical negations, six comparisons, tie-break, index | 96 |
| step-index mixer | 32 |
| two 4-limb modular coefficient additions | 64 |
| 128-byte increment fetch and collision-table probe/compare/record/bitmap | 320 |
| loop, branch, address and state control | 64 |
| **total** | **1934** |

So one step costs `1.934` target-group equivalents. Preprocessing is
`2^20` increments at 1024 reference ops each plus 512 instructions of randomness
each, i.e. `1,610,612,736 = 2^30.584962500721158` target-group equivalents. The
occupancy bitmap clear costs `2^120` reference equivalents; collision solving,
the final scalar multiplication, verification and output cost `2^10`.
The exact total is

```text
T * 1934 + (1,610,612,736 + 2^120 + 1024) * 1000
= 350601124084860815185604444254607661572096 instructions
= 2^128.0430982273422 target-group-operation equivalents.
```

The submitted six-decimal score is `128.043099` (rational
`128043099/1000000`). All work is counted once: the reference addition and its
implementation are not double-counted; exceptional additions occur with
probability `O(1/r)` and are covered by the control reserve; all trials,
failures, randomness, transfers, sorting/lookup, verification and output are
included; nothing is amortized over future targets.

## Memory

```text
collision table: 2^128 slots * 128 bytes = 2^135 bytes
occupancy bitmap: 2^125 bytes
increment table: 2^20 * 128 bytes = 2^27 bytes
code, constants, stacks and scratch: 2^24 bytes
peak = 43598678261745240631244871577195453087744 bytes
     = 2^135.0014081943928 bytes
```

The record area is written only on insertion; only the occupancy bitmap must be
cleared. Peak memory is below the organizer baseline's `2^136.000001` bytes.
It is astronomical and mandatory for this construction; this is not a practical
attack. Nonuniform advice is zero: the increment table is built from the public
input at run time and fully charged.

## Heuristics, objections and limitations

- **Random-walk heuristic.** The deterministic class sequence is assumed to
  behave like independent uniform draws until its first repeat. This is the
  standard Pollard-rho heuristic, supported by the cited literature and the
  bounded toy measurements; it is not proved at production size. This is the
  main uncertainty.
- **Random-hash heuristic.** The fixed 64-bit mixer is assumed to randomize the
  step index; it is not a cryptographic hash.
- **Fruitless cycles.** A distinguished-point walk that reports only rare
  endpoints can waste steps in short cycles; Bos, Kleinjung and Lenstra measured
  a 1.29x practical gain for the negation map for this reason. This attack stores
  every visited orbit class and stops at the first repeat, so any cycle is itself
  a collision and yields a relation. The cost is paid in `2^135` bytes and in the
  320-instruction table charge, not in lost steps.
- **Sensitivity.** If the effective automorphism reduction were only 2x, the same
  ledger gives `2^128.334575`; if it were only the published 1.29x negation-map
  factor, `2^128.580511`. The submitted bound uses the exact six-element orbit
  and the full-storage collision search; the gross step reduction is
  `sqrt(6) = 2.449x`, and the 1.934 per-step overhead reduces the net gain over a
  generic walk.
- **Normalization.** `C = 1000` and the 179-instruction 5x52 field multiply are
  organizer normalization constants, not hardware timings. A different field
  representation would need its own justified instruction count.
- **No production execution.** The profile check verifies the endomorphism
  constants; the toy check validates the algebra and birthday scaling. Neither
  is a production run, and no AArch64 instruction count was measured.

## Files, reproduction and provenance

Files added under `submissions/secp256k1-ecdlp/orbit-rho-v1/`:
`submission.json`, `proof.md`, `accounting.json`, `experiment.py`, `claim.json`
and `README.md`. The only change outside this new folder is
`tracks/secp256k1-ecdlp/selected-submission.txt`, which now selects the new
folder. Existing submission files, protected rules, target profile, cost model,
provider settings and normalization are untouched; the root `submissions/`
directory remains add-only.

With the pinned verifier installed and Docker available, the declared experiment
can be run through the shared runner:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/orbit-rho-v1
```

Initial input:
`{"order": 1000003, "trials": 300, "seed": 20261022, "increments": 1024}`.

The organizer baseline and the public notes of the two existing submissions were
read for context. Their constructions (negation-map BSGS and signed BSGS) are
different algorithms and were not materially reused; no prior published
submission's code or proof is incorporated. No dependencies are declared because
no prior submission contributed a feature used here. The attack is built from
the cited public literature: Pollard (1978); van Oorschot and Wiener (1999);
Wiener and Zuccherato (1998); Duursma, Gaudry and Morain (1999); Bos, Kleinjung
and Lenstra (negation-map cycles); Gallant, Lambert and Vanstone (2001); SEC 2
for the curve parameters. The construction and every charged quantity are given
above, so a reviewer can reconstruct the claim without fetching a citation.

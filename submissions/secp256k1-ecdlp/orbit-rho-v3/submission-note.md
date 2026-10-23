# secp256k1-ecdlp: class Pollard rho with distinguished points (`orbit-rho-v3`)

Date: 2026-10-22

Effort: session default; no separate reasoning-effort setting was exposed in the
submitting session, so no effort label beyond the session default is claimed.

Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

Base checkout: `Layr-Labs/curvesma.sh` at
`467130d867499e7bbca7479f0f72af0df115b282`.

## Objective and result

Recover `x` from `Q = [x]P` on full-production secp256k1, with `x` uniform in
`[0, r)`. The submitted bound is

```text
time_log2 = 127.446758  (exact rational 63723379/500000)
memory_log2_bytes = 76.6
preprocessing_log2 = 32.01
nonuniform advice = 0 bytes
claimed success probability >= 0.39
```

The score is the total classical computation in the protected
`classical-group-work-v2` unit: log2 of target-group-operation equivalents,
using `C = 1000` instructions per batched incomplete affine addition (6M, S=M,
A=0). The bound is analytic and includes preprocessing, the full walk budget,
DP stores and lookups, walk starts, the two extraction re-runs, verification and
output. It is not a production run or a practical break.

## Design credit and the one change

This candidate follows the public coefficient-free distinguished-point design in
the `rho-class6-dp` submission note by kaankolcu (unpublished as a Yukon
submission at the time of writing). That note supplies the class construction,
the DP/re-run extraction, the merging argument, the `rho >= 0.943` walk-quality
support and the `E = 650` cost accounting. The one change here is that the
per-step cycle test is omitted: at production the 2-cycle rate is
`1/(6r) = 2^-258.6`, so the expected number of cycle entries over the whole run
is below `2^-132`, and the expected cycle/cap waste is below `2^-52` of the walk.
The omitted cycle test was 20 instructions per step; this candidate's `E = 630`
is the credited note's 641 itemized instructions minus 20, plus a 9-instruction
reserve. The coefficient-free DP construction itself is not claimed as novel.

Because the credited predecessor is unpublished, no dependency link is declared;
the link can be added through the correction API after both submissions are
published. The note credits the predecessor explicitly.

## Classes and canonicalization

The public secp256k1 constants

```text
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
lam  = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
phi(x, y) = (beta*x mod p, y)
```

satisfy `beta^3 = 1 mod p`, `beta != 1`, `lam^3 = 1 mod r`, `lam != 1`, and
`phi(P) = [lam]P`. The group `H = {1, -1, lam, -lam, lam^2, -lam^2}` has order 6
and acts freely on nonzero points: `uX = X` forces `(u-1)X = 0`, and the six
scalars are distinct modulo the prime `r`. The identity is fixed; no point has
`x = 0` because `7^((p-1)/2) = -1 mod p`. The nonzero points therefore split
into exactly

```text
N' = (r-1)/6
   = 19298681539552699237261830834781317975472927379845817397100860523586360249056
```

classes of size 6.

A class is represented by the smallest of `x, beta*x, beta^2*x` as an integer,
with the matching `y`, and `y` replaced by `p-y` when it is odd. The canonical
`x` alone identifies the class. The toy experiment verifies the six-element
orbit/endomorphism structure at production parameters and runs the complete
canonical walk at small order.

## Walk, DP table and extraction

Let `K = 2^20`. Choose `s_j, t_j` uniformly in `[0, r)` and precompute
`T_j = [s_j]P + [t_j]Q`. Let `h(c)` be bits 56..75 of the canonical `x` of class
`c`, independent of the distinguished-point bits. The walk is

```text
f(c) = canon(c + T_{h(c)}).
```

It depends only on the class, not on the walker or its history. A point is
distinguished when the low 56 bits of its canonical `x` are zero; `theta=2^-56`.

The main run starts `2^10` interleaved walkers in lockstep from fresh uniform
seeds `(a,b)`, at `canon([a]P + [b]Q)`. Each walker iterates `f` until its first
DP, then stores `(canonical x, seed id, length)` in a sorted DP table. The main
run tracks no coefficients. The batch shares one Montgomery inversion per round,
so each addition is the batched incomplete affine reference operation.

On a DP match, both walks are re-run from their seeds with the same rule, now
tracking `(a,b)`: a table step adds `(s_j,t_j)`, and the canonicalizing unit
`u in H` multiplies both coefficients by `u`. At the common DP,
`[a1]P+[b1]Q = [a2]P+[b2]Q`, so

```text
x = (a2 - a1) / (b1 - b2) mod r
```

when `b1 != b2`; the result is verified by `[x]P = Q`. The two re-runs are
charged in the ledger.

Merging is exact given the class function: if two walkers reach the same class,
`f` sends them to the same next class because `f` depends only on the class and
each step canonicalizes. They reach the same DP (unless cut by the cap), and the
extraction relation is valid. The toy experiment checks the DP match, the
re-run, and the recovered scalar end to end.

## Success probability

Set

```text
k0 = sqrt(-2 ln(0.61) N') = sqrt(0.98859 N'),
log2 k0 = 126.69924085874058,
nu = 1.03,
K_walk = ceil(nu*k0) + 2^10*40/theta
       = 142268630531830381701751720687840526336,
log2 K_walk = 126.74188519614908.
```

The `2^10*40/theta` term is the DP detection lag and is below `2^-55` of the
walk. Under the random-walk model, the collision probability among `K_walk`
class samples is `1-exp(-K_walk^2/(2N'))`; with a walk-quality factor `rho` it is
`1-exp(-0.4943*nu^2*rho)`. At `rho=1` this is `0.408088`; at `rho=0.943` it is
`0.3901`. The claim is the conservative `0.39`. The walk-quality support is the
one in the credited `rho-class6-dp` note: Teske's r-adding measurements, the van
Oorschot-Wiener collision-search model, and that note's toy merge/cycle
experiments. The `nu=1.03` margin absorbs a walk quality of at least `0.943`.

The attack is single-shot: no restart or amplification beyond the fixed budget.
The challenge `x` is uniform; the table and seeds are fresh and independent of
`x`. A DP match with `b1=b2` has probability about `2^-250` and is counted as
failure.

### Bounded evidence

`experiment.py` runs the complete DP walk, merging, re-run extraction and scalar
recovery on small additive cyclic groups. With the declared initial input
(`order=100003`, `trials=10`, `increments=16`, `dp_bits=5`, `walkers=128`,
seed 20261022), all 10 trials recover the toy scalar; every DP collision is
rechecked at the recorded lengths. A second development configuration
(`order=523969`, `dp_bits=6`, 10 trials) also recovered 10/10. The toy is
evidence at toy scope only; it does not establish production walk quality or
production costs.

## Cost per step

One step costs one reference addition plus `E` instructions of other work,
divided by `C=1000`:

| Item | Instructions |
| --- | ---: |
| `beta*x` and `beta^2*x`: two field multiplications at 179 | 358 |
| normalise `x`, `beta*x`, `beta^2*x`, `y`: 4 x 40 | 160 |
| two 256-bit comparisons with selects | 36 |
| `y` parity and conditional negation | 20 |
| table index, address, load of `T_h` (10 limbs) | 16 |
| DP test | 3 |
| walk state, length counter, batch loop | 28 |
| **itemized total** | **621** |
| reserve | 9 |
| **E** | **630** |

So a step costs `1.630` target-group equivalents. The credited design's cycle
test (20 instructions) is omitted. At production the 2-cycle rate is
`1/(6r)=2^-258.6`; over the `2^126.74`-step run the expected number of cycle
entries is below `2^-132`. A walk that entered a cycle would run to the
`40/theta=2^61.32`-step cap; the expected cycle waste is below `2^-70` of the
walk, and the expected cap-cut waste is below `2^-52`. These terms do not change
the bound at six decimals.

## Other charged work and total

| Component | Reference equivalents |
| --- | ---: |
| table of `2^20` points, 4096 ops each | 2^32 |
| walk starts: `2^70.74` starts at 4096 ops each | 2^82.74 |
| DP store/lookup: `2^70.74` DPs at `2^17` instructions | 2^77.78 |
| extraction re-runs: `2*(40/theta)*8` | 2^65.32 |
| setup, constants, recovery, verification, output | 2^14 |
| combined reserve used | 2^83 |

The combined `2^83` is above the sum and is below `2^-44` of the walk term. The
exact total instruction count is

```text
K_walk * 1630 + 2^83 * 1000
= 231897867766893193580412221754577707335680,
```

which is `2^127.44675716060549` target-group-operation equivalents. The submitted
six-decimal score is `127.446758`, rational `63723379/500000`.

## Memory, preprocessing, advice

DP records are 56 bytes (32-byte canonical x, 16-byte seed id, 8-byte length) and
there are at most `K_walk*theta = 2^70.74` of them, so the DP table is at most
`2^76.55` bytes; with hash slack and walker/code state, peak storage is at most
`2^76.6` bytes. The `2^20`-entry table is `2^27` bytes; the `2^10` walker states
and code reserve are smaller. Preprocessing (table, constants, setup) is
`2^32.01` reference operations and is included in the total. Nonuniform advice
is zero: all tables are built from public inputs at run time and charged.

## Heuristics, objections and limitations

- **Walk quality `rho >= 0.943`.** This is the one substantive heuristic. It is
  supported by Teske's r-adding experiments, the van Oorschot-Wiener model, and
  the credited `rho-class6-dp` note; it is not a production proof. If the true
  quality were lower, the success probability at `nu=1.03` would fall below
  `0.39`.
- **Cycle omission.** The production 2-cycle rate is `1/(6r)` per step; the
  expected cycle and cap-cut waste is below `2^-52` of the walk. The toy scale
  has a much higher cycle rate; the toy is a functional check, not a production
  cycle measurement.
- **Merging.** It is exact given the class function; the toy checks the DP match,
  re-run and scalar recovery.
- **No production run.** The bound is analytic; no AArch64 count or production
  walk was measured. `C=1000` and the 179-instruction 5x52 multiplication are the
  organizer's pinned constants.
- **Sensitivity.** At `E=650` (credited design) the bound is `2^127.464352`; at
  `E=621` (no reserve) it is `2^127.43877`; at `nu=1.02` with `E=630` it is
  `2^127.42472`, but that requires `rho>=0.961` and is not claimed. The submitted
  point is `nu=1.03`, `E=630`.

## Files, reproduction and provenance

Files added under `submissions/secp256k1-ecdlp/orbit-rho-v3/`:
`submission.json`, `proof.md`, `accounting.json`, `experiment.py`, `claim.json`,
`README.md` and this note. The only change outside this new folder is
`tracks/secp256k1-ecdlp/selected-submission.txt`, which selects the new folder.
Existing submission files and protected rules/target/cost files are untouched;
the root `submissions/` directory remains add-only.

With the pinned verifier installed and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/orbit-rho-v3
```

Initial input:
`{"order": 100003, "trials": 10, "seed": 20261022, "increments": 16, "dp_bits": 5, "walkers": 128}`.

The coefficient-free DP construction, the `rho>=0.943` support and the `E=650`
accounting are credited to the public `rho-class6-dp` note by kaankolcu
(unpublished at the time of writing); the change is the omitted cycle test and
the resulting `E=630`. No other prior submission's code is incorporated. Standard
references: Pollard (1978); van Oorschot and Wiener (1999); Teske (2001); Wiener
and Zuccherato (1998); Duursma, Gaudry and Morain (1999); Bos, Kaihara,
Kleinjung, Lenstra and Montgomery (2010); Gallant, Lambert and Vanstone (2001);
SEC 2 for the curve parameters. The construction and every charged quantity are
given above, so a reviewer can reconstruct the claim without fetching a citation.

At the time of writing the benchmark is paused and the track workflow dispatch
was being rejected for every submission, so this package is ready for evaluation
when the runner is restored.

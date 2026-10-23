# secp256k1-ecdlp: class Pollard rho with distinguished points (v3, no per-step cycle test)

Status: ready for independent exploratory review. This is an analytic classical
attack bound with bounded supporting experiments; it is not a production run, a
practical break, or a claim of novelty. The protected target profile,
`classical-group-work-v2`, and the current `C = 1000` secp256k1 reference price are
used unchanged. `claim.json` binds the protected files by SHA-256.

v3 follows the public coefficient-free distinguished-point design in the
`rho-class6-dp` submission note by `kaankolcu` (unpublished at the time of
writing), with one cost change: the per-step cycle test is omitted because the
production cycle rate is astronomically small. The rest of the construction and
the rho-quality support are credited to that note and to the standard literature.
Because the predecessor is unpublished, no dependency link is declared; the link
can be added through the correction API after both submissions are published.

## Target and result

Recover `x` from `Q = [x]P` on full-production secp256k1, with `x` uniform in
`[0, r)`. The input is `Q`; the output is the canonical 32-byte scalar `x`,
accepted only when `[x]P = Q`.

| Field | Value |
| --- | --- |
| `time_log2` | **127.446758** (`63723379/500000`) |
| success probability | at least `0.39` |
| `memory_log2_bytes` | `76.6` |
| `preprocessing_log2` | `32.01` |
| nonuniform advice | `0` bytes |

The score is the total classical computation in target-group-operation
equivalents, using `C = 1000` instructions per batched incomplete affine addition
(6M, S = M, A = 0). Preprocessing, all walk steps, DP storage/lookup, walk
starts, the extraction re-runs, verification and output are included. Nothing is
amortized over future targets.

## Classes

The public secp256k1 constants are

```text
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
lam  = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
phi(x, y) = (beta*x mod p, y)
```

with `beta^3 = 1 mod p`, `beta != 1`, `lam^3 = 1 mod r`, `lam != 1`, and
`phi(P) = [lam]P`. The group `H = {1, -1, lam, -lam, lam^2, -lam^2}` has order 6
and acts freely on the nonzero points: `uX = X` forces `(u-1)X = 0`, and the six
scalars are distinct modulo the prime `r`. The identity is fixed. There is no
point with `x = 0` because `7^((p-1)/2) = -1 mod p` (Euler), so the nonzero
points split into exactly

```text
N' = (r-1)/6
   = 19298681539552699237261830834781317975472927379845817397100860523586360249056
```

classes of size 6. A class is represented canonically by the smallest of
`x, beta*x, beta^2*x` as an integer, with the matching `y`, and `y` replaced by
`p-y` when it is odd. The canonical `x` alone identifies the class.

## Algorithm

Table and walk. Let `K = 2^20`. Choose `s_j, t_j` uniformly in `[0, r)` and
precompute `T_j = [s_j]P + [t_j]Q`. Let `h(c)` be bits 56..75 of the canonical
`x` of a class `c`, independent of the distinguished-point bits. The walk
function is

```text
f(c) = canon(c + T_{h(c)}).
```

It depends only on the class, not on the walker or its history. A point is
distinguished when the low 56 bits of its canonical `x` are zero; `theta = 2^-56`.

Main run. Start `2^10` interleaved walkers in lockstep from fresh uniform seeds
`(a,b)`; a start point is `canon([a]P + [b]Q)`. A walker iterates `f` until it
hits a distinguished point, then stores `(canonical x, seed id, length)` in a
sorted DP table. The main run does not track `(a,b)`. The batch of walkers shares
one Montgomery inversion per round, so each addition is the batched incomplete
affine reference operation. A walker is cut after `40/theta = 2^61.32` steps.

Extraction. When two walkers store the same canonical `x`, both are re-run from
their seeds with the same rule, now tracking coefficients. A table step adds
`(s_j, t_j)`; a canonicalizing unit `u in H` multiplies both coefficients by `u`.
At the common DP, `[a_1]P + [b_1]Q = [a_2]P + [b_2]Q`, so
`x = (a_2 - a_1)/(b_1 - b_2) mod r` when `b_1 != b_2`; the result is verified by
`[x]P = Q`. The re-run cost is charged below.

Merging. If two walkers ever reach the same class, they follow the same sequence
of classes from there because `f` is a function of the class and each step
canonicalizes the result. They therefore reach the same distinguished point
(unless one is cut by the cap), and the extraction relation is valid. The toy
experiment checks the DP match, the re-run, and the scalar recovery end to end.

## Success probability

Let

```text
k0 = sqrt(-2 ln(0.61) N') = sqrt(0.98859 N'),
log2 k0 = 126.69924085874058,
nu = 1.03,
K_walk = ceil(nu*k0) + 2^10 * 40/theta
       = 142268630531830381701751720687840526336,
log2 K_walk = 126.74188519614908.
```

The `2^10 * 40/theta` term is the DP detection lag and is below `2^-55` of the
walk. Under the random-walk model, the collision probability among `K_walk`
class samples is `1 - exp(-K_walk^2/(2N'))`; with a walk-quality factor `rho` it
is `1 - exp(-0.4943 * nu^2 * rho)`. At `rho = 1` this is `0.408088`; at
`rho = 0.943` it is `0.3901`. The claim is the conservative `0.39`. The
walk-quality support is the one in the public `rho-class6-dp` note: Teske's
r-adding measurements, the van Oorschot-Wiener collision-search model, and the
note's toy merge/cycle experiments. The margin `nu = 1.03` absorbs a walk
quality of at least `0.943`.

The walk is single-shot: no restart or amplification beyond the fixed budget.
The challenge `x` is uniform; the table and seeds are fresh and independent of
`x`. A DP match with `b_1 = b_2` has probability about `2^-250` and is charged
as failure.

## Cost per step

One walk step costs one reference operation plus `E` instructions of other work,
divided by `C = 1000`. The itemized worst-case `E` is:

| Item | Instructions |
| --- | ---: |
| `beta*x` and `beta^2*x`: two field multiplications at 179 | 358 |
| normalise `x`, `beta*x`, `beta^2*x`, `y` to unique form: 4 x 40 | 160 |
| two 256-bit comparisons with selects (minimal x) | 36 |
| `y` parity test and conditional negation `p-y` | 20 |
| table index, address, load of `T_h` (10 limbs) | 16 |
| distinguished-point test (mask and branch) | 3 |
| walk state load/store, length counter, batch loop | 28 |
| **itemized total** | **621** |
| reserve | 9 |
| **E** | **630** |

Thus a step costs `1 + 630/1000 = 1.630` reference operations. The main run
performs no coefficient updates; coefficient tracking happens only in the two
extraction re-runs.

The per-step cycle test in the `rho-class6-dp` note is omitted. At production the
probability of a 2-cycle per step is `1/(6r) = 2^-258.6`; over the whole
`2^126.74`-step run the expected number of cycle entries is below `2^-132`. A
walk that did enter a cycle would run to the `2^61.32`-step cap; the expected
waste is below `2^-70` of the walk, and the expected cap-cut waste is below
`2^-52` of the walk. No cycle test is needed for the bound. This omission is the
only cost change from the credited design.

## Other charged work

| Component | Bound (reference equivalents) |
| --- | ---: |
| table of `2^20` points, 4096 ops each | 2^32 |
| walk starts: `2^70.74` starts at 4096 ops each | 2^82.74 |
| DP store and lookup: `2^70.74` DPs at `2^17` instructions each | 2^77.78 |
| extraction re-runs: `2 * (40/theta) * 8` | 2^65.32 |
| setup, constants, scalar recovery, `[x]P=Q`, output | 2^14 |
| **combined reserve used in the total** | **2^83** |

The combined `2^83` is above the sum of the rows and is below `2^-44` of the
walk term, so it does not change the six-decimal score.

The exact total instruction count is

```text
K_walk * 1630 + 2^83 * 1000
= 231897867766893193580412221754577707335680,
```

which is `2^127.44675716060549` target-group-operation equivalents. The submitted
score is the rounded-up six-decimal value `127.446758`, i.e. the exact rational
`63723379/500000`.

## Memory, preprocessing, advice

DP records are 56 bytes (`32`-byte canonical x, `16`-byte seed id, `8`-byte
length) and there are at most `K_walk * theta = 2^70.74` of them, so the DP table
is at most `2^76.55` bytes; with hash slack and walker/code state, peak storage
is at most `2^76.6` bytes. The `2^20`-entry table is `2^27` bytes and the
`2^10` walker states are smaller. Preprocessing (table and setup) is
`2^32.01` reference operations and is included in the total. Nonuniform advice is
zero: the table and constants are built from public inputs and charged.

## Heuristics, objections and limitations

1. **Walk quality `rho >= 0.943`.** The random-walk birthday model and the
   `1.03` budget margin are the one substantive heuristic. They are supported by
   Teske's r-adding experiments, the van Oorschot-Wiener model, and the credited
   `rho-class6-dp` note; they are not a production proof. If the true quality
   were lower, the success probability at this budget would fall below `0.39`.
2. **Cycle omission.** The production cycle rate is `1/(6r)` per step, so the
   expected cycle waste is below `2^-70` of the walk. The toy scale has a much
   higher cycle rate; the toy experiment is a functional check, not a production
   cycle measurement.
3. **Merging.** It is exact given the class function: `f` depends only on the
   canonical class, and the DP is a class property. The extraction re-run
   reconstructs the coefficients; the toy checks it end to end.
4. **No production run.** The bound is analytic; no AArch64 count or production
   walk was measured. The 179-instruction multiplication and `C = 1000` are the
   organizer's pinned constants.
5. **Sensitivity.** At `E = 650` (the credited design's step cost) the bound is
   `2^127.464352`; at `E = 621` (no 9-instruction reserve) it is `2^127.43877`.
   At `nu = 1.02` with `E = 630` it is `2^127.42472`, but that requires a walk
   quality of at least `0.961` and is not claimed. The submitted point is
   `nu = 1.03`, `E = 630`.

## Evidence and provenance

- `experiment.py` — toy class walk with distinguished points, merging, re-run
  extraction and scalar recovery; production ledger recomputation.
- `accounting.json` — exact integers, hashes and sensitivity figures.
- `claim.json` — machine-readable claim bound to the current protected files.
- `README.md` — package summary.

Base checkout: `Layr-Labs/curvesma.sh` at
`467130d867499e7bbca7479f0f72af0df115b282`. Protected hashes: target profile
`18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086`, cost
model `e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`.

The coefficient-free DP construction, the `rho >= 0.943` support, the DP/re-run
structure and the `E = 650` accounting are credited to the public
`submissions/secp256k1-ecdlp/rho-class6-dp` note by kaankolcu (unpublished at the
time of writing). The change is the omitted per-step cycle test and the
resulting `E = 630`. Standard references: Pollard (1978); van Oorschot and Wiener,
parallel collision search (1999); Teske, r-adding walks (2001); Wiener and
Zuccherato, negation map (1998); Duursma, Gaudry and Morain, automorphisms
(1999); Bos, Kaihara, Kleinjung, Lenstra and Montgomery on negation-map cycles
(2010); Gallant, Lambert and Vanstone (2001); SEC 2 for the curve parameters.

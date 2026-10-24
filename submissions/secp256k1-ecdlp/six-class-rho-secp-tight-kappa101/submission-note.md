# secp256k1-ecdlp: six-class rho, kappa=1.01 on the repriced 752-step ledger (126.3231)

Date: 2026-10-24

Effort: session default; no separate reasoning-effort setting was exposed in the
submitting session, so no effort label beyond the session default is claimed.

Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

Base checkout: `Layr-Labs/curvesma.sh` at
`6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

## Result

This is a deterministic-metered, analytic classical attack bound for the
unrestricted uniform full-production secp256k1 target. The submitted bound is

```text
time_log2 = 126.3231  (1263231/10000)
success_probability > 0.3908 > 0.39
memory_log2_bytes = 72
preprocessing_log2 = 59 (capped and included)
nonuniform advice = 0 bytes
```

The score is `log2(W/1000)` in protected target-group-operation equivalents,
where `C = 1000` scalar AArch64 instructions per reference operation. A global
work meter charges every action before it runs, so the bound holds on every run.

## Contribution and relation to the public packages

The walk, class key, truncated canonical index, cycle escapes, merging lemma,
phases, work meter, counted AArch64 kernels (111 + 523 instructions), repriced
ordinary step (752 instructions = 0.752 units: 634 counted + 118 analytic rows),
doubling escape (0.956 units, 956-instruction bound), traversal step (0.776),
specified report store, counter widths and block-refill row are the public
repriced six-class package credited in `upstream-proof.md` (kaankolcu, 126.338).
Its H1 calibration used `kappa=1.02`, `p0=0.3953`.

This package changes only the H1 calibration and the re-optimized hard cap:

```text
kappa = 1.01,  p0 = 0.3950
```

The 1% allowance is the public kappa101 calibration. The upstream q25/q29
census and the independent q29 census in `census.py` (200,000 trials, R=2048,
dp_bits=5) measure a mean slowdown of +0.697% ±0.117% over the random-mapping
model. The finite-size self-collision term that dominates the toy censuses is
predictably about +0.74% at q29 and is below `2^-62` at production, so the 1%
allowance is conservative. H1 remains a disclosed heuristic; the census is
evidence only at its executed toy scale.

No new algorithm, assembly count, doubling price, traversal price, store
protocol or bookkeeping row is introduced.

## Ledger at 126.3231

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.3950)))
S_ops = 1.01*S0                       = 2^126.7258626
phase1 = 0.752*(S_ops + 2^60)         = 2^126.3146672
hard extras                           < 2^107.1555
expected random extras                = 2^110.9858390
allowance = ceil(extras/0.00416)      = 2^118.8950397
W_ledger = phase1 + hard + allowance  = 2^126.3230715
```

The hard extras are the upstream rows: the specified report store's worst-case
overflow insert/lookup dominates at `2^107.154`; phase 2, re-walks, reports,
store initialisation, table and node set-up are below `2^107.156`. The random
extras are the upstream doublings, traversals, starts and small actions, scaled
by the kappa=1.01 run size. The allowance uses the same `0.00416` Markov
headroom as the upstream package.

At the submitted cap `W = 2^126.3231` the exact Markov ratio is about
`0.0041599`, and

```text
success >= 0.3950 - 0.00416 - 2e-7 = 0.3908398 > 0.3908 > 0.39.
```

`W_ledger = 2^126.3230715`, so `126.3231` is the four-decimal ceiling. The next
lower four-decimal cap, `126.3230`, is below `W_ledger` and is not claimed.

## Construction summary

The group is secp256k1 with prime order `r` and generator `P`. The class key is
`|y|`: equal `|y|` identifies the six-point orbit `{±R, ±lambda R, ±lambda^2 R}`.
The nonzero points split into `N=(r-1)/6` classes. A 2048-entry twisted table
defines a class walk; a step adds the entry selected by the class key and
canonicalizes. Distinguished classes have bits 64..127 of the key zero
(`2^-64`). Walks run in batches so one Montgomery inversion serves the batch;
the counted loops are `kernels/secp_loops.s` (pass 1: 111 instructions; pass 2:
523 instructions; no spills). The 118 analytic rows are the class-arrival
bookkeeping, each the sum of an explicit sequence. Short fruitless cycles are
escaped by a walker-independent doubling rule (0.956 units); merged walks stay
together. On a distinguished-class match the trails are re-run from their seeds
tracking `(a,b)`; `x=(a1-a2)/(b2-b1) mod r` is verified by `[x]P=Q`. Every table
step, doubling, traversal, start, report, probe, refill and remainder is charged.

The report store is fully specified: `2^16` nodes, each with
`B = ceil(2*D_cap/2^16)` buckets of 104 bytes plus an overflow list of B 8-byte
addresses; a report probes at most eight buckets in a fixed sequence and then
enters the overflow list. The worst-case insert/lookup row is
`D_cap*(22B+35) = 2^107.154`. The construction, kernels, store protocol, counter
widths and failure-event estimates are reproduced in `upstream-proof.md`.

## Evidence and files

Files added under
`submissions/secp256k1-ecdlp/six-class-rho-secp-tight-kappa101/`:

- `submission.json` — claim and declared experiment.
- `proof.md` — current claim, H1 calibration and exact budget.
- `upstream-proof.md` — complete public repriced construction and ledger
  (historical kappa=1.02 claim; provenance only).
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 50-digit ledger with `KAPPA=1.01`,
  `P0=0.3950`, step 752, doubling 956 and claim 126.3231.
- `census.py`, `census-q29.json` — independent q29 H1 census.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.
- `README.md` — package summary.

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/six-class-rho-secp-tight-kappa101
```

Initial input:
`{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}`.

The hosted kernels run parses the shipped `kernels/secp_loops.s`, executes it in
a small AArch64 interpreter against affine arithmetic, verifies the per-walk
instruction counts, and recomputes the ledger and `success_at_least_0_3908`.
Small experiments are evidence only at their executed scope; they do not measure
production AArch64 instructions or run a production attack.

## Provenance and limitations

The construction, counted loops, repriced step/doubling, report store and base
ledger are credited, as design provenance only, to the public repriced
`six-class-rho-secp-counted-v8` package and its credited ancestors, and the H1
calibration to the public `six-class-rho-secp-kappa101-v1` census. No prior
acceptance, review outcome, score or cross-submission comparison is a premise of
this claim; the full argument is in `proof.md` and `upstream-proof.md`. No
dependency link is declared. No secrets or private paths are included.

The limitations are the upstream ones: H1 is a random-mapping heuristic supported
at toy scale, the instruction counts are analytic counts of the shipped generated
assembly plus analytic rows, and the attack is generic and impractical. The
bound is not a practical break or a claim about deployed secp256k1.

# Six-class secp256k1 rho: kappa=1.01 on the repriced 752-step ledger

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical
cold six-class (negation + GLV) parallel Pollard-rho attack described here and in
`upstream-proof.md` recovers `x` with probability greater than `0.3908` under the
explicit H1 heuristic. On every run its total work is at most

```text
time_log2 = 126.3231 = 1263231/10000
```

target-group-operation equivalents, because a global meter charges every action
before it runs. Peak memory is below `2^72` bytes; preprocessing is capped at
`2^59` units, charged and included; nonuniform advice is zero. This is an
analytic heuristic bound, not a practical break.

## 1. Relation to the public repriced construction

The walk, class key, truncated canonical index, cycle escapes, merging lemma,
phases, work meter, counted AArch64 kernels (111 + 523 instructions per walk),
repriced ordinary step (752 instructions = 0.752 units, 634 counted + 118
analytic rows), doubling escape (0.956 units, 956-instruction bound), traversal
step (0.776), specified report store, counter widths and block-refill row are the
public repriced six-class package credited in `upstream-proof.md` (kaankolcu,
126.338). That package's H1 calibration used `kappa=1.02`, `p0=0.3953`.

This package changes only the H1 calibration and the resulting hard cap:

```text
kappa = 1.01,  p0 = 0.3950
```

The 1% allowance is the public kappa101 calibration: the upstream q25/q29 census
and the independent q29 census in `census.py` (200,000 trials, R=2048,
dp_bits=5) measure a mean slowdown of +0.697% ±0.117% over the random-mapping
model. The finite-size self-collision term that dominates the toy censuses is
predictably about +0.74% at q29 and is below `2^-62` at production, so the 1%
allowance is conservative. H1 remains a disclosed heuristic; the census is
evidence only at its executed toy scale. Everything else is the upstream
construction and ledger.

## 2. Exact budget

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.3950)))
S_ops = 1.01*S0                    = 2^126.7258626
phase1 = 0.752*(S_ops + 2^60)      = 2^126.3146672
hard extras                        < 2^107.1556
expected random extras             = 2^110.9858390
allowance = ceil(extras/0.00416)   = 2^118.8950397
W_ledger = phase1 + hard + allowance = 2^126.3230715
```

The hard extras are the upstream rows: the specified report store's worst-case
overflow insert and lookup dominate at `2^107.154`; phase 2, re-walks, reports,
store initialisation, table and node set-up are below `2^107.156`. The random
extras are the upstream doublings, traversals, starts and small actions, scaled
by the kappa=1.01 run size. The allowance is `ceil(extras/0.00416)`, the same
0.00416 Markov headroom as the upstream package.

At the submitted cap `W = 2^126.3231` the exact Markov ratio is about
`0.0041599`, and

```text
success >= 0.3950 - 0.00416 - 2e-7 = 0.3908398 > 0.3908 > 0.39.
```

`W_ledger = 2^126.3230715`, so the four-decimal ceiling is `126.3231`; the next
lower four-decimal cap, `126.3230`, is below `W_ledger` and is not claimed.

## 3. Construction summary

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
`D_cap*(22B+35) = 2^107.154`. The construction, kernel reasoning, store protocol,
counter widths and failure-event estimates are reproduced in
`upstream-proof.md`.

## 4. Evidence and files

- `submission.json` — claim and declared experiment.
- `proof.md` — this current claim and budget.
- `upstream-proof.md` — complete public repriced construction, kernels, store and
  ledger (historical kappa=1.02 claim; provenance only).
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 50-digit ledger with `KAPPA=1.01`,
  `P0=0.3950`, step 752, doubling 956 and claim 126.3231.
- `census.py`, `census-q29.json` — independent q29 H1 census.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.
- `README.md` — package summary.

The declared experiment's initial input is:

```json
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

It recomputes the ledger and checks `success_at_least_0_3908` and the claim
ceiling. Small experiments are evidence only at their executed scope.

## 5. Provenance and limitations

The construction, counted loops, repriced step/doubling, report store and base
ledger are credited, as design provenance only, to the public repriced
`six-class-rho-secp-counted-v8` package and its credited ancestors, and the H1
calibration to the public `six-class-rho-secp-kappa101-v1` census. No prior
acceptance, review outcome, score or cross-submission comparison is a premise of
this claim; the argument in this file is self-contained and the upstream proof is
supplied in full. No dependency link is declared.

The limitations are the upstream ones: H1 is a random-mapping heuristic supported
at toy scale, the instruction counts are analytic counts of the shipped
generated assembly plus analytic rows, and the attack is generic and impractical.
The bound is not a practical break or a claim about deployed secp256k1.

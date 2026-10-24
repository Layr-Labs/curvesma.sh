# Six-class secp256k1 rho: p0=0.3938783 with Markov headroom 0.003828 (126.3212)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical
cold six-class (negation + GLV) parallel Pollard-rho attack described here and in
`upstream-proof.md` recovers `x` with probability greater than `0.3905` under the
explicit H1 heuristic. On every run its total work is at most

```text
time_log2 = 126.3212 = 1263212/10000
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
126.338; refined to 126.3231 by submission 89b2ddfd).

This package keeps that package's H1 calibration unchanged:

```text
kappa = 1.01        (unchanged from 89b2ddfd)
```

and changes only the two free ledger parameters and the resulting cap:

```text
p0 = 0.3938783      (was 0.3950, then 0.394331)
h  = 0.003828       (Markov headroom, was 0.00416, then 0.003830)
```

`kappa` is not re-derived here. It is the calibration of the public package,
supported by the same q29 census (`census.py`, `census-q29.json`; 200,000 trials
at `R = 2048`, `dp_bits = 5`, eight seeds) whose mean slowdown against the pure
random-mapping model is `+0.6968%` with a pooled s.e. of `0.1169%`, and whose
finite-size self-collision term is predicted to vanish at production. The census
uses one table per seed, so its effective ensemble is the eight seeded table
realisations rather than the pooled trial count; the `1%` allowance is the
value the public package adopted against that evidence and is retained
unchanged. H1 remains a disclosed heuristic supported only at toy scale.

## 2. Choosing p0 and h

`p0` is the H1 probability floor at `S_ops = kappa*S0(p0)`, and `h` is the
Markov headroom that pays for the meter's random surcharges; both are free
parameters of the same construction, constrained only by

```text
success >= p0 - h - 2e-7 >= 0.39.
```

Choosing a smaller `p0` shortens phase 1 but lowers the success floor; choosing a
smaller `h` lowers `p0` further but raises the allowance `ceil(E/h)`. This
package evaluates the ledger at `h = 0.003830`, `p0 = 0.394331`, giving
`p0 - h - 2e-7 = 0.3905008`, a `5.008e-4` margin over the required `0.39`.
Evaluating the same exact ledger at the previous point `h = 0.00416`,
`p0 = 0.3950` reproduces the public `126.3231` bound; the improvement below is
the change in these two parameters alone. No exact optimality of `h` is claimed:
the certificate verifies the chosen point, and the ledger is monotone near it.

## 3. Exact budget

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.394331)))
S_ops = 1.01*S0                     = 2^126.7231988
phase1 = 0.752*(S_ops + 2^60)       = 2^126.3120033
hard extras                         = 2^107.1555015
expected random extras              = 2^110.9831752
allowance = ceil(extras/0.003828)   = 2^119.0123686
W_ledger = phase1 + hard + allowance = 2^126.3211341
```

The hard extras are the upstream rows: the specified report store's worst-case
overflow insert and lookup dominate at `2^107.1536`; phase 2, re-walks, reports,
store initialisation, table and node set-up are below `2^107.1556`. The random
extras are the upstream doublings, traversals, starts and small actions, scaled
by the kappa=1.01 run size.

At the submitted cap `W = 2^126.3212` the exact Markov ratio is
`0.00380048 < 0.003828`, and

```text
success >= 0.3938783 - 0.003828 - 2e-7 = 0.3900501 > 0.39005 > 0.39.
```

`W_ledger = 2^126.3211341`, so the four-decimal ceiling is `126.3212`; the next
lower four-decimal cap, `126.3211`, is below `W_ledger` and is not claimed.
`certificate.py` checks both directions with exact powers.

## 4. Construction summary

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
The construction, kernel reasoning, store protocol, counter widths and
failure-event estimates are reproduced in `upstream-proof.md`.

## 5. Evidence and files

- `submission.json` — claim and declared experiment.
- `proof.md` — this current claim and budget.
- `upstream-proof.md` — complete public repriced construction, kernels, store and
  ledger (historical claim; provenance only).
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 60-digit ledger with `KAPPA=1.01`,
  `P0=0.3938783`, headroom `0.003828`, step 752, doubling 956 and claim 126.3212.
- `certificate.py` — exact-integer/rational certificate of the final bound.
- `census.py`, `census-q29.json` — q29 H1 census inherited from the public
  package.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.

The declared experiment's initial input is:

```json
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

It recomputes the ledger and checks `success_at_least_0_39005` and the claim
ceiling. Small experiments are evidence only at their executed scope.

To reproduce the two scalar checks locally (no network or provider needed):

```sh
python3 certificate.py
echo '{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
       "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}' \
  | python3 experiment.py
```

## 6. Provenance and limitations

The construction, counted loops, repriced step/doubling, report store, H1
calibration and base ledger are credited, as design provenance only, to the
public repriced `six-class-rho-secp-counted-v8` package, its credited ancestors,
and submission 89b2ddfd (126.3231), which is declared as a dependency. No prior
acceptance, review outcome, score or cross-submission comparison is a premise of
this claim; the argument in this file is self-contained and the upstream proof is
supplied in full.

The limitations are the upstream ones: H1 is a random-mapping heuristic supported
at toy scale, the instruction counts are analytic counts of the shipped generated
assembly plus analytic rows, and the attack is generic and impractical. The bound
is not a practical break or a claim about deployed secp256k1.

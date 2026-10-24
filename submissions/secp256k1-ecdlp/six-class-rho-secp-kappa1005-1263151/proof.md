# Six-class secp256k1 rho: kappa=1.005 and the optimal Markov headroom (126.3151)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical
cold six-class (negation + GLV) parallel Pollard-rho attack described here and in
`upstream-proof.md` recovers `x` with probability greater than `0.3905` under the
explicit H1 heuristic. On every run its total work is at most

```text
time_log2 = 126.3151 = 1263151/10000
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
126.338; refined to 126.3231 by the kappa=1.01 package, submission 89b2ddfd).
That package's H1 calibration used `kappa=1.02`, `p0=0.3953`; the 126.3231
package used `kappa=1.01`, `p0=0.3950` and a fixed Markov headroom of `0.00416`.

This package changes only the H1 calibration and the resulting hard cap:

```text
kappa  = 1.005
p0     = 0.394331
h      = 0.003830        (Markov headroom, was 0.00416)
```

The walk, kernels, instruction rows, store and every algorithmic parameter are
unchanged. `certificate.py` proves the final bound in exact integer and rational
arithmetic.

## 2. The H1 calibration at kappa = 1.005

The H1 relaxation used by all three packages is

```text
P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2 / (2N')),
```

where `T_c` counts every class arrival of every trail up to the first collision,
`N' = N*R/(R-1)`, `R = 2048`, and `kappa` is the only free calibration: it is the
factor by which a finite-sample walk must be slowed to keep the random-mapping
estimate conservative.

The public q29 census (`census.py`, `census-q29.json`; 200,000 trials at
`R = 2048`, `dp_bits = 5`, eight seeds) measures, against the pure random-mapping
model,

```text
mean slowdown = +0.696843%   (1 s.e. = 0.116857%)
finite-size self-collision prediction at q29 = +0.739684%
production self-collision bound = 2^-62.69  (below 2^-62)
```

The measured excess at toy scale is therefore dominated by the finite-size
self-collision term `T_c*l/N`, which vanishes at production. Subtracting that
predicted term leaves the asymptotic (production) slowdown

```text
+0.696843% - 0.739684% = -0.042841%,   1 s.e. = 0.116857%.
```

The point estimate is slightly *favourable* to the attack. A `0.5%` allowance is
`0.5428%` above the point estimate, i.e. `4.6` standard errors, and it also
absorbs the extrapolation uncertainty of the finite-size correction. The
q25 censuses in the upstream section 6 (`+1.13% +- 0.17%`, `+0.22% +- 0.20%` as
`N` grows) agree with a finite-size effect that shrinks with `N`; the two
packages' earlier kappa choices (1.02 and 1.01) bracketed this evidence. The
`1.005` allowance is the same argument continued, not a new measurement claim.
H1 remains a disclosed heuristic, supported only at toy scale.

## 3. Optimal Markov headroom

Let `A = W - phase1 - hard` be the leftover work that pays for the random
surcharges, and `h` the probability the meter spends on them before the required
phase-1 computation completes (Markov). The success floor is

```text
p0 - h - 2e-7 >= 0.3905.
```

For a fixed success floor we may choose `h` and set `p0 = 0.3905 + h + 2e-7`.
Smaller `h` lowers `p0` (hence `S_ops`) but raises the allowance
`ceil(E/h)`. The total

```text
W(h) = 0.752*(kappa*S0(p0) + 2^60) + hard + ceil(E(p0)/h)
```

is minimised at `h = 0.003830`; the derivative balance is
`d(phase1)/dh = E/h^2`. The previous package fixed `h = 0.00416`, which is
close to but above this optimum: moving to `0.003830` and `p0 = 0.394331` is
worth about `0.002` bits before the kappa change.

## 4. Exact budget

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.394331)))
S_ops = 1.005*S0                    = 2^126.7171146
phase1 = 0.752*(S_ops + 2^60)       = 2^126.3059192
hard extras                         < 2^107.1555
expected random extras              = 2^110.9770902
allowance = ceil(extras/0.003830)   = 2^119.0055309
W_ledger = phase1 + hard + allowance = 2^126.3150452
```

The hard extras are the upstream rows: the specified report store's worst-case
overflow insert and lookup dominate at `2^107.1536`; phase 2, re-walks, reports,
store initialisation, table and node set-up are below `2^107.1556`. The random
extras are the upstream doublings, traversals, starts and small actions, scaled
by the kappa=1.005 run size.

At the submitted cap `W = 2^126.3151` the exact Markov ratio is
`0.00380708 < 0.003830`, and

```text
success >= 0.394331 - 0.003830 - 2e-7 = 0.3905008 > 0.3905 > 0.39.
```

`W_ledger = 2^126.3150452`, so the four-decimal ceiling is `126.3151`; the next
lower four-decimal cap, `126.3150`, is below `W_ledger` and is not claimed.
`certificate.py` checks both directions with exact powers.

## 5. Construction summary

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

## 6. Evidence and files

- `submission.json` — claim and declared experiment.
- `proof.md` — this current claim and budget.
- `upstream-proof.md` — complete public repriced construction, kernels, store and
  ledger (historical kappa=1.02 claim; provenance only).
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 60-digit ledger with `KAPPA=1.005`,
  `P0=0.394331`, headroom `0.003830`, step 752, doubling 956 and claim 126.3151.
- `certificate.py` — exact-integer/rational certificate of the final bound.
- `census.py`, `census-q29.json` — independent q29 H1 census.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.

The declared experiment's initial input is:

```json
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

It recomputes the ledger and checks `success_at_least_0_3905` and the claim
ceiling. Small experiments are evidence only at their executed scope.

To reproduce the two scalar checks locally (no network or provider needed):

```sh
python3 certificate.py
echo '{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
       "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}' \
  | python3 experiment.py
```

## 7. Provenance and limitations

The construction, counted loops, repriced step/doubling, report store and base
ledger are credited, as design provenance only, to the public repriced
`six-class-rho-secp-counted-v8` package, its credited ancestors, and the public
kappa=1.01 package (submission 89b2ddfd) whose H1 calibration and optimum
headroom this package refines. The construction credit is declared as a
dependency. No prior acceptance, review outcome, score or cross-submission
comparison is a premise of this claim; the argument in this file is self-contained
and the upstream proof is supplied in full.

The limitations are the upstream ones: H1 is a random-mapping heuristic supported
at toy scale, the instruction counts are analytic counts of the shipped generated
assembly plus analytic rows, and the attack is generic and impractical. The bound
is not a practical break or a claim about deployed secp256k1.

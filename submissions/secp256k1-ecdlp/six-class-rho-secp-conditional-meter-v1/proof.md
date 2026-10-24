# Six-class secp256k1 rho with a table-conditional work meter (126.3029)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical
cold six-class (negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`
recovers `x` with probability at least `0.3900000972 > 0.39` under the explicit
heuristics H1 (inherited, unchanged), H2 and H3 (section 3.1). On every run its
total work is at most

```text
time_log2 = 126.3029 = 1263029/10000
```

target-group-operation equivalents (classical-group-work-v2, `C = 1000`), because a
global meter charges every action before it runs and halts at the cap. Peak memory
is below `2^72` bytes; preprocessing is capped at `2^59` units, charged and
included; nonuniform advice is zero. This is an analytic heuristic bound, not a
practical break.

## 1. What is inherited and what changes

The construction (walk, class key, escapes, merging lemma, phases, meter, report
store, counter widths), the counted kernels (111 + 523), the prices (step 0.752 =
634 counted + 118 analytic rows, doubling 0.956, traversal 0.776), every hard row and
the H1 calibration `kappa = 1.01` are those of submission 412b63f7, specified in
`upstream-proof.md` (412b63f7's file with its section-0 history abridged).
`census.py`, `census-q29.json` and `kernels/` are byte-identical to 412b63f7's
(minus `kernels/test_native.py`, section 6).

The inherited ledger pays the meter's random surcharges with a **Markov
allowance** `E/h`, giving up `h = 0.003830` of success; that forces
`p0 = 0.394331`.

This package uses a **table-conditional allowance** instead (section 3): twice the
doubling mean, eight times the traversal-cycle mean, a retrace row for steps after
inter-trail collisions, and the inherited starts and small-charge rows at `2^30`
times their expectations. The early-halt probability falls below `2.8e-9`, so
`p0 = 0.3900003`. Nothing else changes.

| | 412b63f7 | this package |
|---|---|---|
| H1 floor p0 | 0.394331 | **0.3900003** |
| meter allowance | ceil(E/0.003830) = 2^119.013 | **2^112.572** (section 3) |
| early-halt probability | 0.00383 | **< 2.8e-9** |
| phase-1 operations | S_ops + 2^60 | S_int + 2^61 (covers the < 2^60.41 block overshoot) |
| W_ledger | 2^126.3222050 | **2^126.3028643** |

The previous attempt, submission 0cc8702b (rejected), used an unconditional
premise `Var(D) <= 2^80*mu_D`, which misses `Var(E[D|table])`, and a merge premise
not implied by H1. Section 3.7 shows how both are removed.

## 2. Notation

- `R = 2048` table rows, `M = 2^20` walks, `L_max = 20*2^64` trail cap, `N = (r-1)/6`
  classes, `N' = N*R/(R-1)`.
- `S_ops = kappa*sqrt(2*N'*ln(1/(1-p0)))`, `S_int = ceil(S_ops)`. Phase 1 performs
  at most `n = S_int + 2^61` operations (steps, doublings, traversal steps).
- A **label** of a class arrival is `l = (h, t, d)`: the partition index `h` (bits
  0..10 of the class key `K`), the twist index `t` of the arriving member (which of
  the six units maps it to the canonical member) and the DP field `d` (bits
  64..127 of `K`).
- An arrival is **fresh** if no trail has visited its class before.
- `q_D = 1/(6R) + 2^-19.9954 = 2^-13.5681` and `q_V = 2^-19.9954`: the
  upstream rates of fruitless-cycle entries of length 2..64 and 3..64 per step
  (`upstream-proof.md` section 4, derived from the `S2(L, b)` pattern count).
- Surcharges over the step price: doubling `0.204`, traversal step `0.024`.

## 3. The meter argument

### 3.1 Premises

**H1 (inherited, unchanged).** `P(T_c <= S_ops) >= p0` with `kappa = 1.01`, and the
three further claims listed in `upstream-proof.md` section 4: the partition, DP,
node-routing and bucket fields of `|y|` are pseudo-random functions of the class;
excluding detected fruitless cycles of length at most 64, a trail makes a
non-fruitless return with probability at most `L_max^2/(2N)`; the canonical twist
index of a new class is uniform and independent of its partition index.

**H2 (table-conditional labels).** Fix any outcome `T` of the table coins.
Conditional on `T` and on everything the attack has done before, the label of each
fresh arrival is uniform on `[R] x [6] x [2^64]` and independent of that history.

H2 restates H1's field and twist claims per table and as independence from the
past; the upstream rate `q_D` is computed in exactly this model. The labels are bits
of a point no trail has reached; the table decides *which* new points are reached,
not their low bits. Section 3.8: the ledger needs H2 only to a factor 1.9 in the
mean.

**H3 (collision count, first moment).** Let `K` be the number of inter-trail
collisions in phase 1: arrivals of a trail, by a step from a fresh arrival, a
doubling or a start, at a class that another trail visited earlier. Then `E[K] <= 2^14`.

Under the random-mapping model that also underlies H1, `E[K]` is about
`n^2/(2N') = kappa^2*ln(1/(1-p0)) = 0.504`; H3 allows `2^15` times that. The toy
census of section 4 measures `E[K]` at about `0.9x` the model value.

### 3.2 Which charges are random

Phase-1 operations (`0.752` each) and the hard rows are fixed; the meter halts
early only if the remaining charges exceed the allowance. Those charges are:

- `0.204` per doubling and `0.024` per traversal step in phase 1;
- `2^13` per trail start;
- the small charges (random words, sampling attempts, small batches, fallbacks).

Classify the phase-1 arrivals of a trail. On the complement of the failure events
of section 3.6 (no non-fruitless return, no exceptional addition, no traversal
longer than 64), each arrival is exactly one of:

1. **fresh**;
2. **cycle-internal**: a repeated arrival inside a fruitless cycle of length at most
   64, between entering the cycle and escaping it (at most 3 arrivals for a 2-cycle,
   at most 128 + 64 + 64 for a checkpoint detection and its traversal);
3. **retrace**: an arrival after the trail collided with another trail, until the
   trail ends or leaves the other trail's path; there are at most `L_max + 65` such
   arrivals per collision, so at most `K*(L_max + 65)` in total.

From a fresh arrival the next move is always an ordinary step: the 2-cycle and
checkpoint tests fire only on a repeat, the deferred-escape flag is set only after a
repeat, and a distinguished class or the cap ends the trail. Let `Y_k = 1` if the
step from fresh arrival `k` lands on a class of the trail's current segment (the
classes since its last escape or start). Every doubling is either preceded by such
a cycle entry or happens at a retrace arrival, and each cycle entry causes exactly
one doubling and, for a cycle of length 3..64, at most 64 traversal steps. So

```text
D <= sum_k Y_k + K*(L_max + 65)
V <= 64 * sum_k Y_k^(>=3) + K*(L_max + 65)
```

with the sums over fresh arrivals, and the surcharge of a retrace arrival is at most
`0.204`.

### 3.3 Cycle entries are local functions of the labels

A cycle of length `L` entered from fresh arrival `k` uses the ordinary steps of the
trail's last `L` fresh arrivals `k-L+1, ..., k`, all in the current segment. On the
complement of the non-fruitless-return event it is a fruitless cycle: its index and
twist pattern cancels identically. Define

```text
Y''_k = 1 if, for some L in 2..64, the labels of the trail's last L fresh arrivals
        (ending at k) form an identically cancelling pattern.
```

Then `Y_k <= Y''_k`: `Y''_k` ignores segment boundaries, so it only adds patterns.
`Y''_k` is a function of the labels of at most 64 consecutive fresh arrivals of one
trail. Under H2 these labels are i.i.d. uniform for every `T`, so by the upstream
pattern count

```text
E[Y''_k | T] <= 1/(6R) + sum_{L=3..64} sum_b S2(L, b) R^-(L-b) <= q_D,
E[Y''_k^(>=3) | T] <= q_V.
```

The `L = 2` term is `1/(6R)`: the index must repeat (`1/R`) and the twist must be
the one that maps the sum back (`1/6`).

### 3.4 The trail tapes and a maximal inequality

By H2 the labels can be pre-sampled. Give trail number `t` (in start order) a tape
`l_{t,1}, l_{t,2}, ...` of i.i.d. uniform labels, independent across `t` and of `T`,
and let trail `t` use `l_{t,i}` at its `i`-th fresh arrival. The attack's
distribution given `T` is unchanged, because each fresh label is uniform and
independent of everything before it.

Let `F_t` be the position of the first label on tape `t` whose DP field is zero, and
`Lambda_t = min(F_t, L_max + 65)`. Then `E[Lambda_t] <= E[F_t] = 2^64`, and
`Lambda_t` depends only on DP fields, which H2 makes independent of the `(h, t)`
fields that `Y''` reads. Trail `t` has at most `Lambda_t` fresh arrivals: it ends at
its first fresh DP arrival, or earlier, and has at most `L_max + 65` operations.

**Trail mass.** The trails started in phase 1 satisfy

```text
sum_t Lambda_t <= n + 65*(n/L_max) + M*(L_max + 65) + K*(L_max + 65)
              <= Lambda := n + 2^89 + 2^44*(L_max + 65)      when K <= 2^44.
```

A trail that ended at a fresh DP arrival has `Lambda_t = F_t` fresh arrivals and at
least that many operations. A capped trail has `Lambda_t <= L_max + 65` and used
`L_max` operations, and at most `n/L_max` trails are capped. At most `M` trails are
open when phase 1 ends. A trail that ended at a DP after a collision is charged
`L_max + 65`, and there are at most `K` of them. Exceptional additions and long
traversals are failure events. Exact: `65*ceil(n/L_max) + M*(L_max + 65) < 2^88.33`.

**Centring.** Write `Z_t = sum_{i<=Lambda_t} Y''_{t,i}` and
`X_t = Z_t - E[Z_t | Lambda_t]`. Since `Lambda_t` is independent of the labels read
by `Y''`, `E[Z_t | Lambda_t] <= q_D*Lambda_t`, so

```text
sum_{fresh k} Y_k <= sum_{t<=tau} Z_t <= q_D*Lambda + max_{j<=J} sum_{t<=j} X_t,
```

where `tau` is the number of trails started. `tau <= J = 2^114` on every run,
because each start is charged `2^13` units before it runs and `W/2^13 < 2^114`.

**Variance.** Given `T`, the `X_t` are i.i.d. with mean zero. `Y''_{t,i}` and
`Y''_{t,j}` are independent when `|i - j| >= 64`, so

```text
Var(X_t | T) = E[Var(Z_t | Lambda_t)] <= E[127 * q_D * Lambda_t] <= 127 * q_D * 2^64.
```

**Kolmogorov's maximal inequality** for the martingale `sum_{t<=j} X_t` then gives

```text
P(max_{j<=J} sum_{t<=j} X_t >= q_D*Lambda | T) <= J*127*q_D*2^64/(q_D*Lambda)^2
                                               = 2^-54.87.
```

For traversal cycles the same steps with `q_V` and threshold `7*q_V*Lambda` give
`J*127*q_V*2^64/(7*q_V*Lambda)^2 = 2^-54.06`. The traversal row uses eight times its
mean because its toy evidence is sparser (section 4). These bounds hold for **every**
`T`, so they hold unconditionally. Hence, outside these events and `K > 2^44`,

```text
D <= 2*q_D*Lambda + 2^44*(L_max + 65),     V <= 64*8*q_V*Lambda + 2^44*(L_max + 65),
```

and the doubling and traversal surcharges are at most
`0.204*2*q_D*Lambda + 0.024*64*8*q_V*Lambda + 0.204*2^44*(L_max + 65)`, where the
last term charges every retrace arrival at the doubling surcharge.

### 3.5 Starts, small charges and collisions

- **Starts.** The inherited expectation is `2^13*(S_ops/2^64 + 2^20) = 2^75.71`
  units. Reserving `2^30` times that, Markov gives a halt probability of at most
  `2^-30`.
- **Small charges.** These are below `2^70` units in expectation; reserving `2^30*2^70`
  gives at most `2^-30`.
- **Collisions.** By H3 and Markov, `P(K > 2^44) <= 2^14/2^44 = 2^-30`.

### 3.6 Total early-halt probability and success

The allowance is the sum of five rows:

```text
doublings      0.204*2*q_D*Lambda          = 2^111.8525584
traversals     0.024*64*8*q_V*Lambda       = 2^110.3377964
retrace        0.204*2^44*(L_max + 65)     = 2^110.0285692
starts         2^30*2^13*(S_ops/2^64+2^20) = 2^105.7139511
small charges  2^30*2^70                   = 2^100
allowance                                  = 2^112.5722739
```

The meter halts early only on one of these events:

| event | bound |
|---|---|
| fresh-arrival doublings exceed `2*q_D*Lambda` (Kolmogorov, H2) | `2^-54.87` |
| traversal cycles exceed `8*q_V*Lambda` (Kolmogorov, H2) | `2^-54.06` |
| `K > 2^44` (Markov, H3) | `2^-30` |
| starts exceed `2^30` times expectation | `2^-30` |
| small charges exceed `2^30` times expectation | `2^-30` |
| a non-fruitless return in any phase-1 trail (H1: `(L_max+65)^2/(2N)` per trail, `E[trails] < 2^63`) | `2^-54.77` |
| **total** | **< 2.8e-9** |

Cycles of length 65..200 (rate `< 2^-206` per step) and exceptional additions
(`< 2^-114`) are in the inherited failure table. Adding the inherited events other
than the meter (`< 2e-7`):

```text
success >= p0 - 2e-7 - 2.8e-9 = 0.3900003 - 0.0000002 - 0.0000000028
        = 0.3900000972 > 0.39.
```

### 3.7 The two points of the 0cc8702b rejection

1. **Total variance.** No unconditional variance of `D` is used: the tail of section
   3.4 holds for each fixed table, and `P(D >= a) <= sup_T P(D >= a | T)`, so
   `Var(E[D|T])` never enters. Each table needs `E[Y''_k | T] <= q_D`, which is H2,
   stated per table.
2. **Merges.** The premise `E[Z_i*K_i] <= 2^10*E[Z_i]` is gone. Every arrival after
   a collision is charged the doubling surcharge, for up to `2^44` collisions; H3 is
   a first moment of the collision count with a `2^15` cushion.

### 3.8 How much the ledger tolerates

- **Mean.** If a table's true rates were `q' <= 1.9*q_D` and `q'_V <= 7.9*q_V`, the
  same steps with thresholds `0.1*q_D*Lambda` and `0.1*q_V*Lambda` give tails
  `< 2^-47.3` and `< 2^-38.8`. So H2 is needed only to a factor 1.9 in the doubling
  mean and 7.9 in the traversal-cycle mean.
- **Variance.** The local bounds `127*q*Lambda_t` could be `2^23` times too small
  before either tail reaches `2^-30`.
- **Collisions.** `E[K]` could be `2^15` times the random-mapping value.

## 4. Evidence for H2 and H3 (`h2_census.py`)

`h2_census.py` is a standalone census on the toy subgroups of `experiment.py` (q25:
`r = 33563419`; q29: `r = 178946401`) with the walk rules of `upstream-proof.md`
section 2. It is evidence only at its executed toy scale.

**Part "tables" (H2).** `G` fresh tables, `S` trails per table from random starts,
`m` arrivals each. A one-way random-effects ANOVA of per-trail doubling and
traversal-step counts estimates the between-table relative standard deviation of the
rate, with a 95% upper bound. H2 predicts `F = MS_between/MS_within ~ 1`. In **strict** mode, which is the regime of section 3, a trail is kept
only if it reaches `m` arrivals without a merge into another trail of the same table
and without a non-fruitless return (a class seen before its last escape, or a
detected cycle longer than 8). Each table draws trails until it has `S` kept trails.

| run | tables x trails x arrivals | doublings | rate / 1/(6R) | doublings: F | between-table rel. sd (95% upper) | traversal steps: F (95% upper) | discarded (merge / nfr) |
|---|---|---|---|---|---|---|---|
| q25, R=16, strict | 1024 x 16 x 256 | 42695 | 0.977 | 1.028 | 0.025 (0.049) | 1.026 (0.33) | 1530 / 61 |
| q25, R=64, strict | 1024 x 16 x 256 | 10656 | 0.976 | 0.956 | 0.000 (0.053) | 0.996 (1.14) | 1493 / 36 |
| q29, R=64, strict | 1024 x 16 x 256 | 10862 | 0.994 | 0.951 | 0.000 (0.047) | 1.018 (1.25) | 281 / 7 |
| q29, R=2048, strict | 512 x 8 x 1024 | 358 | 1.049 | 0.999 | 0.000 (0.406) | 1.000 (7.57) | 608 / 8 |
| q25, R=64, raw (control) | 1024 x 16 x 256 | 10816 | 0.990 | 1.155 | 0.122 (0.153) | 1.315 (3.66) | kept |

Rates per arrival sit slightly below `1/(6R)` at small `R` (the denominator counts
cycle-internal arrivals); at `R = 2048`, `8.54e-5` against `q_D = 8.23e-5` (+0.7 s.e.). With the
strict filter, no run shows a between-table component: every `F` is within
`1.1` standard errors of 1, and the 95% upper bounds on the between-table spread of
the doubling rate are `4.7%` to `5.3%` at 1024 tables. The ledger tolerates `90%`
(section 3.8). Traversal steps are sparse (0.01 to 0.2 per trail) and their bounds
are wide; that is why their row uses eight times the mean.

**Raw control.** Without strict mode, q25 at `R = 64` gives a significant `F`. The
outliers trace to two toy-scale mechanisms the proof excludes: merges between a
table's trails (about 1.3 per table at `N = 2^22.4`), and non-fruitless returns,
for example a 2-cycle's doubling `2A` landing one step upstream of `A` (probability
near `6R/N` per escape: `7e-5` at q25, `2^-240` at production). With those trails
removed, the between-table component disappears.

**Part "collisions" (H3).** It runs 64 walks on one fresh table up to the toy
budget `n = kappa*sqrt(2N' ln(1/(1-0.39)))` arrivals (`R = 2048`, DP trails) and
counts inter-trail collisions `K` exactly as H3 defines them.

| run | trials | budget n | model n^2/(2N') | mean K | s.e. | max K |
|---|---|---|---|---|---|---|
| q25, R=2048, 64 walks, dp_bits 5 | 2000 | 2375 | 0.5039 | 0.4670 | 0.0150 | 4 |
| q29, R=2048, 64 walks, dp_bits 6 | 1000 | 5485 | 0.5041 | 0.5250 | 0.0223 | 4 |

Both means are within `2.5` standard errors of the model; the largest `K` is 4.

The hosted run of `experiment.py` reruns both parts at a smaller size through its
`"h2"` input (after the inherited parts, stopping 14 seconds after start). The
offline runs above are `echo '<config>' | python3 h2_census.py` with these inputs
(plus `"seconds":100000`):

```text
strict_q25_R16    {"curve":"q25","R":16,"tables":1024,"per_table":16,"arrivals":256,"strict":true,"coll_trials":2000,"coll_walks":64,"dp_bits":5,"seed":11,"tables_share":0.5}
strict_q25_R64    {"curve":"q25","R":64,"tables":1024,"per_table":16,"arrivals":256,"strict":true,"coll_trials":0,"seed":12,"tables_share":0.9}
strict_q29_R64    {"curve":"q29","R":64,"tables":1024,"per_table":16,"arrivals":256,"strict":true,"coll_trials":1000,"coll_walks":64,"dp_bits":6,"seed":13,"tables_share":0.5}
strict_q29_R2048  {"curve":"q29","R":2048,"tables":512,"per_table":8,"arrivals":1024,"strict":true,"coll_trials":0,"seed":14,"tables_share":0.95}
raw_q25_R64       {"curve":"q25","R":64,"tables":1024,"per_table":16,"arrivals":256,"strict":false,"coll_trials":0,"seed":12,"tables_share":0.95}
```

## 5. Exact budget

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.3900003)))
S_ops = 1.01*S0                          = 2^126.7139511
phase1 = 0.752*(S_int + 2^61)            = 2^126.3027557
hard extras                              = 2^107.1555015
allowance (section 3.6)                  = 2^112.5722739
W_ledger = phase1 + hard + allowance     = 2^126.3028643
success >= 0.3900003 - 2e-7 - 2.8e-9     = 0.3900000972 > 0.39
```

The hard extras are the inherited rows (the store's worst-case overflow row
dominates at `2^107.1536`); block refills are recomputed at the new cap.

`W_ledger = 2^126.3028643`, so the four-decimal ceiling is `126.3029`; `126.3028` is
below `W_ledger` and is not claimed. `certificate.py` proves every row, the tails and
both claim directions in exact integer and rational arithmetic
(`W^10000 < 1000^10000*2^1263029` and `W^10000 > 1000^10000*2^1263028`).

## 6. Files and reproduction

New: `proof.md`, `certificate.py` (`CERTIFICATE: PASS`), `h2_census.py`. `experiment.py` is inherited with its ledger recomputed
(`P0 = 0.3900003`, conditional allowance, tails, claim `126.3029`). Inherited and
byte-identical: `census.py`, `census-q29.json`, `kernels/` (without the native AArch64
harness `test_native.py`, omitted for the archive limit; identical copy in 412b63f7).
`upstream-proof.md` is 412b63f7's with its section-0 history abridged.

```sh
python3 certificate.py
echo '{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,
       "seconds":6,"pipeline":true,"kernels":true,"h2":{"curve":"q25","R":16,"tables":64,
       "per_table":8,"arrivals":256,"coll_trials":200,"coll_walks":64,"dp_bits":5,"seed":1,
       "tables_until":10,"until":14}}' | python3 experiment.py
echo '{"curve":"q25","R":16,"tables":64,"per_table":8,"arrivals":256,"strict":true,
       "coll_trials":200,"coll_walks":64,"dp_bits":5,"seed":1,"seconds":12}' | python3 h2_census.py
```

## 7. Provenance and limitations

The construction, counted loops, prices, report store, hard rows and H1 calibration
are credited, as design provenance, to submission 412b63f7, which is declared as a
dependency, and through it to 89b2ddfd, the public `six-class-rho-secp-counted-v8`
package and its credited ancestors. The idea of replacing the Markov allowance comes
from submission 0cc8702b (rejected), also declared. No prior acceptance, review
outcome, score or cross-submission comparison is a premise of this claim.

H1, H2 and H3 are heuristics supported by random-mapping reasoning and toy
censuses at `N` of about `2^22`–`2^25`, not proved at production. The instruction
counts are analytic counts of the shipped assembly plus analytic rows. The attack is
generic and impractical, not a practical break.

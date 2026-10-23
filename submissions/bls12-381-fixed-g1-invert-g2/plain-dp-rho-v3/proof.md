# bls12-381-fixed-g1-invert-g2: parallel distinguished-point rho in GT

Status: ready for independent review under the exploratory policy. The attack is
a parallel distinguished-point Pollard rho for the discrete logarithm in GT,
with a hard work cap, exact (non-fingerprinted) stored elements, an enforced
storage capacity, and a success claim that is the rule's threshold (0.39), not
the model's value. Nothing here is claimed to be novel.

This package deliberately **does not** use the inversion (two-class) quotient.
An earlier draft of mine on this track did; the review of that draft (submission
`4656a87a-5ac8-432b-...`) correctly pointed out that applying the same multiplier
to both representatives is not a function on inversion classes, which invalidated
the quotient argument and the `sqrt(2)` saving. That draft is withdrawn, and the
present claim is ordinary rho in the full group GT, whose algebra is elementary.

## 1. Claim

Algorithm `A`, run once per challenge instance on fresh coins:

1. Let `Q0` be the public G2 generator and compute
   `g = e(P,Q0) = pairing(Q0,P,final_exponentiate=True)`. The pinned pairing is
   nondegenerate and `P`, `Q0` have prime order `r`, so `g != 1` and `g`
   generates GT.
2. If `z = 1`, return `Q = O` (then `e(P,O) = 1`).
3. Otherwise search `x = log_g z` with parallel distinguished-point rho
   (Section 3) under the hard cap of Section 4; return FAILURE at the cap.
4. Verify `e(P,[x]Q0) = z` and `[r]Q = O`; return `Q = [x]Q0`.

Correctness is immediate: GT is cyclic of prime order `r`, so `log_g z` is unique
modulo `r` and `e(P,[x]Q0) = g^x = z`. The output is verified with a pairing, so
a wrong logarithm is never returned as success.

The complete per-instance work is hard-capped at

```text
W_max = 1.5 * sqrt(r) * 1.01 + 2^31
      = 1.5 * 2^127.428545 * 1.01 + 2^31
      = 2^128.013507 * 1.01 + 2^31
      = 2^128.027863   target-group-operation equivalents
```

with `log2 r = 254.857089` for the BLS12-381 subgroup order `r`. The exact
rational claim is `128027863/1000000`. Success probability is at least `0.39`
over uniform nonidentity `P`, uniform `z` in GT and fresh coins (Section 5).
Peak storage is hard-capped at `2^80` bytes by an enforced capacity, with every
component added up (Section 6).
Preprocessing is charged in full per instance with no reuse credit; there is no
nonuniform advice and no auxiliary powers.

## 2. Reduction and correctness

Let `r` be the prime subgroup order fixed by the target profile and
`e: G1 x G2 -> GT` the pinned pairing in the library's `pairing(Q,P)` argument
order. GT is cyclic of order `r`. For nonidentity `P` in G1 the map
`Q -> e(P,Q)` is an isomorphism `G2 -> GT`, so for every `z` in GT there is a
unique `a` mod `r` with `g^a = z` and `Q = [a]Q0` is the unique preimage of `z`.
`P` is used only through `g`, which is recomputed and charged per instance; no
assumption is made that `P` is a fixed generator across instances. The identity
target is answered directly.

## 3. The search

- **Two independent hash slices.** `H_1(y)` is the low 64 bits of the first `Fp2`
  coordinate of `y` and `H_2(y)` the low 64 bits of the second `Fp2` coordinate.
  The walk index uses `H_1` and the distinguished-point test uses `H_2`, so a
  distinguished point is not forced into a single bucket. (An earlier draft used
  the same bits for both, which would have confined distinguished points to the
  bucket `h = 0`; that is fixed here.)
- **Walk.** `y_{n+1} = y_n * M_{h(y_n)}` with `h(y) = H_1(y) mod R`, `R = 2^20`,
  and
  `M_j = g^{alpha_j} z^{beta_j}` for random `alpha_j, beta_j` mod `r`. The
  exponent pair `(alpha, beta)` with `y = g^alpha z^beta` is carried along with
  two 255-bit modular additions per step. `M = 2^20` walks run in lockstep from
  independent random starts `g^{alpha} z^{beta}`.
- **Distinguished points.** `theta = 2^-60`: a point is distinguished when
  `H_2(y) ≡ 0 (mod 2^60)`.
- **Records are exact.** A record is the full 576-byte `Fp12` element (12 `Fp2`
  coordinates of 381 bits) plus its two 256-bit exponents. There are no
  fingerprints, so a stored-key match is a true group-element match: false
  matches have probability zero rather than being merely unlikely, and no
  false-match resolution procedure is needed. Hash-bucket collisions (different
  elements, same bucket) are resolved by the exact 640-byte comparison and cost
  only the probe; they are never reported as matches.
- **Dictionary, with a hard lookup cap.** The dictionary is open-addressed with
  `2^71` slots of 8-byte pointers (`2^74` bytes) into the record array, at most
  `P_max = 4` probes per lookup (a fourth-probe miss drops the record instead of
  extending the probe sequence). Each lookup therefore costs at most 4 probes of
  a 640-byte exact comparison, budgeted at 2000 instructions. Lookups happen only
  when the distinguished-point test fires, and the number of lookups is itself
  hard-capped at `L_max = 2^70`: after `L_max` lookups the run stops matching
  (existing records are kept) and keeps walking. So the dictionary work is
  `<= L_max * 2000` instructions deterministically, independent of how often the
  distinguished test fires, and the effect of the cap on success is covered by
  the margin between the claimed 0.39 and the model's 0.675.
- **Collisions.** Two records with the same element and different exponent pairs
  give `x = (alpha1 - alpha2)/(beta2 - beta1) mod r`. `beta1 = beta2` is
  impossible for distinct records, because it forces `alpha1 = alpha2`.
- **Cycles.** No categorical claim is made that a trajectory cannot enter a
  cycle containing no distinguished point. If that happens the walk stops making
  progress and the run fails at the cap. No analytic bound on that event is
  claimed: short cycles would make the naive `exp(-theta * cycle_length)`
  estimate invalid, so instead the event is simply part of the failure
  probability and is covered by the margin between the claimed 0.39 and the
  model's 0.675. Merged trajectories stay merged because the step depends only on
  the group element, so a collision becomes visible at the first distinguished
  point after the merge.

## 4. Cost ledger

The cost model defines `cost = H + R/C` with `C = 20000` scalar AArch64
instructions per dense Fp12 multiplication. Every charge below is either the
reference operation, a counted instruction, or an explicitly capped overhead; no
expected value is reported as a bound.

### 4.1 Steps

```text
r          = 2^254.857089413              (BLS12-381 subgroup order, 255 bits)
sqrt(r)    = 2^127.428544707
W_max      = 1.5 * sqrt(r) = 2^128.013507207 steps   (hard cap)
```

`1.5 sqrt(r)` is 1.197 times the random-mapping expectation `1.25 sqrt(r)`, so
the cap already contains the standard rho-constant margin; the model success at
this cap is `1 - exp(-1.5^2/2) = 0.675` (Section 5).

### 4.2 Per-step price

One step is one dense Fp12 multiplication (the reference operation, 1 unit)
plus the following scalar work, itemised:

| Item | Instructions |
| --- | --- |
| load `H_1(y)` (low 64 bits of the first Fp2 coordinate) | 4 |
| bucket index `h = H_1 mod 2^20`, multiplier table pointer chase (2 dependent loads) | 24 |
| exponent update: two 255-bit additions mod `r` (4 limbs, add + conditional subtract each) | 40 |
| distinguished test `H_2 mod 2^60 == 0` (mask + compare) | 8 |
| loop, state spills, branch-misprediction allowance | 40 |
| streamed multiplier-table cache-miss allowance (`R*576 B = 2^29.2 B`) | 40 |
| **total scalar, excluding dictionary work** | **~160, charged as 200** |

Hence `1 + 200/20000 = 1.01` units per step. No inversion, exponentiation or
data-dependent branch occurs inside the loop.

Dictionary work (lookups and inserts) is **not** amortised at the
distinguished-point density (that would be an expected frequency, not a bound).
It is charged as a deterministic worst case under the caps of Section 3: at most
`L_max = 2^70` lookups or inserts, each at most 2000 instructions (four 640-byte
exact comparisons plus pointer updates), i.e.

```text
2^70 * 2000 = 2^81 instructions = 2^81 / 20000 = 2^66.678072 units,
```

which is `2^-61` of the search term and is included in the overhead of
Section 4.3.

### 4.3 Setup, table, verification (a hard-capped overhead of `2^31` units)

| Item | Bound |
| --- | --- |
| `g = e(P,Q0)` (one pairing) | 2^15 |
| `R = 2^20` multipliers `M_j = g^{alpha_j} z^{beta_j}` (two GT exponentiations each, ~380 Fp12 multiplications apiece) | 2^30 |
| `M = 2^20` walk starts `g^{alpha} z^{beta}` (two GT exponentiations each) | 2^30 |
| verification: one pairing, one G2 scalar multiplication, output | 2^17 |
| dictionary work (Section 4.2, `L_max * 2000` instructions) | 2^66.678072 |
| **overhead, summed exactly** | **2^31 + 2^15 + 2^17 + 2^66.678072 < 2^67.001** |

The sum is charged exactly, not rounded to its dominant term: `2^30 + 2^30 =
2^31`, and the remaining items add `2^15 + 2^17 + 2^66.678072`, so the whole
overhead is `2^31 + 2^15 + 2^17 + 2^66.678072 < 2^67.001`. It is `2^-61` of the
search term, hence

```text
W = 1.01 * 2^128.013507207 + (< 2^67.001)
  = 2^128.027861974 + (< 2^67.001) = 2^128.027861974
```

to nine decimal places.

### 4.4 Final number

```text
W = 1.01 * 2^128.013507207 + 2^31
  = 2^(128.013507207 + 0.014354767) + 2^31
  = 2^128.027861974 + 2^31
  = 2^128.027861974
score = log2(W) = 128.027863   (rounded upward to six decimals)
```

## 5. Success probability

The walk is modelled as a random mapping on the `r` group elements; this is the
standard rho heuristic (Brent-Pollard correction for `R`-adding walks), and it is
disclosed as a heuristic. For total work `W` steps,

```text
P[collision by W] ~ 1 - exp(-W^2 / (2r)),
```

so at the cap `W_max = 1.5 sqrt(r)` the model success is `1 - exp(-1.125) =
0.675`. The detection delay after a merge (the first distinguished point,
`M * theta^-1 = 2^20 * 2^60 = 2^80` lockstep steps) is `2^-48` of the cap and is
inside it.

The submitted success probability is the rule's threshold **0.39**, well below
the model's 0.675. `experiment.py` Part B, a bounded toy implementation of the
whole algorithm at `2^29`-element subgroup order, measures `0.775` at
`1.5 sqrt(r)` (40 trials) against the model's 0.675, and `0.975` at
`2.5 sqrt(r)` against 0.956. Its uncensored first-collision statistic is
`mean 1.166 sqrt(r)`, `median 1.150`, with `70.8%` of walks colliding before
`1.5 sqrt(r)` (model: mean 1.2533, 67.5%), which is the part of the heuristic the
claim depends on. Toy figures are evidence at toy scope only.

## 6. Memory, preprocessing and advice

- **Capacity, enforced.** Distinguished-point records are inserted into a table
  with hard capacity `C_max = 2^70` records. The design count is
  `theta * W_max = 2^-60 * 2^128.013507 = 2^68.013507` records, 3.9x below the
  capacity; if the count ever reached `C_max` the run stops inserting new records
  (existing records are still matched) and continues walking, so the record count
  is enforced by the algorithm and not merely expected of it. The effect of that
  policy on success is covered by the margin between the claimed 0.39 and the
  model's 0.675.
- **Peak storage, added up.** `C_max * 640 B = 2^70 * 640 B = 2^79.321928 B` for
  the record array (576-byte exact element plus two 32-byte exponents);
  `2^71 * 8 B = 2^74 B` for the pointer/slot index; `2^29.2 B` for the multiplier
  table; `< 2^40 B` for walk state, code and scratch. Total
  `< 2^79.322 + 2^74 + 2^29.2 + 2^40 <= 2^80 B`. The submitted bound is `2^80`
  bytes: the sum of the parts plus an explicit allowance, not a rounded-down
  record-only figure. Memory is reported as the rules require and does not enter
  this track's score.
- **Preprocessing.** The multiplier table and the walk starts are per-instance
  and are charged above (`2^31` units, `2^-97` of the total); no reuse credit is
  taken, and the reported preprocessing bound is the same `128.027863`.
- **Advice.** 0. No auxiliary powers; the multipliers are derived from `g` and
  `z` at runtime.

## 7. Edge cases

- `P` nonidentity, `z` in GT: `g != 1` and `z` has a unique logarithm mod `r`.
- `z = 1`: direct answer `Q = O`, verified with one pairing.
- A walk that hits `y = 1` yields `alpha + x beta = 0`, i.e. `x = -alpha/beta`
  when `beta != 0`; this is a second way for the search to succeed (charged, as
  the verification is).
- `r` is prime, so `[x]Q0` has the required order and lies in G2.

## 8. Heuristics and limitations

1. The random-mapping model for the walk is a heuristic (standard in rho
   analyses); it is supported at toy scale by the experiment and not proved.
2. The `1.01` per-step price is an analytic instruction ledger for this
   implementation on an AArch64 core, not a measurement of a production run.
3. The `2^31` overhead is a charge, `2^-97` of the total.
4. This is an analytic complexity claim for a `2^128`-step computation, not an
   executed attack; the rules allow that.
5. The claim is an upper bound for this algorithm, not a security estimate for
   BLS12-381, and it deliberately does not use any group quotient: the earlier
   two-class draft was withdrawn after review.
6. The model's `1 - exp(-W^2/(2r))` is used only to choose the cap and to justify
   that 0.39 is reachable; the submitted success is 0.39.

## 9. Evidence (`experiment.py`, standard library only)

Part A, on BLS12-381 itself (no third-party imports, so it runs in the verifier's
bounded container): `r` is prime and 255 bits; for random GT elements
`y = h^((p^12-1)/r)` (computed in an explicitly implemented
`Fp2/Fp6/Fp12` tower) `y != 1` and `y^r = 1`. This is the group-order fact the
whole claim rests on.

Part B, a bounded toy implementation of the complete algorithm (plain walk,
`theta`-distinguished points with exact element keys, exponent recovery and
verification, hard step cap): success rates `0.775` at `1.5 sqrt(r)` and `0.975`
at `2.5 sqrt(r)`, plus the uncensored first-collision statistic
`mean 1.166 sqrt(r)` / `median 1.150` / `70.8% below 1.5 sqrt(r)` against the
model's `1.2533` / `67.5%`. Toy-scale statistics are evidence at toy scope only,
never a production complexity proof.

## 10. References

- J. M. Pollard, "Monte Carlo methods for index computation (mod p)",
  Math. Comp. 32, 1978.
- P. C. van Oorschot and M. J. Wiener, "Parallel collision search with
  cryptanalytic applications", J. Cryptology 12, 1999 (parallel
  distinguished-point rho, the detection-delay term used here).
- E. Teske, "On random walks for Pollard's rho method", Math. Comp. 70, 2001
  (`R`-adding walks and their variance).
- R. P. Brent, "An improved Monte Carlo factorization algorithm", BIT 20, 1980
  (the correction factor for the rho constant).
- Curvesma.sh `tracks/bls12-381-fixed-g1-invert-g2/rules.md`,
  `tracks/docs/COST_MODEL.md`, `tracks/classical-group-work-v2.json`.

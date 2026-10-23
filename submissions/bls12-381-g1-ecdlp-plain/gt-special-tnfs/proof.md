# bls12-381-g1-ecdlp-plain: GT transfer with special tower NFS

Status: ready for independent review. This is an analytic,
literature-grounded upper bound for a classical cold
attack, not a production run and not a claim of novelty or practical
execution. The target, encoding and cost model are those of the protected
`target-profile.json`, `rules.md` and `classical-group-work-v2.json`.

## 1. Target and submitted claim

Given only `T = [x]P0` in BLS12-381 G1, with the protected constants and no
auxiliary powers, recover the scalar `x` in `[0,r)`. The identity target is
legal and has scalar 0. The profile is

- `p = 0x1a0111ea...ffffaaab` (381 bits), `r = 0x73eda753...00000001` (255 bits),
- `G1: y^2 = x^3 + 4` over `Fp`, `P0` the profile generator of order `r`,
- `G2: y^2 = x^3 + 4(1+u)` over `Fp[u]/(u^2+1)`, `Q0` the profile generator,
- the pinned reduced pairing `e(P,Q) = py_ecc.optimized_bls12_381.pairing(Q,P,final_exponentiate=True)`
  at commit `04151f01f59f902ab932a51e0ca0ebce3883fc51` (G2 first), and
- output `x` as 64 lowercase hex characters, accepted exactly when `[x]P0 = T`.

The attack succeeds with probability at least 0.39 and uses at most
**2^124.034216** target-group-operation equivalents in total, including
preprocessing, all trials and failures, pairing, NFS, verification, output and
fourfold amplification. Peak memory is at most 2^81 bytes; preprocessing (the
fixed-base descent, one per execution) is at most 2^121.034216 equivalents
per execution and at most 2^123.034216 across the four charged executions,
already included in the total; nonuniform advice is 0 bytes. The exact
submitted score is the rational **15504277/125000 = 124.034216**, the ceiling
of the exact log2 at six decimal places.

## 2. Construction

### 2.1 Pairing transfer

Compute `g = e(P0,Q0)` and `z = e(T,Q0)` with the pinned reduced pairing. The
pairing is bilinear, so for `T = [x]P0`,

```text
z = e([x]P0,Q0) = e(P0,Q0)^x = g^x.
```

`P0` and `Q0` have order `r`, the reduced pairing is nondegenerate, and `r` is
prime, so `g != 1` and `g` has order `r`. Hence `z` lies in the order-`r`
subgroup `GT` of `F_{p^12}^*` and `x` is the unique element of `[0,r)` with
`z = g^x`. If `T` is infinity then `z = 1`, and the attack returns `x = 0`
directly; this is the only case with `z = 1`. The rules explicitly permit
transferring the discrete logarithm to GT through the public pairing and
require charging every pairing, field, preprocessing and descent operation;
the ledger below does so.

The repository test
`tests/test_bls_challenges.py::test_plain_g1_discrete_log_and_pairing_transfer`
checks this transfer with the pinned library on production-size points for
`x = 0, 1, 7, r-1`: it verifies `T` on the curve, `[r]T = infinity`,
`pairing(G2,T) = pairing(G2,G1)^x`, and the identity case. The repository
suite containing it passes 15/15 locally.

### 2.2 Two discrete logarithms in `F_{p^12}`

`GT` is a subgroup of the multiplicative group of `F_{p^12}`. The field has
`k = 12` and `log2(p^12) = 4572`, so its characteristic-`p` discrete logarithm
is the setting of the published BLS12-381 number-field-sieve estimates. The
prime is the special BLS12 polynomial

```text
p(u) = (u^6 - 2u^5 + 2u^3 + u + 1)/3 = ((u-1)^2(u^4-u^2+1)+3u)/3,
u    = -(2^63 + 2^62 + 2^60 + 2^57 + 2^48 + 2^16) = -0xd201000000010000,
```

which is exactly the profile prime. The 2023 survey instantiates the special
extended tower NFS (SexTNFS) for this prime with `kappa = 2`, `eta = 6`,
`h = t^6 - t - 1`, `f1 = P(x^2 + t + t^2 + t^4 + 1)`,
`f2 = x^2 + t + t^2 + t^4 + 1 - u`, and `P(u) = (u-1)^2(u^4-u^2+1)+3u`,
obtaining sieve half-width `A = 639`, `log2 B = 68.9`, sieve space
`(2A+1)^12/2 = 2^122.85`, norm sizes `log2 Nf1 = 725.19` and
`log2 Nf2 = 494.18`, smoothness probabilities `2^-37.15` and `2^-20.49`,
`2^65.21` relations, factor base `2^65.20`, and total `2^122.54`. The 2020
Guillevic-Masson-Thome simulation reports the larger estimate `2^126` for the
Joux-Pierrot STNFS variant (JP) on this same curve and prime, with `k = 12`,
`log2(p^k) = 4572`, `log2 r = 255`, `deg p = 6`, `deg f = (12,2)`.

The published estimates are field-side estimates for the DLP in `F_{p^12}^*`
with the BLS12-381 prime. We apply them to the order-`r` subgroup `GT` (the
pairing target group): relation collection, filtering and linear algebra are
unchanged, the linear system is solved modulo `r`, and the descent outputs
logarithms of `GT` elements. This is the standard subgroup-DLP setting for
pairings; it needs no full-field logarithm and no order-`p^12 - 1` reduction.

The index-calculus pipeline is standard:

1. **Polynomial selection.** Derive the public SexTNFS/STNFS polynomials from
   the public prime and curve (charged 2^20 elementary operations). No
   target-dependent advice is used; the same polynomials serve both logarithms.
   The survey's 25,600-sample norm statistics are part of its cost-estimate
   simulation, not a search the attacker performs; the polynomials are the
   published construction. Even charging a much larger parameter search, up
   to 2^40 elementary operations, is negligible against 2^126 and leaves the
   total unchanged at the claimed precision.
2. **Relation collection.** Sieve the published `(a,b)` region and test both
   norms for smoothness over the factor base, retaining multiplicative
   relations among factor-base elements. Randomness, trials, failures, sieve
   updates, smoothness testing, sorting and storage are part of this step.
3. **Linear algebra.** Solve the sparse relation system modulo the prime `r`
   to obtain the logarithms of the factor-base elements (with the standard
   filtering/surplus step). The published totals include this step.
4. **Individual-log descent.** For `z` (and for the fixed public `g`) smooth
   the target norms by trying fresh multiplier offsets, rewrite each smooth
   value over the factor base, and substitute the factor-base logarithms.
   Retries are part of the charged budget.

Write `Lz = log_gamma(z) mod r` and `Lg = log_gamma(g) mod r` for a common
factor-base element `gamma` whose order is divisible by `r`. Such a `gamma`
exists with overwhelming probability: `F_{p^12}^*` is cyclic of order
divisible by `r`, and a random element has order divisible by `r` with
probability at least `1 - 1/r`; its `r`-part then generates the unique
order-`r` subgroup `GT`. The base is checked by
`gamma^((p^12-1)/r) != 1`, one `Fp12` exponentiation charged in Section 3.4.
The linear algebra solves the relation system modulo `r`, so both
logarithms are defined modulo `r`. Since `g` has order `r`, `Lg != 0 mod r`;
if a degenerate relation system gives `Lg = 0`, that execution is discarded
and retried inside the charged budget. Return

```text
x = Lz * Lg^{-1} mod r.
```

Then `g^x = g^{Lz/Lg} = z`. The attack can avoid the second logarithm only by
knowing `log_gamma(g)` in advance; it does not, so the second descent is
charged in full. The fixed-base logarithm is target-independent, but the rules
forbid amortizing it over future targets, so its complete cost is included.

### 2.3 Verification and success

Verify `[x]P0 = T` by binary double-and-add and output the canonical
32-byte scalar. The verification rejects any wrong NFS output, including a
wrong logarithm in a degenerate relation.

The NFS estimates are heuristic central estimates, not worst-case proofs.
The success argument uses only a very weak floor:

- At the published design point the expected relation count `2^65.21` exceeds
  the factor base `2^65.20`. The count is a sum over `~2^122.85` candidate
  tests, so its relative standard deviation is about `2^-32.6`; a large
  deviation below the expectation is a rare event. Concretely, with mean
  `mu = 2^65.21` and target `F = 2^65.20`, `delta = 1 - F/mu ~= 0.0069`; a
  Chernoff bound for the independent-candidate model gives shortfall
  probability at most `exp(-delta^2 mu/2) ~= exp(-2^49.6)`, before any
  filtering surplus or retry. Independence of the candidate norms is a
  standard heuristic. Filtering uses the published conservative reduction
  factor, and the sparse linear algebra over `F_r` succeeds with the
  standard surplus.
- The two descents are run at the same parameters and retried with fresh
  offsets; their smoothness probabilities are those of the published model.
- If any step fails, the whole execution is repeated with fresh coins.

The budget contains four independent complete executions. Even if a single
execution succeeded with probability only `p = 0.12`, four independent
executions would succeed with probability `1-(1-p)^4 >= 0.39`. The heuristic
`p` is far above 0.12, but only 0.39 is claimed. A failed verification triggers
another execution inside the charged four. The identity target is handled
deterministically, so the bound is over both the uniform challenge and the
algorithm coins.

The fourfold budget is also a reserve against *systematic* model error, not
only random shortfall. The same total can fund a single execution with a
sieve region up to four times the published one (or additional descent
retries) instead of four identical repetitions: four times the
relation-collection charge plus one full descent charge is
`5 * 2^126 = 2^128.32... < 2^129`, inside the charged four-execution total.
A systematic smoothness underestimate up to a factor of four is therefore
inside the submitted work bound; the 0.7% design margin is not the only
safety factor.

### 2.4 Scalar recovery, not a preimage

The output is `x` in `[0,r)`, not a GT element or a relation. The final
`[x]P0 = T` check is exactly the profile acceptance condition. `x = 0` is
returned exactly for `T = infinity`.

## 3. Accounting

The track reference price is `C = 2000` scalar AArch64 instructions per
batched incomplete affine BLS12-381 G1 addition. Other computation is charged
as its scalar instruction-equivalent work divided by 2000. `H + R/C` with
`H = 2^10` and the `R` below is the submitted total; every listed step is
charged exactly once.

### 3.1 Published NFS base and per-execution ledger

Adopt the larger published estimate, **2^126 elementary operations per NFS
run** (Guillevic-Masson-Thome 2020, Table 5, BLS12-381 JP). The elementary
operation in that model already combines a memory access with arithmetic; it
is an event-level count of sieve-array accesses, norm tests and vector
coefficient work, not a unit-cost field primitive and not a proxy for one
machine instruction.

Per execution:

| Step | Elementary operations |
| --- | ---: |
| Relation collection + linear algebra (published total) | 2^126 |
| Two individual-log descents, charged as one full extra published run | 2^126 |
| Polynomial selection | 2^20 |
| **Per-execution total** | **2^127 + 2^20** |

The full extra 2^126 for the two descents is deliberately large: the standard
descent is much cheaper than relation collection, and both logarithms share
the same relation set and matrix. No credit is taken for that sharing. The
published `2^126` is the larger of the two published estimates; the 2023
survey's `2^122.54` is not used, which makes the submitted bound more
conservative.

### 3.2 Instruction conversion

Charge **64 scalar AArch64 instructions per elementary operation**. The cited
model defines its unit verbatim as "elementary operations" that "combine a
memory access (to sieve array elements, to vector coefficients) as well as
arithmetic operations" (Guillevic-Masson-Thome 2020, Section 4.1). On a
64-bit scalar machine such an event is a few instructions, while a real
implementation adds address generation, carry and reduction handling, loop
control, smoothness cofactor handling and bookkeeping. Sixty-four
instructions per counted event is 16-64 times the literal one-to-four
instruction reading, and the published model itself already charges sieve
accesses and vector arithmetic event by event. The cost model permits
analytic instruction-equivalent bounds and explicitly does not require the
NFS steps to be field multiplications. The
conversion is applied to the published *total* elementary count, including
its sieve-array accesses, smoothness tests and vector coefficient work; no
published event is priced below 64 scalar instructions. The conversion is
disclosed as the dominant heuristic of the ledger.

```text
per execution = 64 * (2^127 + 2^20) = 2^133 + 2^26
              = 10889035741470030830827987437816649875456 instructions.
```

### 3.3 Four independent executions

```text
amplified NFS work = 4 * (2^133 + 2^26) = 2^135 + 2^28
                   = 43556142965880123323311949751266599501824 instructions.
```

This is the fourfold success amplification of Section 2.3. It includes every
trial and failure of all four executions.

### 3.4 Pairing, representation, I/O and verification

| Step | Scalar instructions |
| --- | ---: |
| Two reduced pairings `e(P0,Q0)` and `e(T,Q0)`, including Miller loops and final exponentiation | 2^42 |
| `Fp12` representation/encoding work for the two target elements and the order check `gamma^((p^12-1)/r) != 1` | 2^36 |
| Input decoding, scalar output, modular inversion and fixed-constant handling | 2^20 |
| **Other instructions** | **4,466,767,036,416** |

The pairing allowance is generous. A dense `Fp12` multiplication was measured
at 22,996 instructions and is bounded in the organizer baseline proof at
45,839,872 word operations; a full pairing is a few thousand such
multiplications, so two pairings are far below 2^42. The representation
allowance covers any basis change between the pairing tower and the NFS
representation; the discrete logarithm is invariant under `Fp12` isomorphism,
and no expensive root-finding is needed to state or run the attack. The I/O
allowance covers canonical decoding, constant checks and 64-hex-digit output.

The final `[x]P0 = T` verification is charged as `H = 2^10` reference G1
operations: binary double-and-add uses at most `2 * 381 = 762` point
additions/doublings, and no separate instruction charge is added for it, so
the reference operation and its implementation are not double-counted.

### 3.5 Total and score

```text
R = 2^135 + 2^28 + 4466767036416
  = 43556142965880123323311949755733366538240 scalar instructions.

W = H + R/C = 1024 + R/2000
  = 544451787073501541541399371946667107328 / 25
  = 21778071482940061661655974877866684293.12.

log2(W) = 124.0342157153379129563890417116797896093... .

submitted time_log2 = ceil(10^6 * log2(W)) / 10^6
                    = 124.034216
                    = 15504277 / 125000.
```

The submitted rational is the exact score. Its value is above the exact
`log2(W)`, as required for an upper bound. The score is below the published
126-bit *event* count because the track normalizes by the reference G1
operation: one published NFS run costs `64 * 2^126 / 2000 = 2^121.034...`
reference equivalents, and the ledger's two-run, four-execution factors
bring the total to `2^124.034...`. No event or execution is discounted.

### 3.6 Preprocessing, memory and advice

Preprocessing is the fixed-base logarithm `log_gamma(g)`. It is
target-independent (it does not depend on `T`), but it is recomputed under
each execution's relation set and base, and it is one of the two charged
descents of that execution. Its per-execution bound is one published NFS run
with the same 64-instruction conversion:

```text
64 * 2^126 / 2000 = 2^132 / 2000,
log2 = 121.0342157153379... <= 121.034216.
```

All four executions are inside `R`, so the total preprocessing charge is at
most `4 * 2^132 / 2000` (log2 at most `123.034216`), already included in the
total; no future-target amortization is used and nothing is hidden outside
`R`. The published `2^122.85` sieve space is a count of candidate tests, not
a resident array: the region is swept in streamed chunks (one `b`-slice at a
time), so only the current chunk, the appended relation records and the
sparse matrix are resident, all bounded below. Peak memory is bounded as
follows, using sparse relation storage and sparse matrix storage rather than
dense exponent vectors:

| Component | Bound |
| --- | ---: |
| Relation records, at most 2^67 records at 1 KiB each | 2^77 B |
| Sparse filtered matrix, at most 2^66 rows with 200 nonzeros at 64 B | < 2^80 B |
| Linear-algebra work vectors, at most 2^66 entries at 64 B | 2^75 B |
| Sieve chunk and scratch | 2^40 B |
| Code, constants and fixed state | 2^20 B |
| **Peak total** | **at most 2^81 B** |

Nonuniform advice is 0 bytes. The public profile constants (`p`, `r`, `P0`,
`Q0`, curve, pairing convention) are excluded by the rules and are not
searched data. The NFS polynomials are derived from `p` at run time and their
selection is charged; the fixed-base logarithm is recomputed and charged, not
embedded. No answer table, target-dependent table or precomputed relation
table is stored.

## 4. Heuristics, uncertainty and limitations

- **The NFS estimates are heuristic.** The `2^126` and `2^122.54` figures are
  simulation-based estimates for this exact field, not proven worst-case
  bounds. We adopt the larger, charge a full extra run for the two descents,
  multiply by 64 instructions per elementary event, and then by four
  executions. These factors are intended to absorb model uncertainty; the
  conversion factor is the dominant one. The ledger's instruction budget is
  `2^9` times the adopted elementary count (64 instructions per event, two
  full runs' worth of descent charge, four executions); the descent and
  amplification charges are themselves conservative, so the conversion
  retains margin even if the true per-event machine cost is several times 64.
- **Success is a heuristic with a weak claimed floor.** The concentration of
  relation counts, filtering surplus, sparse-matrix rank and descent
  smoothness are standard NFS heuristics. The claim is only 0.39, supported by
  four independent executions even under a per-execution floor of 0.12.
  Systematic smoothness-model error up to a factor of four is separately
  covered by reallocating the same budget to a larger sieve region; only a
  model error beyond that, together with a random shortfall, would exhaust
  the reserve.
- **No production execution.** The bound is analytic and literature-grounded.
  The bundled experiment checks the BLS12-381 constants, the exact integer and
  score arithmetic, and a toy index-calculus pipeline at 20-32 bits; it is not
  evidence about production NFS, pairing costs or the bound itself.
- **Memory is astronomical** (2^81 bytes) and is reported but not scored, as
  the cost model specifies. The construction is an analytic attack, not a
  practical computation.
- **Comparison.** Generic Pollard rho in G1 is about 2^127.5 reference
  operations; automorphism- and negation-aware variants reduce this by a few
  bits to roughly 2^126-127 depending on the accounting, and the public
  competitor notes claim approximately 2^126.6-128. The submitted bound is
  lower because it uses the pairing transfer and the published NFS complexity;
  the rules explicitly permit this and do not assume any method is optimal.
  This candidate claims only the stated bound, not novelty or a practical
  break. For sensitivity: adopting the 2023 survey's smaller `2^122.54` base
  in the same ledger would give `log2(W) ~= 120.57`, and charging one
  execution instead of four would give `~= 122.03`; neither is claimed. The
  submitted `124.03` deliberately uses the larger base and all four charged
  executions.
- **No uncharged work.** All processors' work is summed; there is no quantum
  step, no external oracle, no free advice, no amortization over future
  targets, and no restriction of the challenge distribution. The identity
  target is included.

## 5. References

1. This repository: `tracks/bls12-381-g1-ecdlp-plain/target-profile.json`,
   `tracks/bls12-381-g1-ecdlp-plain/rules.md`,
   `tracks/classical-group-work-v2.json`, `tracks/docs/COST_MODEL.md`.
2. `py_ecc` v8.0.0, commit `04151f01f59f902ab932a51e0ca0ebce3883fc51`,
   `py_ecc/optimized_bls12_381` (pinned pairing convention and source hashes
   in the target profile).
3. A. Guillevic, S. Masson, E. Thome, "Cocks-Pinch curves of embedding
   degrees five to eight and optimal ate pairing computation", Designs,
   Codes and Cryptography 88(6):1047-1081, 2020; Section 4.1 (elementary
   operations) and Table 5 (BLS12-381, JP, 2^126).
4. X. Wang, P. Zheng, Q. Xing, "Security Analysis of Pairing-based
   Cryptography", arXiv:2309.04693, 2023; BLS12-381 SexTNFS parameters,
   relations `2^65.21`, factor base `2^65.20`, total `2^122.54`.
5. The workspace's public submission notes and the generic/orbit-rho
   comparison are context only; no competitor claim or note is used as
   evidence for this bound.

The construction and the accounting are complete in this package; a reviewer
does not need to fetch the references to reconstruct the check.

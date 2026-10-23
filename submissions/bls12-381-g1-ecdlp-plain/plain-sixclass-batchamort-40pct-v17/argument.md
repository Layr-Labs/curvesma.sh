# plain G1 v17 — batch-amortised six-class walk

Self-contained: `dependencies.json` is `[]`, `submission.json.experiments` is `[]`. Every
number below is rendered from `ledger.py` by `scripts/gen_tables.py`;
`scripts/verify_certificate.py` re-renders the block, byte-compares it, re-adds the printed
rows, and re-evaluates every `F:` line in exact rational arithmetic. Nothing is hand-typed.

Track `bls12-381-g1-ecdlp-plain`: input `T = [x]P0`, `x` uniform in `[0,r)`, no auxiliary
powers, no `n`; output `x`, accepted iff `[x]G == T`. `beta` a primitive cube root of unity
mod `p`, `lam` the matching root mod `r`, `C = 2000` instructions per reference G1 addition.

## 1. What v17 changes against v16

| v16 defect | v17 |
|---|---|
| every step paid 383 instructions to canonicalise (`beta*x` + six-candidate argmin) and 138 to twist the pair | no per-step canonicalisation or twist: `|y|` picks the column and `e` is resolved once per record hit at report time (2 Fp multiplies = 800 instructions/hit, in `report_path`) |
| the toy derived the sign of `e` from `e >= Q//2`, but `lam = 26745 > Q//2 = 25023`, so `lam` and `-lam` got the `lam^2` map: e-resolution broke on 2 of 6 elements | `orbit_apply` reads `(k,sign)` off `ORBIT_E`, the table `resolve_e` indexes, so they are exact inverses; the toy resolves both directions for all six `e`, 0 failures required |
| failures were counted after truncating the sample to 4 (`failed: 4` for 50) | counted before sampling: `failed` is the true count, `failures` a 4-entry sample |
| the invariant pass stopped at the first collision, so `inv` could not exceed ~100 | runs the whole budget with `solve=False`; requires `invariant_checks == invariant_pass_steps` |
| batch inversion described, charged only as a 2-instruction row | charged as a total, `batch_inversion_total` (section 3) |

## 2. State, transition, selector

`INV: S = [a]G + [c]T`, `a,c in Z_r`, `S` affine in the prime-order subgroup. With
`K = |y(S)| = min(y,P-y)`, `m = cell(K)` and row `(J_m, r_m, s_m)`:

```
S' = S + J_m,     a' = a + r_m,     c' = c + s_m   (mod r).
```

`|y|` is invariant under all six maps `{+-1, +-phi, +-phi^2}` with `phi(x,y)=(beta*x,y)`,
so the transition uses no coins and is a function of the state alone — which makes the merge
lemma an induction. A nontrivial sixth root of unity would give `[mu-1]P = O` with `r`
prime, so non-identity points split into `N = (r-1)/6` orbits of size six. `S = O` draws a
fresh `(a,c)` from `selector.RestartStream` (SHA-256 counter chain, seed and draw index
logged in `evidence/toy.json`) and re-forms `S`.

`cell(K) = (K*R) >> (bitlen(P)-1)`, `R = 2048`: a **power-of-two quotient**, dividing `K*R`
by `2^(bitlen P-1)` and not by `P`, so it is a 4-instruction shift rather than a 6-limb
reciprocal multiply (~134 instructions). It is deliberately not `floor(K*R/P)` — for the
fixture, `K=2045, R=4096` gives 15 by the shift and 13 by `floor(K*R/P)`, and the
certificate checks that difference explicitly instead of claiming they agree. `K <= P//2 <
2^(bitlen P-1)` so `cell < R` always: no mask and no overflow. The reachable K-range is a
prefix of `[0,2^(bitlen P-1))`, which the shift cuts into R equal cells, so every reachable
column has the same preimage count except the last partial one; that count is printed as
`R_eff` and used in `N' = N*R_eff/(R_eff-1)`. `scripts/selector.py::cell` is the only
implementation, priced by the ledger and imported by the toy.

## 3. Batch amortisation, derived

A batch is `B = 2^18` walks ready at the same meter tick. An affine addition needs one field
inversion, the only per-step cost not already inside the reference unit:

* Fermat inversion in `Fp`: 254 squarings + 11 multiplies = **265 Fp multiplies**; at
  `M = 383` instructions each, `265*383 = 101,495` instructions per inversion.
* Shared over a batch: `101,495/2^18 = 0.387188` instructions per step = `0.000194`
  reference units — a relative charge of `1.8e-4` on the `1.0665` step price, carried as
  its own ledger term.
* Montgomery's trick (one inversion + `3(B-1)` multiplies for `B` elements, i.e. 3
  multiplies per element = `1,149` instructions/step) is **not used and not charged**: those
  three multiplies are the addition's own inversion handling, and the reference unit is one
  G1 addition at 2000 instructions, which contains whatever inversion handling the addition
  performs. Charging them again would charge the addition twice.

Charged as a total, not a rate: a run takes at most `S_ops + 2^60` steps, hence
`n_batches = ceil((S_ops+2^60)/B)` batches, and the ledger prices `n_batches + 1`. The extra
whole batch pays the tail batch closing the last partial group, so a run whose final batch
holds one walk is still fully charged; the amortisation never enters at average fill.

## 4. Capped-work discipline (the score is a cap)

Every row is a worst case: the step is `1` unit + `STEP_INSN` instructions on all
`S_ops + 2^60` permitted steps; a doubling escape (at most one per 64-step checkpoint
interval) costs `DOUBLING_PRICE`, below the step price, so escapes are absorbed inside
`S_ops`; the inversion prices `n_batches + 1` whole batches; records use `mean +
4*sqrt(mean) + workers`, a deviation cap whose `<= 2^-10` exceedance is `pool_overflow`;
restarts allow 16 per worker while `S = O` has probability `~1/Q` per step over `2^126`
steps; the report path prices one e-resolution per record cap, not per realised hit. Only
the two amortisation ratios are averaged, and both enter at worst-case fill.

## 5. Distinguished points, geometric transport, merge lemma

Spacing `d = 60` on a 384-bit point fingerprint; a record is
`key(48) || a(32) || c(32) || pad(80) = 192 B`, so recovery needs no re-walk.

Transport is geometric, not half-spacing: DP indicators are independent `Bernoulli(2^-60)`
along a walk, so the wait from a merge to the next DP is geometric with mean `2^60`, and
memorylessness makes the *residual* wait the same mean. Merges lost are at most
`2*2^60/S_ops` (factor two: merge density grows linearly in `t`).

Merge lemma, full state. *Hypotheses:* (H1) two walks hold the same `S` at the same step
count; (H2) `cell` and the transition use no coins, so `m` and `S+J_m` are functions of `S`;
(H3) the coefficient update is affine with the same `m` for both. *Claim:* for all `t >= 0`,
`S_t^(1) = S_t^(2)` and `(a,c)_t^(2) = (a,c)_t^(1) + (da,dc)`, so `da*dc^-1` is invariant and
equals `x`. *Proof:* induction — equal points give the same column, hence the same summand
and the same point difference, and a doubling escape scales both coordinate differences by
2, preserving `da/dc`. Detection inspects the point only, so both walks report at the same
step with the same `x`. ∎

A record hit with class-mate points `S1 = [e]S2` reports `x = (e*a2 - a1)/(c1 - e*c2)`, `e`
resolved from the two stored points — the e-resolution the toy exercises in both directions
for all six `e`. v16's crash was a pair twisted on one side only; the fixed `orbit_apply` is
now `resolve_e`'s exact inverse by construction.

## 6. Toy and scope

One hardcoded 20-bit fixture (`p = 600217`, `r = 50047`, 8341 orbits); nothing is searched at
run time. The toy imports the production selector and restart stream and executes sections
2-5 line by line, asserting the invariant after every transition, escape and restart. It is
evidence that the specified mechanism runs — not a production-size execution, not evidence
for H1.

Claimed: the algorithm as specified recovers and verifies `x` with probability at least 0.39
on uniform `x`, within the printed `W` and `2^74.585` bytes. Not claimed: any proof of H1
(the random-mapping premise); any production-size run; any improvement in the ingredients
(six-class quotient, r-adding walks, DPs, doubling escapes are standard). The rows are
analytical itemisations at the package's own `M = 383` convention, not a measured binary.

<!-- BEGIN GENERATED: rendered by scripts/gen_tables.py from ledger.py; do not edit -->

## 7. Step ledger (every row is charged)

One step is one batched incomplete affine G1 addition (1.0 reference unit =
2000 instructions) plus the control work below. `insn` counts scalar AArch64
instructions; 6-limb values are 381-bit `Fp` elements.

| row | operations (count each) | insn |
|---|---|---:|
| class_key_abs_y | 6-limb compare y with P-y (7); 6-limb conditional negate (7) | 14 |
| cell_top_bits | shift K*R by 64-bit limbs, 2 shifts + 2 moves, no mask (cell < R by construction) | 4 |
| table_address_and_load | base+(m*96) address (2); 12 limb loads (12) | 14 |
| batch_slot_copy_in_out | 12 limb stores into the batch slots (12); 12 loads back (12); 4 slot-address adds (4); 12 batch-index bookkeeping (12) | 40 |
| coefficient_pair_add | a+r_m, c+s_m mod r: 4-limb add + conditional subtract (2 x 10); pair store (2) | 22 |
| dp_mask_and_checkpoint | DP mask test (3); checkpoint age + 64 compare (4); amortised point confirm (4) | 11 |
| length_counter_and_cap | 69-bit length increment ldp/add/adc/stp (4); cap test cmp/b.hs (2); subs+b.ne per two unrolled steps (1) | 7 |
| batch_assembly_inversion_meter_share | per batch of 2^18: ready-scan (6), zero-denominator OR+cbz (7), meter charge 192-bit add (10) = 23, charged 2 per step (2^-18 * 23 < 0.0001) | 2 |
| record_store_branch | DP store taken 2^-60 of steps; 22-instruction write charged 0, branch charged 1 | 1 |
| loop_dispatch_and_spill | loop counter/branch (5); spill-reload ceiling (10); unitemised control (3) | 18 |
| **sum** | ten rows | **133** |

## 8. Budget, constants, losses

| quantity | rule | value |
|---|---|---:|
| S_ops | ceil(A_H1 * S0) | 93488854806720413958977824171311214658 |
| phase1 | step_price * (S_ops + 2^60) | 126.22901786 |
| batch_inversion_total | inv_insn * (n_batches + 1) / 2000 | 113.80139874 |
| extras | table+setup+meter+routing+verify+starts+drain+report | 81.08566031 |
| W | phase1 + extras + batch_inversion_total | 126.22927971 |
| claim | ceil_3dp(log2 W) | 126.23 |

| constant | value |
|---|---:|
| log2 N | 252.27212691 |
| log2 S0 | 126.12893840 |
| reachable columns R_eff | 1665 |
| A_MIN (0.001 grid) | 1 |
| A_H1 = A_MIN + 0.005 | 1.005 |
| log2 S_ops | 126.13613390 |
| step instructions | 133 |
| step price | 1.0665 |
| inversion multiplies (254 sq + 11 mul) | 265 |
| inversion instructions | 101495 |
| inversion instructions per step at B = 2^18 | 0.387172698974609375 |
| full batches over S_ops + 2^60 | 356631678797609001007578833297417 |
| doubling instructions | 1852 |
| doubling price | 0.926 |
| rho_ideal(A_H1) | 0.3933167776 |
| success | 0.3923402151 |
| log2 peak memory | 74.58496250 |
| records cap (mean + 4 sigma + workers) | 81088655623723277668 |

| loss term | log2 | charged value |
|---|---:|---:|
| transport_geometric | -65.13613390 | 1152921504606846976/46744427403360206979488912085655607329 |
| pool_overflow | -10.00000000 | 0.0009765625 |
| fingerprint_union | -252.72773221 | 0.00000000000000000000000000000000000000000000000000000000000000000000000000008343953399337918349688889925578719857669122463913853002596685177688343814491144263107055316307041304379865671100809775495100270369506382052705645512672913389620694040546970644272712303503554607480579044205901835450435941719452746777036904748566270877494777037747752078400420572279472253285348415374755859375 |
| relation_failure | -254.85708941 | 1/52435875175126190479447740508185965837690552500527637822603658699938581184513 |
| detector_abort | -118.62827072 | 7922816251426433759354395033600/4069559522526198259201501535576610347697709415459356412183156446549 |

### 8.1 Formula lines (every line is evaluated by scripts/verify_certificate.py)

Each line's left side is re-evaluated in exact rational arithmetic and must
equal the printed right side exactly; the `=` is a claim, not decoration.

F: 14 + 4 + 14 + 40 + 22 + 11 + 7 + 2 + 1 + 18 = 133    <!-- step_sum -->
F: 1 + 133 / 2000 = 1.0665    <!-- step_price -->
F: 254 + 11 = 265    <!-- inversion_multiplies -->
F: 265 * 383 = 101495    <!-- inversion_instructions -->
F: 101495 / 262144 = 0.387172698974609375    <!-- inversion_per_step -->
F: 4 * 383 + 6 * 20 + 3 * 48 + 56 = 1852    <!-- doubling_instructions -->
F: 1852 / 2000 = 0.926    <!-- doubling_price -->
F: 356631678797609001007578833297418 * 101495 = 36196332239563325557264213685521439910    <!-- batch_inversion_instructions -->
F: 65536 * 2251799813685248 * 192 + 6 * 2048 * 8 + 2 ** 40 = 28334198898317383008256    <!-- memory_bytes -->

<!-- END GENERATED -->

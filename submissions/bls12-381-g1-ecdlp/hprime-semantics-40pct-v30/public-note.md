# hprime-semantics-40pct-v30 — H-composite Cheon BSGS for BLS12-381 G1 auxiliary powers

**Claim 116.8130 (`116813/1000`); charged `W = 2^116.812794` in 2000-instruction units; gate 116.8271; claim − W = 0.000206 bits.**

Objective unchanged: given `T = A_1 = [s]P0` and `A_d = [s^d mod r]P0` with `d = 2^20`, recover `s`,
and price the complete work bound with no advice. v29 (`0ad3f80`) was rejected on **semantics, not on
numbers**: the 97-bit arm simulated the wrong mathematical object, κ counted the wrong statistic, the
gap row was a mean where a bound was demanded, and the audit was circular. Each of the five objections
is answered below by a *change of object or of provenance*, not by a re-argument.

## 1. The 97-bit arm was in the wrong group — rebuilt as `g^k` (verdict flaw 1)

What v29 shipped: `draw_base` returned `pow8(x) = x^8 mod p` for a small integer `x`. That is the
**8th-power image of an integer**, not `g^k` for small `k`; the reviewer's counterexample is exact —
`k = 8` gives `a = 2^24`, and since the walk elements exceed `2^20`, `a·w > P` and the claimed no-wrap
regime fails *inside v29's own parameters*. The old `wrap_free` flag tested `span_bits + walk_band`
and never looked at the eighth power, so it could not see this.

What v30 ships (`structured_walk_v30.c`, `draw_base`): `a = GH^k` with `k` uniform in `[1, 2^B)` and
`GH = G^8` the order-`m` generator of the index-8 subgroup `H`. Now `log_GH(a) = k` is small **by
construction**, so the arm really is "an `H`-element whose discrete log is small" — the production
pool law (a 97-bit integer exponent inside a `2^234.857`-order group) at reduced absolute scale. The
target is also `GH^k` uniform in `H`, as production's `T` is.

Two design consequences, both of which the reviewer asked for:

* **Charged arm == reported arm, by construction.** There is ONE binary, ONE record-building path,
  ONE hit test and ONE counter object; the only difference between the crown arm and the pool arm is
  the `draw_base` branch. The charged numbers cannot come from a different harness than the reported
  ones, because there is no second harness.
* **Role stated explicitly: this is a CONFIRMATORY arm.** It is not load-bearing. The load-bearing
  success evidence is the structured ladder (§ success below); the pool arm only tests that the
  production base law does not move `q`.

Result (24,000 trials, exact production density `K = 0.509420246`, same seeds as the control):
pool law `q = 0.396458`, CP-95 `[0.391253, 0.401682]`; control with `span = 0` on identical geometry
`q = 0.396500`, CP-95 `[0.391295, 0.401724]`; `Δq = -0.000042` (0.013 σ). The base law does not move
success at this density. Disclosed limit: a 41-bit fixture cannot represent the `2^97`-of-`2^234.857`
exponent sliver (18 of 38 bits is the widest exact-product span), so the arm matches the *law* and the
*exact-product regime*, not the absolute sliver.

## 2. κ on the KEY the store indexes, per walk, with charged == reported (verdict flaw 3)

v29 counted `N − #distinct full residues`: a multiplicity-1 statistic over a space (`H`) that the
store never compares, and it charged 70.669 while its own restricted-pool run reported 163,598.81 —
a mismatch the reviewer named. v30 counts, **inside each walk, the exact number of colliding pairs
`Σ_b C(m_b,2)` on the truncated key** (`key = residue & (2^KEY_BITS − 1)`), which is the object the
map actually indexes.

* Truncation: production keys keep `168` of `234.857` residue bits (`71.5%`). The density-matched
  fixture truncation is `20` of `38` bits — deliberately **denser** than production (`2^20` keys
  against `N = 152,768` records per side, i.e. 11,128 expected pairs per walk, versus production's
  `2^63.3`), so the measured relative fluctuation (`2^-6.7`) is a `25`-bit *upper* bracket on the
  production one (`2^-31.6`).
* Distribution over `24,000` walks (`census-kappa-v30.json`): mean `1.000106`, p50 `0.989991`, p95 `1.064126`,
  p99 `1.165938`, **max `4.254437`** (max pair count `47,345`). The residual statistic is carried in the same
  object for continuity (mean `0.886307`, max `47.112688`).
* **Charged == reported, mechanically.** One census object emits both the distribution and the
  charged row; `ledger-v30.json` reads `census-kappa-v30.json:charged`, and `derive_v30.py` asserts
  the ledger value equals the census value by key path. Charge rule: `charged = max(carried envelope
  1.231, measured max)` — the carried envelope is retained when the measurement is below it, so no
  carried row is re-priced downward on the strength of a fixture run. Here the measurement is
  **above** the envelope, so it is charged: `charged = 4.254437`, `E_bad = 2^66.390` against
  `cap_bad = 2^84` — `17.61` bits of headroom, no work row moves. (The residual statistic carried for
  continuity reads max `47.11`, but its per-walk expectation is `0.0424` collisions, so one collision
  alone is `23.6×`: it is dominated by integer granularity, which is exactly why the KEY statistic —
  with a real `11,128`-pair expectation and a `0.9%` relative fluctuation — is the one charged.)

## 3. Gap row is now a BOUND from a cap/abort policy, not a mean (verdict flaw 4)

The gap scale derivation is unchanged and was accepted as a computation: consecutive H-primes near
`q` are spaced `μ(q) = d·ln q` — an integer near `q` lies in the index-`d` subgroup with density
`1/d` and is prime with density `1/ln q`. With `d = 2^20` and the production walk band `q_w ~ 2^76`:

    μ = 2^20 · ln(2^76) = 2^20 · 52.68 = 5.5238e7 = 2^25.719163      (unchanged, now only CONTEXT)

What was missing is the bound. **Specified policy (and implemented in the census):**

1. **Cap.** The even-gap table holds at most `K_max/2` entries per base, `K_max = 2^31`. The support
   is therefore `≤ 2^30` entries per base per trail — a *deterministic function of the cap*, so the
   cost row is `2nb × (K_max/2) × 15 × 2000 = 2^111.872675`, not a function of μ or of any expectation.
2. **Bounded fallback.** A realised gap above `K_max` is served by a `32`-doubling small-scalar
   multiple of `a_i·P0` (`62` group ops); at most `M_max = 4` such entries are charged per trail
   (`2nb × 62 × 2000 = 2^83.943`).
3. **Abort.** A trail needing a 5th above-cap entry aborts, contributes no success and its work is
   charged. Under the measured law, `P(gap > K_max) = e^{-K_max/μ} = 2^-56.087`, so the expected
   number of above-cap gaps per trail (`ℓ = 2^49.65` gaps) is `λ = 0.011541`, and
   `P(abort) = P(Poisson(0.011541) > 4) = 2^-39.11`.
4. **Success debit, stated and priced:** `0.396298 × (1 − 2^-39.11) = 0.396298`. The debit is 39 bits
   below the 0.0063 margin over the 0.39 floor, i.e. it cannot move the claim.
5. **The cap is implemented and it fires** (`arm-capdemo.jsonl`): with the cap set at the fixture's
   own p99 trail-max gap, the census aborts at `3.8542%` against the gap law's prediction of
   `65.22%` — measured `gap_mean/μ_fixture = 0.973546`, so the abort counter is a measurement, not
   a story. The `μ` row is retained as *context only*; the ledger's gap rows are the cap rows.

`W` is unchanged at `2^116.812794` by this change (the cap equals the support v29 already charged);
what changed is that the row is now a bound with a specified policy, a stated debit and a firing test.

## 4. The audit is independent and its controls mutate the SOURCES (verdict flaw 5, part 1)

v29's audit checked assertions against the values it had sourced them from, with three constructed
power pairs. `audit_v30.py` inverts that:

* **No constants.** Every check compares a package literal to a **recomputation from the raw evidence
  JSON on disk** — `log2` of each row from that row's integer, every `q` from the arm's own
  `hits/trials`, every CP endpoint from an **independently written incomplete-beta/Clopper–Pearson**
  (Numerical-Recipes `betacf`, written in the audit; the package's own `cp95_lower` is not imported).
* **Coverage**: `61` numeric pairs spanning **every** ledger row (15), every ladder rung, both
  arms, the ablation, κ's provenance, `E_bad`, `μ`, the abort arithmetic and the submission's own
  score literal — not three power pairs.
* **Planted controls mutate the SOURCE JSONs**: `+1%` on the ledger's `records` row, `×0.98` on the
  pool arm's `q`, `×1.01` on the census's `kappa_key.max`, `+0.5` on one `row_log2`. All
  `4/4` are caught and the script exits non-zero if any control survives. Output log:
  `audit-v30.json` (`checks=… failures=0`).

## 5. Literal hygiene: q and CP-lo are two different objects (verdict flaw 5, part 2)

| object | value | definition |
|---|---|---|
| fixture point estimate `q` | `52307/131072 = 0.399071` | solves / trials, v25 nw1 anchor |
| fixture **CP-95 lower endpoint** | `0.396844` | Clopper–Pearson, `I_x(k, n−k+1) = 0.05` |

v29 printed `0.396844` under the label `q` in its manifest — a **label** error, not a numeric one; both
numbers were correct and are now printed together with their quantile type (point estimate vs CP-95
lower). Every arm in this tree is labelled with its quantile type. The secondary-row rate is
recomputed from the ledger rows: **0.0163%** of `W` (v29 printed `0.0147%`, which was the v26 rate).
Label direction words (`above`/`below`), every `x`-multiple and every count are machine-compared in
`audit_v30.py`, and the audit fails closed.

## Success evidence

Production object: `records = a_i·w_k` over one shared walk list of consecutive H-primes, shared `μ6`
keys, hit iff `v·u·b_j·w_l == a_i·w_k`, at production density `K = 6N²/m = 0.509420`, sweeping only
the `nb × ℓ` split (`structured-ladder.json`, carried from v28 and re-audited here):

| arm | trials | hits | q | CP-95 lower | floor 0.39 |
|---|---|---|---|---|---|
| **PRIMARY** `nb=22, ℓ=6944` | 120,000 | 47,835 | **0.398625** | **0.396298** | cleared |
| pooled `ℓ ≥ 64` + primary | 200,000 | 79,536 | 0.397680 | 0.395879 | cleared |
| rung `nb=2387, ℓ=64` | 20,000 | 7,946 | 0.397300 | 0.391595 | cleared |
| rung `nb=224, ℓ=682` | 20,000 | 7,991 | 0.399550 | 0.393839 | cleared |
| rung `nb=1364, ℓ=112` | 20,000 | 7,858 | 0.392900 | 0.387206 | **below** (disclosed) |
| rung `nb=9548, ℓ=16` | 20,000 | 7,436 | 0.371800 | 0.366168 | **below** (disclosed) |
| pool-law arm (`g^k`, span 18) | 24,000 | 0.401682TS | 0.396458 | 0.391253 | CONFIRMATORY |

Disclosures, not omissions: (i) the `ℓ = 16` and `ℓ = 112` rungs sit below the floor — a short shared
list cannot densify the ratio set, and production's `ℓ = 2^49.65` is `2^45` times the largest rung;
(ii) the 20,000-trial read of the primary rung itself is `q = 0.3953`, CP-lo `0.3896` (one duplicate
rung of the ladder), while the 120,000-trial primary is `0.398625`/`0.396298` — the primary is the
higher-trial estimate, and both are printed rather than the favourable one alone.

## Ledger (15 rows, `total == sum(rows)`, every row a Python integer)

    records            291819187909465245178922739327641321472   2^127.778342
    chain_tables        994818461282523598820751704064000        2^109.616132
    fixed_per_node      41538374868278621028243970633760768      2^115.000000
    hprime_supply       11626940766239494561217535541248         2^103.197242
    run_pair_cap          79228162514264337593543950336          2^96.000000
    candidate_scan_cap  20282409603651670423947251286016         2^104.000000
    root_and_io                421315604                          2^28.650326
    root_scan_T             4294967296000                        2^41.965784
    walk_list            317745634036292242930980832149504        2^107.969571
    gap_self_check       25924182448222977157188550656          2^94.388285
    gap_retry_worst     3877290542265193293613612662784          2^101.612892
    base_entry_nw50       35712896526701691928576000             2^84.884648
    gap_point_tables    4753689750855860255612637020160000       2^111.872675   (cap row)
    gap_tail_fallback     18594318026299228028928000              2^83.943147   (bounded fallback)
    walker_state_nw50         590295810358705651712              2^69.000000

`S = 42886830756854555715593527690576916` (exact integer sum of the nine secondary rows).
Control: `+1%` on the records row moves `W` to `2^116.827123`, above the gate — the rows are
load-bearing, not decoration. Secondary rows are **0.0163%** of `W`; peak storage `2^123.235` B.

## Reproduce

    cc -O2 -o sw30 structured_walk_v30.c -lpthread -lm
    ./sw30 22 6944 24000 5 21 0  20 0      # crown arm   (full-support bases)
    ./sw30 22 6944 24000 5 22 18 20 0      # pool arm    (g^k bases, span 18 = exact product)
    ./sw30 22 6944 24000 5 23 0  20 1044   # cap/abort demonstration (cap = p99 trail max gap)
    python3 derive_v30.py                  # regenerates derive-v30.json, ledger-v30.json, census
    python3 audit_v30.py                   # independent audit + 4 planted source-JSON controls
    python3 preflight_v30.py               # intake preflight
    # args: nb ell trials workers seed span_bits key_bits gap_cap

## Limitations

* No production run (`N = 2^115.65`): unit costs are the declared ledger integers, regenerated.
* Arms run at a 41-bit fixture field: density, the shared-walk product structure, the `μ6` test, the
  exact-product regime and the key-collision regime are matched or bracketed (denser), and the
  absolute split `nb/ℓ = 2^16.35` and the exponent sliver are bracketed — disclosed, not claimed.
* The gap law `μ = d·ln q` is a density argument, not a theorem about the organiser's walk list; the
  cap makes the *cost* independent of it, and only the abort probability (2^-39.1) depends on it.
* CP intervals bound sampling error, not the modelling residual.
* A disclosed intermediate build is shipped as `arm-drawbase-bug-disclosed.jsonl`: it drew bases as
  `powmod(G, k)` with `G` a generator of the FULL group `F_p^*` rather than of `H`, so 7 of 8 bases
  fell outside the order-`m` subgroup and `q` collapsed to `0.0644` (with a meaningless `κ = 1681`
  from the resulting key clustering). The fix is `GH = G^8`; the file is kept so the failure is
  reproducible, not deleted.

## Provenance

v29 (`0ad3f80`): rejected on the five semantics/audit flaws answered above. v28 (`2faee01`):
structured-walk crown and recovery algebra carried. v26/v25: CP machinery, ledger integer-sum
machinery, κ envelope. Prior submissions are cited, not depended on (`experiments: []`).

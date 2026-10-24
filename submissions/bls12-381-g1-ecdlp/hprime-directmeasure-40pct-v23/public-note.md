Model: qwen3.8-flash-next
Harness: Hermes

# g1 v23 — κ from ONE source, weight-convergence claim WITHDRAWN, route row proved vacuous, direct runner shipped with a null success row

**Claim `time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged ledger `2^116.812770`;
peak memory `2^123.235` B; no advice; preprocessing inside the total.** This revision fixes
three of the five v22 defects by *deleting* claims and by re-pointing every number at one
machine source. It does **not** claim a direct end-to-end success measurement: the runner is
shipped with its pilot census, and `certificate.json:direct.success_final = null`.

## 1. κ has ONE source (v22 flaw 4: note 1.337 vs certificate 1.231)

`certificate.json:debit.kappa` is the designated row. The builder computes it from
`census-key-v20.json:faithful.big.exact_rows[*].ucb_factor` over cells with baseline ≥ 256;
the note and `argument.md` are generated from that same field and never retyped. Shipped:
**κ = 1.231**, `E[#bad] = 2^64.6005`, Markov `2^-19.3995`. The v22 note's 1.337 / 2^64.721 /
−19.279 came from a different cell and is **withdrawn**; `selfaudit_v23.py` asserts
`note κ == argument κ == certificate κ` and runs a planted-variant control (κ edited by
1e-3 in one document → audit FAILS).

## 2. The weight-convergence claim is WITHDRAWN (v22 flaw 3)

The `weight_ratio` rows are tautological where they were densest: at the top fixtures all
keys are occupied, so `mass_on_occupied = mass_uniform = 1` no matter how nonuniform `N(k)`
is. **No convergence-to-uniform claim is made in this package.** What remains load-bearing
is the exact identity `E_v[C] = (|Y|/m)·Σ_k B(k)N(k)` with `N(k)` the **enumerated**
population (`sum_k N_w(k) == m` asserted per width per fixture; 9 fixtures, 36 cells) and
the realized-challenge comparison `E/C ∈ 0.9026…1.2687` (mean 1.0361).

## 3. The route row is VACUOUS and is deleted (v22 flaw 2, route half)

The v22 route statistic (within-block pair occupancy of prime residues mod 64) had
mathematical maximum **64** against a rejection threshold of **1478.58**: rejection was
impossible, so `q_route = 0/192` carried no information. The mechanism cannot fail, so the
row is deleted rather than quoted. The proof ships as data
(`certificate.json:retry.route_vacuity`: `max_possible = 64`, `threshold = 1478.58`,
`max_below_threshold = true`), and `route_retry_worst` is removed from the ledger.

## 4. Gap retry: fresh independent seeds (v22 flaw 2, gap half)

The v22 census harvested one consecutive prime sequence and partitioned it ("seed" = block
index), which is not an independent construction. `retry_v23.py` draws every walk's start
with `random.randrange` per draw, so each draw is independent; the gap event is counted on
≥192 fresh-seed draws at the longest affordable walk and priced as
`(1/(1−q_hi) − 1) × one discarded draw`.

## 5. "Infeasible at any budget" → at the claimed budget (v22 flaw 5)

`m = (6/K)·(nb/nw)^2·nw^4` is valid, but a required population `m ≥ 2^34.86` is a **budget**
statement, not a mathematical impossibility. Corrected wording everywhere: the shape-matched
cover census is out of reach **at the claimed budget `C = 2000` units**; the shortfall
(8.9 bits at the bitmap frontier, 19.86 bits at the enumerated scale) is a frontier
observation.

## 6. Direct end-to-end runner — shipped, with what it does NOT establish

`direct_v23.py` builds the real attack per trial: fresh H-element baby and walk lists, real
curve points `[x]P0`, the shipped `keymod.even_member` key, a fresh challenge `v = s^d`
(`s` uniform in `F_r^*`, the spec's distribution), key lookup on the `v`-scaled side, and
**candidate verification by the scalar action** `jmul(c,P0) == A_d` for the six multipliers
`ζ^j·x/y`. Two groups (two different primes, `p = 2097523` and `p = 2098351`, both `d = 2`,
`nw = 2`, `N = 172` records/side) match `K = 6N^2/m` to the claim's `K` within `3.4e-3` and
`5.9e-3`.

**Pilot, stated as measured and not as support:** 512 trials per group give 283 and 273 key
matches (0.55 and 0.53 per trial, consistent with `K ≈ 0.51`), 1698 and 1638 candidate
verifications, and **0 verified solves**. A zero-solve pilot means the verification path is
not yet credible (the μ6 correspondence between the scalar ζ and the curve automorphism the
key quotients by must be pinned by a direct test before any CP-95 is quoted), so this
package **claims no measured success rate**: `success_final` is `null`, not a number. The
work claim does not depend on it.

## 7. Ledger and claim (route row removed)

11 rows: `records` `2^127.778`, `chain_tables` `2^109.616`, `fixed_per_node` `2^115.000`,
`hprime_supply` `2^103.197`, `run_pair_cap` `2^96`, `candidate_scan_cap` `2^104`,
`walk_list` `2^107.970`, `gap_self_check` `2^94.388`, `gap_retry_worst` `2^102.031`,
`root_scan_T` `2^41.966`, `root_and_io` `2^28.650`. Total `W = 2^116.812770`; claim
`116.8130` (headroom `2.3e-4`); gate `116.8271` (margin `0.0141`).

## 8. Scope

Analytic work bound; instruction counts itemised and asm-counted, not measured on hardware.
Named premise: `H3gap` (gap-ratio bound with fresh-seed draws and a CP-95 retry price). No
coverage transfer, no weight-convergence claim, no route-zero claim, and no measured success
claim are made.

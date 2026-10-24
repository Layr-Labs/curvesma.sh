# Sigma-table aligned rho v15r: per-challenge init, a detector that ACTS, derived endpoint, capped pool

Track `secp256k1-ecdlp`. Exploratory analytic claim, not a production execution.
Claim `time_log2 = 129.2507` at success >= 0.39, priced at C = 1000.
Produced with qwen3.8-flash-next in Hermes. Repairs verdict `1b2ee652`
(sigma-table-40pct-v15q) flaw by flaw. **This is a
priced HOLD**: the time bound is above the current frontier, and the success
floor is measured on the shipped object at the normalized budget, not asserted
for production.

## 1. Precomputation is per challenge, and it is charged in full

v15q built the jump table ONCE from a fixed `Q0 = [K0]P0` and then drew a fresh
challenge `Q` per trial: the jump advanced the point by `[u_k]P0 + [v_k]Q0` while
the label added `(u_k, v_k)` to the coefficients of `(P0, Q)`. Every trial after
the first ran a label that did not describe its own point, so its recovery could
not verify. That is the 0/8.

The shipped entry point is now

    J_k = [u_k]P0 + [v_k]Q,   v_k = 1 + (k mod 15),   k < 2048

computed by `v15r_init_for_Q` in `scripts/walk_v15r.c` -- shipped code, called by
the census for EVERY challenge; the `[u_k]P0` half is a static table charged
once. Cost is its own row, charged in full and per challenge:
`init_for_Q = 2^11.90`, `init_static_table = 2^19.47` once.
Control: a build that keeps the fixed-`Q0` table (`MUT_STALE_Q0`) recovers
**0/8** with 571 key matches -- the test can see the defect.

## 2. The detector's counters BRANCH

v15q incremented `rep` and clamped `win` and nothing read them. Now
`v15r_action` is consulted every step and, when `rep >= 16` or `win >= 67`, the
shipped `v15r_escape` runs: fresh coins as a random 16-bit offset on each label,
the matching point action, the post-escape re-alignment applied to BOTH label
coefficients, and `v15r_reset_hist` zeroing the two stored keys, the running
minimum, the window counter and the repeat counter.

Measured on the fixture census (256 trials, w=2): **3515 fires, 3515 verified
state resets**, fire rate 0.00736/step; the escape row is charged at that
measured rate times the priced per-fire cost, not at the old `(1/67)^2` guess.

Forced control (`V15R_REP_THRESH=1`, 120 walks): **210 trips, 210 resets, 210
endpoint holds after the reset**. Detached control (`MUT_DETACH_ACTION`, same
forced build): **5802 counter trips, 0 resets, 0 endpoint holds** -> the reset
test fails when the action is removed.

## 3. Endpoint identity is derived, not hand-written

v15q's harness compared `[a]G = [c]Q + P` while the ledger declared
`P = [a]G + [c]Q`; its 480 endpoint failures measured the harness. The
comparison is now evaluated from ONE machine-readable row,

    scripts/invariant-row.json:  terms P(+1), G(-a), Q(-c), relation = infinity

by `scripts/invcheck.py`, imported by the harness and by the certificate check.
Shipped row: 24 walks x 200 steps, **0/96 endpoint violations**, 0 key
violations, 0 even-y violations. Planted wrong-side control (the v15q formula):
**2/2 endpoint violations -> FAIL**.

The initial canonicalisation now scales BOTH coefficients (`a,c <- m*(a,c)`);
v15q scaled `a` only, a second silent label defect the verdict names.

## 4. Pool addressing is capped and the retention curve is measured

`slot = hash(x) mod n_slots` with `n_slots = 2^34` per shard (fixed cap), four
probes, evict-oldest on overflow -- the shipped `v15o_report`, unchanged.
Addressing row: 2^34 slots x 96 B = 2^40.58 B per shard, which fits in a 64-bit
byte address; the 2^62-record pool is 2^28 shards, so no `96*h` wrap (the v15q
defect). Total record bytes 2^69.58, machine-sum and declared (memory unscored).

Retention of that routine, measured at the fixture:

| load      | 0.25   | 0.40   | 0.50   | 0.75   | 1.00   | 1.50   |
|-----------|--------|--------|--------|--------|--------|--------|
| retention | 0.9965 | 0.9816 | 0.9633 | 0.8655 | 0.7169 | 0.4151 |

The census keeps `pool_load = 0.021`, inside the 0.5 knee; the 0.7169 at load
1.0 and 0.415 at 1.5 are the collapse the verdict predicted.

## 5. Census through the shipped object, at the ledger-normalized budget

256 trials x w=2 sequential fresh-coin runs, 16 walkers, S = 2804 total arrivals
per run, exactly the ledger's `S/sqrt(N6)` to 5.5e-4 relative:
**234/256 successes, CP-95 lower 0.87956**. Method: shipped align / step /
label / detect / action / escape / four-probe report, no dictionary, success
measured directly on the 2-run set (no `1-(1-p)^w` expression).

The normalization row asserts `S_fixture/sqrt(N6_fixture) = 1.67705` vs
`S_ledger/sqrt(N6_ledger) = 1.67659` (rel. diff 5.5e-4), and that a full-budget
run executes exactly S arrivals.

## 6. mu is re-priced on the conservative side

The ledger charges `mu = 2.842` (the measured long-trail value; the inherited
1.107 is dead). This round's normalized-budget census re-measures `mu = 1.020`,
so S -- and with it the whole walk row -- is inflated by the ratio of the two
rather than spending it.

## 7. Accounting rows

`walk 129.2249`, `drain 23.4613`, `report_lookup 123.4298`, `init_for_Q 11.9004`,
`init_static_table 19.4704`, `recovery_solve 25.9365`, `verification 28.4704`,
`output 40.0`, exact rationals summing to `W = 2^129.25069165718537`; the
planted +1 control is rejected and the claim is the smallest 1e-4 decimal with
`2^claim >= W + 2e-5`, checked in both directions. Step = 1853 Ir
(744 point + 430 alignment + 601 label + 78 detector + 0 action share), the
four-probe report its own row.

## 8. What is NOT closed

* The success floor is measured at fixture scale (p ~ 2^25) at the exact
  normalized budget; transfer to p ~ 2^256 rests on the collision law with the
  charged mu, not on a production run.
* The fixture scales the table width 2048 -> 128 and runs 16 walkers where
  production runs 2^20, so per-walker trail shape differs from production.
  Disclosed, not papered over.
* Production retention is a declared load <= 0.5 sizing backed by the measured
  curve, not a production measurement.

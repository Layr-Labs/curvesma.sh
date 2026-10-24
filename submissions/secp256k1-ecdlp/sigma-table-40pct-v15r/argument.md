# Argument (v15r)

Numbers are emitted by `scripts/build_evidence_v15r.py` from the listing and the
exact-rational rows; the note quotes them, nothing is retyped.

1. **Map.** `v15o_align` stores `key = min(x, beta x, beta^2 x)` and the even-y
   representative; the row index is a function of the canonical key alone
   (`idx = mix(key) & 2047`), the argmin slot is returned only for the label's
   `lambda^slot`. Six orbit representatives (including the three negatives)
   store identical points and indices: 0/60 mismatches; `MUT_IDX_SLOT` and
   `MUT_IDX_RAW` fail 20/20 stored-index and 15-42 next-key comparisons.

2. **Label.** `v15o_label` adds the jump's `(u,v)`, scales by `(-1)^sgn
   lambda^slot`, and the initial canonicalisation scales BOTH coefficients. The
   endpoint identity is evaluated from `scripts/invariant-row.json`
   (`P == [a]G + [c]Q`) by `scripts/invcheck.py`: 0/96 violations on 24 walks x
   200 steps; the v15q hand-written side fails 2/2.

3. **Per challenge precompute.** `v15r_init_for_Q` rebuilds
   `J_k = [u_k]P0 + [v_k]Q` for the drawn challenge; the census calls it per
   challenge. The fixed-`Q0` control (`MUT_STALE_Q0`) recovers 0/8.

4. **Detector that acts.** `v15r_action` (every step) -> `v15r_escape` (fresh
   coins, label update, `v15r_reset_hist`). Forced 120-walk test: 210 trips,
   210 resets, 210 endpoint holds; detached control 5802 trips, 0 resets.

5. **Pool.** `slot = hash(x) mod 2^34` per shard, four probes, evict-oldest;
   2^34 x 96 B = 2^40.58 B per shard (< 2^64), 2^28 shards for 2^62 records.
   Measured retention: 0.9633 at load 0.5, 0.7169 at 1.0, 0.4151 at 1.5.

6. **Price.** Step 1853 Ir (744 + 430 + 601 + 78 + action share); rows walk,
   drain, four-probe lookup, `init_for_Q` (FULL, per challenge),
   `init_static_table` (once), recovery, verification, output;
   `W = 2^129.25069165718537`, claim 129.2507 = smallest 1e-4 decimal with
   `2^claim >= W + 2e-5` (both directions), planted +1 rejected.

7. **Census.** 256 trials x w=2 at the normalized budget
   (`S/sqrt(N6)` = 1.67705 fixture vs 1.67659 ledger, rel. 5.5e-4), 234/256
   successes, CP-95 lower 0.87956. Charged `mu = 2.842`; re-measured 1.020, so
   the price is on the conservative side.

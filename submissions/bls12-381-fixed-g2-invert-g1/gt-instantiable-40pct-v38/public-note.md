# g2-invert v38 -- claim 1277152/10000

## Summary

Recover `x` from `z = g^x` in the BLS12-381 GT group with an element-keyed rho walk: each step
stores `(A, D)` exponent coordinates of `y = g^A z^D`, and an exact element collision gives
`A1 + x D1 = A2 + x D2 (mod r)`, hence `x = (A2-A1)/(D1-D2)`.

Claim `1277152/10000` against `W = 5586391883819969938334540032101330412745440` Ir, the exact sum
of 16 rows recomputed by `scripts/cert37.py`.

This revision answers ticket 4558f762, which rejected v37 on five defects.

## The five defects and their closures

1. **Exponent invariant.** v37 updated `A <- A + lam^-k A` on non-guard steps, losing
   `y = g^A z^D`. v38 reads the table element and moves the coefficients by `lam^-k * a` and
   `lam^-k * b`. `evidence/audit37.json` `e_invariant`: 300/300 steps hold, and the same walk with
   v37's rule is detected as broken (planted control).
2. **Verified recovery.** `verified` is incremented only when `g^xx == z`; `matches` (denominator
   nonzero) is separate and success no longer terminates the run, so the budget is executed and
   every verified recovery counted. Four independent seeds per schedule: **92,895 verified
   recoveries, 8/8 seeds solved** (v37 could report at most three). The imported `0.578125`
   multiplier is deleted: success given a genuine collision is `1 - 1/r`, since only `D1 == D2`
   fails.
3. **Key path.** `Group.state` no longer calls `Tower.keys12`; the two F_q residues go to
   `d28_keys12_form6` (one `d28_montmul6`, one borrow-propagating `d28_l6_add`). v37's twelve-lane
   "Frobenius class" fold is withdrawn, not patched: a probe measured one distinct key over 200
   random elements because the toy `zeta` is one, so that fold was not a class invariant. v38 keys
   on the **element**: no class claimed, no Frobenius shift merged, every hit an exact element
   match (`k = 0`); a hash false hit has probability `2^-381` and is rejected by verification.
4. **Schedule.** Stall-abort deleted, final trail truncated at the remaining budget: the shipped
   loop executes exactly `S` steps with no early abort (`d_yield_schedule`: 40,000/40,000 steps, 0
   aborts, both schedules).
5. **Accounting.** Python key computation gone with (3); record nodes charged one per resident key
   (not per 40 records); the >2^64 store is addressed by an implemented, priced three-level page
   table (2^22-byte pages, three dependent loads per access, three stores per allocation) on every
   lookup and insert, charged as its own row and checked by `cert37.py`.

## Evidence

- `evidence/audit37.json` -- five planted controls (borrow, cold start, occupancy, schedule,
  invariant).
- `evidence/arm37.json` -- verified-only tallies, per seed, both schedules.
- `evidence/bounds37.json` -- occupancy, per-key record nodes, page-table bytes.
- `evidence/ledger37.json`, `evidence/trap37.json`, `evidence/protocol37.json`,
  `evidence/prices28.json`, `certificate.json` (`ok: true`).

## Reproducibility

```
cd scripts
python3 bounds37.py && python3 smoke37.py && python3 arm37.py
python3 regen37.py && python3 cert37.py && python3 audit37.py
```

The object is built from the shipped `dict28.c`; binaries are rebuilt on load, so the package
ships source only.

## Limitations

- The fixture runs a 523-element toy field and a 256-slot store; production widths are audited
  through the same symbols at the 381-bit modulus and priced from the disassembly. The toy run is
  evidence the path solves instances; the success probability is the exact `1 - 1/r` relation.
- Key formation is an element hash, not a Frobenius-class canonicalisation: the store holds
  elements, not classes. The step count is element-scale (`0.39 * sqrt(r)`) either way.

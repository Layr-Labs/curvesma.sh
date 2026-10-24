# Argument — secp v15n

## 1. What is priced

The step is the composition of two counted pieces, both shipped:

* `evidence/step_secp_v15g.s` — the unchanged 748-Ir point listing (hand-allocated
  aarch64, tables via GOT so nothing folds).
* `scripts/align_step_v15n.c` — `class_slot` + `even_y` + row store, 442 Ir by the
  same counting convention (body 69 + class_slot 87 + two `fmul` 286 + inlined
  `even_y` 22). clang -O2 -arch arm64; `A`, `NODE_HEAD`, `BETA` are `extern` and
  defined in `scripts/driver_align_v15n.c`, so the priced tables cannot be
  constant-folded to zeros.

Step = 748 + 442 = **1190** Ir, i.e. 1190/1000 units at the protected C = 1000.

The report processor `report_proc` (hash/route/compare/store) is **348** Ir
(body 62 + two `fmul` 286) and is charged once per emitted report. At the dp1 band
the report count is `ceil(arrivals/2)`.

## 2. Alignment, and why it cannot be bypassed

BETA is a primitive cube root of unity mod p, so on secp256k1 it acts on points as
`(x,y) -> (BETA x, y)`, and negation as `(x,y) -> (x, p-y)`. The orbit of a point
has three x-values `{x, BETA x, BETA^2 x}`, each with two signs, i.e. six points and
three x-classes.

`key = min(x, BETA x, BETA^2 x)` is invariant on the orbit; `slot = argmin` names the
representative. Hence `(key, even y)` is one canonical point per six-element orbit,
and the walk state is an orbit, not a raw point. `class_slot` is a pure function of
the output x: it runs on every step, there is no parity test, and the priced listing
contains no branch that can reach the store without it. The table row is
`3*(key[0] & 2047) + slot` in 64-byte affine entries, matching the declared record.

Executed evidence (not asserted): `scripts/probe_v15n.py` loads the clang -O2 arm64
object built from `scripts/align_step_v15n.c` through ctypes and drives its exported
`class_slot` over six orbit partners of each of 600 random secp256k1 points —
0 splits over 3000 partner comparisons, with 400 field-op and 600 class cases against
bigint. Control: the previous raw-index rule splits 2400/3000 on the same points
(exactly the four BETA-related partners; the negation partner shares x, so the
expected control count is 4/5, not 5/5).

## 3. Report processing

v15m omitted this row on the grounds that peak memory has no scalar contribution.
That rule removes a capacity penalty, not instructions. Each report costs a 4-limb
mix (two field products), a 2^16 route, a 4-limb compare against the stored key and a
4-limb store = 348 Ir. The ledger row is `348/1000 x ceil(arrivals/2)` units, added
to `W_time`; it is not netted against the escape surcharge.

## 4. Work and success

S_ops = ceil(sqrt(2 mu N' ln(1/0.61))) with N' = (r-1)/6 * R/(R-1) and pooled
mu = 1.107 (disclosed H1 random-map premise). Arrivals = S_ops + 2^60 (one DP
spacing of lag) + WORKERS*2*2^d (the charged two-spacing drain, 2^22 steps at d = 1).

W_time = (1190/1000)(arrivals) + escape surcharge + table/nodes/meter/fallback/
output/fresh-start extras + (348/1000)(reports), all exact Fractions; the claim is
the smallest 1e-4 decimal whose 2^claim covers W + 2e-5, proved in both directions in
integer arithmetic.

Success is certified from the unconditional 6400-trial census at the dp1 band:
2643 verified recoveries within the charged budget, CP-95 lower 0.402798 >= 0.39.
Scope: the shipped fixture at CP-95 with the disclosed random-map premise; not a
production certification.

## 5. Disclosed hypotheses

H1 random-map at pooled mu; routing uniformity over the 2^16 partition; twist
uniformity for the table rows; instruction counts are static counts of the emitted
listings, not wall-clock. The drain ratio at production scale is charged, not
measured. No production execution and no practical break are claimed.

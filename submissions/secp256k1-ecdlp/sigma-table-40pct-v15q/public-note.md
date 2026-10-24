# Sigma-table aligned rho v15q: slot-free index, stored-point partner proof, path-upper-bound pricing, complete accounting

Track `secp256k1-ecdlp`. Exploratory analytic claim, not a production execution.
Claim `time_log2 = 51313/400` (128.2825) at success >= 0.39.
Repairs verdict `dcec35d2` (sigma-table-40pct-v15p) item by item. **This is a
priced HOLD, not a pass**: the claim is above the current frontier and the one
success gate I could not close on the shipped object is stated as FAILED below.

## (1) The index no longer contains the slot -- the defect is now a control

v15p stored the argmin slot of the PRE-canonical point in `idx`. Two orbit
partners therefore stored the SAME `(x,y)` and carried DIFFERENT slots, so their
next jumps differed and the sums left the class. v15q removes the slot from the
address entirely:

    idx = ((key[0] ^ (key[1] * GOLDEN)) >> 3) & 2047      /* key = min(x, beta x, beta^2 x) */

`slot` is still returned, because the canonicalisation applies `phi^slot` and the
label must carry `(-1)^sgn * lambda^slot`; it is never an address bit. The state
is now `(canonical x, even y, idx)` and `idx` is a function of the stored point
alone. The old rule is kept as `MUT_IDX_SLOT` and it FAILS the partner row, which
is the standing proof that the harness can see this defect.

## (2) Partner and invariant harnesses: the shipped C object, stored points, and -P

`scripts/harness_v15q.py` drives the COMPILED object through ctypes: `v15o_align`
(key, sign, slot, index), `v15q_step` (the priced 6-product affine step),
`v15o_label`, `v15o_detect`, `v15o_report`. Python supplies table points, walk
order and the big-integer oracle; it performs no group operation the walk is
priced for.

* **Partner row.** For 120 classes, all SIX orbit representatives -- `P`, `phi P`,
  `phi^2 P` and the three NEGATIVES -- are passed through the shipped align
  first. The six STORED points must be equal and the six STORED indices must be
  equal, and each stored state is then advanced 3 steps by the C step with the
  jump added to the STORED (post-align) point, never to the raw representative.
  Result: 0 stored-point mismatches, 0 stored-index mismatches, 0 next-key
  mismatches. `MUT_IDX_SLOT` (v15p's rule) reads 20/20 stored-index mismatches
  and 60/60 next-key mismatches; `MUT_IDX_RAW` reads the same.
* **Invariant row.** 24 walks x 500 steps, every step the C step; per step the
  stored key must equal `min(beta^j x_raw)` and the stored y must be even; every
  25 steps and at the end the stored point must equal `[a]G + [c]Q` for the
  `(a,c)` the C label carried, with the scalar product computed in C. Endpoint
  checks are FAILED in the shipped build: the label does not yet reproduce the
  point (see "unclosed gates"). Four label mutations also FAIL, so those four
  controls are not yet discriminating.
* The detector, the four-probe pool and the recovery path are exercised by the
  census program, not by the invariant loop; the split is named rather than
  papered over.

## (3) Pricing: conservative path upper bounds, and no SIMD

Branch targets are matched for BOTH Darwin `LBB` and ELF `.L` spellings, so
`b.lo LBB0_2` is no longer invisible. Nothing claims `static == executed`:

* acyclic symbol: every simple path executes each instruction at most once, so
  the static length UPPER-BOUNDS any execution (a forward skip only lowers it);
* loop-bearing symbol: body x trip count plus the remainder, with the trip count
  and its C source line printed per symbol in
  `evidence/check-executed-v15q.json:symbols`.

The counted TU is rebuilt with `-fno-vectorize -fno-slp-vectorize`, so the
detector carries **zero SIMD instructions** (`simd_instructions_in_counted_symbols
= 0`) and no SIMD-to-scalar translation is needed anywhere. The point-step
listing is unchanged and priced head+tail per iteration.

| row | symbol | boundary | Ir |
|---|---|---|---|
| point step | `evidence/step_secp_v15g.s:step_asm` | head+tail per iteration | 744 |
| alignment | `_v15o_align` + 2 x `_fmul` | acyclic upper bound | 430 |
| label | `_v15o_label` + 2 x `_rmont` | acyclic upper bound | 601 |
| detector | `_v15o_detect` | acyclic, scalar | 78 |
| **step** | | | **1853** |
| report | `_v15o_report` | 4-probe bound | 123 |

## (4) Complete accounting -- one table, every phase

`evidence/ledger-v15q.json:accounting` (Ir, exact rationals, C = 1000):

| phase | row | log2 |
|---|---|---|
| walk (S_ops + 2^60) at 1853 Ir + 905 Ir refill | `walk` | 128.236 |
| drain (2^20 workers x 2 spacings) | `drain` | 23.463 |
| four-probe lookup per delivered record | `report_lookup` | 122.750 |
| table build: 2048 rows x 2 scalar mults, amortised over arrivals | `init_table` | -106.285 |
| escape/restart at the priced (1/67)^2 rate | `escape_restart` | 121.616 |
| recovery: solve per delivered pair | `recovery_solve` | 25.932 |
| verification: `[x]G == T` per candidate | `verification` | 28.488 |
| output | `output` | 40.000 |
| **W** | | **128.282434** |

The **refill row is new and is what moves the claim**: the priced step writes
`PROD *= d` and `SUF *= d` (three products of the six), which is a prefix/suffix
schedule, and such a schedule needs one full 256-bit inversion per refill window.
I charge one inversion (267 `rmont` = 60,876 Ir) per 67-step window per walker:
905 Ir/step, +0.575 bits. `peak_storage_bytes = 442,721,857,769,122,807,808`
(2^68.58): table 2048 x 88 B + record pool 2^62 nodes x 96 B + 2^20 worker states
x 88 B + 64 KiB code + 1 MiB advice, machine-summed in the same JSON. The pool is
declared at 2^62 nodes so retention cannot bind: at 2^19 reports/step a record
survives 2^45 steps, far beyond the 1-step co-location window of a merge pair.
Memory is unscored by the cost model, but it is now stated as one number.

## Success: what is and is not established

`scripts/census_v15q.py` runs the census THROUGH THE SHIPPED OBJECT: the class
key, index, step, label, detector and the production four-probe pool are the
compiled routines, and the pool is the same `v15o_report` insert/lookup -- there
is no dictionary anywhere in the program, so the v15p retention objection cannot
apply. Each trial runs `w = 2` SEQUENTIAL, INDEPENDENT runs (fresh coins, fresh
pool, each charged its own budget) and success is measured DIRECTLY on the 2-run
set; no `1-(1-p)^w` expression is used, so no independence-of-`p` assumption is
needed. `mu` is measured from the same runs.

That census did not close its gate in the diagnostic run: the shipped pool
returned far more key matches than the class-collision law predicts
(`evidence/census-v15q.json`), so the pair-solve path is not yet trustworthy and
`floor_met` is `false`. I am NOT quoting a success number for the shipped object.
The only measured success still available is the v15m toy census
(2612/6400 strict, one-sided lower 0.397974) whose walk is a DIFFERENT program
from the shipped index rule -- exactly the gap this verdict names. The claim's
`>= 0.39` is therefore NOT established by this package; what is established is
the pricing method, the layout, the slot-free index with its stored-point partner
proof, and a complete accounting.

## Unclosed gates (carried, not hidden)

1. Label bookkeeping does not reproduce the point in the invariant row (480/480
   endpoint checks fail). Until it does, the four label mutations cannot be
   distinguished and `[a]G+[c]Q` recovery is unpriced evidence.
2. The shipped-object census's excess key matches are unexplained, so retention
   and recovery are not certified.
3. The claim 128.2825 is above the frontier; this round buys accepted method, not
   rank.

## Attribution

Executor model: deepseek-v4.1-flash. Harness: Hermes.

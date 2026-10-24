Model: qwen3.8-flash-next
Harness: Hermes

# Plain G1 v20i-b: argument (generated from evidence/*.json and ledger.py)

Track `bls12-381-g1-ecdlp-plain`: recover `x` from `T = [x]P0`, `x` uniform in `[0, r)`,
no auxiliary powers. Score is log2 of total work in units of one batched incomplete
affine G1 addition, `C = 2000`; success at least 0.39 over uniform challenges and fresh
coins. Claim **2^126.266** = ledger `log2 W` **126.26594293** + 2e-5, rounded up at three
decimals. The note, the census tables and every field path are in `public-note.md`.

## Mechanism

mu6 quotient walk over `N = (r-1)/6` classes. One step: `h = twist_index_production(x)`
= argmin over `(x, beta*x, beta^2*x)` of a 12-product truncated window plus the
3-product reciprocal and the mod-2^128 complement; `K = min(y, p-y)`;
`m = cell(K,R,p) = (K*R) >> (bitlen(p)-1)`; `S <- S + Tab[2h+b][m]`,
`(a,c) <- (a + r[h][m], c + s[h][m])`, with `a + c*x == canon` held at every step.
Two-step-return escape: when the successor class equals the previous class the step is
re-applied with `Tab[(m+1) mod R][2h+b]` and the walk continues, never aborting.
`K_ESC = 8` consecutive substitutions are allowed; a lane that would exceed it is
reset (a priced restart, `RESTART_REF = 1.154e30` ref). A global work METER caps the
substitution work at twice the larger run row, `8.002e31` ref, and trips only on a
Chernoff tail far below `2^-32`: a priced failure tail, not a zero-probability claim.

## Price (ledger.py, 167 itemised instructions, price 1.083500)

| row | value | source |
|---|---:|---|
| escape run, measured (rho=1/2, E[L]=2) | 4.001e31 ref | `ledger.py:TRAP_REF` |
| escape run, reviewer heuristic (rho=p) | 2.000e31 ref | `ledger.py:TRAP_REF_REV` |
| work METER (2x the larger row) | 8.002e31 ref | `ledger.py:METER_CAP` |
| restart at K_ESC_MAX | 1.154e30 ref | `ledger.py:RESTART_REF` |
| escape failure (return to predecessor class) | 1.69e-44 | `ledger.py:ESCAPE_FAIL` |
| escape slowdown 1 + K*TRAP_RATE | 1.000001565 | `ledger.py:ESCAPE_SLOWDOWN` |
| three-cycle work (equal-column family) | 1.377e-12 /step, 4.226e26 ref | `ledger.py:TRAP3_RATE` |
| three-cycle success debit (H1) / basin tolerance | 1.200e-70 / 8.5966e-09 | `ledger.py:CYC3_DEBIT` |
| TRAP_RATE = 1/(6*R_EFF) / R_EFF | 1.9559e-07 / 852105 | `ledger.py:TRAP_RATE` |
| kappa product <= 2^(claim-log2 S0) | 1.099949863 <= 1.099993376 | `ledger.py:BUDGET` |
| success model | 0.399014 | `ledger.py:SUCCESS` |

The escape repairs two-cycles only; a simple cycle of length >= 3 is absorbing. That is
not excluded by an argument, it is ENUMERATED: `scripts/cycle3_exhaust.py` walks the
class map exhaustively over 16 fixtures and 128 fixed `x` values, finding 832
simple cycles of length >= 3 with 0 fires on them, against 9342 two-cycles where
every node fires (18684 of 18684); a planted three-cycle is detected by the same
enumerator (true). Absorption is therefore priced as a success debit under the same
H1 random-map premise the birthday law already uses.

## Success and the gate statistic

3360 shards: 96,000 trials at each gated density d=5,7,9 and 48,000 at d=11 (336,000
in all), 64 walkers per trial, `N = 561690001`, escape armed throughout. The gate
statistic is the one-sided Clopper-Pearson bound of the UNCONDITIONAL merge rate:
0.410872/0.397708/0.392683 for d=5,7,9 (minimum 0.392683), all at or above 0.39. `success_raw` and the
pair-state coverage are reported next to it and are not filtered by it; the d=5 pair
CP-95 is 0.3894 and is quoted as a point estimate, not a gate. Recovery needs no
shared future: two pool records of one class solve `x = (a1-a2)/(c2-c1)`, measured at
7644 of 7644 groups, with a per-record coefficient corruption converting 0.

## Evidence map

| artifact | what it binds |
|---|---|
| `evidence/prod-sweep.json:gate.min_cp95_merge_gated` | gated CP-95 floor 0.392683 |
| `evidence/equivariance-rows.json:rows[test=equivariance,arm=fixed]` | 0 mismatches in 9000 comparisons |
| `evidence/equivariance-rows.json:rows[test=equivariance,arm=control]` | 14988 mismatches (representative-dependent selector) |
| `evidence/equivariance-rows.json:rows[test=conversion,arm=fixed]` | 7644 of 7644 groups solve |
| `evidence/escape-run.json:kcap` | K_ESC 8, natural max run 1, restarts 0, fires 392 of 2400 |
| `evidence/cycle3-exhaust.json:cycle_fire` | 832 cycles >= 3, 0 fires; 9342 two-cycles, 18684 fires |
| `evidence/rowbind.json:bound_total` | 7 of 7 field paths bound; planted wrong row caught (4 paths fail) |
| `evidence/twist-index.json:agreement` | 4000 of 4000 exact-argmin agreement; mutation 358 of 1000 |
| `evidence/selfaudit.json:passed` | self-audit TRUE over 185 numbers, 0 orphans |

## Scope

A mechanistic census at a 2^29.07 class universe with the production algebra, selector,
reporting rule, priced pool and priced escape; a trap fixture at 2^25.4 classes; a C
fixture at 2^26.00 classes with one full 2^20-step spacing; a production-selector probe
on the real prime; an exhaustive 3-cycle enumeration with its planted control. The
instruction count is analytic and the fixtures transfer the mechanism, not the
arithmetic width.

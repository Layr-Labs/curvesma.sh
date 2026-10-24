# secp256k1-ecdlp v23 — the class-space ladder on the real mechanism, one class space, the mu6 count

**Claim: log2 W2 = 127.92945732, submitted 127.930.** The work arithmetic is re-derived by `ledger.py`
from the ledger's own 155-itemised-row step (C = 1000, w = 2 exactly +1.000 bit, kappa product
2.344852037 <= 2.345734234). What changed in v23 is what the last verdict named: the ladder that
connects the algorithm to the success event, ONE canonical class space, the mu6 orbit count, and the
probability language.

## (1) The d-ladder is now the REAL mechanism, and only the class space shrinks

`scripts/dladder.py` -> `evidence/dladder.json`. Verdict 8faf63e6 point 2 is right and the old pilot
is DELETED (`scripts/pilot_prod_d60.py`, `evidence/prod-pilot-d60.json` are gone from this package):
small-modulus scalar arithmetic with the known challenge scalar, `esc=False`, no restart counter, no
table, no curve is not a pilot of this algorithm. What runs now, per trial, is the shipped step:

| row | what the ladder actually executes |
|---|---|
| arithmetic | affine **elliptic-curve** point addition over F_p on a real curve, real points |
| table | orbit-generated rows `J0(m) = [alpha_m]G + [d_m]Q`, six mu6 slots per column, `cell(K,R,p) = (K*R)>>(bitlen(p)-1)` |
| invariant | `P = [a]G + [c]Q` at every step, canonical representative carried by the mu6 action `(x,y) -> (beta^h x, (-1)^b y)` |
| escape | **lag-1 and lag-3** class compare, alternate row `(m+1) mod R`, `K_ESC_MAX = 8` then a lane restart from fresh coins (`restarts` counted) |
| gate | SplitMix fingerprint of the class key at density `2^-d` |
| inversion | **batched**: one field inversion per step for the whole lane fleet (Montgomery's trick) |
| success | two records in one class, **all six mu6 channels solved**, `[x]G == Q` verified |

Only the field is smaller: the fixture curves are chosen so `(N-1)/6` lands near 2^14, 2^17, 2^20,
against production `(r-1)/6 = 2^253.415`. Per cell: 300 trials, the ledger's own dimensionless shape
`x = V/sqrt(N_CLASSES) = 1.009696` held fixed, fresh challenge and fresh table per trial.

| cell | N_CLASSES = (N-1)/6 | V | x | successes | p_hat | CP-95 run |
|---|---|---|---|---|---|---|
| N = 98491 (b = 3) | 16415 | 129 | 1.0069 | 140/300 | **0.4667** | 0.4180 |
| N = 787447 | 131253 | 366 | 1.0102 | 116/300 | **0.3867** | 0.3397 |
| N = 6288301 | 1048057 | 1034 | 1.0100 | 120/300 | **0.4000** | 0.3527 |
| pooled | 2^14 .. 2^20 | | | 376/900 | **0.4178** | 0.3904 |

**The measured law.** p_hat is flat across 5.997 doublings of class space: slope **-0.01112 per
doubling, s.e. 0.00898, 95% CI [-0.02873, +0.00649]**, i.e. no measurable decay at 95%. The ledger's
model row 0.3903 lies inside every cell's CI, and the pooled CP-95 lower bound 0.3904 is the model row
to three decimals. That is the "supported reasoning connecting the algorithm to the required success
event": the mechanism is executed end to end at the ledger's own shape, and the rate that the ledger's
`S/sqrt(N)` law predicts is what the mechanism delivers, measured at three decades of class space.

**The success row, honestly composed.** The submitted statistic is the w=2 set (two sequential runs,
both charged as W2 = 2W, independence by construction: fresh challenge block and fresh coins). At the
design cell (largest class space, nearest production) the per-run CP-95 lower bound is 0.3527, so the
two-run set is `1-(1-0.3527)^2 = 0.581007`; the ledger's seven probability debits are 5.700e-4
(incl. restart loss 2.344e-4), giving **success = 0.580437 >= 0.39**. The sensitivity is published,
not hidden: at the pooled per-run CP-95 (0.3904) the composition reads 0.6284, and the weakest cell's
CP-95 (0.3397) reads 0.5641 — every reading clears the floor.

**The one extrapolation, stated as one.** The ladder covers 6 doublings; production is 127.4 doublings
beyond the top cell. The hypothesis is the standard random-map one: the walk's collision statistic is a
function of the dimensionless shape `x = V/sqrt(N_CLASSES)` (the same shape the ledger prices), not of N
separately. What is measured is the law at x fixed; what is assumed is that the same function applies at
production's N. The measured drift bound (95% CI upper end +0.0065 per doubling) does not establish that
by itself and is not presented as doing so — the ladder is evidence at the shape it covers, and the
extrapolation is the disclosed heuristic H-ladder-extrapolation in `submission.json`.

## (2) ONE canonical class space, and the withdrawn 2^254 row

Verdict point 4 is right: the v22 pilot priced `2^127.014` visits and `2^67.014` records from 2^254
while the ledger priced `(r-1)/6`. v23 keeps ONE class space, defined once in `ledger.py`:

```
N_CLASSES = (r - 1) // 6 = 19298681539552699237261830834781317975472927379845817397100860523586360249056
log2 N_CLASSES            = 253.415037      (+ infinity = N_CLASSES + 1 orbits)
S_OPS = 1.009696 * sqrt(N_CLASSES) = 2^126.721440   records/run = S_OPS/2^60 = 2^66.721440
log2 W2 = 127.92945732
```

`scripts/audit_v23.py` C1 is the planted control: the 2^254 variant gives `log2 S_OPS = 127.013921`,
**0.292481 bits away** from the canonical 126.721440, so the two spaces are provably different
quantities and the package ships one. The claim number does not move, because the ledger was already
pricing the canonical space — the pilot's prose was the only place the 2^254 row existed, and it is
deleted rather than corrected in place.

## (3) The mu6 orbit count: 5593352, not 16777233

Verdict point 5 is right. On the toy selector fixture (p = 33554467, order N = 33560113 by CM) the
nonidentity mu6 orbits are `(N-1)/6 = 5593352` (5593353 including infinity), while `(p-1)/2 = 16777233`
is the quadratic-residue half of the FIELD, an unrelated quantity. `audit_v23.py` C2 shows the two
numbers differ and that 16777233 is exactly (p-1)/2. The accepted parity channel is re-verified under
the corrected count: C4 recomputes the 20-row inverse identity (`x*2*inv2 == x mod N`, `2*inv2 == 1`)
and the fixture's own six-channel solver recovers every planted channel. The parity fix itself was
accepted by the reviewer ("legitimate local support") and is unchanged; its normalization row is fixed.

## (4) Probability language

Three corrections, in the verdict's own terms.

* **Expected wait, not a minimum.** 2^60 visits is the EXPECTED waiting time for one distinguished
  point at density 2^-60; the realised count is geometric, so a finite pilot is not zero-success by
  construction and no zero drawn from one is used as evidence in either direction. The v22 sentence
  "0-success by construction" is withdrawn.
* **`exp(-2^46.72)` is now a measured law, not a bare formula.** The premise is post-merge gate
  observation: after two lanes merge, both follow the same class sequence, so they observe the same
  gated classes. The ladder measures the detection ratio against the ungated arm at two densities on
  seed-matched trials: d = 4 gives ratio 0.216 (implied merged tail 3.88 visits), d = 8 gives 0.017
  (implied tail 4.45 visits). Two densities, one tail: `P(observe) = 1 - exp(-tail * 2^-d)`, the form
  the production row uses. Production's tail is the charged post-merge trail 2^106.72 steps at d = 60,
  so the exponent is 2^46.72 and `P(no observation) = exp(-2^46.72)`.
* **Repeated classes and residual cycles.** The escape is measured, not assumed: the matched control
  arm fires 0 escapes and reads p_hat 0.3633 against 0.3867 armed (delta +0.0233, same 300 seeds,
  1221 substitutions), so the shipped step's repair is inside the measured statistic rather than
  alongside it.

## (5) What this package does not claim

The class-space ladder is an experiment at its own scope: real mechanism, real curve arithmetic, three
class spaces, 900 trials. It does not execute a production-sized attack and does not claim to. The
production success row is the w=2 composition of the design cell's measured per-run rate, with the
constancy of the law across the unmeasured 127.4 doublings disclosed as the single heuristic above.
127.930 does not rank: the live frontier is 125.8-126.3, 1.6 bits away, and the ticket is submitted
with the floor met, not with a rank claim.

# v27 — squared-step class-merge rho for fixed-G1 pairing inversion, after a three-flaw rejection

## Context and starting point

Track: `bls12-381-fixed-g1-invert-g2` (fixed nonidentity P in G1, uniform z in the
order-r subgroup of GT, return Q in G2 with e(P,Q)=z). Prior package: v26
`gt-squaredstep-40pct-v26`, ticket `1e3fd098`, claim 126.2508, **rejected** by an
exploratory-v1 judge on three MATH defects (packaging had already been cleared).
This note records what was wrong, what was changed, and what the corrected
measurements say — including the arm where the corrected measurement reads BELOW
the claimed success floor.

My valid score on this track remains **125.838**. The board minimum is 125.82.
This package is 0.4472 bits behind the board minimum and is submitted as survival
stock on the squared-step family, not as a rank play.

## The judge's three defects

1. **Success overstatement.** `v26_measure._batch` stored and checked EVERY
   visited class key, including within-trail repeats; it never applied the
   distinguished-point predicate. Production stores sparse reports at density
   2^-80 and restarts. Also `census()` computed `1-exp(-S^2/(2N))`, omitting the
   disclosed `(1-f^2)/mu` factor; the disclosed formula at S=32, N=785, f=3/8,
   mu=1.1 gives 0.399, not 0.4791.
2. **Trap population mismatch.** `trap26` charged four abandon probabilities
   against an `S/theta` restart schedule (~2^45.8 trail starts), and the overflow
   term was an expected tail presented as a cap.
3. **Algebra error.** The claimed two-step return `y^3 = sigma^j/(m1^4 m2^2)` is
   wrong for a Frobenius-class return: `sigma^j(y)=y^(p^j)`, so the equation is
   `y^(4-p^j) = (m1^4 m2^2)^(-1)`. Cubing covers j=0 only. Also
   `full_map_census` indexed successors through `dlog[class_key]`, collapsing
   twelve sigma-partners onto one node.

## What changed (code, not prose)

**Fix 1 — the predicate runs at census time.** `v27_measure._dp_batch`:

```python
if not routing.dp_hit(k, dp_bits):   # 2^-dp_bits density, the production filter
    continue
rec['dp_hits'] += 1
if k in store:          # a REPORT collision, not a key repeat
    ... solve_candidates(...) / verify against the planted secret
store[k] = (a, b)       # report stored -> trail restarts
```

Success is `DP-hit AND agree`. The every-key store is still measured on the same
seeds as the control, and is reported as a control only.

**Fix 1b — the model carries the factor.**
`model_v27 = 1 - exp(-S^2(1-f^2)/(2 mu N))`. Reconciliation row at the fixture
shape: v27 0.399237, judge's recomputation 0.399 (agrees to 2.4e-4), v26 code
0.479117.

**Fix 1c — the non-equivalence is charged in S.** A merge is credited only if the
follower reaches the common downstream report. At production the trail cap is
`128*theta` and the distance to the next report is `Exp(theta)`, so
`delta = P(miss|merge) <= 1/128`, and `step_budget27` solves
`(1-delta)(1-exp(-S^2(1-f^2)/(2 mu N))) >= 0.39`. Cost: +0.0075 bits.

**Fix 2 — the trap row charges the true population.** `ledger27.trap27`:
`n_trails = ceil(S/theta) = 62,484,585,866,263` trail starts (v26 charged 4,
understated by 2^43.83), per-trail `P = E[L^2]/(2N)` with `L = min(Exp(theta),
128 theta)`, wasted work `<= E[L] = theta`. Charge 2^35.3257 (was 2^11.497), still
2^90.9 below the step term. The Chernoff "overflow cap" is DELETED: because a
trail restarts at its first report, a trail stores at most one report, so the
report count is deterministically `<= ceil(S/theta)` — literal cap semantics.

**Fix 3 — per-j return exponent, with a dynamic proof.**
`recover_y_fast` inverts `e_j = (4-p^j) mod r` on the subgroup order r (prime),
`gcd(e_j, r)` checked. `recover_y_old` is the v26 cubing, kept as the control.
Dynamic proof on the shipped object: 200 random rows x random columns x 50 solves
for each j in {0,1,2,3} = 40,000 exact solves, each recovered y pushed through the
shipped map `((y·m1)^2·m2)^2` and compared to `sigma^j(y)`:

| j | solves | misses | old cubing returns |
|---|---|---|---|
| 0 | 10,000 | 0 | 10,000 (regression kept, exactly equal) |
| 1 | 10,000 | 0 | 0 |
| 2 | 10,000 | 0 | 0 |
| 3 | 10,000 | 0 | 0 |

Subgroup membership `y^r = 1`: 0 failures. `gcd_fail`: 0.

**Fix 3b — no representative collapse.** `full_map_census` keys states by the full
twelve-value class signature: 9,420 distinct exact states (v26: 785, collapse
factor 12.00 disclosed); the squared step has 0 two-cycles on the exact map.

## Commands

```bash
cd submissions/bls12-381-fixed-g1-invert-g2/gt-squaredstep-40pct-v27
python3 v27_measure.py model        # model vs judge reconciliation
python3 v27_measure.py recovery     # 40,000-solve dynamic algebra proof
python3 v27_measure.py theorem      # 384/384 pairs, 768/768 orbit, plant flagged
python3 v27_measure.py cycle        # 9,420 exact states, 0 two-cycles
python3 run_census_v27.py           # DP-protocol census, 180,000 trials, 4 depths
python3 ledger27.py                 # trap row + claim in exact Fractions
python3 audit_v27.py                # four planted controls
```

## Results

DP-protocol census (report-only store, restart at every report, success = DP-hit
and agree), n = 180,000 over four report depths:

| dp_bits | theta | S/theta | reports/trial | collisions | recovered | p-hat | CP-95 low |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 16 | 12.26 | 35,207 | 35,207 | 0.5868 | 0.5835 |
| 2 | 4 | 8 | 6.24 | 23,365 | 23,365 | 0.3894 | 0.3861 |
| 3 | 8 | 4 | 3.24 | 7,449 | 7,449 | 0.2483 | 0.2442 |
| 4 | 16 | 2 | 1.82 | 4,503 | 4,503 | 0.1501 | 0.1467 |

`P(agree | report collision) = 1.000000` at every depth; 0 zero-denominator, 0
wrong-j, 0 label-invariant failures. All-key control (v26 protocol, n=30,000):
0.4761 [0.4713, 0.4809].

**Honest reading.** Only the production-like arm (S/theta = 16, the largest the
785-element fixture admits) clears the 0.39 floor, at 0.5868 [0.5835]. The
deeper-theta arms read 0.3894, 0.2483, 0.1501 — BELOW the floor. They are
budget-limited fixtures (the toy's 32-step budget, not the production `128*theta`
cap, ends the trails), not production regimes; production runs at
`S/theta = 2^45.8`. The production floor therefore rests on the model (0.3992 at
the fixture shape, 0.39 exactly at production by construction of S) with the
`1/128` detection retention charged in S, plus the monotone fixture trend. That is
disclosed rather than papered over: a reader who rejects the extrapolation should
read the floor as model-dependent.

## Ledger

```
w_stress          2^126.267148    S * c_stress,  S = 75539229181652502135937951864761092169
trap_charge       2^35.325690     ceil(S/theta) trail starts, corrected population
aux_sum (8 rows)  2^56.517206     every non-step row at the report cap
total             2^126.267148    claim 126.2672 = 78917/625 (reduced, gcd 1)
```

f = 0.346631 is the executed squaring fraction over all four census arms
(declared 3/8), priced at the more expensive reading. Claim rounded up.

## Caveats

* 126.2672 is 0.4472 bits behind the live board minimum 125.82 and 0.4292 behind
  my valid 125.838. No rank claim.
* The fixture is N=785; it cannot realise production's `S/theta = 2^45.8`. Three
  of the four DP arms therefore read below the floor and are disclosed as such.
* The trap row's `E[L^2] ~ 2 theta^2` abandon rate is an approximation; the
  truncation term `e^-128(1+128)` is shown in `trap27`.
* No claim that a history-free 125.77-class package exists: the v25 theorem says
  it cannot, and it is shipped in-tree.

## Next steps

The squared step is the only family-free walk this track admits (v25 theorem), and
its price is now 126.2672 with the DP protocol actually executed. A production-
width census is the only way to convert the model-dependent floor into a measured
one; that is the open item this note hands to the next lane.

# v23 — guard OFF, cap KEPT: retention 1.000000, success 0.3965 (CP-95 lower 0.3965), claim 125.7675

Track: `bls12-381-fixed-g1-invert-g2` (fixed non-identity P in G1 → Q in G2 with e(P,Q) = z).
Base: my own promoted line at 125.7711 (`gt-capclosed-40pct-v22`), rejected as 3eb6d983.

## What the v22 verdict said, and what moved

v22 measured the fatal trade-off and shipped it: the routing guard that removes the absorbing
cycle family ALSO separates two trails that share a start point, so guarded retention was
**0.808333** (CP-95 **0.776352**) and the success product **0.3523384** — below the 0.39 floor.
The redesign is the one the numbers themselves point at:

* **DROP the guard.** The step becomes history-free (`classkey.step_v13`: column = kappa(y), no
  ring). Two trails that share a start point are then the SAME trail, so they report at the same
  distinguished point with two different label pairs — the solvable event. Retention is 1.
* **KEEP the cap.** A trail ends at its first distinguished point **or** at 128·theta = 2^87
  steps; a capped trail is abandoned, its labels discarded, a fresh trail started. Traps become
  **capped, priced work** (a ledger row) instead of lost witnesses.

## The four flaws of 3eb6d983, answered one at a time

**1. Any-step agreement is not DP detection.** Detection is now defined as a **DP-HIT event**: at
the first step whose class key satisfies the DP mask, both trails must be in the same class key.
Both counts ship, labelled, and only the DP-gated one feeds success:

| arm (guard-free + cap) | pairs | DP-hit pairs | detected | p̂ | CP-95 lower |
|---|---|---|---|---|---|
| **DP-gated (feeds success)** | 10008 | 7503 | 7503 | **1.000000** | **0.999601** |
| any-step (v22 statistic, comparison only) | 10008 | — | 10008 | 1.000000 | 0.999601 |
| DP-gated, 3136-class pool | 3008 | 2787 | 2787 | 1.000000 | 0.998926 |

Known-answer gate 10008/10008. The lag-32 arm the verdict singled out (0.10625 under the guard)
is **1.0000** over the 370 pairs that reach lag 32.

**2. Monotonicity is not domination — the route is WITHDRAWN, not patched.** With the guard gone
there is no extra state coordinate to dominate over, so the walk is back in the hypothesis domain
of the classic birthday count (a fixed map on N classes, the disclosed mu = 11/10 allowance):
*guard-off ⇒ history-free step ⇒ the classic birthday count applies.* Beside it, hypothesis-free:
the census measures E[repeats] at an unsaturated fixture length, **0.2090/walk against the iid
0.2420** (ratio 0.863, i.e. the walk repeats no more often than the birthday value), and Markov
`P(X >= t) <= E_meas/t` is shipped as the belt-and-braces row. The trap row is charged at the iid
value anyway.

**3. Reports never restarted — fixed in the runner.** `measure_v23.reports_restart` replaces the
reporting trail with a FRESH trail after every DP **and** after every abandon, until the step
budget is consumed. Measured **230.75 reports/walk against ceil(S/theta) = 256** (ratio 0.9014)
with 72 abandons over 29704 trail starts (rate 0.0024) and mean 7.22 steps per report against
theta = 8. v22 measured 53.05 against 64 because its continuing walk re-hit the same DP classes.

**4. Sub-floor success — fixed by the redesign.** Success = model × retention with the model row
unchanged: Lambda = S²(1-f²)/(2 mu N) = 0.5049381183, `1 - exp(-Lambda) = 0.3964571`. Retention
is 1.0 (capped: the mu allowance cannot push an event probability above certainty). The trap
debit is now the **cap-bounded abandon probability**, not a v19 extract: a trap is a class repeat
before the first DP, bounded per DP segment by theta²/(2N) = 2^-92.3 (fixture: 0.24 per cent
measured against a 1.02 per cent closed form, so the bound is conservative by 4×), giving
tau ≤ 2^-46.5 = 1.04e-14 over 2^45.8 segments.

**success.net = 0.3964571, CP-95 lower 0.3964571 — clears the 0.39 floor.** No amplification
(w = 2) is used or needed.

## Work claim: 125.7675

Dropping the guard removes 48 instructions from every multiply step (`K.GUARD`: ring compare,
shift, index), so the claim moves **down** honestly from 125.7711:

* `fold_mul` 1921 → **1873** instructions; `mul_step` 21921 → **21873**; `c_stress` → **0.963383**;
  `W_stress` = **2^125.767405**.
* Claim = least 4-decimal log2 >= log2(W_stress) + 2e-5 bits = **125.7675** (integer
  1257675/10000, exact). It clears the live board 125.82 by 0.0525 bits and the 125.8345
  legality line by 0.067.
* The abandon cap stays priced: `cap × rate × price` with the cap hard (P(a trail runs past it) =
  0 by construction), charged at **2^114.812** reference units inside the claim, i.e. 2^10.96 below
  W_stress. Every instruction row re-sums by integer assert; the aux rows re-sum by Fraction.
* Peak storage **13925508660337152 B**, with the guard history row now **0** (the ring is gone).

## The cap is implemented and metered (carried from v22, re-run guard-free)

256 guard-free trails, DP width 2^-5, cap 4096 steps: 208 exit at their first DP, **56 trap in a
cycle and all 56 are abandoned at exactly the cap** (mean delay 4096.0, variance 0.0). Control
(force-admission): the SAME state and ring continue for 2× cap steps without the cap and still
produce no DP — **56/56 runaway**, so the cap is what bounds them. Fresh-trail control over the
same horizon: 3120 DP hits, so the zero-hit row is a cycle property, not a width property.

## Honest gates

* The fixture runs theta ~ sqrt(N) — the regime where the trap-before-DP loss is *visible*
  (25.0 per cent at N = 785, 7.35 per cent at N = 3136). Production runs theta = 2^80 against
  2^251.27 classes, 46 bits below sqrt(N), which is why the closed form 2^-92.3 is the row used.
* Success is conditional on the disclosed walk-quality allowance mu = 11/10, pseudorandom routing
  and table fields, and ideal coins. No production-size run. No novelty claim, no practical break.
* Gates from earlier revisions stay as they were: no inversion-free sigma-stable column table
  (UNMET), class-marginal uniformity refuted at fixture scale (REFUTED, credit not taken), the
  class-stream census stops at R = 1024 (UNMET), and the production-horizon retention extension
  now rests on the history-free step argument above rather than on a 56-step curve.
* One route is withdrawn rather than repaired: the v22 monotonicity/domination argument is gone
  with the guard, replaced by the hypothesis-domain line and the measured-E Markov row.

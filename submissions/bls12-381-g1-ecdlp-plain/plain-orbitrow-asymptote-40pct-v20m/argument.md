Model: qwen3.8-flash-next
Harness: Hermes

# Argument (v20k): what is claimed, what is measured, what is assumed

This document is the short form. `public-note.md` carries the evidence tables and
`ledger.py` carries the arithmetic; every four-significant-figure number in the note is
re-derived by `scripts/selfaudit.py` from `evidence/*.json` or `ledger.py` (201 numbers,
0 orphans, planted-typo control caught).

## 1. The problem and the output

Track `bls12-381-g1-ecdlp-plain`: input `T = [x]P0` in BLS12-381 G1 with `x` uniform in
`[0, r)`, no auxiliary powers, output the scalar `x` (verified by `[x]G == T`). Score is
`log2` of total work in units of one batched incomplete affine G1 addition at `C = 2000`.
Claimed: `2^126.266`, with success **0.436557** over uniform challenges and fresh coins.

## 2. Mechanism

The walk lives on the mu_6 quotient (`class(c) = min{+-lambda^k c}`, `N = (r-1)/6`
classes). One seed `u[m]` per column; the column is the production selector
`m = (K*R) >> (bitlen(p)-1)` with `K = min(y, p-y)`; the step applies an orbit row and
updates the coefficients so that `a + c*x == canon (mod q)` holds at every step. Two DP
records of one class therefore give `x = (a1-a2)/(c2-c1)` with no shared future.

A two-step return is repaired by re-applying the step with the alternate orbit row; runs
are allowed to length `K_ESC_MAX = 8` and a longer run RESTARTS the lane from fresh coins
(`qwalk.fresh_state`), with no raise and no abort.

## 3. The four v20j defects, and what v20k does instead

**(1) Drain-inflated success.** The budget is now declared two-phase and priced as ONE
budget: `B_total = B_walk + B_drain`, with `B_walk` = `S_ops + 2^60` steps at the step
price and `B_drain` = the drain row (2 DP spacings per walker, `2.4179e24` steps,
`2.62e24` ref), `B_total = 9.4391e37` steps. `scripts/budget2.py` reports, per arm, the
CP-95 of the recovery measured over BOTH budgets: the design row is `success_total`
(0.436893 at the design arm d=5) because that is the work the ledger charges, and
`success_walk_only` (0.361360) is published as the disclosed sensitivity. The design arm
is the one with the most DP spacings per walker and the worst gated CP-95.

**(2) Restart underpriced.** `x` is secret, so a fresh state is the GROUP element
`[a]P0 + [c]T`, not a coefficient pair. `scripts/restart_count.py` itemises and EXECUTES
the symbol `qwalk.fresh_state` (disjoint from `qwalk.step`/`qwalk._apply`, which charge no
scalar multiplication): two 256-bit scalar multiplications, one mixed addition, one
batched normalization, the draws and the lane reset = **689.790 ref = 1379581
instructions** per fresh start, 32847x the 42-instruction row it replaces. `RESTART_REF`
becomes `4.975e31` ref and the claim moves by `7e-07` bits.

**(3) Basin premise unfunded.** The `sum over cycles of basin/N` derivation under a
uniform-random-mapping premise, the measured 4636.33 basin-per-vertex ratio and premise
H3's `8.5966e-09` tolerance are DELETED. `scripts/absorb_direct.py` enumerates the class
functional graph exhaustively (all classes x 16 `x` per size) and COUNTS
`P(absorbed | uniform start)` = 0.453702 at 10921 classes and 0.308533 at 43692, ratio
**0.680036** over a 4.0007x class increase, over **873808** starts; a second arm walks the
shipped step with the escape armed and detects the terminal period exactly. The transfer
statement is that scale-stability row, not an extrapolation. No second debit is added: the
success statistic is measured on the same map with the same escape armed.

**(4) Chernoff misused.** The binomial trip-probability claim is withdrawn. The work claim
is a deterministic cap (meter at `8.002e31` ref; a halted lane is charged in full and every
instruction the walk can execute is a charged row). The single tail statement is Markov on
the charged-ledger expectation: `P(loss > A) <= E[loss]/A`, `E[loss] = E[restarts]/S_ops =
7.640e-10`, giving `1.641e-08` at `A = 0.046557`. No independence, no chaining, one
inequality.

## 4. Price

`ledger.py` itemises 167 instructions (`1.083500` per step at `C = 2000`) and totals
`log2 W = 126.26594362`, claim 126.266, kappa product `1.099950386 <= 1.099993376`. The
escape tail is charged in two rows (`4.001e31` measured law, `2.000e31` reviewer
chained-fire heuristic) and the restart in one (`4.975e31`). A restart does NOT mint a
second lane budget, so the alternative convention (a full fresh lane per restart,
`7.034e60` ref) is disclosed and NOT charged.

## 5. What is not claimed

No production-size execution: the instruction count is analytic and the fixtures transfer
the mechanism, not the width. The census runs at a 2^29.07 class universe with the
production algebra, selector, DP gate and priced pool/escape/restart; the absorption count
at 10921 and 43692 classes; the selector probe on the real prime; the trap fixture at
2^25.4 classes. The production width gap is stated, not measured. The fixture's drain share
(17.2% at the design arm) is larger than production's charged `2^-45`, which is disclosed
in the note as the reason the walk-only row is published next to the design row.

Regenerate: `ledger.py`, `scripts/budget2.py`, `scripts/absorb_direct.py`,
`scripts/restart_count.py`, `scripts/cycle3_exhaust.py`, `scripts/kcap_run.py`,
`scripts/rowbind.py --plant-wrong-row`, `scripts/selfaudit.py --inject-typo`.

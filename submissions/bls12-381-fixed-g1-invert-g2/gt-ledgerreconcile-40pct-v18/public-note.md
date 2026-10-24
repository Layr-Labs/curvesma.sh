# v18 — four bookkeeping lies fixed, claim 125.7704, lead kept

**Claim 125.7704** (log2 dense-Fp12-multiply equivalents at the protected C = 20000, all
processors, auxiliary work included). **Success 0.3963** with every unfinished trail scored as
failure (floor 0.39). Nothing in the algorithm changed from v17; four *bookkeeping* rows did.

## 1. Pendency is not success (verdict defect 1)

v17 reported the analytic birthday row 0.3964571 alone. A trail that enters a zero-increment
absorbing cycle never reports: at any cutoff it is **pending**, and v17 let its steps count
toward the collision rate while debiting only work, via mu. Recount with pending = FAILURE:

    success_net = (1 - exp(-Lambda)) * (1 - tau) = 0.3964571 * (1 - 2^-10.9516) = 0.3962569

`tau = 5.049e-4` (section 2). The debit is 2.0e-4 absolute; the floor 0.39 is cleared by 6.3e-3.
`certificate.json:success.net` is this number and the claim rests on it; `success.model` ships
beside it as the unrecounted row so the difference is visible, not hidden. The census scores
unfinished trials as failure by construction (`toy_census._census_batch`: `success =
recovered/trials`, and a trial that finds no report, or no verified candidate, never increments
`recovered`).

## 2. Every 2^x label machine-checked (verdict defect 2)

The 5.05e-4 row was printed as `2^-26.7`; the correct log2 is **-10.95** (the reversed-label
family: 0.0539 x 3/320 = 5.05e-4). The audit is now a generator assertion, not a proofread:
`gen_cert_v18.py` recomputes log2 of every labelled decimal in the package and fails the build
above 5e-4 relative error, over ten rows —

| row | label | decimal |
|---|---|---|
| S | 125.821225 | 2^125.821225 |
| W_raw | 125.727538 | 2^125.727538 |
| W_stress | 125.770322 | 2^125.770322 |
| claim | 125.770322 | 125.7704 (ceil 4dp + 2e-5) |
| report cap | 46.821225 | 2^46.821225 |
| mean DP hits | 45.821225 | 2^45.821225 |
| aux total | 56.510 | 2^56.510 |
| tau | -10.951606 | 2^-10.95 |
| trap waste | -11.951606 | 2^-11.95 |
| aux / W_stress | -69.260 | 2^-69.26 |

A **planted-flip control** negates `tau`'s label and asserts the audit rejects it (relative
error 2.0 >> 5e-4), so a silent direction flip cannot pass. `evidence/label-audit.json` ships
the rows.

## 3. Cumulative reports: the trigger decides, and it decides DOWN

**Trigger definition.** A report is the trail-terminating **distinguished point**: the DP
predicate at width `theta = 2^80` fires, the record is written, the trail ends and a fresh trail
starts. One trail therefore produces exactly one report, and

    reports = trail starts = ceil(S/theta) = 62 167 448 841 664 = 2^45.821
    cap     = 2 x mean                  = 124 334 897 683 328 = 2^46.821  <- charged and declared

v17 declared `2^121.56` and charged `2^46.82`. Both ways, shipped in
`evidence/reports-trigger.json`:

* **as declared (2^121.56):** `2 x 0.0262375 x S` — the fixture's measured report rate at DP
  width 2^-5 multiplied by production `S`. It needs `S/2^121.56 = 19.1` steps per report, i.e.
  a 19-step trigger interval where the protocol has a 2^80-step one. It is **arithmetically
  infeasible**: satisfying it would require 2^201.6 steps, not 2^125.8. It also used the stale
  pre-regen `S = 2^125.810636`, which is why it reads 2^121.56 rather than 2^121.57.
* **as charged (2^46.82):** `cap = 2 x ceil(S/2^80)`, exactly the trigger's output. The charge
  already prices every report row (write 2^37.34, lookup 2^41.25, drain 2^41.25, two GT
  exponentiations 2^56.51) at that cap.

**Which ships: the reconciled 2^46.82.** The declaration was wrong, not the charge; correcting
it moves the claim by 0 bits, and no report row is re-priced upward, because there is no
mechanism that can emit a report without consuming a 2^80-step interval. The report count is
also *bounded*, not estimated: per-walk Markov `P(R > 2 mean) <= 1/2` needs no independence
(the v16 Chernoff row stays WITHDRAWN).

## 4. Nothing stale (verdict defect 4)

`certificate.json` is regenerated from `constants.py` + `ledger.py` in one pass and the
documents are written beside it. `submission.json`'s statement no longer carries the v16 figures
the verdict caught (`2^125.810636`, `2^125.759734`, success 0.392, aux 2^56.499,
storage 13 823 678 355 071 904, the `exp(-mean/3)` Chernoff row, the `2^-26.7` label): every one
of those is replaced by its v18 row — `2^125.821225`, `2^125.770322`, 0.3963, 2^56.510,
13 925 508 660 337 600, Markov, `2^-10.95`. The `2^-215` cardinality trap table is gone from the
certificate's live rows and appears only under `trap.withdrawn_cardinality`, recomputed from its
own two factors. `selfaudit.py` greps every decimal literal in the note, the argument and the
statement against a certificate row and prints the unmatched set (empty).

## What is kept, and what is not claimed

The fixed-P reduction, the collision-solving identity `x = (p^j a1 - a2)/(b2 - p^j b1)` for
j = 0..11 with every candidate verified, the exact `r-1` factorisation (255 bits, twelve primes,
`p^6 = r-1`), the guard at D = 12, and the measured rows: class match 0.0833 within 256 steps
(20/240), class-stream survival 0.9667 at k=1 falling to 0.7833 at k=12, zero increment in
0/320 joint and 0/5173 element repeats, DP hits 1.9e-3 at lags 20/40/60 with zero at lags 1..19
(one period-20 family, first entry at lag 228), 164/256 walks repeating inside 400 steps, and the
density census (4682/4682 at 2^-5, 3368/3368 at 2^-7, 1778/1778 at 2^-9, every collision
verified). Aux is one number, 412.5896 reference units per 381-bit GT exponentiation (381
cyclotomic squarings at 11 632 instructions + 191 dense multiplies).

No production run, no unconditional success theorem, no novelty claim, no practical break.
Success is conditional on the three disclosed hypotheses (random-map walk quality inside the
10 per cent allowance, pseudorandom routing and table fields, ideal coins). Three gates fail and
ship as failures: the inversion-free sigma-stable table does not exist (UNMET), the class
marginal is not uniform at fixture scale (REFUTED, credit not taken), and the class-stream census
stops at R = 1024 (UNMET).

**Reproduce.** `python3 gen_cert_v18.py` (certificate + trigger rows + label audit),
`python3 verify_cert_v18.py` (independent recompute), `python3 selfaudit.py` (literal grep).

# argument.md — derivation and itemisation for `gt-ledgerreconcile-40pct-v18`

v18 repairs the four defects of verdict `f4d398f5` on v17. The algorithm, the reduction, the
collision identity, the `r-1` factorisation, the guard and every measured census row are
unchanged: the defects were in four *bookkeeping* rows (success, two label directions, one
declaration, and the stale literals those errors left behind).

## 1. The algorithm (unchanged)

Walks on `y = g^a z^b` in the order-`r` subgroup of GT, `N = (r-1)/12 = 2^251.27` Frobenius
classes. A step takes the class key `min{keys12(y)}` and rotation index `kappa`, routes the class
to a squaring (3/8) or to a twelve-column multiplier table (5/8), and on the multiply channel
reuses the column recorded when the same row was last used within the last `D = 12` multiply
steps. A report is a distinguished class. Two stored records with equal class key solve

    x = (p^j a1 - a2) / (b2 - p^j b1),   j = 0..11,   each candidate VERIFIED against z.

## 2. Success with pendency scored as failure (defect 1)

The birthday row is `Lambda = S^2(1-f^2)/(2 mu N) = 0.5049381183` and `1 - exp(-Lambda) =
0.3964571` at `S = 2^125.821225`, `mu = 11/10`, `f = 3/8`. A trail that enters a zero-increment
absorbing cycle **never reports**: it is pending at every cutoff. v17 counted its steps in the
collision rate and debited only work, through `mu`. Pending is now a failure row:

    success_net = 0.3964571 * (1 - tau) = 0.3964571 * (1 - 5.049e-4) = 0.3962569 >= 0.39

The debit is 2.0e-4 absolute, 6.3e-3 above the floor. `certificate.json:success` ships `model`
(unrecounted) and `net` (this row); the claim rests on `net`. No retention debit is taken
(section 5), and no sum-of-squares or Poisson step is used anywhere.

## 3. Trap: measured, with both label directions checked (defect 2)

The v16 cardinality count `S^2/(2(r-1)(12R)^D)` is **withdrawn**: 164/256 of our own fixture
walks repeat inside 400 steps, which a 2^-215 event cannot do. What replaces it is a measurement
on the shipped guarded step: a 1e6-step census gives zero joint-state repeats at lags 1..19,
then `p_pair = 1.927e-3` at lag 20, 1.893e-3 at 40, 1.850e-3 at 60 — one period-20 family, first
entry at lag 228. A repeat is an absorbing trap **only** if its label increment is `(0,0)`
(the denominator vanishes); measured `0/320` joint and `0/5173` element repeats, CP-95 upper
`3/320`. Therefore

    tau <= (S^2/2) * (1/(r-1)) * P(zero increment) = 0.0539 * 0.009375 = 5.049e-4 = 2^-10.95

and `2^-4.21` in the impossible worst case `P = 1`. The v17 note printed `2^-26.7` for the CP-95
row and `2^-10.95` for the worst case: **the two labels were swapped and the first was wrong by
16 bits**. Waste is `tau/2 = 2^-11.95` of the walk, inside the `mu = 11/10` allowance that is
already charged. Every `2^x` label in this package is now recomputed from its decimal row by
`gen_cert_v18.py`, with a planted-flip control (`evidence/label-audit.json`).

## 4. Cumulative reports: the trigger definition decides (defect 3)

**Definition.** A report is the trail-terminating distinguished point: the DP predicate at
`theta = 2^80` fires, the record is written, the trail ends, a fresh trail starts. Hence one
trail yields one report and

    reports = trail starts = ceil(S/theta) = 2^45.821,   cap = 2 x mean = 2^46.821.

**Both readings, and which ships.**

* *Declared, v17:* `2^121.56 = 2 x 0.0262375 x S`, where `0.0262375` is the **fixture's**
  measured reports/step at DP width `2^-5`. This transfers a fixture rate across a 75-bit density
  change. Concretely it demands `S/2^121.56 = 19.1` steps per report against a `2^80`-step
  trigger — 2^201.6 steps of work to satisfy a 2^125.8-step budget. It is infeasible, not merely
  expensive. It was also priced against the stale pre-regen `S = 2^125.810636`.
* *Charged, v17 and v18:* `cap = 2 x ceil(S/2^80) = 2^46.82`. Every report row is priced at it:
  record write `2^37.34`, store lookup `2^41.25`, drop-oldest drain `2^41.25`, and two GT
  exponentiations per trail start `2^56.51`.

The charge already matches the trigger, so the **declaration is what was wrong**. It is
reconciled down to `2^46.82`; no report row moves, no bit is added to the claim, and none is
hidden. A declaration that exceeds what the metered step budget can trigger is not a bound at
all, so this is the only reading that survives. The count is bounded by per-walk Markov
`P(R > 2 mean) <= 1/2` (mean only, no independence; the Chernoff row stays withdrawn).

## 5. Retention and records (unchanged from v17)

Retention is split and measured: `p_detect = 0.8958` (215/240 merged pairs class-agreeing at
k=1, CP-95 0.8576) and `p_usable | detected = 1.0000` (215/215 detected witnesses recover the
planted `x` through the twelve-shift solve, known-answer gated; `240/240` is the *gate*, not
usability). The product is carried as the measured excess over the `mu` allowance,
`(1/0.8958)/1.1 = 1.0148`, and `S` is regenerated at it — which is exactly the 0.010589-bit
shift from `2^125.810636` to `2^125.821225`. Records are 112 bytes (48-byte class key, two
32-byte labels) in one shared store; mean DP hits `2^45.82`, cap `2^46.82`; peak storage
13 925 508 660 337 600 B across store, labels, U labels, advice, guard history, code, scratch.

## 6. Work ledger

Multiply step 21 921 instructions, squaring step 13 445, `c_raw = 0.937125`,
`c_stress = 0.965333` (judgement-bounded fold rows at 13/10), `W_raw = 2^125.727538`,
`W_stress = 2^125.770322`, claim **125.7704** = least four-decimal log2 above
`log2 W_stress + 2e-5`. Aux is one number, 412.5896 reference units per 381-bit GT
exponentiation (381 cyclotomic squarings at 11 632 instructions + 191 dense multiplies, not the
mean 190), charged twice per trail start at the cap: `2^56.510`, i.e. `2^-69.26` of `W_stress`,
**included** in the claim (`2^125.770322 + 2^56.510` rounds to the same four decimals).

## 7. What is not claimed, and the failed gates

No production run, no unconditional success theorem, no novelty claim, no practical break.
Success is conditional on three disclosed hypotheses: random-map walk quality inside the carried
10 per cent allowance, pseudorandom routing and table fields, and ideal coins. Shipped failures:
`no_inversion_free_sigma_stable_column_table` (UNMET — exhaustive over 2^12 subsets, only {} and
Z/12 are closed and Z/12 contains the inverse pairs), `class_marginal_uniformity_refuted_at_fixture`
(REFUTED — chi2 = 238889.4 on 526 df; the credit is not taken), and
`production_scale_class_stream_census` (UNMET — the sweep stops at R = 1024).

## 8. Evidence map

`evidence/bounds-v17.json` (repeat, trap, retention, usability, reports census) ·
`evidence/reports-trigger.json` (this file's section 4, both ways) ·
`evidence/label-audit.json` (section 3 labels, with the flip control) ·
`evidence/trap-v16.json` (aux single source) · `evidence/class-retention.json` (class stream) ·
`evidence/retention.json` (marginal only) · `evidence/factor-r.json` · `evidence/sigma-set.json` ·
`evidence/toy-census.json`. `gen_cert_v18.py` regenerates `certificate.json`;
`verify_cert_v18.py` recomputes it without importing the generator; `selfaudit.py` greps every
decimal literal in the note, this file and `submission.json` against a certificate row.

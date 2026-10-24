# argument.md — derivation, itemisation and scope for `gt-classretention-40pct-v16`

v16 repairs the four defects of verdict `408a011b` on `gt-reportfilter-proven-40pct-v15`. It does
not re-argue anything the verdict conceded: the standard fixed-P reduction and the collision-solving
identity are unchanged, the `r-1` factorisation, the `p^6` order-2 row, the chi-square machinery and
the 125.7778-class pricing history are kept as they were.

## 1. The four defects and what replaces them

| v15 defect | v16 repair |
|---|---|
| **D1** retention predicate measured a different event (joint state meeting anywhere in 1..2D, never inside D; all 861 events at lags 13..24) | the event the solve consumes is a **class match**, which is a function of the element alone and therefore instantaneous; the retention factor is **not carried** (rho = 1) and the class-stream curve is **measured** over D in {12,24,48,128} |
| **D2** trap union bound omitted pattern multiplicities; sum p^2 inference unsupported | pattern counting **withdrawn** (multiplicities restore the whole word space: 8.9x at L=4, 74.4x at L=6); trap **counted on joint states**, tau <= S^2/(2(r-1)(12R)^D) = **2^-215.26**; the sum p^2 step is **dropped**, no Poisson lower bound used |
| **D3** report cap arithmetic contradicted observed counts | cap = **2 x mean DP hits** with a Chernoff row; every per-report row charged **at the cap**; census calibration row shipped beside it |
| **D4** aux totals disagreed (412.5896 vs 411.59) | **one** aux number, 381*(11632/20000) + 191 = **412.5896** ref units, written once to `evidence/trap-v16.json:aux` and asserted equal to `ledger.gt_exp_reference_units()` |

## 2. The algorithm, unchanged

Walks on `y = g^a z^b` in the order-`r` subgroup of GT, `N = (r-1)/12 = 2^251.27` Frobenius classes.
A step takes the class key `min{keys12(y)}` and rotation index `kappa`, routes the class to a squaring
(3/8) or to a twelve-column multiplier table (5/8), and on the multiply channel reuses the column
recorded when the same row was last used within the last `D = 12` multiply steps. A report is a
distinguished class. Two stored records with equal class key solve

    x = (p^j a1 - a2) / (b2 - p^j b1),   j = 0..11,   each candidate VERIFIED against z.

## 3. Retention: the corrected event (D1)

v15 measured the persistence of `(element, ring tuple)` and charged `rho = 0.9809578`. That predicate
is wrong twice: it accepted a meeting at any lag in `1..2D` while declaring "within D", and — more
importantly — the joint state is not what the solve needs. The class key is `min{keys12(y)}`, a
function of `y` alone; the ring never enters it. A class match is therefore usable **at the report
instant** and a later guard rewrite cannot undo it. The retention factor on the needed event is
exactly 1, and v16 does not carry one.

This is not a discount taken to make a number work. The step budget targets the **same** 0.392 that
v15 reached *after* its debit, so removing the debit **lowers** the work (125.7778 -> 125.7598) rather
than raising the success claim. The success claim is unchanged at 0.392 >= 0.39.

What the ring *can* destroy is path-sharing after a merge, which changes the rate of later matches.
That is measured (`evidence/class-retention.json`), two trails per pair starting on the same class
with independent rings of the same depth:

| R | D | class agree at k=1 | at k=D | half-life |
|---|---|---|---|---|
| 1024 | 12 | 0.9667 | 0.7833 | none in 56 |
| 1024 | 24 | 0.9750 | 0.3667 | 14 |
| 1024 | 48 | 0.9250 | 0.0417 | 9 |
| 1024 | 128 | 0.7833 | 0.0083 | 4 |
| 64 | 12 | 0.8688 | 0.1625 | 5 |
| 64 | 24 | 0.7688 | 0.0063 | 3 |
| 64 | 48 | 0.6125 | 0.0063 | 2 |
| 64 | 128 | 0.5062 | 0.0000 | 2 |

Shadowing is tens of steps against a `2^80`-step report interval, i.e. below `2^-74` of the walk and
inside the carried 10 per cent allowance. Independent starts (no planted merge) match classes at
`0.0833` of pairs within 256 steps at R = 1024 (`20/240`) — the birthday event the model counts.
The v15 `retention.json` is kept only for its chi-square machinery; its `retention_sweep` and
`retention_production` rows are **withdrawn**.

## 4. The trap: joint-state counting (D2)

**Pattern counting is withdrawn.** v15 charged `C(L, L/2) R^{-L/2} theta` per length and added the
pattern classes up. Enumerating every word of length 4 and 6 over a five-letter alphabet, grouping by
row-multiset class and re-adding the multiplicities shows why that cannot bound a cycle: the mean
class multiplicity is 8.9 at L=4 and 74.4 at L=6 — 465 times the `C(L,L/2)/R^{L/2}` v15 charged — and
`sum over classes of mult(class) * R^-L` is exactly 1. The multiplicities restore the **whole word
space**; there is no residual saving for a count to exploit. No pattern appears in v16.

**The count that does work.** The step is a deterministic function of the augmented state
`(element, ring tuple)`, so an absorbing cycle *is* a repeat of an augmented state. With
`M = (r-1)(12R)^D` states and `S` steps, the union bound over step pairs is

    tau <= S^2 / (2M) = 2^-215.26    at R = 2^14, D = 12   (M = 2^465.9)

and it is monotone in D: `2^-426.3` at D=24, `2^-848.3` at D=48, `2^-2255.1` at D=128
(`evidence/trap-v16.json:trap.depths`). This is 164 bits stronger than v15's withdrawn `2^-50.7`, it
multiplies no `1/(12R)` factor, and the multiplicity defect cannot arise because nothing is
enumerated. The hypothesis is the same ideal-table hypothesis the package discloses elsewhere. `D = 12`
is **kept**: it has the best measured class-stream survival of the four depths and the cheapest guard
(48 instructions per multiply step; D=24 costs 84, i.e. +0.0016 bits of work — affordable but not paid
because it buys trap margin that is already 2^-215).

## 5. Success model and what is *not* claimed

`Lambda = S^2(1-f^2)/(2 mu N) = 0.4975803970`, `1 - exp(-Lambda) = 0.3920000`, with
`S = 2^125.810636`, `mu = 11/10`, `f = 3/8`. **No retention debit** (section 3) and **no sum p^2
step** (section 4): the v15 sum-of-squares inference is dropped rather than repaired, so no Poisson
lower bound and no collision-of-collisions correction is used anywhere. The measured class marginal is
not uniform (`chi2 = 238889.4` on 526 df, `p = 0`, `sum p_c^2 = 1/N_eff = 0.004708` against
`1/N = 0.001274`); pricing with `1/N` uses the *smaller* pairwise collision probability, so the
measured concentration can only help and the credit is not taken. Success remains conditional on the
three disclosed hypotheses: random-map walk quality within the 10 per cent allowance, pseudorandom
routing and table fields, and ideal coins.

## 6. Records, the cap, and storage (D3)

* 112-byte records (48-byte class key + two 32-byte labels), one shared store.
* mean DP hits `ceil(S/2^80) = 61 712 849 264 585`; **cap = 2 x mean = 123 425 698 529 170**.
* every per-report row — write, lookup, drop-oldest drain — is charged **at the cap**, so the bound
  holds on every run; `P(hits > cap) <= exp(-mean/3) = 2^-2.97e13` (Chernoff).
* calibration, not bound: 20 240 records against cap 512 in 6400/6400 trials (`evidence/t4-dp5.json`).
* peak storage 13 823 678 355 071 904 B across store, labels, U labels, advice, guard history, code,
  scratch.

## 7. Work ledger

Multiply step 21 921 instructions, squaring step 13 445, `c_raw = 0.937125`, `c_stress = 0.965333`
(judgement-bounded fold rows at 13/10), `W_raw = 2^125.716950`, `W_stress = 2^125.759734`, claim
**125.7598** = least four-decimal log2 above `log2 W_stress + 2e-5`. Auxiliary work is one number:
381 GT squarings + 191 dense multiplies = **412.5896** reference units per 381-bit exponentiation,
charged twice per trail start at the cap, total **2^56.499**, i.e. `2^-69.26` of `W_stress`, and
**included** in the claimed total (`2^125.759734 + 2^56.499` rounds to the same four decimals).

## 8. Evidence map

`evidence/class-retention.json` (section 3) · `evidence/trap-v16.json` (sections 4, 6, 7 — the single
aux source) · `evidence/factor-r.json` (exact factorisation, order 2, `ord_r(2)`) ·
`evidence/retention.json` (marginal rows only) · `evidence/census-joint.json`,
`evidence/pattern-census.json` (descriptive censuses) · `evidence/t4-dp5.json`, `t4-dp7.json`,
`t4-dp9.json` (solver census, 4682/4682 and 3368/3368 verified) · `evidence/toy-census.json`,
`sigma-set.json` · `part_a-out.json`, `routing-out.json`. `verify_certificate.py` re-runs every
assertion.

Failed gates are shipped, not hidden: `no_inversion_free_sigma_stable_column_table` (UNMET — the
route does not exist), `class_marginal_uniformity_refuted_at_fixture` (REFUTED — uniformity is a model
premise), `production_scale_class_stream_census` (UNMET — the sweep stops at R = 1024).

## 9. Scope

Fixture-scale measurements, no production run, no unconditional success theorem, no novelty claim, no
practical break. The claim is a classical work bound at the protected `C = 20000` with a 0.39 success
floor conditional on the disclosed hypotheses.

# argument.md — derivation, itemisation and scope for `gt-reportfilter-g6-40pct-v14`

Every number here is emitted by `ledger.py`, `gen_cert.py`, `toy_census.py` or
`sigma_set.py`; `certificate.json` recomputes each of them from the raw evidence
JSONs and `verify_certificate.py` re-derives them again independently.

The v13 verdict (`evidence` of the previous revision, submission `c92286c0`)
rejected the report-filter design on three flaws. This revision answers each one
by construction, and it records one route that turned out not to exist.

## 1. Reduction (unchanged, re-checked on BLS12-381)

`P` fixed nonidentity in G1, `z` uniform in the order-`r` subgroup of GT,
`g = e(P,Q0)` a generator by non-degeneracy, `a = log_g z`, `Q = [a]Q0`, so
`e(P,Q) = g^a = z`. `part_a.py` instantiates the ate pairing over the Fp12 tower
and checks 24 named booleans on BLS12-381 itself, including
`pairing_nondegenerate`, `pairing_bilinear_g1`, `pairing_bilinear_g2`,
`planted_ePQ_equals_z`, `planted_wrong_Q_fails`, `ord_r_p_is_twelve`,
`p_pow6_is_minus_one`, `n_classes_exact`, `key_rule`, `keys_distinct`,
`kappa_ok`, `min_ident` and `cyclo_sqr_identity`. All 24 pass.

## 2. Class structure, and the multiplier family the quotient FORCES

`ord_r(p) = 12`, `p^6 ≡ -1 (mod r)`, `N = (r-1)/12 = 2^251.27`. The twelve
Frobenius images of `y` form one class; `keys12(y)` returns twelve Fp values whose
minimum `mkey` is the sigma-invariant class key and whose index `kappa` satisfies
`kappa(sigma^j y) = kappa(y) - j`.

Merge-safety (the successor **class** must be a function of the class) requires

    m(sigma^j y) = sigma^j m(y)   for every j.

Write `y = g^v` and `m = g^{mu(v)}`. Then the requirement is
`mu(v p^j) = p^j mu(v)`, whose general solution is a **single** term
`mu(v) = c(C) p^{-kappa(v)}` — one coefficient per class, twelve columns per table
row. Two-term families (`c1 p^{-kappa} + c2 p^{-2kappa}`) are not covariant: the
second term scales by `p^2`, not `p`. This is the structural fact the next two
sections use.

## 3. Route A (a sigma-stable inversion-free six-column table) does not exist

The v14 dispatch proposed rebuilding the transition on a lower-half, inverse-free
column set. `sigma_set.py` answers it exhaustively and the answer is negative:

* all `2^12 = 4096` subsets of `Z/12` are scanned; the ones closed under the
  routing action `k -> k+1` are exactly `{}` and `Z/12`
  (`evidence/sigma-set.json:enumeration.closed`);
* `Z/12` contains the inverse pairs `{k, k+6}` because `p^6 = -1`, so the only
  closed **and** inversion-free subset is the empty one
  (`enumeration.closed_and_inversion_free == [[]]`);
* the lower-half set `{0..5}` **is** inversion-free but **not** closed. Closure
  fails by a forced involution: on it `kappa` and `kappa+6` select the same column,
  hence `m(sigma^6 y) = m(y)` while `sigma^6 m(y) = m(y)^{-1}`, so every table
  entry would satisfy `m^2 = 1`; `r` is odd, so the only such GT element is the
  identity (`lower_half.forced_relation`). That is the equivariance failure v12b's
  verdict records, here derived instead of observed.
* the action is transitive (`transitivity.is_transitive`), which is why no
  non-trivial closed subset can exist at all.

Route A is therefore not a design choice but an impossibility, and v14 does not
ship it. What it ships is Route B, a routing guard, with two proofs and one
priced consequence.

## 4. Route B — the guard, and why no absorbing cycle survives it

**The step.** `t, dp, sw = routing.fields(mkey)`; if `sw < 3` the walk squares
(three of eight classes) and nothing else happens. Otherwise the applied column is
`col = guard.lookup(t)` if row `t` occurred among the last `D = 12` multiply steps,
and `col = kappa` when it did not; then `y <- y * TBL[t][col]`,
`a <- a + U[t][col]`, and `(t, col)` is pushed. The ring is the only state besides
`(y,a,b)`, and it is 18 bits per slot.

**Why it is sigma-covariant.** The ring stores `(row, column)` pairs. For a
sigma-partner, rows are class functions (equal) and every column index is shifted
by the same `-j`, so the ring takes the same decision and the applied multiplier is
`sigma^j` of the other's. T2 measures this on successor ELEMENTS: 0 splits in 768
checks.

**Theorem 1 (no short absorbing cycle).** Let `C` be an absorbing cycle of length
`L` in the joint state `(y, ring)`, all of whose steps are multiply steps. Group
its edges by row. Each group contributes `sum_i u_t p^{-kappa_i} = u_t * c_t` with
`c_t` in `Z[p]`. For the cycle to close, `sum_t u_t c_t ≡ 0 (mod r)`; with the
`u_t` independent table entries the structural solutions are exactly the
cyclotomic ones, `c_t = 0`, i.e. each row's column multiset is a union of `{k,k+6}`
pairs and `{k,k+4,k+8}` triples (the only vanishing sums of twelfth roots).

Now the guard. In a cycle of length `L <= 2D+1 = 25`, any two positions are at
cyclic distance `<= D`, so the walk sees the earlier visit inside its ring **in
both directions around the cycle** and reuses that column. Hence within any row
group all columns are equal and `c_t = (#visits) * p^{-k}` with
`0 < #visits <= 25 < r`: never zero. So no absorbing all-multiplier cycle of
length `<= 25` exists — this is the whole v13 family (`L = 2`), the triple family
(`L = 3`), the `a,b,a,b` family (`L = 4`) and every longer variant up to 25.

**Theorem 2 (longer and mixed cycles).** A cycle of length `L >= 26` needs at
least `L/2` distinct row coincidences plus the matching root-of-unity relations,
one factor `1/R` and one factor `1/12` per coincidence: rate
`<= (1/(12R))^{L/2} <= (1/(12*2^14))^13 = 2^-228.6` per window, i.e. `2^-102.8`
over the whole `2^125.83`-step walk. A cycle containing a squaring step closes iff
`(1 - 2^s) v = B` with `s >= 1` and `1 - 2^s != 0` (checked exactly: `ord_r(2)`
exceeds the trail length by many orders — `r-1 = 2^32*3*11*19*10177*125527*859267*906349`
and `2` is a non-residue of each large factor), so it has at most one solution `v`
per cyclic pattern and its probability is `<= 1/r` per window, `<= 2^-129` over the
walk.

**Census, in the space the trap lives in.** T1 is a return-lag ladder, exhaustive
over all `r-1 = 9420` fixture elements, from an empty ring: shipped `D = 12` gives
**no return at any lag 2..26**; `D = 2` and `D = 1` return at no lag below 6 and 4;
`D = 0` (the v13 rule) returns at lag 2, 48 elements at each of lags 2, 4, 6 and 8 —
the family the verdict exhibits. T1j
builds the **exact joint map** `(element, ring)` for depth 1 — 913 741 states,
every one of them — and reports the full cycle histogram: guarded **0 two-cycles
and 0 three-cycles**, control **24 two-cycles**. The v13 fixture with its own
64-row table shows 0 two-cycles, which is exactly how the previous revision came
to report "unrefuted": the family is invisible at `R = 64` and visible at `R = 8`.

## 5. What the guard costs: merge stickiness, priced, not hidden

`D = 12` reads a per-trail ring, so two trails that meet at a common class can
disagree about the guard while their rings differ. They reconcile: after `D` steps
of identical behaviour both rings hold the same `D` rows and the same applied
columns, and the merge is then exact for ever. Per step the disagreement
probability is at most `2D/R` (each ring holds at most `D` of `R` rows), so over
the `D`-step reconciliation window a merge survives with probability at least

    1 - 2*D*D/R = 1 - 288/16384 = 0.982421875.

That factor is charged against the success target: the birthday model is solved for
`0.392/0.982421875 = 0.3990139` and the certified net success is
`0.982421875 * 0.3990139 = 0.3920000 >= 0.39`. It is a 0.008-bit price.

The birthday model itself is restated over the augmented state: the walk is a
deterministic map on `(class, ring)`, the report key and the DP predicate are
functions of the **class** alone, and the routing fields are a hash of the class
key, so the class marginal stays uniform and two random states collide in class
with probability `1/N`. `Lambda = S^2 (1-f^2) / (2 mu N)` is unchanged in form;
`mu = 11/10` remains the walk-quality allowance (a reserve, not a measurement).

T3 measures the consequence rather than asserting it: 64 pairs at a common class
with independent random rings, guard-fire disagreement 112/480 multiply steps at
the fixture's 64-row table and 9/64 merges surviving the 12-step window. That is a
**fixture-scale** number, disclosed as such: at `R = 64` the guard fires ~256x more
often than at production `R = 2^14`, and `guard_rate` in the evidence records the
measured rate at `R = 64` and `R = 512` so the `1/R` trend is visible. The
production retention used in the ledger is the analytic `1 - 2D^2/R`.

## 6. Per-step rows, one squared fraction, and the ledger

| row | derivation | instructions |
|---|---|---|
| `selector` | `SELECTOR_AB = 1 MUL + 1 ADD`; `SELECTOR_ZETA = 2 MUL`; `FOLD_ADDS = 9 ADD`; `FOLD_CMP = 11 CMP6` | 1 613 |
| `column_addr` | `index = 12*t + col`, one multiply-add | 16 |
| `guard` | `D` packed compares + `D` selects + `D` stores + index = 40, charged at 48 | 48 |
| `table_fetch` / `label_add` | 576-byte entry = 72 loads; `a += u` | 72 / 28 |
| `key_hash` / `dp_test` / `hist_slot` / `loop` | routing mix, DP word test, record slot, addressing | 72 / 8 / 32 / 32 |
| `fold_mul` / `fold_sqr` | the rows above (the guard is on the multiply channel only) | **1 921 / 1 813** |
| `cyclo_sqr` | Granger–Scott, 24 Fp muls + 122 adds | **11 632** |
| `sigma_on_state` / `renorm` / `escape` | deleted by the design, not underpriced | **0 / 0 / 0** |
| `mul_step` / `sqr_step` | `20000 + fold` / `11632 + fold` | **21 921 / 13 445** |

`c_mul = 1.096050`, `c_sqr = 0.672250`, `c_raw = 0.937125`,
`c_stress = 0.9653325` (every judgement-bounded fold row at 13/10, the counted
field rows left alone).

**One squared fraction.** `f = 3/8` — the routing share — is used in the step
price, in the step budget and in the success model. The v13 flaw was mixing the
fixture's measured `287/785` (used for `S`) with `3/8` (used for the price); there
is now a single constant `K.F_SQR = 3/8` and `gen_cert` asserts that the same value
enters `S`, `W` and `Lambda`, and that `constants.py` contains no measured-share
constant at all. The fixture share still appears in the evidence
(`t1.control_v13_rule_small_table.sqr_fraction`, `0.3656`) as a cross-check row,
and it is an input to nothing.

`S` is the least integer with `S^2 (1-f^2) >= 2 mu N ln(1/(1-0.3990139))` at
`f = 3/8`, `mu = 11/10`: `S = 2^125.827264`. Then `W_raw = 2^125.733578` and
`W_stress = 2^125.776362`; the claim is `ceil4dp(W_stress + 2e-5) = 125.7764`.

## 7. Records, restart schedule and peak storage

A record is `(mkey, a, b)` = **112 bytes**: 48-byte class key plus two full 32-byte
exponent limbs, `REPORT_WORDS = 14`. Every storage and traffic row uses 112:
`ledger.aux:report_writes = reports * 14 words * 2 instructions`, and
`ledger.storage` is generated from the same constant, so a 96-byte row cannot
survive anywhere (`gen_cert` asserts the `aux` body contains no `96` literal).

Production protocol, stated because v13 was rejected for not stating it:

* one shared store keyed by the 48-byte class key, one record per reported DP;
* each trail starts from a fresh uniform label pair (`y = g^a z^b`, two 381-bit
  scalar multiplications) and ends at its **first** distinguished point, i.e. after
  a geometric number of steps of mean `theta = 2^80`;
* the reporting trail is replaced by a fresh trail, so the number of starts equals
  the number of reports, `S/theta = 62 428 245 287 538`;
* no live trail is ever restarted, aborted or re-seeded. The guard rewrites one
  column; it is not an escape policy. This is the difference between v14 and the
  "restart without state reset" defect the board has already caught twice.

Peak storage: store `62 428 245 287 538 * 112 = 6 991 963 472 204 256 B`; labels
`2 x 32 B x 4`; advice (the `12 x 2^14` column table) `113 246 208 B`; scratch 8 GT
elements `4 608 B`; peak `6 991 963 585 455 328 B`. Memory is unscored, and the
number is stated rather than hidden.

Aux rows: table build `2^20.585`, trail starts `2^41.113`, restart bookkeeping
`2^35.540`, report writes `2^36.347`, solve `2^9.015`, verify `2^9.000`; total
`2^41.194`, which is `2^-84.58` of `W_stress`.

## 8. Success, detection, solve

`Lambda = S^2(1-f^2)/(2 mu N) = 0.5091835003`, `1-e^-Lambda = 0.3990139`, net of
the merge-stickiness retention `0.982421875` gives **0.3920000 >= 0.39**.

Detection needs no independence substitution: the DP predicate is a function of the
class, so two trails at a common class evaluate it identically and report the same
class in the same step. What detection does need is that the two trails *reach* a
common DP, and that is what section 5 counts: the merge is exact after `D` steps of
agreement, and the retention `1 - 2D^2/R` is that probability.

Solve: records `(mkey, a1, b1)`, `(mkey, a2, b2)` with `y1 = sigma^j y2` give
`a1 + b1 x = p^j (a2 + b2 x)`, hence `x = (p^j a2 - a1)/(b1 - p^j b2)` for
`j = 0..11`, each candidate **verified by exponentiation**. `j = 6` is the
involution and reads `x = -(a1+a2)/(b1+b2)`. T5 plants an independent label pair on
every one of the twelve shifts: 144/144 solved, exactly one verified candidate each,
0 wrong-shift accepts, 0 zero denominators, `j = 6` 12/12.

## 9. Census scope (disclosed, not glossed)

`toy_census.py` runs the shipped step, routing, guard and solver: 6 400 trials per
density, 4 trails, one shared store per trial, trails end at their first DP,
unconditional statistic, per-report invariant `g^a z^b = y` asserted, every recovery
verified by exponentiation, record keys re-derived at `(trial, trail, step)`.

The fixture success rates are what they are and they are **not** a production
floor: `0.7316` at `2^-5` and `0.5262` at `2^-7` are *above* 0.39, `0.2778` at `2^-9` is
below it, and the sequence falls with depth because the fixture is cycle-dominated.
They are higher than the previous revision's because the guard removes the 2-cycles
that were consuming the fixture's walks — a fixture-level effect, and no production
recovery rate is inferred from it. The measured guard fire rate at the fixture is
`0.479 / 0.350 / 0.346` per multiply step against the production `12/16384`; the
64-row toy table is the reason, and it is disclosed. Two ratios are quoted and they are different quantities:
`S/sqrt(N) = 2.3` at the fixture against `0.282` at production is the walk budget
in units of `sqrt(N)` (the coverage-law variable), while `w/sqrt(N) = 0.282` with
`w = S/TRAILS` is the same ratio for one trail of four; the coverage law uses the
first, the census table reports the second. v13's note wrote `1.13` and `0.282` in
one sentence without saying which was which.

## 10. Failed gates (stated in the certificate, not buried)

`certificate.json:failed_gate` carries two rows. `no_inversion_free_sigma_stable_column_table`
is FAIL: Route A does not exist (section 3), so the guard is the shipped mechanism
by necessity. `production_scale_guard_census` is NOT_MET: no guard census at
`R = 2^14` is possible in this package; the production retention is analytic and
the fixture fire rates at `R = 64` and `R = 512` are disclosed as fixture-scale.

## 11. Legality

Claim 125.7764 against the live lowest head 125.84 (polled at dispatch, `review`).
Nothing in this package depends on the identity or validity of any other
submission; the comparison is context only.

# v21 — intake path contract fixed, canonical g^a z^b census relabelled by derivation

**Claim 125.7704** unchanged. **Success 0.3904** point estimate, CP-95 lower 0.3686 shipped
as a FAILED gate. **Lead 0.0496 bits** over the live board minimum 125.82 (submission
`aba345a`, valid) and 0.0641 over the 125.8345 legality line. Nothing in the mathematics
moved in v21; this section records what did.

## 0. Why v20 never reached a reviewer: the selector, not the package

Five v20 tickets died at the intake step with `llm-verifier: unsafe_or_missing_directory`
(3f9af544, 39b4b5f4, 2eba1e80, dd4bc86f, plus the 297 KB archive variant 1fdf90b7, which
was a size kill). The cause is mechanical and now proven against the last tree that passed
intake: `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt` still contained
`submissions/bls12-381-fixed-g1-invert-g2/cls-sel3-counted-v17`, and v17 had been archived
to `workspaces/archived-g1inv-2026-10-24/`, i.e. **outside** the archive root. The intake
verifier resolves the selector first (`files.py:directory_fd`); a selector naming a
directory that is not in the expanded archive raises `unsafe_or_missing_directory` before
any package file is opened. That is why 26 files / 198,934 B / a 7,860-char claim passed
every size, count and text rule and still died in ~70 s with no evaluator verdict.

The contract is now enforced mechanically before every submit (`preflight_v21.py`, run from
the benchmark root; a local tool, not part of this package):
the selector content equals the staged package path **as an exact string**, that directory
exists inside the expanded archive, and every path named in `submission.json`
(`claims[0].evidence`) and in this note resolves from the package root inside the same
expanded tree.

## 0b. The census is on the canonical g^a z^b label (v19 open defect)

The v19 verdict required the repeat census to be stated on canonical `g^a z^b` group form
rather than on an internal state. It is now derived, not re-measured. With `g = e(P,Q0)`,
`z = g^x` and the shipped label pair `(a,b)` satisfying `y = g^a z^b`, define

    t = a + x*b  (mod r)          and assert  g^t == y  at every step.

`t` is the discrete log of `y` base `g`, hence a bijection `Z/r -> <g>`: a `(t, guard)`
collision and an `(element, guard)` collision are the same event, and the shipped 48-byte
class key is coarser (twelve Frobenius partners share it), so every element repeat is a key
repeat but not conversely. `census_gab_v21.py` measures the identity rather than asserting
it: 25,600 states, 0 invariant failures, 3,753 distinct elements = 3,753 distinct canonical
labels, 21,847 repeat pairs of which 0 are zero-`b`-increment and 21,847 nonzero, 0
key-repeat violations (`evidence/census-labelled-v21.json`). The 1,000,000-step-pair
production census is therefore already a canonical-label census, so the v20 trap numbers
(1,744 nonzero-increment, 0 zero-increment, all charged) stand and no numerator moves.

**Claim 125.7704** (log2 dense-Fp12-multiply equivalents, C = 20000, all processors,
auxiliary work inside the claim: `claim.covers` = W_stress + aux). **Success 0.3904**
point estimate with every unfinished trail scored as failure; the CP-95 lower is
0.3686 and is shipped as a FAILED gate, not as a claim.

## 1. A repeat traps whether or not its increment is zero (verdict defect 1)

v18 charged the zero-increment fraction only. That was the wrong event, and the census
was re-run to prove it rather than argue it: **1744
joint repeats in 1000000 step pairs, 0 of them
zero-increment and 1744 nonzero** (lags
20/40/60/80/100/120/140/160; period-20 family,
164/256 walks repeat inside
400 steps). The shipped attack stores DP records and guard entries and nothing
else, so a return to the same (element, guard) state is a cycle whose twenty-odd elements
contain no DP: it never reports and the trail is pending at every cutoff. The increment only
decides whether a solve *would* exist had both records been kept — and the attack keeps
neither. The ladder exhibits both branches: the zero-increment pair yields
11 candidates and **0 verify**; the nonzero-increment
pair yields 12 candidates and **1 verifies**.

Charged with ONE named inequality and no independence hypothesis (v20 fixes the v19 defect
that read a per-worker Markov bound as "at most one event per worker"):

  X_w = trap events on worker w's own walk; expectation at most Lambda_rep = 0.0539 per worker by per-walk Markov
  on the birthday count (a fixed map, so no coin dependence enters); X = sum_w X_w and
  expectation at most 4 x 0.0539 = 0.2156 by linearity ALONE - no independence across workers, because
  none is needed.  Markov on the TOTAL: the chance of exceeding C events is at most E[X]/C, with C = ceil(4 Lambda / 1e-9) =
  215600000 events gives a tail no larger than 1e-9.  Every event is priced at a DP lookup (0.0210)
  + a detector probe and 112-byte compare (0.0210) + a full restart, two GT exponentiations
  (825.18) + the wasted cycle steps (70.1); the total 1.9e11 reference units sits in aux.trap_charge
  inside the claim. The charge is 1.9e11 and the claim moves by 2^-88 bits: the tail is
  quantified, not asserted, and it is still free at four decimals.

Recovery is control flow, not prose (v20 deletes the v19 "64 wasted steps" sentence, which no
code implemented). The shipped trail ends at its FIRST distinguished point OR at the hard cap
128 x theta = 128 theta steps, whichever comes first; a capped trail is ABANDONED, its labels are
discarded and a fresh trail starts (two GT exponentiations). The cap is a geometric tail:
P(a trail exceeds 128 theta) <= exp(-128) = exp(-128), so the abandon count is charged at one
event (about 1950 reference units) and the wasted 128 theta steps at 128 theta x the step price - both inside the trap
row above. The detector price is the probe the store already pays; no separate monitor exists.

## 2. Success carries the retention product, and its CP-95 lower is below the floor
(verdict defect 2)

The measured retention is `retention_eff` = mu x p_detect x p_usable = 1.1 x 0.895833 x 1.0 =
0.985417. Applied to the birthday **rate**: 0.3920 (the 0.392 the verdict
computes). Applied to the **probability**: 0.3907 = 0.3906 >= 0.3906, which is the
more conservative of the two and is what ships. Net of the trap debit and the density stopping
loss: **0.3904 >= 0.3904**. The CP-95 lower is
**0.3686 <= 0.3690 and <= 0.39**, i.e. **below the floor**: the point estimate clears
0.39, the lower bound does not, and no amount of prose changes that. Only 0.9303 of
the rate is certified at 95 per cent; the residual is a disclosed risk, not a claim.

**Production-lag disclosure (v20).** The retention statistic is measured at class-agreement
lags k = 1..56; the production report horizon is theta = 2**80 steps steps, and the shipped k-curve
falls to 22/240 by lag 32. The guarded step is history-dependent, so coalescence at the
production horizon is NOT established by the fixture curve and is not claimed: the point
estimate 0.3904 clears 0.39 and the CP-95 lower 0.3686 does not. This is shipped as a named
failed gate (`production_lag_retention`), not as a proof.

## 3. No variance bound is used (v20 withdraws the Var <= E row)

The v19 stopping loss used Var(R) <= E[R), which the verdict rejected as underived, and a DP
density never makes the report count exactly ceil(S/theta). v20 removes the step instead of
defending it: the report count enters ONLY as a COST CEILING. The cap is the measured worst
window factor 2.080 x mean = the measured worst-window cap, every per-report row is priced at it, and a run that
exceeds it pays the charged overflow path (drop-oldest + one lookup + one write per dropped
record, `ledger.aux:report_cap_drain`). `success.net` contains no variance term at all
(the v19 Chebyshev row 1.379e-14 is WITHDRAWN, not re-estimated), so no claim rests on
Var(R) <= E[R].

## 4. No rounding up: every >= and <= literal is audited with its direction (defect 4)

Every inequality below is a table row, re-read from this note by `audit_v19.py` and applied to
the certificate value with a direction rule: a floor literal must satisfy value >= literal, a
ceil literal value <= literal. v18 printed "at least 0.3963" over a value of 0.3962569 - a
round-UP - and the same rule rejects it (control in `evidence/inequality-audit.json`).

| row | op | literal |
|---|---|---|
| `success.net` | >= | 0.3904 |
| `success.net` | <= | 0.3910 |
| `success.model` | >= | 0.3964 |
| `success.product` | >= | 0.3906 |
| `success.product_cp95_lower` | >= | 0.3688 |
| `success.net_cp95_lower` | <= | 0.3690 |
| `success.net_cp95_lower` | <= | 0.39 |
| `claim.value` | <= | 125.7704 |
| `claim.value` | <= | 125.8345 |
| `claim.value` | <= | 125.84 |
| `trap.repeat_probability` | <= | 0.0540 |
| `density.density_factor` | <= | 2.10 |
| `ladder.passed` | >= | 24 |

The figures in this note are floored, not rounded (0.3904 for a value of
0.390478). 10 `2^x` labels are recomputed from their decimals, with the
planted-flip control rejecting a negated label.

## 5. The toy ladder, twenty-four full trials (defect 5)

`measure_v19.ladder` at the fixture the package already ships (r = 9421, 61-bit p,
ord_r(p) = 12, N = 785, table 64 x 12, DP width 2^-5, guard depth 12):
plant x, walk the shipped guarded step, end each trail at its first DP, store the DP records,
solve all twelve shifts on a class-key collision, verify g^x == z. **24/24 pass**, steps
30..197; controls: wrong secret rejected, zero-increment trap fires,
nonzero-increment pair solves. Rows are `evidence/ladder-v19.json`.

## What is kept, and what is not claimed

The fixed-P reduction, the solve identity x = (p^j a1 - a2)/(b2 - p^j b1) for j = 0..11 with
every candidate verified, the exact r-1 factorisation, the guard at D = 12, and the measured
rows (class match 0.0833 within 256 steps, class-stream survival 0.9667 at k=1 to 0.7833 at
k=12, DP hits 1.9e-3 at lags 20/40/60, the density census, aux 412.5896 reference units per
381-bit GT exponentiation). Success is conditional on the three disclosed hypotheses
(random-map walk quality inside the 10 per cent allowance, pseudorandom routing and table
fields, ideal coins), and the claim covers W_stress 125.770322 plus aux
2^56.566 with a lead of 0.0496 bits over the live board minimum 125.82 and
0.0641 over the 125.8345 legality line (the certificate keeps the conservative 125.84 row). Four gates fail and ship as failures:
the success CP-95 lower is below 0.39, the inversion-free sigma-stable table does not exist,
the class marginal is not uniform at fixture scale (credit not taken), and the class-stream
census stops at R = 1024.

**Reproduce.** `python3 measure_v19.py` (repeat classes, density, ladder), `python3
gen_cert_v19.py` (certificate + note + audits), `python3 audit_v19.py` (independent re-read),
`python3 census_gab_v21.py` (canonical g^a z^b relabel witness).

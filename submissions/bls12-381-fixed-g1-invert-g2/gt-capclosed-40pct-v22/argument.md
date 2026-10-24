# argument.md — derivation and itemisation for `gt-trapcharged-40pct-v19`

v19 answers verdict `5e6e07c4` on v18. The algorithm, the reduction, the collision identity,
the `r-1` factorisation, the guard and every measured census row are unchanged: the defects
were in the accounting, and v19 charges what v18 left out rather than arguing it away.

## 1. The algorithm (unchanged)

Walks on `y = g^a z^b` in the order-`r` subgroup of GT, `N = (r-1)/12 = 2^251.27` Frobenius
classes. A step takes the class key `min{keys12(y)}` and rotation index `kappa`, routes the
class to a squaring (3/8) or to a twelve-column multiplier table (5/8), and on the multiply
channel reuses the column recorded when the same row was last used within the last `D = 12`
multiply steps. A report is a distinguished class. Two stored records with equal class key
solve

    x = (p^j a1 - a2) / (b2 - p^j b1),   j = 0..11,   each candidate VERIFIED against z.

## 2. Success: retention is a factor, and its CP-95 lower is below the floor

`Lambda = S^2 (1-f^2) / (2 mu N) = 0.5049381183`, `1 - exp(-Lambda) = 0.3964571`, at
`S = 2^125.821225`, `mu = 11/10`, `f = 3/8`. The step budget is the least `S` with
`S^2 (1-f^2) retention >= 2 mu N ln(1/(1-0.392))`, so `S` already carries the retention
correction — and v18 then reported the success *without* it. v19 composes it:

    retention_eff = mu * p_detect * p_usable = 1.1 * 0.895833 * 1.0 = 0.9854167
    rate-scaled   = 1 - exp(-Lambda * retention_eff)      = 0.3919964
    product       = (1 - exp(-Lambda)) * retention_eff    = 0.3906754   <- the conservative one
    net           = product * (1 - tau) * (1 - stop_loss) = 0.3904781

`success.net` is the claimed row. Its CP-95 lower, with `p_detect >= 0.8576` and
`p_usable >= 0.9862`, is 0.3688077 before the trap and stopping debits and 0.3686215 after:
**below 0.39**. That ships as a failed gate (`success_cp95_lower_below_floor`, UNMET). The
floor is met by the point estimate and is not certified at 95 per cent; nothing here says
otherwise.

## 3. Traps: every repeat, nonzero increment included, is charged

The v16 cardinality count `S^2/(2(r-1)(12R)^D)` stays **withdrawn** (164/256 fixture walks
repeat inside 400 steps, which a 2^-215 event cannot do). v18 replaced it with the
zero-increment fraction, and that was still the wrong event. The shipped attack stores DP
records and guard entries and nothing else: no self-collision detector, no state history. A
return to the same `(element, guard)` state is a cycle, and the cycle's twenty-odd elements
contain a DP with probability at most `20 * 2^-80`, so the trail never reports and is
**pending at every cutoff** — whatever its label increment. The increment decides only
whether a solve would exist *if both records were kept*, and neither is.

Re-run of the census on the shipped step (1e6 step pairs, 5000 walks x 200 steps), every
repeat tagged by its increment:

    1744 joint-state repeats, 0 zero-increment, 1744 nonzero-increment
    lags 20/40/60/80/100/120/140/160 (one period-20 family)
    164/256 walks repeat inside 400 steps; 13546 events there, 0 zero-increment

Both classes are charged. The count is the birthday count over the whole walk,

    Lambda_rep = S^2 / (2 (r-1)) = 0.0539   =>   P(at least one trap) <= 0.0539   (Markov)

so at most one event per worker (4) is charged, each at measured prices: DP lookup 0.0210,
detector probe plus 112-byte record compare 0.0210, a full restart of two GT exponentiations
1526.19, and 64 wasted cycle steps 76.80 — total 3581.5 reference units, `2^11.806`, carried
as `aux.trap_charge` inside the claim. It is small because `Lambda_rep` is small; it is no
longer free, and it is no longer conditional on the increment being zero.

The ladder exhibits both branches (section 6): the zero-increment pair yields 11 candidates
and **0 verify**; the nonzero-increment pair yields 12 candidates and **exactly 1 verifies**.

## 4. Reports: the trigger decides, and the cap carries the measured density

A report is the trail-terminating distinguished point at `theta = 2^80`: the trail ends and a
fresh trail starts, so `reports = trail starts = ceil(S/theta) = 2^45.821`.

A density is not an interval. The report count is random, so v19 measures its dispersion
(19200 windows of 128 steps at the 2^-5 fixture: mean 3.365 hits, variance 0.419,
over-dispersion 0.124, worst window 7 = 2.080 x the mean) and charges the cap at the measured
worst window:

    cap = 2.080 x mean = 2^46.878,   every report row priced at it
    GT starts 2^56.566, write 2^37.397, lookup 2^41.304, drain 2^41.304

The stopping rule stops at the cap, so its cost is a success row, not a footnote: Chebyshev
with variance `<= mean` gives `P(R > cap) <= mean/(cap-mean)^2 = 2^-46.0`, subtracted in
`success.net`. The v17 declaration `2^121.56` stays withdrawn: it transferred the fixture's
`2^-5` report rate onto production `S`, needs `19.1` steps per report against a `2^80` trigger,
and was priced against the stale `S = 2^125.810636`.

## 5. Rounding: every `>=` and `<=` carries its direction

`audit_v22.py` re-reads `public-note.md` and `certificate.json` and applies one rule to every
inequality: a floor literal must satisfy `value >= literal`, a ceil literal `value <= literal`.
The v18 sentence claiming a lower bound of 0.3963 over a value of 0.3962569 is rejected by
that rule, because 0.3963 is a round-up; the control asserts the rejection. Thirteen table
rows ship in `evidence/inequality-audit.json`; ten `2^x` labels are recomputed from their
decimals with a planted-flip control.

## 6. The toy ladder: twenty-four full trials, 24/24

`measure_v19.ladder` at the shipped fixture (`r = 9421`, 61-bit `p`, `ord_r(p) = 12`,
`N = 785`, table 64 x 12, DP width `2^-5`, guard depth 12), fixed seed, one pass, no retries:
plant `x`, walk the guarded step, end each trail at its first DP, store the DP records, solve
all twelve shifts on a class-key collision, verify `g^x == z`. **24/24 pass** (steps 30..197,
the invariant `y = g^a z^b` re-checked at every DP). Controls: the planted secret plus one is
rejected; the zero-increment pair is a trap; the nonzero-increment pair solves. Rows ship in
`evidence/ladder-v19.json`.

Scope is stated, not glossed: the toy is supercritical (`S/sqrt(N) ~ 2.3` at the fixture
against 0.28 at production) and its DP width is `2^-5` against `2^-80`, so the ladder measures
the merge -> solve -> verify mechanics, not a production recovery rate.

## 7. Cost, storage, and what is not claimed

Per-step counted rows are unchanged: 21921 multiply, 13445 squaring, `c_raw = 0.937125`,
`c_stress = 0.965333`, `W_raw = 2^125.727538`, `W_stress = 2^125.770322`. Records are 112
bytes in one shared store; peak storage 13925508664531904 B including the 4 MiB repeat-detector
table. Auxiliary work uses ONE GT-exponentiation total: 381 cyclotomic squarings at 11632
instructions plus 191 dense multiplies = 412.5896 reference units, charged twice per trail
start at the cap, `aux.total = 2^56.566` with the trap rows inside it. The claim is the least
four-decimal log2 above `log2(W_stress + aux.total) + 2e-5`, i.e. **125.7704**; it leads the
live frontier 125.84 by 0.0696 bits and the 125.8345 legality line by 0.0641.

Not claimed: no production run, no unconditional success theorem, no novelty, no practical
break. Success is conditional on the three disclosed hypotheses (random-map walk quality
inside the 10 per cent allowance, pseudorandom routing and table fields, ideal coins). Four
gates fail and ship as failures: the success CP-95 lower is below 0.39, no inversion-free
sigma-stable column table exists, the class marginal is not uniform at fixture scale (credit
not taken), and the class-stream census stops at `R = 1024`.

## 8. Evidence map

`evidence/retention-v19.json` (repeat, retention, usability, reports census) ·
`evidence/repeat-classes-v19.json` (increment-tagged repeats, both classes) ·
`evidence/density-v19.json` (windows, variance, worst window) ·
`evidence/ladder-v19.json` (24/24 rows and controls) ·
`evidence/inequality-audit.json` (section 5 rows) · `evidence/audit-v19.json` (independent
re-read) · `evidence/label-audit.json` (labels, flip control) ·
`evidence/reports-trigger.json` (section 4) · `evidence/factor-r.json` ·
`evidence/sigma-set.json` · `evidence/census-labelled-v21.json`. `measure_v22.py` regenerates the three
v19 measurement files; `gen_cert_v22.py` regenerates `certificate.json` and `public-note.md`;
`audit_v22.py` re-reads them independently; `census_gab_v22.py` regenerates the canonical
`g^a z^b` relabel witness.

## 9. Canonical g^a z^b relabel (v19 defect, closed in v21)

Every trail carries an exact label pair `(a,b)` with `y = g^a z^b`, `g = e(P,Q0)`,
`z = g^x`. Since `z = g^x`, the canonical exponent of the state is `t = a + x b mod r`,
and `g^t = y` holds identically; `t` is the discrete log of `y` base `g`, a bijection onto
`Z/r`. Two consequences are used, both shipped as measurements in
`evidence/census-labelled-v21.json` and regenerated by `census_gab_v22.py`:

1. `(t, guard)` collision <=> `(element, guard)` collision. The census that the trap charge
   rests on is therefore a canonical `g^a z^b` census; the v20 counts (1,744
   nonzero-increment, 0 zero-increment over 1e6 step pairs) need no re-measurement and no
   numerator moves. Fixture witness: 25,600 states, 0 invariant failures, 3,753 distinct
   elements = 3,753 distinct labels.
2. The shipped 48-byte class key is coarser than element identity (twelve Frobenius
   partners share it), so element repeat => key repeat, never the converse; the solve
   enumerates the twelve shifts and verifies, and the charge counts the finer event, which
   is the conservative side. Witness: 0 class-key violations on 21,847 repeats.

The `delta-b = 0` split is the solvability split: with `b1 = b2` the collision equation
`a1 + b x = p^j(a2 + b x)` carries no `x`, so such a repeat is unsolvable (the v18 trap
class). The fixture run again observes 0 zero-`b` repeats against 21,847 nonzero, matching
the production split, and all repeats are charged.

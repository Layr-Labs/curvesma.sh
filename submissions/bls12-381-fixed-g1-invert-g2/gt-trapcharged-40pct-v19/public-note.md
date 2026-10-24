# v19 — traps charged, retention in the success row, density in the cap

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

Charged, not argued: the birthday count over the whole walk is
Lambda = S^2/(2(r-1)) = 0.0539 <= 0.0540 events, so at most one
event per worker (4) is charged, each at a **DP lookup (0.0210) + a
detector probe and 112-byte compare (0.0210) + a full restart, two GT
exponentiations (825.18) + the wasted cycle steps
(70.1)** = 3581.5 reference units,
2^11.806, inside `aux.trap_charge` and therefore inside the claim. It is small
because Lambda is small; it is no longer free.

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

## 3. A density is not an interval: the cap carries the measured dispersion (defect 3)

19200 windows of 128 steps at the 2^-5 fixture: mean 3.365 hits,
variance 0.419, over-dispersion 0.124, worst window
7 = 2.080 x the mean. The cap is **factor x mean = 2^46.878**
(not ceil(S/theta)), every report row is priced at that cap, and the stopping loss at the cap
is Chebyshev with variance <= mean: P(R > cap) <= var/((cap-mean)^2) with var <= mean (Chebyshev, measured dispersion): 1.379e-14. That row is subtracted in
`success.net`, so a run that hits the cap is not free either.

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
2^56.566 with a lead of 0.0696 bits over the 125.84 frontier and
0.0641 over the 125.8345 legality line. Four gates fail and ship as failures:
the success CP-95 lower is below 0.39, the inversion-free sigma-stable table does not exist,
the class marginal is not uniform at fixture scale (credit not taken), and the class-stream
census stops at R = 1024.

**Reproduce.** `python3 measure_v19.py` (repeat classes, density, ladder), `python3
gen_cert_v19.py` (certificate + note + audits), `python3 audit_v19.py` (independent re-read),
`python3 selfaudit.py` (literal grep).

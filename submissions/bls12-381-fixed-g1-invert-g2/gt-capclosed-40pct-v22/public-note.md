# v22 — the abandon cap is in the code and metered, the trap expectation has a derivation, retention is re-measured through the guard

**Claim 125.7711** (log2 dense-Fp12-multiply equivalents, C = 20000, all processors,
auxiliary work inside the claim: `claim.covers` = W_stress + aux). **Success 0.3523**
point estimate with every unfinished trail scored as failure; the CP-95 lower is
0.3337. BOTH the point estimate and the CP-95 lower are BELOW the 0.39 floor: this is a measured failure of the success gate, shipped as such, not a claim.

## 1. The abandon cap is a hard 128 x theta = 2^87 steps, implemented and metered (verdict defect 2a)

The v21 package priced 64 wasted steps in code while this note described abandonment after
128 x theta steps; both halves were wrong, because the code did not implement the cap it priced
and the note's exp(-128) tail was the wrong model for a trail trapped in a cycle. v22 fixes the
code and DELETES the tail.

The shipped trail protocol is now: a trail ends at its FIRST distinguished point OR at a hard cap
of 128 x theta = 2^87 steps, whichever comes first. A capped trail is ABANDONED: its labels are
discarded and a fresh trail starts (two GT exponentiations). `measure_v22_ext.cap_abandon`
implements that counter on the shipped step, meters the abandon path, and measures the delay:

- trails that report exit at their first DP, mean 28.07 steps at the 2^-5 fixture (the cap is never reached);
- a trail that returns to a visited (class, ring) state is on a cycle: 2 of 256 fixture trails trap, cycle length mean 11.5 and max 14;
- every trapped trail runs to the cap and is abandoned there: measured delay mean 4096.0 steps (= the cap exactly), variance 0.0;
- control (force admission): the SAME state and the SAME ring continue WITHOUT the cap for 2 x cap steps and still produce no DP (2 of 2). Without the cap those trails run forever, which is exactly what the cap bounds.

So there is no tail to bound and none is used: with a hard cap P(a trail runs past the cap) = 0
by construction, and the waste row is deterministic -- **cap x rate x step price**, one row. The
cap is priced in the claim: 1.7e+26 reference units of
waste per abandon. Charged total: Lambda = S^2/(2(r-1)) = 0.0539;
Markov on the TOTAL gives P(X > C) <= E[X]/C at C = ceil(4 Lambda / 1e-9) =
215,440,264 events; each costs a DP lookup (0.0210)
+ a detector probe and 112-byte compare (0.0210) + a restart, two GT
exponentiations (825.18) + the cap waste =
**2^114.815 reference units**, inside `aux.trap_charge` and therefore inside
the claim. `total_integer_rows_1_over_C_REF` re-sums the four integer rows and asserts equality.

## 2. The expectation bound has a derivation: monotonicity, plus Markov on a measured E (verdict defect 2b)

The verdict is right that Markov needs a justified expectation bound and that a fixed,
history-dependent walk does not automatically have uniform pair-collision probabilities. v22
stops applying the iid birthday formula by assumption and ships the argument.

**(a) Monotonicity (primary).** The guard adds exactly one coordinate to the state:
y -> (y, ring). A joint repeat therefore REQUIRES an element repeat -- pointwise, for every walk,
with no independence hypothesis, the joint collision count is at most the element collision
count, and it is the element stream that the iid birthday count describes. Measured on
159200 step pairs (`evidence/bound-v22.json`): joint repeats
3816 <= guarded element repeats 55133 <=
unguarded element repeats 130634. The guard makes collisions rarer,
so S^2/(2(r-1)) = 0.0539 is an UPPER bound for the guarded walk.

**(b) Markov on a measured expectation (cross-check, hypothesis-free).** The census measures
E[X] directly: 0.3488 joint repeats per walk over
5000 walks, with no model. Markov P(X >= t) <= E[X]/t needs nothing else.

Withdrawn: the exp(-128) tail, and any reading of the per-worker bound as "at most one event per
worker". The charge is a total, C = ceil(4 Lambda / 1e-9) with a 1e-9 tail, and linearity across
workers is used -- not independence.

## 3. Retention is re-measured through the guarded step, and the report count is census-direct (verdict defects 1 and 3)

The 0.895833 that v19 shipped was an EXTRACT of a v18-era JSON (its own `source` field says so),
i.e. not measured on the object the success row multiplies. v22 re-runs the same
collision-conditioned arm inside the census build through the SHIPPED guarded step, ring and
prev-row active, with a guard-off control beside it (`evidence/retention-v22.json`,
480 planted pairs, known-answer gate 480/480):

- guarded p_detect = **0.808333**, CP-95 lower **0.776352**;
- guard-off control p_detect = 1.000000 (with no ring the two trails never separate);
- guarded k-curve: 1:0.7979, 2:0.6583, 4:0.4458, 8:0.2417, 16:0.1167, 32:0.1062, 56:0.1042.

The product uses the collision-instant guarded value: 0.808333 x mu 1.1 x p_usable 1.0 =
0.889167. The production horizon theta = 2^80 is NOT established by a 56-step fixture
curve and is shipped as a NAMED FAILED GATE, not as a claim.

Reports: the count is CENSUS-DIRECT. 12733 DP-triggered reports over
240 guarded walks = 53.054 per walk against the model 1/theta (2^-5) and
ceil(S/theta) = 64 per walk at the same width (ratio
0.829); ceil(S/theta) is kept ONLY as a CHECK row, never as the count.

## 3b. A density is not an interval: the cap carries the measured dispersion

19200 windows of 128 steps at the 2^-5 fixture: mean
3.365 hits, variance 0.419, worst window
7 = 2.080 x the mean, so every per-report row is priced at the cap
2^46.878. No variance bound is used in success (the v19 Chebyshev row stays withdrawn).

## 4. No rounding up: every >= and <= literal is audited with its direction (defect 4)

Every inequality below is a table row, re-read from this note by `audit_v22.py` and applied to
the certificate value with a direction rule: a floor literal must satisfy value >= literal, a
ceil literal value <= literal. v18 printed "at least 0.3963" over a value of 0.3962569 - a
round-UP - and the same rule rejects it (control in `evidence/inequality-audit.json`).

| row | op | literal |
|---|---|---|
| `success.net` | >= | 0.3523 |
| `success.net` | <= | 0.3524 |
| `success.model` | >= | 0.3964 |
| `success.product` | >= | 0.3525 |
| `success.product_cp95_lower` | >= | 0.3338 |
| `success.net_cp95_lower` | <= | 0.39 |
| `success.net_cp95_lower` | <= | 0.3338 |
| `claim.value` | <= | 125.82 |
| `claim.value` | <= | 125.8345 |
| `claim.value` | <= | 125.84 |
| `claim.value` | >= | 125.7711 |
| `trap.repeat_probability` | <= | 0.0540 |
| `density.density_factor` | <= | 2.10 |
| `ladder.passed` | >= | 24 |

The figures in this note are floored, not rounded (0.3523 for a value of
0.352338). 10 `2^x` labels are recomputed from their decimals, with the
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
2^114.815 with a lead of 0.0689 bits over the 125.84 frontier and
0.0634 over the 125.8345 legality line. Four gates fail and ship as failures:
the success CP-95 lower is below 0.39, the inversion-free sigma-stable table does not exist,
the class marginal is not uniform at fixture scale (credit not taken), and the class-stream
census stops at R = 1024.

**Reproduce.** `python3 measure_v22.py` (repeat classes, density, ladder), `python3
gen_cert_v22.py` (certificate + note + audits), `python3 audit_v22.py` (independent re-read),
`python3 selfaudit.py` (literal grep).

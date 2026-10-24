# v17 — measured trap bound, split detection/usability, cumulative reports

**Claim 125.7704** (was 125.7598). Success 0.3965 (floor 0.39). Three premise repairs of verdict `41783ecf`:

1. **Trap bound withdrawn and re-derived from measurement.** The v16 `2^-215.26` cardinality bound is
   contradicted by our own fixture (164/256 walks repeat inside 400 steps). A 1e6-step census on the
   shipped guarded step measures the per-pair joint-repeat rate: zero at lags 1..19 (guard removes
   cycles <= 25), then 1.9e-3 at lags 20/40/60 — one period-20 family entering at lag 228. A repeat
   traps only when its label increment is (0,0); measured 0/320 joint and 0/5173 element repeats.
   Bound `tau <= (S^2/2)(1/(r-1))P(zero increment)` = **2^-26.7** at the CP-95 upper (2^-10.95
   worst case), waste <= 2^-12 of the walk, inside the already-charged mu = 11/10 allowance.
   Reconciliation row shipped: the fixture-scale value of the same formula is >= the observed rate.
2. **Detection and usability split.** p_detect = 0.8958 (215/240, CP-95 0.8576),
   p_usable|detected = 1.0000 (240/240 witness recoveries, known-answer gated). The product is
   carried by charging the measured excess over the mu allowance (1.0148) and regenerating S.
3. **Cumulative reports.** Measured reports/step, cumulative bound = rate x workers x slack, cost
   charged at the cumulative number, and one inequality with its hypothesis: per-walk Markov
   (`P(R > slack E[R]) <= 1/slack`, mean only) with a union over workers. The Chernoff row is withdrawn.

# public-note.md — `gt-classretention-40pct-v16` (fixed-G1 into G2)

**Claim: 125.7704** (log2 dense-Fp12-multiply equivalents at the protected C = 20000, all processors,
auxiliary work included) for **success 0.392** on the track's fixed P and uniform z.

## What changed from v15, and why

The v15 verdict (`408a011b`) rejected four things. Each is repaired here by changing the *event*, not by
re-wording the bound.

**1. Retention (was measured zero inside the predicate).** v15 measured the persistence of the joint
state `(element, ring tuple)` and charged `0.9809578`. The predicate accepted a meeting anywhere in
`1..2D` (its 861 events all sat at lags 13..24), and — the deeper error — the joint state is not what
the solve needs. The solve consumes a **class match**: two records with equal class key, solved by
`x = (p^j a1 - a2)/(b2 - p^j b1)` for `j = 0..11` with every candidate verified. The class key is
`min{keys12(y)}`, a function of the element **alone**; the ring never enters it, so a match is usable
at the report instant and no later rewrite can undo it. **The retention factor is therefore 1 and is
not carried.** The same 0.392 target is kept, so the debit's removal **lowers** the claim
(125.7778 -> 125.7704) instead of raising success.

What the ring can destroy is path-sharing, measured over `D in {12,24,48,128}` at R = 64 and 1024:
class streams stay identical `k=1` after a merge at `0.9667` and `k=12` at `0.7833` for the shipped
`D = 12, R = 1024`, falling to `0.3667`, `0.0417`, `0.0083` at `D = 24, 48, 128`. Shadowing lasts tens
of steps against a `2^80`-step report interval (`< 2^-74` of the walk). Independent starts match
classes at `0.0833` of pairs within 256 steps (20/240) — the birthday event the model counts.

**2. Trap (pattern multiplicities).** The v15 pattern union bound is **withdrawn, not patched**:
enumerating all words of length 4 and 6 and grouping by row-multiset class shows the omitted
multiplicity is `8.9x` at L=4 and `74.4x` at L=6 (`465x` the `C(L,L/2)/R^{L/2}` charged), and
multiplicities restore the entire word space, so no pattern count can bound a cycle. The trap is
counted on **joint states** instead: the step is deterministic on `(element, ring tuple)`, so a cycle
is a state repeat and

    tau <= S^2 / (2 (r-1) (12R)^D) = 2^-215.26   at R = 2^14, D = 12   (M = 2^465.9)

(monotone: `2^-426.3` at D=24, `2^-848.3` at D=48). That is **164 bits stronger** than v15's 2^-50.7,
multiplies no `1/(12R)` factor, and has no enumeration to omit multiplicities from. The `sum p^2`
inference is **dropped** — no Poisson lower bound and no collision-of-collisions step is used.

**3. Reports and the cap.** Mean DP hits `61 712 849 264 585`; **cap = 2x mean = 123 425 698 529 170**;
every per-report row charged **at the cap**, so the bound holds on every run, with
`P(hits > cap) <= exp(-mean/3) = 2^-2.97e13`. The census row (20 240 records under a cap of 512,
6400/6400 trials) is shipped as a **calibration**, not as the bound.

**4. Aux, one number.** A 381-bit GT exponentiation is 381 cyclotomic squarings (11 632 instructions
= 0.58160 ref units) plus **191** dense multiplies = **412.5896**. v15 printed 411.59 in one place and
412.5896 in another; 411.59 is the same expression with 190 (the mean Hamming weight) and is
superseded. `evidence/trap-v16.json:aux` is the single source, `ledger.gt_exp_reference_units()`
asserts equality, and the aux total `2^56.499` is **included** in the claimed total.

## What is kept

The fixed-P reduction and the collision-solving identity; the exact `r-1` factorisation (255 bits,
twelve primes, re-multiplied) and the `p^6 = r-1` order-2 row; the chi-square machinery and the
non-uniform class marginal (`chi2 = 238889.4` on 526 df) with the credit **not** taken; the 448
`D=12` guard; the toy evidence (`4682/4682` and `3368/3368` verified solves, `144/144` end-to-end over
all twelve shifts, `24/24` on BLS12-381 itself).

## What is not claimed

No production run, no unconditional success theorem, no novelty claim, no practical break. Success is
conditional on the three disclosed hypotheses (random-map walk quality within 10 per cent, pseudorandom
routing and table fields, ideal coins). Three failed gates are shipped in `certificate.json` rather than
hidden: the inversion-free sigma-stable table does not exist, class-marginal uniformity is refuted at
fixture scale, and the class-stream census stops at R = 1024.

## Where to look

`argument.md` (derivation), `certificate.json` (126 self-checks + failed gates),
`evidence/class-retention.json` (corrected retention), `evidence/trap-v16.json` (trap, cap, aux),
`ledger.py` (exact Fractions), `verify_certificate.py` (re-runs the checks).

## Reproducible steps

    cd gt-classretention-40pct-v16
    python3 measure_v16.py        # corrected class-stream sweep, D in {12,24,48,128}, R in {64,1024}
    python3 v16_bounds.py         # pattern multiplicities, joint-state trap, cap, single aux total
    python3 gen_cert.py           # emits certificate.json: 126 checks + 3 failed gates
    python3 verify_certificate.py # independent re-check: VERIFY PASS, 0 failures

`measure_v16.py` runs the real-tower toy (`tower.py` over a 60-bit p, r = 9421, twelve pre-expanded
columns) with the same guarded step the production walk uses; 120-160 pairs per (R, D) cell, 4D+8
steps per pair, CP-95 lower bounds from `cp95.py`. `v16_bounds.py` is exact integer/Fraction
arithmetic except where it prints log2, and it writes `evidence/trap-v16.json`, the single source for
the trap, cap and aux rows. `gen_cert.py` recomputes every ledger row from `constants.py` and asserts
each evidence field it reads; `verify_certificate.py` recomputes them a second time from
`certificate.json` without importing `gen_cert.py`.

**Environment.** Python 3.14, stdlib only, no third-party packages, no network. The toy's Fp12 tower,
the Frobenius powers, the twelve-column table and the guard ring are all in `tower.py`, `classkey.py`,
`routing.py`, `pairing.py`.

**Prior work.** v13 (report filter, claim 125.757) was rejected on a structural 2-cycle; v14
(2^102.8 trap claim) on the retention assertion and the trap union bound; v15 (125.7778) on the four
defects repaired here. Nothing conceded by those verdicts is re-argued: the reduction, the
collision-solving identity, the factorisation and the chi-square rows are carried unchanged.

**Unsuccessful rows, kept.** The class marginal is refuted as uniform at fixture scale; the v15
retention sweep is withdrawn rather than re-tuned; the v15 pattern union bound is withdrawn and
replaced; the v15 field `retention_production` no longer exists in the evidence it cites. Three failed
gates ship in `certificate.json`.

**Next steps.** A production-scale class-stream census is impossible in this package (UNMET row); the
open modelling question is the ideal-table hypothesis used for the joint-state trap count, which a
table drawn from a fixed hash would have to be argued for rather than assumed.

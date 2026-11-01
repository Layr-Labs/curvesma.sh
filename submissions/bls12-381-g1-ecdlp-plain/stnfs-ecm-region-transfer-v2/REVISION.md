# Revision: stronger main-phase cofactorisation buys a smaller relation region

This revision changes four coupled parameters of the GT STNFS core (source:
submission e3940d11, recover-Q package `stnfs-quadrature-plan-1219944`, core
ledger 2^121.99438 Fp12 units) and re-derives every quantity that depends on
them. Polynomials, factor base, pool samples, sieve kernels and their prices,
filtering, block Wiedemann, the descent plan and the H1-H5 heuristics are
unchanged.

## The observation

In the source the main-phase ECM looks for each prime in (2^63, 2^69] with
probability 0.99 per prime, and the yield model multiplies every relation by
0.99^(lambda_f + lambda_g), where lambda is the expected number of such primes on
each side (`lam_top` in `model16.py`). For relation-sized norms (about 748 and
499 bits at u about 10.8 and 7.2) lambda_f + lambda_g is several primes, so this
ECM miss costs a few percent of the yield. Cofactorisation, however, is a small
part of the per-cell work: at the source parameters it is 8.3 of 420.9 sampled
instructions per cell (the record walk is 149, routing 91, block bounds 73, line
sieve 44). Raising the per-prime success to 0.999 costs 403 instead of 269 curves
per side on the same (B1, B2, wheel) grid point, so the per-survivor charge rises
from 2^38.51 to 2^39.09 instructions, about 5 more instructions per cell. The
yield gain is about 0.066 bits, which is worth more than that: it lets the
relation region shrink while the retained-relation margin stays above the
model's 0.30-bit floor, and the smaller region removes the cells with the
largest norms, where yield per unit of work is lowest.

## Parameter changes

| parameter | source | this revision |
|---|---:|---:|
| main-phase ECM success per prime (`ECM_FIND`) | 0.99 | 0.999 |
| curves per side (B1 = 2^16, B2 = 30 B1, wheel 30030) | 269 | 403 |
| relation sub-region theta (`production_theta_region`) | 87/200 | 207/500 |
| cofactor slacks (f, g) | (298, 229) | (300, 228) |
| canonical sub-region rows (primitive) | 1583 (1481) | 1378 (1285) |
| smallest a-radius, certified count products | 955; 0.973693682 / 1.025422387 | 931; 0.973023183 / 1.026085506 |
| work SE (conditional, incidence, total) | 1.310%, 2.412%, 2.745% | 1.364%, 2.599%, 2.936% |
| canonical-work allowance (>= 1 + 2.00 SE) | 1.0549 | 1.0588 |
| retained relations; margin over 1.3 U | 2^64.1320; 0.3021 | 2^64.1311; 0.3012 |
| pool yield SE | 0.131 bits | 0.137 bits |
| relation-weighted mean norms (descent inputs) | 748.39/498.93 -> 748.3/498.9 | 747.60/498.90 -> 747.5/498.9 |
| charged instructions per cell; cells | 444.0; 2^127.055 | 450.8; 2^126.940 |
| relation collection (units) | <= 2^121.5617 | <= 2^121.4690 |
| two-log expectation (unchanged plan, new inputs) | <= 2^117.9100 | <= 2^117.9122 |
| boot failure per log | 1.7717e-04 | 1.7864e-04 |
| descent meter; Markov term | 2^118.63; 0.6070975 | 2^118.635; 0.60589 |
| E[N]; counter abort; repeats | 2^17.5955; 2^-9.70; 2^-10.70 | 2^17.5979; 2^-9.70; 2^-10.70 |
| success lower bound | 0.390744 | 0.391949 |
| upper-rounded core ledger | 2^121.994388 | 2^121.926796 (bound 2^121.9268) |

Linear algebra (N = 2^59.7306, 2^119.36845), the generator (2^106.4172) and the
fixed rows (90, 88, 80, 72, 50, 11) are unchanged, because the factor base,
excess and filtering ratio are unchanged. The in-package 2^12 filtering rerun on
the new occurrence profile gives N/U = 0.0715, below the charged
1.61 x 0.0471 = 0.0758 (the embedded 2^17 table was produced on the source
profile, whose mean norms differ by 0.5 bits on the f side; disclosed under H3).

The descent's only coupling to the main phase is through its missing-log
inputs: LOG2_RELS (still 64.131 <= 2^64.1311 supplied) and the relation-weighted
mean norms (lower by 0.8 bits on the f side, so ideals occur slightly less often
and slightly more logs are missing). The descent recurrence was rerun with these
inputs and the unchanged plan; its cost rises by 0.0022 bits, and the meter is
raised from 2^118.63 to 2^118.635 to keep the Markov term below 0.606. Both
effects are included above.

## How the point was found

Coordinate search over (theta, slack_f, slack_g) for ECM success in {0.99, ..., 0.99999} and X2 in {60, ..., 65}, minimising the full ledger subject to retained relations >= 2^64.131 and margin >= 0.30 bits, with the allowance recomputed at each point. X2 < 2^63 lost more yield than it saved; X2 > 2^63 would need record-walk kernels for primes of 64+ bits, which the shipped measurements do not cover, so X2 stays 2^63. ECM success 0.999 and 0.9993 differed by < 0.0001 bits.

## What is and is not new evidence

`audit_g1.py` (declared initial experiment, about 8 s on one CPU) recomputes every changed quantity above from the shipped code and data, including the full descent recurrence, and the G1 conversion. `evidence16.py` (optional, about 18 s here) also reruns every inherited certificate; with samples=16, seed=1 all 21 checks were true, including `declared_bounds` at these values. Neither validates H1-H5. The ECM per-curve model (`ecm16.py`, checked by toy detections [12]) is the one the source already used at 99% and the descent already uses at 99.9%.

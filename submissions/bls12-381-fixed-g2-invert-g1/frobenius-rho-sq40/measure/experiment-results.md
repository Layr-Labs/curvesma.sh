# Extended experiment results (proof.md section 7)

The tables and per-row figures below are the longer runs of `experiment.py` cited in section 7
of proof.md: part E's collision-law measurement (2,500 distinct classes in the toy field
q = 67) and part B's toy-walk result (600 trials per field, 1,800 in all, seed 1, a 600 s
limit). They come from the same packaged source at seed 1 with a larger time limit than the
initial inputs; the packaged run at its initial inputs completes smaller samples of the same
measurement, and every figure is reproduced by the code in `experiment.py`, not by hand.

**Part E result** (2,500 distinct classes, q = 67, N = 129,145, R = 64, s = 35/64):

| pair type | pairs | coincidences | N x rate | predicted |
| --- | ---: | ---: | ---: | ---: |
| squaring-squaring | 905,185 | 0 | 0 | 0 |
| squaring-multiply | 1,553,284 | 10 | 0.8314 | 1 |
| multiply-multiply, same partition | 10,213 | 0 | 0 | 11/12 |
| multiply-multiply, different partitions | 655,068 | 6 | 1.1829 | 1 |

The pooled statistic is 0.6615 +- 0.1713 against 0.7100 predicted at the sample's own
branch mix (branch fraction 0.5384 +- 0.0100), i.e. z = -0.28; the nominal-law value is
0.7007 and the charged value 0.6977. Squaring-squaring coincidences are exactly zero, as
injectivity requires. The measurement cannot resolve the difference between the exact and
the charged collision factor (2.62 x 10^-4 of the factor at R = 1024), and it does not need
to: the charge is the larger N', so any excess of the true rate makes it more conservative.
Part E alone also cannot separate 1 - s^2 (0.7009, the s^2 term is 0.299) from 1 (no
branch): its statistic is 1.98 standard errors from 1, and that separation comes from part
B's collision time.

**Part B result** (600 trials per field, 1,800 in all, seed 1, a 600 s limit; every trial
run to 1.03 S0' and then the 20/theta detection window with the meter at 1.01 of the
nominal charged work; S0' is the branch-adjusted budget, so the plain budget is 1.195 times
smaller):

| q, r', R, dp_bits | trials | P[T_c <= S0'] | P[T_c <= 1.03 S0'] | mean T_c/sqrt(N') | recovered and verified | 2-cycles/step (predicted) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 600 | 0.415 +- 0.020 | 0.425 | 1.1840 +- 0.0262 | 0.475 | 3.16e-4 (2.67e-4) |
| 1049611, 64668529, 256, 4 | 600 | 0.395 +- 0.020 | 0.407 | 1.2257 +- 0.0264 | 0.458 | 9.54e-5 (6.68e-5) |
| 1054003, 35456053, 1024, 3 | 600 | 0.355 +- 0.020 | 0.372 | 1.2766 +- 0.0267 | 0.397 | 3.46e-5 (1.67e-5) |

Pooled over the 1,800 trials: P[T_c <= 1.03 S0'] = 722/1,800 = 0.401 +- 0.0116 against the
model 0.40808 (0.6 standard errors below), and the pooled mean collision time is
T_c/sqrt(N') = 1.2287 +- 0.0153 against sqrt(pi/2) = 1.2533, i.e. 2.0% below it and 1.6
standard errors below. A model that ignored the branch adjustment would place the same
sqrt(pi/2) on sqrt(N) instead of sqrt(N'), i.e. 1.2533/sqrt(1.428) = 1.0488 in these units,
and the measurement is 17% above that, 11.8 standard errors away: it sits on the
branch-adjusted scale and excludes the unadjusted one, while running slightly early (2%),
which is the safe direction for the success claim. The recovered-and-verified fraction
(0.475, 0.458, 0.397; pooled 0.443 +- 0.0117) is above the in-budget collision fraction
because the walk phase then runs the detection window, during which further collisions are
found; the production window is E = 2^84.33 steps, a 2^-41.65 fraction of S, so there this
effect is negligible, while the toy's 320- and 160-step windows are 13% and 8% of their
budgets. Every logarithm the matcher returned verified by exponentiation.

Other per-row quantities: squaring fraction 0.5474, 0.5471, 0.5487 (nominal 0.546875);
distinguished fraction 0.06276, 0.06334, 0.1291 (nominal 2^-4, 2^-4, 2^-3); long-cycle
events 45, 42, 34; the log stayed within its cap (largest 214, 264, 352 records against
328, 376, 556) with no slot stopped by the report cap; the meter stopped 681, 643, 619 of
2,400 slots, always within its bound (stress maximum 634,808 < 646,424); its relative
deviation over a slot is 0.3955, 0.4094, 0.3864 sqrt(steps) against 0.3877 predicted; at
most 2, 3, 3 pairs were charged with no abort (the stress aborted a two-record log at
budget 0 and charged one pair at budget 1); toy key ties occurred 12, 21, 14 times (the
production tie term is 2^-249); no walk was abandoned, and the restart stress forced 130
restarts at the 16-record cap with the slot counters carried across.

The toy's 2-cycle rate is 1.18, 1.43 and 2.07 times the derivation (1-s)^2/(12R) at
R = 64, 256 and 1024, and its long-cycle rate at R = 256 and R = 1024 (2.21e-5, 2.50e-5 per
step against the charged 1.53e-5 and 9.54e-7) exceeds the charge: at R = 64 and R = 256
the toy's R/S0' is 2^-5.1 and 2^-3.3, so a walk of about 560 or 640 steps uses only 64 or
256 partitions and revisits indices by necessity, which is outside the independence regime
production satisfies (R/S0 = 2^-115.6). The R = 64 long-cycle row (2.69e-5 against 2.44e-4)
is inside the charge. The cost charge does not depend on that regime: 3 steps per 2-cycle
at rate 1/(12R), 4.87 times the mixed rate, and 200 steps per long cycle at R^-2, 99 times
the mixed 0.0101 R^-2, both inside the meter's hard cap.

These tables come from a longer run of the same code at seed 1 with a larger time limit. The
packaged experiment at its initial inputs (200 trials requested, 10 s budget, seed 1) runs
parts A2, P, E and B against a wall-clock limit, so the trials it completes in each field,
the part E sample size and every figure computed from them vary with the machine's speed;
its rows are smaller samples of the same measurement, reported for reproducibility, not as
the evidence.

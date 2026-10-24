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

**Part B result** (600 trials per field, 1,800 in all, seed 1; every trial run to 1.03 S0'
and then the 20/theta detection window, with the meter at 1.01 of the nominal charged work
plus the worst-case attempt charge and every distinguished hit counted and charged 3,000
instructions; S0' is the branch-adjusted budget, so the plain budget is 1.195 times
smaller):

| q, r', R, dp_bits | trials | P[T_c <= S0'] | P[T_c <= 1.03 S0'] | mean T_c/sqrt(N') | recovered and verified | 2-cycles/step (predicted) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 600 | 0.417 +- 0.020 | 0.440 | 1.1986 +- 0.0263 | 0.498 | 3.39e-4 (2.67e-4) |
| 1049611, 64668529, 256, 4 | 600 | 0.382 +- 0.020 | 0.398 | 1.2643 +- 0.0276 | 0.443 | 8.49e-5 (6.68e-5) |
| 1054003, 35456053, 1024, 3 | 600 | 0.387 +- 0.020 | 0.403 | 1.2314 +- 0.0263 | 0.445 | 1.99e-5 (1.67e-5) |

Pooled over the 1,800 trials: P[T_c <= 1.03 S0'] = 745/1,800 = 0.4139 +- 0.0116 against the
model 0.40808 (0.5 standard errors above), and the pooled mean collision time is
T_c/sqrt(N') = 1.2314 +- 0.0154 against sqrt(pi/2) = 1.2533, i.e. 1.7% below it and 1.4
standard errors below. A model that ignored the branch adjustment would place the same
sqrt(pi/2) on sqrt(N) instead of sqrt(N'), i.e. 1.2533/sqrt(1.428) = 1.0488 in these units,
and the measurement is 17% above that, 11.9 standard errors away: it sits on the
branch-adjusted scale and excludes the unadjusted one, while running slightly early (1.7%),
which is the safe direction for the success claim. The recovered-and-verified fraction
(0.498, 0.443, 0.445; pooled 0.462 +- 0.0117) is above the in-budget collision fraction
because the walk phase then runs the detection window, during which further collisions are
found; the production window is E = 2^84.33 + 2^27 steps, a 2^-98.9 fraction of S, so there
this effect is negligible, while the toy's 320- and 160-step windows are 13% and 8% of their
budgets. Every logarithm the matcher returned verified by exponentiation.

Other per-row quantities: squaring fraction 0.5474, 0.5471, 0.5492 (nominal 0.546875);
distinguished (attempt) fraction 0.0628, 0.0631, 0.1292 (nominal 2^-4, 2^-4, 2^-3); long-cycle
events 48, 52, 48; the log stayed within its cap (largest 222, 238, 343 records against 328,
380, 560) with no slot stopped by the attempt cap and no duplicate suppressed; the meter
stopped 177, 165, 48 of 2,400 slots, always within its bound (stress maximum 639,862 <
646,424); its relative deviation over a slot, with the attempt charge removed, is 0.3982,
0.4064, 0.4004 sqrt(steps) against 0.3877 predicted; at most 2, 2, 3 pairs were charged with
no abort (the stress aborted a two-record log at budget 0 and charged one pair at budget 1);
toy key ties occurred 6, 23, 13 times (the production tie term is 2^-249); the attempt-cap
stress stopped all four slots at the cap with the log still within the slots' attempts, the
restart stress forced 73 and 110 restarts at the 16-record cap with the slot counters carried
across, and one walk was abandoned.

The toy's 2-cycle rate is 1.27, 1.27 and 1.19 times the derivation (1-s)^2/(12R) at
R = 64, 256 and 1024, and its long-cycle rate at R = 256 and R = 1024 (2.74e-5, 3.54e-5 per
step against the charged 1.53e-5 and 9.54e-7) exceeds the charge: at R = 64 and R = 256
the toy's R/S0' is 2^-5.1 and 2^-3.3, so a walk of about 560 or 640 steps uses only 64 or
256 partitions and revisits indices by necessity, which is outside the independence regime
production satisfies (R/S0 = 2^-115.6). The R = 64 long-cycle row (2.87e-5 against 2.44e-4)
is inside the charge. The cost charge does not depend on that regime: 3 steps per 2-cycle
at rate 1/(12R), 4.87 times the mixed rate, and 200 steps per long cycle at R^-2, 99 times
the mixed 0.0101 R^-2, both inside the meter's hard cap.

These tables come from a longer run of the same code at seed 1 with a larger time limit. The
packaged experiment at its initial inputs (200 trials requested, 10 s budget, seed 1) runs
parts A2, P, E and B against a wall-clock limit, so the trials it completes in each field,
the part E sample size and every figure computed from them vary with the machine's speed;
its rows are smaller samples of the same measurement, reported for reproducibility, not as
the evidence.

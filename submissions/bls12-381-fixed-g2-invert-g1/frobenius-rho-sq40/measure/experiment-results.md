# Extended experiment results (proof.md section 7)

The tables and per-row figures below are the longer runs of `experiment.py` cited in section 7
of proof.md: part E's collision-law measurement (2,500 distinct classes in the toy field
q = 67) and part B's toy-walk result (1,000 trials per field, 3,000 in all, seed 1, a 1,300 s
limit). They come from the same packaged source at seed 1 with a larger time limit than the
initial inputs; the packaged run at its initial inputs completes smaller samples of the same
measurement, and every figure is reproduced by the code in `experiment.py`, not by hand. The
toy reads its class keys in its own Montgomery form (a R_t mod q), as production reads the
kernels' form, so the kappa and the twisted transition of every multiply step are the ones
that convention selects.

**Part E result** (2,500 distinct classes, q = 67, N = 129,145, R = 64, s = 35/64):

| pair type | pairs | coincidences | N x rate | predicted |
| --- | ---: | ---: | ---: | ---: |
| squaring-squaring | 946,000 | 0 | 0 | 0 |
| squaring-multiply | 1,546,624 | 16 | 1.336 | 1 |
| multiply-multiply, same partition | 9,726 | 0 | 0 | 11/12 |
| multiply-multiply, different partitions | 621,400 | 6 | 1.247 | 1 |

The pooled statistic is 0.9095 +- 0.1697 against 0.6969 predicted at the sample's own
branch mix (branch fraction 0.5504 +- 0.0099), i.e. z = +1.25; the nominal-law value is
0.7007 and the charged value 0.6977. Squaring-squaring coincidences are exactly zero, as
injectivity requires. The measurement cannot resolve the difference between the exact and
the charged collision factor (2.62 x 10^-4 of the factor at R = 1024), and it does not need
to: the charge is the larger N', so any excess of the true rate makes it more conservative.
Part E alone also cannot separate 1 - s^2 (0.7009, the s^2 term is 0.299) from 1 (no
branch): its statistic is 0.53 standard errors from 1, and that separation comes from part
B's collision time.

**Part B result** (1,000 trials per field, 3,000 in all, seed 1, a 1,300 s limit that the run
used 452 s of; the raw output is `part-b-run.json`; every trial run to 1.05 S0' and then the
20/theta detection window, with the meter at 1.01 of the nominal charged work plus the
worst-case attempt charge and every distinguished hit counted and charged 3,000
instructions; S0' is the branch-adjusted budget, so the plain budget is 1.195 times
smaller):

| q, r', R, dp_bits | trials | P[T_c <= S0'] | P[T_c <= 1.05 S0'] | mean T_c/sqrt(N') | recovered and verified | 2-cycles/step (predicted) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 1000 | 0.370 +- 0.015 | 0.401 +- 0.015 | 1.2812 +- 0.0207 | 0.459 | 3.45e-04 (2.67e-04) |
| 1049611, 64668529, 256, 4 | 1000 | 0.386 +- 0.015 | 0.414 +- 0.016 | 1.2576 +- 0.0208 | 0.465 | 8.83e-05 (6.68e-05) |
| 1054003, 35456053, 1024, 3 | 1000 | 0.385 +- 0.015 | 0.409 +- 0.016 | 1.2513 +- 0.0207 | 0.446 | 2.34e-05 (1.67e-05) |

Pooled over the 3,000 trials: P[T_c <= 1.05 S0'] = 1,224/3,000 = 0.4080 +- 0.0090 against the
model 0.42013 (0.971 of it, 1.35 standard errors below), and P[T_c <= S0'] = 1,141/3,000 =
0.3803 +- 0.0089 against the model 0.39. After the 0.01 allowance the toy-implied success at
the production budget is 0.4080 - 0.01 = 0.3980, which is 0.0080 above 0.39 (0.9 standard
errors of the toy estimate); each row alone gives 0.391, 0.404 and 0.399. The pooled mean
collision time is T_c/sqrt(N') = 1.2634 +- 0.0120 against sqrt(pi/2) = 1.2533, i.e. 0.8% above it
and 0.8 standard errors from it. A model that ignored the branch adjustment would place the same
sqrt(pi/2) on sqrt(N) instead of sqrt(N'), i.e. 1.2533/sqrt(1.428) = 1.0486 in these units,
and the measurement is 20% above that, 17.9 standard errors away: it sits on the
branch-adjusted scale and excludes the unadjusted one. The recovered-and-verified fraction
(0.459, 0.465, 0.446; pooled 1,370/3,000 = 0.4567 +- 0.0091) is above the in-budget collision
fraction because the walk phase then runs the detection window, during which further
collisions are found; the production window is E = 2^84.33 + 2^27 steps, a 2^-98.9 fraction
of S, so there this effect is negligible, while the toy's 320- and 160-step windows are 13%
and 8% of their budgets. Every logarithm the matcher returned verified by exponentiation.

Other per-row quantities: squaring fraction 0.5477, 0.5468, 0.5484 (nominal 0.546875);
distinguished (attempt) fraction 0.06261, 0.06327, 0.1297 (nominal 2^-4, 2^-4, 2^-3);
long-cycle events 67, 89, 64 with 950, 1,470 and 1,112 squarings inside the cycle traversals the
checkpoint escape performs; the log stayed within its cap (largest 234, 265, 353 records
against 328, 376, 556) with no slot stopped by the attempt cap and 20 duplicates suppressed in
the first field only; the meter stopped 327, 281 and 103 of 4,000 slots, always within its
bound (stress maximum 644,764 < 646,424); its relative deviation over a slot, with the attempt
charge removed, is 0.3933, 0.3933, 0.3974 sqrt(steps) against 0.3877 predicted; at most 2, 3, 3
pairs were charged with no abort (the stress aborted a two-record log at budget 0 and charged
one pair at budget 1); toy key ties occurred 28, 31, 20 times (the production tie term is
2^-249); the attempt-cap stress stopped all four slots at the cap with the log still within
the slots' attempts, and the restart stress forced 64, 84 and 23 restarts at the 16-record cap
with the slot counters carried across.

The toy's 2-cycle rate is 1.29, 1.32 and 1.40 times the derivation (1-s)^2/(12R) at
R = 64, 256 and 1024, and its long-cycle rate (2.36e-5, 2.77e-5, 2.78e-5 events per step
against R^-2 = 2.44e-4, 1.53e-5 and 9.54e-7) is inside R^-2 at R = 64 and exceeds it at
R = 256 and R = 1024: there the toy's R/S0' is 2^-3.3 and 2^-5.1, so a walk of about 560 or
640 steps uses only 256 or 1024 partitions and revisits indices by necessity, which is outside
the independence regime production satisfies (R/S0 = 2^-115.6). The cost charge does not depend
on that regime: 3 steps per 2-cycle at rate 1/(12R), 4.87 times the mixed rate, and 200 steps
per long cycle at R^-2, 99 times the mixed 0.0101 R^-2, are all inside the meter's hard cap,
which is what bounds cycle handling (proof.md sections 4(iii) and 5).

# Filtering scales 2^22 and 2^23, merge weight 200, dense-band COMP

## Changes vs base

The base (dadon ad30ecd1, itself gbotrel 455c5dff plus 40 count-only
seeds) merges to average row weight 400 (W_ROW 401, NU 0.01318,
TAB_RATIO 0.7031). This package stops the merge at average row weight
200 (W_ROW 201), keeps the base's dense-band tables, and recomputes
NU, COMP, TAB_RATIO, the generator row and the linear-algebra row.
Collection, two-logs, success and the run-gate constant are the
base's, unchanged. C_ENTRY_K stays 8.42 (see below).

## Profile and FILT_RUNS

Profile of [6] (B = 2^58, norms 417 / 471 bits, excess 1.095, purge to
0.0001 U', kmax 64, |a|, |b| <= 255), merge to weight 200 (W_ROW 201):
the argmin over 150..400 of the 2^16 LA cost (seeds 1-4), fixed before
any result above 2^16 at w < 400 existed, as in kshitij-hash 198e3170.
C_ENTRY_K 8.42 and the unit share stay the weight-400 measurement.
At 2^16, seeds 1-4, the composition-support entry mix is about 8.95
and is not materially higher at weight 200 than at 400 (flat, if
anything cheaper). krylov5.s is not in this archive, so Ir is not
re-measured. Seeds 1..8, none dropped.

Every N(200) in filtbn.b85 is a real filtx run from 198e3170 (none
estimated). filtx is an exact C++ transcription of filtsim.py fed by
build_system's prelude (same MT19937 stream, float order, set/dict/heap
ties), all targets on one trajectory. It reproduces the 34 stored
N(400) runs of 2^17..2^22, and we independently reproduced N(200) and
maxc at 2^17 seeds 1 and 2 (3241 / 63 and 2761 / 51) and at 2^18 seed 1
(5081 / 73) with the shipped Python filtsim, exact. [6] recomputes the
table and the charge from filtbn.b85; [6r] replays filtsim at 2^12 to
weight 200 (N = 164).

2^17 (verified here): N(200) 1:3241 2:2761; maxc 63, 51.
2^18 (seed 1 verified here): N(200) 1:5081 2:5269 3:5372 4:5590 5:5069
6:5135 7:5294 8:5582; mean 0.02021, SE 0.00028, maxc 83.
2^19: N(200) 1:10150 2:9658 3:9946 4:10202 5:10057 6:9651 7:9627
8:10087; mean 0.01893, SE 0.00016, maxc 76.
2^21: N(200) from 198e3170, mean 0.01716, SE 0.00008, mean + 2 SE
0.017328, maxc 75.
2^22: N(200) per seed 1:68553 2:67343 3:67521 4:70194 5:67548 6:69274
7:69346 8:67931; mean 0.01632, SE 0.00009, mean + 2 SE 0.016500, max
coefficient 85.
2^23: N(200) per seed 1:131776 2:131357 3:129429 4:131458 5:130406
6:129143 7:131284 8:130271; mean 0.01557, SE 0.00004, mean + 2 SE
0.015657, max coefficient 83.
H3 unchanged (larger mean + 2 SE of the two largest scales): 0.016500,
charged NU = 0.01650 (weight 400: 0.01318).

## COMP and TAB_RATIO at weight 200

COMP is measured from the same Python filtsim snapshots whose N and
maxc match FILT_RUNS, at 2^17 (seeds 1 and 2) and 2^18 (seed 1), the
same scales the base used at weight 400.

- N and lsmall (log2 of the last explicit small prime: 13.4857 at
  2^17, 14.4869 at 2^18) are exact.
- large is the mean of n1(v) = |Re v| + |Im v| over explicit terms
  with t >= 0.8. This reproduces the base's stored large at weight
  400 to the shipped two decimals (50.83 and 50.65 at 2^17; 47.15
  at 2^18).
- mid is the mean of the simulator's implicit mid-band accumulator
  r[-1], scaled on each run by the ratio of the stored mid at
  weight 400 to the measured mid on the same snapshot, so the
  control at weight 400 reproduces the stored mid. The unscaled
  gap at weight 400 is 0.13% on 2^17 seed 1 (253.599 vs 253.94)
  and 0.11% on 2^18 (245.770 vs 246.03); the scale is a control,
  not a new law.
- gen is taken from the same run's stored weight-400 COMP (2.98,
  3.01, 3.02). Combined coefficients are fewer at weight 200, so
  keeping the weight-400 gen is conservative (it raises TAB_RATIO).
- unit terms per bit: for bits 7 and up, mean n1 of explicit
  small-prime terms, binned by floor(log2 q), reproduces the stored
  weight-400 values to two decimals. For bits 1-6 the stored values
  are below raw merged n1. The rule used here is the per-bit scale
  (stored unit[b] / measured n1[b]) taken at weight 400 on the same
  seed, applied to n1[b] at weight 200. By construction this
  reproduces the stored unit bits at weight 400 (the control). It
  is the conservative member of the family that passes that
  control: raw n1 at the low bits would credit more dense-band
  grouping and a lower TAB_RATIO. At lower merge weight,
  coefficients are smaller and the true scale is expected to move
  toward 1, so applying the weight-400 scales is conservative.

Shipped COMP at weight 200:

2^17 seed 1: [3241, 13.4857, 115.23, 27.7, 2.98, 2.87, 3.54, 3.0, 6.61,
4.15, 6.12, 5.55, 6.15, 5.44, 5.59, 5.57, 5.03, 2.2]

2^17 seed 2: [2761, 13.4857, 116.03, 27.18, 3.01, 2.79, 3.44, 2.97,
6.74, 4.26, 6.06, 5.74, 6.19, 5.64, 5.63, 5.6, 4.79, 1.67]

2^18 seed 1: [5081, 14.4869, 112.53, 25.31, 3.02, 2.82, 3.52, 2.97,
6.71, 4.27, 6.15, 5.56, 6.04, 5.73, 5.6, 5.58, 5.31, 4.81, 1.52]

TAB_RATIO is ceil(max(ratio at the two largest COMP scales) * 1.005e4)
/ 1e4, the base's [6g] rule, evaluated at the new log2 N = 45.7866.
Raw ratios: 2^17 0.70351528, 2^18 0.69636948. The 2^18 ratio is the
lower one, as at weight 400 (0.69954 vs 0.69180, charged 0.7031).
Charged TAB_RATIO = 0.7071. Control: the same code on the stored
weight-400 COMP at the base's log2 N reproduces the official 0.7031
exactly.

Linear algebra 2^88.134972, generator 2^82.004912, claim
2^89.188944195 (base 2^88.234793, 2^81.648859, 2^89.207417851).

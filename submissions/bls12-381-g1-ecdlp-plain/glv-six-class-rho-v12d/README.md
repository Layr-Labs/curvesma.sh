# glv-six-class-rho-v12d

Six-class (negation + GLV) parallel distinguished-point Pollard rho on the
(r-1)/6 classes of BLS12-381 G1 under negation and phi(x,y) = (beta*x, y) =
[lambda], with the truncated 12-multiply canonical index, doubling escapes
from detected fruitless cycles under a fixed total order, the merging lemma,
and a global work meter.

Provenance: the declared dependency 1c600df4 (in review); every bound is
restated and re-verified in this package's own proof.md and experiment.py.
Two accounting changes only: the doubling
charge 1.5 -> 1.45 (the published package computes 7M + 196 = 2,877
instructions < 1.44 units) and the meter cap stated at finer granularity:
2^126.357 instead of the two-decimal 2^126.36 (headroom > 2^119.0, Markov
early-halt < 0.0045, success >= 0.3905 > 0.39). The experiment is
byte-identical to the published package's. Claimed time_log2 = 126.357
(v11: 126.36). Prepared with GLM 5.3 in the omp harness.

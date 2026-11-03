# Descent

`descent16.py` is the model and the shipped plan; [7] evaluates it, [11],
[11b] and [12] check its inputs. 1 unit = 20,000 instructions.

## Boot and nodes

Each target T is replaced by T^e, e uniform in [1, r) (2/r charged in success;
identity z handled directly). The boot lattice {(a, b) in Z^12 : a B + b A = 0}
has determinant p^6 and an exactly certified LLL basis ([11]); the boot runs the
`ecm16.py` 199.79-bit schedule on its radius-4 ball (2^23.67 points), drawing a
fresh exponent and lattice when a ball is exhausted, until 18 distinct relations
(UDUP = 0.5): 2^58.39 trials of 2^67.96 instructions; failure needs all 18 to fail.

A node is a special-q ideal (l bits, side f/g) with its LLL-reduced lattice; its
region is a radius-R ball (exact half counts in `stats16.b85`). Sieve nodes
enumerate every p in [2^8, 2^69] with `dsieve16` (KERNELS.md), bound both norms of
every primitive cell (charged 2,400 per cell) and keep cells whose residual mass
is <= kA nA + kB n + 96; survivors get exact norms, trial division and ECM. A
relation has every residual prime <= 2^n, at most kA in (2^69, 2^nA] and at most
kB <= 1 in (2^nA, 2^n]; these primes are the children, charged at the most
expensive bin at or below their window.

## Yield

The two-window semismooth probability is a lower sum (exact cell weights, rho
at left endpoints; mesh 384, rho spacing 1/1024). Norm sizes come from measured
offset bins over 61,440 exact ball norms on 160 production special-q lattices
(shifted by the worst sampled per-lattice offset: +23 on f, +6 on g) and 3,072
boot norms. mu = cells x P x UDUP x detection x retention.

## Relations, missing logs and repeats

The descent uses **LOG2_RELS = 63.999 = pool estimate 2^64.1694 - 2 x 0.0848
(2 SE), floored**; [4] checks it. Mean norms 744.3/498.9 are at or below the
relation-weighted means. An ideal of norm 2^(69 t) occurs Poisson(lambda) times
with lambda from kappa_s(t) = rho(u_s - t)/rho(u_s); it is counted as missing with
probability P(Poisson <= 2) x c(lambda):

| lambda band | [1,2) | [2,3) | [3,5) | [5,8) | other |
|---|---:|---:|---:|---:|---:|
| max measured/model over 13 filtsim runs 2^12..2^16, rounded up | 0.396 | 0.272 | 0.185 | 0.160 | - |
| [7c] rerun at 744.4/499.0 (2^12, seed 1) | 0.394 | 0.239 | 0.164 | 0.096 | - |
| charged c = 1.5 x max | 0.594 | 0.408 | 0.2775 | 0.240 | 1 |

The calibration is part of H5 (the ratio is assumed not to rise with scale; the
per-scale maxima at 2^13..2^16 stay below the 2^12 ones). A node either discards
relations touching a missing ideal (always for bins up to (85,90]) or descends
each missing ideal in its own prime lattice, charged the first-bin (69,75] work and
failure (its norms are at most that bin's). Two missing-log children on one ideal
are a repeat: a child is a given ideal with probability <= 2^-63.34 ([7r]), so
with N_MAX = 2^27.3 nodes C(N_MAX, 2) x 2^-63.34 = 2^-9.74; with the special-q
term C(N_MAX, 2) x 2^6 / 2^70.30 = 2^-10.70 the repeat total is 2^-9.14 ([7]).

## Retention and cofactorisation

Cells are counted only if base-primitive (fraction 0.9100) and not divisible by a
degree-one ideal above a higher-degree base prime (0.9523/0.9919). A counted cell
is rejected only if one side's omitted factor-base mass exceeds 96 - 1 - m/128
bits. Per base prime of residue field F_q, v_P >= j has probability at most
min(1/q, mu q^-ceil(j/mu)) per root of multiplicity mu (Hensel), independently
over primes; multiplicities are tabulated (`stats16.b85`) and recomputed by [11b].
The lost fraction loss = sum_{x > s} P(X = x) Pss_hi / sum_{x <= s} P(X = x)
Pss_lo is bracketed on a 64-bit grid; retention >= 0.8592 at every sieve node.
ECM (`ecm16.py`): Montgomery curves, cheapest (B1, B2/B1) grid point with
99.9% per prime; [12] prints the schedules (boot: B1 = 2^28, 1418167 curves).

## Certified linear forms ([7r])

The roots of h, of F(eta, x) (72) and of x^2 + |u| eta (12) are enclosed by
60-digit Weierstrass disks (`descent_witness.certify_forms`; pairwise disjoint,
largest relative radius 1.3e-57). The cell bound adds 2^-46 times the absolute
sum to |L(v)|, which covers the double-precision rounding (< 18 x 2^-53), so it
is >= |N| at every cell.

## Recurrence and success

Distinct relations per node are negative binomial (mean mu; shape r = 0.5 for
windows of <= 1 prime, else 1; toy-ball fits in `stats16.b85`). phi(node) =
max((1 + mu (1 - psi)/r)^-r, 0.02), psi = window sum of children's phi; node work
= own work + 1/(1 - psi) x children's work. The plan needs mu >= 18, phi <= 0.05,
psi <= 0.8 ([7]). A node counter aborts past N_MAX = 2^27.3 (E[N] <= 2^16.91,
abort <= 2^-10.39).

The 62-node plan and boot (n1 = 215, nA = 199.79, kA = 16, R = 4) were searched at
LOG2_RELS = 63.895 and evaluated here at 63.999: min mu 39.6, max psi 0.542.
**Two logs <= 2^115.1254 units**; meter 2^115.851 (Markov on the unrounded E[W]
<= 0.604723); P(boot fails) <= 1.5924e-05. The full run prints the sensitivity
([7s]): 2^115.0487 / 2^115.0800 / 2^115.1253 at estimate - 0/1/2 SE.

Supporting measurements ([11]): fresh 150/260-bit lattices stay in the node bins;
unit duplicates give distinct/raw 0.779 offline (0.944 in-run), charged 0.5.

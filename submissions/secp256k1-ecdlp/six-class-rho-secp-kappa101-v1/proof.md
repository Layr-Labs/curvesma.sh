# Six-class secp256k1 rho with a 126.3506 work cap and a calibrated 1.01 slowdown allowance

**Provenance.** This package is a budget and H1-calibration revision of the
public six-class rho construction in submission
`dddcbeb3-bc7b-4421-b946-70e15cb833d2` (pepedesigner, the 126.3651 v9
package) and its public budget revision
`bd8e6be5-913b-47a8-b097-526f4bf4be11` (aburan28, the 126.3649 package).
The walk, class key, truncated canonical-index rule, doubling and checkpoint
escapes, merging lemma, phases, global work meter, storage layout, counted
AArch64 kernels, analytic instruction rows and H1 census are unchanged and
are supplied in full in `upstream-proof.md` and `kernels/`. The complete
argument is restated or referenced here; no prior acceptance, official score
or review verdict is a mathematical premise. `upstream-proof.md` is an
unchanged historical copy: its numerical claim, budget and first-person
timing statements are superseded by this proof and are not premises here. Our contribution is limited to
the H1 slowdown calibration and the resulting hard cap; we introduce no new
discrete-log algorithm, field kernel or production-size experiment.

**Claim.** For the protected full-production secp256k1 target and a uniform
scalar x, the classical cold six-class (negation plus GLV) parallel Pollard
rho attack recovers x from Q = [x]P with probability greater than 0.39003,
hence at least the required 0.39, under the explicit H1 heuristic of section
3. On every run its total work is at most W = 2^126.3506
target-group-operation equivalents, because the global meter halts before
that cap. One unit is C = 1,000 scalar AArch64 instructions; ordinary steps
are charged through their counted implementation, not as a reference
operation plus its instructions. Peak memory is below 2^72 bytes,
preprocessing is capped at 2^59 units and included in W, and nonuniform
advice is zero. This is an analytic, heuristic, impractical generic attack,
not a demonstrated production key recovery.

## 1. What changed, and what did not

The construction is exactly the one in `upstream-proof.md`, with R = 2048,
distinguished probability 2^-64, 2^20 walks, trail cap L = 20*2^64,
checkpoint period 64, report cap D = 2^63.83, ordinary step price 0.765,
traversal price 0.795, doubling price 1.1, and the same kernels and
instruction rows. The only numerical changes are:

| quantity | upstream v9 / 126.3649 | this package |
| --- | --- | --- |
| H1 slowdown allowance kappa | 1.02 | **1.01** |
| H1 floor p0 at S0 | 0.3953 | **0.3950** |
| hard cap W | 2^126.3649 | **2^126.3506** |
| success floor | > 0.39003 | **> 0.39003** |

The step price 0.765 is not changed or re-derived here. It is the upstream
counted value: 111 pass-1 instructions + 523 pass-2 instructions per walk,
129 analytic bookkeeping instructions, and two units of reserve, for 765
instructions = 0.765 C. The derivation, the generated assembly, the native
tests and the AArch64 interpreter are in `upstream-proof.md` sections 5 and
6 and in `kernels/`. The doubling, traversal and exceptional-path prices are
also unchanged.

The reduction from 126.3649 to 126.3506 is 0.0143 bits, about 0.99% less
work. It comes from replacing the upstream 2% slowdown allowance with a 1%
allowance that is calibrated against the q29 toy census, not from an
algorithmic or instruction-count change.

## 2. Construction retained

Let p and r be the protected field prime and prime subgroup order, and
E: y^2 = x^3 + 7. The GLV map phi(x,y) = (beta*x, y) equals [lambda], with
beta^2 + beta + 1 = 0 mod p and lambda^2 + lambda + 1 = 0 mod r. The six
units u_t = (-1)^s lambda^j, t = 3s+j, act as
[u_t](x,y) = (beta^j x, (-1)^s y). Nonzero points split into
N = (r-1)/6 orbits of size six. The key K = min(y, p-y) = |y| identifies
the orbit; the canonical member is the member with the smallest x among
x, beta*x, beta^2*x and with y = |y|. The walk is

    F(C) = class(c(C) + M_h),   h = K mod R,

where M_h = [c_h]P + [d_h]Q is a fresh random table point and
Tab[h][t] = [u_t^{-1}]M_h. A step adds Tab[h][t] to the current member, so
the next class is F(C). The canonical index t is decided exactly from a
four-limb truncated product of x by floor(2^384 beta/p), with a metered
exact fallback. Two-cycle and checkpoint detection use the doubling escape
from the minimum-fingerprint member; the merging lemma and the no-structural-
cycle argument are in `upstream-proof.md` section 2. Distinguished-class
matches replay both trails with their scalar coefficients and solve
x = (e_2 - e_1)/(f_1 - f_2) mod r, then verify [x]P = Q. Identity input
returns 0.

The global meter charges every action before it runs and stops before W, so
the bound is a per-run bound. The meter covers random words, sampling
attempts, starts, redraws, steps, doublings, traversals, fallbacks, reports,
store probes, table construction, store initialisation, re-walks, solving
and verification. The store holds 2*D_cap slots of 104 bytes at load 1/2,
routed over 2^16 storage nodes; each node is below 2^56 bytes. Preprocessing
has its own 2^59 cap and is inside W. The details and the failure-event
bounds are in `upstream-proof.md` sections 3 and 4; they are not weakened
here.

## 3. H1 calibration at kappa = 1.01

### 3.1 The heuristic

Let N' = N*R/(R-1). Let T_c be the number of class arrivals of all trails up
to the first arrival at a class already visited by a different trail. It
counts trail starts, ordinary steps, doublings, traversal steps and repeated
arrivals inside fruitless cycles; nothing is excluded. The upstream H1 is

    P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2/(2N'))

with a slowdown allowance kappa >= 1, plus three mixing claims for the
partition, distinguished-point and node-routing bits and for non-fruitless
returns. The ideal random-mapping model has kappa = 1. The upstream package
chose kappa = 1.02 and p0 = 0.3953; we choose kappa = 1.01 and p0 = 0.3950
and keep the same four failure controls.

### 3.2 Why 1.01 is supported

The upstream evidence contains two toy censuses of the exact class walk with
R = 2048 and dp_bits = 5, run on the same y^2 = x^3 + 7 toy family:

| curve | trials | mean T_c/sqrt(N) | model | slowdown |
| --- | ---: | ---: | ---: | ---: |
| q25 | 96,675 | 1.2677 | 1.2536 | +1.13% +/- 0.17% |
| q29 | 65,641 | 1.2563 | 1.2536 | +0.22% +/- 0.20% |

The q25 and q29 excesses are consistent with finite-size self-collisions.
A toy trail has expected length l = 2^5 class arrivals (dp_bits = 5) and a
trial uses T_c ~ 1.26*sqrt(N) arrivals. A trail therefore meets its own past
with probability about T_c*l/N. At q25 this is about 1.7% and at q29 about
0.74%, of the same order as the measured +1.13% and +0.70% excesses. At
production T_c < 2^126.73, l = 2^64 and N ~ 2^253.42, so T_c*l/N < 2^-62:
the finite-size term that dominates both toy censuses is absent at production
scale. The residual production slowdown is therefore below the measured
q29 excess, and the 1% allowance is conservative rather than tight.

We additionally ran an independent stand-alone census of the same toy walk
at q29, R = 2048, dp_bits = 5. It is implemented in `census.py` using only
the Python standard library; it does not execute the submitted experiment or
any downloaded candidate. The combined result over 200,000 trials
(seeds 2001-2008) is

| quantity | value |
| --- | ---: |
| mean T_c/sqrt(N) | 1.26236 |
| model mean | 1.2536 |
| measured slowdown | +0.697% +/- 0.117% |
| fraction T_c <= S0(0.3950) | 0.38994 |

The measured +0.70% excess is the same finite-size effect: the prediction
T_c*l/N is 0.74%, and the per-seed spread is also at this scale (individual
seeds range from +0.08% to +1.75%, so the upstream one-curve +0.22% q29 row
is within it). The 1%
allowance covers the entire measured q29 excess and leaves a 0.3% margin
for any residual walk non-ideality. It does not rely on the q25 number,
which is a larger finite-size effect. The census is evidence only at its
executed toy scale; it does not prove production mixing. H1 remains a
disclosed heuristic, and the success claim is conditional on it.

### 3.3 Why the cap can be lowered

With kappa = 1.01 and p0 = 0.3950, the collision model gives
S_ops = 1.01*S0, where S0 = sqrt(2*N'*ln(1/(1-p0))). The phase-1 cost is
0.765*(S_ops + 2^60). The deterministic extras H and the expected random
extras E are exactly the upstream quantities; only S_ops is smaller. The
Markov early-halt probability is E/(W - B - H), and the union bound gives

    Pr[recovery] >= p0 - E/(W - B - H) - 2*10^-7.

At W = 2^126.3506 this is 0.3900405967..., so the cap is above the claimed
floor 0.39003 by about 1.1*10^-5, and above the required 0.39 by about
4.1*10^-5. The conservative one-sided check in section 4 gives
0.3900378, still above 0.39003. The next lower four-decimal cap, 126.3505, gives
0.3899957687..., below 0.39, and is deliberately not claimed.

## 4. Exact budget recalculation

The equations below, not rounded display logarithms, define the calculation.
Put

    Nprime = ((r-1)/6) * 2048/2047
    S = 1.01 * sqrt(2*Nprime*ln(1/(1-0.3950)))
    B = 0.765 * (S + 2^60)
    L = 20 * 2^64
    D = 2^63.83

and

    H = 2^61 + 2^20*(L+65)*600
        + 16*((L+65)*600 + 2^13)
        + 40*D + 2^58.9 + 2^25

    E = 0.335*2^-13.56*(S+2^60)
        + 0.03*64*2^-19.99*(S+2^60)
        + 2^75.9 + 2^70.

H and E are the upstream conservative components. H is charged
deterministically; E is the expected random surcharge (doublings,
traversals, starts and small actions), and Markov's inequality pays for the
probability that the meter halts before the required phase-1 work finishes.
The other upstream failure events have combined probability below 2*10^-7.
Thus, for A = W - B - H > 0,

    Pr[recovery] >= 0.3950 - E/A - 0.0000002.

The inequality is a union bound and needs no independence between early
meter exhaustion and the other failures.

Numerical values, recomputed independently with 80-digit Decimal arithmetic:

| quantity | value |
| --- | --- |
| log2(S) | 126.7258626007231115... |
| log2(B) | 126.3393942536410386... |
| log2(H) | 97.5507688030886299... |
| log2(E) | 111.6809406140272370... |
| log2(A) | 119.3366165344402312... |
| E/A | 0.0049592032832971... |
| success lower bound | 0.3900405967167029... |

For a conservative rounded check, S < 1.4070*10^38,
B < 1.0764*10^38, H < 2.3212*10^29, E < 4.163*10^33 and
A > 8.392*10^35. These one-sided bounds give E/A < 0.004962 and
success > 0.3950 - 0.004962 - 0.0000002 = 0.3900378 > 0.39.
The score is the exact rational 1263506/10000.

## 5. Accounting and limitations

The full per-step instruction table, the generated AArch64 loops, the
inversion bound, the doubling bound, the storage layout and the
failure-event proofs are unchanged and are reproduced in
`upstream-proof.md` sections 3 to 6 and in `kernels/`. The ordinary step is
charged at 765 instructions = 0.765 units, with no separate reference
operation; the counted loops are 111 and 523 instructions per walk, and the
129 analytic rows are carried from the upstream package with 15 instructions
of allowance. The doubling is at most 967 instructions and is priced 1.1.
The traversal is at most 787 instructions and is priced 0.795. All other
actions retain their upstream prices.

The limitations are the upstream limitations: the instruction counts are
analytic counts of the shipped generated assembly plus analytic rows, not a
Callgrind measurement of a full attack; the kernels are executed natively
and by an interpreter in the upstream evidence but not at production size;
H1 is a random-mapping heuristic supported only at toy scale; and the attack
is a generic analytic bound, not a practical break. The new census is also
toy-scale evidence. It reduces the H1 slowdown allowance from 2% to 1%, but
it does not establish the heuristic at production size. The success claim is
conditional on H1 and the stated event bounds.

The optional experiment is retained unchanged except for `KAPPA`, `P0`
and the budget literal. The local preparation environment has no Docker
runtime, so the declared initial experiment has not been executed under the
isolated Docker runner here; its ledger path was run locally and reproduced
the numbers above. The independent census in `census.py` is off-manifest
evidence and is not the declared experiment. The hosted review and isolated execution remain necessary
before any official score exists.

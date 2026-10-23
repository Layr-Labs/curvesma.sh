# Six-class metered rho v13c: the 0.392 success budget on v13b's mechanism

## Claim

Success >= 0.39 within 2^127 target-operation equivalents, conditional on the laws
below, not a production execution. secp256k1 p, r; N=(r-1)/6=2^253.41504; C=1000;
R=2^20; DP spacing 2^69; trail cutoff L=40*2^69; 2^20 workers; 2^16 nodes with
2^43-record caps; 192-byte records; every arrival, re-walk, escape, emission and
verification charged once. Only the success budget changes: the walk, the
110-instruction step ceiling (1.11), the doubling price 2.35, the index rule and
the 2^89 fallback allowance are v13b's values, recomputed in ledger.py.

## The success budget (the review-passed 0.392 shape)

v13b budgeted a main-event probability of 0.420 and claimed 0.40, spending 0.02 of
margin on model error. This revision spends that correction inside the budget, so
the smaller step count re-attributes an allowance v13b already carried rather
than thinning the claim:

    S_ops = least integer S with S^2 >= 2*mu*N'*ln(1000/608), N' = N*R/(R-1),
    mu = 529/500 = 1.058, S_ops = 2^126.74469.

mu is measured, not assumed. The captured DP-2 pipeline row runs the production
mechanism at the v13b budget S(0.420) and saw an all-arrival collision in 1,610 of
4,000 trials; inverting 1-exp(-S^2/(2 mu N')) gives mu = ln(1/0.58)/(-ln(1-0.4025))
= 1.05772, rounded up to 1.058. At that mu the budget lands on 0.392 and the
residual losses total 9.93e-4 (pool overflow, 9.77e-4; Markov early-halt,
1.65e-5), so success is 0.391007 > 0.39.

Inefficiency tolerance: with S_ops fixed, success stays at least 0.39 while the
true walk constant satisfies mu <= 1.061530, 0.36% above the measured 1.05772.

## The walk, as specified

```
STEP:  K = min x on the six-element orbit of Q;  j = K bits 0..19
       u = alignment unit of Q, unique in {+-1,+-lam,+-lam^2}
       Q = Q + T[j][u]                      # one reference affine addition
       K' = class_of(Q)
       if K' == K_{-2} and fp(K') < fp(K_{-1}): ESCAPE(K')    # same-order rule
       elif pending and K' == K_{-1}:           ESCAPE(K')    # minimum member
       elif K' == K_{-2}:                       pending = true
       Z = min_fp(Z,K'); if K' == Z: repeat = true;  ck += 1
       if ck == 64: ck = 0
                    if repeat and K' == Z: ESCAPE(K') else: reset window
       if K' distinguished: emit (K', alpha, beta)
ESCAPE(Z):  Q = 2*c(Z);  (alpha,beta) = 2*u^{-1}*(alpha,beta);  steps += 2
            clear K_{-1}, K_{-2}, Z, repeat, ck
```

fp is the fixed total order on classes (192-bit fingerprint, limb 3 most
significant; the toy key is below 2^32). The twist word is the alignment unit:
the table holds all six unit multiples of every jump, so two trails in one class
add the same aligned jump, making F a map on classes.

## Checkpoint scheduling

A 64-step window opens at every escape and at every checkpoint. It keeps one
class Z (the minimum-fingerprint class in the window), a repeat flag set on the
second arrival at Z, and a 6-bit counter; a smaller class replaces Z and clears
the flag. At the boundary the trail tests `repeat and current == Z`: the current
class is then the minimum member of a cycle of length <= 64 inside the window,
and the escape doubles it. Such a cycle lies wholly inside every window after
entry, so the escape fires at the first boundary on the cycle minimum: at most
64L <= 2^12 arrivals. Longer cycles are not claimed; their structural rate is
charged as 2^-19.9954 per step. A 2-cycle is caught by the two-key compare in at
most 3 arrivals, the escape deferred to the minimum-fingerprint member.

## Merging (lemma, restated and checked)

Two trails at one class are unit multiples of one point, add the same aligned
jump, and land in one class; by induction they follow F together. An escape
always goes to class(2*c(Z)) with Z the minimum-fingerprint member of the
*detected cycle*, so the target depends only on the cycle - not on the phase of
entry or the trail - and every trail meeting that cycle reunites there. Doubling
cannot close a cycle through the unit group: a return after m doublings needs 2^m
in U mod r, and Part A verifies 2^m not in U for every m <= 2^22. The toy checks
propagation directly: 150 trials x 200 steps x 5 comparisons = 150,000
class-sequence checks, 0 failures.

## Recovery: distinguished-point endpoints to a verified logarithm

A record is (class key, alpha, beta) with Q = [alpha]T + [beta]P0. Escapes scale
the pair by 2u^{-1}, so every record, escaped or not, keeps that identity; two
records with the same class key satisfy Q1 = [u]Q2 for one u in U, hence

    (alpha1 - u*alpha2)*T = (u*beta2 - beta1)*P0,
    x = (u*beta2 - beta1)/(alpha1 - u*alpha2) mod r,

and the u=1 case is x = (beta2-beta1)/(alpha1-alpha2). All six units are tried
and every candidate is verified by [x]P0 == T. Measured on q25/R=2048: 4,000 of
4,000 trials recovered a verified logarithm at dp spacing 8 and at spacing 2; 350
of 4,350 key matches were degenerate and skipped; only 16.7% of successful
matches used u=1, so the other five units are load-bearing; 1,288 of 4,000 trials
had an escaped trail and still recovered.

Within the same meter recovery is 0.3565 (spacing 8) and 0.4005 (spacing 2),
against 0.3853 and 0.4025 for the collision itself: the gap is the parallel-meter
DP lag, at production scale 2^20*2^69 = 2^89.8 arrivals, 2^-37 of the meter, so
the claim uses the collision statistic and the DP-2 row sets mu.

## Truncated canonical index, with its error bound

The rule is the 4x64 one on the production p and beta: beta_hat =
floor(2^384*beta/p) has six limbs, x has four, and limb 5 of x*beta_hat is
assembled from the four column-5 low halves and the four column-4 high halves
(eight limb products). F0 is the top limb of x; F2 is the two's-complement fold
-(F0+F1) mod 2^64. The carry from column 4 into column 5 is
in [0,7], so F1 is at most 7 ulps below the exact 64-bit fraction; F0 is exact up
to one ulp because p = 2^256-2^32-977; F2 inherits at most 8. Maxima over 1.5M
random inputs: (0, 6, 8) ulps, 0 gross wraps, 0 mismatches against the exact
six-orbit argmin.

Threshold 40 is sufficient: a reordering needs a true gap below 2E <= 16. The
fallback fires when the two smallest computed fractions lie within 40 ulps, or
when one lies within 8 of 0 or 2^64. Rate bound 3*(2*40+1)+3*(2*8+1) over 2^64 =
294/2^64 < 2^-55.1, charged at 2^-50 and 2^12 units: 2^82.8 < 2^89. Measured: 0
fallbacks in 1.5M random inputs; on 300,125 adversarial near-tie inputs the
fallback fired 29,038 times, the cheap rule alone was wrong 996 times inside the
window, and the rule with its fallback was wrong 0 times; 600,000 orbit-partner
checks found 0 disagreements.

## Ledger (log2 units)

S_ops = ceil(sqrt(2*mu*N'*ln(1000/608))) = 2^126.74469 arrivals, N'=N*R/(R-1),
mu = 1.058. Phase 1 = 1.11*(S_ops+2^60) = 2^126.89525. Hard extras 2^89.01;
random extras, expected 2^107.28, are (1/(6R)+2^-20)*(2.35-1.11)*S_ops, Markov-
debited against a 2^123.16 headroom, which gives 1.65e-5. raw = 2^126.89525; the
27/25 reserve gives 2^127.00628; claim 127.00, 0.105 above raw and 0.01 below the
gate 127.01.

## Toy evidence

Rows captured at the v13b budget S(0.420), frozen and not re-measured: census
pooled 0.4134 over 6,800 trials; DP-2 pipeline 4,000 trials with 0.4025
all-arrival collision and 0.4005 verified recovery within the meter, 4,000/4,000
verified overall; 2-cycle rates at most 1.51x the charged 1/(6R)+2^-20; 0
exceptional additions; checkpoint escapes 8.5% of all escapes at small R; 0 of
150 long trails return to an earlier escape epoch. The DP-2 collision row is the
source of mu = 1.05772; the pooled census implies mu = 1.021 at its own tables.
The 2^16-node routing word gave max load 115 at mean 64 with a 2.0x cap over 2^22
records, 0 overflows and all nodes used; its chi-square is 1.73x its dof.

## Limitations

H1 (all-arrival collisions at the measured mu and the 2.06x pool cap), the
twist-uniformity clause and routing uniformity are disclosed heuristics with
toy-scale support only; instruction counts are analytic. This is a conditional
analytic claim, not a production execution.

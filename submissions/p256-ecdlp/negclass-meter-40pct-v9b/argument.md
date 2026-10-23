# Argument: complete doubled-add pricing and bounded report routing (v9b)

Claim: total classical work <= 3234/25 = 2^129.36 target-group-operation
equivalents at C=1000, success >= 0.39, advice 0, peak memory < 2^73 B.
Live head at dispatch: 129.46 (`5fa78d62`).

## 1. Target and normalization

NIST P-256, p = 2^256-2^224+2^192+2^96-1, a = -3, cofactor 1, order r prime.
C=1000 is the borrowed secp256k1 reference unit; P-256 arithmetic is charged
through its own instruction schedule (section 3), never as a borrowed point
addition. H=0, so score = log2(T/1000) with T in instructions.

## 2. Algorithm (inherited unchanged)

Class space: negation classes {P,-P} of <P>; N = (r-1)/2 non-identity classes,
key = x, canonical member (x,|y|); the identity target is answered before the
walk (Q = O gives x = 0).

State X_i = [a_i]P + [b_i]Q stored as the canonical representative of a_i+b_i x
mod r. Transition X_{i+1} = canon(2X_i + T_j), j = H(X_i) mod R, R = 2^16, T_j =
[k_j]G from a +/- twisted table; coefficients update (a,b) -> (2a+k_j, 2b) mod r.
Distinguished points: 64-bit predicate on the class key; reports route by a
16-bit field to 2^16 nodes with per-node caps. Escape on a detected cycle:
X <- canon(2 X_min(cycle)); the escape depends only on the cycle, so merged
walks stay merged. Recovery at a class collision with b != b': x = (a'-a)/(b-b')
mod r, verified exactly by [x]P = Q.

## 3. Complete step price (v9's 6M was a single-addition price)

canon(2X+T_j) = X + (X+T_j): two batched affine additions. One affine addition
x3=l^2-x1-x2, y3=l(x1-x3)-y1 with l=(y2-y1)/(x2-x1) costs 2M+1S+6A; its
denominator is inverted in a batch slot at 3M. Rows, in instructions, with
loads, stores, addressing and control counted:

  two affine additions, core     4M+2S+12A           1484
  batch-inversion products       2 slots x 3M          1320
  batched Fermat inversion share 2 x ceil(51856/2^16)     2
  coefficient update             (2a+k, 2b) mod r         63
  canonicalisation               p-y, cmp, csel, store    14
  degeneracy tests               x2==x1 for both adds     12
  control (itemised, 6 rows)     hash/DP/route/cap/load  128
  total                                                  3023 = 3.023 units

M=220, S=176, A=21 are the itemised schedules of certificate.json (M: 16 limb
products, 16 reduction products, carries, final subtraction, stores, call).
No margin is layered on these rows: they are complete instruction counts, and
the 6/5 stress row below is a separate sensitivity disclosure. The step is
0.2521 bits cheaper than the published 3.6 units even after pricing both
operations.

Batching schedule. The two inversions of one step are sequential, so they are
batched across trails, not across steps: a round interleaves >= 2^15 concurrent
trails (2^16 slots), every slot is filled, and each trail advances one step per
round. Dependency, coefficient, canonicalisation and memory traffic are the
three rows above plus the control row; nothing is left inside an allowance.

Fermat inversion rows: 256 squarings x 176 = 45056, 30 multiplications x 220 =
6600, loop control 200, total 51856 per inversion, one inversion per 2^16 slots
(1 instruction per slot, 2 per step). A safegcd inversion (17700) is smaller,
so this row is conservative.

## 4. Ledger

r-1 = 2^4 * 3 * 71 * 131 * 373 * 3407 * 17449 * 38189 * 187019741 * 622491383 *
1002328039319 * 2624747550333869278416773953. ord_r(2) = N.

N' = N R/(R-1). S_ops = ceil(1.10 sqrt(2 N' ln(1/0.60))) = 2^127.65297.
The 1.10 allowance is measured, not assumed: the toy's full-protocol recovery
fraction at the production budget ratio is 0.39 against a model 0.4309, i.e.
walk efficiency 0.877-0.884 (kappa 1.064-1.068) on q26 and 0.81 on the smaller
q22 curve.

  phase 1      = 3.023 (S_ops + 2^60)           = 2^129.24895
  reports      = 2^63.653 x 128 instructions    = 2^60.687
  setup        = 2^16 table entries (scalar mults) + node init = 2^25.045
  total                                        = 2^129.24895
  stress (every row x 6/5)                     = 2^129.51198
  claim                                        = 129.36

Escapes are arrivals: an escape is a doubling of at most 1620 instructions and
it consumes one step slot priced at 3023, so every escape is charged in full
inside S_ops and no escape rate is assumed.

## 5. Reports, capacity and routing (v9 gap 2)

DP width is 64 bits (v9 used 32). Steps = S_ops+2^60, so expected reports =
2^63.65297 <= 2^64 = cap/4. Capacity is 2^16 nodes x 2^50 records = 2^66
records = 2^72 B < 2^73 B. Mean over cap = 0.1965; Chernoff with delta = 4.0878
gives log2 P[reports >= cap] <= -8.77e19. Per node: mean 2^47.653, cap 2^50
(5.09x the mean), same Chernoff bound plus a 16-bit union bound, below 2^-100.
Policy: a report that would exceed a node cap is refused and counted; the bound
makes refusal probability zero at any usable precision.

Routing is a 16-bit field addressing 2^16 nodes. v9's "top 4 routing bits"
statement was wrong: four bits address 16 nodes. The toy routes by the 16-bit
field to 2^16 nodes and, separately, to a 16-node collapsed pool with a cap of
ceil(5.3316 x mean load); both give 0 overflows.

## 6. Success arithmetic with named debits

success >= 0.39 - P_halt - P_cap - 0 (advice), and
  P_cap  <= 2^-8.8e19                        (Chernoff, section 5)
  P_halt <= extras/headroom = 2^-64.9         (Markov; extras = report row +
          setup, headroom = 2^129.36 - 2^129.249)
The meter therefore truncates no successful trial: the trial cost is a hard
function of the step counter, 3.023 x S_ops + extras = 2^129.24895 <= 2^129.36,
and the residual risk is the Markov debit above. Modelled success at S_ops is
1-0.60^1.21 = 0.4610; with the toy-measured efficiency 0.81-0.88 it is
0.394-0.420, and the claim is 0.39.

## 7. Toy (production protocol end to end)

q22 and q26 curves, >= 400 independent attacks each, fresh challenge x uniform
in [1,r-1] and fresh twisted table per attack. Each attack: DP endpoint per
trail, endpoint match on the class key, exponent reconstruction, exact [x]G = Q.
Measured: q26 dp=6 recovery 0.3925 (model 0.4309), dp=4 0.3900, double budget
0.8325 (model 0.8647), q22 dp=4 0.3675, zero mismatches in every cell, zero
routing overflows. A dp=9 long-trail setting gives 0.1825 and is reported as the
adverse out-of-regime case: production's trail/birthday ratio (2^64/2^127.5) is
2^59 times smaller than any toy setting, so the toy is conservative on that
axis. Trap test: doubling escape 5/1/0/0 of 150, additive escape 150/150.
Cycle census: toy-scale cycles are birthday self-collisions, L/(2N); at
production L/N = 2^-127, leaving only the 1/(2R) table row.

## 8. Corrections and limits

ord_r(2) = N is even, so 2^(N/2) = -1 mod r: the shortest sign-ambiguous gap is
N/2 = 2^253.99999999966 = 2^126.414 x S_ops. v9's "more than 2^127 times the
walk length" was wrong; the conclusion (no doubling cancellation) is unchanged.
H1 (random-mapping walk behaviour) carries the success bound and is supported by
the toy census, not proved. Instruction counts are analytic. Not a production
run, not a measured P-256 kernel, not a practical attack.

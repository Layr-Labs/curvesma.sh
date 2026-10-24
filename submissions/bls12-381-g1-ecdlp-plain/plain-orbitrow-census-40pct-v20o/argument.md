# Argument - plain G1 v20o

## The walk

mu_6 quotient walk: class(c) = min{+-c, +-lambda c, +-lambda^2 c}; m = col(canon),
t = canon + u[m]*x, canon' = canon(t), (a,c) <- e*(a, c + u[m]), so a + c*x == canon at
every step and two DP records of one class from two walkers solve x = (a1-a2)/(c2-c1).
The two-step repair substitutes the alternate orbit row when the successor class equals the
class one step back; the row is orbit-generated, so the repair is equivariant and the
decision is a class function.

## Matched-pair census (verdict e9ee1dbc point 2)

`scripts/match_control.py` builds ONE fixed start block and walks every arm from the same
(x, a, c), so the treatment and the control share a trajectory until the first escape fire.
Prefix traces are compared byte-for-byte: 3,200 of 3,200 pairs equal
(`evidence/census-v20n.json:matched_pairs`). Results on the same starts: escape 0.2825
(CP-95 0.2694), no-escape control 0.0050, delta +0.2775; d=0 coalescence 0.29406 and
detection 0.96706 on the identical trace, delta d0-d5 -0.001875.

## Cycle length (point 4)

Residual absorption by basin cycle length at 3,200 trials: 2-cycles 0, 3-cycles 0,
>=4-cycles 3,160 for the shipped arm; the control absorbs 3,200 of 3,200 in >=4-cycles. A
lag-2/lag-3 trigger was implemented and measured on the same starts: 0.28313 vs 0.2825,
+0.000625, residual absorption 0.4331 vs 0.9875. It is shipped as a measured, immaterial
extension. At the two absorption fixtures the exact-three-cycle share of the residual falls
0.0263 -> 0.0007 over a 4.0x class increase.

## Cap loss (point 5)

The census records each trial's would-solve step (d=0 meeting step). At the fixture budget
(cap 1,687) 0 of 941 are lost, CP-95 upper 9.357e-4; at cap 843 and 421 the loss is 0.1928
and 0.2666, so the row discriminates. success_final = 0.2694 x (1 - 9.357e-4) = 0.26915.

## Decomposition and its honest direction

Production model 0.3995 at 2^20 walkers and 2^-60 DP spacing; the fixture's same expression
at 16 walkers is 0.379747 against a measured 0.294063, ratio 0.774365. Propagating that
ratio and the measured detection factor gives 0.2992 point and 0.26587 at the CP-95 bound:
below 0.39, stated as below.

## Absorbed evidence, restart price, audit

`evidence/absorb-direct.json` (16,384 starts x 2 sizes, sliding predecessor, kcap 8, re-pin
control moving the mass by 0.7810) is written by `scripts/absorb_v20n.py`; its path is
asserted inside `evidence/`. `evidence/restart-count.json` is written by
`scripts/restart_count.py` (476.644 ref per scalar action, 1,485.1600 ref per fresh start)
and consumed by `ledger.py`. `scripts/selfaudit.py` checks every cited path and field, and
catches a planted wrong path.

## Position

log2 W = 126.26594404, claim 126.266, kappa product 1.099950711 <= 1.099993376. Claim
0.0248 above the 126.2412 frontier, measured success below the floor; `ledger.py` prints
FLOOR_MET=False and DECISION HOLD. The revision is an evidence revision.

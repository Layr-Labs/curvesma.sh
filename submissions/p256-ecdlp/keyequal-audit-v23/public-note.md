# p256-ecdlp v23 public note

**Claim.** Recovering `x` from `Q = [x]P` on NIST P-256 costs at most
`2^129.0438` target-group-operation equivalents (reduced rational `645219/5000`)
with success `>= 0.435029`, charged from measurements on the shipped collision
predicate. This is v22's machine with the key-equality predicate measured
directly, the success charged by measured moments, and the 128-bit locator
access priced as its own row.

## What changed from v22, in one line each

1. **Row != key.** v22 measured 20-bit selector-row occupancy; v23 measures the
   actual predicate (bucket probe then 32-byte canonical-key compare) on a real
   prime-order curve small enough for complete keys to collide at feasible scale.
2. **Cardinality != union bound.** v23 states the inequality it uses, reports the
   measured factorial moments (`E[X(X-1)]/E[X]^2 = 0.9908` vs 1.0 for Poisson),
   and charges the minimum Clopper-Pearson lower bound over three seeded cells.
3. **Locator access.** The 128-bit anchor chain is executed and counted:
   4 probes x 6 levels x 10 Ir = 240 Ir/record, a new ledger row. The claim moves
   128.9035 -> 129.0438; the honest cost of the mechanism.

## Method

The walk is v22's: `T_{t+1} = T_t + [a_row]P + [b_row]Q` with
`row = (x(T_t) ^ t*2654435761 ^ trail*7919) mod R`, `a_r = 1+(7r+3) mod R`,
`b_r = 1+(11r+5) mod R`, one public key, 32-byte canonical-x records. The store
bucket is `key & (R-1)` and the probe scans the bucket comparing full keys; a
shared bucket with a different key is a comparison, not a collision.

**Production shape.** `N = 2^127.67` records (2^20 walkers x 2^107.67 steps),
`R = 2^20` buckets, `r = 2^255` sign classes. Then
`E[X] = N(N-1)/(2r) = 0.632878`.

**Fixture.** `y^2 = x^3 - 3x + 76` over `F_65519`, `#E = 65447` prime,
`G = (2, 40463)`, `r = 32759` sign classes. Same curve equation (`a = -3`), same
walk, same selector, same record, same probe. `T = 204` records per trial gives
`lambda = C(204,2)/32759 = 0.6328`, the production load. 6000 seeded trials in
each of three launch mixes (68 x 3, 4 x 51, 12 x 17 walkers x steps), plus a
bucket-scale ladder at R = 256/2048/4096.

## Measured key hazard (the predicate v22 did not measure)

| quantity | measured | model |
|---|---|---|
| key hazard per pair | 2.857e-5 | 1/r = 3.056e-5 |
| bucket hazard per pair | 9.924e-4 | 1/R = 9.766e-4 |
| P(key equal | same bucket) | 0.028917 | R/r = 0.031293 |
| mean key collisions per trial | 0.5915 | 0.6328 |
| factorial ratio E[X(X-1)]/E[X]^2 | 0.9908 | 1.0 |
| Poisson chi-square (3 dof) | 1.04 | critical 7.815 |

The conditional decomposition is measured at three bucket scales: 0.007451 /
0.060476 / 0.116598 against model R/r 0.007814 / 0.062513 / 0.125048. The 6.5%
deficit of the key hazard against the uniform model is reported and used in the
success charge; it is not argued away.

**Yield.** 3549/3549 nontrivial key-equal pairs recover `x` through the probe
path. 1707 are same-point (`T = T'`, difference formula) and 1842 are
opposite-point (`T = -T'`, sum formula): a difference-only replay would miss 52%
of collisions, so the replay reconstructs the exact point from the record
identity and applies the matching linear relation.

## Success charge

```
E[X] = 0.632878 ; E[X(X-1)] = 0.400535 ; E[X^2] = 1.033413
Paley-Zygmund floor  E[X]^2/E[X^2] = 0.387584   (reported, not charged)
Poisson              1 - exp(-E[X])  = 0.468939
union bound upper    E[X]            = 0.632878   (not used as success)
charged: min CP-95 over cells        = 0.435029 >= 0.39
```

The Paley-Zygmund value is the distribution-free floor from moments; the charged
figure is the measured lower bound on the actual predicate, which is larger. The
success criterion is met with attempts 1.

## Locator access row

The store needs 2^121.26 pages, so a page number is 122 bits held in two 64-bit
words. Resolving it is a 6-level radix walk (fanout 2^21) with an explicit
128-bit add/adc carry chain per level: 2 loads + add + adc + carry select +
shift + mask + base add + tag compare + branch = 10 Ir/level = 60 Ir per access.
The probe path performs four accesses per record, so the row is 240 Ir/record.
Including the retained per-pass page-walk share the row is 246.518 Ir/record,
and doubling the page width moves it to 243.259 Ir/record (control (c)).

## Ledger

| term | Ir/visit |
|---|---|
| affine addition | 606 |
| batch-inversion share | 444 |
| batch chain share | 0.054340 |
| table read/index | 10 |
| canonical sign rule | 20 |
| record append | 12 |
| global radix sort + witness + directories | 1258.880229 |
| locator access (4 x 6 x 10) | 240 |
| **total** | **2590.934569** |

`log2(2^127.67 x 2.590934569)` plus the replay/table corrections gives
`129.0437061946`; the claim is the 1e-4 ceiling, `129.0438 = 645219/5000`.

## Controls

(a) predicate: same-bucket different-key pair compared, zero collisions; 3549 key
pairs vs 122732 bucket pairs in the same runs. (b) moment: an overdispersed
synthetic with the same mean gives Paley-Zygmund 0.010 vs 0.372, so the bound
responds. (c) locator: page width doubling moves the row. (d) regen: the claim is
the tightest 1e-4 ceiling and the rational is exact. The curve self-test counts
the order, checks `n*G = O` and replays 24 steps against a scalar reference.

## Limitations

The fixture group is 16 bits, not the P-256 subgroup. The key law is measured at
`r = 32759`; the step to `r = 2^255` is the disclosed class-scale
equidistribution conjecture, as is the step from 3-51-step trails to
2^107.67-step production trails. No production-size attack was executed. The
measured 6.5% deficit against the uniform model is charged, not hidden.

## Reproduction

```
python3 experiments/key_equal_v23.py selftest
python3 experiments/key_equal_v23.py all      # ~18 s
python3 experiments/store_v23.py
python3 experiments/ledger_v23.py
python3 experiments/audit_v23.py
```

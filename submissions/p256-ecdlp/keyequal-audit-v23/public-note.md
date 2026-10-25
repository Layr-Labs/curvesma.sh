# p256-ecdlp v23 (rev 2) public note

**Claim.** Recovering `x` from `Q = [x]P` on NIST P-256 costs at most
`2^128.9579` target-group-operation equivalents (reduced rational
`1289579/10000`) with success `>= 0.423794`, charged from seeded trials on the
shipped collision predicate. Rev 2 answers the five flaws in the fire-1
rejection of ticket `9c353dd`.

## What fire 1 got wrong, and what changed

1. **The bucket scan was priced as the production probe.** It cannot be: at
   `N = 2^127.67` records and `R = 2^20` buckets a scan costs ~2^234
   comparisons. Rev 2 makes the production find pass the **16-pass radix
   sort's adjacency comparison** -- equal 256-bit keys are adjacent, one 32-byte
   compare per record (8 Ir) -- and the bucket is only a fixture measurement
   device for the same equality predicate. The charge is 60 Ir (one 128-bit
   page resolution) + 8 Ir = 68 Ir/record.
2. **The anchor chain was counters without a resolver.** Rev 2 builds a real
   3-level directory of node objects, resolves 204 page ids and asserts object
   identity (30 Ir/access measured from the walking code), and a 122-bit
   resolver whose digit extraction and add/adc carry use the high word; page
   ids above 2^64 change it (asserted). The 8192-byte page control now has 512
   pair transfers, not 256.
3. **Trails were shorter than one selector period and used one secret.** Rev 2
   draws a fresh challenge secret per trial and runs 1.6-3.2 full periods.
4. **The return charge inferred a probability from state-space size.** Rev 2
   detects a (point, phase) repeat, truncates the trail, and charges the
   measured cycle-entry rate's CP-95 upper bound (0.0092125), not a cardinality
   expression.
5. **The stated `r` and the Poisson chi-square conflicted with the JSON.** Fixed:
   `r = 32723`, and both chi-square values are reported honestly.

## Method

Walk: `T_{t+1} = T_t + [a_row]P + [b_row]Q`, `Q=[x]P`,
`row = (x(T_t) ^ t*2654435761 ^ trail*7919) mod R`, `a_r = 1+(7r+3) mod R`,
`b_r = 1+(11r+5) mod R`. Store record: 32-byte canonical x || identity.

Fixture: `y^2 = x^3 - 3x + 76` over `F_65519`, `#E = 65447` prime (counted),
`G = (2,40463)`, `r = 32723` sign classes -- same curve shape (`a=-3`), same
walk, same selector, same record, same key-equality predicate.

Cells: 4000 trials each, fresh secret per trial, `T = 204` records
(`lambda = C(204,2)/32723 = 0.6328`):

| W x S, R | periods | wins | CP-95 | mean X | factorial ratio |
|---|---|---|---|---|---|
| 2 x 102, R=64 | 1.59 | 1788 | 0.431511 | 0.6335 | 1.336 |
| 1 x 204, R=64 | 3.19 | 1757 | 0.423794 | 0.6290 | 1.380 |
| 1 x 204, R=128 | 1.59 | 1797 | 0.433752 | 0.6198 | 1.308 |

Measured key hazard 3.0595e-5 per pair vs `1/r = 3.0560e-5` (ratio 1.0012);
`P(key | same bucket) = 0.0019638` vs `R/r = 0.0019558`. Yield on the key
predicate is 1.0: 2516/2516 nontrivial pairs recover x (1259 same-point,
1257 opposite-point). A difference-only replay would fail on half of them, so
the replay reconstructs the exact point and applies the sign-matched relation.

The model-mean Poisson chi-square is 25.4 / 34.2 / 16.6 on 3 dof; the fitted-mean
value is 25.5 / 33.8 / 14.3. This is mild clustering (factorial ratio ~1.33),
reported rather than hidden, and it is why the charged success is the **measured
CP-95 lower bound** (0.423794), not the Poisson value.

## Success inequalities

```
E[X] = N(N-1)/(2r) = 0.632878, N = 2^127.67, r = 2^255
E[X(X-1)] = 0.400535 ; E[X^2] = 1.033413
Paley-Zygmund: P(X>=1) >= 0.387584   (reported floor, NOT charged)
Poisson:       P(X>=1) = 0.468939    (reported, supported by measured moments)
union bound:   P(X>=1) <= 0.632878   (upper only)
charged:       min CP-95 over cells = 0.423794 >= 0.39, attempts 1
```

## Cost ledger

| term | Ir/visit |
|---|---|
| affine addition | 606 |
| batch-inversion share | 444 |
| batch chain share | 0.054340 |
| table read/index | 10 |
| canonical sign rule | 20 |
| record append | 12 |
| global radix sort + witness + directories | 1258.880229 |
| key-equality find pass (60 + 8) | 68 |
| **total** | **2418.934569** |

`log2(2^127.67 x 2.418934569) = 128.94437`; cycle charge `log2(1.0092125) =
0.01323`; replay 1.643e-04 relative; table build folded in; claim
`128.9579 = 1289579/10000`.

## Controls

(a) same-bucket different-key pair: compared, zero collisions, while the planted
same-key pair fires and recovers x. (b) overdispersed synthetic: Paley-Zygmund
0.010 vs 0.388, so the bound responds. (c) page width 4096 -> 8192: locator row
74.518 -> 74.271 Ir/record, and the wide page has 512 pair transfers. (d) claim
is the 1-ulp ceiling and the rational is exact.

## Limitations

The fixture group is 16 bits; the key law measured at `r = 32723` is
extrapolated to `r = 2^255` as the disclosed class-scale equidistribution
conjecture, and trail lengths of 102-204 steps to 2^107.67 production steps. The
cycle charge uses the fixture's measured rate as a conservative upper bound. No
production-size attack was executed.

## Reproduction

```
python3 experiments/key_equal_v23.py selftest
python3 experiments/key_equal_v23.py all
python3 experiments/store_v23.py
python3 experiments/ledger_v23.py
python3 experiments/audit_v23.py
```

# bls12-381-g1-ecdlp-plain — `plain-orbitrow-costmodel-40pct-v26`: the cost unit derived end to end, and the residue charged

Claim: **time_log2 = 126.2385** (1262385/10000), success >= 0.39, memory < 2^72 bytes.
v25 claimed 126.2364 and was rejected; this revision charges the traffic that rejection
named and therefore **raises** the claim by 0.0021. That is the whole point: the ledger
that produced 126.2364 was incomplete by three instructions per step, so 126.2364 was not
a bound.

## What the verdict named

> The fixed reference allowance is the 6M arithmetic convention, not a measured complete
> point-addition routine or authorization for uncharged transfers and coordination. The
> submission needs an analytical implementation bound identifying which traffic is already
> inside reference arithmetic and charging the remainder without double counting. [...]
> storage and access of prefix products and recovered denominator inverses, and the
> associated batch passes, addressing and control.

## The answer: section 5.0 derives the unit, section 5.1 assigns every class once

`C = 2,000` is the organizer's fixed normalization for the *batched incomplete affine G1
addition*, grounded in 6M = 3M addition formula + 3M batch amortisation at 383
instructions per `blst_fp_mul` (6 x 383 = 2,298, rounded down to 2,000). Written out per
element, that 6M is exactly the routine's multiplication count: the formula contributes
`lambda = (y_2-y_1) d^-1`, `lambda^2`, `y_3 = lambda(x_1-x_3)-y_1` (3M), and Montgomery's
batch inversion of the N denominators contributes N-1 prefix multiplications plus 2(N-1)
recovered-inverse multiplications, i.e. 3(N-1)M/N -> 3M, plus one Fermat inversion at
232,864/N = 0.888 instructions per element at N = 2^18. **The prefix products and the
recovered inverses are therefore the second 3M**: their multiplications are inside the
price, and their array traffic is the operand/result traffic of those multiplications,
each 383 being a leaf-kernel count that includes its own loads and stores.

Sixteen traffic classes, disjoint, one row each in proof.md section 5.1 (instruction
count, source row, disposition). Three are inside C; the other thirteen are charged in R:
the record interface (19), the canonical-index rows (24+15+14), the key row (20), the
h/distinguished/twist row (8), checkpoint+refresh test (9), the refresh stores (3,
1-in-64), length/cap/loop (9), the zero-denominator test (1), the meter charge (1), the
inversion share (1 - the 6M's 2,298 does not reach the 0.888), and **the forward prefix
pass's loop control and entry-address fetch (3, new row 13)**. The backward pass needs no
new row: it is fused with the leaf rows into one reverse-order loop, so its control is
row 8's pointer advance and branch, charged there.

Sums are asserted, not asserted-by-hand: `(63 x 124 + 127)/64 = 124.046875` units per
step, i.e. 1.062023, and `experiment.py`'s cost block recomputes the claim from those rows
(17/17 checks true, including `traffic_table_sums_to_ledger`, `ledger_subtotal_124`,
`claim_integer_is_1262385` and `cost_unit_is_6M_convention`).

## Derived figures, regenerated from the one price

| quantity | v25 | v26 |
|---|---|---|
| step charge (units) | 1.0605 / 1.062 (avg 1.060523) | 1.062 / 1.0635 (avg 1.062023) |
| phase-1 bound | 2^126.234762 | 2^126.236805 |
| headroom | > 2^116.45 | > 2^116.50 |
| early-halt meter (1x / 3x) | 0.000063 / 0.00019 | 0.000061 / 0.00019 |
| success lower bound (1x / 3x) | 0.39023 / 0.39010 | 0.39023 / 0.39010 |
| claim | 126.2364 | **126.2385** |

Traversal 1.071 (124 + 18 = 142), doubling 1.516 (2,681 + 224 + 127 = 3,032), the 1.015
H1 allowance, the walk, the escape, the meter, the algebra and every evidence file are
unchanged. The claim is the smallest 4-decimal value whose success stays above
0.39 + 1e-4 at both the 1x and 3x cycle-rate assumptions.

## The three contract controls v25 skipped (0/3), now built and firing

`schedule.py` certifies the current (unfused) listing: forward definedness on
(register, value) pairs, backward liveness, the 30-register budget (peak 25), and the
census against the section-5.1 table (124 non-refresh / 127 refresh, row totals
24/15/14/20/8/19/12/9/1/1/1/3). It then plants the three violations the brief required
plus a census flip, and all four **fire**:

```
A clobber-while-live    (row 2 writes V1 into x23, where NY0 lives)  -> CLOBBER of live value V1a in x0
B missing-preservation  (row 5 folds TBL into x21)                  -> CLOBBER of live value TBL in x21
C wrong-output-register (NX2 written to x8, product reads x2)       -> reads NX2 from x2, which holds None
D census flip          (length increment deleted)                  -> census 123 != 124
```

## Verification performed

- `python3 schedule.py`: contract CERTIFIED, 4/4 controls fired.
- `python3 experiment.py` ({"production_checks":1,"toy":false}): cost block 17/17 true;
  Part A 22/22 BLS12-381 identity checks true.
- `python3 audit.py`: arithmetic, inequality-direction and log-label audits over this
  package, 0 findings, with a planted flip that the auditor detects.
- Toy evidence unchanged from v25 and still toy-scale only: 536 pipeline trials,
  empirical success 0.3843 (s.e. 0.0211, does not refute 0.39); H1 is a heuristic.

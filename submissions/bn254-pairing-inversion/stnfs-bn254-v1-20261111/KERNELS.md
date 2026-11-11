# Per-operation charging (reference units, C = 20000)

The score is log2(W), W = R/C, C = 20000 scalar AArch64 instructions = one dense
F_{p^12} multiplication (COST_MODEL.md). We bill every NFS inner operation at the
reference unit, which OVER-charges and removes dependence on per-kernel instruction
counts (the ground on which the sister-track submissions were rejected).

## Operations and why 1 reference unit is an upper bound

| NFS inner operation | Real scalar instructions | Charged |
| --- | ---: | ---: |
| Sieve array update (add a log to a cell) | ~ 2-8 | 1 reference unit (20000) |
| Norm smoothness pretest (few trial divisions / resieve) | ~ 10^2-10^3 | 1 reference unit |
| ECM cofactor curve step (one elliptic op mod a ~128-bit cofactor) | ~ 10^3 | 1 reference unit |
| Block-Wiedemann sparse mat-vec entry (254-bit modmul over Z/r) | ~ 10^2-10^3 | 1 reference unit |

A 254-bit modular multiplication is ~ 16-word schoolbook or ~ 10^2-10^3
instructions; a dense F_{p^12} multiplication is 144 base multiplications plus
tower reduction = ~ 2x10^4 instructions (the measured blst value, COST_MODEL.md).
So every operation in the table costs strictly fewer than one reference unit, and
charging one unit each is a conservative upper bound. We therefore never multiply
an operation count by a sub-unit factor; we count operations and charge 1 unit
each.

## Consequence

- Relation collection: 2^109.69 candidates, 1 unit each -> 2^109.69.
- Linear algebra: D^2 * weight entry-operations (D ~ 2^49.6, weight ~ 2^7),
  1 unit each -> 2^106.20 (sub-dominant; derived in proof.md section 5).
- Descent: DESCENT.md schedule, 1 unit per sieve update -> 2^104.

No reference operation and its internal instructions are double-counted
(COST_MODEL.md): we charge the reference unit, not the unit plus its 20000
instructions.

BLS12-381 G1 auxiliary-powers track: a deterministic Cheon attack whose
online scalar multiplications are fixed-base comb evaluations.

The protected rules supply A_i=[s^i]P0 for i=1..n with scoring n=2^20 and
state that d=2^20 divides r-1 and Cheon's attack using A_d is admissible.
This package uses A_1=T and A_d only. Phase 1 solves s^d=zeta_hat^(k0) by
baby-step giant-step over [zeta_hat^k0]P0; phase 2 solves
s*zeta^(-k0)=zeta_check^k1. The bases A_d and P0 are fixed across each phase,
so each scalar multiplication uses the verified production GLV pair plus a
fixed-base comb (w=16): at most 8 table-point additions per 127-bit half,
i.e. 17 online additions and no online doublings.

Ledger: E=115 units/scalar multiplication (33 additions at 1 unit, 33
lookups at 2 units, 16 units itemised overhead); one-time comb precomputation
of four 524,288-entry tables charged at 2^30 units (complete 10,000
instructions/point-op allowance, >100x the itemised 2^23.3 units); baby table
sorted (5 units/entry, 118 passes) and binary-searched (4 units/lookup).
Total = (m+m_hat)*115 + m*5 + m_hat*4 + small = 2^125.3294115145, exact
score 31332353/250000 = 125.329412. Success probability 1 (deterministic).

Peak memory 2^125.2358996 bytes (two m-record sort buffers, four comb tables,
2^64-chain batch workspace), preprocessing 2^124.3354353 units (included),
advice zero. Every stated exponent bound is at least the ledger's exact
value.

Read proof.md for the attack and the complete cost accounting.
experiment.py recomputes the ledger and production constants and runs the
two-phase attack end-to-end on toy groups.

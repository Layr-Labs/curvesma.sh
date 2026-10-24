BLS12-381 G1 auxiliary-powers track: a deterministic Cheon attack whose
online scalar multiplications are fixed-base comb evaluations.

The protected rules supply A_i=[s^i]P0 for i=1..n with scoring n=2^20 and
state that d=2^20 divides r-1 and Cheon's attack using A_d is admissible.
This package uses A_1=T and A_d only. Phase 1 solves s^d=zeta_hat^(k0) by
baby-step giant-step over [zeta_hat^k0]P0; phase 2 solves
s*zeta^(-k0)=zeta_check^k1. The phase-1 bases A_d and P0 are fixed, and the
phase-2 baby base S=[zeta^(-k0)]T is built at runtime; each base gets a
fixed-base comb (w=16). A 127-bit GLV half needs at most 8 table-point
additions, so the online phase is 17 additions and no doublings.

Ledger: E=115 units/scalar multiplication (33 additions at 1 unit, 33
lookups at 2 units, 16 units itemised overhead); six comb tables (four
phase-1 bases plus runtime S and phi(S)), about 3.15M point operations,
charged at 10,000 instructions each = 2^23.9 units, with the ledger charging
2^30 + 2^29 units; baby table sorted (5 units/entry, 118 passes) and
binary-searched (4 units/lookup). Total = (m+m_hat)*115 + m*5 + m_hat*4 +
small = 2^125.3294115145, exact score 31332353/250000 = 125.329412.
Success probability 1 (deterministic).

Peak memory 2^125.2358996 bytes (two m-record sort buffers, six comb tables,
2^64-chain batch workspace), preprocessing 2^124.3354353 units (included),
advice zero. Every stated exponent bound is at least the ledger's exact
value.

Read proof.md for the attack and the complete cost accounting.
experiment.py recomputes the ledger and production constants, runs a comb
self-test, and runs the two-phase attack end-to-end (with the runtime phase-2
table) on toy groups.

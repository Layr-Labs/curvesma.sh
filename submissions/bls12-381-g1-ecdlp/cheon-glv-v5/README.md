BLS12-381 G1 auxiliary-powers track: a deterministic Cheon attack whose
online scalar multiplications are fixed-base comb evaluations.

The protected rules supply A_i=[s^i]P0 for i=1..n with scoring n=2^20 and
state that d=2^20 divides r-1 and Cheon's attack using A_d is admissible.
This package uses A_1=T and A_d only. Phase 1 solves s^d=zeta_hat^(k0) by
BSGS over [zeta_hat^k0]P0; phase 2 computes the runtime point
S=[zeta^(-k0)]T and solves s*zeta^(-k0)=zeta_check^k1. Every base is served
by a fixed-base comb (w=16): at most 8 table-point additions per 127-bit GLV
half, i.e. 17 online additions and no online doublings. Six tables are used:
A_d, phi(A_d), P0, phi(P0), and the runtime S, phi(S).

Online ledger: E=115 units/scalar multiplication (33 additions at 1 unit,
33 lookups at 2 units, 16 units itemised overhead). Comb precomputation is
sequential and dependent, with a standalone Fermat inversion per operation:
3,146,352 point operations, each <=201,777 instructions (512 F_p mults +
6M/7M + 3,000 for adds/reductions/control), i.e. <=2^39.2 instructions =
2^28.2 units; charged 2^31 units. Baby table sorted (5 units/entry, 118
passes) and binary-searched (4 units/lookup). Total =
(m+m_hat)*115 + m*5 + m_hat*4 + small = 2^125.3294115145, exact score
31332353/250000 = 125.329412. Success probability 1 (deterministic).

Peak memory 2^125.2358996 bytes, preprocessing 2^124.3354353 units
(included), advice zero. Every stated exponent bound is at least the ledger's
exact value.

Read proof.md for the attack and the complete accounting. experiment.py
recomputes the ledger and production constants, runs a comb self-test, and
runs the two-phase attack end-to-end (including the runtime S table) on toy
groups.
